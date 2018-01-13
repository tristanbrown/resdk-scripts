#!/usr/bin/env python3

"""
A script to migrate reads data from DictyExpress to Genialis' BCM Server
"""

import sys
import os
import config
import time
import resdk
from resdk.resources import Collection
from genesis import Genesis
from requests.exceptions import RequestException

def main(args=None):
    """Usage: python migrate_de.py username pw
    """
    # Login to Genesis
    try:
        user = config.de_login['user']
        pw = config.de_login['pw']
    except KeyError:
        user = input('DictyExpress username: ')
        pw = input('DictyExpress password: ')

    gen = Genesis(email=user, password=pw)

    # Login to BCM
    try:
        user = config.bcm_login['user']
        pw = config.bcm_login['pw']
    except KeyError:
        user = input('BCM username: ')
        pw = input('BCM password: ')

    res = resdk.Resolwe(user, pw, config.bcm_login['addr'])
    resdk.start_logging()

    # Create a Migrate object that will transfer a data file as a Sample object
    M = Migrate(gen, res)
    collection, datalist = M.datalist(config.de_data)
    for data in datalist:
        sample = Sample(collection, data)
        M.migrate(sample)

class Migrate():
    """
    """
    def __init__(self, server1, server2):
        self.gen = server1
        self.res = server2
        self.download = True # Set to True when running script for real. 
        self.upload = True # Set to True when running script for real.

    def datalist(self, gen_data):
        """Uses the parameters in gen_data to generate the list of data on the 
        DE server to be migrated. Returns the collection name and data list. 
        """
        data = self.gen.project_data(gen_data['col'])
        collection = self.gen.projects()[gen_data['col']]
        return (collection,
            [d for d in data if d.type.startswith(gen_data['type'])]
        )
    def download_file(self, read_file, n):
        """Downloads the DE reads file, 1 or 2, specified by n.
        """
        if n == 1:
            path = read_file.path
        elif n == 2:
            path = read_file.path2
        
        downloaded = (
            os.path.exists(path) and
            os.stat(path).st_size > 0
            )

        if not downloaded:
            print(
                'Downloading Read {0} for sample: {1}'.format(
                    n, read_file.name
                    )
                )

            time1 = time.time()
            with open(path, 'wb') as f:
                if n == 1:
                    f.write(read_file.fastq())
                elif n == 2:
                    f.write(read_file.fastq2())
            time2 = time.time()

            print(
                'Read {0} downloaded in {1} seconds'.format(
                    n, int(time2-time1)
                    )
                )

    def migrate(self, read_file):
        """Downloads a DE data file and uploads it to the collection 
        specified.
        """
        time0 = time.time()

        # Create the collection, if necessary
        coll = read_file.collection
        if coll:
            try:
                collection = self.res.collection.get(name=coll)
            except LookupError:
                collection = Collection(resolwe=self.res)
                collection.name = coll
                collection.save()
        
        # Check if the data already exists on BCM
        try:
            reads = collection.data.get(name=read_file.file)
            print('File already uploaded: {}'.format(read_file.file))
            self.download = False
            self.upload = False
        except LookupError:
            self.download = True
            self.upload = True

        # Download the data files from DictyExpress to temporary storage
        # Takes ~9 min per GB

        if self.download:
            self.download_file(read_file, 1)
            self.download_file(read_file, 2)

        # Upload the reads files to BCM
        if self.upload:
            print('Uploading reads for sample: {}'.format(read_file.name))
            time4 = time.time()
            reads = self.res.run('upload-fastq-paired',
                            input={
                                'src1': read_file.path,
                                'src2': read_file.path2})

            if coll:
                collection.add_data(reads)
            time5 = time.time()
            print('Reads uploaded in {} seconds'.format(int(time5-time4)))

        # Annotate the reads
        reads.descriptor_schema = 'reads'
        reads.descriptor = read_file.reads_annotation
        reads.save()

        # Get sample object
        main_sample = reads.sample

        # rename a main_sample
        main_sample.name = read_file.name
        main_sample.slug = read_file.slug
        main_sample.save()

        # provide sample annotation
        main_sample.descriptor_schema = 'sample'
        main_sample.descriptor = read_file.sample_annotation
        main_sample.save()

        # confirm that the sample is annotated and attach it to the collection
        main_sample.confirm_is_annotated()
        if coll:
            collection.add_samples(main_sample)

        # Delete the reads files from temporary storage
        try:
            os.remove(read_file.path)
            os.remove(read_file.path2)
        except FileNotFoundError:
            pass

        # Done
        time6 = time.time()
        print('Migration of {0} completed in {1} seconds.'.format(
                read_file.name,
                int(time6-time0)))

class Sample():
    """Create a sample-like object from a Genesis object. 
    """
    def __init__(self, coll, genobject):
        self.object = genobject
        
        an = genobject.annotation
        
        self.name = genobject.name
        self.slug = genobject.name
        self.collection = coll.name
        self.file = an['output.fastq']['value']['file']
        self.file2 = an['output.fastq2']['value']['file']
        self.path = os.path.join(config.storage, self.file)
        self.path2 = os.path.join(config.storage, self.file2)

        self.reads_annotation = {
            'experiment_type': 'Chemical mutagenesis',
            'protocols': {
                'growth_protocol': an['var.sample.growth']['value'],
                'treatment_protocol': an['var.sample.treatment']['value'],
                'extract_protocol': 'phenol/chloroform',
                'library_prep': 'DNA shearing',
                'fragmentation_method': an['var.seqrun.fragmenting']['value']
            },
            'reads_info': {
                'instrument_type': an['var.seqrun.adapter']['value'],
                'facility': an['var.seqrun.center']['value']
            }
        }

        self.sample_annotation = {
            'sample': {
                'annotator': an['var.experimenter']['value'],
                'organism': 'Dictyostelium discoideum',
                'source': an['var.experiment']['value'],
                'strain': an['var.sample.strain']['value'],
                'genotype': an['var.sample.genotype']['value'],
                'molecule': 'genomic DNA',
                'optional_char': ['{0}:{1}'.format(
                        'screen', an['var.experiment']['value']
                )]
            }
        }
    
    def fastq(self):
        """Start the download request for the first fastq file.
        """
        return self.object.download('output.fastq').content

    def fastq2(self):
        """Start the download request for the second fastq file.
        """
        return self.object.download('output.fastq2').content

if __name__ == "__main__":
    main(sys.argv[1:])