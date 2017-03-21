#!/usr/bin/env python2
import argparse
import os
import csv
import resdk
from resdk.resources import Collection


parser = argparse.ArgumentParser(description='Upload raw data.')

parser.add_argument('-sample_sheet', type=str, help='Sample sheet', required=True)
parser.add_argument('-username', type=str, help='Username', required=True)
parser.add_argument('-password', type=str, help='Password', required=True)
parser.add_argument('-URL', type=str, help='URL', required=True)

args = parser.parse_args()

res = resdk.Resolwe(args.username, args.password, args.URL)
resdk.start_logging()

samples = {}

with open(args.sample_sheet, 'rb') as sample_sheet:
    sample_reader = csv.reader(sample_sheet, delimiter='\t')
    header = next(sample_reader)
    for row in sample_reader:
        samples[row[0]] = {col: '' for col in header}
        for i, column in enumerate(row):
            if i == 0:
                continue  # skip sample name
            samples[row[0]] [header[i]] = column


for s in samples:

    collection_name = samples[s]['COLLECTION']
    try:
        collection = res.collection.get(name=collection_name)
    except LookupError:
        collection = Collection(resolwe=res)
        collection.name = collection_name
        collection.save()

    if samples[s]['PAIRED'] == '0':
        src = samples[s]['FASTQ_R1'].split(",")
        print('Uploading data for sample: {}'.format(s))
        reads = res.run('upload-fastq-single',
                        input={'src': src},
                        collections=[collection.id])

    elif samples[s]['PAIRED'] == '1':
        src1 = samples[s]['FASTQ_R1'].strip(',').split(',')
        src2 = samples[s]['FASTQ_R2'].strip(',').split(',')

        print('Uploading data for sample: {}'.format(s))
        reads = res.run('upload-fastq-paired',
                        input={
                            'src1': src1,
                            'src2': src2},
                        collections=[collection.id])

    else:
        raise KeyError('Invalid PAIRED option.')

    reads_annotation = {
        'experiment_type': samples[s]['SEQ_TYPE'],
        'protocols': {
            'extract_protocol': samples[s]['EXTRACTION_PROTOCOL'],
            'library_prep': samples[s]['LIBRARY_CONSTRUCTION_PROTOCOL'],
	    'treatment_protocol': samples[s]['TREATMENT_PROTOCOL'],
	    'growth_protocol': samples[s]['GROWTH_PROTOCOL'],
        }
    }

    reads.descriptor_schema = 'reads'
    reads.descriptor = reads_annotation
    reads.save()

    # Get sample object
    main_sample = reads.sample

    # rename a main_sample
    main_sample.name = s
    main_sample.slug = s
    main_sample.save()

    # provide sample annotation
    sample_annotation = {
        'sample': {
            'annotator': samples[s]['ANNOTATOR'],
            'source': samples[s]['SOURCE'],
            'organism': samples[s]['ORGANISM'],
            'strain': samples[s]['STRAIN'],
            'genotype': samples[s]['GENOTYPE'],
            'molecule': samples[s]['MOLECULE'],
            'optional_char': []
        }
    }

    if samples[s]['LIBRARY_STRATEGY']:
        sample_annotation['sample']['optional_char'].append(
            'LIBRARY_STRATEGY:{}'.format(samples[s]['LIBRARY_STRATEGY']))

    if samples[s]['TISSUE']:
        sample_annotation['sample']['optional_char'].append(
            'TISSUE:{}'.format(samples[s]['TISSUE']))

    if samples[s]['AGE']:
        sample_annotation['sample']['optional_char'].append(
            'AGE:{}'.format(samples[s]['AGE']))

     if samples[s]['OTHER_CHAR_1']:
        sample_annotation['sample']['optional_char'].append(
            'OTHER_CHAR_1:{}'.format(samples[s]['OTHER_CHAR_1']))

     if samples[s]['OTHER_CHAR_2']:
        sample_annotation['sample']['optional_char'].append(
            'OTHER_CHAR_2:{}'.format(samples[s]['OTHER_CHAR_2']))

    main_sample.descriptor = sample_annotation
    main_sample.save()

    # confirm that the sample is annotated and attach it to the collection
    main_sample.confirm_is_annotated()
    collection.add_samples(main_sample)

