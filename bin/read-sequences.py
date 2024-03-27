#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import pandas as pd
from Bio import SeqIO

# Function to parse dates from sequence labels
def parse_date_from_label(label):
    match = re.search(r'(\d{4}-\d{2}-\d{2})', label)
    if match:
        return match.group(1)
    return 'N/A'

# Function to parse sequences from FASTA files
def parse_sequences(files):
    sequences = []
    unknown_labels = []
    for file in files:
        if file.name == '<stdin>':
            print("(reading from standard input)", file=sys.stderr)
        for record_idx, record in enumerate(SeqIO.parse(file, 'fasta')):
            date = parse_date_from_label(record.id)
            if date == 'N/A':
                unknown_labels.append(record.id)
            sequence = str(record.seq)
            sequences.append({
                'sequenceLabel': record.id,
                'sequenceCharacters': sequence,
                'sequenceDateISO': date,
                'sequenceDateYear': date.split('-')[0] if date != 'N/A' else 'N/A',
                'sequenceDateMonth': date.split('-')[1] if date != 'N/A' else 'N/A',
                'sequenceDateDay': date.split('-')[2] if date != 'N/A' else 'N/A',
                'sequenceLengthMin': len(sequence),
                'sequenceLengthMax': len(sequence),
                'sequenceSourcePath': file.name if file.name != '<stdin>' else 'standard input',
                'sequenceSourceOffset': record_idx,
            })
    return sequences, unknown_labels

# Function to output summaries based on arguments
def output_summaries(df, summarize_prefix, summarize_format):
    base_filename = f"{summarize_prefix}.{{}}"
    summary_functions = {
        'yearly': ['sequenceDateYear'],
        'monthly': ['sequenceDateYear', 'sequenceDateMonth'],
        'daily': ['sequenceDateYear', 'sequenceDateMonth', 'sequenceDateDay']
    }
    for summary_name, group_columns in summary_functions.items():
        summary_df = df.groupby(group_columns).size().reset_index(name='NumSequences')
        filename = base_filename.format(summary_name)
        if summarize_format == 'json':
            summary_df.to_json(filename, orient='records', lines=True)
        elif summarize_format == 'csv':
            summary_df.to_csv(filename, index=False)
        elif summarize_format == 'tsv':
            summary_df.to_csv(filename, index=False, sep='\t')

# Main function
def main():
    parser = argparse.ArgumentParser(description="Reads sequences from FASTA files and outputs them along with extracted dates.")
    parser.add_argument('paths', nargs='*', default=['-'], help="Paths to FASTA files")
    parser.add_argument('-o', '--output-filepath', type=argparse.FileType('w'), default=sys.stdout, help="Output filepath")
    parser.add_argument('-f', '--output-format', choices=['json', 'csv', 'tsv'], default='json', help="Output format")
    parser.add_argument('-s', '--summarize', type=str, help="Prefix for summary files")
    parser.add_argument('--summarize-format', choices=['json', 'csv', 'tsv'], help="Output format for summary files")
    args = parser.parse_args()

    # Open input files
    input_files = [open(file, 'r') if file != '-' else sys.stdin for file in args.paths]

    # Parse sequences
    sequences, unknown_labels = parse_sequences(input_files)
    df = pd.DataFrame(sequences)

    # Output
    if args.output_format == 'json':
        df.to_json(args.output_filepath, orient='records', lines=True)
    elif args.output_format == 'csv':
        df.to_csv(args.output_filepath, index=False)
    elif args.output_format == 'tsv':
        df.to_csv(args.output_filepath, index=False, sep='\t')

    # Summarize
    if args.summarize:
        summary_format = args.summarize_format if args.summarize_format else args.output_format
        output_summaries(df, args.summarize, summary_format)

    # Report to standard error
    print(f"Total sequences read: {len(sequences)}", file=sys.stderr)
    print(f"Total sequences parsed properly: {len(df) - len(unknown_labels)}", file=sys.stderr)
    print(f"Total unknown sequences: {len(unknown_labels)}", file=sys.stderr)
    if unknown_labels:
        print("Unknown sequence labels:", file=sys.stderr)
        for label in unknown_labels:
            print(label, file=sys.stderr)

    # Close input files
    for file in input_files:
        if file is not sys.stdin:
            file.close


if __name__ == "__main__":
    main()

