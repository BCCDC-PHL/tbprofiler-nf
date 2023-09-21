#!/usr/bin/env python3

import argparse
import csv
import sys

def main(args):
    bed = []
    current_bed_record = None
    with open(args.input, 'r') as f:
        header = f.readline()
        for line in f:
            line_split = line.strip().split('\t')
            chrom = line_split[0]
            pos = int(line_split[1]) - 1
            depth = int(line_split[3])
            if depth < args.threshold:
                if current_bed_record is None:
                    current_bed_record = {
                        'chrom': chrom,
                        'start': pos,
                        'end': pos,
                    }
            else:
                if current_bed_record is not None:
                    current_bed_record['end'] = pos
                    bed.append(current_bed_record)
                    current_bed_record = None

    writer = csv.DictWriter(sys.stdout, fieldnames=['chrom', 'start', 'end'], delimiter='\t')

    for record in bed:
        writer.writerow(record)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize depth of coverage')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-t', '--threshold', type=int, default=10, help='Threshold for depth of coverage')
    args = parser.parse_args()
    main(args)
