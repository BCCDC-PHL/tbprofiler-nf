#!/usr/bin/env python3

import argparse

def main(args):
    exit(-1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check outputs')
    parser.add_argument('--pipeline-outdir', type=str, help='Path to the pipeline output directory')
    args = parser.parse_args()
    main(args)
