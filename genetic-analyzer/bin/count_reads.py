#!/usr/bin/env python3

import argparse
import os
import genetic_analyzer as ga

def main():
    parser = argparse.ArgumentParser(description="count_reads")
    parser.add_argument("-settings",     dest="settings",    required=True, help="settings.txt")
    parser.add_argument("-samples",      dest="samples",     required=True, help="sample IDs, comma separated")
    parser.add_argument("-feature",      dest="feature",     required=True, help="e.g. exon")
    parser.add_argument("-attribute",    dest="attribute",   required=True, help="e.g. gene_id")
    parser.add_argument("-strand_opt",   dest="strand_opt",  required=True, help="yes/no/reverse")
    args = parser.parse_args()

    # Load settings and determine output path override
    settings = ga.load_settings(args.settings)
    samples = args.samples.split(',')
    output_override = os.getcwd() + os.sep  # ensures trailing slash

    for sample in samples:
        # Override the output path to point to the current working directory
        settings[sample]['output_path'] = output_override

        # Run analysis steps
        ga.count_reads(settings, sample, feature=args.feature, attribute=args.attribute, strand_opt=args.strand_opt)
        ga.gene_lengths(settings, sample, feature=args.feature, attribute=args.attribute)
        ga.mapped_reads(settings, sample)

if __name__ == "__main__":
    main()

