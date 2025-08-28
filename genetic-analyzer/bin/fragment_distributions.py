#!/usr/bin/env python3

# Copyright (C) 2017 by Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# Released under MIT license (see LICENSE.txt)

import argparse
import os
import genetic_analyzer as ga

def main():
    parser = argparse.ArgumentParser(description="fragment_distributions")
    parser.add_argument("-settings", dest="settings", required=True, help="settings.txt")
    parser.add_argument("-samples",  dest="samples",  required=True, help="1,2")
    args = parser.parse_args()

    samples  = args.samples.split(',')
    settings = ga.load_settings(args.settings)

    # Override output paths to current working directory (Nextflow work dir)
    output_dir = os.getcwd() + os.sep
    for s in samples:
        if s in settings:
            settings[s]['output_path'] = output_dir

    for s in samples:
        ga.fragment_length_dists(settings, s, reads_to_sample=1000000)

if __name__ == "__main__":
    main()

