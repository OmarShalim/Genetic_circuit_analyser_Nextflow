#!/usr/bin/env python

# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Supporting modules
import argparse
import os
import genetic_analyzer as ga

def main():
    # Parse command line inputs
    parser = argparse.ArgumentParser(description="transcription_profile")
    parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
    parser.add_argument("-samples",   dest="samples",   required=True,  help="tube_1", metavar="string")
    parser.add_argument("-chroms",    dest="chroms",    required=True,  help="0x58v50", metavar="string")
    args = parser.parse_args()

    # Setup output to Nextflow's work dir
    output_dir = os.getcwd() + os.sep
    settings = ga.load_settings(args.settings)

    # Force all outputs to current directory
    for s in settings.keys():
        settings[s]['output_path'] = output_dir
    settings['None']['output_path'] = output_dir

    samples = args.samples.split(',')
    chroms = args.chroms.split(',')

    for s in samples:
        ga.create_profiles(settings, s, chroms, scale=10e-5)

if __name__ == "__main__":
    main()

