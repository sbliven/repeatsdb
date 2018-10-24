#!/usr/bin/python
"""
Fetch info about a single RepeatDB repeat

@author Spencer Bliven <spencer.bliven@gmail.com>
"""

import sys
import argparse
import logging
from repeatsdb import RepeatsDBRegion


def main(args=None):
    parser = argparse.ArgumentParser(description='Fetch info about a single RepeatDB repeat')
    parser.add_argument("pdbChain", help="PDB ID and chain (eg 1g3jA)")
    parser.add_argument("outfile", type=argparse.FileType('w'),
                        default=sys.stdout, nargs="?", help="Output file (default: stdout)")
    parser.add_argument("-n", "--repeat-region", default=0, type=int, help="Repeat region number to use")
    parser.add_argument("-v", "--verbose", help="Long messages",
                        dest="verbose", default=False, action="store_true")
    args = parser.parse_args(args)

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG if args.verbose else logging.WARN)

    try:
        regions = RepeatsDBRegion.fetchRepeatRegionsJSON(args.pdbChain)
        repeat = RepeatsDBRegion(json=regions[args.repeat_region])
        args.outfile.write(str(repeat))
        args.outfile.write("\n")
    except BrokenPipeError:
        pass


if __name__ == "__main__":
    main()
