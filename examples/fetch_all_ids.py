#!/usr/bin/python
"""
Fetch current list of all RepeatDB structures

@author Spencer Bliven <spencer.bliven@gmail.com>
"""

import sys
import os
import argparse
import logging
from repeatsdb import RepeatsDBSearch

def getRepeatsDBChains():
    """Fetch a list of all reviewed proteins from RepeatsDB

    Returns: (list of str)
    """
    return [result["id"][0] for result
            in RepeatsDBSearch().search(query="repDB_source:Reviewed",params="id")
           ]

def main(args=None):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("outfile", type=argparse.FileType('w'), 
                        default=sys.stdout, nargs="?", help="Output file (default: stdout)")
    parser.add_argument("-v", "--verbose", help="Long messages",
        dest="verbose", default=False, action="store_true")
    args = parser.parse_args(args)

    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG if args.verbose else logging.WARN)
    
    try:
        for pdbChain in getRepeatsDBChains():
            args.outfile.write(pdbChain)
            args.outfile.write("\n")
    except BrokenPipeError: pass

if __name__ == "__main__":
    main()
