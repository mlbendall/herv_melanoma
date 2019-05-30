#! /usr/bin/env python
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import argparse

def main(input, output, sdict):
    gtf = [l.strip('\n').split('\t') for l in input if not l.startswith('#')]
    if sdict is None:
        gtf.sort(key=lambda x:int(x[3]))  
        gtf.sort(key=lambda x:x[0])
    else:
        chroms = []
        dfile = (l.strip('\n').split('\t') for l in sdict)
        for r in dfile:
            if r[0] == '@SQ':
                d = {_[:2]:_[3:] for _ in r[1:]}
                chroms.append(d['SN'])
        chromd = {v:i for i,v in enumerate(chroms)}

        gtf = [g for g in gtf if g[0] in chromd]
        gtf.sort(key=lambda x:int(x[3]))
        gtf.sort(key=lambda x:chromd[x[0]])
    
    print('\n'.join('\t'.join(g) for g in gtf), file=output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sort GTF.')
    parser.add_argument('--sdict',
                        type=argparse.FileType('r'),
                        help='''Sequence dictionary (i.e. from picard 
                                CreateSequenceDictionary) for chromosome order'''
    )
    parser.add_argument('input',
                        type=argparse.FileType('r'),
                        nargs='?',
                        default=sys.stdin,
                        help='GTF file to be sorted'
    )
    parser.add_argument('output',
                        type=argparse.FileType('w'),
                        nargs='?',
                        default=sys.stdout,
                        help='Sorted GTF file'
    )
    args = parser.parse_args()
    try:
        main(args.input, args.output, args.sdict)
    except BrokenPipeError:
        pass
