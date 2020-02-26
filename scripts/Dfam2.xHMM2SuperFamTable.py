#!/usr/bin/env python3

import sys

args = sys.argv

if '-h' in args:
    print('''
        usage:
            DfamHMM2SuperFamTable.py < <Dfam.hmm>

        description:
            Takes Dfam.hmm on stdin and returns a 2-column list
            on stdout with (1) element name and (2) superfamily

        ''')
    sys.exit()


name = None
superfamily = None
for line in sys.stdin:
    if line.startswith('NAME'):
        name = line.strip().split(' ')[-1]
    elif line.startswith('CC'):
        if 'SubType:' in line:
            superfamily = line.strip().split(' ')[-1]
            if superfamily == 'SubType:':
                superfamily = 'Unknown'
            if name == None or superfamily == None:
                sys.exit('Problem: a record does not have superfamily or name entry.\nNAME found: {0},Superfamily found:{1}'.format(name, superfamily))
            else:
                print(name, superfamily, sep='\t')
