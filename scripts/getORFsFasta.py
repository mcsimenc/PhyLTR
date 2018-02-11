#!/usr/bin/env python3

import sys
import os

def getORFsFASTA(elementsFl, orffasta, outfasta):
	'''
	Write a new fasta with the orfs from orffasta for elements in elements.
	'''
	elements = set()
	with open(elementsFl, 'r') as inFl:
		for line in inFl:
			elements.add(line.strip())
	if os.path.isfile(outfasta):
		os.remove(outfasta)
	with open(orffasta, 'r') as inFl:
		for line in inFl:
			if line.startswith('>'):
				el = line.strip()[1:].split('.')[0]
				WRITE = False
				if el in elements:
					WRITE = True
			if WRITE:
				with open(outfasta, 'a') as outFl:
					outFl.write(line)

def help():
	print('''
	usage:
		getORFsFasta.py -elements <path> -orffasta <path> -outfasta <path>

		-elements		A one-column list of LTR RTs, like LTR_retrotransposon3
	''', file=sys.stderr)


args = sys.argv

if '-elements' not in args or '-orffasta' not in args or '-outfasta' not in args or len(args) < 7:
	help()
	sys.exit()


elementsFl = args[args.index('-elements')+1]
orffasta = args[args.index('-orffasta')+1]
outfasta = args[args.index('-outfasta')+1]


getORFsFASTA(elementsFl, orffasta, outfasta)

print('Done.', file=sys.stderr)
