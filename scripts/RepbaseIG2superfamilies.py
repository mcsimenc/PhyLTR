#!/usr/bin/env python3

import sys

if '-h' in sys.argv:
	print('''
	usage:
		RepbaseIG2superfamilies.py < <Repbase_extract.IG>
	
	description:

		This script is meant to be used to generate an
		ID-superfamily mapping for Repbase entries for use
		in the PhyLTR pipeline. Repbase is often updated so
		this is necessary to keep PhyLTR's classification
		abilities up-to-date.

	instructions:

		Download all ERV and LTR retrotransposon entries
		from Repbase in IG format. Repbase downloads are
		restricted to those with accounts. Accounts are
		often given freely upon request:

		http://www.girinst.org/repbase/update/browse.php

		Run this script and combine the output using append
		(>>). i.e.

		RepbaseIG2superfamilies.py < Repbase_LTR_RT.IG > Repbase.SF
		RepbaseIG2superfamilies.py < Repbase_ERV.IG >> Repbase.SF

		Then place Repbase.SF in the PhyLTR directory

		RepeatDatabases/Repbase/
		''')
	sys.exit()

name = None
superfmaily = None
GET = False
for line in sys.stdin:
	if line.startswith(';ID'):
		name = line.strip().split()[1]
		GET = True
	elif line.startswith(';KW') and GET:
		superfamily = line.strip().split()[1][:-1]
		if superfamily.startswith('Endo'):
			superfamily = 'ERV'
		elif superfamily.startswith('BEL'):
			superfamily = 'Pao'
		elif superfamily == 'LT':
			superfamily = 'Unknown'
		print(name, superfamily, sep='\t')
		GET = False
		name = None
		superfamily = None

# Example IG file
#;ID   ALTR2       DNA   ; MAM   ; 321 BP
#;XX
#;DE   Retroviral long terminal repeat.
#;XX
#;AC   .
#;XX
#;DT   09-OCT-1997 (Rel. 2.09, Created)
#;DT   18-FEB-2011 (Rel. 16.02, Last updated, Version 3)
#;XX
#;KW   ERV1; Endogenous Retrovirus; Transposable Element; LTR;
#;KW   Long terminal repeat; AMER1; ALTR2.
#;XX
#;NM   ALTR2.
#;XX
#;OS   Bos taurus
#;XX
#;OC   Bos taurus
#;OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
#;OC   Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia;
#;OC   Pecora; Bovidae; Bovinae; Bos.
#;XX
#;RN   [1]  (bases 1 to 321)
#;RA   Smit,A.F.
#;RT   ALTR2.
#;RL   Direct Submission to Repbase Update (30-NOV-1996)
#;XX
#;CC   ~81% identical to consensus.
#;CC   Putative LTR. AATAAA signal present, possibly 4 bp
#;CC   duplication.x.
#;XX
#;DR   [1] (Consensus)
#;XX
#;SQ   Sequence 321 BP; 90 A; 84 C; 50 G; 92 T; 5 other;
#ALTR2
