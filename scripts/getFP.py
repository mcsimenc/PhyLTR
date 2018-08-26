#!/usr/bin/env python3

import sys

class GFF3_line:
	'''
	Attributes:
			field0, ... , field8 - string
			attributes - dictionary

	Methods:    str()		prints gff line
		    refreshAttrStr()	updates gff line attributes (needed if changes to attributes were made)

	kwargs is a dictionary
	'''
	def __init__(self, line=None, **kwargs):

		if line == None:
			(self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes_str) = [None]*9
			self.attributes_order = []
			self.attributes = {}
			
		else:
			(self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes_str) = line.strip().split('\t')
			self.start = int(self.start)
			self.end = int(self.end)
			attributes_list = self.attributes_str.split(';')
			self.attributes_order = [ attr.split('=')[0] for attr in attributes_list ]
			self.attributes = { attr.split('=')[0]:attr.split('=')[1] for attr in attributes_list }

		self.line_number = None

		if 'line_number' in kwargs:	
			self.line_number = kwargs['line_number']

		# rename the name attribute so it conforms to GFF3 specifications, where Name is a reserved attribute key. The version of LTRDigest makes attribute key name
		if 'name' in self.attributes:
			self.attributes['Name'] = self.attributes.pop('name')
			self.attributes_order[self.attributes_order.index('name')] = 'Name'

	def __repr__(self): # for when str() or repr() are called

		return '\t'.join( [ str(self.seqid), str(self.source), str(self.type), str(self.start), str(self.end), str(self.score), str(self.strand), str(self.phase), str(self.attributes_str) ] )

	def refreshAttrStr(self, attrOrder=None):
		'''
		If the attributes have been changed this needs to be called
		before the changes are reflected on a str() call on this object.

		If an attribute has been added it should have also been added to self.attributes_order
		'''
		self.attributes_str = ';'.join([ '='.join([ attr, self.attributes[attr] ]) for attr in self.attributes_order ])


args = sys.argv

if len(args) < 3:
	print('''
	usage:
		getFP.py <noFP.gff> <withFP.gff>

	Prints features from withFP.gff that aren't in
	noFP.gff.
	''')
	sys.exit()



TP = set()
with open(args[1]) as inFl:
	for line in inFl:
		if line.startswith('#'):
			continue
		gl = GFF3_line(line)
		if gl.type == 'LTR_retrotransposon':
			TP.add(gl.attributes['ID'])

with open(args[2]) as inFl:
	for line in inFl:
		if line.startswith('#'):
			continue
		gl = GFF3_line(line)
		if gl.type=='repeat_region':
			print('###')
		if 'Parent' in gl.attributes:
			if 'LTR_retrotransposon' in gl.attributes['Parent']:
				if gl.attributes['Parent'] not in TP:
					print(line, end='')
			elif 'repeat_region' in gl.attributes['Parent']:
				if 'LTR_retrotransposon{0}'.format(gl.attributes['Parent'].lstrip('repeat_region')) not in TP:
					print(line, end='')
			else:
				print('error: Parent', file=sys.stderr)
				sys.exit(line)

		elif 'ID' in gl.attributes:
			if 'LTR_retrotransposon' in gl.attributes['ID']:
				if gl.attributes['ID'] not in TP:
					print(line, end='')
			elif 'repeat_region' in gl.attributes['ID']:
				if 'LTR_retrotransposon{0}'.format(gl.attributes['ID'].lstrip('repeat_region')) not in TP:
					print(line, end='')
			else:
				print('error: ID', file=sys.stderr)
				sys.exit(line)
		else:
			print('error: neither', file=sys.stderr)
			sys.exit(line)
			
			
#Azfi_s4673	LTRdigest	protein_match	1921	2086	5.8e-06	+	.	Parent=LTR_retrotransposon26200;reading_frame=0;name=Retrotrans_gag
#Azfi_s4673	LTRdigest	protein_match	2999	3185	6.2e-07	+	.	Parent=LTR_retrotransposon26200;reading_frame=1;name=gag_pre-integrs
#Azfi_s4673	LTRdigest	protein_match	3233	3566	8.5e-17	+	.	Parent=LTR_retrotransposon26200;reading_frame=1;name=rve
#Azfi_s4673	LTRdigest	protein_match	4449	4596	1.6e-17	+	.	Parent=LTR_retrotransposon26200;reading_frame=2;name=RVT_2
#Azfi_s4673	LTRdigest	protein_match	4660	4867	3.5e-28	+	.	Parent=LTR_retrotransposon26200;reading_frame=0;name=RVT_2
#Azfi_s4673	LTRdigest	protein_match	4883	5132	1.6e-20	+	.	Parent=LTR_retrotransposon26200;reading_frame=1;name=RVT_2
#Azfi_s4673	getorf	ORF	5572	5872	.	+	.	ID=LTR_retrotransposon26200.ORF.01;Parent=LTR_retrotransposon26200;translated_seq=GRMVQTNSEGLGVPIMDPIRILCDNMSSIYLARNPVFHTRTKHIEVHYHFIRERVQSGEIDLQHVSTNLQVADIFTKALGIDKLGQFASGLGLTPSALPA
#Azfi_s4673	LTRharvest	long_terminal_repeat	5971	6115	.	+	.	Parent=LTR_retrotransposon26200
#Azfi_s4673	LTRharvest	target_site_duplication	6116	6119	.	+	.	Parent=repeat_region26200
####
#Azfi_s4685	LTRharvest	repeat_region	148	5434	.	+	.	ID=repeat_region26201
#Azfi_s4685	LTRharvest	target_site_duplication	148	151	.	+	.	Parent=repeat_region26201
#Azfi_s4685	LTRharvest	LTR_retrotransposon	152	5430	.	+	.	ID=LTR_retrotransposon26201;Parent=repeat_region26201;ltr_similarity=85.19;seq_number=3822;dfamClassification=None;repbaseClassification=Gypsy-4_PPa-I
#Azfi_s4685	LTRharvest	long_terminal_repeat	152	286	.	+	.	Parent=LTR_retrotransposon26201
#Azfi_s4685	LTRdigest	protein_match	395	938	1.2e-09	+	.	Parent=LTR_retrotransposon26201;reading_frame=0;name=RVT_1
#Azfi_s4685	getorf	ORF	2005	2332	.	+	.	ID=LTR_retrotransposon26201.ORF.01;Parent=LTR_retrotransposon26201;translated_seq=EAEGEVLVAGTVPDVHRFVTTCESCQMHSIVRPGRASPYVPSDHPLQVDGRSGDDADGGRTDAVSGPSAGRPDEPGGRPSPSEQDDRGRMSVPNRGCDMPVRVRREDRG
#Azfi_s4685	getorf	ORF	2555	3197	.	+	.	ID=LTR_retrotransposon26201.ORF.02;Parent=LTR_retrotransposon26201;translated_seq=PGTCRRSSCSVRSRSCRWNERSPRGHGRLEDEMSREELLAARIRQLERRPEDVETAAEKIRMARTRTKPGSTERTDSDRRRSRKVTGCSSTTAASTTSTEQRGSSRGGGSGRTRYKRQRQRDVSSGGTRRYENGDTGGRKEDQSLQEAARRRDPVSPGQQRRPIRGRGRNRWRRVKAAGFLRANLIWRMPVDARLGGADVVKKTVRGGIEHKAW
#Azfi_s4685	getorf	ORF	3286	3619	.	+	.	ID=LTR_retrotransposon26201.ORF.03;Parent=LTR_retrotransposon26201;translated_seq=MPNGLLTGYKRRGRRRQRYLARRGRSAAEPDSELLSLFLFLLLLHRVPVQVSRRLVLFVERPECVLQLLRPVLLFLLRRNKQGSHLFAALFALLVSFPCPIAGCNPCLVIC
#Azfi_s4685	getorf	ORF	3737	4481	.	+	.	ID=LTR_retrotransposon26201.ORF.04;Parent=LTR_retrotransposon26201;translated_seq=SVFGAIGYPSVFFLSLLFFGTRPPRSVSVHTRCVITPTTGPGDNRVERPGRGEDGLRRSRLPELRGTGGRQTGGETRTVRDLATIRDPVRASTKRSRRIRYAKIPQEGRVPSWDGNRRPRTAEGPRRTSPRQQGGADATGGVGNRKEGVSTTPLWTCQAEEAVDGRGTGKAGRGANDAPGRHGAGRRVGREPLFAATLGPSDHGNPGKDRDVQEPVVALVDHGSEINLMSMDFYKKGKWPINTSTGGR
#Azfi_s4685	getorf	ORF	4742	5069	.	+	.	ID=LTR_retrotransposon26201.ORF.05;Parent=LTR_retrotransposon26201;translated_seq=PRLPRFPKTERIFREGCQGVGLRQSGHPGKPHTCTIEPKKTCTTGRNFCFSDSGACTIVGIGGPREGEGGGEASSVRSSRPGRGNQRSRGGDLLRRHLGNAGQGSRRRP
#Azfi_s4685	LTRharvest	long_terminal_repeat	5315	5430	.	+	.	Parent=LTR_retrotransposon26201
#Azfi_s4685	LTRharvest	target_site_duplication	5431	5434	.	+	.	Parent=repeat_region26201
####
#Azfi_s4686	LTRharvest	repeat_region	3881	5200	.	?	.	ID=repeat_region26202
#Azfi_s4686	LTRharvest	target_site_duplication	3881	3885	.	?	.	Parent=repeat_region26202
#Azfi_s4686	LTRharvest	LTR_retrotransposon	3886	5195	.	?	.	ID=LTR_retrotransposon26202;Parent=repeat_region26202;ltr_similarity=92.97;seq_number=3823;dfamClassification=None;repbaseClassification=None
#Azfi_s4686	LTRharvest	long_terminal_repeat	3886	4013	.	?	.	Parent=LTR_retrotransposon26202
#Azfi_s4686	LTRharvest	long_terminal_repeat	5071	5195	.	?	.	Parent=LTR_retrotransposon26202
#Azfi_s4686	LTRharvest	target_site_duplication	5196	5200	.	?	.	Parent=repeat_region26202
####
#Azfi_s4687	LTRharvest	repeat_region	2348	4391	.	?	.	ID=repeat_region26203
#Azfi_s4687	LTRharvest	target_site_duplication	2348	2364	.	?	.	Parent=repeat_region26203
#Azfi_s4687	LTRharvest	LTR_retrotransposon	2365	4374	.	?	.	ID=LTR_retrotransposon26203;Parent=repeat_region26203;ltr_similarity=97.48;seq_number=3824;dfamClassification=None;repbaseClassification=None
#Azfi_s4687	LTRharvest	long_terminal_repeat	2365	2721	.	?	.	Parent=LTR_retrotransposon26203
#Azfi_s4687	LTRharvest	long_terminal_repeat	4020	4374	.	?	.	Parent=LTR_retrotransposon26203
#Azfi_s4687	LTRharvest	target_site_duplication	4375	4391	.	?	.	Parent=repeat_region26203
####
#Azfi_s4711	LTRharvest	repeat_region	20	6582	.	-	.	ID=repeat_region26204
#Azfi_s4711	LTRharvest	target_site_duplication	20	23	.	-	.	Parent=repeat_region26204
#Azfi_s4711	LTRharvest	LTR_retrotransposon	24	6578	.	-	.	ID=LTR_retrotransposon26204;Parent=repeat_region26204;ltr_similarity=94.25;seq_number=3828;dfamClassification=None;repbaseClassification=Gypsy-4_PPa-I
#Azfi_s4711	LTRharvest	long_terminal_repeat	24	249	.	-	.	Parent=LTR_retrotransposon26204
#Azfi_s4711	LTRdigest	protein_match	1235	1418	7.9e-08	-	.	Parent=LTR_retrotransposon26204;reading_frame=0;name=rve
#Azfi_s4711	getorf	ORF	1681	2296	.	-	.	ID=LTR_retrotransposon26204.ORF.01;Parent=LTR_retrotransposon26204;translated_seq=TRTLGYRFGGKGRQDYLIGTEVIIETDCLPILGMVSGCATPDLAMLMDSVHQVPGPGNPTHFRERQRHGRHAFEARFDDEGGMVSEDEEVGVDFFEAAYVTTDGTSTPALNDFDESKYDGEWLQIGGFLRTMTPDASWTKDEANRIRKKAYRFFLRDGYLWKHPKKRNGVPLRVVAMKETGGAAEGLPRQPMGRTPWDMGPRSRN
#Azfi_s4711	getorf	ORF	2326	2899	.	-	.	ID=LTR_retrotransposon26204.ORF.02;Parent=LTR_retrotransposon26204;translated_seq=TTSRSRDARWKRRTRRSGTDAGDSWRRTSTTARRCCKARGCEAHLLRGEIGLRAIGIMVVGHLCGPYGRKPSPAKVEAISAMKPTALGHGGTKFLGACAFYHIWIPHYAHVAERCTVAEERAKIRVAEHTESVRKLKEALAAAPALRKGLREGYPGVHNGRYESDRDRMGRQPRRRGRHSVSDPIRREGPQ
#Azfi_s4711	getorf	ORF	3092	3641	.	-	.	ID=LTR_retrotransposon26204.ORF.03;Parent=LTR_retrotransposon26204;translated_seq=GQLRAPLPEGSESRIRGVASDPSLRDPACIGHTFTDGTLRELKIGGGGFLLPAEETGFGECSEGTGKPSHSHRKRSVRRPDDRRTDGDIHGTACAVESQADSGPRAHIPKLMELLRQKVDMGILEPSSAPYSNRWFTVPRKTDVALHSRPAAGQQSDHSERRDRTDHRRVRGGLRRKVDLLGR
#Azfi_s4711	getorf	ORF	4126	4729	.	-	.	ID=LTR_retrotransposon26204.ORF.04;Parent=LTR_retrotransposon26204;translated_seq=KNGSSTAKWSCHLREVLGIAKKEFHDSIVDLVKRKRLSTETEPERPVEVRTTHIDDMALEDEWAESHYSRPHWARATTETPVKIGDVQEPVVALVDHGSEINLMSMDFYKKGKWPINTKHGWKIRAATRATEELHGACPNVRVKIGDVEIDQHFFVQETSSHPVILGNRTSRRREWRPRCWTTVRPTRGSKARTGAIRSNF
#Azfi_s4711	getorf	ORF	4905	5589	.	-	.	ID=LTR_retrotransposon26204.ORF.05;Parent=LTR_retrotransposon26204;translated_seq=RACGRFDKRREWADNGRPGRARSQRGSLNRRRRCEGKRPGWLETGSASIRGSTGGAALEELTRAIKDLQIAQARREGGEPARDRRGSRKQVHVEVGHIRKDCGDFAEALRNRVVYLWEGRVHASDTRRALNRNDGRGGMKRLMEEAAARHAETVHYSASAGIGSGQRGPKGEQTGFWPTMLEGLAGARLKKDEADRAEKRVREITGWSDPVEEKTGFVEAAARTTRLW
#Azfi_s4711	getorf	ORF	5614	6100	.	-	.	ID=LTR_retrotransposon26204.ORF.06;Parent=LTR_retrotransposon26204;translated_seq=TREGTGALDRAIEAVVGTIGRFSGKDATKYLASYGAEMLMRDIPEERRLAGFPRVAMPSIHAEVLEVRAESRTWEEFEGRLLEKYGLDDALRLSKRNFMEWVESPGREETRRYSSGSLRNTLRAFRRLTERSWIRAGSYCSSSRWTSGIGIKWAPCWKPKTG
#Azfi_s4711	LTRharvest	long_terminal_repeat	6355	6578	.	-	.	Parent=LTR_retrotransposon26204
#Azfi_s4711	LTRharvest	target_site_duplication	6579	6582	.	-	.	Parent=repeat_region26204
