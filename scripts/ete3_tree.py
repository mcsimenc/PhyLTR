#!/usr/bin/env python3

import sys
from ete3 import Tree, TreeStyle, SeqMotifFace, CircleFace, faces, add_face_to_node


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


def hex_to_RGB(hex):
        ''' "#FFFFFF" -> [255,255,255] '''
# Pass 16 to the integer function for change of base
        return [int(hex[i:i+2], 16) for i in range(1,6,2)]


def RGB_to_hex(RGB):
        ''' [255,255,255] -> "#FFFFFF" '''
        # Components need to be integers for hex to make sense
        RGB = [int(x) for x in RGB]
        return "#"+"".join(["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in RGB])


def color_dict(gradient):
        ''' Takes in a list of RGB sub-lists and returns dictionary of
          colors in RGB and hex form for use in a graphing function
          defined later on '''
        return {"hex":[RGB_to_hex(RGB) for RGB in gradient], "r":[RGB[0] for RGB in gradient], "g":[RGB[1] for RGB in gradient], "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
        ''' returns a gradient list of (n) colors between
          two hex colors. start_hex and finish_hex
          should be the full six-digit color string,
          inlcuding the number sign ("#FFFFFF") '''
        # Starting and ending colors in RGB form
        s = hex_to_RGB(start_hex)
        f = hex_to_RGB(finish_hex)
        # Initilize a list of the output colors with the starting color
        RGB_list = [s]
        # Calcuate a color at each evenly spaced value of t from 1 to n
        for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
                curr_vector = [int(s[j] + (float(t)/(n-1))*(f[j]-s[j])) for j in range(3)]
                # Add it to our list of output colors
                RGB_list.append(curr_vector)

        return color_dict(RGB_list)

def polylinear_gradient(colors, n):
        ''' returns a list of colors forming linear gradients between
            all sequential pairs of colors. "n" specifies the total
            number of desired output colors '''
        # The number of colors per individual linear gradient
        n_out = int(float(n) / (len(colors) - 1))
        # returns dictionary defined by color_dict()
        gradient_dict = linear_gradient(colors[0], colors[1], n_out)

        if len(colors) > 1:
                for col in range(1, len(colors) - 1):
                        next = linear_gradient(colors[col], colors[col+1], n_out)
                        for k in ("hex", "r", "g", "b"):
                        # Exclude first point to avoid duplicates
                                gradient_dict[k] += next[k][1:]

        return gradient_dict

def gff2strands(gff=None, attr='ID', Type='LTR_retrotransposon', trim=True):
	strandDct = {}
	with open(gff, 'r') as inFl:
		for line in inFl:
			if line.startswith('#'):
				continue
			if '\t{0}\t'.format(Type) in line:
				gffLine = GFF3_line(line)
				A = gffLine.attributes[attr]
				if trim:
					A = A[len(Type):]
				strandDct[A] = gffLine.strand
	return strandDct

args = sys.argv

if len(args) < 5:

	print('''
		usage:
			ete3_tree.py -t <newick_tree> -d <divergences> -g <gff> [-reroot [<node>|auto]] [-nlabel] [-lflabel] [-annot] [-ultrametric]

		-lflabel
			show leaf labels

		-annot
			show leaf annotations (blast hits and domain architectures

		-reroot
			if <node>, the tree is rerooted at that node. If 'auto', then the filename for the input tree is parsed to get
			the name of the outgroup used and the tree is rerooted on the outgroup.	

		-ultrametric
			Run ete3's convert_to_ultrametric() on the tree
		-orfhits <path>		List of ORFs with hits to nr
		''', file=sys.stderr)

	sys.exit()

categories = {
'asp':['AP2','Asp','Asp_protease_2','Asp_protease','dUTPase_2','dUTPase','Exc','Exo_endo_phos_2','FB_lectin','Foamy_virus_ENV'],

'various':['ATHILA','BCNT','DBD_Tnp_Mut','DDE_Tnp_IS1595','DDE_Tnp_ISL3','DpnI','HTH_psq','HTH_Tnp_1','HTH_Tnp_Tc3_2','HTH_Tnp_Tc5','IN_DBD_C','Inhibitor_I34','Gypsy','Herpes_ORF11','Intron_maturas2','LEDGF','Maelstrom','Maff2','Mu-transpos_C','N-Term_TEN','Nup153','Nup_retrotrp_bd','PEN-2','Peptidase_A17','Peptidase_A2B','Peptidase_A2E','Phage_Cox','Phage_GPA','PHINT_rpt','Piwi','RAG1','RdRP_5','RE_AlwI','RE_SacI','Retro_M','Slu7','SNF5','SPP','SQAPI','SSV1_ORF_D-335','Sulfolobus_pRN','TagA','Telomerase_RBD','Thg1','TLV_coat','Tn916-Xis','Tnp_P_element_C','Tnp_zf-ribbon_2','Transposase_22','Transposase_28','Transp_Tc5_C','TYA','Vif','WCCH','Y1_Tnp','Yuri_gagarin','Zea_mays_MuDR'],

'duf':['DUF1725','DUF1759','DUF3158','DUF3173','DUF3258','DUF3435','DUF3701','DUF3806','DUF390','DUF4102','DUF4219','DUF4413'],

'gag':['gag-asp_proteas','Gag_MA','Gag_p10','Gag_p12','Gag_p17','Gag_p19','Gag_p24','Gag_p30','Gag_p6','gag_pre-integrs','Gag_spuma', 'Retrotran_gag_2','Retrotran_gag_3','Retrotrans_gag'],

'integrase':['Integrase_1','Integrase_AP2','Integrase_DNA','Integrase_Zn','Phage_integ_N','Phage_Integr_2','Phage_integr_3','Phage_integrase','Phage_int_SAM_1','Phage_int_SAM_2','Phage_int_SAM_3','Phage_int_SAM_4','Phage_int_SAM_5','rve_2','rve_3','rve'],

'recombinase':['Recombinase'],

'rnaseh':['RHSP','RNase_H2_suC','RNase_H2-Ydr279','RNaseH_C','RNase_H','RNaseH_like'],

'rvp':['RVP_2','RVP'],

'rvt':['RVT_1','RVT_2','RVT_3','RVT_connect','RVT_N','RVT_thumb'],

'zf':['zf-C2H2','zf-CCHC_2','zf-CCHC_3','zf-CCHC_4','zf-CCHC_5','zf-CCHC_6','zf-CCHC','zf-H2C2','zf-H3C2','zf-RVT']
}

tree_flpath = args[args.index('-t') + 1]
treeName = tree_flpath.split('/')[-1]
divergences_flpath = args[args.index('-d') + 1]
divergences = {}
divergencesCorrected = {}
#domains_dct = {}
IGCdct = {}
classifDct = {}

if '-reroot' in args:
	REROOT = True
	reroot_at = args[args.index('-reroot') + 1]
	if reroot_at == 'auto':
		#reroot_at = tree_flpath.split('.')[-1].split('_')[-1][15:] # Copia_0.bootstrapped.pathd8_ultrametric.outgroup_LTR_retrotransposon58
		reroot_at = tree_flpath.split('.')[-2].split('_')[-1] # Copia_cluster_0.bootstrapped.pathd8_ultrametric.outgroup_69.newick
else:
	REROOT = False

if '-orfhits' in args:
	orfhitsfl = args[args.index('-orfhits')+1]
	orfhits = { line.strip() for line in open(orfhitsfl, 'r') }
	
	
if '-lflabel' in args:
	LEAF_LABELS = True
else:
	LEAF_LABELS = False

if '-annot' in args:
	ANNOT = True
else:
	ANNOT = False
if '-ultrametric' in args:
	ULTRAMETRIC = True
else:
	ULTRAMETRIC = False

LTRRTs = {}
LTRRTlengths = {}
gffFlPth = args[args.index('-g')+1]
strandDct = gff2strands(gff=gffFlPth, attr='ID', Type='LTR_retrotransposon', trim=False)
with open(gffFlPth, 'r') as gffFl:
	for line in gffFl:
		if line.startswith('#'):
			continue
		gffLine = GFF3_line(line)
		feat = gffLine.type
		start = int(gffLine.start)
		end = int(gffLine.end)
		length = end - start + 1

		# Get get span of whole element
		if feat == 'repeat_region':
			el = 'LTR_retrotransposon{0}'.format(gffLine.attributes['ID'].split('repeat_region')[-1])
			LTRRTlengths[el] = {'start':start, 'end':end, 'length':length }
			continue
			
		# Get span of LTRs
		if feat == 'long_terminal_repeat':
			el = gffLine.attributes['Parent']
			start = start - LTRRTlengths[el]['start'] + 1
			end = end - LTRRTlengths[el]['start'] + 1
			if el in LTRRTs:
				if feat in LTRRTs[el]:
					LTRRTs[el][feat]['1'] = (start, end)
				else:
					LTRRTs[el] = {feat:{'0':(start, end)}}
			else:
				LTRRTs[el] = {feat:{'0':(start, end)}}

		# Get span of domains
		elif feat == 'protein_match':
			el = gffLine.attributes['Parent']
			domain = gffLine.attributes['Name']
			start = start - LTRRTlengths[el]['start'] + 1
			end = end - LTRRTlengths[el]['start'] + 1
			if el in LTRRTs:
				if domain in LTRRTs[el]:
					ct = 0
					while domain in LTRRTs[el]:
						ct += 1
						domain = '{0}_{1}'.format(domain, ct)
					LTRRTs[el][domain] = (start, end)
				else:
					LTRRTs[el][domain] = (start, end)
			else:
				LTRRTs[el] = {domain:(start, end)}
		elif feat == 'ORF':
			el = gffLine.attributes['Parent']
			orfid = gffLine.attributes['ID']
			if orfid in orfhits:
				domain = 'ORFHIT'
			else:
				domain = 'ORF'
			start = start - LTRRTlengths[el]['start'] + 1
			end = end - LTRRTlengths[el]['start'] + 1
			if el in LTRRTs:
				if domain in LTRRTs[el]:
					ct = 0
					while domain in LTRRTs[el]:
						ct += 1
						domain = '{0}_{1}'.format(domain, ct)
					LTRRTs[el][domain] = (start, end)
				else:
					LTRRTs[el][domain] = (start, end)
			else:
				LTRRTs[el] = {domain:(start, end)}

with open(divergences_flpath) as in_fl:
	for line in in_fl:
		if line.startswith('LTR'):

			#elementName	classification	MCLinflationValue	cluster	clusterSize	model	divergence	correctedDivergence	IntraelementGeneConversion
#name	full_alignment_length	effective_alignment_length	non_identities	proportion_non_identities	est_subs_per_site_baseml_HKY85_model	bp_between_ltrs	domain_architecture	best_blast_hit
			#fields = line.strip().split('\t')
			rt_name, classification, I, clust, clustSize, model, div, divc, IGC = line.strip().split('\t')
			#domains = fields[7]
			divergences[rt_name] = div
			divergencesCorrected[rt_name] = divc
			#domains_dct[rt_name] = domains
			IGCdct[rt_name] = IGC
			classifDct[rt_name] = classification

#num_colors = int('{0:.4f}'.format(cutoff).replace('.', ''))
num_colors = 20000
color_gradient = polylinear_gradient(['#FAFF00', '#FF1800', '#001EFF', '#000000'], num_colors) # yellow, red, blue, black
t = Tree(tree_flpath)
#for leaf in t.get_leaf_names():
#    el = 'LTR_retrotransposon{0}'.format(leaf)
#    Motifs = []
#    for feat in LTRRTs[el]:
#        if feat == 'long_terminal_repeat':
#            for LTR in LTRRTs[el][feat]:
#                motif = [ LTRRTs[el][feat][LTR][0], LTRRTs[el][feat][LTR][1], "[]", None, "black", "black" ]
#                Motifs.append(motif)
#        elif feat == 'protein_match':
#            [ LTRRTs[el][feat][0], LTRRTs[el][feat][1], "[]", None, "black", "black" ]
#            Motifs.append(motif)
#
#    seqFace = SeqMotifFace(seq=None, motifs=Motifs, gap_format="line")
#    (t & leaf).add_face(seqFace, 0, aligned)
#
    # start, end, shape, width, fgcolor, bgcolor
    #[ start, end, "[]", None, "black", "black" ]



greatest_div = {'element':None, 'value':0}
for node in t:
	node_name = str(node).split('-')[-1]
	rt_name = 'LTR_retrotransposon{0}'.format(node_name.split('_')[0])
	if rt_name in divergences:
		if greatest_div['value'] < float(divergences[rt_name]):
			greatest_div['value'] = float(divergences[rt_name])
			greatest_div['element'] = rt_name
	#	print(int(divergences[rt_name][:6].replace('.', '')))
	#	print(len(color_gradient['hex']))
		col = color_gradient['hex'][int(divergences[rt_name][:6].replace('.', ''))]
	else:
		col = '#417849'
	if ANNOT:
		try:
			classifFace = faces.TextFace(classifDct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		except KeyError:
			classifFace = faces.TextFace('ERROR', fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		try:
			IGCface = faces.TextFace(IGCdct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')
		except KeyError:
			IGCface = faces.TextFace('ERROR', fsize = 8, penwidth = 10, fgcolor = 'Black')
		#(t & node_name).add_face(domain_face, 0, 'aligned')
		(t & node_name).add_face(IGCface, 0, 'aligned')
		#(t & node_name).add_face(blast_hit_face, 1, 'branch-right')
		(t & node_name).add_face(classifFace, 1, 'branch-right')

	cf = CircleFace(radius=4, color=col)
	(t & node_name).add_face(cf, 0, 'branch-right')

        #el = 'LTR_retrotransposon{0}'.format(leaf)
	el = 'LTR_retrotransposon{0}'.format(node_name)
	Motifs = []
	for feat in LTRRTs[el]:
		if feat == 'long_terminal_repeat':
			for LTR in LTRRTs[el][feat]:
				# Add box for LTRs
				start = LTRRTs[el][feat][LTR][0]
				end = LTRRTs[el][feat][LTR][1]
				if strandDct[el] == '-': # invert the graphic
					start = int(LTRRTlengths[el]['length']) - start
					end = int(LTRRTlengths[el]['length']) - end
				motif = [ start//50, end//50, "[]", None, 8, "black", "black", None ]# start, end, typ, w, h, fg, bg, name
					#motif = [ LTRRTs[el][feat][LTR][0]//100, LTRRTs[el][feat][LTR][1]//100, "[]", None, 4, "black", "black", None ]
					#motif = [ LTRRTs[el][feat][LTR][0]//10, LTRRTs[el][feat][LTR][1]//10, "[]", None, 8, "black", "black", None ]
					#motif = [ LTRRTs[el][feat][LTR][0], LTRRTs[el][feat][LTR][1], "[]", None, 1, "black", "black", None ]
				Motifs.append(motif)
		else: # domains
			#motif = [ LTRRTs[el][feat][0], LTRRTs[el][feat][1], "[]", None, 1, "blue", "blue", None ]

			
			if feat in categories['zf'] or len([ cat for cat in categories['zf'] if feat.startswith(cat) ]) > 0:
				domainColor = "Pink" # Zinc finger is pink
			elif feat in categories['rvt'] or len([ cat for cat in categories['rvt'] if feat.startswith(cat) ]) > 0:
				domainColor = "Red" # Reverse transcriptase is red
			elif feat in categories['rvp'] or len([ cat for cat in categories['rvp'] if feat.startswith(cat) ]) > 0:
				domainColor = "DeepPink" # Retroviral protein is deeppink
			elif feat in categories['rnaseh'] or len([ cat for cat in categories['rnaseh'] if feat.startswith(cat) ]) > 0:
				domainColor = "DarkOrange" #RNaseH is Orange
			elif feat in categories['recombinase'] or len([ cat for cat in categories['recombinase'] if feat.startswith(cat) ]) > 0:
				domainColor = "Yellow" # Recombinase is yellow
			elif feat in categories['integrase'] or len([ cat for cat in categories['integrase'] if feat.startswith(cat) ]) > 0:
				domainColor = "cyan" # Integrase is Cyan
			elif feat in categories['gag'] or len([ cat for cat in categories['gag'] if feat.startswith(cat) ]) > 0:
				domainColor = "purple" # Gag is purple
			elif feat in categories['duf'] or len([ cat for cat in categories['duf'] if feat.startswith(cat) ]) > 0:
				domainColor = "ForestGreen" 
			elif feat in categories['asp'] or len([ cat for cat in categories['asp'] if feat.startswith(cat) ]) > 0:
				domainColor = "Lime"
			elif feat in categories['various'] or len([ cat for cat in categories['various'] if feat.startswith(cat) ]) > 0:
				domainColor = "silver"
			elif 'ORFHIT' in feat:
				domainColor = "LightSeaGreen"
			elif 'ORF' in feat:
				domainColor =  "BurlyWood"
			else:
				sys.exit('feat not in categories: {0}'.format(feat))

			start = LTRRTs[el][feat][0]
			end = LTRRTs[el][feat][1]
			if strandDct[el] == '-': # invert the graphic
				start = int(LTRRTlengths[el]['length']) - start
				end = int(LTRRTlengths[el]['length']) - end
			if domainColor == "silver":
				motif = [ start//50, end//50, "[]", None, 8, domainColor, domainColor, "arial|6|black|{0}".format(feat) ]
			else:
				#motif = [ LTRRTs[el][feat][0], LTRRTs[el][feat][1], "[]", None, 4, domainColor, domainColor, None ]
				motif = [ start//50, end//50, "[]", None, 8, domainColor, domainColor, None ]
				#motif = [ LTRRTs[el][feat][0]//100, LTRRTs[el][feat][1]//100, "[]", None, 4, domainColor, domainColor, None ]
			Motifs.append(motif)

	#box_motifs = [
	#	# seq.start, seq.end, shape, width, height, fgcolor, bgcolor
	#	[0,  5, "[]", None, 10, "black", "rgradient:blue", "arial|8|white|10"],
	#	[10, 25, "[]", None, 10, "black", "rgradient:ref", "arial|8|white|10"],
	#	[30, 45, "[]", None, 10, "black", "rgradient:orange", "arial|8|white|20"],
	#	[50, 65, "[]", None, 10, "black", "rgradient:pink", "arial|8|white|20"],
	#	[70, 85, "[]", None, 10, "black", "rgradient:green", "arial|8|white|20"],
	#	[90, 105, "[]", None, 10, "black", "rgradient:brown", "arial|8|white|20"],
	#	[110, 125, "[]", None, 10, "black", "rgradient:yellow", "arial|8|white|20"],
	#]

	#seqFace = SeqMotifFace(seq=None, motifs=box_motifs, gap_format="line")
	seqFace = SeqMotifFace(seq=None, motifs=Motifs, gap_format="line")
	(t & node_name).add_face(seqFace, 0, 'aligned')




if REROOT:
	t.set_outgroup( t & reroot_at )
else:
	# Auto-reroot on taxon with highest divergence corrected
	t.set_outgroup( t & greatest_div['element'].lstrip('LTR_retrotransposon') )

ts = TreeStyle()

if LEAF_LABELS:
	ts.show_leaf_name = True
else:
	ts.show_leaf_name = False

ts.show_branch_support = True

#ts.mode = "c"
#ts.arc_start = -180 # 0 degrees = 3 o'clock
#ts.arc_span = 180

if ULTRAMETRIC:
	t.convert_to_ultrametric()
#t.render("{0}_phylo_uncorrectedDivergences.png".format(treeName), w=10, units="in", tree_style=ts)
t.render("{0}_phylo_uncorrectedDivergences.png".format(treeName), w=35, units="in", tree_style=ts)


#
#
#
#
#
#
#
#

t = Tree(tree_flpath)

greatest_divc = {'element':None, 'value':0}
for node in t:
	node_name = str(node).split('-')[-1]
	rt_name = 'LTR_retrotransposon{0}'.format(node_name.split('_')[0])
	if rt_name in divergences:
		if greatest_divc['value'] < float(divergencesCorrected[rt_name]):
			greatest_divc['value'] = float(divergencesCorrected[rt_name])
			greatest_divc['element'] = rt_name
		col = color_gradient['hex'][int(divergencesCorrected[rt_name][:6].replace('.', ''))]
	else:
		col = '#417849'
	if ANNOT:
#        descFace.margin_top = 10
#        descFace.margin_bottom = 10
		#blast_hit_face = faces.TextFace(blast_hit_dct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		try:
			classifFace = faces.TextFace(classifDct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		except KeyError:
			classifFace = faces.TextFace('ERROR', fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		#domain_face = faces.TextFace(domains_dct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')
		try:
			IGCface = faces.TextFace(IGCdct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')
		except KeyError:
			IGCface = faces.TextFace('ERROR', fsize = 8, penwidth = 10, fgcolor = 'Black')
		#(t & node_name).add_face(domain_face, 0, 'aligned')
		(t & node_name).add_face(IGCface, 0, 'aligned')
		#(t & node_name).add_face(blast_hit_face, 1, 'branch-right')
		(t & node_name).add_face(classifFace, 1, 'branch-right')

	cf = CircleFace(radius=4, color=col)
	(t & node_name).add_face(cf, 0, 'branch-right')

        #el = 'LTR_retrotransposon{0}'.format(leaf)
	el = 'LTR_retrotransposon{0}'.format(node_name)
	Motifs = []
	for feat in LTRRTs[el]:
		if feat == 'long_terminal_repeat':
			for LTR in LTRRTs[el][feat]:
				# Add box for LTRs
				start = LTRRTs[el][feat][LTR][0]
				end = LTRRTs[el][feat][LTR][1]
				if strandDct[el] == '-': # invert the graphic
					start = int(LTRRTlengths[el]['length']) - start
					end = int(LTRRTlengths[el]['length']) - end
				motif = [ start//50, end//50, "[]", None, 8, "black", "black", None ]# start, end, typ, w, h, fg, bg, name
				#start, end, typ, w, h, fg, bg, name
				#motif = [ LTRRTs[el][feat][LTR][0], LTRRTs[el][feat][LTR][1], "[]", None, 4, "black", "black", None ]
				#motif = [ LTRRTs[el][feat][LTR][0]//100, LTRRTs[el][feat][LTR][1]//100, "[]", None, 4, "black", "black", None ]
				#motif = [ LTRRTs[el][feat][LTR][0]//50, LTRRTs[el][feat][LTR][1]//50, "[]", None, 8, "black", "black", None ]
				Motifs.append(motif)
		else: # domains
			#motif = [ LTRRTs[el][feat][0], LTRRTs[el][feat][1], "[]", None, 1, "blue", "blue", None ]

			
			if feat in categories['zf'] or len([ cat for cat in categories['zf'] if feat.startswith(cat) ]) > 0:
				domainColor = "Pink" # Zinc finger is pink
			elif feat in categories['rvt'] or len([ cat for cat in categories['rvt'] if feat.startswith(cat) ]) > 0:
				domainColor = "Red" # Reverse transcriptase is red
			elif feat in categories['rvp'] or len([ cat for cat in categories['rvp'] if feat.startswith(cat) ]) > 0:
				domainColor = "DeepPink" # Retroviral protein is deeppink
			elif feat in categories['rnaseh'] or len([ cat for cat in categories['rnaseh'] if feat.startswith(cat) ]) > 0:
				domainColor = "DarkOrange" #RNaseH is Orange
			elif feat in categories['recombinase'] or len([ cat for cat in categories['recombinase'] if feat.startswith(cat) ]) > 0:
				domainColor = "Yellow" # Recombinase is yellow
			elif feat in categories['integrase'] or len([ cat for cat in categories['integrase'] if feat.startswith(cat) ]) > 0:
				domainColor = "cyan" # Integrase is Cyan
			elif feat in categories['gag'] or len([ cat for cat in categories['gag'] if feat.startswith(cat) ]) > 0:
				domainColor = "purple" # Gag is purple
			elif feat in categories['duf'] or len([ cat for cat in categories['duf'] if feat.startswith(cat) ]) > 0:
				domainColor = "ForestGreen" 
			elif feat in categories['asp'] or len([ cat for cat in categories['asp'] if feat.startswith(cat) ]) > 0:
				domainColor = "Lime"
			elif feat in categories['various'] or len([ cat for cat in categories['various'] if feat.startswith(cat) ]) > 0:
				domainColor = "silver"
			elif 'ORFHIT' in feat:
				domainColor = "LightSeaGreen"
			elif 'ORF' in feat:
				domainColor =  "BurlyWood"
			else:
				sys.exit('feat not in categories: {0}'.format(feat))

			start = LTRRTs[el][feat][0]
			end = LTRRTs[el][feat][1]
			if strandDct[el] == '-': # invert the graphic
				start = int(LTRRTlengths[el]['length']) - start
				end = int(LTRRTlengths[el]['length']) - end
			if domainColor == "silver":
				motif = [ start//50, end//50, "[]", None, 8, domainColor, domainColor, "arial|6|black|{0}".format(feat) ]
			else:
				#motif = [ LTRRTs[el][feat][0], LTRRTs[el][feat][1], "[]", None, 4, domainColor, domainColor, None ]
				motif = [ start//50, end//50, "[]", None, 8, domainColor, domainColor, None ]
				#motif = [ LTRRTs[el][feat][0]//100, LTRRTs[el][feat][1]//100, "[]", None, 4, domainColor, domainColor, None ]
			Motifs.append(motif)

	#box_motifs = [
	#	# seq.start, seq.end, shape, width, height, fgcolor, bgcolor
	#	[0,  5, "[]", None, 10, "black", "rgradient:blue", "arial|8|white|10"],
	#	[10, 25, "[]", None, 10, "black", "rgradient:ref", "arial|8|white|10"],
	#	[30, 45, "[]", None, 10, "black", "rgradient:orange", "arial|8|white|20"],
	#	[50, 65, "[]", None, 10, "black", "rgradient:pink", "arial|8|white|20"],
	#	[70, 85, "[]", None, 10, "black", "rgradient:green", "arial|8|white|20"],
	#	[90, 105, "[]", None, 10, "black", "rgradient:brown", "arial|8|white|20"],
	#	[110, 125, "[]", None, 10, "black", "rgradient:yellow", "arial|8|white|20"],
	#]

	#seqFace = SeqMotifFace(seq=None, motifs=box_motifs, gap_format="line")
	seqFace = SeqMotifFace(seq=None, motifs=Motifs, gap_format="line")
	(t & node_name).add_face(seqFace, 0, 'aligned')


if REROOT:
	t.set_outgroup( t & reroot_at )
else:
	t.set_outgroup( t & greatest_divc['element'].lstrip('LTR_retrotransposon') )

ts = TreeStyle()

if LEAF_LABELS:
	ts.show_leaf_name = True
else:
	ts.show_leaf_name = False

ts.show_branch_support = True

#ts.mode = "c"
#ts.arc_start = -180 # 0 degrees = 3 o'clock
#ts.arc_span = 180

if ULTRAMETRIC:
	t.convert_to_ultrametric()
#t.show(tree_style=ts)
#t.render("{0}_phylo_correctedDivergences.png".format(treeName), h=40, w=40, units="in", tree_style=ts)
#t.render("{0}_phylo_correctedDivergences.png".format(treeName),w=15, units='in', tree_style=ts)
t.render("{0}_phylo_correctedDivergences.png".format(treeName),w=35, units='in', tree_style=ts)
#"arial|8|white|hello"
