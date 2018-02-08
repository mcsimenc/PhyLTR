#!/usr/bin/env python3

import sys
from ete3 import Tree, TreeStyle, SeqMotifFace, CircleFace, faces


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

args = sys.argv

if len(args) < 5:

	print('''
		usage:
			ete3_tree.py -t <newick_tree> -d <divergences> [-reroot <node>] [-nlabel] [-lflabel] [-annot]

		-lflabel
			show leaf labels

		-annot
			show leaf annotations (blast hits and domain architectures
		''', file=sys.stderr)

	sys.exit()

if '-reroot' in args:
	REROOT = True
	reroot_at = args[args.index('-reroot') + 1]
else:
	REROOT = False
	
if '-lflabel' in args:
	LEAF_LABELS = True
else:
	LEAF_LABELS = False

if '-annot' in args:
	ANNOT = True
else:
	ANNOT = False

	
tree_flpath = args[args.index('-t') + 1]
divergences_flpath = args[args.index('-d') + 1]
divergences = {}
divergencesCorrected = {}
#domains_dct = {}
IGCdct = {}
classifDct = {}

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
ts = TreeStyle()

if LEAF_LABELS:
	ts.show_leaf_name = True
else:
	ts.show_leaf_name = False

for node in t:
	node_name = str(node).split('-')[-1]
	rt_name = 'LTR_retrotransposon{0}'.format(node_name.split('_')[0])
#	print(int(divergences[rt_name][:6].replace('.', '')))
#	print(len(color_gradient['hex']))
	col = color_gradient['hex'][int(divergences[rt_name][:6].replace('.', ''))]
	if ANNOT:
#        descFace.margin_top = 10
#        descFace.margin_bottom = 10
		#blast_hit_face = faces.TextFace(blast_hit_dct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		classifFace = faces.TextFace(classifDct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		#domain_face = faces.TextFace(domains_dct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')
		IGCface = faces.TextFace(IGCdct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')

	cf = CircleFace(radius=4, color=col)
	(t & node_name).add_face(cf, 0, 'branch-right')
	if ANNOT:
		#(t & node_name).add_face(domain_face, 0, 'aligned')
		(t & node_name).add_face(IGCface, 0, 'aligned')
		#(t & node_name).add_face(blast_hit_face, 1, 'branch-right')
		(t & node_name).add_face(classifFace, 1, 'branch-right')

if REROOT:
	t.set_outgroup( t & reroot_at )

#t.show(tree_style=ts)
t.render("phylo_uncorrectedDivergences.png", w=10, units="in", tree_style=ts)



for node in t:
	node_name = str(node).split('-')[-1]
	rt_name = 'LTR_retrotransposon{0}'.format(node_name.split('_')[0])
	col = color_gradient['hex'][int(divergencesCorrected[rt_name][:6].replace('.', ''))]
	if ANNOT:
#        descFace.margin_top = 10
#        descFace.margin_bottom = 10
		#blast_hit_face = faces.TextFace(blast_hit_dct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		classifFace = faces.TextFace(classifDct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
		#domain_face = faces.TextFace(domains_dct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')
		IGCface = faces.TextFace(IGCdct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')

	cf = CircleFace(radius=4, color=col)
	(t & node_name).add_face(cf, 0, 'branch-right')
	if ANNOT:
		#(t & node_name).add_face(domain_face, 0, 'aligned')
		(t & node_name).add_face(IGCface, 0, 'aligned')
		#(t & node_name).add_face(blast_hit_face, 1, 'branch-right')
		(t & node_name).add_face(classifFace, 1, 'branch-right')

if REROOT:
	t.set_outgroup( t & reroot_at )

#t.show(tree_style=ts)
t.render("phylo_correctedDivergences.png", w=10, units="in", tree_style=ts)
