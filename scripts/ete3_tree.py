#!/usr/bin/env python3

import sys
from ete3 import (Tree, TreeStyle, SeqMotifFace, CircleFace, faces, 
                                                              add_face_to_node)


def help():
    print('''
    Usage:
    ------------
    ete3_tree.py -t <newick_tree> \\
                 -d <divergences> \\
                 -g <gff> \\
                 [options] 

    Description:
    ------------
    Renders phylogenetic trees for long terminal repeat retrotransposon 
    phylogenies with LTR-RT diagrams and other annotations.

    Options:
    ------------
    -lflabel
        Show element IDs (integers) as leaf labels.

    -classif
        Add the superfamily classification for each element above each 
        LTR RT diagram.

    -geneconv
        Add the word 'Yes' or 'No' to the immediately to the right of 
        each leaf depending on whether intra-element gene conversion 
        tracts were detected between the LTRs of that element.

    -reroot <int>|auto
        Roots tree to one of the elements. Two options are possible for 
        -reroot: 'auto', or an <int> corresponding to the element name 
        (i.e. LTR_retrotransposon<int>) to position as the earliest 
        diverging lineage. Only use -reroot auto if the newick filename 
        contains the outgroup number in the name of the file as
        formatted by PhyLTR. 

    -ultrametric
        Draw tree after applying ete3's convert_to_ultrametric() 
        function. In my experience the trees don't look ultrametric.

    -orfhits <path>
        A file containing List of ORF IDs 
        (e.g. LTR_retrotransposon1224.ORF.08) to color teal. For 
        example, a file wth a list of IDs for ORFs that had blast hits 
        in a database.

    -transcribed <path>
        A file containing a list of element IDs 
        (e.g. LTR_retrotransposon123) to mark with a green asterisk.
        Instead of an asterisk, a 'T' is shown if -classif is used.
    ''', file=sys.stderr)
    sys.exit()


class GFF3_line:
    """A class to represet GFF3 lines and allow modification of the
    values of its fields.

    Attributes:
    ------------
    field0, ..., field8   strings containing the values of each field
                          in a GFF3 line
    attributes            a dictionary containing the key-value pairs
                          in the GFF3 line 9th field

    Methods:    
    ------------
    str()               Outputs GFF3 line
    repr()              Outputs GFF3 line
    refreshAttrStr()    This needs to be called if changes were made to
                        any of the attributes. It refreshes
    """

    def __init__(self, line):
        """GFF3_line is initialized to contain the fields of the GFF3
        line provided as attributes. The attributes are kept in a
        dictionary and a the keys are ordered in a list to preserve
        the order of attributes upon getting the sring representation
        of the line from the GFF3_line object.
        """
        (self.seqid, 
         self.source, 
         self.type, 
         self.start, 
         self.end, 
         self.score, 
         self.strand, 
         self.phase, 
         self.attributes_str) = line.strip().split('\t')
        # preserve attribute order as a list of keys (attributes_order)
        attributes_list = self.attributes_str.split(';')
        self.attributes_order = [attr.split('=')[0] for attr in 
                                                               attributes_list]
        # store attribute keys and their values in a dictionary
        self.attributes = {attr.split('=')[0]:attr.split('=')[1] for attr in 
                                                               attributes_list}
        # rename the name attribute key to Name so it conforms to the
        # GFF3 specification, where Name is a reserved attribute key
        if 'name' in self.attributes:
            self.attributes['Name'] = self.attributes.pop('name')
            self.attributes_order[self.attributes_order.index('name')] = 'Name'

    def __repr__(self):
        """Output for overloaded functions str() and repr()"""
        return '\t'.join([str(self.seqid), 
                          str(self.source), 
                          str(self.type), 
                          str(self.start), 
                          str(self.end), 
                          str(self.score), 
                          str(self.strand), 
                          str(self.phase), 
                          str(self.attributes_str)])

    def refreshAttrStr(self):
        """If the attributes dictionary or attributes_order has been 
        altered this should be called to update attributes_str.
        """
        self.attributes_str = ';'.join(['='.join(
             [attr, self.attributes[attr]]) for attr in self.attributes_order])


def hex_to_RGB(hex):
        """"Converts hexadecimal triplet to RGB format
        e.g.: #FFFFFF" -> [255,255,255]
        """
        # Pass 16 to the integer function for change of base
        return [int(hex[i:i + 2], 16) for i in range(1, 6, 2)]


def RGB_to_hex(RGB):
        """Converts RGB to hexadecimal format 
        e.g.: [255,255,255] -> "#FFFFFF"
        """
        # Components need to be integers for hex to make sense
        RGB = [int(x) for x in RGB]
        return '#'+''.join(['0{0:x}'.format(v) if v < 16 
                                          else '{0:x}'.format(v) for v in RGB])


def color_dict(gradient):
        """Takes in a list of RGB sub-lists and returns dictionary of
        colors in RGB and hexadecimal format
        """
        return {'hex':[RGB_to_hex(RGB) for RGB in gradient], 
                  'r':[RGB[0] for RGB in gradient], 
                  'g':[RGB[1] for RGB in gradient], 
                  'b':[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
        """Returns a gradient list of (n) colors between
        two hex colors. start_hex and finish_hex
        should be the full six-digit color string,
        inlcuding the number sign ("#FFFFFF")
        """
        # Starting and ending colors in RGB form
        s = hex_to_RGB(start_hex)
        f = hex_to_RGB(finish_hex)
        # Initilize a list of the output colors with the starting color
        RGB_list = [s]
        # Calcuate a color at each evenly spaced value of t from 1 to n
        for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
                curr_vector = [int(s[j] + (float(t)/(n-1))*(f[j]-s[j])) 
                                                             for j in range(3)]
                # Add it to our list of output colors
                RGB_list.append(curr_vector)
        return color_dict(RGB_list)


def polylinear_gradient(colors, n):
        """Returns a list of colors forming linear gradients between
        all sequential pairs of colors. "n" specifies the total
        number of desired output colors
        """
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
    """Reads a GFF3 file and outputs a dictionary containing a mapping
    of feature ID (default) of type LTR_retrotransposon (default) to 
    strand. Trim removes the first n characters from the attribute value
    where n is the number of characters in the string Type.
    """
    strandDct = {}
    with open(gff, 'r') as inFl:
        for line in inFl:
            # skip commented lines
            if line.startswith('#'):
                continue
            # read only lines of feature type Type
            if '\t{0}\t'.format(Type) in line:
                gffLine = GFF3_line(line)
                Attr = gffLine.attributes[attr]
                # remove prefix from the attribute value
                if trim:
                    Attr = Attr[len(Type):]
                strandDct[Attr] = gffLine.strand
    return strandDct


# output help information if command line arguments are missing
args = sys.argv
if (len(args) < 5
      or '-t' not in args
      or '-g' not in args
      or '-d' not in args):
    help()
# create a map of each type of domain to the domain names present in
# the pHMM. Categories are: aspartyl protease, various-which are drawn
# in grey with the pHMM name written in black, DUFs, gag, integrase,
# recombinase, RNaseH, retroviral protein, reverse transcriptase, and
# zinc finger
categories = {
                'asp':['AP2',
                       'Asp',
                       'Asp_protease_2',
                       'Asp_protease',
                       'dUTPase_2',
                       'dUTPase',
                       'Exc',
                       'Exo_endo_phos_2',
                       'FB_lectin',
                       'Foamy_virus_ENV'],

                'various':['ATHILA',
                           'BCNT',
                           'DBD_Tnp_Mut',
                           'DDE_Tnp_IS1595',
                           'DDE_Tnp_ISL3',
                           'DpnI',
                           'HTH_psq',
                           'HTH_Tnp_1',
                           'HTH_Tnp_Tc3_2',
                           'HTH_Tnp_Tc5',
                           'IN_DBD_C',
                           'Inhibitor_I34',
                           'Gypsy',
                           'Herpes_ORF11',
                           'Intron_maturas2',
                           'LEDGF',
                           'Maelstrom',
                           'Maff2',
                           'Mu-transpos_C',
                           'N-Term_TEN',
                           'Nup153',
                           'Nup_retrotrp_bd',
                           'PEN-2',
                           'Peptidase_A17',
                           'Peptidase_A2B',
                           'Peptidase_A2E',
                           'Phage_Cox',
                           'Phage_GPA',
                           'PHINT_rpt',
                           'Piwi',
                           'RAG1',
                           'RdRP_5',
                           'RE_AlwI',
                           'RE_SacI',
                           'Retro_M',
                           'Slu7',
                           'SNF5',
                           'SPP',
                           'SQAPI',
                           'SSV1_ORF_D-335',
                           'Sulfolobus_pRN',
                           'TagA',
                           'Telomerase_RBD',
                           'Thg1',
                           'TLV_coat',
                           'Tn916-Xis',
                           'Tnp_P_element_C',
                           'Tnp_zf-ribbon_2',
                           'Transposase_22',
                           'Transposase_28',
                           'Transp_Tc5_C',
                           'TYA',
                           'Vif',
                           'WCCH',
                           'Y1_Tnp',
                           'Yuri_gagarin',
                           'Zea_mays_MuDR'],

                'duf':['DUF1725',
                       'DUF1759',
                       'DUF3158',
                       'DUF3173',
                       'DUF3258',
                       'DUF3435',
                       'DUF3701',
                       'DUF3806',
                       'DUF390',
                       'DUF4102',
                       'DUF4219',
                       'DUF4413'],

                'gag':['gag-asp_proteas',
                       'Gag_MA',
                       'Gag_p10',
                       'Gag_p12',
                       'Gag_p17',
                       'Gag_p19',
                       'Gag_p24',
                       'Gag_p30',
                       'Gag_p6',
                       'gag_pre-integrs',
                       'Gag_spuma',
                       'Retrotran_gag_2',
                       'Retrotran_gag_3',
                       'Retrotrans_gag'],

                 'integrase':['Integrase_1',
                              'Integrase_AP2',
                              'Integrase_DNA',
                              'Integrase_Zn',
                              'Phage_integ_N',
                              'Phage_Integr_2',
                              'Phage_integr_3',
                              'Phage_integrase',
                              'Phage_int_SAM_1',
                              'Phage_int_SAM_2',
                              'Phage_int_SAM_3',
                              'Phage_int_SAM_4',
                              'Phage_int_SAM_5',
                              'rve_2',
                              'rve_3',
                              'rve'],

                'recombinase':['Recombinase'],

                'rnaseh':['RHSP',
                          'RNase_H2_suC',
                          'RNase_H2-Ydr279',
                          'RNaseH_C',
                          'RNase_H',
                          'RNaseH_like'],

                'rvp':['RVP_2',
                       'RVP'],

                'rvt':['RVT_1',
                       'RVT_2',
                       'RVT_3',
                       'RVT_connect',
                       'RVT_N',
                       'RVT_thumb'],

                'zf':['zf-C2H2',
                      'zf-CCHC_2',
                      'zf-CCHC_3',
                      'zf-CCHC_4',
                      'zf-CCHC_5',
                      'zf-CCHC_6',
                      'zf-CCHC',
                      'zf-H2C2',
                      'zf-H3C2',
                      'zf-RVT']
            }
# read command line arguments and set booleans which control which parts
# of the script are used
tree_flpath = args[args.index('-t') + 1]
treeName = tree_flpath.split('/')[-1]
divergences_flpath = args[args.index('-d') + 1]
ANYANNOT = False
if '-reroot' in args:
    REROOT = True
    reroot_at = args[args.index('-reroot') + 1]
    if reroot_at == 'auto':
        reroot_at = tree_flpath.split('.')[-2].split('_')[-1]
else:
    REROOT = False
# read this list of ORFs with sequence homology to a sequence in a
# database other than the pHMMs searched by LTRdigest so it can be
# colored seagreen in the rendered diagram
ORFHITS = False
if '-orfhits' in args:
    orfhitsfl = args[args.index('-orfhits')+1]
    orfhits = { line.strip() for line in open(orfhitsfl, 'r') }
    ORFHITS = True
    ANYANNOT = True
TRANSCRIBED = False
if '-transcribed' in args:
    transcribedfl = args[args.index('-transcribed')+1]
    transcribed = { line.strip():'*' for line in open(transcribedfl, 'r') }
    TRANSCRIBED = True
    ANYANNOT = True
if '-lflabel' in args:
    LEAF_LABELS = True
    ANYANNOT = True
else:
    LEAF_LABELS = False
if '-classif' in args:
    ANNOT = True
    ANYANNOT = True
else:
    ANNOT = False
if '-geneconv' in args:
    GCLABEL = True
    ANYANNOT = True
else:
    GCLABEL = False
if '-ultrametric' in args:
    ULTRAMETRIC = True
else:
    ULTRAMETRIC = False
# dictionaries to store information about each element
divergences = {}
divergencesCorrected = {}
IGCdct = {}
classifDct = {}
LTRRTs = {}
LTRRTlengths = {}
gffFlPth = args[args.index('-g')+1]
strandDct = gff2strands(gff=gffFlPth, 
                        attr='ID', 
                        Type='LTR_retrotransposon', 
                        trim=False)
# parse the gff file for information about element architecture
# and domain types, for drawing diagrams
with open(gffFlPth, 'r') as gffFl:
    for line in gffFl:
        # skip commented lines in the gff
        if line.startswith('#'):
            continue
        # read each gff line and save feature type, start, end
        gffLine = GFF3_line(line)
        feat = gffLine.type
        start = int(gffLine.start)
        end = int(gffLine.end)
        # calculate length of entire element
        length = end - start + 1
        # if feature is repeat_region, record the start and end
        # coordinates and the length of the element for drawing diagrams
        if feat == 'repeat_region':
            el = 'LTR_retrotransposon{0}'.format(
                           gffLine.attributes['ID'].split('repeat_region')[-1])
            LTRRTlengths[el] = {'start':start, 'end':end, 'length':length }
            continue
        # record the start and end coordinates of the two LTRs relative
        # to the beginning of the element for drawing diagrams
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
        # record the name and position of the domains relative to the
        # start of the element for drawing diagrams
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
        # record the name and position of the orfs relative to the start
        # of the element for drawing diagrams
        elif feat == 'ORF':
            el = gffLine.attributes['Parent']
            orfid = gffLine.attributes['ID']
            domain = 'ORF'
            # label orfs with sequence homology to a sequence in a
            # database other than the pHMMs searched by LTRdigest so it
            # can be colored seagreen in the rendered diagram
            if ORFHITS:
                if orfid in orfhits:
                    domain = 'ORFHIT'
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
# read divergences file
with open(divergences_flpath) as in_fl:
    for line in in_fl:
        if line.startswith('LTR'):
            # record name, divergence, and scaled divergence (IGC)
            (rt_name,
                classification, 
                I, 
                clust, 
                clustSize, 
                model, 
                div, 
                divc, 
                IGC) = line.strip().split('\t')
            divergences[rt_name] = div
            divergencesCorrected[rt_name] = divc
            IGCdct[rt_name] = IGC
            classifDct[rt_name] = classification
# generate color gradient for representing divergence values as colored
# circles at tips of leaves in the rendered tree
num_colors = 20000
# yellow, red, blue, black
color_gradient = polylinear_gradient(['#FAFF00', 
                                      '#FF1800', 
                                      '#001EFF', 
                                      '#000000'], 
                                      num_colors)
# load the newick tree
t = Tree(tree_flpath)
# for marking which elements did not have information in the divergences
# file for coloring the circles white
NOLTRDIVERGENCES = False
# record the greatest divergence value for automatically setting the
# outgroup as the element with the most divergent LTRs (estimating the
# branch containing the oldest element as the first split in the tree)
greatest_div = {'element':None, 'value':0}
# scale bootstrap values to percentages
for node in t.traverse():
    node.support = node.support * 100
# assign the colors for the node circles based on divergence
for node in t:
    node_name = str(node).split('-')[-1]
    rt_name = 'LTR_retrotransposon{0}'.format(node_name.split('_')[0])
    if rt_name in divergences:
        if greatest_div['value'] < float(divergences[rt_name]):
            greatest_div['value'] = float(divergences[rt_name])
            greatest_div['element'] = rt_name
        col = color_gradient['hex'][int(
                                    divergences[rt_name][:6].replace('.', ''))]
    else:
        NOLTRDIVERGENCES = True
    if ANYANNOT:
        if ANNOT:
            # add superfamily labels to the rendering
            try:
                classifFace = faces.TextFace(classifDct[rt_name], 
                                             fsize = 8, 
                                             fgcolor = 'DarkBlue', 
                                             penwidth = 8)
            except KeyError:
                classifFace = faces.TextFace('?', 
                                            fsize = 8, 
                                            fgcolor = 'DarkBlue', 
                                            penwidth = 8)
            (t & node_name).add_face(classifFace, 0, 'aligned')
        # add a "yes" or "no" text label depending on whether or not 
        # intra-element gene conversion between LTRs was found for this
        # element
        if GCLABEL:
            GCAVAIL = False
            try:
                IGCface = faces.TextFace(IGCdct[rt_name], 
                                         fsize = 8, 
                                         penwidth = 10, 
                                         fgcolor = 'Black')
                GCAVAIL = True
            except KeyError:
                pass
            if GCAVAIL:
                (t & node_name).add_face(IGCface, 1, 'branch-right')
        # add a green asterisk if the element had evidence of
        # transcription
        if TRANSCRIBED:
            try:
                transcribedFace = faces.TextFace(transcribed[rt_name],
                                                 fsize = 18, 
                                                 fgcolor='Green', 
                                                 penwidth = 18)
            except KeyError:
                transcribedFace = faces.TextFace('', 
                                                 fsize = 8, 
                                                 fgcolor = 'Green', 
                                                 penwidth = 12)
            (t & node_name).add_face(transcribedFace, 2, 'branch-right')
    else:
        classifFace = faces.TextFace(' ', 
                                     fsize = 8, 
                                     fgcolor = 'DarkBlue', 
                                     penwidth = 8)
        (t & node_name).add_face(classifFace, 1, 'branch-right')
    # if not ltr divergence was available for this element draw a white
    # circle, otherwise draw a colored circle
    if NOLTRDIVERGENCES:
        cf = CircleFace(radius=4, color='white')
        (t & node_name).add_face(cf, 0, 'branch-right')
    else:
        cf = CircleFace(radius=4, color=col)
        (t & node_name).add_face(cf, 0, 'branch-right')
    # prepare to draw element diagrams which include LTRs, domains, and
    # ORFs. record information about the length and position of each
    el = 'LTR_retrotransposon{0}'.format(node_name)
    Motifs = []
    for feat in LTRRTs[el]:
        # record information about the length and position of LTRs
        if feat == 'long_terminal_repeat':
            for LTR in LTRRTs[el][feat]:
                start = LTRRTs[el][feat][LTR][0]
                end = LTRRTs[el][feat][LTR][1]
                # invert the diagram if the element is on the - strand
                if strandDct[el] == '-':
                    start = int(LTRRTlengths[el]['length']) - start
                    end = int(LTRRTlengths[el]['length']) - end
                # scale the position by 1/50 to make them a similar size
                # to the tree. the values expected for the motif are:
                # start, end, shape, w, h, fg, bg, name
                motif = [start // 50, end // 50, "[]", None, 8, "black", 
                                                                 "black", None]
                Motifs.append(motif)
        # draw domains with predetermined colors
        else:
            # Zinc finger is pink
            if (feat in categories['zf'] 
                or len([ cat for cat in categories['zf'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "Pink"
            # Reverse transcriptase is red
            elif (feat in categories['rvt'] 
                  or len([ cat for cat in categories['rvt'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "Red"
            # Retroviral protein is deeppink
            elif (feat in categories['rvp'] 
                  or len([ cat for cat in categories['rvp'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "DeepPink"
            # RNaseH is orange
            elif (feat in categories['rnaseh'] 
                  or len([ cat for cat in categories['rnaseh'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "DarkOrange"
            # Recombinase is yellow
            elif (feat in categories['recombinase'] 
                  or len([ cat for cat in categories['recombinase'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "Lime"
            # Integrase is cyan
            elif (feat in categories['integrase'] 
                 or len([ cat for cat in categories['integrase'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "cyan"
            # Gag is purple
            elif (feat in categories['gag'] 
                 or len([ cat for cat in categories['gag'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "purple" 
            # DUFs are forestgreen
            elif (feat in categories['duf'] 
                  or len([ cat for cat in categories['duf'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "ForestGreen" 
            # Aspartyl protease is yellow
            elif (feat in categories['asp'] 
                  or len([ cat for cat in categories['asp'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "Yellow"
            # "Various" domains are silver with the pHMM name of the
            # domain written in black
            elif (feat in categories['various'] 
                  or len([ cat for cat in categories['various'] 
                                               if feat.startswith(cat) ])) > 0:
                domainColor = "silver"
            # ORFs with sequence homology to a sequence in a database
            # other than the pHMMs searched by LTRdigest are
            # LightSeaGreen
            elif 'ORFHIT' in feat:
                domainColor = "LightSeaGreen"
            # Other ORFs are Burlywood (tan-ish)
            elif 'ORF' in feat:
                domainColor =  "BurlyWood"
            else:
                sys.exit('feat not in categories: {0}'.format(feat))
            start = LTRRTs[el][feat][0]
            end = LTRRTs[el][feat][1]
            # invert the location of the domain if it is on the - strand
            if strandDct[el] == '-':
                start = int(LTRRTlengths[el]['length']) - start
                end = int(LTRRTlengths[el]['length']) - end
            # if the domain is "Various" add the name from the pHMM of
            # the domain in black writing on the domain. scale length
            # of domain by 1/50 so it is a good size relative to the
            # tree
            if domainColor == "silver":
                motif = [start // 50, end // 50, "[]", None, 8, domainColor, 
                                 domainColor, "arial|6|black|{0}".format(feat)]
            else:
                motif = [start // 50, end // 50, "[]", None, 8, domainColor, 
                                                             domainColor, None]
            Motifs.append(motif)
    # draw domains and LTRs
    seqFace = SeqMotifFace(seq=None, motifs=Motifs, gap_format="line")
    (t & node_name).add_face(seqFace, 0, position='aligned')
# set the outgroup as the specified element given by -reroot, otherwise
# set the outgroup as the element with the greatest LTR divergence
if REROOT:
    t.set_outgroup(t & reroot_at)
else:
    if greatest_div['element'] == None:
        pass
    else:
        t.set_outgroup( t & greatest_div['element'].lstrip(
                                                       'LTR_retrotransposon') )
ts = TreeStyle()
# draw element numbers on the leaves of the tree
if LEAF_LABELS:
    ts.show_leaf_name = True
else:
    ts.show_leaf_name = False
ts.show_branch_support = True
# convert tree to an ultrametric tree using ete3's process. it does not
# appear as an ultrametric tree in the rendering
if ULTRAMETRIC:
    t.convert_to_ultrametric()
# render the output tree and diagrams
t.render("{0}_phylo_uncorrectedDivergences.png".format(treeName), w=35,
                                                     units="in", tree_style=ts)
# repeat the process again using gene conversion-scaled ltr divergences
t = Tree(tree_flpath)
NOLTRDIVERGENCES = False
greatest_divc = {'element':None, 'value':0}
for node in t.traverse():
    node.support = node.support * 100
for node in t:
    node_name = str(node).split('-')[-1]
    rt_name = 'LTR_retrotransposon{0}'.format(node_name.split('_')[0])
    if rt_name in divergences:
        if greatest_divc['value'] < float(divergencesCorrected[rt_name]):
            greatest_divc['value'] = float(divergencesCorrected[rt_name])
            greatest_divc['element'] = rt_name
        col = color_gradient['hex'][int(divergencesCorrected[rt_name][:6].replace('.', ''))]
    else:
        NOLTRDIVERGENCES = True
    if ANYANNOT:
        if ANNOT:
            try:
                classifFace = faces.TextFace(classifDct[rt_name], fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
            except KeyError:
                classifFace = faces.TextFace('?', fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
            (t & node_name).add_face(classifFace, 0, 'aligned')
        if GCLABEL:
            GCAVAIL = False
            try:
                IGCface = faces.TextFace(IGCdct[rt_name], fsize = 8, penwidth = 10, fgcolor = 'Black')
                GCAVAIL = True
            except KeyError:
                pass
            if GCAVAIL:
                (t & node_name).add_face(IGCface, 1, 'branch-right')
        if TRANSCRIBED:
            try:
                transcribedFace = faces.TextFace(transcribed[rt_name], fsize = 18, fgcolor = 'Green', penwidth = 18)
            except KeyError:
                transcribedFace = faces.TextFace('', fsize = 8, fgcolor = 'Green', penwidth = 12)
            (t & node_name).add_face(transcribedFace, 2, 'branch-right')
    else:
        classifFace = faces.TextFace(' ', fsize = 8, fgcolor = 'DarkBlue', penwidth = 8)
        (t & node_name).add_face(classifFace, 1, 'branch-right')
    if NOLTRDIVERGENCES:
        cf = CircleFace(radius=4, color='white')
        (t & node_name).add_face(cf, 0, 'branch-right')
    else:
        cf = CircleFace(radius=4, color=col)
        (t & node_name).add_face(cf, 0, 'branch-right')
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
                motif = [ start//50, end//50, "[]", None, 8, "black", "black", None ]# start, end, shape, w, h, fg, bg, name
                Motifs.append(motif)
        else: # domains
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
                motif = [ start//50, end//50, "[]", None, 8, domainColor, domainColor, None ]
            Motifs.append(motif)
    seqFace = SeqMotifFace(seq=None, motifs=Motifs, gap_format="line")
    (t & node_name).add_face(seqFace, 0, position='aligned')
if REROOT:
    t.set_outgroup( t & reroot_at )
else:
    if greatest_div['element'] == None: # This happens when divergences are not obtained for a given cluster/superfamily (e.g. DIRS)
        pass
    else:
        t.set_outgroup( t & greatest_div['element'].lstrip('LTR_retrotransposon') )
ts = TreeStyle()
if LEAF_LABELS:
    ts.show_leaf_name = True
else:
    ts.show_leaf_name = False

ts.show_branch_support = True
if ULTRAMETRIC:
    t.convert_to_ultrametric()
t.render("{0}_phylo_correctedDivergences.png".format(treeName),w=35, units='in', tree_style=ts)
