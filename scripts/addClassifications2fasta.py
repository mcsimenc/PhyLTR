#!/usr/bin/env python3

import sys


def help():
    print('''
    Usage:
    ------------
    addClassifications2fasta.py -gff <path>  \\
                                -fasta <path> \\
                                -dfam <path> \\
                                -repbase <path> > output.fasta

    Description: 
    ------------
    This is intended to be a post-processing script for PhyLTR.py output
    which adds classification (superfamily) names to the sequence
    headers in the file given by -fasta, e.g

        Original header:
            >LTR_retrotransposon1000

        New header:
            >LTR_retrotransposon1000_ERV1

    See below for locations of the specific files from PhyLTR.py output.

    Options:
    ------------
    -gff     <path> PhyLTR.output/GFF_output/*LTRdigestClassifiedNoFP.gff

    -fasta   <path> PhyLTR.output/FASTA_output/LTRdigest_LTR_retrotransposons.fasta

    -dfam    <path> PhyLTR/RepeatDatabases/Dfam/Dfam_ERV_LTR.SF

    -repbase <path> PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.SF
        ''', file=sys.stderr)
    sys.exit(0)


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


def read2dict(filePath):
    """Expects the path to a two column tab-delimited as input. Reads 
    the first two columns in the file given by filePath into a 
    dictionary {column1:column2} and returns the dictionary.
    """
    outDict = {}
    with open(filePath) as fl:
        for line in fl:
            key, val = line.strip().split()
            outDict[key] = val
    return outDict


# print help information if not eanough arguments are present
args = sys.argv
if (len(args) < 9
      or '-gff' not in args
      or '-fasta' not in args
      or '-repbase' not in args
      or '-dfam' not in args):
    help()
# read command line arguments
gffFilePath = args[args.index('-gff') + 1]
fastaFilePath = args[args.index('-fasta') + 1]
repbaseFilePath = args[args.index('-repbase') + 1]
dfamFilePath = args[args.index('-dfam') + 1]
# read repbase file in to a dictionary with element name keys and
# superfamily names as values
repbaseDict = read2dict(repbaseFilePath)
# read dfam file in to a dictionary with element name keys and
# superfamily names as values
dfamDict = read2dict(dfamFilePath)
# read the values for the dfam and repbase homologous results from the
# input gff into a dictionary
gffDict = {}
with open(gffFilePath) as gffFile:
    for line in gffFile:
        if not line.startswith('#'):
            gffLine = GFF3_line(line)
            if gffLine.type == 'LTR_retrotransposon':
                dfamHit = gffLine.attributes['dfamClassification']
                repbaseHit = gffLine.attributes['repbaseClassification']
                ID = gffLine.attributes['ID']
                # allow dfam hits to take precedence over repbase hits
                # if different superfamilies are found for each
                if repbaseHit != 'None':
                    superfamily = repbaseDict[repbaseHit]
                if dfamHit != 'None':
                    superfamily = dfamDict[dfamHit]
                gffDict[ID] = superfamily
# read fasta file, modify sequence headers, and output to stdout
with open(fastaFilePath) as fastaFile:
    for line in fastaFile:
        if line.startswith('>'):
            oldHeader = line.strip()[1:]
            superfamily = gffDict[oldHeader]
            newHeader = '>' + oldHeader + '_' + superfamily
            print(newHeader)
        else:
            print(line.strip())
