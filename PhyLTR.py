#!/usr/bin/env python3
import sys
import re
import os
import time
import os.path
import subprocess
import random
from shutil import copyfile, copytree, rmtree
from datetime import datetime
from math import ceil
from Bio import SeqIO, AlignIO
from multiprocessing import Pool, Manager
from copy import copy, deepcopy
from inspect import currentframe, getframeinfo


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
    addAttr()           Adds new GFF3 attribute
    delAttr()           Delete a GFF3 attribute
    refreshAttrStr()    This needs to be called if changes were made to
                        any of the attributes. It refreshes
    """
    def __init__(self, line=None, **kwargs):
        """GFF3_line is initialized to contain the fields of the GFF3
        line provided as attributes. The attributes are kept in a
        dictionary and a the keys are ordered in a list to preserve
        the order of attributes upon getting the sring representation
        of the line from the GFF3_line object.
        """
        if line == None:
            (self.seqid, self.source, self.type, self.start, self.end, self.score, 
                self.strand, self.phase, self.attributes_str) = [None]*9
            self.attributes_order = []
            self.attributes = {}
            
        else:
            (self.seqid, 
             self.source, 
             self.type, 
             self.start, 
             self.end, 
             self.score, 
             self.strand, 
             self.phase, 
             self.attributes_str) = line.strip().split('\t')
            self.start = int(self.start)
            self.end = int(self.end)
            assert self.start <= self.end
            self.coords = (self.start, self.end)
            self.length = self.end - self.start + 1
            attributes_list = self.attributes_str.split(';')
            self.attributes_order = [attr.split('=')[0] for attr in attributes_list]
            self.attributes = {attr.split('=')[0]:attr.split('=')[1] for attr in attributes_list}

        self.line_number = None

        if 'line_number' in kwargs:    
            self.line_number = kwargs['line_number']

        # rename the name attribute so it conforms to GFF3 specifications, 
        # where Name is a reserved attribute key. The version of LTRDigest 
        # makes attribute key name
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

    def addAttr(self, attr_key, attr_val, replace=False, attr_pos=0):
        """ds attribute, default is at the start of the list.
        Default behavior is to add attr_val to a growing
        comma-sep list if attr_key exists. Set replace=True to
        replace existing attr_val.
        """
        if attr_key in self.attributes:
            if replace:
                delAttr(attr_key)
                self.attributes[attr_key] = attr_val
                self.attributes_order.insert(attr_pos, attr_key)
                self.refreshAttrStr()
            else: # grow list
                self.attributes[attr_key] = '{0},{1}'.format(self.attributes[attr_key], attr_val)
                self.refreshAttrStr()
        else:
            self.attributes[attr_key] = attr_val
            self.attributes_order.insert(attr_pos, attr_key)
            self.refreshAttrStr()

    def delAttr(self, attr_key):
        """Deletes attribute."""
        del self.attributes[attr_key]
        self.attributes_order.pop(self.attributes_order.index(attr_key))
        self.refreshAttrStr()


def rename_fasta_seq_headers(in_flpath, 
                             in_header_pattern, 
                             header_to_name_map, 
                             out_flpath):
    """in_fasta_filepath - the fasta file with headers to be renamed.
    in_header_pattern - a regex that will match the part of the header 
                        that corresponds to keys in header_to_name_map
    header_to_name_map - a dict with keys as part of in_fasta_filepath 
                         seq headers that match in_header_pattern and 
                         values as the abbreviated LTR-RT ID (e.g. 24_1, 
                         which corresponds to LTR_retrotransposon24)

    This function depends on the re and Bio.SeqIO modules
    """
    fasta_file = list(SeqIO.parse(in_flpath, 'fasta'))
    in_header_pattern = re.compile(in_header_pattern)

    for seq in fasta_file:
        seq_id = re.search(in_header_pattern, seq.id).group(1)
        seq.id = header_to_name_map[seq_id]
        seq.description = ''
    SeqIO.write(fasta_file, out_flpath, 'fasta')


def write_ltrs_gff3(data):
    """Writes items from a list data to a file (intended to be a GFF3 file 
    in PhyLTR.py). For use with -ltrs option in PhyLTR.py
    """
    if len(data) < 3:
        print('Less than 2 LTRs...?  ...GFF3 not written. DATA:\n{0}'.format(
                                             '\n'.join(data)), file=sys.stderr)

    elif len(data) > 3:
        print('More than 2 LTRs...?  ...GFF3 not written. DATA:\n{0}'.format(
                                             '\n'.join(data)), file=sys.stderr)

    else:
        # data[-1] should have the output filepath while data[0:2] are 
        # the GFF3 lines
        with open(data[-1], 'w') as out_fl: 
            out_fl.write(''.join(data[0:2]))


def call_process(call_str):
    """Just runs subprocess.call. For simplifying parallel execution of 
    PhyLTR.py using the multiprocessing module.
    """
    subprocess.call(call_str, shell=True)


def count_end_gaps(aln):
    """Input is a two-element list of seqs of identical length 
    (a pairwise alignment).

    Returns the number of bases to subtract from the alignment length. 
    The reason being that since this is meant for alignments between 
    two LTRs, LTR boundaries may be off and part of LTR may have been 
    left in the genome. If LTR boundaries were accurate it's possible 
    there would be no end gaps. The length is used in PhyLTR.py for 
    estimating sequence divergence.

    The number of gaps from the left end in the seq with the most 
    gaps on the left end plus the number from the right end

    Algorithm: walk inward from each seq end one char at a time and 
    count how many gaps. Stop a walk when the char is not the gap char. 
    When all walks are non-gap chars, return the sum of the highest 
    left and highest right count.
    """
    # make sure input seqs are the same length
    if len(aln[0]) != len(aln[1]): 
        raise ValueError('Sequences in alignment have different length') 


    # Keys are all same in these dicts. The number from the key can be 
    # stripped and used as the index with which to access the sequence.
    # keeps track of how many gaps on each end of each sequence
    end_gaps = {'L0':0, 'L1':0, 'R0':0, 'R1':0}
    # keeps track of the current character in each walk
    chars = {'L0':'-', 'L1':'-', 'R0':'-', 'R1':'-'} 
    # keeps track of the current position of each walk
    pos = {'L0':0, 'L1':0, 'R0':-1, 'R1':-1}
    while chars['L0'] == ('-' 
                          or chars['L1'] == '-' 
                          or chars['R0'] == '-' 
                          or chars['R1'] == '-'):
        for k in end_gaps:
            if chars[k] == '-':
                # Strip num from key k to get index of sequence
                seq = int(k[1]) 
                # get next position in current walk
                chars[k] = aln[seq][pos[k]]
                # add one to count if this char is a gap char
                if chars[k] == '-':
                    end_gaps[k] += 1
            # get next index to continue walk from left side
            if 'L' in k:
                pos[k] += 1
            # get next indext to continue walk from right side
            if 'R' in k:
                pos[k] -= 1
    # return gap count
    return ( max([end_gaps['L0'],end_gaps['L1']]) 
             + max([end_gaps['R0'],end_gaps['R1']]) )


def fastas2supermatrix(**kwargs):
    """Needs biopython installed

    kwargs expected:
    input_dir='path to dir containing fasta files to be concatenated'
    input_dir is expected to have only fasta files to be concatenated

    1. Check that lengths of all seqs in each input fasta alignment are 
       the same, report seq lengths and quit if they're different.
    2. For each fasta file input, add each sequence identifier as a key 
       in a dictionary with the value of the key a dictionary with the
       name of the fasta file the key and the value the sequence.
    3. Make a set of unique sequence ids.
    4. Build a new dicitonary with keys as sequence identifiers and
       values are concatenated sequences, with Ns or gaps if no 
       sequence is available for that identifier for a given fasta file 
       input.
    5. Write new fasta of concatenated sequences.
    """
    # get names of input files. dir with input files passed to this 
    # function as kwarg input_dir.
    input_fastas = [fl for fl in os.listdir(kwargs['input_dir']) 
                      if not fl.startswith('.')]
    # will become a nested dict holding sequences from each fasta file 
    # for each sequence header
    seqs_dct = {}
    # will hold the length of each alignment
    seq_lengths = {}

    # process each fasta seq
    for fasta in input_fastas:
        # read in alignment
        fasta_contents = list(SeqIO.parse('{0}/{1}'.format(
                                                     kwargs['input_dir'],
                                                     fasta), 'fasta')) 
        # make sure seqs in the alignment are the same length. If they 
        # aren't then print seq ids and lengths to stderr and quit
        all_seqs_lengths = {seq.id:len(seq) for seq in fasta_contents}
        if len(set(all_seqs_lengths.values())) != 1:
            print('Sequences in file: {0} not all of equal length'.format(
                                                                        fasta), 
                                                               file=sys.stderr)
            print('\n'.join(
                        ['sequence\tlength'] 
                      + ['{0}\t{1}'.format(seq, all_seqs_lengths[seq]) 
                            for seq in all_seqs_lengths]), file=sys.stderr)
            sys.exit()
        seq_lengths[fasta] = len(fasta_contents[0])
        # Add sequences to seqs_dct
        for seq in fasta_contents:
            if seq.id in seqs_dct:
                # check for duplicate seq id and stop program if any 
                # duplicates are found
                if fasta in seqs_dct[seq.id]:
                    print('Duplicate sequence identifier in {0}: {1}'.format(
                                                                fasta, seq.id),
                                                              file=sys.stderr)
                    sys.exit()
                # if no duplicate add sequence to nested dict with key
                # as fasta file name
                seqs_dct[seq.id][fasta] = str(seq.seq)
            else:
                # add sequence to nested dict with key as fasta file 
                # name 
                seqs_dct[seq.id] = { fasta:str(seq.seq) } 

    # print concatenated fasta supermatrix
    if 'output_fl' in kwargs:
        out_fl = open(kwargs['output_fl'], 'a')

    for seq in seqs_dct:
        output_seq = ''
        for fasta in input_fastas:
            if fasta not in seqs_dct[seq]:
                output_seq += 'N' * seq_lengths[fasta]
            else:
                output_seq += seqs_dct[seq][fasta]
        
        if 'output_fl' in kwargs:
            out_fl.write('>{0}\n'.format(seq))
            out_fl.write('{0}\n'.format(output_seq))

        else:
            print('>{0}'.format(seq))
            print(output_seq)

    if 'output_fl' in kwargs:
        out_fl.close()


def bedtoolsid2attr(gff_flpath,
                    attr='ID',
                    strand=False, 
                    lstrip=None):
    """This takes a GFF and creates a map of the current attribute 
    specified by attr, to the expected names that will be output by
    bedtools getfasta (a tool for extracting sequences from a FASTA
    which correspond to features in a GFF3.

    Also as a standalone script:

        bedtoolsid2attr.py -gff <gff> [-attr <str>] [-strand]
    """
    STRAND = strand
    attr_pat = re.compile('{0}=(.+?)(;|$)'.format(attr))
    map_dct = {}
    with open(gff_flpath) as in_fl:
        for line in in_fl:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                seqid = fields[0]
                start = str(int(fields[3]) - 1)
                end = fields[4]
                strand = fields[6]
                attr_value = re.search(attr_pat, fields[8]).group(1)
                
                if not lstrip == None:
                    attr_value = attr_value.lstrip(lstrip)

                if STRAND:
                    map_dct['{0}:{1}-{2}({3})'.format(seqid,
                                                      start, 
                                                      end, 
                                                      strand)
                                                      ] = attr_value

                else:
                    map_dct['{0}:{1}-{2}'.format(seqid, 
                                                 start, 
                                                 end)
                                                 ] = attr_value

    return map_dct


def rename_fasta_seq_headers(in_flpath, 
                             in_header_pattern, 
                             header_to_name_map, 
                             out_flpath):
    """I/O for changing fasta sequence headers. Uses regex matching and 
    dictionary associating current with new name.

    in_fasta_filepath - the fasta file with headers to be renamed.
    in_header_pattern - a regex that will match the part of the header 
                        that corresponds to keys in header_to_name_map
    header_to_name_map - a dict with keys as part of in_fasta_filepath 
                         seq headers that match in_header_pattern and 
                         values as the abbreviated LTR-RT ID (e.g. 24_1,
                         which corresponds to LTR_retrotransposon24)

    This function depends on the re and Bio.SeqIO modules
    """
    fasta_file = list(SeqIO.parse(in_flpath, 'fasta'))
    in_header_pattern = re.compile(in_header_pattern)

    for seq in fasta_file:
        seq_id = re.search(in_header_pattern, seq.id).group(1)
        seq.id = header_to_name_map[seq_id]
        seq.description = ''
    
    SeqIO.write(fasta_file, out_flpath, 'fasta')


def ChangeFastaHeaders(inputFastaPath, inputGFFpath, attribute='ID'):
    """Creates map of bedtools getfasta-style features, reads in 
    inputFasta, writes newFasta, deletes inpuFasta, renames newFasta 
    as inputFasta
    """
    bedtoolsIDmap = bedtoolsid2attr(inputGFFpath, attr=attribute)
    newFasta = '{0}.new.tmp'.format(inputFastaPath)
    header_pattern='(.+?:\d+?-\d+?)(?:$|\D)'
    rename_fasta_seq_headers(inputFastaPath,
                             header_pattern, 
                             bedtoolsIDmap, 
                             newFasta)
    os.replace(newFasta, inputFastaPath)


def ChangeFastaHeadersMultiprocessing(bundle):
    """Creates map of bedtools getfasta-style features, reads in 
    inputFasta, writes newFasta, deletes inpuFasta, renames 
    newFasta as inputFasta
    """
    inputFastaPath = bundle[0]
    inputGFFpath = bundle[1]
    attribute=bundle[2]

    bedtoolsIDmap = bedtoolsid2attr(inputGFFpath, attr=attribute)
    newFasta = '{0}.new.tmp'.format(inputFastaPath)
    header_pattern='(.+?:\d+?-\d+?)(?:$|\D)'
    rename_fasta_seq_headers(inputFastaPath, 
                             header_pattern, 
                             bedtoolsIDmap, 
                             newFasta)
    os.replace(newFasta, inputFastaPath)


def mergeCoords(A,B):
    """
    Takes two tuples and outputs two tuples, which will be identical 
    if the original overlap otherwise will be the originals

        A = (a1, a2), B = (b1, b2) and a1<=b1, a1<=a2, b1<=b2

        case 1: a2<=b1 ---> output originals
        case 2: b1<a2 && b2>a2 ---> output (a1, b2)
        case 3: b2<=a2 ---> output A
    """
    assert (min(A) <= min(B)), (
        'tuples given to mergeCoords in wrong order: A={0}, B={1}').format(A,B)

    if min(B) >= max(A):
        return ((A,B), 0)
    elif min(B) < max(A) and max(B) > max(A):
        output = (min(A),max(B))
        return ((output, output), 1)
    elif max(B) <= max(A):
        return ((A,A), 2)
    else:
        raise Exception(
            'Unexpected result from mergeCoords(A,B) using A={0}, B={1}'.format(
                                                                          A,B))


def append2logfile(directory, logfilename, content):
    """Appends string 'content' to directory/logfilename"""
    with open('{0}/{1}'.format(directory, logfilename), 'a') as logfile:
        logtime = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
        logfile.write('{0}\n'.format(logtime))
        logfile.write('{0}\n\n'.format(content))


def write2summary(text):
    """Write text to a hardcoded summary file at output_dir/summary"""
    scriptpath = os.path.realpath(__file__)
    lineno = getframeinfo(currentframe()).lineno + 1
    append2logfile(paths['output_top_dir'], 
                   mainlogfile, 
                   'line {1} in {2}\nWriting to summary file at {0}'.format(
                                        '{0}/summary'.format(
                                                      paths['output_top_dir']),
                                        lineno,
                                        scriptpath))
    with open('{0}/summary'.format(paths['output_top_dir']), 'a') as summaryFl:
        summaryFl.write('{0}\n'.format(text))


def makecall(call, stdout=None, stderr=None, stdin=None):
    """Handles running subprocess.call. Used when making calls without 
    multiprocessing
    """
    if stdout == None and stderr == None and stdin == None:
        subprocess.call(call)
    elif stdout != None:
        with open(stdout, 'w') as outfl:
            if stderr != None:
                with open(stderr, 'w') as errfl:
                    if stdin == None:
                        subprocess.call(call, stdout=outfl, stderr=errfl)
                    else:
                        with open(stdin, 'r') as inFl:
                            subprocess.call(call, 
                                            stdin=inFl,
                                            stdout=outfl,
                                            stderr=errfl)
            elif stderr == None:
                if stdin == None:
                    subprocess.call(call, stdout=outfl)
                else:
                    with open(stdin, 'r') as inFl:
                        subprocess.call(call, stdin=inFl, stdout=outfl)

    elif stderr != None and stdout == None:
        with open(stderr, 'w') as errfl:
            if stdin == None:
                subprocess.call(call, stderr=errfl)
            else:
                with open(stdin, 'r') as inFl:
                    subprocess.call(call, stdin=inFl,  stderr=errfl)
        
    elif stdin != None and stderr == None and stdout == None:
        with open(stdin, 'r') as inFl:
            subprocess.call(call, stdin=inFl)


def makecallMultiprocessing(callBundle):
    """Handles running subprocess.call when using multiprocessing.
    Used to iterate over calls with Pool().map
    """
    call = callBundle[0]
    stdout = callBundle[1]
    stderr = callBundle[2]
    stdin = callBundle[3]
    if stdout == None and stderr == None and stdin == None:
        subprocess.call(call)
    elif stdout != None:
        with open(stdout, 'w') as outfl:
            if stderr != None:
                with open(stderr, 'w') as errfl:
                    if stdin == None:
                        subprocess.call(call, stdout=outfl, stderr=errfl)
                    else:
                        with open(stdin, 'r') as inFl:
                            subprocess.call(call, 
                                            stdin=inFl, 
                                            stdout=outfl, 
                                            stderr=errfl)
            elif stderr == None:
                if stdin == None:
                    subprocess.call(call, stdout=outfl)
                else:
                    with open(stdin, 'r') as inFl:
                        subprocess.call(call, stdin=inFl, stdout=outfl)
    elif stderr != None and stdout == None:
        with open(stderr, 'w') as errfl:
            if stdin == None:
                subprocess.call(call, stderr=errfl)
            else:
                with open(stdin, 'r') as inFl:
                    subprocess.call(call, stdin=inFl,  stderr=errfl)
    elif stdin != None and stderr == None and stdout == None:
        with open(stdin, 'r') as inFl:
            subprocess.call(call, stdin=inFl)


def MakeDir(pathsname, path):
    """Makes a directory path and stores it in the global dict paths 
    under the key pathsname and writes to logfile.
    """
    global paths

    paths[pathsname] = path
    if not os.path.exists(path): 
        os.makedirs(path)
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 1
        append2logfile(paths['output_top_dir'], 
                       mainlogfile, (
                'line {2} in {3}\nCreated dir for {0}:\n{1}').format(pathsname, 
                                                                    path, 
                                                                    lineno, 
                                                                    scriptpath))


def addStrandToGFF(strandDct, GFFpth):
    """Updates strand field for element with ? as strand based on Dfam 
    and Repbase results. Provide a dictionary with LTR RT # (e.g. 4 for 
    LTR_retrotransposon4) as keys and strand as values. A new GFF will 
    be written, and the old one removed.
    """
    newgff = '{0}.updatingstrandinprocess'.format(GFFpth)
    if os.path.exists(newgff):
        os.remove(newgff)

    with open(GFFpth) as inFl:
        for line in inFl:
            if not line.startswith('#'):
                try:
                    gffLine = GFF3_line(line)
                except:
                    print(line)
                    sys.exit()
                if gffLine.strand == '?':
                    if 'Parent' in gffLine.attributes:
                        elNum = re.search('\d*$', 
                                          gffLine.attributes['Parent']
                                          ).group(0)
                    elif 'ID' in gffLine.attributes:
                        elNum = re.search('\d*$', 
                                          gffLine.attributes['ID']
                                          ).group(0)
                    else:
                        print(('WARNING:\taddStrandToGFF() found the following '
                                'GFF file to contain the following line where '
                                'strand is unknown and the attributes lack ID or '
                                'Parent keys\n{0}\n{1}').format(GFFpth, 
                                                               line.strip()), 
                                file=sys.stderr)

                    if elNum in strandDct:
                        if elNum == '?':
                            continue
                        gffLine.strand = strandDct[elNum]
                        with open(newgff, 'a') as outGFFfl:
                            outGFFfl.write(str(gffLine)+'\n')
                    else:
                        with open(newgff, 'a') as outGFFfl:
                            outGFFfl.write(line)
                else:
                    with open(newgff, 'a') as outGFFfl:
                        outGFFfl.write(line)
            else:
                with open(newgff, 'a') as outGFFfl:
                    outGFFfl.write(line)
    os.rename(newgff, GFFpth)


def RemoveNonLTRretrotransposons(LTRdigestGFFfl, 
                                 annotAttr2DbDict, 
                                 outputFlName, 
                                 REPORTCONFLICTS=True, 
                                 KEEPCONFLICTS=False, 
                                 KEEPNOCLASSIFICATION=False, 
                                 logFilePth='conflictingAnnotations.log'):
    """Removes non-LTR retrotransposons from LTRdigest GFF3 input based 
    on classifications from e.g. Repbase or Dfam

    Parses LTRdigestGFFfl and checks attributes specified in 
    annotAttr2DbDict for annotations and looks for them in the Db 
    specified in annotAttr2DbDict. If found, that element is output. 
    If conflicting annotations are found, REPORTCONFLICTS=True writes 
    conflicts to logfile at logFilePth; if KEEPCONFLICTS=True,
    the element is retained in the output GFF.

    Structure of annotAttr2DbDict should be:

        { attr1:DbFlName1, ... }

    where attr are keys in the attributes field of the input GFF3 that 
    contain classification information derived from the database of 
    which DbFlName should be the subset of classifications in the 
    database that are LTR retrotransposons.

    If an attribute's value is "None", it is treated as "no 
    information". If all attributes specified in annotAttr2DbDict are 
    "None" and KEEPNOCLASSIFICATION=True, that element is output.
    """
    LTR_retrotransposon_GFF_lines = {}
    FOUNDNONLTR = False
    Db_LTR_retrotransposon_features = {}
    Db_files = list(annotAttr2DbDict.values())

    for dbFlPth in Db_files:
        with open(dbFlPth) as dbFl:
            for line in dbFl:
                flName = dbFlPth.split('/')[-1]
                if flName in Db_LTR_retrotransposon_features:
                    Db_LTR_retrotransposon_features[flName].add(line.strip())
                else:
                    Db_LTR_retrotransposon_features[flName] = set(
                                                                [line.strip()])

    logfile = open(logFilePth, 'a')
    logfile.write(('#Conflicting annotations (e.g. Dfam->Copia Repbase->Gypsy) '
                   'are reported. If an annotation is Unknown it does not count '
                   'as a conflict\n'))
    with open(LTRdigestGFFfl) as gffFl:
        for line in gffFl:
            if line.startswith('#'):
                continue
            gffLine = GFF3_line(line)
            # First line in LTR retrotransposon block in LTRdigest GFF3
            if gffLine.type == 'repeat_region': 
                FOUNDNONLTR = False
                if gffLine.type in LTR_retrotransposon_GFF_lines:
                    sys.exit('Line\n{0}\nout of order'.format(line.strip()))
                # add repeat_region line to output ('true positive') set
                LTR_retrotransposon_GFF_lines[
                                   gffLine.attributes['ID']] = [ line.strip() ] 
                continue
            elif gffLine.type == 'LTR_retrotransposon':
                ## Check if LTR
                NoClassification = []
                LTRmatching = {}
                # Check database true positives list# Check database 
                # true positives list
                for attr in annotAttr2DbDict: 
                    db = annotAttr2DbDict[attr].split('/')[-1]
                    annot = gffLine.attributes[attr]
                    if not annot == 'None':
                        if (annot in Db_LTR_retrotransposon_features[db] 
                        or annot == 'Unknown'):
                            # is homologous to LTR retrotransposon in 
                            # given db
                            LTRmatching[attr] = True
                        else:
                            # homologous to a non-LTR retrotransposon 
                            # in given db. Flag for removal as 'false 
                            # positive'
                            LTRmatching[attr] = False
                    else:
                        # not homologous to a known LTR RT in given db.
                        NoClassification.append('None')

                LTRmatching_values = list(LTRmatching.values())
                if True in LTRmatching_values and False in LTRmatching_values:

                    # Found conflicting annotations: LTR and Non-LTR
                    if REPORTCONFLICTS:
                        logfile.write(
        '{0}\thas conflicting annotations LTR and Non-LTR\t{1}\n'.format(
            gffLine.attributes['Parent'],
            '\t'.join(['{0}={1}'.format(
                        attr, gffLine.attributes[attr]) 
                        for attr in LTRmatching]
                      )))
                    if KEEPCONFLICTS:
                        LTR_retrotransposon_GFF_lines[
                             gffLine.attributes['Parent']].append(line.strip())
                    else:
                        FOUNDNONLTR = True
                        del(LTR_retrotransposon_GFF_lines[
                                                 gffLine.attributes['Parent']])
                        continue

                # no demonstrated homology to LTR RT in any db
                elif 'None' in NoClassification and LTRmatching_values == []:
                        if KEEPNOCLASSIFICATION:
                            LTR_retrotransposon_GFF_lines[
                             gffLine.attributes['Parent']].append(line.strip())
                            continue

                        del(LTR_retrotransposon_GFF_lines[
                                                 gffLine.attributes['Parent']])
                        FOUNDNONLTR = True
                        continue

                # no demonstratd homology to LTR RT in any db
                elif NoClassification == [] and LTRmatching_values == []:
                    logfile.write('{0}\thas no classification\n'.format(
                                    gffLine.attributes['Parent']))
                    FOUNDNONLTR = True
                    del(LTR_retrotransposon_GFF_lines[
                                                 gffLine.attributes['Parent']])
                    continue

                # homology to non-LTR retrotransposon in a db
                elif False in LTRmatching_values:

                    del(LTR_retrotransposon_GFF_lines[
                                                 gffLine.attributes['Parent']])
                    FOUNDNONLTR = True
                    continue

                else:
                    LTR_retrotransposon_GFF_lines[
                             gffLine.attributes['Parent']].append(line.strip())
                    continue
            else:
                if FOUNDNONLTR:
                    continue

                parent = gffLine.attributes['Parent']

                if 'LTR_retrotransposon' in parent:
                    parent = 'repeat_region{0}'.format(parent.lstrip(
                                                        'LTR_retrotransposon'))
                if not parent in LTR_retrotransposon_GFF_lines:
                    sys.exit(('Line\n{0}\nParent attribute not in LTR '
                              'retrotransposon dictionary').format(line.strip()))

                LTR_retrotransposon_GFF_lines[parent].append(line.strip())

    with open(outputFlName, 'w') as outFl:
        for element in sorted(list(LTR_retrotransposon_GFF_lines.keys())):
            if LTR_retrotransposon_GFF_lines[element] != []:
                for line in LTR_retrotransposon_GFF_lines[element]:
                    outFl.write(line.strip()+'\n')
                outFl.write('###\n')

    logfile.close()


def writeLTRretrotransposonInternalRegions(inputGFFpth, 
                                           outputGFFpth, 
                                           elementSet=None, 
                                           truncateParent=False):
    """Requires Class GFF3_line
    Writes GFF3 for region between two LTRs from a LTRharvest-type file
    Only for elements in elementSet if provided, if elementSet == None,
    all elements are written. if truncateParent=True, Parent attribute 
    has 'LTR_retrotranspson' trimmed from it
    """
    with open(inputGFFpth, 'r') as inGFF:
        currentNewElement = GFF3_line()
        for line in inGFF:
            if '\tlong_terminal_repeat\t' in line:
                gffLine = GFF3_line(line)
                if elementSet == None or (elementSet != None 
                                                    and gffLine.attributes[
                                                                      'Parent'] 
                                                        in elementSet):
                        if not currentNewElement.start == None:
                            if truncateParent == True:
                                gffLine.attributes[
                                              'Parent'] = gffLine.attributes[
                                                        'Parent'].lstrip(
                                                         'LTR_retrotransposon')
                            currentNewElement.end = gffLine.start

                            if (currentNewElement.end 
                                - currentNewElement.start 
                                + 1) <= 0:
                                currentNewElement = GFF3_line()
                                continue

                            assert currentNewElement.attributes['Parent'] == gffLine.attributes['Parent'], (
                                    'GFF long_terminal_repeats may be out of '
                                    'order. Check near {0} or {1} in {2}').format(
                                        currentNewElement.attributes['Parent'],
                                        gffLine.attributes['Parent'],
                                        inputGFFpth)

                            currentNewElement.seqid = gffLine.seqid
                            currentNewElement.source = 'PhyLTR'
                            currentNewElement.type = \
                                           'LTR_retrotransposon_InternalRegion'
                            currentNewElement.score = '.'
                            currentNewElement.strand = gffLine.strand
                            currentNewElement.phase = '.'

                            with open(outputGFFpth, 'a') as outGFFfl:
                                outGFFfl.write('{0}\n'.format(
                                                       str(currentNewElement)))
                                currentNewElement = GFF3_line()
                        else:
                            currentNewElement.start = gffLine.end
                            if truncateParent == True:
                                currentNewElement.attributes[
                                            'Parent'] = gffLine.attributes[
                                                        'Parent'].lstrip(
                                                         'LTR_retrotransposon')
                            else:
                                currentNewElement.attributes[
                                       'Parent'] = gffLine.attributes['Parent']
                            currentNewElement.attributes_order.append('Parent')
                            currentNewElement.refreshAttrStr()
                                

def writeLTRretrotransposonGFF(inputGFFpth, 
                               outputGFFpth, 
                               elementSet=None, 
                               REPEATREGION=False, 
                               truncateParent=True):
    """Requires Class GFF3_line
    Writes GFF3 for LTR_retrotransposon type features from a 
    LTRharvest-type file. Only for elements in elementSet if provided,
    if elementSet == None, all elements are written. If 
    truncateParent=True, Parent attribute has 'LTR_retrotranspson' 
    trimmed from it. If REPEATREGION==True, extract entire repeat region
    """
    global paths

    if os.path.isfile(outputGFFpth):
        os.remove(outputGFFpth)

    scriptpath = os.path.realpath(__file__)
    lineno = getframeinfo(currentframe()).lineno + 1
    if not elementSet == None:
        append2logfile(paths['output_top_dir'], 
                       mainlogfile, 
                       (('line {3} in {4}\n '
                        'Writing LTR_retrotransposon features:\n '
                        '{0}\nfrom:\n{1}\nto:\n{2}')).format(
                                            ','.join(sorted(list(elementSet))),
                                            inputGFFpth, 
                                            outputGFFpth, 
                                            lineno, 
                                            scriptpath))
    else:
        append2logfile(paths['output_top_dir'], 
                       mainlogfile, 
                       ('line {2} in {3}\nWriting LTR_retrotransposon features:\n '
                       'all elements\nfrom:\n{0}\nto:\n{1}').format(inputGFFpth, 
                                                                  outputGFFpth, 
                                                                  lineno, 
                                                                  scriptpath))
    if REPEATREGION:
        feat = '\trepeat_region\t'
    else:
        feat = '\tLTR_retrotransposon\t'
    with open(inputGFFpth, 'r') as inGFF:
        currentNewElement = GFF3_line()
        for line in inGFF:
            if feat in line:
                gffLine = GFF3_line(line)
                if elementSet == None or (elementSet != None 
                                                    and gffLine.attributes[
                                                          'ID'] in elementSet):
                    with open(outputGFFpth, 'a') as outFl:
                        if truncateParent:
                            gffLine.attributes['ID'] = gffLine.attributes[
                                                                     'ID'][19:]
                            gffLine.refreshAttrStr()
                        outFl.write(str(gffLine)+'\n')


def writeLTRsGFF(inputGFFpth, outputGFFpth, elementSet=None):
    """Requires Class GFF3_line
    Writes one GFF3 for each pair of LTRs from a LTRharvest-type file.
    Only for elements in elementSet if provided, if elementSet == None,
    all elements are written. If truncateParent=True, Parent attribute 
    has 'LTR_retrotranspson' trimmed from it
    """
    global paths

    if os.path.isfile(outputGFFpth):
        os.remove(outputGFFpth)

    scriptpath = os.path.realpath(__file__)
    lineno = getframeinfo(currentframe()).lineno + 1
    append2logfile(paths['output_top_dir'], 
                   mainlogfile, 
                   ('line {3} in {4}\nWriting long_terminal_repeat features:\n '
                    '{0}\nfrom:\n{1}\nto:\n{2}').format(
                                            ','.join(sorted(list(elementSet))),
                                            inputGFFpth,
                                            outputGFFpth, 
                                            lineno, 
                                            scriptpath))
    with open(inputGFFpth, 'r') as inGFF:
        currentNewElement = GFF3_line()
        LTR_counts = {}
        for line in inGFF:
            if '\tlong_terminal_repeat\t' in line:
                gffLine = GFF3_line(line)
                if elementSet == None or (elementSet != None
                                              and gffLine.attributes['Parent'] 
                                                in elementSet):
                    with open(outputGFFpth, 'a') as outFl:
                        if gffLine.attributes['Parent'] in LTR_counts:
                            if LTR_counts[gffLine.attributes['Parent']] > 2:
                                sys.exit(('writeLTRsGFF found element {0} to '
                                      'contain more that 2 LTRs in\n{1}').format(
                                                  gffLine.attributes['Parent'], 
                                                  inputGFFpth))
                            gffLine.attributes['ID'] = gffLine.attributes[
                                                               'Parent'] + '_R'
                            LTR_counts[gffLine.attributes['Parent']] += 1
                        else:
                            gffLine.attributes['ID'] = gffLine.attributes[
                                                               'Parent'] + '_L'
                            LTR_counts[gffLine.attributes['Parent']] = 1
                        gffLine.attributes_order = ['ID', 'Parent']
                        gffLine.refreshAttrStr()
                        outFl.write(str(gffLine)+'\n')
                            

def ltrharvest():
    """Runs LTRharvest. LTRharvest options can be specified on the 
    command line. See phyltr -h for defaults and phyltr -help for explanation.
    """
    global paths
    global filenames

    if LTRHARVEST:
        # If this is in paths this step has been completed. Skip
        if not 'inputFastaSuffixArray' in paths:
            MakeDir('suffixerator_dir', '{0}/suffixerator'.format(
                                                      paths['output_top_dir']))
            paths['inputFastaSuffixArray'] = '{0}/{1}.index'.format(
                                                     paths['suffixerator_dir'], 
                                                     paths['inputFasta'])
            gt_suffixerator_call_string = ('gt suffixerator -db {1} '
                                                           '-indexname {0} '
                                                           '-dna -tis -suf -lcp '
                                                           '-des -ssp '
                                                           '1>suffixerator.stdout '
                                            '2>suffixerator.stderr').format(
                                               paths['inputFastaSuffixArray'], 
                                               paths['inputFasta'])
            gt_suffixerator_call = [ executables['genometools'], 
                                     'suffixerator', 
                                     '-db',
                                     paths['inputFasta'], 
                                     '-indexname', 
                                     paths['inputFastaSuffixArray'], 
                                    '-dna', 
                                    '-tis', 
                                    '-suf', 
                                    '-lcp', 
                                    '-des', 
                                    '-ssp' ]
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 1
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('line {2} in {3}\nBegan creating suffix array for  '
                           '{0} using the call:\n{1}').format(
                                                      paths['inputFasta'], 
                                                      gt_suffixerator_call_string, 
                                                      lineno, scriptpath))

            # Run suffixerator
            makecall(gt_suffixerator_call, 
                     '{0}/suffixerator.stdout'.format(
                                                    paths['suffixerator_dir']), 
                     '{0}/suffixerator.stderr'.format(
                                                    paths['suffixerator_dir']))
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 1
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, ('line {0} in {1}\n '
                           'Finished gt suffixerator').format(lineno, 
                                                            scriptpath) )
            paths['suffixeratorInputFastaCopy'] = '{0}/{1}'.format(
                                                     paths['suffixerator_dir'], 
                                                     filenames['inputFasta'])
            copyfile(paths['inputFasta'], paths['suffixeratorInputFastaCopy'])
            # Add suffix array path to status file (for resuming later)
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('inputFastaSuffixArray\t{0}\n'.format(
                                               paths['inputFastaSuffixArray']))
        # If this is in paths this step has been completed, skip
        if not 'LTRharvestGFF' in paths: 
            # Make dir for LTRharvest
            MakeDir('ltrharvest_dir', 
                    '{0}/LTRharvest'.format(paths['output_top_dir']))
            paths['LTRharvestGFF'] = '{0}/{1}.ltrharvest.out.gff'.format(
                                                       paths['ltrharvest_dir'], 
                                                       filenames['inputFasta'])

            # Run LTRharvest
            gt_ltrharvest_call = [executables['genometools'], 
                                  'ltrharvest', 
                                  '-similar', str(ltrharvest_similar), 
                                  '-index', paths['inputFastaSuffixArray'], 
                                  '-gff3', paths['LTRharvestGFF'], 
                                  '-seqids', 'yes', 
                                  '-v', 'yes', 
                                  '-mintsd', str(ltrharvest_mintsd), 
                                  '-maxtsd', str(ltrharvest_maxtsd), 
                                  '-xdrop', str(ltrharvest_xdrop), 
                                  '-mat', str(ltrharvest_mat), 
                                  '-mis', str(ltrharvest_mis), 
                                  '-ins', str(ltrharvest_ins), 
                                  '-del', str(ltrharvest_del), 
                                  '-minlenltr', str(ltrharvest_minlenltr), 
                                  '-maxlenltr', str(ltrharvest_maxlenltr), 
                                  '-mindistltr', str(ltrharvest_mindistltr), 
                                  '-maxdistltr', str(ltrharvest_maxdistltr), 
                                  '-vic', str(ltrharvest_vic)]
            gt_ltrharvest_call_string = '{0} {1}'.format(' '.join(
                                                           gt_ltrharvest_call),  
                                            ('1>ltrharvest.stdout '
                                             '2>ltrharvest.stderr').format(
                                            paths['inputFastaSuffixArray'], 
                                            paths['LTRharvestGFF'],
                                            executables['genometools']))
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 1
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('line {2} in {3}\nBegan running LTRharvest on {0} '
                           'using the call:\n{1}').format(paths['inputFasta'], 
                                                       gt_ltrharvest_call_string, 
                                                       lineno, scriptpath))
            makecall(gt_ltrharvest_call, 
                     '{0}/ltrharvest.stdout'.format(paths['ltrharvest_dir']), 
                     '{0}/ltrharvest.stderr'.format(paths['ltrharvest_dir']))
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 1
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'line {0} in {1}\nFinished gt ltrharvest'.format(
                                                           lineno, scriptpath))

            # Need to sort GFF3, sometimes it's not sorted like 
            # LTRdigest needs it sorted
            gt_sort_call = [executables['genometools'], 'gff3', '-sort', 
                                          '-retainids', paths['LTRharvestGFF']]
            gt_sort_call_string = ('gt gff3 -sort -retainids {0} > {0}.sorted '
                                    '2>{0}.gff3sort.err').format(
                                                        paths['LTRharvestGFF'])
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 1
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('line {1} in {2}\nBegan sorting LTRharvest:\n '
                            '{0}').format(gt_ltrharvest_call_string, 
                                        lineno, 
                                        scriptpath))
            makecall(gt_sort_call,  stdout='{0}.sorted'.format(
                                             paths['LTRharvestGFF']), 
                                             stderr='{0}.gff3sort.err'.format(
                                                       paths['LTRharvestGFF']))
            os.rename('{0}.sorted'.format(paths['LTRharvestGFF']), 
                                          paths['LTRharvestGFF'])
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, ('Below log entry is from line '
                                        '{0} in {1}').format(lineno, 
                                                           scriptpath))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'line {0} in {1}\nFinished sorting GFF3'.format(
                                                                   lineno, 
                                                                   scriptpath))

            # Add LTRharvest GFF3 path to status file (for resuming later)
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('LTRharvestGFF\t{0}\n'.format(
                                                       paths['LTRharvestGFF']))
            paths['CurrentGFF'] = paths['LTRharvestGFF']


def ltrdigest():
    """
    Runs LTRdigest. Domains to search for can be be given as a multi 
    HMM with the command line flag: --ltrdigest_hmms
    """
    global paths
    global filenames

    # Identify parts of element internal regions with evidence of 
    # homology to LTR RT protein coding domains
    if LTRDIGEST:
        os.environ['PATH'] = '{0}:{1}'.format(executables['hmmer'], 
                                              os.environ['PATH'])
        # If this is in paths this step has been completed. Skip it now.
        if not 'LTRdigestGFF' in paths:
            if not 'suffixeratorInputFastaCopy' in paths:
                paths['suffixeratorInputFastaCopy'] = paths['inputFasta']

            MakeDir('ltrdigest_dir', '{0}/LTRdigest'.format(
                                                      paths['output_top_dir']))
            paths['LTRdigestOutputPrefix'] = '{0}/{1}.LTRdigest'.format(
                                                        paths['ltrdigest_dir'], 
                                                        paths['inputFasta'])
            paths['LTRdigestGFF'] = '{0}.gff'.format(
                                                paths['LTRdigestOutputPrefix'])
            filenames['LTRdigestGFF'] = '{0}.LTRdigest.gff'.format(
                                                       filenames['inputFasta'])

            gt_ltrdigest_call = [executables['genometools'], 
                                '-j', str(procs), 
                                'ltrdigest', 
                                '-matchdescstart', 
                                '-outfileprefix', paths['LTRdigestOutputPrefix'], 
                                '-hmms', '{0}'.format(paths['LTRdigestHMMs']), 
                                '-seqfile', paths['suffixeratorInputFastaCopy']]

            gt_ltrdigest_call_string = ('{0} -j {1} ltrdigest -matchdescstart '
                                        '-outfileprefix {2} -hmms {3} -seqfile '
                                        '{4} < {5} > {6}').format(
                                           executables['genometools'], 
                                           procs, 
                                           paths['LTRdigestOutputPrefix'], 
                                           paths['LTRdigestHMMs'], 
                                           paths['suffixeratorInputFastaCopy'], 
                                           paths['CurrentGFF'], 
                                           '{0}.gff'.format(
                                               paths['LTRdigestOutputPrefix']))

            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'Below log entry is from line {0} in {1}'.format(
                                                                   lineno, 
                                                                   scriptpath))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Began running LTRdigest on {0} using the call:\n '
                           '{1}').format(paths['inputFasta'], 
                                       gt_ltrdigest_call_string))
            makecall(gt_ltrdigest_call, '{0}.gff'.format(
                                               paths['LTRdigestOutputPrefix']), 
                                               '{0}/ltrdigest.stderr'.format(
                                                       paths['ltrdigest_dir']), 
                                                       paths['CurrentGFF'])
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'Below log entry is from line {0} in {1}'.format(
                                                                   lineno, 
                                                                   scriptpath))
            append2logfile(paths['output_top_dir'], mainlogfile, 
                                                      'Finished gt ltrdigest' )

            # Add LTRdigest GFF3 path to status file (for resuming later)
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('LTRdigestGFF\t{0}\n'.format(
                                                        paths['LTRdigestGFF']))
            paths['CurrentGFF'] = paths['LTRdigestGFF']

            # Remove suffixerator-generated files
            if not KEEP_UNUSED_FILES:
                rmtree(paths['suffixerator_dir'])
                

def Overlaps(j, k):
    """Inputs, j and k, are tuples/lists with a start and a end 
    coordinate as in: (start, end). If they overlap True is returned, 
    if they do not overlap, False is returned.
    """
    j = sorted([int(i) for i in j])
    k = sorted([int(i) for i in k])
    jk = sorted([j,k], key=lambda x:x[0])
    if jk[0][1] >= jk[1][0]:
        return True
    else:
        return False


def bestORFs(fasta, outdir, gff, minLen=300):
    """Finds all ORFs in fasta using EMBOSS getorf and writes the best 
    ones to a GFF3 and protein FASTA. The best ORFs are the set of 
    non-overlapping ORFs containing the longest ORF out of sets on that
    strand,given by records in gff, or, if strand is unknown, the best 
    set out of all.

    The coordinates output by getorf are 1-based

    Only ORFs with nucleotide sequences longer than minLen are kept. 
    (default 300 bp = 80 aa)
    """
    # the name of the gff file which is the final output of this
    # function
    outgff = '{0}/{1}.orfs.gff'.format(outdir, fasta.split('/')[-1])
    # remove output file if it exists
    if os.path.isfile(outgff):
        os.remove(outgff)
    # read all LTR_retrotransposon features from the gff and record
    # strandedness for all elements
    strands = {}
    with open(gff, 'r') as gffFl:
        for line in gffFl:
            if '\tLTR_retrotransposon\t' in line:
                gffLine = GFF3_line(line)
                el = gffLine.attributes['ID']
                strand = gffLine.strand
                strands[el] = strand
    # run EMBOSS getorf
    outseq = '{0}/{1}.orfs'.format(outdir, fasta.split('/')[-1])
    if os.path.isfile(outseq):
        os.remove(outseq)
    getorf_call = [ executables['getorf'], '-sequence', fasta, 
                                                            '-outseq', outseq ]
    makecall(getorf_call)
    # read orf sequences into list of nonredundant orfs
    orfs = list(SeqIO.parse(outseq, 'fasta'))
    nonredundant_orfs = list()
    unique_orfs = set()
    for seqio in orfs:
        if not seqio.description in unique_orfs:
            unique_orfs.add(seqio.description)
            nonredundant_orfs.append(seqio)
    orfs_ordered_lengths = {}
    orfs_ordered_coords = {}
    orfs_seqs_dct = {}
    orfs_coords = {}
    coords2lenkey = {}
    for orf in nonredundant_orfs:
        element = '_'.join(orf.id.split('_')[:-1])
        strand_gff = strands[element]
        orfnum = orf.id.split('_')[-1]
        desc = orf.description
        seq = str(orf.seq)
        if 'REVERSE' in desc:
            # don't consider reverse strand ORFs if the element is 
            # already assigned to the forward strand.
            if strand_gff == '+':
                continue
            strand = '-'
            # coords from getorf are 1-based, need to -1 from start
            start, end = sorted([int(desc.split()[1][1:]), 
                                 int(desc.split()[3][:-1])]) 
            start -= 1
            length = end - start + 1
            if length < minLen:
                continue
        else:
            # don't consider forward strand ORFs if the element is 
            # already assigned to the reverse strand.
            if strand_gff == '-':
                continue
            strand = '+'
            # coords from getorf are 1-based, need to -1 from start
            start, end = sorted([int(desc.split()[1][1:]), 
                                 int(desc.split()[3][:-1])])
            start -= 1
            length = end - start + 1
            if length < minLen:
                continue

        if element in orfs_coords:
            if strand in orfs_coords[element]:
                orfs_coords[element][strand][orfnum] = (start, end)
            else:
                orfs_coords[element][strand] = {orfnum:(start,end)}
        else:
            orfs_coords[element] = {strand:{orfnum:(start,end)}}

        if element in coords2lenkey:
            if strand in coords2lenkey[element]:
                coords2lenkey[element][strand][(start, end)] = (orfnum, length)
            else:
                coords2lenkey[element][strand] = {(start,end):(orfnum, length)}
        else:
            coords2lenkey[element] = {strand:{(start,end):(orfnum, length)}}

        if element in orfs_ordered_lengths:
            if strand in orfs_ordered_lengths[element]:
                orfs_ordered_lengths[element][strand].append((orfnum, length))
            else:
                orfs_ordered_lengths[element][strand] = [(orfnum, length)]
        else:
            orfs_ordered_lengths[element] = {strand:[(orfnum, length)]}

        if element in orfs_ordered_coords:
            if strand in orfs_ordered_coords[element]:
                orfs_ordered_coords[element][strand].append((start, end))
            else:
                orfs_ordered_coords[element][strand] = [(start, end)]
        else:
            orfs_ordered_coords[element] = {strand:[(start, end)]}

        if element in orfs_seqs_dct:
            if strand in orfs_seqs_dct[element]:
                orfs_seqs_dct[element][strand][orfnum] = seq
            else:
                orfs_seqs_dct[element][strand] = {orfnum:seq}
        else:
            orfs_seqs_dct[element] = {strand:{orfnum:seq}}

    # For each element, find ORFs for its strand or both strands if 
    # strand = ?
    # Order ORFs names in list from element/strand longest to shortest.
    # Order another ORF name list with coordinate position of starts, 
    # smallest to largest
    # Make a dict with the sequences to use after selecting which orfs 
    # to keep
    sorted_elements = sorted(list(orfs_ordered_lengths.keys()))
    for element in sorted_elements:
        strand = strands[element]
        if strand == '?' or strand == '.':
            strand = ['+', '-']
        else:
            strand = [strand]
        best_orf_sets = {'+':None, '-':None}
        for s in strand:
            i = 0
            if s not in orfs_ordered_lengths[element]:
                continue
            orfs_ordered_lengths[element][s].sort(reverse=True, 
                                                  key=lambda x:x[1])
            orfs_ordered_coords[element][s].sort(key=lambda x:x[0])
            # orfs_ordered_length is a list that is modified. i gets 
            # incremented
            while len(orfs_ordered_lengths[element][s]) > i+1:
                orfnum = orfs_ordered_lengths[element][s][i][0]
                # coords of the current orf
                coord = orfs_coords[element][s][orfnum] 
                # current largest orf
                j = orfs_ordered_coords[element][s].index(coord)
                # check for overlaps with next in proximity
                k = j+1
                # Compare j with successively further away orfs until 
                # a non-overlap is reached
                if not k > len(orfs_ordered_coords[element][s])-1:
                    # corresponding occurence in lengths dict for j, 
                    # the current longest ORF
                    J = orfs_ordered_lengths[element][s].index(
                                                    coords2lenkey[element][s][
                                                      orfs_ordered_coords[
                                                               element][s][j]]) 
                    # corresponding occurence in lengths dict for k, the 
                    # current ORF closest to j if moving toward position 0
                    K = orfs_ordered_lengths[element][s].index(coords2lenkey[
                                              element][s][orfs_ordered_coords[
                                                               element][s][k]])
                    while Overlaps( orfs_ordered_coords[element][s][j], 
                                          orfs_ordered_coords[element][s][k] ):
                        # k overlaps j. remove k. because it is shorter than j.
                        coord_removed = orfs_ordered_coords[element][s][k]
                        orfs_ordered_coords[element][s] = [ 
                                   item for item in orfs_ordered_coords[
                                    element][s] if not item == coord_removed ]
                        lenkey = coords2lenkey[element][s][coord_removed]
                        orfs_ordered_lengths[element][s].remove(lenkey)
                        # current largest orf
                        j = orfs_ordered_coords[element][s].index(coord) 
                        # check for overlaps with next in proximity
                        k = j+1
                        if k > len(orfs_ordered_coords[element][s])-1:
                            break
                        # corresponding occurence in lengths dict for j, 
                        # the current longest ORF
                        J = orfs_ordered_lengths[element][s].index(
                                           coords2lenkey[element][s][
                                           orfs_ordered_coords[element][s][j]]) 
                        # corresponding occurence in lengths dict for k, 
                        # the current ORF closest to j if moving toward 
                        # position 0
                        K = orfs_ordered_lengths[element][s].index(
                                           coords2lenkey[element][s][
                                           orfs_ordered_coords[element][s][k]]) 
                # Compare j with successively further away orfs until a
                # non-overlap is reached
                # current largest orf
                j = orfs_ordered_coords[element][s].index(coord) 
                # check for overlaps with previous in proximity
                m = j-1 
                if m > 0:
                    # corresponding occurence in lengths dict for j, 
                    # the current longest ORF
                    J = orfs_ordered_lengths[element][s].index(
                                        coords2lenkey[element][s][
                                        orfs_ordered_coords[element][s][j]])
                    # corresponding occurence in lengths dict for m, 
                    # the current ORF closest to j if moving toward 
                    # position 0
                    M = orfs_ordered_lengths[element][s].index(
                                            coords2lenkey[element][s][
                                            orfs_ordered_coords[element][s][m]])
                    while Overlaps( orfs_ordered_coords[element][s][j], 
                                    orfs_ordered_coords[element][s][m] ):
                        coord_removed = orfs_ordered_coords[element][s][m]
                        orfs_ordered_coords[element][s] = [
                            item for item in orfs_ordered_coords[element][s] if
                                                     not item == coord_removed]
                        lenkey = coords2lenkey[element][s][coord_removed]
                        orfs_ordered_lengths[element][s].remove(lenkey)
                        # current largest orf
                        j = orfs_ordered_coords[element][s].index(coord) 
                        # check for overlaps with previous in proximity
                        m = j-1 
                        if m < 0:
                            break
                        # corresponding occurence in lengths dict for j, 
                        # the current longest ORF
                        J = orfs_ordered_lengths[element][s].index(
                                        coords2lenkey[element][s][
                                        orfs_ordered_coords[element][s][j]])
                        # corresponding occurence in lengths dict for m, 
                        # the current ORF closest to j if moving toward 
                        # position 0
                        M = orfs_ordered_lengths[element][s].index(
                                        coords2lenkey[element][s][
                                        orfs_ordered_coords[element][s][m]])
                i += 1
            orf_ids = [ p[0] for p in orfs_ordered_lengths[element][s] ]
            best_orfs = [[element, 
                          s, 
                          orf_id, 
                          orfs_seqs_dct[element][s][orf_id], 
                          orfs_coords[element][s][orf_id]]  for 
                                                             orf_id in orf_ids]
            best_orf_sets[s] = best_orfs
        if best_orf_sets['+'] == None and best_orf_sets['-'] == None:
            sys.exit('bestORFs() did not populate best_orf_sets')
        if best_orf_sets['+'] == None:
            best_orfs = best_orf_sets['-']
            strand_used = '-'
        elif best_orf_sets['-'] == None:
            best_orfs = best_orf_sets['+']
            strand_used = '+'
        else:
            lengths_plus = sum([i[4][1]-i[4][0]+1 for i in best_orf_sets['+']])
            lengths_minus = sum([i[4][1]-i[4][0]+1 for i in best_orf_sets['-']])
            if lengths_plus > lengths_minus:
                best_orfs = best_orf_sets['+']
                strand_used = '+'
            else:
                best_orfs = best_orf_sets['-']
                strand_used = '-'
        best_orfs.sort(key=lambda x:x[4][0])
        for orf in best_orfs:
            element, s, orf_id, seq, coords = orf
            start, end = coords
            with open(outgff, 'a') as outFl:
                outFl.write(('{0}\tgetorf\tORF\t{1}\t{2}\t.\t{3}\t.\tParent={4};'
                                        'translated_seq={5}\n').format(element, 
                                                                     start, 
                                                                     end, 
                                                                     strand_used, 
                                                                     element, 
                                                                     seq))


def addORFs(maingff, orfgff, newgff):
    """
    Inserts ORFs into existing LTRdigest/LTRharvest GFF. Expects Orfs 
    were obtained from EMBOSS getorf on output from 
    writeLTRretrotransposonInternalRegions()
    Existing features take precedence, and if ORFs overlap existing 
    features, those ORFs are not included in the final output, newgff.
    """
    # Read orf gff store lines in lists in dict with parent as key
    orfs = {}
    if newgff.endswith('.gff'):
        orffasta = '{0}/{1}'.format('/'.join(newgff.split('/')[:-1]), 
                           '{0}.ORFs.fasta'.format(newgff.split('/')[-1][:-4]))
    else:
        orffasta = '{0}/{1}'.format('/'.join(newgff.split('/')[:-1]), 
                                '{0}.ORFs.fasta'.format(newgff.split('/')[-1]))
    append2logfile(paths['output_top_dir'], 
                   mainlogfile, 
                   'Incorporating ORFs in {0} with GFF {1} in to {2}'.format(
                                                      orfgff, maingff, newgff))
    # read gff containing just ORF features output by bestORFs()
    with open(orfgff, 'r') as inFl:
        for line in inFl:
            if not line.startswith('#'):
                gffLine = GFF3_line(line)
                parent = gffLine.attributes['Parent']
                if parent in orfs:
                    orfs[parent].append(gffLine)
                else:
                    orfs[parent] = [gffLine]
    # Read in main gff
    GFFLines = []
    NewGFFLines = []
    with open(maingff, 'r') as inFl:
        for line in inFl:
            if line.startswith('#'):
                continue
            else:
                gffLine = GFF3_line(line)
                GFFLines.append(gffLine)
    with open(orffasta, 'w') as outFl:
        internalparts = []
        el = None
        firstLTRend = None
        for i in range(len(GFFLines)):
            gl = GFFLines[i]
            if gl.type == 'repeat_region':
                internalparts = []
                el = 'LTR_retrotransposon' + gl.attributes['ID'][13:]
                NewGFFLines.append(gl)
            elif gl.type == 'target_site_duplication':
                NewGFFLines.append(gl)
            elif gl.type == 'LTR_retrotransposon':
                NewGFFLines.append(gl)
            elif gl.type == 'long_terminal_repeat':
                # This is the second LTR
                if firstLTRend != None:
                    if el in orfs:
                        orf_ct = 0
                        with open('ERR','a') as errFl:
                            errFl.write('GOING IN\n')
                        for orf in orfs[el]:
                            orf.start = firstLTRend + int(orf.start)
                            orf.end = firstLTRend + int(orf.end)
                            # Change the scaffold name
                            orf.seqid = gl.seqid
                            OVERLAP = False
                            for part in internalparts:
                                if Overlaps([orf.start, orf.end], 
                                                       [part.start, part.end]):
                                    OVERLAP = True
                                    break
                            if not OVERLAP:
                                orf_ct += 1
                                orf.attributes['ID'] = '{0}.ORF.{1:02d}'.format(
                                              orf.attributes['Parent'], orf_ct)
                                orf.attributes_order.insert(0, 'ID')
                                orf.refreshAttrStr()
                                outFl.write('>{0}\n{1}\n'.format(
                                             orf.attributes['ID'], 
                                             orf.attributes['translated_seq']))
                                internalparts.append(orf)
                    internalparts.sort(key=lambda x:int(x.start))
                    NewGFFLines += internalparts
                    NewGFFLines.append(gl)
                    firstLTRend = None
                # This is the first LTR
                elif firstLTRend == None:
                    firstLTRend = int(gl.end)
                    NewGFFLines.append(gl)
            else:
                internalparts.append(gl)
    
    with open(newgff, 'w') as outFl:
        outFl.write('##gff-version 3\n')
        for gl in NewGFFLines:
            if gl.type == 'repeat_region':
                outFl.write('###\n')
            outFl.write('{0}\n'.format(str(gl)))


def AnnotateORFs(minLen):
    """
    Uses bestORFs() and addORFs() to add ORFs of length > minLen
    to the GFF3 if they don't overlap existing features.
    """
    global paths

    if not checkStatusFl('WithORFsGFF'):
        MakeDir('ORFsDir', '{0}/AnnotateORFs'.format(paths['output_top_dir']))
        internalGFF = '{0}/internals.gff'.format(paths['ORFsDir'])
        internalFASTA = '{0}/internals.fasta'.format(paths['ORFsDir'])
        writeLTRretrotransposonInternalRegions(paths['CurrentGFF'], 
                                               internalGFF, 
                                               elementSet=None, 
                                               truncateParent=False)
        getfasta_call = [executables['bedtools'], 'getfasta', '-fi', 
                                      paths['inputFasta'], '-bed', internalGFF]
        makecall(getfasta_call, internalFASTA)
        ChangeFastaHeaders(internalFASTA, internalGFF, attribute='Parent')
        # run EMBOSS getorf
        bestORFs(fasta=internalFASTA, outdir=paths['ORFsDir'], 
                                        gff=paths['CurrentGFF'], minLen=minLen)
        # name of the file output by bestORFs() in the previous call.
        # this file contains gff lines for only the ORF features
        orfgff = '{0}/{1}.orfs.gff'.format(paths['ORFsDir'], 
                                           internalFASTA.split('/')[-1])
        # name of the file to be output by addORFs() in the subsequent
        # call, which will contain the LTRdigest-output gff plus
        # ORF features
        withorfsgff='{0}/FullWithORFs_gt_{1}bp.gff'.format(paths['ORFsDir'], 
                                                           minLen)
        addORFs(maingff=paths['CurrentGFF'], orfgff=orfgff, newgff=withorfsgff)
        paths['WithORFsGFF'] = '{0}/{1}.withORFs_gt_{2}bp.gff'.format(
                                        paths['GFFOutputDir'], 
                                        '.'.join(paths['CurrentGFF'].split(
                                             '/')[-1].split('.')[:-1]), minLen)
        copyfile(withorfsgff, paths['WithORFsGFF'])
        with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
            statusFlAppend.write('WithORFsGFF\t{0}\n'.format(
                                                         paths['WithORFsGFF']))
        paths['CurrentGFF'] = paths['WithORFsGFF']


def classify_by_homology(KEEPCONFLICTS=False, 
                         KEEPNOCLASSIFICATION=False, 
                         repbase_tblastx_evalue=1e-5, 
                         nhmmer_reporting_evalue=5e-2, 
                         nhmmer_inclusion_evalue=1e-2):
    """Parses results from tblastx->Repbase and/or nhmmer->Dfam,
    assigning the superfamily annotation from the highest scoring hit
    in the databases as the classification of the element queried.

    "False positives" are defined as elements in the LTRharvest results
    which do not have significant homology to any element in either the
    Repbase or Dfam databases. A GFF3 file is written which contains
    only and all elements which do have significant homology to a
    sequence in Repbase and/or Dfam and was thus assigned a
    classification.
    """
    global paths
    global filenames

    # Extract LTR_retrotransposon sequences for classification using 
    # homology
    if CLASSIFYDFAM or CLASSIFYREPBASE:
        # If this is in paths this step has been completed. Skip
        if not 'LTRharvest_LTR_retrotransposons_fasta' in paths: 
            paths['LTRharvest_LTR_retrotransposons_GFF'] = (
                        '{0}/LTRharvest_LTR_retrotransposons.gff').format(
                                                         paths['GFFOutputDir'])
            paths['LTRharvest_LTR_retrotransposons_fasta'] = \
                        '{0}/LTRharvest_LTR_retrotransposons.fasta'.format(
                                                       paths['FastaOutputDir'])
            append2logfile(paths['output_top_dir'], 
                   mainlogfile, 
                   'Began extracting LTR_retrotransposons from LTRharvest GFF')
            # Write GFF for just LTR_retrotransposon features
            with open(paths['LTRharvestGFF'], 'r') as harvestFl:
                for line in harvestFl:
                    if not line.startswith('#'):
                        gffLine = GFF3_line(line)
                        if gffLine.type == 'LTR_retrotransposon':
                            with open(
                                paths['LTRharvest_LTR_retrotransposons_GFF'], 
                                   'a') as LTRharvest_LTR_retrotransposons_GFF:
                                LTRharvest_LTR_retrotransposons_GFF.write(
                                                  '{0}\n'.format(str(gffLine)))

            append2logfile(paths['output_top_dir'], 
                mainlogfile, 
                'Finished extracting LTR_retrotransposons from LTRharvest GFF')
            getfasta_ltrretrotransposons_call = [executables['bedtools'], 
                                                'getfasta', 
                                                '-fi', paths['inputFasta'], 
                                                '-s', 
                                                '-bed', '{0}'.format(
                                 paths['LTRharvest_LTR_retrotransposons_GFF'])]
            getfasta_ltrretrotransposons_call_string = \
                'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(
                              paths['inputFasta'], 
                              paths['LTRharvest_LTR_retrotransposons_GFF'], 
                              paths['LTRharvest_LTR_retrotransposons_fasta'], 
                              '{0}/bedtools_getfasta.stderr'.format(
                                                      paths['FastaOutputDir']))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Began extracting LTR_retrotransposon sequences '
                              'from LTRharvest GFF:\n{0}').format(
                                     getfasta_ltrretrotransposons_call_string))
            makecall(getfasta_ltrretrotransposons_call, 
                     paths['LTRharvest_LTR_retrotransposons_fasta'], 
                    '{0}/bedtools_getfasta.stderr'.format(
                                                      paths['FastaOutputDir']))
            append2logfile(paths['output_top_dir'], 
                          mainlogfile, 
                         ('Finished extracting LTR_retrotransposon sequences '
                            'from LTRharvest GFF'))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Changing FASTA headers from bedtools getfasta-style '
                              'to LTR_retrotransposon ID'))
            ChangeFastaHeaders(paths['LTRharvest_LTR_retrotransposons_fasta'], 
                               paths['LTRharvest_LTR_retrotransposons_GFF'])
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Done changing FASTA headers from bedtools '
                                'getfasta-style to LTR_retrotransposon ID'))
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write(
                        'LTRharvest_LTR_retrotransposons_fasta\t{0}\n'.format(
                               paths['LTRharvest_LTR_retrotransposons_fasta']))
    # Find evidence of homology to repeats in Dfam using nhmmer
    if CLASSIFYDFAM: 
        # If this is in paths this step has been completed. Skip
        if not 'DfamTable' in paths: 
            # make Dfam classification output dir
            MakeDir('DfamClassificationDir', '{0}/DfamClassification'.format(
                                                      paths['output_top_dir']))
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('DfamClassificationDir\t{0}\n'.format(
                                               paths['DfamClassificationDir']))
            # run hmmsearch of LTR_retrotransposon features from 
            # LTRdigest or LTRharvest on Dfam
            paths['nhmmer_DfamHits_table'] = (
                                    '{0}/{1}.nhmmer_DfamHits.table').format(
                                                paths['DfamClassificationDir'], 
                                                filenames['inputFasta'])
            nhmmer_dfam_call = ['{0}/nhmmer'.format(executables['hmmer']), 
                                '--tblout', paths['nhmmer_DfamHits_table'], 
                                '--incE', str(nhmmer_inclusion_evalue), 
                                '-E', str(nhmmer_reporting_evalue), 
                                '--cpu', str(procs), 
                                paths['DfamDB'], 
                                paths['LTRharvest_LTR_retrotransposons_fasta']]
            nhmmer_dfam_call_string = '{0} {1}'.format(
                           ' '.join(nhmmer_dfam_call), 
                           '1>/dev/null 2>{0}.nhmmer_DfamHits.stderr'.format(
                               '{0}/{1}'.format(paths['DfamClassificationDir'], 
                                                filenames['inputFasta'])))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'Began nhmmer of Dfam:\n{0}'.format(
                                                      nhmmer_dfam_call_string))
            makecall(nhmmer_dfam_call, 
                    '/dev/null', 
                    '{0}.nhmmer_DfamHits.stderr'.format(
                               '{0}/{1}'.format(paths['DfamClassificationDir'], 
                                               filenames['inputFasta'])))
            append2logfile(paths['output_top_dir'], 
                                        mainlogfile, 'Finished nhmmer of Dfam')
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write('DfamTable\t{0}\n'.format(
                                               paths['nhmmer_DfamHits_table']))
            paths['DfamTable'] = paths['nhmmer_DfamHits_table']

        # add Dfam classifications to GFF
        # If this is in paths this step has been completed. Skip
        if not 'GFFwithDfamClassification' in paths:
            # Extract best hits for each query seq
            paths['DfamResultsTableParsed'] = (
                    '{0}/{1}.LTR_retrotransposon_DfamBestHits.tab').format(
                                                paths['DfamClassificationDir'], 
                                                filenames['inputFasta'])
            dfam_results_parse_call_string = ('{0}/nhmmer_table2columns.py < '
                       '{1} > {2} 2>{3}/nhmmer_table2columns.py.stderr').format(
                                               paths['scriptsDir'], 
                                               paths['DfamTable'], 
                                               paths['DfamResultsTableParsed'], 
                                               paths['DfamClassificationDir'])
            dfam_results_parse_call = ['{0}/nhmmer_table2columns.py'.format(
                                                         paths['scriptsDir'])]
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Began extracting best hits from  nhmmer on Dfam '
                         'results:\n{0}').format(dfam_results_parse_call_string))
            makecall(dfam_results_parse_call, 
                     paths['DfamResultsTableParsed'], 
                     '{0}/nhmmer_table2columns.py.stderr'.format(
                                               paths['DfamClassificationDir']), 
                                               stdin=paths['DfamTable'])
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Finished extracting best hits from  nhmmer on Dfam '
                                                                      'results'))
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('DfamResultsTableParsed\t{0}\n'.format(
                                              paths['DfamResultsTableParsed']))

            # Add best hits to GFF
            paths['GFFwithDfamClassification'] = (
                 '{0}/{1}.LTRdigest.withDfam.gff').format(paths['GFFOutputDir'], 
                                                         filenames['inputFasta'])
            add_dfam_hits_to_ltrdigest_gff_call_string = ('{0}/gffAddAttr.py '
                        '-gff {1} -attr dfamClassification -map {2} -mapKey ID '
                        '-restrictType LTR_retrotransposon -replaceIfNone > {3} '
                         '2>{4}/gffAddAttr.py.DfamHits.stderr').format(
                                             paths['scriptsDir'], 
                                             paths['CurrentGFF'], 
                                             paths['DfamResultsTableParsed'], 
                                             paths['GFFwithDfamClassification'], 
                                             paths['GFFOutputDir'])
            add_dfam_hits_to_ltrdigest_gff_call = ['{0}/gffAddAttr.py'.format(
                                                          paths['scriptsDir']), 
                                       '-gff', paths['CurrentGFF'], 
                                       '-attr', 'dfamClassification', 
                                       '-map', paths['DfamResultsTableParsed'], 
                                       '-mapKey', 'ID', 
                                       '-restrictType', 'LTR_retrotransposon', 
                                       '-replaceIfNone']
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Began adding best hits from nhmmer on Dfam results '
                           'to LTRdigest GFF:\n{0}').format(
                                   add_dfam_hits_to_ltrdigest_gff_call_string))
            makecall(add_dfam_hits_to_ltrdigest_gff_call, 
                     paths['GFFwithDfamClassification'], 
                     '{0}/gffAddAttr.py.DfamHits.stderr'.format(
                                                        paths['GFFOutputDir']))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Finished adding best hits from nhmmer on Dfam '
                                                     'results to LTRdigest GFF'))

            # Add LTRdigest GFF3 with Dfam classifications path to 
            # status file (for resuming later)
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write('GFFwithDfamClassification\t{0}\n'.format(
                                           paths['GFFwithDfamClassification']))
            paths['CurrentGFF'] = paths['GFFwithDfamClassification']
    # Find evidence of homology to repeats in Repbase using tblastx
    if CLASSIFYREPBASE: 
        os.environ['BLASTDB'] = '{0}/RepeatDatabases/Repbase/:{0}/{1}'.format(
                                      paths['selfDir'],paths['FastaOutputDir'])
        # If this is in paths this step has been completed. Skip
        if not 'RepbaseTable' in paths: 
            # make Repbase annotation dir
            MakeDir('RepbaseClassificationDir', 
                    '{0}/RepbaseClassification'.format(paths['output_top_dir']))
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write('RepbaseClassificationDir\t{0}\n'.format(
                                            paths['RepbaseClassificationDir']))

            # run tblastx of LTR_retrotransposon features from LTRdigest 
            # or LTRharvest on Repbase
            paths['tblastx_RepbaseHits_table'] = (
                                '{0}/{1}.tblastx_Repbase.tab').format(
                                             paths['RepbaseClassificationDir'], 
                                             filenames['inputFasta'])
            tblastx_repbase_call = ['{0}/tblastx'.format(executables['blast']), 
                                '-db', 'Repbase_ERV_LTR.fasta', 
                                '-query', 
                                paths['LTRharvest_LTR_retrotransposons_fasta'], 
                                '-evalue', str(repbase_tblastx_evalue), 
                                '-outfmt', ('7 qseqid sseqid pident length '
                                                'mismatch gapopen qstart qend '
                                                'sstart send evalue bitscore '
                                                'sstrand'), 
                                '-num_threads', str(procs), 
                                '-max_hsps', '25']
            tblastx_repbase_call_string = ('{0} -db {1} -query {2} -evalue {5} '
                    '-outfmt "7 qseqid sseqid pident length mismatch gapopen '
                    'qstart qend sstart send evalue bitscore sstrand" '
                    '-num_threads {3} -max_hsps 25 1>{4} 2>{4}.stderr').format(
                                    '{0}/tblastx'.format(executables['blast']), 
                                    'Repbase_ERV_LTR.fasta', 
                                    paths['LTRharvest_LTR_retrotransposons_fasta'],
                                    procs, 
                                    paths['tblastx_RepbaseHits_table'], 
                                    repbase_tblastx_evalue)
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 'Began tblastx of Repbase:\n{0}'.format(
                                                  tblastx_repbase_call_string))
            makecall(tblastx_repbase_call, 
                     paths['tblastx_RepbaseHits_table'], 
                     '{0}.stderr'.format(paths['tblastx_RepbaseHits_table']))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'Finished tblastx of Repbase')
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write('RepbaseTable\t{0}\n'.format(
                                           paths['tblastx_RepbaseHits_table']))
            paths['RepbaseTable'] = paths['tblastx_RepbaseHits_table']

        # add Repbase classifications to GFF
        # If this is in paths this step has been completed. Skip
        if not 'GFFwithRepbaseClassification' in paths:
            # Extract best hits for each query seq
            paths['RepbaseResultsTableParsed'] = (
                    '{0}/{1}.LTR_retrotransposon_RepbaseBestHits.tab').format(
                                             paths['RepbaseClassificationDir'], 
                                             filenames['inputFasta'])
            repbase_results_parse_call_string = ('{0}/best_blast_hit.py < {1} > '
                                    '{2} 2>{3}/best_blast_hit.py.stderr').format(
                                            paths['scriptsDir'], 
                                            paths['RepbaseTable'], 
                                            paths['RepbaseResultsTableParsed'], 
                                            paths['RepbaseClassificationDir'])
            repbase_results_parse_call = ['{0}/best_blast_hit.py'.format(
                                                          paths['scriptsDir'])]
            append2logfile(paths['output_top_dir'], mainlogfile, 
                                    ('Began extracting best hits from  tblastx '
                                        'on Repbase results:\n{0}').format(
                                            repbase_results_parse_call_string))
            makecall(repbase_results_parse_call, 
                     paths['RepbaseResultsTableParsed'], 
                     '{0}/best_blast_hit.py.stderr'.format(
                                            paths['RepbaseClassificationDir']), 
                                            stdin=paths['RepbaseTable'])
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write('RepbaseResultsTableParsed\t{0}\n'.format(
                                           paths['RepbaseResultsTableParsed']))
            append2logfile(paths['output_top_dir'], 
              mainlogfile, 
              'Finished extracting best hits from  tblastx on Repbase results')

            # Add best hits to GFF
            if CLASSIFYDFAM:
                gff_for_repbase_classification = paths[
                                                   'GFFwithDfamClassification']
                paths['GFFwithRepbaseClassification'] = (
                        '{0}/{1}.LTRdigest.withDfam.withRepbase.gff').format(
                                                        paths['GFFOutputDir'], 
                                                        filenames['inputFasta'])
            else:
                gff_for_repbase_classification = paths['LTRdigestGFF']
                paths['GFFwithRepbaseClassification'] = (
                                '{0}/{1}.LTRdigest.withRepbase.gff').format(
                                                        paths['GFFOutputDir'], 
                                                        filenames['inputFasta'])

            add_repbase_hits_to_ltrdigest_gff_call_string = ('{0}/gffAddAttr.py '
                     '-gff {1} -attr repbaseClassification -map {2} -mapKey ID '
                     '-restrictType LTR_retrotransposon -replaceIfNone > {3} '
                     '2>{4}/gffAddAttr.py.RepbaseHits.stderr').format(
                                    paths['scriptsDir'], 
                                    gff_for_repbase_classification, 
                                    paths['RepbaseResultsTableParsed'], 
                                    paths['GFFwithRepbaseClassification'], 
                                    paths['GFFOutputDir'])
            add_repbase_hits_to_ltrdigest_gff_call = ['{0}/gffAddAttr.py'.format(
                                    paths['scriptsDir']), 
                                    '-gff', gff_for_repbase_classification, 
                                    '-attr', 'repbaseClassification', 
                                    '-map', paths['RepbaseResultsTableParsed'], 
                                    '-mapKey', 'ID', 
                                    '-restrictType', 'LTR_retrotransposon', 
                                    '-replaceIfNone']
            append2logfile(paths['output_top_dir'], mainlogfile, 
                            ('Began adding best hits from tblastx on Repbase '
                            'results to LTRdigest GFF:\n{0}').format(
                                add_repbase_hits_to_ltrdigest_gff_call_string))
            makecall(add_repbase_hits_to_ltrdigest_gff_call, 
                        paths['GFFwithRepbaseClassification'], 
                        '{0}/gffAddAttr.py.RepbaseHits.stderr'.format(
                                                        paths['GFFOutputDir']))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Finished adding best hits from tblastx on Repbase '
                                                     'results to LTRdigest GFF'))

            # Add LTRdigest GFF3 with Repbase classifications path to 
            # status file (for resuming later)
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write(
                                'GFFwithRepbaseClassification\t{0}\n'.format(
                                        paths['GFFwithRepbaseClassification']))
            paths['CurrentGFF'] = paths['GFFwithRepbaseClassification']
    # Remove false positives from LTRdigest GFF3
    if CLASSIFYDFAM or CLASSIFYREPBASE:
        # If this is in paths this step has been completed. Skip
        if not 'LTRdigestClassifiedNoFP' in paths:
            if CLASSIFYREPBASE:
                # Will have both Dfam and Repbase classifications if 
                # CLASSIFYDFAM==True also
                gff_classified = paths['GFFwithRepbaseClassification'] 
            else:
                gff_classified = paths['GFFwithDfamClassification']

            paths['LTRdigestClassifiedNoFP'] = (
                                '{0}/{1}.LTRdigestClassifiedNoFP.gff').format(
                                                        paths['GFFOutputDir'], 
                                                        filenames['inputFasta'])
            TruePositiveLTRclassificationsDct = {
                        'dfamClassification':paths['DfamTruePosLTRlist'], 
                        'repbaseClassification':paths['RepbaseTruePosLTRlist']}
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Began removing false positives from LTRdigest GFF '
                                                        'with classifications.'))
            # only elements are preseved if they have a LTR-R 
            # classification, or an unknown classification.
            RemoveNonLTRretrotransposons(gff_classified, 
                      TruePositiveLTRclassificationsDct, 
                      outputFlName=paths['LTRdigestClassifiedNoFP'], 
                      REPORTCONFLICTS=True, 
                      KEEPCONFLICTS=KEEPCONFLICTS, 
                      KEEPNOCLASSIFICATION=KEEPNOCLASSIFICATION, 
                      logFilePth='{0}/RemoveNonLTRretrotransposons.log'.format(
                                                        paths['GFFOutputDir']))
            paths['CurrentGFF'] = paths['LTRdigestClassifiedNoFP']

            if CLASSIFYDFAM:
                strands = {}
                with open(paths['DfamResultsTableParsed'], 'r') as inFl:
                    for line in inFl:
                        el, hit, strand = line.strip().split()
                        el = el.lstrip('LTR_retrotransposon')
                        strands[el] = strand
                addStrandToGFF(strands, paths['CurrentGFF'])

            if CLASSIFYREPBASE:
                strands = {}
                with open(paths['RepbaseResultsTableParsed'], 'r') as inFl:
                    for line in inFl:
                        #el, hit, strand = line.strip().split()
                        el, hit = line.strip().split()
                        el = el.lstrip('LTR_retrotransposon')
                        #strands[el] = strand
                #addStrandToGFF(strands, paths['CurrentGFF'])

            append2logfile(paths['output_top_dir'], mainlogfile, 
             'Update strandedness in GFF3 based on Dfam and/or Repbase results')
            append2logfile(paths['output_top_dir'], mainlogfile, 
                ('Finished removing false positives from LTRdigest GFF with '
                                    'classifications. TP file at:\n{0}').format(
                                             paths['LTRdigestClassifiedNoFP']))
            # Add LTRdigest GFF3 with FP removed path to status file 
            # (for resuming later)
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write('LTRdigestClassifiedNoFP\t{0}\n'.format(
                                             paths['LTRdigestClassifiedNoFP']))

            # Remove large tblastx output and Dfam output. best hits are kept
            if not KEEP_UNUSED_FILES:
                #rmtree(paths['RepbaseTable'])
                os.remove(paths['RepbaseTable'])
                #rmtree(paths['DfamTable'])
                os.remove(paths['DfamTable'])


def shortClassif(ElNames=False):
    """- Opens Dfam and Repbase and associates short, class level names
      for elements with the elements in this annotation which have 
      significant homology to them
    - Opens GFF3 and assigns classification based on attributes. 
      Conflicting attributes (e.g. Dfam=Gypsy, Repbase=Copia) are 
      resolved by using the Dfam classification.
    - Returns a dictionary with assignments as keys and a list of 
      elements as values.
    - ElNames=True means this function will return a dict with LTR RT 
       names as keys.
    - ElNames=False means this function will return a dict with classifs 
      as keys and sets of LTR RT names as values.
    """
    DfamNames = {}
    with open(paths['DfamShortNames']) as DfamFl:
        for line in DfamFl:
            full, short = line.strip().split()
            DfamNames[full] = short
    RepbaseNames = {}
    with open(paths['RepbaseShortNames']) as RepbaseFl:
        for line in RepbaseFl:
            full, short = line.strip().split()
            RepbaseNames[full] = short
    
    ElementNames = {}

    with open(paths['CurrentGFF']) as gffFl:
        for line in gffFl:
            if not line.startswith('#'):
                if '\tLTR_retrotransposon\t' in line:
                    gffLine = GFF3_line(line)
                    el = gffLine.attributes['ID']
                    for attr in gffLine.attributes:
                        annot = gffLine.attributes[attr]

                        # There is an annotation for this Dfam hit in 
                        # the paths['DfamShortNames'] file
                        if 'dfam' in attr:
                            if annot in DfamNames:
                                if el in ElementNames:
                                    # Use Repbase classification if it's 
                                    # not Unknown and Dfam is Unknown
                                    if (ElementNames[el] != 'Unknown' 
                                             and DfamNames[annot]) == 'Unknown': 
                                        continue
                                # A Dfam annotation will take precedence 
                                # over a Repbase annotation
                                ElementNames[el] = DfamNames[annot] 
                            else: 
                                # There is not an annotation for this 
                                # Dfam hit in the paths['DfamShortNames'] 
                                # file
                                if el not in ElementNames:
                                    ElementNames[el] = 'Unknown'

                        elif 'repbase' in attr:
                            if el in ElementNames:
                                # Use Dfam classification if Repbase 
                                # classification is unknown, unless it's 
                                # 'Unknown'
                                if ElementNames[el] != 'Unknown': 
                                    continue
                            # Repbase annotation available
                            if annot in RepbaseNames: 
                                ElementNames[el] = RepbaseNames[annot]
                            # No Repbase annotation available, mark as 
                            # Unknown if there is not already a Dfam 
                            # annotation
                            else: 
                                if el not in ElementNames:
                                    ElementNames[el] = 'Unknown'

    paths['GFFByClassification'] = '{0}/ByClassification'.format(
                                                         paths['GFFOutputDir'])
    if not checkStatusFl('GFFByClassification'):
        if os.path.exists(paths['GFFByClassification']):
            rmtree(paths['GFFByClassification'])
        MakeDir('GFFByClassification', paths['GFFByClassification'])

        with open(paths['CurrentGFF']) as gffFl:
            for line in gffFl:
                if not line.startswith('#'):
                    gffLine = GFF3_line(line)
                    if not 'ID' in gffLine.attributes or gffLine.type == 'ORF':
                        el = gffLine.attributes['Parent']
                    else:
                        el = gffLine.attributes['ID']
                    if 'repeat_region' in el:
                        el = 'LTR_retrotransposon{0}'.format(
                                                  el.split('repeat_region')[1])

                    # Write GFFs for each classification
                    with open('{0}/{1}.LTR_RTs.{2}.gff'.format(
                                              paths['GFFByClassification'], 
                                              filenames['inputFasta'], 
                                              ElementNames[el]), 'a') as outFl:
                        if gffLine.type == 'repeat_region':
                            outFl.write('###\n')
                        outFl.write(line)

        with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
            statusFlAppend.write("{0}\t{1}\n".format('GFFByClassification', 
                                                 paths['GFFByClassification']))
    # return 1
    if ElNames:
        return ElementNames

    Classifications = {}
    for el, clasif in ElementNames.items():
        if clasif in Classifications:
            Classifications[clasif].add(el)
        else:
            Classifications[clasif] = set([el])
    # return 2
    return Classifications


def checkStatusFl(key):
    """The 'status' file contains paths to files and directories to help 
    the program resume at points it has completed.
    """
    with open('{0}/status'.format(paths['output_top_dir']), 'r') as statusFl:
        for line in statusFl:
            if line.strip().split('\t')[0] == key:
                return True
        return False


def graph2groups(G):
    """Returns lists of each element in each connected
    component in G as dictionary. DFS recursive algorithm.
    Adapted from: 
    stackoverflow.com/questions/21078445/find-connected-components-in-a-graph
    """
    def dfs(node1):

        nonlocal visited
        nonlocal groups
        nonlocal group
        nonlocal G

        for node2 in G[node1]:
            if not node2 in visited:
                visited.add(node2)
                groups[group].append(node2)
                dfs(node2)
    visited = set()
    group = 0
    groups = {}
    for node in G:
        if node not in visited:
            visited.add(node)
            group += 1
            groups[group] = [node]
            dfs(node)
    return groups


def WickerFam(pId=80, percAln=80, minLen=80, use_ltrs=True, use_internal=True):
    """Creates a classification of elements using the Wicker et al. 2007 
    protocol:

    '80% sequence similarity or more in at least 80% of the aligned 
    sequence.

    Two elements belong to the same family if they:
        1. share 80% (or more) sequence identity
        2. in at least 80% of their coding or internal domain, or 
           within their terminal repeat regions, or in both.
    
    To avoid misclassification of short and possibly random stretches 
    of homology, we recommend analysing only segments of longer than 
    80 bp. TEs that are smaller than 80 bp require specialized analyses.

    The terminal repeat regions and other non-coding regions are the 
    fastest evolving parts of TEs. Therefore, they offer the most 
    specificity in defining families. Allowing the 80-80-80 rule for 
    DNA sequence identity in either the internal domain or in the 
    terminal regions, or both, also addresses the problem caused by 
    frequent TE truncations. In some cases only terminal repeats and 
    non-coding regions may be present, whereas in other cases only 
    parts of the coding region but no terminal repeats may be 
    available for analysis.
    
    In some cases, it may be necessary to add the subfamily taxon, 
    depending on the population structure of a TE family. Subfamilies 
    can be populations of non-autonomous deletion derivatives or 
    distinct subpopulations in large families that can be clearly 
    segregated. The similarity threshold can differ between 
    subfamilies, depending on the number and homogeneity of elements
    described. However, such distinctions are matters for TE specialists 
    and should not be a burden for annotators (see the wikiPoson 
    web site for further discussion). Importantly, the term should
    not be used for groupings above the family level."


    LTRs and internal regions will compared separately

    I. LTRs:
        1. Extract all LTRs and run blastn all by all (runblast())
        2. Parse blast output for pairs that satisfy the parameters
    
    II. Internals:
        1. Extract all internals and run blastn all by all (runblast())
        2. Parse blast output for pairs that satisfy the parameters
    
    III. Assign
        1. Write family designations (Like mcl output)
    """
    global paths
    
    WICKERLTRS = use_ltrs
    WICKERINTERNAL = use_internal
    MakeDir('WickerFamDir', '{0}/WickerFamDir'.format(paths['output_top_dir']))
    time.sleep(.3)
    MakeDir('WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, percAln, 
                                                                       minLen), 
            '{0}/{1}_pId_{2}_percAln_{3}_minLen'.format(paths['WickerFamDir'], 
                                                         pId, percAln, minLen))
    OutDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, 
                                                                        percAln, 
                                                                        minLen)]
    OutDirKey = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, 
                                                                     percAln, 
                                                                     minLen)

    if not checkStatusFl(OutDirKey):
        with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
            statusFlAppend.write("{0}\t{1}\n".format(
                  'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(pId, 
                                                                       percAln, 
                                                                       minLen), 
                  paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                                     pId, 
                                                                     percAln, 
                                                                     minLen)]))

    MakeDir('WickerFam_Cluster_dir', '{0}/Clusters'.format(OutDir))
    OutDir = paths['WickerFam_Cluster_dir']
    for classif in clusters_by_classif:

        if not checkStatusFl(
                     'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                               pId, percAln, minLen, classif)):


            MakeDir('WickerFam_{0}_dir'.format(classif), '{0}/{1}'.format(
                                                              OutDir, classif))
            cOutDir = paths['WickerFam_{0}_dir'.format(classif)]
            # Dictionary representation of a graph that will hold the 
            # blast results
            G = {} 

            if WICKERINTERNAL:
                # Extract internal regions of selected elements
                paths['InternalelementsGFF'] = '{0}/internals.gff'.format(
                                                                       cOutDir)
                paths['InternalelementsFasta'] = '{0}/internals.fasta'.format(
                                                                       cOutDir)
                writeLTRretrotransposonInternalRegions(paths['CurrentGFF'], 
                                             paths['InternalelementsGFF'], 
                                             elementSet=set(
                                                 clusters_by_classif[classif]), 
                                                          truncateParent=False)
                getfasta_call = [executables['bedtools'], 'getfasta', 
                                '-fi', paths['inputFasta'], 
                                '-s', 
                                '-bed', paths['InternalelementsGFF']]
                getfasta_call_string = ('bedtools getfasta -fi {0} -s -bed {1} '
                                        '> {2} 2> {3}').format(paths['inputFasta'], 
                                                  paths['InternalelementsGFF'], 
                                                  paths['InternalelementsFasta'], 
                        '{0}/bedtools_getfasta_InternalRegions.stderr'.format(
                                                                       OutDir))
                makecall(getfasta_call, 
                         paths['InternalelementsFasta'], 
                         '{0}/bedtools_getfasta_InternalRegions.stderr'.format(
                                                                      cOutDir))
                ChangeFastaHeaders(paths['InternalelementsFasta'], 
                                   paths['InternalelementsGFF'], 
                                   attribute='Parent')

                # Get internals seq lengths
                internal_seq_lengths = {s.id:len(s) for s in list(
                                    SeqIO.parse(paths['InternalelementsFasta'], 
                                    'fasta'))}
                elements = list(internal_seq_lengths.keys())

                # blast all by all internals
                paths['Internals_{0}_selfBlastnOut'.format(classif)] = (
                                  '{0}/Internals_selfBlastn.tab').format(cOutDir)
                runblast(query=paths['InternalelementsFasta'], 
                               subject=paths['InternalelementsFasta'], 
                               out=paths['Internals_{0}_selfBlastnOut'.format(
                                                                     classif)], 
                               evalue='1e-5', 
                               outfmt='7', 
                               percid=pId, 
                               blast='blastn', 
                               procs=procs)

                # Parse blast hits to define families
                InternalsBlastPth = paths['Internals_{0}_selfBlastnOut'.format(
                                                                      classif)]
                with open(InternalsBlastPth, 'r') as blastFl:
                    internal_aln_lens = {}
                    for line in blastFl:
                        if line.startswith('#'):
                            continue

                        (query, subj, percent_id, alignment_len, mismatches, 
                                            gap_opens, q_start, q_end, s_start, 
                                            s_end, E_value, 
                                            bit_score) = line.strip().split('\t')
                        # 1. Ignore self-hits
                        if query == subj:
                            continue
                        # 2. Ignore any alignment less than minLen
                        if int(alignment_len) < minLen:
                            continue
                        # 3. Build dict of longest aln len between two 
                        # elements
                        aln_pair = frozenset([query, subj])
                        if aln_pair in internal_aln_lens:
                            if int(alignment_len) > internal_aln_lens[aln_pair]:
                                internal_aln_lens[aln_pair] = int(alignment_len) 
                        else:
                            internal_aln_lens[aln_pair] = int(alignment_len) 

                    # 4. Ignore any alignment with 
                    # alnLn[i]/seqLen[i][j] < percAln for each seq in 
                    # each alignment
                    for aln_pair in internal_aln_lens:
                        aln_pair_lst = list(aln_pair)
                        el1 =  aln_pair_lst[0]
                        el2 =  aln_pair_lst[1]
                        el1_aln_ratio = internal_aln_lens[aln_pair] \
                                         /internal_seq_lengths[el1]*100
                        el2_aln_ratio = internal_aln_lens[aln_pair] \
                                         /internal_seq_lengths[el2]*100
                        if el1_aln_ratio < percAln or el2_aln_ratio < percAln:
                            continue
                        # 4. Build graph 
                        if el1 in G:
                            G[el1].add(el2)
                        else:
                            G[el1] = set([el2])
                        if el2 in G:
                            G[el2].add(el1)
                        else:
                            G[el2] = set([el1])

            if WICKERLTRS:
                # Extract all LTRs from CurrentGFF separated by classif
                paths['LTRs_{0}_GFF'.format(classif)] = '{0}/LTRs_{1}.gff'.format(
                                                              cOutDir, classif)
                paths['LTRs_{0}_Fasta'.format(classif)] = '{0}/LTRs_{1}.fasta'.format(
                                                              cOutDir, classif)
                writeLTRsGFF(paths['CurrentGFF'], paths['LTRs_{0}_GFF'.format(
                                                                     classif)], 
                                       elementSet=clusters_by_classif[classif])
                getfasta_call = [executables['bedtools'], 'getfasta', 
                                 '-fi', paths['inputFasta'], 
                                 '-s', 
                                 '-bed', paths['LTRs_{0}_GFF'.format(classif)]]
                getfasta_call_string = ('bedtools getfasta -fi {0} -s -bed {1} '
                                    '> {2} 2> {3}').format(paths['inputFasta'],
                                         paths['LTRs_{0}_GFF'.format(classif)], 
                                       paths['LTRs_{0}_Fasta'.format(classif)], 
                       '{0}/bedtools_getfasta_InternalRegions.stderr'.format(
                                                                      cOutDir))
                makecall(getfasta_call, stdout=paths['LTRs_{0}_Fasta'.format(
                                                                     classif)], 
                  stderr='{0}/bedtools_getfasta_InternalRegions.stderr'.format(
                                                                      cOutDir))
                ChangeFastaHeaders(paths['LTRs_{0}_Fasta'.format(classif)], 
                         paths['LTRs_{0}_GFF'.format(classif)], attribute='ID')

                # Get LTRs seq lengths
                ltrs_seq_lengths = {s.id:len(s) for s in list(SeqIO.parse(
                                                paths['LTRs_{0}_Fasta'.format(
                                                                     classif)], 
                                                'fasta'))}
                element_combined_ltr_lengths = {}
                for ltr in ltrs_seq_lengths:
                    el = ltr[:-2] # trim off the _L or _R
                    if el in element_combined_ltr_lengths:
                        element_combined_ltr_lengths[el] += ltrs_seq_lengths[ltr]
                    else:
                        element_combined_ltr_lengths[el] = ltrs_seq_lengths[ltr]
                elements = list(element_combined_ltr_lengths.keys())
                # blast all by all LTRs
                paths['LTRs_{0}_selfBlastnOut'.format(classif)] = (
                                      '{0}/LTRs_selfBlastn.tab').format(cOutDir)
                runblast(query=paths['LTRs_{0}_Fasta'.format(classif)], 
                         subject=paths['LTRs_{0}_Fasta'.format(classif)], 
                         out=paths['LTRs_{0}_selfBlastnOut'.format(classif)], 
                         outfmt='7', 
                         evalue='1e-5', 
                         percid=pId, 
                         blast='blastn', 
                         procs=procs)

                # Parse blast hits to define families
                LTRsBlastPth = paths['LTRs_{0}_selfBlastnOut'.format(classif)]
                with open(LTRsBlastPth, 'r') as blastFl:
                    ltr_aln_lens = {}
                    for line in blastFl:
                        if line.startswith('#'):
                            continue

                        (query, subj, percent_id, alignment_len, mismatches, 
                           gap_opens, q_start, q_end, s_start, s_end, E_value, 
                           bit_score) = line.strip().split('\t')

                        # 1. Ignore self-hits
                        if query == subj:
                            continue
                        # 2. Ignore any alignment less than minLen
                        if int(alignment_len) < minLen:
                            continue
                        # 3.1 Ignore any alignment with 
                        # alnLn[i]/seqLen[i][j] < percAln for each seq 
                        # in each alignment
                        # 3.2 Build dict of longest aln len between two 
                        # elements
                        aln_pair = frozenset([query, subj])
                        q_el = query[:-2]
                        s_el = subj[:-2]
                        if q_el in ltr_aln_lens:
                            if s_el in ltr_aln_lens[q_el]:
                                if aln_pair in ltr_aln_lens[q_el][s_el]:
                                    if (int(alignment_len) 
                                         > ltr_aln_lens[q_el][s_el][aln_pair]):
                                        ltr_aln_lens[q_el][s_el][aln_pair] = \
                                                             int(alignment_len) 
                                else:
                                    ltr_aln_lens[q_el][s_el][aln_pair] = \
                                                             int(alignment_len) 
                            else:
                                ltr_aln_lens[q_el][s_el] = {aln_pair:int(
                                                                alignment_len)}
                        else:
                            ltr_aln_lens[q_el] = {s_el:{aln_pair:int(
                                                               alignment_len)}}

                        if s_el in ltr_aln_lens:
                            if q_el in ltr_aln_lens[s_el]:
                                if aln_pair in ltr_aln_lens[s_el][q_el]:
                                    if (int(alignment_len) 
                                        > ltr_aln_lens[s_el][q_el][aln_pair]):
                                        ltr_aln_lens[s_el][q_el][aln_pair] = \
                                                             int(alignment_len) 
                                else:
                                    ltr_aln_lens[s_el][q_el][aln_pair] = \
                                                             int(alignment_len) 
                            else:
                                ltr_aln_lens[s_el][q_el] = {aln_pair:int(
                                                                alignment_len)}
                        else:
                            ltr_aln_lens[s_el] = {q_el:{aln_pair:int(
                                                               alignment_len)}}

                    for el1 in ltr_aln_lens:
                        for el2 in ltr_aln_lens[el1]:
                            pairs = ltr_aln_lens[el1][el2]
                            kinds = {}
                            for pair in pairs:
                                pLst = sorted(list(pair))
                                if (pLst[0].endswith('_L') 
                                    and pLst[1].endswith('_L')):
                                    if 'LL' in kinds:
                                        if pairs[pair] > kinds['LL']:
                                            kinds['LL'] = pairs[pair]
                                    else:
                                        kinds['LL'] = pairs[pair]

                                elif (pLst[0].endswith('_R') 
                                      and pLst[1].endswith('_R')):
                                    if 'RR' in kinds:
                                        if pairs[pair] > kinds['RR']:
                                            kinds['RR'] = pairs[pair]
                                    else:
                                        kinds['RR'] = pairs[pair]
                                elif (pLst[0].endswith('_L') 
                                      and pLst[1].endswith('_R')):
                                    if 'LR' in kinds:
                                        if pairs[pair] > kinds['LR']:
                                            kinds['LR'] = pairs[pair]
                                    else:
                                        kinds['LR'] = pairs[pair]
                                elif (pLst[0].endswith('_R') 
                                      and pLst[1].endswith('_L')):
                                    if 'RL' in kinds:
                                        if pairs[pair] > kinds['RL']:
                                            kinds['RL'] = pairs[pair]
                                    else:
                                        kinds['RL'] = pairs[pair]
                                else:
                                    sys.exit(
                                'Uknown element name LTR suffixes: {0}'.format( 
                                                                         pair))

                            len_both_ltrs = [element_combined_ltr_lengths[el1], 
                                             element_combined_ltr_lengths[el2]]
                            PASS = False
                            for l in len_both_ltrs:
                                ratio1 = 0
                                ratio2 = 0
                                ratio = 0
                                if 'LL' in kinds and 'RR' in kinds:
                                    ratio1 = ((kinds['LL'] + kinds['RR'])/l)*100
                                if 'LR' in kinds and 'RL' in kinds:
                                    ratio2 = ((kinds['LR'] + kinds['RL'])/l)*100
                                if not ratio1 == 0 and ratio2 == 0:
                                    ratio = max([ratio1, ratio2])
                                if ratio >= percAln:
                                    PASS = True
                            if PASS:
                                # 4. Build graph 
                                if el1 in G:
                                    G[el1].add(el2)
                                else:
                                    G[el1] = set([el2])
                                if el2 in G:
                                    G[el2].add(el1)
                                else:
                                    G[el2] = set([el1])

            # 5. Report nodes in each connected component as each family
            if G != {}:
                wicker_groups = graph2groups(G)
                group_num = max(list(wicker_groups.keys()))
                non_singletons = set()
                for g in wicker_groups:
                    for el in wicker_groups[g]:
                        non_singletons.add(el)

                for el in elements:
                    if el not in non_singletons:
                        group_num += 1
                        wicker_groups[group_num] = [el]

                paths['Wicker_{0}'.format(classif)] = (
                               '{0}/wicker_groups_{1}').format(cOutDir, classif)
                paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                              pId, percAln, minLen, classif)] = \
                                            paths['Wicker_{0}'.format(classif)]
                group_lens = sorted([(g,len(wicker_groups[g])) for 
                          g in wicker_groups], key=lambda x:x[1], reverse=True)
                # Write families
                with open(paths['Wicker_{0}'.format(classif)], 'w') as outFl:
                    for group in group_lens:
                        group = group[0]
                        outFl.write('\t'.join(wicker_groups[group])+'\n')
                if not checkStatusFl(
                      'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                               pId, percAln, minLen, classif)):
                    with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                        statusFlAppend.write("{0}\t{1}\n".format(
                      'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                                pId, percAln, minLen, classif), 
                                          paths['Wicker_{0}'.format(classif)]))
            else:
                # Write families: all singletons
                paths['Wicker_{0}'.format(classif)] = \
                               '{0}/wicker_groups_{1}'.format(cOutDir, classif)
                paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                              pId, percAln, minLen, classif)] = \
                                            paths['Wicker_{0}'.format(classif)]
                if os.path.isfile(paths['Wicker_{0}'.format(classif)]):
                    os.remove(paths['Wicker_{0}'.format(classif)])
                for el in elements:
                    with open(paths['Wicker_{0}'.format(classif)], 'a') as outFl:
                        outFl.write('{0}\n'.format(el))
                
                if not checkStatusFl(
                      'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                               pId, percAln, minLen, classif)):
                    with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                        statusFlAppend.write("{0}\t{1}\n".format(
                      'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                    pId, percAln, minLen, classif), paths[
                                                'Wicker_{0}'.format(classif)]))


def MCL(I=6, minClustSize=30, CombineIfTooFew=False):
    """CombineIfTooFew: Not available via command line flag.
    If there are less than minClustSize elements for a given 
    classification or for all classifications combined, they are put 
    into one cluster and MCL is not used.

    This function performs an all-by-all blast of elements within each
    classification and uses the results to separate the elements into
    clusters with MCL, following the protocol here:
        https://micans.org/mcl/man/clmprotocols.html#blast
    """
    global paths
    global filenames

    if USEMCL:
        if not 'LTRdigest_LTR_retrotransposons_fasta' in paths:

            MakeDir('MCLdir', '{0}/MCL'.format(paths['output_top_dir'], I))
            paths['LTRdigest_LTR_retrotransposons_GFF'] = (
                              '{0}/LTRdigest_LTR_retrotransposons.gff').format(
                                                         paths['GFFOutputDir'])
            paths['LTRdigest_LTR_retrotransposons_fasta'] = (
                            '{0}/LTRdigest_LTR_retrotransposons.fasta').format(
                                                       paths['FastaOutputDir'])
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'Began extracting LTR_retrotransposons from {0}'.format(
                                                          paths['CurrentGFF']))

            # Write GFF for just LTR_retrotransposon features
            with open(paths['CurrentGFF'], 'r') as gffFl:
                for line in gffFl:
                    if not line.startswith('#'):
                        gffLine = GFF3_line(line)
                        if gffLine.type == 'LTR_retrotransposon':
                            with open(paths['LTRdigest_LTR_retrotransposons_GFF'], 
                                    'a') as LTRdigest_LTR_retrotransposons_GFF:
                                LTRdigest_LTR_retrotransposons_GFF.write(
                                                  '{0}\n'.format(str(gffLine)))

            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           'Finished extracting LTR_retrotransposons from {0}'.format(
                                                          paths['CurrentGFF']))

            # Extract sequences for true positive LTR_retrotransposon features
            getfasta_ltrretrotransposons_call = [executables['bedtools'], 
                                                 'getfasta', 
                                                 '-fi', paths['inputFasta'], 
                                                 '-s', 
                                                 '-bed', '{0}'.format(
                                  paths['LTRdigest_LTR_retrotransposons_GFF'])]
            getfasta_ltrretrotransposons_call_string = \
                    'bedtools getfasta -fi {0} -s -bed {1} > {2} 2> {3}'.format(
                                    paths['inputFasta'], 
                                    paths['LTRdigest_LTR_retrotransposons_GFF'], 
                                    paths['LTRdigest_LTR_retrotransposons_fasta'], 
                                    '{0}/bedtools_getfasta.stderr'.format(
                                                      paths['FastaOutputDir']))
            append2logfile(paths['output_top_dir'], mainlogfile, ('Began extracting '
                                          'LTR_retrotransposon sequences from '
                                          'LTRdigest_TruePositives GFF:\n{0}').format(
                                     getfasta_ltrretrotransposons_call_string))
            makecall(getfasta_ltrretrotransposons_call, 
                     paths['LTRdigest_LTR_retrotransposons_fasta'],
                     '{0}/bedtools_getfasta.stderr'.format(
                                                      paths['FastaOutputDir']))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, 
                           ('Finished extracting LTR_retrotransposon sequences '
                                    'from LTRdigest_TruePositives GFF'))
            append2logfile(paths['output_top_dir'], 
                           mainlogfile, ('Changing FASTA headers from bedtools '
                                     'getfasta-style to LTR_retrotransposon ID'))
            ChangeFastaHeaders(paths['LTRdigest_LTR_retrotransposons_fasta'], 
                               paths['LTRdigest_LTR_retrotransposons_GFF'])
            append2logfile(paths['output_top_dir'], mainlogfile, 
                                    ('Done changing FASTA headers from bedtools '
                                     'getfasta-style to LTR_retrotransposon ID'))

            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write(
                        'LTRdigest_LTR_retrotransposons_fasta\t{0}\n'.format(
                                paths['LTRdigest_LTR_retrotransposons_fasta'])) 

        AllFASTA = list(SeqIO.parse(
                                paths['LTRdigest_LTR_retrotransposons_fasta'], 
                                'fasta'))
        classifFastas = {}

        MakeDir('MCLdir', '{0}/MCL'.format(paths['output_top_dir'], I))
        # If all elements combined are less than the min clust size 
        # specified by the user, then all elements are put into one cluster.
        # NOT ENABLED FOR USER. IT SHOULD WORK THOUGH.
        if CombineIfTooFew:
            if sum([len(clusters_by_classif[c]) for 
                                     c in clusters_by_classif]) < minClustSize:
                classif = 'All'
                allClassifs = '{0}/out.allClust'.format(paths['MCLdir'])
                c = []
                for cls in classifs:
                    c += clusters_by_classif[cls]
                    
                with open(allClassifs, 'w') as outFl:
                    outFl.write('{0}\n'.format('\t'.join(c)))
                
                newRecord = 'MCL_{0}_I{1}'.format(classif, I)
                if not checkStatusFl(newRecord):
                    paths['MCL_{0}_I{1}'.format(classif, I)] = (
                                     '{0}/out.allClust').format(paths['MCLdir'])
                    with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                        statusFlAppend.write('{0}\t{1}\n'.format(
                                             'MCL_{0}_I{1}'.format(classif, I), 
                                             paths['MCL_{0}_I{1}'.format(
                                                                 classif, I)])) 
                        return

        MakeDir('MCL_I{0}'.format(I), '{0}/I{1}'.format(paths['MCLdir'], I))
        if not checkStatusFl('MCL_I{0}'.format(I)):
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format('MCL_I{0}'.format(I), 
                                                         '{0}/I{1}'.format(
                                                          paths['MCLdir'], I))) 
        MakeDir('MCL_I{0}_ClustersDir'.format(I), '{0}/Clusters'.format(
                                                  paths['MCL_I{0}'.format(I)]))
        for classif in classifs:
            if not checkStatusFl('MCL_{0}_I{1}'.format(classif, I)):
                MakeDir('MCL_{0}_I{1}_dir'.format(classif, I), '{0}/{1}'.format(
                             paths['MCL_I{0}_ClustersDir'.format(I)], classif))
                outputPth = paths['MCL_{0}_I{1}_dir'.format(classif, I)]
                if len(clusters_by_classif[classif]) < minClustSize:
                    with open('{0}/{1}_MCL_clusters.I{2}'.format(outputPth, 
                                                    classif, I), 'w') as outFl:
                        outFl.write('{0}\n'.format('\t'.join(
                                                clusters_by_classif[classif])))

                    if 'MCL_{0}_I{1}'.format(classif, I) not in paths:
                        paths['MCL_{0}_I{1}'.format(classif, I)] = (
                                           '{0}/{1}_MCL_clusters.I{2}').format(
                                                         outputPth, classif, I)
                        with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                            statusFlAppend.write('{0}\t{1}\n'.format(
                                     'MCL_{0}_I{1}'.format(classif, I), 
                                     paths['MCL_{0}_I{1}'.format(classif, I)])) 

                if not 'MCL_{0}_abc'.format(classif) in paths:
                    classifFastas[classif] = '{0}/{1}.fasta'.format(
                                                            outputPth, classif)
                    SeqIO.write([rec for rec in AllFASTA if rec.id in 
                                                 clusters_by_classif[classif]], 
                                classifFastas[classif], 'fasta')
                    # make blast db for all-by-all blast
                    makeblastdb_AllByAll_call_string = (
                                '{0}/makeblastdb -in {1} -dbtype nucl').format(
                                  executables['blast'], classifFastas[classif])
                    makeblastdb_AllByAll_call = ['{0}/makeblastdb'.format(
                                                         executables['blast']), 
                                                         '-in', 
                                                         classifFastas[classif], 
                                                         '-dbtype', 'nucl']
                    scriptpath = os.path.realpath(__file__)
                    lineno = getframeinfo(currentframe()).lineno + 2
                    append2logfile(paths['output_top_dir'], 
                                mainlogfile, 
                                'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                   ('Began creating blast db for all-by-all '
                                        'blast of LTR_retrotransposon sequences '
                                        'using the call:\n{0}').format(
                                             makeblastdb_AllByAll_call_string))
                    makecall(makeblastdb_AllByAll_call, 
                            '{0}.makeblastdb.stdout'.format(
                                                        classifFastas[classif]), 
                            '{0}.makeblastdb.stderr'.format(
                                                        classifFastas[classif]))
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                  ('Finished creating blast db for all-by-all '
                         'blast of {0} LTR_retrotransposon sequences.').format(
                                                                      classif))
                    # Perform all-by-all blastn of true positive 
                    # LTR_retrotransposon sequences
                    paths['LTR_retrotransposonAllByAllblastnTable'] = (
        '{0}/{1}.LTRdigest_LTR_retrotransposon_{2}_AllByAll.blastn.tab').format(
                                   outputPth, filenames['inputFasta'], classif)
                    blastn_AllByAll_call = ['{0}/blastn'.format(executables['blast']), 
                                            '-db', classifFastas[classif], 
                                            '-query',  classifFastas[classif], 
                                            '-evalue', '1e-5', 
                                            '-outfmt', '6', 
                                            '-num_threads', str(procs), 
                                            '-max_hsps', '25']
                    blastn_AllByAll_call_string = ('{0} -db {1} -query {1} '
                         '-evalue 1e-5 -outfmt 6 -num_threads {2} -max_hsps 25 '
                           '1>{3} 2>{3}.stderr').format(
                                    '{0}/blastn'.format(executables['blast']),  
                                              classifFastas[classif], procs, 
                               paths['LTR_retrotransposonAllByAllblastnTable'])
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                   ('Began blastn all-by-all blast of '
                                      'LTR_retrotransposon sequences using the '
                                        'call:\n{0}').format(
                                                  blastn_AllByAll_call_string))
                    makecall(blastn_AllByAll_call, 
                             paths['LTR_retrotransposonAllByAllblastnTable'], 
                             '{0}.stderr'.format(
                              paths['LTR_retrotransposonAllByAllblastnTable']))
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                   ('Finished blastn all-by-all blast of {0} '
                               'LTR_retrotransposon sequences.').format(classif))

                    # Perform MCL clustering based on all-by-all blast 
                    # results convert blastn output  to abc format for mcl
                    paths['MCL_{0}_abc'.format(classif)] = (
                            '{0}/LTR_retrotransposon.AllByAllblastn.abc').format(
                                                                     outputPth)
                    blasttable2abc_call = ['cut', '-f', '1,2,11', 
                               paths['LTR_retrotransposonAllByAllblastnTable']]
                    blasttable2abc_call_string = ('cut -f 1,2,11 {0} > {1} '
                                                    '2>{1}.stderr').format(
                               paths['LTR_retrotransposonAllByAllblastnTable'], 
                               paths['MCL_{0}_abc'.format(classif)])
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                   ('Began converting blastn tabular output to '
                              'abc format for MCL using the call:\n{0}').format(
                                                   blasttable2abc_call_string))
                    makecall(blasttable2abc_call, 
                             paths['MCL_{0}_abc'.format(classif)], 
                             '{0}.stderr'.format(
                                         paths['MCL_{0}_abc'.format(classif)]))
                    append2logfile(paths['output_top_dir'], mainlogfile, 
             'Finished converting blastn tabular output to abc format for MCL')

                    with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                        statusFlAppend.write('MCL_{0}_abc\t{1}\n'.format(
                                         classif, 
                                         paths['MCL_{0}_abc'.format(classif)]))

                # create network for mcl
                paths['LTRretrotransposonNetworkmci'] = (
                                '{0}/LTR_retrotransposon_{1}.mci').format(
                                            paths['MCL_{0}_I{1}_dir'.format(
                                                         classif, I)], classif)
                paths['LTRretrotransposonNetworktab'] = (
                                    '{0}/LTR_retrotransposon_{1}.tab').format(
                                            paths['MCL_{0}_I{1}_dir'.format(
                                                         classif, I)], classif)
                mcxload_call = ['{0}/mcxload'.format(executables['mcl']), 
                                '-abc', paths['MCL_{0}_abc'.format(classif)], 
                                '--stream-mirror', 
                                '--stream-neg-log10', 
                                '-stream-tf', 'ceil(200)', 
                                '-o', paths['LTRretrotransposonNetworkmci'], 
                                '-write-tab', paths['LTRretrotransposonNetworktab']]
                mcxload_call_string = ('{0}/mcxload -abc {1} --stream-mirror '
                                         '--stream-neg-log10 -stream-tf ceil(200) '
                                         '-o {2} -write-tab {3} > {4}/mcxload.stdout '
                                         '2>{4}mcxload.stderr').format(
                                                            executables['mcl'], 
                                         paths['MCL_{0}_abc'.format(classif)], 
                                         paths['LTRretrotransposonNetworkmci'], 
                                         paths['LTRretrotransposonNetworktab'], 
                                         paths['MCL_{0}_I{1}_dir'.format(
                                                                  classif, I)])
                scriptpath = os.path.realpath(__file__)
                lineno = getframeinfo(currentframe()).lineno + 2
                append2logfile(paths['output_top_dir'], mainlogfile, 
                            'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
                append2logfile(paths['output_top_dir'], mainlogfile, 
                  'Began creating network for {0} using mcxload:\n{1}'.format(
                                                 classif, mcxload_call_string))
                makecall(mcxload_call, '{0}/mcxload.stdout'.format(
                                            paths['MCL_{0}_I{1}_dir'.format(
                                                                 classif, I)]), 
                         '{0}/mcxload.stderr'.format(
                         paths['MCL_{0}_I{1}_dir'.format(classif, I)]))
                append2logfile(paths['output_top_dir'], 
                               mainlogfile, 
                               ('Finished creating network for {0} LTR RTs '
                                    'using mcxload.').format(classif))
                
                # cluster using MCL
                current_wd = os.getcwd()
                mcl_call = ['{0}/mcl'.format(executables['mcl']), 
                            'LTR_retrotransposon_{0}.mci'.format(classif), 
                            '-I', str(I), 
                            '-use-tab', 
                            'LTR_retrotransposon_{0}.tab'.format(classif), 
                            '-te', str(procs) ]
                mcl_call_string = ('{0}/mcl {1} -I {2} -use-tab {3} -te {4} > '
                            'mcl.stdout 2>mcl.stderr').format(executables['mcl'], 
                              'LTR_retrotransposon_{0}.mci'.format(classif), I, 
                              'LTR_retrotransposons_{0}.tab'.format(classif), 
                              str(procs))
                os.chdir(current_wd)
                append2logfile(paths['output_top_dir'], 
                               mainlogfile, 
                              'Began clustering {0} using mcl:\n{1}'.format(
                                                     classif, mcl_call_string))
                os.chdir(paths['MCL_{0}_I{1}_dir'.format(classif, I)])
                makecall(mcl_call, 'mcl.stdout', 'mcl.stderr')
                os.chdir(current_wd)
                append2logfile(paths['output_top_dir'], 
                               mainlogfile, 
                               'Finished clustering {0} using mcl'.format(
                                                                      classif))
                os.chdir(paths['MCL_{0}_I{1}_dir'.format(classif, I)])
                os.chdir(current_wd)

                for fl in os.listdir(paths['MCL_{0}_I{1}_dir'.format(
                                                                 classif, I)]):
                    if fl.startswith('out'):
                        paths['MCL_{0}_I{1}'.format(classif, I)] = (
                             '{0}/{1}').format(paths['MCL_{0}_I{1}_dir'.format(
                                                              classif, I)], fl)
                        
                with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                    statusFlAppend.write('{0}\t{1}\n'.format(
                                     'MCL_{0}_I{1}'.format(classif, I), 
                                     paths['MCL_{0}_I{1}'.format(classif, I)])) 


def full2flankgff(inGFFpth, outGFFpth, bpflank):
    """Creates GFF3 files containing features of type
    LTR_RT_flank where each feature is the bpflank on either side of
    each feature of type repeat_region (the LTR-RT boundaries) in the
    input GFF3.

    inGFFpth       A LTRharvest/LTRdigest GFF3
    outGFFpth      The new flanking regions from elements in inGFFpth
    bpflank        The length of each flanking region
    """
    # overwrites an old file if present
    with open(outGFFpth, 'w') as outFl:
        outFl.write('##gff-version 3\n')
    with open(inGFFpth, 'r') as inGFFfl:
        outLines = set()
        for line in inGFFfl:
            if not line.startswith('#'):
                gffLine = GFF3_line(line)
                if gffLine.type == 'repeat_region':
                    SKIP_L = False
                    leftFlankGFFline = GFF3_line()
                    leftFlankGFFline.seqid = gffLine.seqid
                    leftFlankGFFline.source = 'full2flankgff'
                    leftFlankGFFline.type = 'LTR_RT_flank'
                    if int(gffLine.start) == 1:
                        SKIP_L = True
                        leftFlankGFFline.start = 1
                    else:
                        leftFlankGFFline.start = int(gffLine.start)-bpflank
                    if int(leftFlankGFFline.start) < 1:
                        leftFlankGFFline.start = 1
                    leftFlankGFFline.end = int(gffLine.start)-1
                    leftFlankGFFline.score = '.'
                    leftFlankGFFline.strand = gffLine.strand
                    leftFlankGFFline.phase = '.'
                    leftFlankGFFline.attributes['ID'] = (
                                    '{0}_leftflank_{1}_bp').format(
                                            gffLine.attributes['ID'], bpflank)
                    leftFlankGFFline.attributes['Parent'] = \
                                                       gffLine.attributes['ID']
                    leftFlankGFFline.attributes_order = ['ID', 'Parent']
                    leftFlankGFFline.refreshAttrStr()
                    if not SKIP_L:
                        outLines.add(str(leftFlankGFFline))
                    
                    rightFlankGFFline = GFF3_line()
                    rightFlankGFFline.seqid = gffLine.seqid
                    rightFlankGFFline.source = 'full2flankgff'
                    rightFlankGFFline.type = 'LTR_RT_flank'
                    rightFlankGFFline.start = int(gffLine.end)+1
                    rightFlankGFFline.end = int(gffLine.end)+bpflank
                    rightFlankGFFline.score = '.'
                    rightFlankGFFline.strand = gffLine.strand
                    rightFlankGFFline.phase = '.'
                    rightFlankGFFline.attributes['ID'] = (
                                    '{0}_rightflank_{1}_bp').format(
                                             gffLine.attributes['ID'], bpflank)
                    rightFlankGFFline.attributes['Parent'] = \
                                                       gffLine.attributes['ID']
                    rightFlankGFFline.attributes_order = ['ID', 'Parent']
                    rightFlankGFFline.refreshAttrStr()
                    outLines.add(str(rightFlankGFFline))
        with open(outGFFpth, 'w') as outFl:
            for line in outLines:
                outFl.write('{0}\n'.format(line)) 


def removeRedundant(fastaPth):
    """Removes all but the first occuring sequence if there are
    redundant sequences in this FASTA file"""
    seqnames = set()
    tmp = '{0}.nonredundant'.format(fastaPth)
    if os.path.isfile(tmp):
        os.remove(tmp)
    with open(tmp, 'w') as outFl:
        with open(fastaPth, 'r') as inFl:
            SKIP = False
            for line in inFl:
                if line.startswith('>'):
                    if line in seqnames:
                        SKIP = True
                        continue
                    else:
                        SKIP = False
                        seqnames.add(line)
                        outFl.write(line)
                else:
                    if SKIP:
                        continue
                    else:
                        outFl.write(line)
    os.rename(tmp, fastaPth)


def getfasta(inGFFpth, fastaRefPth, outFastaPth, headerKey='ID'):
    """Runs bedtools getfasta to extract sequences for features in a
    GFF3 file from the corresponding FASTA file and rename the sequences
    from the bedtools getfasta default to the value of a GFF3 attribute
    key.

    inGFFpth        the GFF3 file for which to extract sequences
    fastaRefPth     the FASTA file from which to extract sequences
    outFastaPth     the path which will become the FASTA file with
                    extracted sequences
    headerKey       the GFF3 attribute key from which the name of the
                    new FASTA sequences will be derived
    """
    call = [executables['bedtools'], 'getfasta', 
            '-fi', fastaRefPth, 
            '-s', 
            '-bed', inGFFpth]
    call_string = 'bedtools getfasta -fi {0} -s -bed {1} > {2}'.format(
                                            fastaRefPth, inGFFpth, outFastaPth)
    makecall(call, stdout=outFastaPth)
    ChangeFastaHeaders(outFastaPth, inGFFpth, attribute=headerKey)
    removeRedundant(outFastaPth)


def runblast(query, 
             subject, 
             out, 
             evalue, 
             outfmt, 
             percid=None, 
             blast='blastn', 
             procs=1):
    """Runs blastn (not tested with other blast programs), creating a
    blast database using makeblastdb.

    query       query sequences in FASTA format
    subject     sequences in FASTA format used to make a blast-
                searchable sequence database
    out         file to receive stdout of the blast call
    evalue      maximum evalue for an alignment to be reported
    outfmt      output format (6 is tabular, 7 is tabular with comments)
    percid      minimum percent identity for alignment to be reported
    blast       blast program to run. Only blastn has been tested
    procs       number of processors to use in parallel for blast run
    """
    # mapping of blast program to database type (the makeblastdb flag)
    dbtypes = {'blastn':'nucl', 
               'blastp':'prot', 
               'blastx':'prot', 
               'tblastx':'nucl', 
               'tblastn':'nucl'}
    assert blast in dbtypes,('runblast() blast param must be one of blastn, '
                             'blastp, blastx, tblastn, tblastx')
    # run makeblastdb to construct a blast database
    makeblastdb_call_string = ('{0}/makeblastdb -in {1} -dbtype {2} '
                                  '-parse_seqids').format(executables['blast'], 
                                                          subject, 
                                                          dbtypes[blast])
    makeblastdb_call = ['{0}/makeblastdb'.format(executables['blast']), 
                                                 '-in', subject, 
                                                 '-dbtype', dbtypes[blast], 
                                                 '-parse_seqids']
    makecall(makeblastdb_call)
    # blast programs other than blastn do not use the -percid flag, so
    # run them without it
    if percid == None or blast != 'blastn':
        blast_call = [ '{0}/{1}'.format(executables['blast'], blast), 
                       '-db', subject, 
                       '-query', query, 
                       '-evalue', str(evalue), 
                       '-outfmt', str(outfmt), 
                       '-num_threads', str(procs)]
        blast_call_string = ('{0}/{1} -db {2} -query {3} -evalue {4} -outfmt {5}'
                     ' -num_threads {6} 1>{7} 2>{8}').format(executables['blast'], 
                                                           blast, subject, query, 
                                                           evalue, outfmt, procs, 
                                                           out, 
                              '{0}/runblast.{1}.err'.format('/'.join(
                                                  out.split('/')[:-1]), blast))
    # blastn is called
    else:
        blast_call = ['{0}/{1}'.format(executables['blast'], blast), 
                                       '-db', subject, 
                                       '-query', query, 
                                       '-evalue', str(evalue), 
                                       '-outfmt', str(outfmt), 
                                       '-num_threads', str(procs), 
                                       '-perc_identity', str(percid)]
        blast_call_string = ('{0}/{1} -db {2} -query {3} -evalue {4} -outfmt {5}'
                    ' -num_threads {6} -perc_identity {7} 1>{8} 2>{9}').format(
                                           executables['blast'], blast, subject, 
                                           query, evalue, outfmt, procs, percid, 
                                           out, 
                              '{0}/runblast.{1}.err'.format('/'.join(
                                                  out.split('/')[:-1]), blast))
    append2logfile(paths['output_top_dir'], mainlogfile, blast_call_string)
    makecall(blast_call, stdout=out, stderr='{0}/runblast.{1}.err'.format(
                                         '/'.join(out.split('/')[:-1]), blast))


def reportpairswithhomologousflanks(blastqueryfasta, 
                                    blastResults, 
                                    outFlPth, 
                                    bpflank, 
                                    perc_len_cutoff):
    """Returns void; writes output file outFlPth, a two-column tab-
    delimited list of pairs of elements whose flanks are significantly
    homologous. All but one of these elements will be discarded from
    this set of elements.
    
    blastqueryfasta     sequences of the flanks of the elements used
                        in this homologous flank query
    blastResults        table of results from a blastn of flanks
                        from elements in blastqueryFasta
    outFlPth            path for the output file containing a list of
                        pairs of elements with significantly homologous
                        flanking sequences
    bpflank             the number of bp flanking each element queried
    perc_len_cutoff     the minimum percent of length of a flanka in the
                        blastn alignment for a pair of elements to be 
                        considered homologous flanks
    """
    # counts the number of flanks for each element in the input fasta
    # because elements near the ends of scaffold might have one flank
    # only and these are treated differently below
    flank_counts = {}
    with open(blastqueryfasta, 'r') as fastafl:
        for line in fastafl:
            if line.startswith('>'):
                el = line[1:].split('_')[1].lstrip('region')
                if el in flank_counts:
                    flank_counts[el] += 1
                else:
                    flank_counts[el] = 1
    # read blast results table and record pairs of elements whose
    # flanks were homologous with a minimum percent of the length of
    # a flanks perc_len_cutoff
    pairs = set()
    with open(blastResults, 'r') as blastfl:
        for line in blastfl:
            if not line.startswith('#'):
                contents = line.strip().split('\t')
                aln_len = int(contents[3])
                rpt1 = contents[0].split('_')[1].lstrip('region')
                rpt2 = contents[1].split('_')[1].lstrip('region')
                if rpt1 == rpt2:
                    continue
                else:
                    if 100*(aln_len/bpflank) >= float(perc_len_cutoff):
                        pairs.add(frozenset([contents[0], contents[1]]))
    # sort the pairs of homologous flanks so the output file is sorted
    # by element name
    pairs = sorted([sorted(list(pair)) for pair in pairs], key=lambda x:x[0])
    # for each pair of flank alignments record only nonredundant pairs
    # in the dictionary matches
    matches = {}
    for pair in pairs:
        # get the number of the repeats in this pair (e.g. if the
        # element name is repeat_region1 the number is 1
        rpt1 = pair[0].split('_')[1].lstrip('region')
        rpt2 = pair[1].split('_')[1].lstrip('region')
        # get the string indicating if the flank in this pair is the
        # left or the right flanks (leftflank, rightflank)
        flank1 = pair[0].split('_')[2]
        flank2 = pair[1].split('_')[2]
        if flank1 == 'leftflank':
            f1 = 'L1'
            F1 = 'L2'
        elif flank1 == 'rightflank':
            f1 = 'R1'
            F1 = 'R2'
        if flank2 == 'leftflank':
            f2 = 'L2'
            F2 = 'L1'
        elif flank2 == 'rightflank':
            f2 = 'R2'
            F2 = 'R1'
        match1 = frozenset([f1,f2])
        match2 = frozenset([F1,F2])
        if rpt1 in matches:
            if rpt2 in matches[rpt1]:
                matches[rpt1][rpt2].add(match1)
            else:
                matches[rpt1][rpt2] = set([match1])
        else:
            matches[rpt1] = {rpt2:set([match1])}

        if rpt2 in matches:
            if rpt1 in matches[rpt2]:
                matches[rpt2][rpt1].add(match2)
            else:
                matches[rpt2][rpt1] = set([match2])
        else:
            matches[rpt2] = {rpt1:set([match2])}
    options = [set([frozenset(['L1', 'L2']), 
                    frozenset(['R1','R2'])]), 
               set([frozenset(['L1', 'R2']), 
                    frozenset(['L2', 'R1'])])]
    # record pairs of elements where both left and right flanks were 
    # found to be homologous
    matchesset =  set()
    for rpt1 in matches:
        for rpt2 in matches[rpt1]:
            # this element only has one flank
            if flank_counts[rpt2] == 1:
                matchesset.add(frozenset([rpt1, rpt2]))
            else:
                if matches[rpt1][rpt2] in options:
                    matchesset.add(frozenset([rpt1, rpt2]))
    # if output file exists, remove it
    if os.path.isfile(outFlPth):
        os.remove(outFlPth)
    # write output file of element pairs whose left and right flanks
    # of length bpflank are homologous
    for match in matchesset:
        match = list(match)
        with open(outFlPth, 'a') as outFl:
            outFl.write('{0}\t{1}\n'.format('LTR_retrotransposon{0}'.format(
                                                                      match[0]), 
                                            'LTR_retrotransposon{0}'.format(
                                                                     match[1])))


def elementsWithHomologousFlanks(ingff, 
                                 infasta, 
                                 outdir, 
                                 bpflank=None, 
                                 outfmt='7', 
                                 percid=None, 
                                 evalue=None, 
                                 perc_len_cutoff=None, 
                                 procs=None):
    """Coordinates discovering element pairs whose left and right 
    homologous flanking regions are homologous. The reason this is done
    is so elements which were duplicated by non-transpositional events
    can be excluded from the phylogenies which are to be used for
    transposition rate analysis.

    ingff           the full LTRdigest+ GFF3
    infasta         the input genome
    outdir          location for outputfile containing list of elements
                    whose flanks are both homologous
    bpflank         the number of bp on each side of each element to
                    consider as flanks
    outfmt          blast output format (7 is commented tabular)
    percid          the minimum percent id in a blastn alignment for
                    a pair of flanks to be considered homologous
    evalue          maximum E-value of a blastn alignment for a pair
                    of flanks to be considered homologous
    perc_len_cutoff minimum percent of the length of a flank for it to
                    be considered for homology with other flanks
    procs           number of processors to use for blastn
    """
    # make input directory and prepare output filenames
    MakeDir('FlankDir', outdir)
    blastout = '{0}/blast_all_by_all_flanks'.format(outdir)
    ingffBasename = ingff.split('/')[-1].rstrip('.gff')
    flankfasta = '{0}/{1}.flanks.fasta'.format(outdir, ingffBasename)
    flankgff = '{0}/{1}.flanks.gff'.format(outdir, ingffBasename)
    # only run blastn on the flanking regions if the blast table does
    # not already exists
    if not os.path.isfile(blastout):
        if os.path.isfile(flankfasta):
            os.remove(flankfasta)
        if os.path.isfile(flankgff):
            os.remove(flankgff)
        # create a gff for all left and right flanks of each element in the
        # input gff (LTRdigest+ output)
        full2flankgff(inGFFpth=ingff, outGFFpth=flankgff, bpflank=bpflank)
        # extract sequences for each element's left and right flank
        getfasta(inGFFpth=flankgff, fastaRefPth=infasta, 
                                        outFastaPth=flankfasta, headerKey='ID')
        # run blastn of each flank against each other flank
        runblast(query=flankfasta, 
                 subject=flankfasta, 
                 out=blastout, 
                 evalue=str(evalue), 
                 outfmt=outfmt, 
                 percid=percid, 
                 blast='blastn', 
                 procs=procs)
    # write pairs of output files whose flanks are homologous
    reportfl = '{0}/LTR_element_pairs_with_homologous_flanks'.format(outdir)
    reportpairswithhomologousflanks(blastqueryfasta=flankfasta, 
                                                    blastResults=blastout, 
                                                    outFlPth=reportfl, 
                                                    bpflank=bpflank, 
                                                    perc_len_cutoff=perc_len_cutoff)


def CleanMafft(mafft_fasta):
    """Removes non-FASTA format text occurring before FASTA lines; some
    MAFFT versions output additional text.
    """
    if os.path.exists(mafft_fasta):
        with open('{0}.fixingmafftdefaultoutput.tmp'.format(
                                                   mafft_fasta), 'w') as outFl:
            with open(mafft_fasta, 'r') as inFl:
                STARTFASTA = False
                for line in inFl:
                    if line.startswith('>'):
                        STARTFASTA = True
                    if STARTFASTA == True:
                        outFl.write(line)
        copyfile('{0}.fixingmafftdefaultoutput.tmp'.format(mafft_fasta), mafft_fasta)
        os.remove('{0}.fixingmafftdefaultoutput.tmp'.format(mafft_fasta))


def aligner(elementList, OutDir, statusFlAlnKey, part):
    """Generic aligner to work with internal and whole elements from 
    LTRharvest GFF3. Called by AutoAlign()
    1. Writes GFF for elements requested (elementList)
    2. Extracts sequences using bedtools getfasta and converts the 
       headers to reflect the GFF3 ID attributes.
    3. Aligns using mafft with given settings.
    4. Trims alignment using trimal -automated1
    5. Adds statsFlAlnKey to the status file.
    """
    global paths

    if part == 'entire':
        ENTIRE = True
        INTERNAL = False
    elif part == 'internal':
        INTERNAL = True
        ENTIRE = True

    dirKey = '{0}-dir'.format(statusFlAlnKey)
    MakeDir(dirKey, OutDir)
    if not checkStatusFl(statusFlAlnKey):

        if INTERNAL:
            paths['InternalelementsGFF'] = '{0}/elements.gff'.format(OutDir)
            paths['InternalelementsFasta'] = '{1}/elements.fasta'.format(OutDir)
            writeLTRretrotransposonInternalRegions(paths['CurrentGFF'],
                                                   paths['InternalelementsGFF'], 
                                                   elementSet=set(elementList), 
                                                   truncateParent=True)
            getfasta_call = [executables['bedtools'], 'getfasta', 
                             '-fi', paths['inputFasta'], 
                             '-s', 
                             '-bed', paths['InternalelementsGFF']]
            getfasta_call_string = ('bedtools getfasta -fi {0} -s -bed {1} > '
                                    '{2} 2> {3}').format(paths['inputFasta'], 
                                                  paths['InternalelementsGFF'], 
                                                  paths['InternalelementsFasta'], 
                 '{0}/bedtools_getfasta_InternalRegions.stderr'.format(OutDir))
            makecall(getfasta_call, 
                     paths['InternalelementsFasta'], 
                 '{0}/bedtools_getfasta_InternalRegions.stderr'.format(OutDir))
            ChangeFastaHeaders(paths['InternalelementsFasta'], 
                               paths['InternalelementsGFF'], 
                               attribute = 'Parent')
            paths['AlnFasta'] = paths['InternalelementsFasta']
            paths['AlnPth'] = '{0}.aln'.format(paths['InternalelementsFasta'])

        elif ENTIRE:

            paths['EntireelementsGFF'] = '{0}/elements.gff'.format(OutDir)
            paths['EntireelementsFasta'] = '{0}/elements.fasta'.format(OutDir)
            writeLTRretrotransposonGFF(inputGFFpth = paths['CurrentGFF'], 
                                       outputGFFpth = paths['EntireelementsGFF'], 
                                       elementSet=set(elementList))
            getfasta_call = [executables['bedtools'], 'getfasta', 
                             '-fi', paths['inputFasta'], 
                             '-s', 
                             '-bed', paths['EntireelementsGFF']]
            getfasta_call_string = ('bedtools getfasta -fi {0} -s -bed {1} > '
                                    '{2} 2> {3}').format(paths['inputFasta'], 
                                                   paths['EntireelementsGFF'], 
                                                   paths['EntireelementsFasta'], 
                          '{0}/bedtools_getfasta_Entire.stderr'.format(OutDir))
            makecall(getfasta_call, paths['EntireelementsFasta'], 
                          '{0}/bedtools_getfasta_Entire.stderr'.format(OutDir))
            ChangeFastaHeaders(paths['EntireelementsFasta'], 
                               paths['EntireelementsGFF'], 
                               attribute='ID')
            paths['AlnFasta'] = paths['EntireelementsFasta']
            paths['AlnPth'] = '{0}.aln'.format(paths['EntireelementsFasta'])

        # Align regions from selected elements
        if len(elementList) <= mafft_smallAln_maxclustsize:
            mafft_call = [executables['mafft'], 
                          '--quiet', 
                          '--retree', '2', 
                          '--thread', str(procs), 
                          '--maxiterate', str(mafft_smallAln_maxiterate), 
                          paths['AlnFasta']]
        elif (len(elementList) > mafft_smallAln_maxclustsize 
               and len(elementList) <= mafft_mediumAln_maxclustsize):
            mafft_call = [executables['mafft'], 
                          '--quiet', 
                          '--retree', '2',
                          '--thread', str(procs), 
                          '--maxiterate', str(mafft_mediumAln_maxiterate), 
                          paths['AlnFasta']]
        elif (len(elementList) > mafft_mediumAln_maxclustsize 
               and len(elementList) <= mafft_largeAln_maxclustsize):
            mafft_call = [executables['mafft'], 
                         '--quiet', 
                         '--retree', '1', 
                         '--thread', str(procs), 
                         paths['AlnFasta']]
        elif len(elementList) > mafft_largeAln_maxclustsize:
            # Write would-be alignment path to status file
            paths['ClusterTrimmedAln'] = '{0}.trimal'.format(paths['AlnPth'])
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                paths[statusFlAlnKey] = paths['ClusterTrimmedAln']
                if os.path.isfile(paths[statusFlAlnKey]):
                    statusFlAppend.write('{0}\t{1}\n'.format(statusFlAlnKey, 
                                                   paths['ClusterTrimmedAln']))
            return
        mafft_call_string = '{0} {1}'.format(' '.join(mafft_call),
                                  ' >{0} 2>{0}.stderr'.format(paths['AlnPth']))
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], 
                       mainlogfile, 
                       'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                   scriptpath))
        append2logfile(paths['output_top_dir'], 
                       mainlogfile, 
                       'Aligning\n{0}'.format(mafft_call_string))
        makecall(mafft_call, paths['AlnPth'], 
                                          '{0}.stderr'.format(paths['AlnPth']))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                          'Cleaning MAFFT output\n{0}'.format(paths['AlnPth']))
        CleanMafft(paths['AlnPth'])
        append2logfile(paths['output_top_dir'], mainlogfile, 'Finished aligning')
        # Trim alignment
        paths['ClusterTrimmedAln'] = '{0}.trimal'.format(paths['AlnPth'])
        trimal_Aln_call =  [executables['trimal'], 
                            '-in', paths['AlnPth'], 
                            '-out', paths['ClusterTrimmedAln'], 
                            '-automated1']
        trimal_Aln_call_string = ('{0} -in {1} -out {2} -automated1 '
                        '>{2}.stdout 2>{2}.stderr').format(executables['trimal'], 
                                                    paths['AlnPth'], 
                                                    paths['ClusterTrimmedAln'])
        append2logfile(paths['output_top_dir'],
                       mainlogfile, 
                       'Began cleaning alignment using TrimAl:\n{0}'.format(
                                                       trimal_Aln_call_string))
        makecall(trimal_Aln_call, 
                 '{0}.stdout'.format(paths['ClusterTrimmedAln'], 
                                     '{0}.stdout'.format(
                                                  paths['ClusterTrimmedAln'])))
        append2logfile(paths['output_top_dir'], 
                       mainlogfile, 
                       'Finished cleaning alignment using TrimAl')
        paths[dirKey] = OutDir
        if os.path.isfile(paths['ClusterTrimmedAln']):
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                paths[statusFlAlnKey] = paths['ClusterTrimmedAln']
                if os.path.isfile(paths[statusFlAlnKey]):
                    statusFlAppend.write('{0}\t{1}\n'.format(statusFlAlnKey, 
                                                   paths['ClusterTrimmedAln']))


def AutoAlign(I=6, 
              part='entire', 
              rmgeneconv=False, 
              minClustSize=4, 
              align='clusters', 
              rmhomologflank=False, 
              clustering_method='WickerFam', 
              WickerParams={'pId':80,'percAln':80,'minLen':80}, 
              auto_outgroup=False, 
              bpflank=None, 
              flank_pId=None, 
              flank_evalue=None, 
              flank_plencutoff=None, 
              combine_and_do_small_clusters=True, 
              LTRSONLY=False):
    """AutoAlign handles aligning whole or internal regions of LTRharvest GFF3 
    LTR RTs that had been clustered using MCL() or WickerFam().

    Options:
    -------------
    I                           Inflation parameter used for MCL
    part                        'entire' or 'internal', referring to the 
                                part of the LTR RTs to align.
    rmgeneconv                  Remove elements within cluster that show 
                                evidence of inter-element gene conversion.
    minClustSize                The minimum cluster size to align. If 
                                combine_and_do_small_clusters=True, 
                                clusters smaller than minClustSize will 
                                be aligned together.
    align                       'clusters' (align clusters separately)
                                or 'classif' (align all clusters for a 
                                given superfamily together)
    rmhomologflank              Identify and remove one of every pair of 
                                elements from a given cluster whose 
                                flanking regions are homologous.
    clustering_method           'WickerFam' or 'MCL'. Clustering needs 
                                to have been done already.
    WickerParams                Parameters used with the WickerFam() 
                                clustering to be aligned.
    auto_outgroup               Automatically pick an ougroup and align 
                                with cluster. This option is not 
                                available for align='classif'
    bpflank                     If rmhomologflank=True, this specifies 
                                how many bp beyond the edge of each 
                                element to test for homology.
    flank_pId                   Minimum percent identity in flank 
                                alignment to consider as evidence of 
                                homology.
    flank_evalue                Minimum E-value of flank alignment to 
                                consider as evidence of homology.
    flank_plencutoff            Minimum percentage of flanking region 
                                required to participate in alignment 
                                to consider as evidence of homology.
    combine_and_do_small_clusters    Align all elements in clusters 
                                     smaller than minClustSize together.
    """
    global paths
    ENTIRE = False
    INTERNAL = False
    AUTO_OUTGROUP = auto_outgroup
    if align == 'clusters':
        ALIGNCLUSTERS = True
        ALIGNCLASSIFS = False
    elif align == 'classifs':
        ALIGNCLASSIFS = True
        ALIGNCLUSTERS = False
    else:
        sys.exit(("AutoAlign() parameter 'align' needs to be either "
                    "'clusters' or 'classifs', and it is {0}").format(align))
    if rmhomologflank:
        REMOVEHOMOLOGOUSFLANK = True
    else:
        REMOVEHOMOLOGOUSFLANK = False
    if rmgeneconv:
        REMOVEGENECONV = True
        gc = 'GCfiltered'
    else:
        REMOVEGENECONV = False
        gc = 'NoGCfiltering'
    if part == 'entire':
        ENTIRE = True
    elif part == 'internal':
        INTERNAL = True
    else:
        sys.exit(('part parameter to AutoAlign() is invalid: {0}. '
                'Must be either "entire" or "internal"').format(part))
    WICKERCLUST=False
    MCLCLUST=False
    if clustering_method == 'WickerFam':
        clusterDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])]
        MakeDir('WickerAlignments', '{0}/Alignments'.format(clusterDir))
        paths['Alignments'] = paths['WickerAlignments']
        WICKERCLUST=True
    elif clustering_method == 'MCL':
        clusterDir = paths['MCL_I{0}'.format(I)]
        MakeDir('MCLAlignments', '{0}/Alignments'.format(clusterDir))
        paths['Alignments'] = paths['MCLAlignments']
        MCLCLUST=True
    else:
        sys.exit(('AutoAlign() parameter clustering_method needs to be either '
                   'WickerFam or MCL, and it is {0}').format(clustering_method))
    if INTERNAL:
        MakeDir('InternalRegionsAlignments', 
                '{0}/InternalRegions'.format(paths['Alignments']))
        paths['AlnDir'] = paths['InternalRegionsAlignments']
    elif ENTIRE:
        MakeDir('WholeElementAlignments', 
                '{0}/WholeElements'.format(paths['Alignments']))
        paths['AlnDir'] = paths['WholeElementAlignments']
    if ALIGNCLASSIFS:
        MakeDir('Superfamilies', '{0}/Superfamilies'.format(paths['AlnDir']))
        if REMOVEHOMOLOGOUSFLANK:
            strHomoflank = 'homoflank'
            MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(
                                                       paths['Superfamilies']))
        else:
            strHomoflank = 'nohomoflank'
            MakeDir('HomoFlankDir', '{0}/AllElements'.format(
                                                       paths['Superfamilies']))
        paths['AlnDir'] = paths['HomoFlankDir']
    elif ALIGNCLUSTERS:
        if REMOVEGENECONV:
            MakeDir('GCDir', '{0}/GeneconversionDisallowed'.format(
                                                              paths['AlnDir']))
        else:
            MakeDir('GCDir', '{0}/NoGCFiltering'.format(paths['AlnDir']))
        if REMOVEHOMOLOGOUSFLANK:
            strHomoflank = 'homoflank'
            MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(
                                                               paths['GCDir']))
        else:
            strHomoflank = 'nohomoflank'
            MakeDir('HomoFlankDir', '{0}/AllElements'.format(paths['GCDir']))
        if AUTO_OUTGROUP:
            strOutgroup = 'withOutgroup'
            MakeDir('OutgroupDir', '{0}/WithOutgroup'.format(
                                                        paths['HomoFlankDir']))
        else:
            strOutgroup = 'noOutgroup'
            MakeDir('OutgroupDir', '{0}/NoOutgroup'.format(
                                                        paths['HomoFlankDir']))
        paths['AlnDir'] = paths['OutgroupDir']
    gcDct = {}
    if REMOVEGENECONV:
        if MCLCLUST:
            gcSummaryPth = 'MCL_I{0}_GENECONV_summary'.format(I)
        elif WICKERCLUST:
            gcSummaryPth = ('Wicker_{0}_pId_{1}_percAln_{2}_minLen'
                            '_GENECONV_summary').format(WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        if os.path.isfile(paths[gcSummaryPth]):
            with open(paths[gcSummaryPth], 'r') as gcFl:
                for line in gcFl:
                    el, clust, classif, g = line.strip().split()
                    el = 'LTR_retrotransposon{0}'.format(el)
                    clust = int(clust)
                    if classif in gcDct:
                        if clust in gcDct[classif]:
                            if not el in gcDct[classif][clust]:
                                gcDct[classif][clust].append(el)
                        else:
                            gcDct[classif][clust] = [el]
                    else:
                        gcDct[classif] = {clust:[el]}
        else:
            print(('modeltest(removegeneconv=True) used '
                   'but {0} not found').format(paths[gcSummaryPth]), 
                                               file=sys.stderr)
    for classif in classifs:
        if LTRSONLY:
            YESLTRS = False
            for SF in LTR_SFs:
                if classif.startswith(SF):
                    YESLTRS = True
                    break
            if not YESLTRS:
                continue
        if clustering_method == 'WickerFam':
            clusterPath = paths[('WickerFamDir_{0}_pId_{1}_'
                                 'percAln_{2}_minLen_{3}').format(
                                                      WickerParams['pId'], 
                                                      WickerParams['percAln'], 
                                                      WickerParams['minLen'], 
                                                      classif)]
        elif clustering_method == 'MCL':
            clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
        smallClusters = []
        MakeDir('AlnDir_{0}'.format(classif), '{0}/{1}'.format(paths['AlnDir'], 
                                                               classif))
        clusters = [clust.split('\t') for clust 
                           in open(clusterPath,'r').read().strip().split('\n')]
        if ALIGNCLASSIFS:
            if WICKERCLUST:
                statusFlKey = ('WickerAln_{0}_pId_{1}_percAln_'
                               '{2}_minLen_{3}_cluster_all.{4}.{5}').format(
                                                       WickerParams['pId'], 
                                                       WickerParams['percAln'], 
                                                       WickerParams['minLen'], 
                                                       classif, 
                                                       strHomoflank)
            if MCLCLUST:
                statusFlKey = 'Aln_{0}_I{1}_cluster_all.{2}.{3}'.format(
                                                                 classif, 
                                                                 I,
                                                                 strHomoflank)
            elementlist = [el for clust in clusters for el in clust]
            # No gene conversion implemented for this step
            if REMOVEHOMOLOGOUSFLANK:
                outDir = paths['AlnDir_{0}'.format(classif)]
                repeat_region_gff = '{0}/repeat_regions.gff'.format(outDir)
                elements = set()
                for el in elementlist:
                    el = el.lstrip('LTR_retrotransposon')
                    el = 'repeat_region{0}'.format(el)
                    elements.add(el)
                with open(paths['CurrentGFF'], 'r') as inFl:
                    with open(repeat_region_gff, 'w') as outFl:
                        for line in inFl:
                            if '\trepeat_region\t' in line:
                                gffLine = GFF3_line(line)
                                if gffLine.attributes['ID'] in elements:
                                    outFl.write(line)
                elementsWithHomologousFlanks(ingff = repeat_region_gff, 
                                             infasta = paths['inputFasta'], 
                                             outdir = outDir, 
                                             bpflank = bpflank, 
                                             outfmt = '7', 
                                             percid = flank_pId, 
                                             evalue = flank_evalue, 
                                             perc_len_cutoff = flank_plencutoff, 
                                             procs = procs)
                if 'LTR_element_pairs_with_homologous_flanks' in os.listdir(
                                                                       outDir):
                    print('Found some homologous flankings', file=sys.stderr)
                    print('{0}/LTR_element_pairs_with_homologous_flanks'.format(
                                                      outDir), file=sys.stderr)
                    # parse file and choose LTR RTs to leave out
                    with open(('{0}/LTR_element_pairs_with_'
                          'homologous_flanks').format(outDir), 'r') as pairsFl:
                        to_remove = [line.strip().split('\t')[0] for line 
                                                                    in pairsFl]
                        write2summary(('These elments were found have '
                                       'homologous flanking sequences with '
                                       'another {0} element and were therefore '
                                       'removed:\n{1}\n').format(classif, 
                                                         '\n'.join(to_remove)))
                        elementlist = [el for el 
                                         in elementlist if not el in to_remove]
                else:
                    print('Found no homologous flankings', file=sys.stderr)
                    print(outDir, file=sys.stderr)
                aligner(elementList = elementlist, 
                        OutDir = outDir, 
                        statusFlAlnKey = statusFlKey.format(classif), part = part)
        elif ALIGNCLUSTERS:
            if AUTO_OUTGROUP:
                # Write header for record of outgroups file, overwriting an old file if it exisits.
                if clustering_method == 'WickerFam':
                    OutgroupSummaryKey = ('WickerOutgroups_{0}_pId_{1}_percAln'
                            '_{2}_minLen_{3}_{4}.{5}.OutgroupSummary').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        gc, 
                                                        strHomoflank)
                elif clustering_method == 'MCL':
                    OutgroupSummaryKey = ('MCLOutgroups_{0}_I{1}_{2}.{3}.'
                                          'OutgroupSummary').format(classif, 
                                                                   I, 
                                                                   gc, 
                                                                   strHomoflank)
                outgroupFile = '{0}/Outgroups'.format(
                                                paths['AlnDir_{0}'.format(
                                                                     classif)])
                paths[OutgroupSummaryKey] = outgroupFile
                if not checkStatusFl(OutgroupSummaryKey):
                    with open(outgroupFile, 'w') as outFl:
                        outFl.write('Alignment\toutgroup\toutgroup_cluster\n')
            for j in range(len(clusters)):
                elementlist = clusters[j]
                if REMOVEGENECONV:
                    if os.path.isfile(paths[gcSummaryPth]):
                        try:
                            if j in gcDct[classif]:
                                elementlist = [el for el in clusters[j] 
                                                if not el in gcDct[classif][j]]
                        except KeyError:
                            pass
                    else:
                        pass # No gene conversion was found at all
                else:
                    pass

                if len(elementlist) < minClustSize: # Put small clusters together and align at the end
                    if combine_and_do_small_clusters:
                        smallClusters += elementlist
                else:
                    MakeDir('AlnDir_{0}_I{1}-cluster{2}_{3}'.format(classif, 
                                                                    I, 
                                                                    j, 
                                                                    gc), 
                                    '{0}/cluster_{1}'.format(
                                       paths['AlnDir_{0}'.format(classif)], j))
                    outDir = paths['AlnDir_{0}_I{1}-cluster{2}_{3}'.format(
                                                                       classif, 
                                                                       I, 
                                                                       j, 
                                                                       gc)]
                    if REMOVEHOMOLOGOUSFLANK:
                        repeat_region_gff = '{0}/repeat_regions.gff'.format(
                                                                        outDir)
                        elements = set()
                        for el in elementlist:
                            el = el.lstrip('LTR_retrotransposon')
                            el = 'repeat_region{0}'.format(el)
                            elements.add(el)
                        with open(paths['CurrentGFF'], 'r') as inFl:
                            with open(repeat_region_gff, 'w') as outFl:
                                for line in inFl:
                                    if '\trepeat_region\t' in line:
                                        gffLine = GFF3_line(line)
                                        if gffLine.attributes['ID'] in elements:
                                            outFl.write(line)
                        elementsWithHomologousFlanks(ingff = repeat_region_gff, 
                                                     infasta = paths['inputFasta'], 
                                                     outdir = outDir, 
                                                     bpflank = bpflank, 
                                                     outfmt = '7', 
                                                     percid = flank_pId, 
                                                     evalue = flank_evalue, 
                                                     perc_len_cutoff = flank_plencutoff, 
                                                     procs = procs)
                        if ('LTR_element_pairs_with_homologous_flanks' 
                                                        in os.listdir(outDir)):
                            write2summary(('Found some homologous flankings '
                             'among cluster {0}, clustering {1}').format(j, 
                                                            clustering_method))
                            # parse file and choose LTR RTs to leave out
                            with open(('{0}/LTR_element_pairs_with_'
                                       'homologous_flanks').format(
                                                      outDir), 'r') as pairsFl:
                                to_remove = [line.strip().split('\t')[0] for 
                                                               line in pairsFl]
                                write2summary(('These elments were found have '
                                               'homologous flanking sequences '
                                               'with another {0} element and '
                                               'were therefore '
                                               'removed:\n{1}\n').format(
                                                         classif, 
                                                         '\n'.join(to_remove)))
                                elementlist = [el for el in elementlist if 
                                                           not el in to_remove]
                        else:
                            write2summary(('Found no homologous flankings '
                                'among cluster {0}, clustering {1}').format(j, 
                                                            clustering_method))
                    if clustering_method == 'WickerFam':
                        statusFlKey = ('WickerAln_{0}_pId_{1}_percAln_{2}_'
                               'minLen_{3}_cluster_{4}_{5}.{6}.{7}').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        j, 
                                                        gc, 
                                                        strHomoflank, 
                                                        strOutgroup)
                    elif clustering_method == 'MCL':
                        statusFlKey = ('Aln_{0}_I{1}_cluster{2}_'
                                       '{3}.{4}.{5}').format(classif, 
                                                             I, 
                                                             j, 
                                                             gc, 
                                                             strHomoflank, 
                                                             strOutgroup)
                    if not statusFlKey in paths:
                        if AUTO_OUTGROUP:
                            # the outgroup shall be a random element 
                            # from cluster k where cluster k is the 
                            # largest of the clusters that is not j if
                            # j is the first cluster then the next 
                            # smallest cluster is k if there is no other 
                            # cluster, no outgroup is used
                            OUTGROUP_POSSIBLE = True
                            if j==0:
                                try:
                                    clusters[j+1]
                                    k = j+1
                                except KeyError:
                                    OUTGROUP_POSSIBLE = False
                                    print(('statusFlKey: {0}\nj: {1}\nclusters:'
                                            ' {2}\nclusters_len: {3}\n').format(
                                                                 statusFlKey, 
                                                                 j, 
                                                                 clusters, 
                                                                 len(clusters)))
                                    print('KeyError, line 2178')
                                    sys.exit()
                            elif j > 0:
                                k = 0
                            if OUTGROUP_POSSIBLE:
                                outgroup = random.choice(clusters[k])
                                clustAndOutgroup = clusters[j] + [outgroup]
                                with open(outgroupFile, 'a') as outFl:
                                    outFl.write('{0}\t{1}\t{2}\n'.format(
                                                                  statusFlKey, 
                                                                  outgroup, 
                                                                  k))
                                elementlist = clustAndOutgroup
                        aligner(elementlist, 
                                OutDir = outDir, 
                                statusFlAlnKey = statusFlKey, 
                                part = part)
            if AUTO_OUTGROUP:
                if not checkStatusFl(OutgroupSummaryKey):
                    with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                        if os.path.isfile(paths[OutgroupSummaryKey]):
                            statusFlAppend.write('{0}\t{1}\n'.format(
                                                   OutgroupSummaryKey, 
                                                   paths[OutgroupSummaryKey]))
            if combine_and_do_small_clusters:
                if len(smallClusters) > 1:
                    MakeDir('AlnDir_{0}_I{1}-cluster{2}_{3}'.format(classif, 
                                                                    I, 
                                                                    'small', 
                                                                    gc), 
                            '{0}/cluster_{1}'.format(paths['AlnDir_{0}'.format(
                                                           classif)], 'small'))
                    outDir = paths['AlnDir_{0}_I{1}-cluster{2}_{3}'.format(
                                                                       classif, 
                                                                       I, 
                                                                       'small', 
                                                                       gc)]
                    if clustering_method == 'WickerFam':
                        statusFlKey = ('WickerAln_{0}_pId_{1}_percAln_{2}_'
                                       'minLen_{3}_clustersmall_{4}').format(
                                                       WickerParams['pId'], 
                                                       WickerParams['percAln'], 
                                                       WickerParams['minLen'], 
                                                       classif, 
                                                       gc)
                    elif clustering_method == 'MCL':
                        statusFlKey = 'Aln_{0}_I{1}_clustersmall_{2}'.format(
                                                                classif, I, gc)
                    if not statusFlKey in paths:
                        aligner(smallClusters, 
                                OutDir = outDir, 
                                statusFlAlnKey = statusFlKey, 
                                part = part)


def geneconvClusters(g = '/g0', 
                     clust = None, 
                     I = 6, 
                     minClustSize = 4, 
                     clustering_method = 'WickerFam', 
                     WickerParams = {'pId':80,'percAln':80,'minLen':80}, 
                     combine_and_do_small_clusters = True, 
                     LTRSONLY = True):
    """Searches each cluster for evidence of inter-element gene
    conversion using GENECONV.
    
    g can be one of /g0, /g1, or /g2
    g is proportional to the tolerance for mismatches in fragments by 
    GENECONV.
    """
    global paths
    global filenames
    if GENECONVCLUSTERS:
        WICKERCLUST = False
        MCLCLUST = False
        if clustering_method == 'WickerFam':
            WICKERCLUST = True
        elif clustering_method == 'MCL':
            MCLCLUST = True
        else:
            sys.exit(('geneconvClusters() parameter clustering_method needs to'
                      ' be either WickerFam or MCL, and it is: {0}').format(
                                                            clustering_method))
        if WICKERCLUST:
            WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                       WickerParams['pId'], 
                                                       WickerParams['percAln'], 
                                                       WickerParams['minLen'])]
            paths['Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(
                         WickerParams['pId'], 
                         WickerParams['percAln'], 
                         WickerParams['minLen'])] = '{0}/GENECONV'.format(
                                                                     WickerDir)
            WickerGCdirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(
                                                      WickerParams['pId'], 
                                                      WickerParams['percAln'], 
                                                      WickerParams['minLen'])
            if not checkStatusFl(WickerGCdirkey):
                with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                    statusFlAppend.write('{0}\t{1}\n'.format(
                                                        WickerGCdirkey, 
                                                        paths[WickerGCdirkey]))
            geneconvOutputDir = WickerGCdirkey
            gcSummaryFl = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_{3}_summary'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'], 
                                                    g[1:])
        elif MCLCLUST:
            MCLdir = paths['MCL_I{0}'.format(I)]
            paths['MCL_I{0}_GENECONVdir'.format(I)] = '{0}/GENECONV'.format(
                                                                        MCLdir)
            if not checkStatusFl('MCL_I{0}_GENECONVdir'.format(I)):
                with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
                    statusFlAppend.write('{0}\t{1}\n'.format(
                                      'MCL_I{0}_GENECONVdir'.format(I), 
                                      paths['MCL_I{0}_GENECONVdir'.format(I)]))
            geneconvOutputDir = 'MCL_I{0}_GENECONVdir'.format(I)
            gcSummaryFl = 'MCL_I{0}_GENECONV_{1}_summary'.format(I, g[1:])
        if not gcSummaryFl in paths:
            geneconv_calls = []
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], mainlogfile, 
                            'Below log entry is from line {0} in {1}'.format(
                                                                   lineno, 
                                                                   scriptpath))
            append2logfile(paths['output_top_dir'], mainlogfile, 
              'Checking directory structure for GENECONV using {0}'.format(g))

            for classif in classifs:
                # Align only clusters from superfamilies with identical 
                # LTRs on transposition (Copia, Gypsy, ERV, BEL/Pao)
                if LTRSONLY:
                    YESLTRS = False
                    for SF in LTR_SFs:
                        if classif.startswith(SF):
                            YESLTRS = True
                            break
                if not YESLTRS:
                    continue

                MakeDir('GENECONV_{0}_dir'.format(classif), 
                           '{0}/{1}'.format(paths[geneconvOutputDir], classif))
                MakeDir('GENECONV_{0}_{1}_dir'.format(classif, g[1:]), 
                       '{0}/{1}'.format(paths['GENECONV_{0}_dir'.format(
                                                             classif)], g[1:]))
                if MCLCLUST:
                    clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
                elif WICKERCLUST:
                    clusterPath = paths[('WickerFamDir_{0}_pId_{1}_percAln'
                                                       '_{2}_minLen_{3}').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif)]
                clusters = [clust.split('\t') for clust in open(
                                   clusterPath,'r').read().strip().split('\n')]
                smalls = []
                for j in range(len(clusters)):
                    if len(clusters[j]) >= minClustSize:
                        if MCLCLUST:
                            alnPth = paths[('Aln_{0}_I{1}_cluster{2}_NoGCfiltering'
                                '.nohomoflank.noOutgroup').format(classif, I, j)]
                        elif WICKERCLUST:
                            alnPth = paths[('WickerAln_{0}_pId_{1}_percAln_{2}'
                                        '_minLen_{3}_cluster_{4}_NoGCfiltering.'
                                        'nohomoflank.noOutgroup').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        j)]
                        if os.path.isfile(alnPth):
                            # if alignment file is non-empty
                            if not os.stat(alnPth).st_size == 0:
                                scriptpath = os.path.realpath(__file__)
                                lineno = getframeinfo(currentframe()).lineno + 2
                                append2logfile(paths['output_top_dir'], 
                                               mainlogfile, 
                                  'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
                                append2logfile(paths['output_top_dir'], 
                                               mainlogfile, 
                               'Preparing to run GENECONV on:\n{0}'.format(alnPth))
                                call = [executables['geneconv'], 
                                        alnPth, 
                                        '/w124', 
                                        g, 
                                        '-include_monosites', 
                                        '-nolog', 
                                        '-Dumptab', 
                                        '-Fancy']
                                geneconv_calls.append(
                                               (call, '/dev/null', None, None))
                        else:
                            continue
                    else:
                        if combine_and_do_small_clusters:
                            smalls += clusters[j]
                if combine_and_do_small_clusters:
                    if smalls != []:
                        if MCLCLUST:
                            if 'Aln_{0}_I{1}_clustersmall_NoGCfiltering'.format(
                                                          classif, I) in paths:
                                alnPth = paths[('Aln_{0}_I{1}_clustersmall_'
                                             'NoGCfiltering').format(classif, I)]
                            else:
                                continue
                        elif WICKERCLUST:
                            if (('WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}'
                                 '_clustersmall_NoGCfiltering').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif) in paths):
                                alnPth = paths[('WickerAln_{0}_pId_{1}_percAln'
                                                '_{2}_minLen_{3}_clustersmall_'
                                                'NoGCfiltering').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif)]
                            else:
                                continue
                        if os.path.isfile(alnPth):
                            # if alignment file is non-empty
                            if not os.stat(alnPth).st_size == 0:
                                scriptpath = os.path.realpath(__file__)
                                lineno = getframeinfo(currentframe()).lineno + 2
                                append2logfile(paths['output_top_dir'], 
                                               mainlogfile, 
                                  'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
                                append2logfile(paths['output_top_dir'], 
                                               mainlogfile, 
                                  'Preparing to run GENECONV on:\n{0}'.format(
                                                                       alnPth))
                                call = [executables['geneconv'], 
                                        alnPth,
                                        '/w124', 
                                        g,
                                        '-include_monosites',
                                        '-nolog',
                                        '-Dumptab',
                                        '-Fancy']
                                geneconv_calls.append((call, '/dev/null', 
                                                '{0}.geneconv.err'.format(
                                                                alnPth), None))
                        else:
                            continue
            if geneconv_calls == []:
                return
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], mainlogfile, 
                  'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
            append2logfile(paths['output_top_dir'], mainlogfile, 
                  'Running GENECONV, example call:\n{0}'.format(
                                              ' '.join(geneconv_calls[0][0])))
            chunk_size = ceil(len(geneconv_calls)/procs)
            with Pool(processes=procs) as p:
                p.map(makecallMultiprocessing, geneconv_calls, 
                                                          chunksize=chunk_size)
            p.join()

            hasEvidenceOfGC = set()
            for classif in classifs:
                # Align only clusters from superfamilies with identical 
                # LTRs on transposition (Copia, Gypsy, ERV, BEL/Pao)
                if LTRSONLY:
                    YESLTRS = False
                    for SF in LTR_SFs:
                        if classif.startswith(SF):
                            YESLTRS = True
                            break
                if not YESLTRS:
                    continue
                if MCLCLUST:
                    clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
                elif WICKERCLUST:
                    clusterPath = paths[('WickerFamDir_{0}_pId_{1}_percAln_{2}'
                                        '_minLen_{3}').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif)]
                clusters = [clust.split('\t') for clust 
                              in open(clusterPath,'r').read().strip().split('\n')]
                for j in range(len(clusters)):
                    if len(clusters[j]) >= minClustSize: 
                        if MCLCLUST:
                            trimalOutput = paths[('Aln_{0}_I{1}_cluster{2}_'
                                'NoGCfiltering.nohomoflank.noOutgroup').format(
                                                     classif, I, j)].split('.')
                        elif WICKERCLUST:
                            trimalOutput = paths[('WickerAln_{0}_pId_{1}_'
                                              'percAln_{2}_minLen_{3}_'
                                              'cluster_{4}_NoGCfiltering.'
                                              'nohomoflank.noOutgroup').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        j)].split('.')
                        # GENECONV runs and creates the file were the alignment is
                        gcOutput = copy(trimalOutput)
                        gcOutput[-1] = 'tab'
                        geneconvOutputPth = '.'.join(gcOutput)
                        if not os.path.isfile(geneconvOutputPth):
                            continue
                        with open(geneconvOutputPth, 'r') as gcFl:
                            for line in gcFl:
                                if line.startswith('GI'):
    ## Add parsing of new format here and make Circos plots.
#    #   Seq       Sim     BC KA    Aligned Offsets         In Seq1            In Seq2        Num  Num  Tot  MisM
#    #   Names    Pvalue   Pvalue   Begin  End   Len    Begin  End   Len   Begin  End   Len   Poly Dif  Difs Pen.
#    GI  S18;S37  0.0000  4.73e-41   235   1605 1371     235   1605 1371    229   1599 1371    213   0  229  None
#    GI  S21;S37  0.0000  3.32e-38   235   1383 1149     235   1383 1149    229   1377 1149    179   0  249  None
#    GI  S34;S37  0.0000  6.18e-37   322   1605 1284     322   1605 1284    316   1599 1284    206   0  218  None
#    GI  S26;S37  0.0000  1.33e-33   426   1383  958     426   1383  958    420   1377  958    157   0  252  None
                                    GI, 
                                    els, 
                                    sim_p_val, 
                                    bc_ka_p_val, 
                                    aln_start, 
                                    aln_end, 
                                    aln_len, 
                                    el1_start, 
                                    el1_end, 
                                    el1_len, 
                                    el2_start, 
                                    el2_end, 
                                    el2_len, 
                                    num_polys, 
                                    num_difs, 
                                    tot_difs, 
                                    mism_pen  = line.strip().split('\t')
                                    element1, element2 = [i[1:] for i in 
                                                                els.split(';')]
                                    if int(tot_difs) < 3:
                                        continue
                                    hasEvidenceOfGC.add((element1, j, 
                                                               classif, g[1:]))
                                    hasEvidenceOfGC.add((element2, j, 
                                                               classif, g[1:]))
                                    line='{0}\t{1}\t{2}\n'.format(
                                                      line.strip(), classif, j)
                                    with open('{0}/{1}_{2}.summary'.format(
                                       paths['GENECONV_{0}_dir'.format(
                                                             classif)],
                                                classif, g[1:]), 'a') as outFl:
                                        outFl.write(line)

                        fName = geneconvOutputPth.split('/')[-1]
                        os.rename(geneconvOutputPth, '{0}/clust{1}_{2}'.format(
                           paths['GENECONV_{0}_{1}_dir'.format(
                                                   classif, g[1:])], j, fName))
                trimalOutput = None
                if combine_and_do_small_clusters:
                    if MCLCLUST:
                        if ('Aln_{0}_I{1}_clustersmall_NoGCfiltering'.format(
                                                          classif, I) in paths):
                            if os.path.isfile(paths[('Aln_{0}_I{1}_clustersmall'
                                          '_NoGCfiltering').format(classif, I)]):
                                trimalOutput = paths[('Aln_{0}_I{1}_clustersmall'
                                 '_NoGCfiltering').format(classif, I)].split('.')
                    elif WICKERCLUST:
                        if ('WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_'
                               'clustersmall_NoGCfiltering').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'],
                                                        WickerParams['minLen'], 
                                                        classif) in paths:
                            # 1 sequence, no alignment
                            if os.path.isfile(paths[('WickerAln_{0}_pId_{1}_'
                                    'percAln_{2}_minLen_{3}_clustersmall_'
                                    'NoGCfiltering').format(
                                                       WickerParams['pId'], 
                                                       WickerParams['percAln'], 
                                                       WickerParams['minLen'], 
                                                       classif)]):
                                trimalOutput = paths[('WickerAln_{0}_pId_{1}_'
                                    'percAln_{2}_minLen_{3}_clustersmall_'
                                    'NoGCfiltering').format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif)].split('.')
                if not trimalOutput == None:
                    # GENECONV runs and creates the file were the alignment is
                    gcOutput = copy(trimalOutput)
                    gcOutput[-1] = 'tab'
                    geneconvOutputPth = '.'.join(gcOutput)
                    if not os.path.isfile(geneconvOutputPth):
                        continue
                    with open(geneconvOutputPth, 'r') as gcFl:
                        for line in gcFl:
                            if line.startswith('GI'):
                                totDifs = int(line.strip().split()[9])
                                if totDifs < 5:
                                    continue
                                element1 = line.strip().split()[1].split(';')[0]
                                element2 = line.strip().split()[1].split(';')[1]
                                hasEvidenceOfGC.add(
                                      (element1.lstrip('S'), j, classif, g[1:]))
                                hasEvidenceOfGC.add(
                                      (element2.lstrip('S'), j, classif, g[1:]))
                                line='{0}\t{1}\t{2}\n'.format(
                                                      line.strip(), classif, j)
                                with open('{0}/{1}_{2}.summary'.format(
                                            paths['GENECONV_{0}_dir'.format(
                                                                classif)], 
                                                                classif, 
                                                                g[1:]), 
                                                                'a') as outFl:
                                    outFl.write(line)
                    fName = geneconvOutputPth.split('/')[-1]
                    os.rename(geneconvOutputPth, 
                              '{0}/clust{1}_{2}'.format(
                                    paths['GENECONV_{0}_{1}_dir'.format(
                                                              classif, g[1:])], 
                                    j, fName))
            for el in sorted(list(
                            hasEvidenceOfGC), key=lambda x:(x[2], x[1], x[3])):
                with open('{0}/ElementsWithEvidenceOfGeneConversion'.format(
                                      paths[geneconvOutputDir]), 'a') as outFl:
                    outFl.write('{0}\t{1}\t{2}\t{3}\n'.format(
                                                   el[0], el[1], el[2], g[1:]))
            paths['GENECONV_clusters_I{0}_summary'.format(I)] = ('{0}/'
                            'ElementsWithEvidenceOfGeneConversion').format(
                                                      paths[geneconvOutputDir])
            if MCLCLUST:
                paths['MCL_I{0}_GENECONV_summary'.format(I)] = ('{0}/'
                                 'ElementsWithEvidenceOfGeneConversion').format(
                                                      paths[geneconvOutputDir])
            elif WICKERCLUST:
                paths['Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_summary'.format(
                                           WickerParams['pId'], 
                                           WickerParams['percAln'], 
                                           WickerParams['minLen'])] = ('{0}/'
                                 'ElementsWithEvidenceOfGeneConversion').format(
                                                      paths[geneconvOutputDir])
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                if MCLCLUST:
                    if not checkStatusFl('MCL_I{0}_GENECONV_{1}_summary'.format(
                                                                    I, g[1:])):
                        statusFlAppend.write('{0}\t{1}\n'.format(
                            'MCL_I{0}_GENECONV_{1}_summary'.format(I, g[1:]), 
                            '{0}/ElementsWithEvidenceOfGeneConversion'.format(
                                                    paths[geneconvOutputDir])))
                    if not checkStatusFl('MCL_I{0}_GENECONV_summary'.format(I)):
                        statusFlAppend.write('{0}\t{1}\n'.format(
                            'MCL_I{0}_GENECONV_summary'.format(I), 
                            '{0}/ElementsWithEvidenceOfGeneConversion'.format(
                                                    paths[geneconvOutputDir])))
                elif WICKERCLUST:
                    if not checkStatusFl(('Wicker_{0}_pId_{1}_percAln_{2}_'
                                          'minLen_GENECONV_{3}_summary').format(
                                                       WickerParams['pId'], 
                                                       WickerParams['percAln'], 
                                                       WickerParams['minLen'], 
                                                       g[1:])):
                        statusFlAppend.write('{0}\t{1}\n'.format(
                            ('Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_'
                             '{3}_summary').format(WickerParams['pId'], 
                                                   WickerParams['percAln'], 
                                                   WickerParams['minLen'], 
                                                   g[1:]), 
                           '{0}/ElementsWithEvidenceOfGeneConversion'.format(
                                                    paths[geneconvOutputDir])))
                    if not checkStatusFl(('Wicker_{0}_pId_{1}_percAln_{2}_'
                            'minLen_GENECONV_summary').format(
                                                        WickerParams['pId'],
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])):
                        statusFlAppend.write('{0}\t{1}\n'.format(
                            ('Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_'
                             'summary').format(WickerParams['pId'], 
                                               WickerParams['percAln'], 
                                               WickerParams['minLen']), 
                            '{0}/ElementsWithEvidenceOfGeneConversion'.format(
                                                    paths[geneconvOutputDir])))
        

def modeltest(iters=1, 
              I=6, 
              removegeneconv=True, 
              part='entire', 
              clustering_method='WickerFam', 
              WickerParams={'pId':80,'percAln':80,'minLen':80}, 
              minClustSize=4, 
              bpflank=None, 
              combine_and_do_small_clusters=True):
    """Runs jModeltest2 on each cluster and writes information for PAUP
    to summary files.
    """
    global paths
    global fileames
    REMOVEGENECONV = removegeneconv
    WICKERCLUST = False
    MCLCLUST = False
    if clustering_method == 'WickerFam':
        WICKERCLUST = True
    elif clustering_method == 'MCL':
        MCLCLUST = True
    else:
        sys.exit(('modeltest() parameter clustering_method needs to be either'
            ' WickerFam or MCL, and it is: {0}').format(clustering_method))
    if REMOVEGENECONV:
        gc = 'GCfiltered'
    else:
        gc = 'NoGCfiltering'
    keySuffix = '{0}.nohomoflank.noOutgroup'.format(gc)
    if LTRDIVERGENCE:
        if WICKERCLUST:
            WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])]
            WickerMTdirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_modeltestDir'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
            paths[WickerMTdirkey] = '{0}/Modeltest'.format(WickerDir)
            paths['ModelTestDir'] = paths[WickerMTdirkey]
            if not checkStatusFl(WickerMTdirkey):
                with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                    statusFlAppend.write('{0}\t{1}\n'.format(
                                        WickerMTdirkey, paths[WickerMTdirkey]))
            if REMOVEGENECONV:
                AutoAlign(I=None, 
                          part='entire', 
                          rmgeneconv=True, 
                          minClustSize=minClustSize, 
                          align='clusters', 
                          rmhomologflank=False, 
                          clustering_method='WickerFam', 
                          WickerParams={'pId':WickerParams['pId'], 
                          'percAln':WickerParams['percAln'], 
                          'minLen':WickerParams['minLen']}, 
                          auto_outgroup=False, 
                          bpflank=bpflank, 
                          combine_and_do_small_clusters=combine_and_do_small_clusters, 
                          flank_pId=flank_pId, 
                          flank_evalue=flank_evalue, 
                          flank_plencutoff=flank_plencutoff, 
                          LTRSONLY=True)
                MakeDir('ModelTestDir_{0}'.format(keySuffix), 
                    '{0}/GeneconversionDisallowed'.format(paths['ModelTestDir']))
            else:
                AutoAlign(I=None, 
                          part='entire', 
                          rmgeneconv=False, 
                          minClustSize=minClustSize, 
                          align='clusters', 
                          rmhomologflank=False, 
                          clustering_method='WickerFam', 
                          WickerParams={'pId':WickerParams['pId'], 
                          'percAln':WickerParams['percAln'], 
                          'minLen':WickerParams['minLen']}, 
                          auto_outgroup=False, 
                          bpflank=bpflank, 
                          combine_and_do_small_clusters=combine_and_do_small_clusters, 
                          flank_pId=flank_pId, 
                          flank_evalue=flank_evalue, 
                          flank_plencutoff=flank_plencutoff, 
                          LTRSONLY=True)
                MakeDir('ModelTestDir_{0}'.format(keySuffix), 
                             '{0}/NoGCFiltering'.format(paths['ModelTestDir']))
            # model testing not done for small clusters
            alignmentsForModeltesting = [pth for pth in paths if pth.startswith(
                        'WickerAln_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])) 
                     and pth.endswith('{0}.nohomoflank.noOutgroup'.format(gc))]
        elif MCLCLUST:
            MCLdir = paths['MCL_I{0}'.format(I)]
            MCL_MT_dirkey = 'MCL_I{0}_modeltestDir'.format(I)
            paths[MCL_MT_dirkey] = '{0}/Modeltest'.format(MCLdir)
            paths['ModelTestDir'] = paths[MCL_MT_dirkey]
            if not checkStatusFl(MCL_MT_dirkey):
                with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                    statusFlAppend.write('{0}\t{1}\n'.format(
                                          MCL_MT_dirkey, paths[MCL_MT_dirkey]))
            if REMOVEGENECONV:
                AutoAlign(I=I, 
                          part=part, 
                          rmgeneconv=True, 
                          minClustSize=minClustSize, 
                          align='clusters', 
                          rmhomologflank=False, 
                          clustering_method='MCL', 
                          WickerParams=None, 
                          auto_outgroup=False, 
                          bpflank=bpflank, 
                          combine_and_do_small_clusters=combine_and_do_small_clusters, 
                          flank_pId=flank_pId, 
                          flank_evalue=flank_evalue, 
                          flank_plencutoff=flank_plencutoff, 
                          LTRSONLY=True)
                MakeDir('ModelTestDir_{0}'.format(keySuffix), 
                  '{0}/GeneconversionDisallowed'.format(paths['ModelTestDir']))
            else:
                AutoAlign(I=I, 
                          part=part, 
                          rmgeneconv=False, 
                          minClustSize=minClustSize, 
                          align='clusters', 
                          rmhomologflank=False, 
                          clustering_method='MCL', 
                          WickerParams=None, 
                          auto_outgroup=False, 
                          bpflank=bpflank, 
                          combine_and_do_small_clusters=combine_and_do_small_clusters, 
                          flank_pId=flank_pId, 
                          flank_evalue=flank_evalue, 
                          flank_plencutoff=flank_plencutoff, 
                          LTRSONLY=True)
                MakeDir('ModelTestDir_{0}'.format(keySuffix), 
                             '{0}/NoGCFiltering'.format(paths['ModelTestDir']))
            alignmentsForModeltesting = [pth for pth in paths if 
                                            pth.startswith('Aln_') and 
                                            pth.endswith('{0}.nohomoflank.noOutgroup'.format(
                                              gc)) and 'I{0}'.format(I) in pth]
        if alignmentsForModeltesting == []:
            with open('{0}/NO_MODEL_TESTING_DONE'.format(
                   paths['ModelTestDir_{0}'.format(keySuffix)]), 'w') as outFl:
                outFl.write('Perhaps there are too few LTR RTs')
        for aln in alignmentsForModeltesting:
            if not os.path.isfile(paths[aln]):
                continue
            if MCLCLUST:
                # e.g. aln = Aln_Other_I8_cluster0_GCfiltered.nohomoflank.noOutgroup
                a = aln.split('_')
                classif = a[1]
                if '-' in a[2]:
                    j = a[2].split('-')[1][7:]
                    suffix = a[3] 
                else:
                    j = a[3][7:]
                    suffix = a[4]
                OutDirKey = 'MCLModelTestDir_{0}_iters_I{1}_{2}_cluster_{3}_{4}'.format(
                                                  iters,I, classif, j, suffix)
                OutSummaryKey = 'MCLModelTestSummary_{0}_iters_I{1}_{2}'.format(
                                                              iters,I, classif)
                if j == 'small':
                    if not combine_and_do_small_clusters:    
                        continue
                    j = 'clustersmall'
            elif WICKERCLUST:
                # e.g. 
                # aln = WickerAln_80_pId_80_percAln_80_minLen_Copia_cluster_0_GCfiltered
                a = aln.split('_')
                classif = a[7]
                j = a[-2] # Cluster
                suffix = a[-1]
                OutDirKey = ('WickerModelTestDir_{0}_iters_{1}_pId_{2}_percAln'
                            '_{3}_minLen_{4}_cluster_{5}_{6}').format(
                                                        iters, 
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        j, 
                                                        suffix)
                OutSummaryKey = ('WickerModelTestSummary_{0}_iters_{1}_pId_{2}'
                                '_percAln_{3}_minLen_{4}').format(
                                                        iters, 
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif)
            MakeDir('ModelTest_{0}_{1}_dir'.format(suffix, classif), 
                    '{0}/{1}'.format(paths['ModelTestDir_{0}'.format(suffix)], 
                                                                      classif))
            sessionDir = paths['ModelTest_{0}_{1}_dir'.format(suffix, classif)]
            # If true it means the model testing already finished for that group
            if not OutDirKey in paths:

                MakeDir(OutDirKey, '{0}/{1}_iters'.format(sessionDir, iters))
                MakeDir('ModelTestClustDir', '{0}/cluster_{1}'.format(
                                                          paths[OutDirKey], j))
                for i in range(int(iters)):
                    MakeDir('ModelTestIterationDir', 
                            '{0}/iter_{1}'.format(
                                         paths['ModelTestClustDir'], str(i+1)))
                    # FastTree
                    paths['Tree'] = '{0}/{1}_I{2}_{3}.tree'.format(
                                 paths['ModelTestIterationDir'], classif, I, j)
                    filenames['Tree'] = '{0}_I{1}_{2}.tree'.format(
                                                                 classif, I, j)
                    fasttree_call = [ executables['fasttree'], '-nt', '-gtr' ]
                    fasttree_call_string = '{0} -nt -gtr <{1} >{2} 2>{2}.stderr'.format(
                             executables['fasttree'], paths[aln],paths['Tree'])
                    scriptpath = os.path.realpath(__file__)
                    lineno = getframeinfo(currentframe()).lineno + 2
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                   'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
                    append2logfile(paths['output_top_dir'],
                                   mainlogfile, 
                                   'Began inferring phylogeny using FastTree:\n{0}'.format(
                                                         fasttree_call_string))
                    makecall(fasttree_call, stdout=paths['Tree'], 
                            stderr='{0}.stderr'.format(paths['Tree']), 
                            stdin=paths[aln])
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                   'Finished inferring phylogeny using FastTree')
                    paths['jModeltest2out_{0}'.format(
                            classif)] = '{0}/{1}.jModelTest2.out'.format(
                             paths['ModelTestIterationDir'], filenames['Tree'])
                    jmodeltestCallString = ('java -jar {0} -d {1} -w -g 4 -f '
                        '-BIC -a -u {2} -o {3} -tr {4} -s 11').format(
                                                 executables['jmodeltest2'], 
                                                 paths[aln], 
                                                 paths['Tree'], 
                                                 paths['jModeltest2out_{0}'.format(
                                                         classif)], str(procs))
                    append2logfile(paths['output_top_dir'], 
                                   mainlogfile, 
                                   'Starting jModeltest2\n{0}'.format(
                                                         jmodeltestCallString))
                    subprocess.call(['java', '-jar', executables['jmodeltest2'], 
                                     '-d', 
                                     paths[aln], 
                                     '-w', 
                                     '-g', 
                                     '4', 
                                     '-f', 
                                     '-BIC', 
                                     '-a', 
                                     '-u', 
                                     paths['Tree'], 
                                     '-o', 
                                     paths['jModeltest2out_{0}'.format(classif)], 
                                     '-tr', 
                                     str(procs), 
                                     '-s', 
                                     '11'])
                                    
                    topDir = paths[OutDirKey]
                    jmt2summaryFlPth = '{0}/Summary.txt'.format(topDir)
                    paths['jModeltest2summary_{0}'.format(classif)] = jmt2summaryFlPth
                    clustSize = len([ 1 for line in open(
                                paths[aln],'r').read().split('\n') 
                                if line.startswith('>') ])
                    with open(paths['jModeltest2out_{0}'.format(classif)], 'r') as jmt2outFl:
                        END = False
                        PAUP = False
                        getPAUP = False
                        COMPLETE_RUN = False
                        paupLines = ''
                        lastline = ''
                        for line in jmt2outFl:
                            if line.startswith('PAUP* Commands Block:'):
                                PAUP = True
                            elif line.startswith('::Best Models::'):
                                END = True
                            if line.startswith('END;') and getPAUP:
                                paupLines += 'Dset distance=ML;\n'
                                paupLines += 'SaveDist format=oneColumn file=dist;\n'
                                paupLines += line
                                getPAUP = False
                            elif line.startswith('[!') and PAUP:
                                paupLines += line
                                PAUP = False
                                getPAUP = True
                            elif getPAUP:
                                paupLines += line
                            if line.startswith('BIC') and END:
                                with open(jmt2summaryFlPth, 'a') as summaryFl:
                                    summaryFl.write('{0}\t{1}\t{2}\n'.format(
                                                   line.strip(), j, clustSize))
                                    summaryFl.write(paupLines)
                                    COMPLETE_RUN = True
            if not checkStatusFl(OutDirKey):
                with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                    statusFlAppend.write('{0}\t{1}\n'.format(
                                    OutDirKey, paths['ModelTestIterationDir']))
                    if 'jModeltest2summary_{0}'.format(
                                            classif) in paths and COMPLETE_RUN:
                        if not checkStatusFl(OutSummaryKey):
                            statusFlAppend.write('{0}\t{1}\n'.format(
                              OutSummaryKey, 
                              paths['jModeltest2summary_{0}'.format(classif)]))
                    else:
                        pass


def align_ltrs(I=6, 
               clustering_method='WickerFam', 
               WickerParams={'pId':80,'percAln':80,'minLen':80}, 
               DONTALIGN=False):
    """Make a GFF3 for every LTR pair from which to use to extract 
    sequences and align each pair of LTRs using MAFFT.

    Elements with Superfamily classifications besides the ones in 
    LTR_SFs are ignored because they might be DIRS or another kind of 
    SF with non-identical LTRs, thus the LTR divergence process would 
    not be accurate.
    """
    global paths
    global filenames
    WICKERCLUST = False
    MCLCLUST = False
    # And set up directory structure for output (LTR pairs GFFs, 
    # FASTAs, and alignments)
    if clustering_method == 'WickerFam':
        WICKERCLUST = True
        key_base = 'WickerFam_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
        WickerDir = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
        MakeDir('Alignments', '{0}/Alignments'.format(paths[WickerDir]))
        AlnKey = '{0}.LTRAlnDir'.format(key_base)
        paths[AlnKey] = '{0}/LTRs'.format(paths['Alignments'])
        MakeDir(AlnKey, paths[AlnKey])
        GFFKey = '{0}.GFFDir'.format(key_base)
        paths[GFFKey] = '{0}/GFFs'.format(paths[WickerDir])
        MakeDir(GFFKey, paths[GFFKey])
        FASTAKey = '{0}.FASTADir'.format(key_base)
        paths[FASTAKey] = '{0}/FASTAs'.format(paths[WickerDir])
        MakeDir(FASTAKey, paths[FASTAKey])
        TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
    elif clustering_method == 'MCL':
        MCLCLUST = True
        key_base = 'MCL_I{0}'.format(I)
        MCLdir = 'MCL_I{0}'.format(I)
        MakeDir('Alignments', '{0}/Alignments'.format(paths[MCLdir]))
        AlnKey = '{0}.LTRAlnDir'.format(key_base)
        paths[AlnKey] = '{0}/LTRs'.format(paths['Alignments'])
        MakeDir(AlnKey, paths[AlnKey])
        GFFKey = '{0}.GFFDir'.format(key_base)
        paths[GFFKey] = '{0}/GFFs'.format(paths[MCLdir])
        MakeDir(GFFKey, paths[GFFKey])
        FASTAKey = '{0}.FASTADir'.format(key_base)
        paths[FASTAKey] = '{0}/FASTAs'.format(paths[MCLdir])
        MakeDir(FASTAKey, paths[FASTAKey])
        TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
    else:
        sys.exit(('align_ltrs() parameter clustering_method needs to be either '
                  'WickerFam or MCL and it is: {0}').format(clustering_method))
    if checkStatusFl('{0}.LTR_divergence_complete'.format(key_base)):
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
                'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                 scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                'ltr_divergence() already completed: {0}'.format(
                        paths['{0}.LTR_divergence_complete'.format(key_base)]))
        return
    if GENECONVLTRS or LTRDIVERGENCE:
        MakeDir('LTRsGFFoutputDir', '{0}/LTRs'.format(paths['GFFOutputDir']))
        MakeDir('LTRsFASTAoutputDir', '{0}/LTRs'.format(paths['FastaOutputDir']))
        MakeDir('Alignments', '{0}/Alignments'.format(paths['output_top_dir']))
        MakeDir('LTRsAlignments_dir', '{0}/LTRs'.format(paths['Alignments']))
        # this will have the GFF lines for each LTR RT's LTRs
        ltrs = {}
        ltrs_getfasta_calls = {}
        ltrs_changefastaheaders_calls = {}
        ltrs_mafft_calls = {}
        ltrs_clean_mafft_output_calls = {}
        ltrs_trimal_calls = {}
        num_pairs = 0
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
                  'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                   scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                    'Parsing LTRs from GFF3:\n{0}'.format(paths['CurrentGFF']))
        with open(paths['CurrentGFF'], 'r') as GFF_fl:
            for line in GFF_fl:
                if '\tlong_terminal_repeat\t' in line:
                    gffLine = GFF3_line(line)
                    elementName = gffLine.attributes['Parent']
                    # Skip element kinds with potentially non-identical 
                    # LTRs upon insertion, e.g. DIRS, Ngaro
                    YESLTR = False
                    for SF in LTR_SFs:
                        if classifs_by_element[elementName].startswith(SF):
                            YESLTR = True
                            break
                    if YESLTR:
                        if elementName in ltrs:
                            ltrs[elementName].append(line)
                            if len(ltrs[elementName]) == 2:
                                LTRsGFFfilepath = '{0}/{1}_LTRs.gff'.format(
                                                            paths[GFFKey], 
                                                            elementName)
                                LTRsFASTAfilepath = '{0}/{1}_LTRs.fasta'.format(
                                                            paths[FASTAKey], 
                                                            elementName)
                                LTRsAlignmentFilepath = '{0}/{1}_LTRs.fasta.aln'.format(
                                                            paths[AlnKey], 
                                                            elementName)
                                LTRsTrimmedAlnFilepath = '{0}.trimmed'.format(
                                                         LTRsAlignmentFilepath)
                                # For writing GFF3
                                ltrs[elementName].append(LTRsGFFfilepath)
                                getfasta_LTRs_call = [executables['bedtools'], 
                                                     'getfasta', 
                                                     '-fi', paths['inputFasta'], 
                                                     '-s', 
                                                     '-bed', LTRsGFFfilepath]
                                ltrs_getfasta_calls[elementName] = (
                                                    getfasta_LTRs_call, 
                                                    LTRsFASTAfilepath, 
                                                    None, None)
                                ltrs_changefastaheaders_calls[elementName] = (
                                                    LTRsFASTAfilepath, 
                                                    LTRsGFFfilepath, 
                                                    'Parent')
                                mafft_LTRs_call = [executables['mafft'], 
                                                  '--quiet', 
                                                  '--globalpair', 
                                                  '--maxiterate', '1000', 
                                                  LTRsFASTAfilepath]
                                ltrs_mafft_calls[elementName] = (mafft_LTRs_call, 
                                                    LTRsAlignmentFilepath, 
                                                    None, None)
                                ltrs_clean_mafft_output_calls[
                                           elementName] = LTRsAlignmentFilepath
                                trimal_LTRs_call = [executables['trimal'], 
                                                   '-in', LTRsAlignmentFilepath, 
                                                   '-out', LTRsTrimmedAlnFilepath, 
                                                   '-automated1']
                                ltrs_trimal_calls[elementName] = (
                                                        trimal_LTRs_call, 
                                                        LTRsTrimmedAlnFilepath, 
                                                        None, None)
                                num_pairs += 1
                        else:
                            ltrs[elementName] = [line]
        chunk_size = ceil(num_pairs/procs)
        # Write to log file here about chunk size and processors used
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], 
                       mainlogfile, 
                       'Below log entry is from line {0} in {1}'.format(
                                                                lineno, 
                                                                scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                        'For align_ltrs(): procs={0} chunk_size={1}'.format(
                                                                  procs,
                                                                  chunk_size))
        if not checkStatusFl(GFFKey):
            append2logfile(paths['output_top_dir'], mainlogfile, 
                            'Writing GFF3s for each LTR pair:\n{0}'.format(
                                                                paths[GFFKey]))
            # Write GFF3s for each LTR pair
            with Pool(processes = procs) as p:
                p.map(write_ltrs_gff3, ltrs.values(), chunksize=chunk_size)
            p.join()
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                # Add LTRs GFF path to status file (for resuming later)
                statusFlAppend.write('{0}\t{1}\n'.format(GFFKey, paths[GFFKey]))
        if not checkStatusFl(FASTAKey):
            append2logfile(paths['output_top_dir'], mainlogfile, 
                            'Extracting sequences for LTR pairs:\n{0}'.format(
                                        list(ltrs_getfasta_calls.values())[0]))
            # Write FASTA for each LTR pair
            with Pool(processes = procs) as p:
                p.map(makecallMultiprocessing, ltrs_getfasta_calls.values(), 
                                                        chunksize = chunk_size)
            p.join()
            with open('{0}/status'.format(paths['output_top_dir']),
                                                        'a') as statusFlAppend:
                # Add LTRs FASTA path to status file (for resuming later)
                statusFlAppend.write('{0}\t{1}\n'.format(FASTAKey, 
                                                         paths[FASTAKey]))
        NewHeaderKey = '{0}.LTRsFASTAnewheaders'.format(key_base)
        if not checkStatusFl(NewHeaderKey):
            append2logfile(paths['output_top_dir'], mainlogfile, 
                'Changing bedtools getfasta default headers to LTR RT names:\n{0}'.format(
                                                              paths[FASTAKey]))
            # Write FASTA for each LTR pair
            with Pool(processes = procs) as p:
                p.map(ChangeFastaHeadersMultiprocessing,
                      ltrs_changefastaheaders_calls.values(), 
                      chunksize = chunk_size)
            p.join()
            paths[NewHeaderKey] = paths[FASTAKey]
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                # Add LTRs FASTA path to status file (for resuming later)
                statusFlAppend.write('{0}\t{1}\n'.format(NewHeaderKey, 
                                                         paths[FASTAKey]))
        # Stop here. Happens when only doing solo LTR search.
        if DONTALIGN:
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], mainlogfile, 
                            'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
            append2logfile(paths['output_top_dir'], mainlogfile, 
             'align_ltrs() stopping after writing FASTAs for each LTR'.format(
                                                             procs,chunk_size))
            return
        if not checkStatusFl(AlnKey):
            append2logfile(paths['output_top_dir'], mainlogfile, 
                    'Making MAFFT alignments for each LTR pair:\n{0}'.format(
                                           list(ltrs_mafft_calls.values())[0]))
            # Do alignment for each LTR pair
            with Pool(processes=procs) as p:
                p.map(makecallMultiprocessing, ltrs_mafft_calls.values(), 
                                                          chunksize=chunk_size)
            p.join()
            # Remove non-fasta format text from alignment if present 
            # (some versions of MAFFT output aln method info with alignments
            with Pool(processes=procs) as p:
                p.map(CleanMafft, ltrs_clean_mafft_output_calls.values(), 
                                                          chunksize=chunk_size)
            p.join()
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                # Add LTRs FASTA path to status file (for resuming later)
                statusFlAppend.write('{0}\t{1}\n'.format(AlnKey, paths[AlnKey]))
        if not checkStatusFl(TrimalKey):
            append2logfile(paths['output_top_dir'], mainlogfile, 
              'Trimming MAFFT LTRs alignment using TrimAl -automated1:\n{0}'.format(
                                         list(ltrs_trimal_calls.values())[0]))
            # Do alignment for each LTR pair
            with Pool(processes=procs) as p:
                p.map(makecallMultiprocessing, ltrs_trimal_calls.values(), 
                                                          chunksize=chunk_size)
            p.join()
            paths[TrimalKey] = paths[AlnKey]
            with open('{0}/status'.format(paths['output_top_dir']), 
                                                        'a') as statusFlAppend:
                # Add LTRs FASTA path to status file (for resuming later)
                statusFlAppend.write('{0}\t{1}\n'.format(TrimalKey, paths[
                                                                   TrimalKey]))
                statusFlAppend.write('{0}.LTR_divergence_complete\t{0}'.format(
                                                                     key_base))


def SoloLTRsearch(I=6, 
                  clustering_method='WickerFam', 
                  WickerParams={'pId':80,'percAln':80,'minLen':80}):
    """Finds LTRs which are homologous to an LTR from a full-length
    LTR-RT but are not associated with any of the full-length LTR-RTs
    and do not overlap other solo LTRs. The blast highest scoring solo
    LTR is retained.
    """
    global paths
    global filenames
    # Set up directory structure for output (LTRs GFF & FASTA)
    if clustering_method == 'WickerFam':
        WICKERCLUST = True
        clustMethodTopDir = '{0}/WickerFamDir/{1}_pId_{2}_percAln_{3}_minLen'.format(
                                                    paths['output_top_dir'], 
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
        key_base = 'WickerFam_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
        ClusterSummaryFl = 'WickerClusterSummary_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
        ClusterMembershipFl = 'WickerClusterMembership_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
    elif clustering_method == 'MCL':
        MCLCLUST = True
        clustMethodTopDir = '{0}/MCL/I{1}'.format(paths['output_top_dir'], I)
        key_base = 'MCL_I{0}'.format(I)
        ClusterSummaryFl = 'MCL_ClusterSummary_I{0}'.format(I)
        ClusterMembershipFl = 'MCL_ClusterMembership_I{0}'.format(I)
    append2logfile(paths['output_top_dir'], mainlogfile,
                             'Beginning SoloLTRsearch(): {0}'.format(key_base))
    OutputDir = 'SoloLTRSearch.{0}'.format(key_base)
    paths[OutputDir] = '{0}/SoloLTRsearch'.format(clustMethodTopDir)
    MakeDir(OutputDir, paths[OutputDir])
    paths['SoloLTRsGFFsDir'] = '{0}/GFFs'.format(paths[OutputDir])
    MakeDir('SoloLTRsGFFsDir', paths['SoloLTRsGFFsDir'])

    LTRsGFF = '{0}.LTRs'.format(key_base)
    paths[LTRsGFF] = '{0}/{1}.gff'.format(paths[OutputDir], LTRsGFF)
    LTRsFASTA = '{0}.LTRs_FASTA'.format(key_base)
    paths[LTRsFASTA] = '{0}/{1}.fasta'.format(paths[OutputDir], LTRsFASTA)
    BLASToutput = '{0}.LTRs.blastn2ref'.format(key_base)
    paths[BLASToutput] = '{0}.blastn.tsv'.format(paths[LTRsFASTA])
    SoloLTRsGFF = '{0}.SoloLTRsGFF'.format(key_base)
    paths[SoloLTRsGFF] = '{0}/{1}.gff'.format(paths['SoloLTRsGFFsDir'], LTRsGFF)
    RepeatRegionsGFF = 'repeat_regions'
    paths['RepeatRegionsGFF'] = '{0}/{1}.gff'.format(paths[OutputDir],
                                                     RepeatRegionsGFF)

    SoloLTRsummary = '{0}.SoloLTRsummary'.format(key_base)
    paths[SoloLTRsummary] = '{0}/{1}.tsv'.format(paths[OutputDir],
                                                 SoloLTRsummary)
    if checkStatusFl(SoloLTRsummary) and checkStatusFl(SoloLTRsGFF):
        append2logfile(paths['output_top_dir'], mainlogfile, 
                  'SoloLTRsearch() already completed for {0}'.format(key_base))
        return
    append2logfile(paths['output_top_dir'], mainlogfile, 
                                'Writing GFFs for SoloLTRsearch(): {0}'.format(
                                                                     key_base))
    # Write repeat_region features for checking overlaps with LTR hits
    writeLTRretrotransposonGFF(paths['CurrentGFF'], 
                               paths['RepeatRegionsGFF'], 
                               elementSet=None, 
                               REPEATREGION=True, 
                               truncateParent=False) 


    lastEl = None
    LTRlengths = {}
    with open(paths[LTRsGFF], 'w') as outFl:
        # Write LTRs GFF
        with open(paths['CurrentGFF'], 'r') as inFl:
            for line in inFl:
                if line.startswith('#'):
                    continue
                gl = GFF3_line(line)
                if gl.type == 'long_terminal_repeat':
                    el = gl.attributes['Parent']
                    if el != lastEl:
                        gl.attributes['ID'] = '{0}.1'.format(el)
                    else:
                        gl.attributes['ID'] = '{0}.2'.format(el)
                    LTRlengths[gl.attributes['ID']] = gl.end - gl.start + 1
                    gl.attributes_order.insert(0, 'ID')
                    gl.refreshAttrStr()
                    outFl.write(str(gl)+'\n')
                    lastEl = el
    # Get LTRs FASTA
    getfasta(paths[LTRsGFF], paths['inputFasta'], paths[LTRsFASTA], 
                                                               headerKey='ID') 
    if not checkStatusFl(BLASToutput):
        append2logfile(paths['output_top_dir'], mainlogfile, 
                    'Running blastn for SoloLTRsearch(): {0}'.format(key_base))
        # Do blastn x input fasta
        runblast(query = paths[LTRsFASTA], 
                 subject = paths['inputFasta'], 
                 out = paths[BLASToutput], 
                 evalue = soloLTRmaxEvalue, 
                 outfmt = '7', 
                 percid = soloLTRminPid, 
                 blast = 'blastn', 
                 procs = procs)
        with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
            statusFlAppend.write('{0}\t{1}\n'.format(BLASToutput, 
                                                           paths[BLASToutput]))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                   'Finished blastn for SoloLTRsearch(): {0}'.format(key_base))
    # parse results and assign LTRs to clusters and write GFF3 for each cluster
    append2logfile(paths['output_top_dir'], mainlogfile, 
       'Parsing repeat_regions GFF3 for SoloLTRsearch(): {0}'.format(
                                                                     key_base))
    # repeat_region features that hits shouldn't overlap
    RR = {}
    i = 0
    with open(paths['RepeatRegionsGFF'], 'r') as inFl:
        for line in inFl:
            if not line.startswith('#'):
                gl = GFF3_line(line)
                if not gl.seqid in RR:
                    # convert coords from 0-based (GFF3) to 1-based for 
                    # easier overlap comparison
                    RR[gl.seqid] = {i:{'coords':(gl.start+1, gl.end)}}
                else:
                    RR[gl.seqid][i] = {'coords':(gl.start+1, gl.end)}
                i+=1
    append2logfile(paths['output_top_dir'], mainlogfile, 
                                 ('Parsing potential solo LTR hits from blastn'
                                 ' output for SoloLTRsearch(): {0}').format(
                                                                     key_base))
    # hits from genomic blast of LTRs
    Hits = {}
    i = 0
    with open(paths[BLASToutput], 'r') as inFl:
        for line in inFl:
            if not line.startswith('#'):
                c = line.strip().split()
                query, subj, pid, alnLen = c[:4]
                LTRlen = int(c[7]) - int(c[6]) + 1
                # percent of LTR length in the alignment
                pLen = LTRlen/LTRlengths[query]*100 
                # skip alignments shorter 
                if pLen < soloLTRminLen:
                    continue
                bit = float(c[-1])
                # subject start and end
                s, e = [int(x) for x in c[8:10]]
                if not subj in Hits:
                    # convert coords from 0-based (GFF3) to 1-based for easier
                    # overlap comparison
                    Hits[subj] = {i:{'coords':(s+1, e), 'bit':bit, 'pLen':pLen, 
                                                                  'LTR':query}}
                else:
                    Hits[subj][i] = {'coords':(s+1, e), 'bit':bit, 'pLen':pLen, 
                                                                   'LTR':query}
                i+=1
    append2logfile(paths['output_top_dir'], mainlogfile, 
                        ('SoloLTRsearch(): {0}\ndiscarding hits that overlap '
                         'repeat_regions').format(key_base))
    # Discard hits that overlap a repeat region
    HitsR1 = {}
    for scaf in Hits:
        if scaf not in RR:
            HitsR1[scaf] = Hits[scaf]
            continue
        KEEP = True
        for i in Hits[scaf]:
            coordsLTR = Hits[scaf][i]['coords']
            LTR = Hits[scaf][i]['LTR']
            for j in RR[scaf]:
                if not KEEP:
                    break
                coordsLTRRT = RR[scaf][j]['coords']
                if Overlaps(coordsLTR, coordsLTRRT):
                    KEEP = False
            if KEEP:
                if scaf in HitsR1:
                    HitsR1[scaf][i] = Hits[scaf][i]
                else:
                    HitsR1[scaf] = {i:Hits[scaf][i]}

    # If hits overlap, keep only the hit with the highest bit score
    append2logfile(paths['output_top_dir'], mainlogfile, 
                        ('SoloLTRsearch(): {0}\nretaining highest scoring hits'
                         ' among hits that overlap').format(key_base))
    HitsR2 = {}
    for scaf in HitsR1:
        # Sort hits by bit score to allow not checking bit score if 
        # Overlaps()==True
        sortedHits = sorted([(H, sorted(list(HitsR1[scaf][H].items()), 
                                key=lambda x:x[0])) for H in HitsR1[scaf]], 
                                         key=lambda y:y[1][1][1], reverse=True)
        # Elements are discarded if they overlap or are added to HitsR2
        Remaining = [ LTRinfo[0] for LTRinfo in sortedHits ]
        while Remaining != []:
            LTR1 = Remaining[0]
            Remove = set()
            for LTR2 in Remaining:
                # skip identities
                if LTR1 == LTR2:
                    continue
                if Overlaps(HitsR1[scaf][LTR1]['coords'], 
                                         HitsR1[scaf][LTR2]['coords']):
                    # discard
                    Remove.add(LTR2)
            if scaf in HitsR2:
                HitsR2[scaf][LTR1] = HitsR1[scaf][LTR1]
            else:
                HitsR2[scaf] = {LTR1:HitsR1[scaf][LTR1]}
            # kept, but remove from Remaining
            Remove.add(LTR1)
            Remaining = [LTR for LTR in Remaining if not LTR in Remove]
            if Remaining == []:
                break
    # write table with cluster and # of solo LTRs
    append2logfile(paths['output_top_dir'], mainlogfile, 
                        ('SoloLTRsearch(): {0}\nparsing cluster membership and'
                         ' summary files').format(key_base))
    ClusterMembership = {}
    ClusterSizes = {}
    with open(paths[ClusterMembershipFl], 'r') as inFl:
        for line in inFl:
            if not line.startswith('element'):
                el, classif, clust = line.strip().split()
                ClusterMembership[el] = (classif, clust)
                if classif in ClusterSizes:
                    if clust in ClusterSizes[classif]:
                        ClusterSizes[classif][clust] += 1
                    else:
                        ClusterSizes[classif][clust] = 1
                else:
                    ClusterSizes[classif] = {clust:1}
    SoloLTRclusterMembership = {}
    GFFoutput = {}
    for scaf in HitsR2:
        GFFoutput[scaf] = {}
        for LTRinfo in HitsR2[scaf]:
            LTRname = HitsR2[scaf][LTRinfo]['LTR']
            el = LTRname.split('.')[0].lstrip('LTR_retrotransposon')
            try:
                classif, clust = ClusterMembership[el]
            except KeyError:
                print('{0} not in ClusterMembership file'.format(el), 
                                                               file=sys.stderr)
            if classif not in GFFoutput[scaf]:
                GFFoutput[scaf][classif] = {clust:{LTRname:HitsR2[scaf][LTRinfo]}}
            else:
                if clust not in GFFoutput[scaf][classif]:
                    GFFoutput[scaf][classif][clust] = {LTRname:HitsR2[scaf][LTRinfo]}
                else:
                    GFFoutput[scaf][classif][clust][LTRname] = HitsR2[scaf][LTRinfo]
            if classif in SoloLTRclusterMembership:
                if clust in SoloLTRclusterMembership[classif]:
                    SoloLTRclusterMembership[classif][clust] += 1
                else:
                    SoloLTRclusterMembership[classif][clust] = 1
            else:
                SoloLTRclusterMembership[classif] = {clust:1}
    if not checkStatusFl(SoloLTRsummary):
        append2logfile(paths['output_top_dir'], mainlogfile, 
                    'SoloLTRsearch(): {0}\nwriting summary output'.format(
                                                                     key_base))
        with open(paths[SoloLTRsummary], 'w') as outFl:
            outFl.write(('classif\tclust\tFullLengthElements\tNumberOfSoloLTRs'
                                                          '\tSolo2FullRatio\n'))
            for classif in sorted(list(SoloLTRclusterMembership.keys())):
                for clust in sorted(list(SoloLTRclusterMembership[classif]), 
                                                                      key=int):
                    solos = SoloLTRclusterMembership[classif][clust]
                    fulls = ClusterSizes[classif][clust]
                    ratio = solos/fulls
                    outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(classif, 
                                                                   clust, 
                                                                   fulls, 
                                                                   solos, 
                                                                   ratio))
        with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
            statusFlAppend.write('{0}\t{1}\n'.format(SoloLTRsummary, paths[SoloLTRsummary]))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                        'SoloLTRsearch(): {0}\nSummary file written'.format(
                                                                     key_base))
    if not checkStatusFl(SoloLTRsGFF):
        append2logfile(paths['output_top_dir'], mainlogfile, 
                        'SoloLTRsearch(): {0}\nWriting GFF3 output'.format(
                                                                     key_base))
        # remove any existing GFFs soas not to double-write when 
        # appending with write() below
        for scaf in GFFoutput:
            for classif in GFFoutput[scaf]:
                if os.path.isfile('{0}/{1}_{2}.SoloLTRs.gff'.format(
                                                    paths['SoloLTRsGFFsDir'], 
                                                    key_base, 
                                                    classif)):
                    os.remove('{0}/{1}_{2}.SoloLTRs.gff'.format(
                                                    paths['SoloLTRsGFFsDir'], 
                                                    key_base, 
                                                    classif))
                clustersOut = 'SoloLTRs{0}'.format(classif)
                if classif in GFFoutput[scaf]:
                    for clust in GFFoutput[scaf][classif]:
                        if clustersOut in paths:
                            if os.path.isfile('{0}/{1}_{2}_cluster_{3}.SoloLTRs.gff'.format(
                                                           paths[clustersOut], 
                                                           key_base, 
                                                           classif, 
                                                           clust)):
                                os.remove('{0}/{1}_{2}_cluster_{3}.SoloLTRs.gff'.format(
                                                           paths[clustersOut], 
                                                           key_base, 
                                                           classif, 
                                                           clust))
        # write GFF files with coordinates and store info for summary file
        with open(paths[SoloLTRsGFF], 'w') as outFl:
            i=0
            outFl.write('##gff-version 3\n')
            for scaf in GFFoutput:
                append2logfile(paths['output_top_dir'], mainlogfile, 
                        'SoloLTRsearch(): {0}\nwriting GFF3 output to\n{1}'.format(
                                                           key_base, 
                                                           paths[SoloLTRsGFF]))
                for classif in GFFoutput[scaf]:
                    with open('{0}/{1}_{2}.SoloLTRs.gff'.format(
                                                paths['SoloLTRsGFFsDir'], 
                                                key_base, 
                                                classif), 'a') as outClassifFl:
                        outClassifFl.write('##gff-version 3\n')
                        clustersOut = 'SoloLTRs{0}'.format(classif)
                        MakeDir(clustersOut, '{0}/{1}_clusters'.format(
                                            paths['SoloLTRsGFFsDir'], classif))
                        for clust in GFFoutput[scaf][classif]:
                            with open('{0}/{1}_{2}_cluster_{3}.SoloLTRs.gff'.format(
                                                  paths[clustersOut], 
                                                  key_base, 
                                                  classif, 
                                                  clust), 'a') as outClusterFl:
                                outClusterFl.write('##gff-version 3\n')
                                for relatedLTR in GFFoutput[scaf][classif][clust]:
                                    i+=1
                                    Info = GFFoutput[scaf][classif][clust][relatedLTR]
                                    closestLivingRelative = Info['LTR'].split(
                                         '.')[0] + '-LTR-' + Info['LTR'].split(
                                                                        '.')[1]
                                    start, end = sorted([int(c) for 
                                                          c in Info['coords']])
                                    score = Info['bit']
                                    outFl.write(('{0}\tPhyLTR\tSoloLTR\t{1}\t{2}'
                                                 '\t{3}\t?\t.\tID=LTR-{4}.{5}.'
                                                 'cluster_{6};relative={7}\n').format(
                                                        scaf, 
                                                        start, 
                                                        end, 
                                                        score, 
                                                        i, 
                                                        classif, 
                                                        clust, 
                                                        closestLivingRelative))
                                    outClassifFl.write(('{0}\tPhyLTR\tSoloLTR'
                                           '\t{1}\t{2}\t{3}\t?\t.\tID=LTR-{4}.'
                                           '{5}.cluster_{6};relative={7}\n').format(
                                                        scaf, 
                                                        start, 
                                                        end, 
                                                        score, 
                                                        i, 
                                                        classif, 
                                                        clust, 
                                                        closestLivingRelative))
                                    outClusterFl.write(('{0}\tPhyLTR\tSoloLTR'
                                           '\t{1}\t{2}\t{3}\t?\t.\tID=LTR-{4}.'
                                           '{5}.cluster_{6};relative={7}\n').format(
                                                        scaf, 
                                                        start, 
                                                        end, 
                                                        score, 
                                                        i, 
                                                        classif, 
                                                        clust, 
                                                        closestLivingRelative))
        with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
            statusFlAppend.write('{0}\t{1}\n'.format(SoloLTRsGFF, paths[SoloLTRsGFF]))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                         'SoloLTRsearch(): {0}\nGFFs written'.format(key_base))


def geneconvLTRs(g='/g0', 
                 I=6, 
                 clustering_method='WickerFam', 
                 WickerParams={'pId':80,'percAln':80,'minLen':80}):
    """Assesses intra-element gene conversion between LTRs
    g can be one of /g0, /g1, or /g2
    g is proportional to the tolerance for mismatches in fragments by 
    geneconv
    """
    global paths
    global filenames
    if clustering_method == 'WickerFam':
        WICKERCLUST = True
        key_base = 'WickerFam_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        WickerDir = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        AlnKey = '{0}.LTRAlnDir'.format(key_base)
        if not checkStatusFl(AlnKey):
            sys.exit('LTR alignment not in status file: {0}'.format(AlnKey))
        GENECONVTopDirKey = '{0}.GENECONV'.format(key_base)
        paths[GENECONVTopDirKey] = '{0}/GENECONV'.format(paths[WickerDir])
        MakeDir(GENECONVTopDirKey, paths[GENECONVTopDirKey])
        GENECONVDirKey = '{0}.GENECONVLTRs'.format(key_base)
        paths[GENECONVDirKey] = '{0}/LTRs'.format(paths[GENECONVTopDirKey])
        MakeDir(GENECONVDirKey, paths[GENECONVDirKey])
        GENECONVgDirKey = 'GENECONVgDir'
        paths[GENECONVgDirKey] = '{0}/{1}'.format(paths[GENECONVDirKey], g[1:])
        MakeDir(GENECONVgDirKey, paths[GENECONVgDirKey])
        SummaryKey = '{0}.GENECONVLTRs.Summary'.format(key_base)
    elif clustering_method == 'MCL':
        MCLCLUST = True
        key_base = 'MCL_I{0}'.format(I)
        MCLdir = 'MCL_I{0}'.format(I)
        AlnKey = '{0}.LTRAlnDir'.format(key_base)
        if not checkStatusFl(AlnKey):
            sys.exit('LTR alignment not in status file: {0}'.format(AlnKey))
        GENECONVTopDirKey = '{0}.GENECONV'.format(key_base)
        paths[GENECONVTopDirKey] = '{0}/GENECONV'.format(paths[MCLdir])
        MakeDir(GENECONVTopDirKey, paths[GENECONVTopDirKey])

        GENECONVDirKey = '{0}.GENECONVLTRs'.format(key_base)
        paths[GENECONVDirKey] = '{0}/LTRs'.format(paths[GENECONVTopDirKey])
        MakeDir(GENECONVDirKey, paths[GENECONVDirKey])

        GENECONVgDirKey = 'GENECONVgDir'
        paths[GENECONVgDirKey] = '{0}/{1}'.format(paths[GENECONVDirKey], g[1:])
        MakeDir(GENECONVgDirKey, paths[GENECONVgDirKey])

        SummaryKey = '{0}.GENECONVLTRs.Summary'.format(key_base)

    else:
        sys.exit(('modeltest() parameter clustering_method needs to be either '
                  'WickerFam or MCL and it is: {0}').format(clustering_method))

    if checkStatusFl(SummaryKey):
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
                        'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                   scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                            'ltr_divergence() already completed: {0}'.format(
                                   paths['{0}.GENECONVLTRs'.format(key_base)]))
        return
    if not checkStatusFl('{0}.GENECONVLTRs.{1}'.format(key_base, g[1:])):
        geneconv_calls = []
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
                        'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                   scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile,
               'Checking directory structure for GENECONV using {0}'.format(g))
        files = [f for f in os.listdir(paths[AlnKey]) if f.endswith('trimmed')]
        num_elements = len(files)
        alnLens = {}
        show = False
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
                        'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                   scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, ('Preparing to run'
             ' GENECONV for finding intraelement gene conversion between LTRs'))
        for f in files:
            flpth = '{0}/{1}'.format(paths[AlnKey], f)
            call = [executables['geneconv'], flpth, '/w123', g, 
                                   '-include_monosites', '-nolog', '-Dumptab']
            geneconv_calls.append((call, '/dev/null', None, None))
            elementName = '_'.join(f.split('_')[:2])
            seqs = list(SeqIO.parse(flpth, 'fasta'))
            try:
                alnLen = len(seqs[0].seq)
            except IndexError:
                continue
            elementName = seqs[0].id
            alnLens[elementName] = alnLen
            show=False
        # Make calls for each ltr pair
        if not checkStatusFl('{0}.GENECONVLTRs.{1}'.format(key_base, g[1:])):
            chunk_size = ceil(len(files)/procs)
            with Pool(processes = procs) as p:
                p.map(makecallMultiprocessing, geneconv_calls, 
                                                        chunksize = chunk_size)
            p.join()
            output = [f for f in os.listdir(paths[AlnKey]) if f.endswith('tab')]
            # Move geneconv files to geneconv dir
            for f in output:
                os.rename('{0}/{1}'.format(paths[AlnKey], f), 
                                   '{0}/{1}'.format(paths[GENECONVgDirKey], f))
            paths['{0}.GENECONVLTRs.{1}'.format(
                                     key_base, g[1:])] = paths[GENECONVgDirKey]
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}.GENECONVLTRs.{1}\t{2}\n'.format(
                                                    key_base, g[1:], 
                                                    paths[GENECONVgDirKey])) 
        append2logfile(paths['output_top_dir'], mainlogfile, 
                                                     'Parsing GENECONV output')
        # Parse geneconv files
        # sig is short for "significant". To contain significant 
        # global inner fragments identified by GENECONV
        sig = [] 
        for f in os.listdir(paths[GENECONVgDirKey]):
            if f.endswith('tab'):
                with open('{0}/{1}'.format(paths[GENECONVgDirKey], f)) as fl:
                    sig += [line for line in fl.read().split('\n') if 
                                                         line.startswith('GI')]
        # Sort by element name, which at sig[i][1] are as: 
        # LTR_retrotransposon1;LTR_retrotransposon1
        sig = sorted(sig, key=lambda x:x[1])
        paths['GENECONVsummary'] = '{0}/GENECONV_summary'.format(
                                                         paths[GENECONVDirKey])
        paths['GENECONV_output'] = '{0}/GENECONVoutput_tab'.format(
                                                         paths[GENECONVDirKey])
        IAGCpositive = set()
        with open(paths['GENECONV_output'], 'w') as outputFl:
            with open(paths['GENECONVsummary'], 'w') as summaryFl:
                summaryFl.write(('# ratio is the ratio of the alignment length '
                   'to the alignment length minus the gene conversion tract\n'))
                try:
                    summaryFl.write(('# {0} elements out of {1}, or {2:.1f}% '
                         'with possible evidence of gene conversion\n').format(
                                                 len(sig), 
                                                 num_elements, 
                                                 ((len(sig)/num_elements)*100)))
                except ZeroDivisionError:
                    summaryFl.write(('# {0} elements out of {1}, or 0% with '
                               'possible evidence of gene conversion\n').format(
                                                       len(sig), num_elements))
                summaryFl.write('# 5% are expected by chance\n')
                summaryFl.write(('#elementName\tsim_p-val\talnLen\tstart\tend'
                          '\ttractLen\tratio_alnLen2alnLen-tractLen\tgScale\n'))
                for line in sig:
                    totDifs = int(line.strip().split()[9])
                    if totDifs < 3:
                        continue
                    outputFl.write(line + '\n')
                    contents = line.split('\t')
                    element = contents[1].split(';')[0]
                    IAGCpositive.add(element)
                    sim_p = float(contents[2])
                    #KA_p = float(contents[3])
                    start=int(contents[4])
                    end=int(contents[5])
                    tractLen = int(contents[6])
                    alnLen = alnLens[element]
                    ratio = alnLen / (alnLen-tractLen)
                    summaryFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(
                                                            element, 
                                                            sim_p, 
                                                            alnLen, 
                                                            start, 
                                                            end, 
                                                            tractLen, 
                                                            ratio, 
                                                            g[1:]))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                'GENECONV output and a summary written to:\n{0}\n{1}'.format(
                                                    paths['GENECONV_output'], 
                                                    paths['GENECONVsummary']))
        paths[SummaryKey] = paths['GENECONVsummary']
        with open('{0}/status'.format(paths['output_top_dir']), 'a') as statusFlAppend:
            statusFlAppend.write('{0}\t{1}\n'.format(SummaryKey,
                                                     paths['GENECONVsummary']))


def ltr_divergence(I=6, 
                   clustering_method='WickerFam', 
                   WickerParams={'pId':80,'percAln':80,'minLen':80}, 
                   iters=1, 
                   model='hky85'):
    """
    Runs PAUP
    iters used with modeltest(). deprecated pretty much
    """
    global paths
    WICKERCLUST = False
    MCLCLUST = False
    if clustering_method == 'WickerFam':
        key_base = 'WickerFam_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        WICKERCLUST = True
        WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])]
        WickerLTRdivDirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_LTRdivDir'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        statusFlKey = '{0}_PAUP_divergence_dir'.format(WickerLTRdivDirkey)
        paths[WickerLTRdivDirkey] = '{0}/LTR_divergence'.format(WickerDir)
        paths['DivergenceTopDir'] = paths[WickerLTRdivDirkey]
        SummaryKey = 'WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_LTR_divergence_summary'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
    elif clustering_method == 'MCL':
        key_base = 'MCL_I{0}'.format(I)
        MCLCLUST = True
        MCLdir = paths['MCL_I{0}'.format(I)]
        MCL_LTRdivDirkey = 'MCL_I{0}_LTRdivDir'.format(I)
        statusFlKey = '{0}_PAUP_divergence_dir'.format(MCL_LTRdivDirkey)
        paths[MCL_LTRdivDirkey] = '{0}/LTR_divergence'.format(MCLdir)
        paths['DivergenceTopDir'] = paths[MCL_LTRdivDirkey]
        SummaryKey = 'MCL_I{0}_LTR_divergence_summary'.format(I)
        TrimalKey = '{0}.LTRs.aln.trimal'.format(key_base)
    else:
        sys.exit(('modeltest() parameter clustering_method needs to be either '
                    'WickerFam or MCL and it is: {0}').format(clustering_method))
    MakeDir('PAUPdivergenceDir', '{0}/PAUP'.format(paths['DivergenceTopDir']))
    MakeDir('PAUPNexusInputDir', '{0}/nexus'.format(
                                         paths['PAUPdivergenceDir'.format(I)]))
    MakeDir('PAUPDivOutDir', '{0}/divergences'.format(
                                         paths['PAUPdivergenceDir'.format(I)]))
    if checkStatusFl(SummaryKey):
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
          'Below log entry is from line {0} in {1}'.format(lineno, scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, 
           'ltr_divergence() already completed: {0}'.format(paths[SummaryKey]))
        return
    if LTRDIVERGENCE:
        paupCalls = []
        modeltestResults = {}
        if not checkStatusFl(statusFlKey):
            for classif in classifs:
                # Align only clusters from superfamilies with 
                # identical LTRs on transposition (Copia, Gypsy, ERV, BEL/Pao)
                YESLTRS = False
                for SF in LTR_SFs:
                    if classif.startswith(SF):
                        YESLTRS = True
                        break
                if not YESLTRS:
                    continue
                # parse model test results for paup block and summary 
                # for best model
                if MCLCLUST:
                    ModeltestSummaryKey = 'MCLModelTestSummary_{0}_iters_I{1}_{2}'.format(
                                                              iters,I, classif)
                elif WICKERCLUST:
                    ModeltestSummaryKey = 'WickerModelTestSummary_{0}_iters_{1}_pId_{2}_percAln_{3}_minLen_{4}'.format(
                                                        iters, 
                                                        WickerParams['pId'],   
                                                        WickerParams['percAln'],   
                                                        WickerParams['minLen'],    
                                                        classif)
                modeltestResults[classif] = {}
                # Create PAUP calls for each cluster
                if MCLCLUST:
                    clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
                elif WICKERCLUST:
                    clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                                        WickerParams['pId'],   
                                                        WickerParams['percAln'],   
                                                        WickerParams['minLen'],    
                                                        classif)]
                clusters = [clust.split('\t') for clust in open(
                                  clusterPath,'r').read().strip().split('\n')]
                PhyLTRalnPths = {'_'.join(f.split('_')[:2]):'{0}/{1}'.format(
                                   paths[TrimalKey], f) for f in os.listdir(
                                    paths[TrimalKey]) if f.endswith('trimmed')} 
                if ModeltestSummaryKey in paths:
                    with open(paths[ModeltestSummaryKey]) as testOutputFl:
                        paupLines = ''
                        method = None
                        model = None
                        clust = None
                        for line in testOutputFl:
                            if line.startswith('BIC'):
                                if not method == None:
                                    modeltestResults[classif][clust] = (
                                                              paupLines, model)
                                if len(line.strip().split()) == 17:
                                    method, model, pa, pc, pg, pt, kappa, titv, 
                                    Ra, Rb, Rc, Rd, Re, Rf, pInv, gamma, 
                                    clust = line.strip().split()
                                else:
                                    method, model, pa, pc, pg, pt, kappa, titv, 
                                    Ra, Rb, Rc, Rd, Re, Rf, pInv, gamma, clust, 
                                    clustSize = line.strip().split()
                                paupLines = ''
                            else:
                                paupLines += line
                        if not method == None:
                            modeltestResults[classif][clust] = (paupLines, model)
                    for j in range(len(clusters)):
                        for el in clusters[j]:
                            if el == '':
                                continue
                            paupBlock = ''
                            seqRec = list(SeqIO.parse(PhyLTRalnPths[el], 'fasta'))
                            paupBlock = """#NEXUS
begin DATA;
DIMENSIONS ntax=2 nchar={0};
FORMAT datatype=dna missing=? gap=-;
MATRIX
'{1}' {3}
'{2}' {4}
;
END;
""".format(len(seqRec[0].seq), seqRec[0].id+'_L', seqRec[1].id+'_R', str(
                                            seqRec[0].seq), str(seqRec[1].seq))

                            # Have modeltest result for this cluster
                            if str(j) in modeltestResults[classif]:
                                model = modeltestResults[classif][str(j)][1]
                                for line in modeltestResults[classif][str(j)][0].split('\n'):
                                    if line.startswith('SaveDist'):
                                        paupBlock += ('SaveDist format=oneColumn '
                                               'file=../divergences/divergence.'
                                               '{0}.{1};\n').format(model, el)
                                    else:
                                        paupBlock += line + '\n'
                            # No model test result, perhaps the cluster 
                            # is too small
                            else:
                                # HKY+85 does not work when there are 
                                # only 2 character states
                                if len(set(str(seqRec[0].seq))) == 2 and len(
                                                set(str(seqRec[1].seq))) == 2:
                                    model = 'JC'
                                    append2logfile(paths['output_top_dir'], 
                                        mainlogfile, 
                                        ('Model testing not done for {0}\tcluster'
                                        ':{1}. Using JC').format(classif, j))
                                    paupBlock += """[!
Likelihood settings from best-fit model (JC) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold and there are only 2 character
states.]

BEGIN PAUP;
Dset distance=jc;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
""".format(model, el)
                                else:
                                    model = 'HKY85'
                                    scriptpath = os.path.realpath(__file__)
                                    lineno = getframeinfo(currentframe()).lineno + 2
                                    append2logfile(paths['output_top_dir'], 
                                        mainlogfile, 
                                        'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
                                    append2logfile(paths['output_top_dir'], 
                                        mainlogfile, 
                                        'Model testing not done for {0}\tcluster:{1}. Using HKY85'.format(
                                                                   classif, j))
                                    paupBlock += """[!
Likelihood settings from best-fit model (HKY) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold.]

BEGIN PAUP;
Dset distance=hky85;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
""".format(model, el)
                            with open('{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(
                                     paths['PAUPNexusInputDir'], seqRec[0].id, 
                                           model, classif, j), 'w') as nexusFl:
                                nexusFl.write(paupBlock)
                            paup_call = [executables['paup'], '-n', 
                                '{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(
                                                   paths['PAUPNexusInputDir'], 
                                                   seqRec[0].id, 
                                                   model, 
                                                   classif, 
                                                   j)]
                            packet = (paup_call, None, None, None)
                            paupCalls.append(packet)
                else:
                    for j in range(len(clusters)):
                        for el in clusters[j]:
                            if el == '':
                                continue
                            paupBlock = ''
                            seqRec = list(SeqIO.parse(PhyLTRalnPths[el], 'fasta'))
                            paupBlock = """#NEXUS
begin DATA;
DIMENSIONS ntax=2 nchar={0};
FORMAT datatype=dna missing=? gap=-;
MATRIX
'{1}' {3}
'{2}' {4}
;
END;
""".format(len(seqRec[0].seq), seqRec[0].id+'_L', seqRec[1].id+'_R', str(
                                            seqRec[0].seq), str(seqRec[1].seq))
                             # HKY+85 does not work when there are only 
                             # 2 character states
                            if (len(set(str(seqRec[0].seq))) == 2 and 
                                        len(set(str(seqRec[1].seq))) == 2):
                                append2logfile(paths['output_top_dir'], 
                                    mainlogfile, 
                                    'Model testing not done for {0}\tcluster:{1}. Using JC'.format(
                                                                   classif, j))
                                model = 'JC'
                                paupBlock += """[!
Likelihood settings from best-fit model (JC) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold and there are only 2 character
states.]

BEGIN PAUP;
Dset distance=jc;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
""".format(model, el)
                            else:
                                append2logfile(paths['output_top_dir'], 
                                    mainlogfile, 
                                    'Model testing not done for {0}\tcluster:{1}. Using HKY85'.format(
                                                                   classif, j))
                                model = 'HKY85'
                                paupBlock += """[!
Likelihood settings from best-fit model (HKY) selected by default
Model testing not done for some reason, possibly the cluster size is
less than the user specified threshold.]

BEGIN PAUP;
Dset distance=hky85;
SaveDist format=oneColumn file=../divergences/divergence.{0}.{1};
END;
""".format(model, el)
                            with open('{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(
                                            paths['PAUPNexusInputDir'], 
                                            seqRec[0].id, 
                                            model, 
                                            classif, 
                                            j), 'w') as nexusFl:
                                nexusFl.write(paupBlock)
                            paup_call = [executables['paup'], '-n', 
                                '{0}/{3}_{4}_{1}_divergence_{2}.nex'.format(
                                            paths['PAUPNexusInputDir'], 
                                            seqRec[0].id, 
                                            model, 
                                            classif, 
                                            j)]
                            packet = (paup_call, None, None, None)
                            paupCalls.append(packet)
            chunk_size = ceil(len(paupCalls)/procs)
            with Pool(processes=procs) as p:
                p.map(makecallMultiprocessing, paupCalls, chunksize=chunk_size)
            p.join()
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format(statusFlKey, 
                                                   paths['PAUPdivergenceDir']))
    gcDct = {}
    if GENECONVLTRS:
        GENECONVSummaryKey = '{0}.GENECONVLTRs.Summary'.format(key_base)
        with open(paths[GENECONVSummaryKey], 'r') as gcFl:
            for line in gcFl:
                if not line.startswith('#'):
                    (el, p, alnLen, start, end, tractLen, ratio, 
                                                  gscale) = line.strip().split()
                    ratio = float(ratio)
                    alnLen = int(alnLen)
                    start = int(start)
                    end = int(end)

                    if not el in gcDct:
                        gcDct[el] = [(start, end, alnLen, gscale)]
                    else:
                        gcDct[el].append((start, end, alnLen, gscale))

        # Average ratio for those elements with multiple predicted GC tracts.
        # Currently all LTR divergence estimates are scaled as if all 
        # GC tracts were g0 (no mismatches)
        for el in gcDct:
            # If True there is only one GC tract for this element.
            if len(gcDct[el]) == 1: 
                tractLen = gcDct[el][0][1] - gcDct[el][0][0]
                alnLen = gcDct[el][0][2]
                gcDct[el] = alnLen / (alnLen - tractLen)
            else:
                alnLen = gcDct[el][0][2]
                coords = []
                for tract in gcDct[el]:
                    coords.append(tract[:2])
                coords = sorted(coords, key=lambda x:x[0])
                c1 = coords[0]
                tracts = []
                for c2 in coords[1:]:
                    c3 = mergeCoords(c1, c2) 
                    if c3[0][0] == c3[0][1]:
                        c1 = c3[0][0]
                    else:
                        tracts.append(c3[0][0])
                        c1 = c3[0][1]
                tracts.append(c1)
                tractLen = 0
                for t in tracts:
                    tractLen += t[1] - t[0] + 1
                gcDct[el] = alnLen / (alnLen - tractLen)
    # Read PAUP output
    clustLenDcts = {}
    clustDct = {}
    for classif in classifs:
        # Align only clusters from superfamilies with identical LTRs 
        # on transposition (Copia, Gypsy, ERV, BEL/Pao)
        YESLTRS = False
        for SF in LTR_SFs:
            if classif.startswith(SF):
                YESLTRS = True
                break
        if not YESLTRS:
            continue
        clustDct[classif] = {}
        if MCLCLUST:
            clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
        elif WICKERCLUST:
            clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                                WickerParams['pId'], 
                                                WickerParams['percAln'], 
                                                WickerParams['minLen'], 
                                                classif)]
        clusters = [clust.split('\t') for clust in open(
                                   clusterPath,'r').read().strip().split('\n')]
        for j in range(len(clusters)):
            for el in clusters[j]:
                clustDct[classif][el] = j

        clustLenDct = {i:len(clusters[i]) for i in range(len(clusters))}
        clustLenDcts[classif] = clustLenDct

    if MCLCLUST:
        paths['DivergenceSummary'] = '{0}/MCL_I{1}.LTR_divergences.tab'.format(
                                                  paths['DivergenceTopDir'], I)
    elif WICKERCLUST:
        paths['DivergenceSummary'] = '{0}/WickerFamDir_{1}_pId_{2}_percAln_{3}_minLen.LTR_divergences.tab'.format(
                                                paths['DivergenceTopDir'], 
                                                WickerParams['pId'], 
                                                WickerParams['percAln'], 
                                                WickerParams['minLen'])
    with open(paths['DivergenceSummary'.format(I)],'w') as divSummaryFl:
        divSummaryFl.write(('elementName\tclassification\tMCLinflationValue\t'
                            'cluster\tclusterSize\tmodel\tdivergence\t'
                            'correctedDivergence\tIntraelementGeneConversion\n'))
    for fname in os.listdir(paths['PAUPDivOutDir']):
        fPth = '{0}/{1}'.format(paths['PAUPDivOutDir'], fname)
        el1, el2, div = open(fPth, 'r').read().strip().split('\t')
        el = el1[:-2]
        div = float(div)
        model = fname.split('.')[1]
        if list(clustDct.keys()) == ['All']:
            classif == 'All'
        else:
            classif = classifs_by_element[el]
        if classif in clustDct:
            if el in clustDct[classif]:
                clust = clustDct[classif][el]
            else:
                continue
        else:
            continue
        clustSize = clustLenDcts[classif][clust]
        if el in gcDct:
            divc = div*gcDct[el]
            GC = 'Yes'
        else:
            divc = div
            GC = 'No'
        if div > 5 or divc > 5:
            print(('({0} LTRs have estimated substitutions per site of {1} and '
                  '{2} (gene conversion-corrected) using {3}. It was excluded '
                  'from the summary table at: {4}').format(el, 
                                                           div, 
                                                           divc, 
                                                           model, 
                                                           paths['DivergenceSummary']), 
                                                       file=sys.stderr)
            continue
        # Combine as new output
        with open(paths['DivergenceSummary'],'a') as divSummaryFl:
            divSummaryFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(
                       el, classif, I, clust, clustSize, model, div, divc, GC))
    if not SummaryKey in paths:
        with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
            statusFlAppend.write('{0}\t{1}\n'.format(SummaryKey, 
                                                   paths['DivergenceSummary']))


def phylo(removegeneconv = True, 
          BOOTSTRAP = True, 
          I = 6, 
          align = 'cluster', 
          removehomologouspair = True, 
          part = 'entire', 
          clustering_method = 'WickerFam', 
          WickerParams = {'pId':80,'percAln':80,'minLen':80}, 
          auto_outgroup = False, 
          bootstrap_reps = 100, 
          minClustSize = 4, 
          convert_to_ultrametric = False, 
          bpflank = None, 
          combine_and_do_small_clusters = True):
    """
    align              one of cluster or classifs
    removegeneconv     If True GENECONV output will be parsed and 
                       elements with evidence of gene conversion will 
                       be removed from cluster prior to aligning for 
                       tree
    BOOTSTRAP=False    is not implemented yet.
    I                  inflation parameter used for MCL clusters desired 
                       to use for this function
    align              'cluster' align and make trees for each WickerFam 
                       or MCL cluster, depending on setting of 
                       clustering_method. 'classif' align and make trees 
                       for each superfamily
    removehomologouspair    Not implementd for align='cluster' yet.
                            If True then the flanking X bp on both sides 
                            of the element are extracted and blastn'd 
                            together. (constraints can be changed, see). 
                            One of each pair of hits is removed prior to 
                            alignment.
    part                'entire' Use whole sequence of element for 
                        trees. Not including TSD. From: 
                        start(LTR1) -> end(LTR2)
                        'inernal' Use sequence between LTRs for trees. 
                        From: end(LTR1) -> start(LTR2)
    clustering_method   'WickerFam' or 'MCL'
    WickerParams        A dictionary of the format 
                        {'pId':A,'percAln':B,'minLen':C} where A, B, C 
                        are numbers specifying the WickerFam() params.
    auto_outgroup       Not available for align='classif'.
                        Automatically roots each cluster's tree with an 
                        outgroup by this process: the outgroup shall be 
                        a random element from cluster k where cluster k 
                        is the largest of the clusters that is not j if 
                        j is the first cluster then the next smallest 
                        cluster is k if there is no other cluster, no 
                        outgroup is used
    """
    global paths
    if not PHYLO:
        return
    append2logfile(paths['output_top_dir'], mainlogfile, 
                     'Start phylo() for {0}'.format(clustering_method))
    AUTO_OUTGROUP = auto_outgroup
    REMOVEHOMOLOGOUSFLANK = removehomologouspair
    REMOVEGENECONV = removegeneconv
    ULTRAMETRIC = convert_to_ultrametric
    if AUTO_OUTGROUP:
        strOutgroup = 'withOutgroup'
    else:
        strOutgroup = 'noOutgroup'
    if REMOVEHOMOLOGOUSFLANK:
        strHomoflank = 'homoflank'
    else:
        strHomoflank = 'nohomoflank'
    WICKERCLUST = False
    MCLCLUST = False
    if clustering_method == 'WickerFam':
        WICKERCLUST = True
        WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])]
        WickerTreesDirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_TreesDir'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
        paths[WickerTreesDirkey] = '{0}/Trees'.format(WickerDir)
        paths['TreesDir'] = paths[WickerTreesDirkey]
        MakeDir(WickerTreesDirkey, paths['TreesDir'])

    elif clustering_method == 'MCL':
        MCLCLUST = True
        MCLdir = paths['MCL_I{0}'.format(I)]
        MCLTreesDirkey = 'MCL_I{0}_TreesDir'.format(I)
        paths[MCLTreesDirkey] = '{0}/Trees'.format(MCLdir)
        paths['TreesDir'] = paths[MCLTreesDirkey]
        MakeDir(MCLTreesDirkey, paths['TreesDir'])
    else:
        sys.exit(('modeltest() parameter clustering_method needs to be either '
                  'WickerFam or MCL, and it is: {0}').format(clustering_method))
    OutPth = None
    alnPthKeys = []
    ENTIRE = False
    INTERNAL = False
    if part == 'entire':
        ENTIRE = True
    elif part == 'internal':
        INTERNAL = True
    else:
        sys.exit(('part parameter to AutoAlign() is invalid: {0}. Must be '
                  'either "entire" or "internal"').format(part))
    if INTERNAL:
        MakeDir('InternalRegionsAlignments', '{0}/InternalRegions'.format(
                                                            paths['TreesDir']))
        paths['RegionDir'] = paths['InternalRegionsAlignments']
    elif ENTIRE:
        MakeDir('WholeElementAlignments', '{0}/WholeElements'.format(
                                                            paths['TreesDir']))
        paths['RegionDir'] = paths['WholeElementAlignments']
    # OUTGROUP is not possible here
    if align == 'classif':
        MakeDir('Superfamilies', '{0}/Superfamilies'.format(paths['RegionDir']))
        if REMOVEHOMOLOGOUSFLANK:
            strHomoflank = 'homoflank'
            MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(
                                                       paths['Superfamilies']))
        else:
            strHomoflank = 'nohomoflank'
            MakeDir('HomoFlankDir', '{0}/AllElements'.format(
                                                       paths['Superfamilies']))
        OutPth = paths['HomoFlankDir']
        if WICKERCLUST:
            AutoAlign(I = None, 
                      part = part, 
                      rmgeneconv = removegeneconv, 
                      minClustSize = minClustSize, 
                      align = 'classif', 
                      rmhomologflank = REMOVEHOMOLOGOUSFLANK, 
                      clustering_method = 'WickerFam', 
                      WickerParams = {'pId':WickerParams['pId'],
                        'percAln':WickerParams['percAln'],'minLen':['minLen']}, 
                      auto_outgroup = AUTO_OUTGROUP, 
                      bpflank = bpflank, 
                      combine_and_do_small_clusters = combine_and_do_small_clusters, 
                      flank_pId = flank_pId, 
                      flank_evalue = flank_evalue, 
                      flank_plencutoff = flank_plencutoff, 
                      LTRSONLY = False)
            alnPthKeys.append('WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_all.{4}.{5}'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'],
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        strHomoflank))
        elif MCLCLUST:
            AutoAlign(I = I, 
                      part = part, 
                      rmgeneconv = removegeneconv, 
                      minClustSize = minClustSize, 
                      align = 'classif', 
                      rmhomologflank = REMOVEHOMOLOGOUSFLANK, 
                      clustering_method = 'MCL', 
                      WickerParams = None, 
                      auto_outgroup = AUTO_OUTGROUP, 
                      bpflank = bpflank, 
                      combine_and_do_small_clusters = combine_and_do_small_clusters, 
                      flank_pId = flank_pId, 
                      flank_evalue = flank_evalue, 
                      flank_plencutoff = flank_plencutoff, 
                      LTRSONLY = False)
            alnPthKeys.append('Aln_{0}_I{1}_cluster_all.{2}.{3}'.format(
                                                     classif, I, strHomoflank))
    elif align == 'cluster':
        if REMOVEGENECONV:
            MakeDir('GCDir', '{0}/GeneconversionDisallowed'.format(
                                                           paths['RegionDir']))
        else:
            MakeDir('GCDir', '{0}/NoGCFiltering'.format(paths['RegionDir']))
        if REMOVEHOMOLOGOUSFLANK:
            strHomoflank = 'homoflank'
            MakeDir('HomoFlankDir', '{0}/NoPairsWithHomologousFlanks'.format(
                                                               paths['GCDir']))
        else:
            strHomoflank = 'nohomoflank'
            MakeDir('HomoFlankDir', '{0}/AllElements'.format(paths['GCDir']))
        if AUTO_OUTGROUP:
            MakeDir('OutgroupDir', '{0}/WithOutgroup'.format(
                                                        paths['HomoFlankDir']))
        else:
            MakeDir('OutgroupDir', '{0}/NoOutgroup'.format(
                                                        paths['HomoFlankDir']))
        OutPth= paths['OutgroupDir']
        if WICKERCLUST:
            append2logfile(paths['output_top_dir'], mainlogfile, 
              'Start AutoAlign() in phylo() for {0}'.format(clustering_method))
            AutoAlign(I = None, 
                      part = part, 
                      rmgeneconv = removegeneconv, 
                      minClustSize = minClustSize, 
                      align = 'clusters', 
                      rmhomologflank = REMOVEHOMOLOGOUSFLANK, 
                      clustering_method = 'WickerFam', 
                      WickerParams = {'pId':WickerParams['pId'],
                                      'percAln':WickerParams['percAln'],
                                      'minLen':WickerParams['minLen']}, 
                      auto_outgroup = AUTO_OUTGROUP, 
                      bpflank = bpflank, 
                      combine_and_do_small_clusters = combine_and_do_small_clusters, 
                      flank_pId = flank_pId, 
                      flank_evalue = flank_evalue, 
                      flank_plencutoff = flank_plencutoff, 
                      LTRSONLY = False)
        elif MCLCLUST:
            append2logfile(paths['output_top_dir'], mainlogfile, 
              'Start AutoAlign() in phylo() for {0}'.format(clustering_method))
            AutoAlign(I = I, 
                      part = part, 
                      rmgeneconv = removegeneconv, 
                      minClustSize = minClustSize, 
                      align = 'clusters', 
                      rmhomologflank = REMOVEHOMOLOGOUSFLANK, 
                      clustering_method = 'MCL', 
                      WickerParams = None, 
                      auto_outgroup = AUTO_OUTGROUP, 
                      bpflank = bpflank, 
                      combine_and_do_small_clusters = combine_and_do_small_clusters, 
                      flank_pId = flank_pId, 
                      flank_evalue = flank_evalue, 
                      flank_plencutoff = flank_plencutoff, 
                      LTRSONLY = False)

    for classif in classifs:
        append2logfile(paths['output_top_dir'], mainlogfile, 
                    'Start processing {0} clusters in phylo() for {1}'.format(
                                                            classif, 
                                                            clustering_method))
        # Read clusters
        if MCLCLUST:
            clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
        elif WICKERCLUST:
            clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                                WickerParams['pId'], 
                                                WickerParams['percAln'], 
                                                WickerParams['minLen'], 
                                                classif)]
        clusters = [clust.split('\t') for clust in open(
                                 clusterPath, 'r').read().strip().split('\n')]
        # Read GENECONV output
        if REMOVEGENECONV:
            append2logfile(paths['output_top_dir'], mainlogfile, 
              'Start processing gene conversion results phylo() for {0}'.format(
                                                            clustering_method))
            gc = 'GCfiltered'
            append2logfile(paths['output_top_dir'], mainlogfile, 
                'Excluding {0} elements with evidence of gene conversion'.format(
                                                                      classif))
            if MCLCLUST:
                gcSummaryPth = 'MCL_I{0}_GENECONV_summary'.format(I)
            elif WICKERCLUST:
                gcSummaryPth = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONV_summary'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
            if gcSummaryPth in paths:
                if os.path.isfile(gcSummaryPth):
                    gcDct = {}
                    with open(paths[gcSummaryPth], 'r') as gcFl:
                        for line in gcFl:
                            el, clust, classif, g = line.strip().split()
                            el = 'LTR_retrotransposon{0}'.format(el)
                            clust = int(clust)
                            if classif in gcDct:
                                if clust in gcDct[classif]:
                                    if not el in gcDct[classif][clust]:
                                        gcDct[classif][clust].append(el)
                                else:
                                    gcDct[classif][clust] = [el]
                            else:
                                gcDct[classif] = {clust:[el]}
            elif gcSummaryPth not in paths:
                print(('phylo(): modeltest(removegeneconv = True) used but {0} '
                       'not in paths or status file. This happens at least '
                       'when there is only 1 element with a classification and '
                       'therefore nothing is aligned').format(gcSummaryPth), 
                                                             file = sys.stderr)
            if os.path.isfile(gcSummaryPth):
                for j in range(len(clusters)):
                    if j in gcDct[classif]:
                        clusters[j] = [el for el in clusters[j] if 
                                                   not el in gcDct[classif][j]]
        else:
            gc = 'NoGCfiltering'
        if align == 'cluster':
            for j in range(len(clusters)):
                if WICKERCLUST:
                    alnPth = 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_cluster_{4}_{5}.{6}.{7}'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        j, 
                                                        gc, 
                                                        strHomoflank, 
                                                        strOutgroup)
                elif MCLCLUST:
                    alnPth = 'Aln_{0}_I{1}_cluster{2}_{3}.{4}.{5}'.format(
                                                        classif, 
                                                        I, 
                                                        j, 
                                                        gc, 
                                                        strHomoflank, 
                                                        strOutgroup)
                if  alnPth in paths:
                    alnPthKeys.append(alnPth)
            # clustersmall will be left out of OUTGROUP. No LTT for clustersmall
            if combine_and_do_small_clusters:
                if WICKERCLUST:
                    alnPth = 'WickerAln_{0}_pId_{1}_percAln_{2}_minLen_{3}_clustersmall_{4}'.format(WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        gc)
                    if alnPth in paths:
                        alnPthKeys.append(alnPth)
                elif MCLCLUST:
                    alnPth = 'Aln_{0}_I{1}_clustersmall_{2}'.format(classif, I, 
                                                                            gc)
                    if alnPth in paths:
                        alnPthKeys.append(alnPth)
    if BOOTSTRAP:
        append2logfile(paths['output_top_dir'], mainlogfile, 
            'Ready to begin bootstrapping for:\n{0}'.format('\n'.join(
                                                                  alnPthKeys)))
        if ULTRAMETRIC:
            if WICKERCLUST:
                bootstrap(alnPthsLst = alnPthKeys, 
                          reps = bootstrap_reps, 
                          OutPth = OutPth, 
                          convert_to_ultrametric = True, 
                          WickerParams = {'pId':WickerParams['pId'],
                                        'percAln':WickerParams['percAln'],
                                        'minLen':WickerParams['minLen']}, 
                          gc = gc, 
                          strHomoflank = strHomoflank, 
                          strOutgroup = strOutgroup, 
                          I = None)
            elif MCLCLUST:
                bootstrap(alnPthsLst = alnPthKeys, 
                          reps = bootstrap_reps, 
                          OutPth = OutPth, 
                          convert_to_ultrametric = True, 
                          WickerParams = None, 
                          gc = gc, 
                          strHomoflank = strHomoflank, 
                          strOutgroup = strOutgroup, 
                          I = I)
        else:
            if WICKERCLUST:
                bootstrap(alnPthsLst = alnPthKeys, 
                          reps = bootstrap_reps, 
                          OutPth = OutPth, 
                          convert_to_ultrametric = False, 
                          WickerParams = {'pId':WickerParams['pId'],
                                          'percAln':WickerParams['percAln'],
                                          'minLen':WickerParams['minLen']}, 
                          gc = gc, 
                          strHomoflank = strHomoflank, 
                          strOutgroup = strOutgroup, 
                          I = None)
            elif MCLCLUST:
                bootstrap(alnPthsLst = alnPthKeys, 
                          reps = bootstrap_reps, 
                          OutPth = OutPth, 
                          convert_to_ultrametric = False, 
                          WickerParams = None, 
                          gc = gc, 
                          strHomoflank = strHomoflank, 
                          strOutgroup = strOutgroup, 
                          I = I)


def SeqbootCall(Call):
    """Runs seqboot"""
    baseDir = os.getcwd()
    os.chdir(Call[1])
    subprocess.call(Call[0], stdin=open('seqboot.conf', 'r'))
    os.chdir(baseDir)


def bootstrap(alnPthsLst, 
              reps, 
              OutPth = None, 
              convert_to_ultrametric = False, 
              WickerParams = None, 
              gc = None, 
              strHomoflank = None, 
              strOutgroup = None, 
              I = None):
    global paths
    append2logfile(paths['output_top_dir'], mainlogfile, 'Start bootstrap()')
    OutgroupSummaryKey = None
    seqbootCalls = []
    ULTRAMETRIC = convert_to_ultrametric
    # Store calls for SEQBOOT
    for AlnFasta in alnPthsLst:
        if os.path.isfile(paths[AlnFasta]):
            # if alignment file is empty
            if os.stat(paths[AlnFasta]).st_size == 0:
                continue
             # no sequences
            if len(list(SeqIO.parse(paths[AlnFasta], 'fasta'))[0]) == 0: 
                continue
        else:
            continue
        append2logfile(paths['output_top_dir'], mainlogfile, 
                  'Preparing to run SEQBOOT for:\n{0}'.format(paths[AlnFasta]))
        AlignIO.convert(in_file=paths[AlnFasta], 
                        out_file='{0}.phylip'.format(paths[AlnFasta]), 
                        in_format='fasta', 
                        out_format='phylip')
        SeqBootInstructionsFlPth = '{0}/seqboot.conf'.format('/'.join(
                                              paths[AlnFasta].split('/')[:-1]))
        with open(SeqBootInstructionsFlPth, 'w') as SeqBootInstructionsFl:
            Phylip = '{0}.phylip'.format(paths[AlnFasta]).split('/')[-1]
            SeqBootInstructionsFl.write('{0}\nR\n{1}\nY\n1\nR\n'.format(
                                                                 Phylip, reps))
        seqbootCalls.append([['{0}/seqboot'.format(executables['phylip'])], 
                                    '/'.join(paths[AlnFasta].split('/')[:-1])])
    if not seqbootCalls == []:
        chunk_size = ceil(len(seqbootCalls)/procs)
        with Pool(processes=procs) as p:
            p.map(SeqbootCall, seqbootCalls, chunksize=chunk_size)
        p.join()
        append2logfile(paths['output_top_dir'], mainlogfile, 
                                                    'Finished running SEQBOOT')
    for AlnFasta in alnPthsLst:
        fasttreeCalls = []
        if os.path.isfile(paths[AlnFasta]):
            # if alignment file is empty
            if os.stat(paths[AlnFasta]).st_size == 0:
                continue
            # no sequences
            if len(list(SeqIO.parse(paths[AlnFasta], 'fasta'))[0]) == 0:
                continue
        else:
            continue
        alignment_length = len(list(SeqIO.parse(paths[AlnFasta], 'fasta'))[0])
        Alns = list(SeqIO.parse(paths[AlnFasta], 'fasta'))
        numAlns = len(Alns)
        alnLen = len(Alns[0].seq)
        pth = '/'.join(paths[AlnFasta].split('/')[:-1])
        clustMethod = paths[AlnFasta].split('/')[1]
        if clustMethod.startswith('Wicker'):
            clustMethod = 'WickerFam'
        if clustMethod == 'MCL':
            settings  = paths[AlnFasta].split('/')[2]
            classif = paths[AlnFasta].split('/')[-3]
            if classif == 'whole_classif':
                classif = paths[AlnFasta].split('/')[-2]
            clust = paths[AlnFasta].split('/')[-2].split('_')[-1]
            if clust.endswith('.fasta.aln.trimal'):
                clust = 'wholeClassif'
        elif clustMethod == 'WickerFam':
            classif = paths[AlnFasta].split('/')[-3]
            clust = paths[AlnFasta].split('/')[-2]
            settings = paths[AlnFasta].split('/')[2]
        else:
            sys.exit(('modeltest() parameter clustering_method needs to be '
                      'either WickerFam or MCL, and it is, from within '
                      'bootstrap(): {0}').format(clustMethod))
        if ULTRAMETRIC:
            if clustMethod == 'WickerFam':
                OutgroupSummaryKey = 'WickerOutgroups_{0}_pId_{1}_percAln_{2}_minLen_{3}_{4}.{5}.OutgroupSummary'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        gc, 
                                                        strHomoflank)
            elif clustMethod == 'MCL':
                OutgroupSummaryKey = 'MCLOutgroups_{0}_I{1}_{2}.{3}.OutgroupSummary'.format(
                                                        classif, 
                                                        I, 
                                                        gc, 
                                                        strHomoflank)
        seqbootOutputPhylip = '{0}/outfile'.format(pth)
        # Split seqboot multi-phylip output
        MakeDir('classifDir', '{0}/{1}'.format(OutPth, classif))
        MakeDir('clustDir', '{0}/{1}'.format(paths['classifDir'], clust))
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
          'Below log entry is from line {0} in {1}'.format(lineno, scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                 'Preparing to run FastTree for:\n{0}'.format(paths[AlnFasta]))
        # Integrate both Wicker and MCL
        runID = '{0}_{1}_Bootstrap_{2}_{3}'.format(clustMethod, settings, 
                                                                classif, clust)
        # Run FastTree for main tree
        if not 'mainTree_{0}'.format(runID) in paths:
            append2logfile(paths['output_top_dir'], mainlogfile, 
                                    'Start inferring phylogeny in bootstrap()')
            mainTree = '{0}/{1}.main_tree'.format(paths['clustDir'], paths[
                                                      AlnFasta].split('/')[-1])
            makecallMultiprocessing(([executables['fasttree'], '-nt', '-gtr'], 
                    mainTree, '{0}/fasttree.err'.format(pth), paths[AlnFasta]))
            paths['mainTree_{0}'.format(runID)] = mainTree
            with open('{0}/status'.format(paths[
                                    'output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format('mainTree_{0}'.format(
                                                             runID), mainTree))
            scriptpath = os.path.realpath(__file__)
            lineno = getframeinfo(currentframe()).lineno + 2
            append2logfile(paths['output_top_dir'], mainlogfile, 
                     'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                   scriptpath))
            append2logfile(paths['output_top_dir'], mainlogfile, 
                                  'Generated main tree:\n{0}'.format(mainTree))
        if not 'Trees_{0}'.format(runID) in paths:
            with open(seqbootOutputPhylip, 'r') as alnIters:
                replicate = 0
                repFlPth = None
                scriptpath = os.path.realpath(__file__)
                lineno = getframeinfo(currentframe()).lineno + 2
                append2logfile(paths['output_top_dir'], mainlogfile, 
                            'Below log entry is from line {0} in {1}'.format(
                                                           lineno, scriptpath))
                append2logfile(paths['output_top_dir'], mainlogfile, 
                                                 'Preparing for bootstrapping')
                for line in alnIters:
                    # new alignment
                    if line.strip().replace(' ', '') == '{0}{1}'.format(
                                                              numAlns, alnLen):
                        if repFlPth != None:
                            iterDir = '{0}/{1}'.format(paths['clustDir'], 
                                                                     replicate)
                            fasttree_call = ([executables['fasttree'], '-nt', '-gtr'], 
                                            '{0}/tree.newick'.format(iterDir), 
                                            '{0}/fasttree.err'.format(iterDir), 
                                            repFlPth)
                            fasttreeCalls.append(fasttree_call)
                        replicate += 1
                        # make dir for bootstrap replicates output
                        MakeDir('iterDir', '{0}/{1}'.format(paths['clustDir'], 
                                                                    replicate))
                        repFlPth = '{0}/alnReplicate_{1}.fasta'.format(
                                                   paths['iterDir'], replicate)
                        with open(repFlPth, 'a') as repFl:
                            repFl.write(line)
                    else:
                        with open(repFlPth, 'a') as repFl:
                            repFl.write(line)
                fasttree_call = ([executables['fasttree'], '-nt', '-gtr'], 
                        '{0}/tree.newick'.format(paths['iterDir']),
                        '{0}/fasttree.err'.format(paths['iterDir']), repFlPth)
                fasttreeCalls.append(fasttree_call)
            chunk_size = ceil(len(fasttreeCalls)/procs)
            with Pool(processes=procs) as p:
                append2logfile(paths['output_top_dir'], mainlogfile, 
                    'Start FastTree bootstrapping for alignments:\n{0}'.format(
                                                                     AlnFasta))
                p.map(makecallMultiprocessing, fasttreeCalls,
                                                        chunksize = chunk_size)
            p.join()
            paths['Trees_{0}'.format(runID)] = '{0}/allReplicates.newick'.format(
                                                             paths['clustDir'])
            for d in [D for D in os.listdir(paths['clustDir']) if 
                        os.path.isdir('{0}/{1}'.format(paths['clustDir'], D))]:
                for f in os.listdir('{0}/{1}'.format(paths['clustDir'], d)):
                    if f.endswith('newick'):
                        tree = open('{0}/{1}/{2}'.format(paths['clustDir'], d, 
                                                        f), 'r').read().strip()
                        with open(paths['Trees_{0}'.format(runID)], 'a') as treeFl:
                            treeFl.write('{0}\n'.format(tree))
            append2logfile(paths['output_top_dir'], mainlogfile, 
                    'All bootstrap replicate trees in this file:\n{0}'.format(
                                             paths['Trees_{0}'.format(runID)]))
            paths[runID] = paths['clustDir']
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format(runID,  paths[runID]))
            with open('{0}/status'.format(paths[
                                    'output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format('Trees_{0}'.format(
                                     runID), paths['Trees_{0}'.format(runID)]))
        # Get bootstrap values
        append2logfile(paths['output_top_dir'], mainlogfile, 
                       'Start CompareToBootstrap.pl for\n{0}'.format(AlnFasta))
        compare2bootstrap_call = [executables['perl'], 
                                  '{0}/{1}'.format( 
                                  paths['scriptsDir'], 
                                  'CompareToBootstrap.pl'), 
                                  '-tree', 
                                  paths['mainTree_{0}'.format(runID)], 
                                  '-boot', 
                                  paths['Trees_{0}'.format(runID)]]
        makecall(compare2bootstrap_call, stdout = '{0}.bootstrapped'.format(
                                          paths['mainTree_{0}'.format(runID)]))
        bootstrapped = '{0}.bootstrapped'.format(paths['mainTree_{0}'.format(
                                                                       runID)])
        # run PATHd8 to convert tree to ultrametric
        if ULTRAMETRIC:
            append2logfile(paths['output_top_dir'], mainlogfile, 
                  'Beginning Ultrametric transformation for {0}'.format(runID))
            with open(bootstrapped, 'r') as bstreefl:
                tree = bstreefl.read().strip()
            outgroup = None
            with open(paths[OutgroupSummaryKey], 'r') as outgroup_file:
                for line in outgroup_file:
                    if line == 'Alignment\toutgroup\toutgroup_cluster\n':
                        continue
                    key, outG, outgroup_clust = line.strip().split('\t')
                    if key == AlnFasta:
                        outgroup = outG[19:]
            if outgroup == None:
                print('Unable to root {0}. No ultrametric tree will be made'.format(
                                                                     AlnFasta))
            else:
                # a random element from the current cluster to tell PATHd8 
                # where the mrca is (outgroupXother_taxon) for determining 
                # relative branch lengths
                other_taxon = random.choice(list(SeqIO.parse(paths[AlnFasta], 
                                                                  'fasta'))).id 
                append2logfile(paths['output_top_dir'], mainlogfile, 
                                      'Beginning PATHd8 for {0}'.format(runID))
                pathd8_file_str = """Sequence length={0};
mrca: {1}, {2}, fixage=1;
{3}
    """.format(alignment_length, outgroup, other_taxon, tree)
                pathd8flpth = '{0}/pathd8.in'.format(paths[runID])
                with open(pathd8flpth, 'w') as outFl:
                    outFl.write(pathd8_file_str)
                paths['PATHd8_output_{0}'.format(
                        runID)] = '{0}.pathd8_ultrametric'.format(bootstrapped)
                pathd8_call = [executables['pathd8'], pathd8flpth, 
                                      paths['PATHd8_output_{0}'.format(runID)]]
                makecall(pathd8_call)
                if os.path.isfile(paths['PATHd8_output_{0}'.format(runID)]):
                    with open(paths['PATHd8_output_{0}'.format(runID)], 'r') as d8Fl:
                        for line in d8Fl:
                            if line.startswith('d8 tree '):
                                d8tree = line[line.index('('):]
                                newLocation = '{0}/{1}_{2}.bootstrapped.pathd8_ultrametric.outgroup_{3}.newick'.format(
                                                         paths['classifDir'], 
                                                         classif, 
                                                         clust, 
                                                         outgroup)
                                if clustMethod == 'WickerFam':
                                    newLocationKey = 'WickerOutgroups_{0}_pId_{1}_percAln_{2}_minLen_{3}_{4}.{5}.outgroup_{6}'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif, 
                                                        gc, 
                                                        strHomoflank, 
                                                        outgroup)
                                elif clustMethod == 'MCL':
                                    newLocationKey = 'MCLOutgroups_{0}_I{1}_{2}.{3}.outgroup_{4}'.format(
                                                        classif, 
                                                        I, 
                                                        gc, 
                                                        strHomoflank, 
                                                        outgroup)
                                with open(newLocation, 'w') as outFl:
                                    outFl.write(d8tree)
                                if not checkStatusFl(newLocationKey):
                                    with open('{0}/status'.format(paths[
                                          'output_top_dir']), 'a') as statusFl:
                                        statusFl.write('{0}\t{1}\n'.format(
                                                  newLocationKey, newLocation))
                append2logfile(paths['output_top_dir'], mainlogfile, 
                                       'Finished PATHd8 for {0}'.format(runID))
        bootstrappedLocation = '{0}/{1}_{2}.bootstrapped.outgroup_{3}.newick'.format(
                                  paths['classifDir'], classif, clust, classif)
        copyfile(bootstrapped, bootstrappedLocation)
        append2logfile(paths['output_top_dir'], mainlogfile, (('Finished bootstrapping. '
                    'Tree with bootstrap support values here:\n{0}')).format(
                                                         bootstrappedLocation))


def clusterSummary():
    """Creates a map of LTR-RT IDs to the superfamily classification
    they were assigned in the classify step.
    """
    clustPaths = [line for line in open('{0}/status'.format(
                    paths['output_top_dir']), 'r').read().strip().split('\n') 
                   if not 'divergence' in line and not line.endswith('minLen') 
                   and (line.startswith('WickerFamDir_') 
                   or re.match('MCL_.+_I\d', line))]
    # To avoid double-appending, remove cluster files if they exist 
    # before adding the new stuff
    for p in clustPaths:
        p = p.strip().split('\t')
        if p[0].startswith('MCL'):
            I = p[0].split('_')[-1][1:]
            MCLclustdir = '/'.join(p[1].split('/')[:-2])
            if os.path.exists('{0}/MCL_I{1}_summary'.format(MCLclustdir, I)):
                os.remove('{0}/MCL_I{1}_summary'.format(MCLclustdir, I))
        elif p[0].startswith('Wicker'):
            Wickerclustdir = '/'.join(p[1].split('/')[:-2])
            params = '_'.join(p[0].split('_')[1:-2])
            if os.path.exists('{0}/Wicker_{1}_summary'.format(Wickerclustdir, 
                                                                      params)):
                os.remove('{0}/Wicker_{1}_summary'.format(Wickerclustdir, 
                                                                       params))
    # Holds whether the header was written already for a particular file
    headers = set()
    for p in clustPaths:
        p = p.strip().split('\t')
        if p[0].startswith('MCL'):
            MCLclustdir = '/'.join(p[1].split('/')[:-2])
            superfamily = p[0].split('_')[1]
            I = p[0].split('_')[-1][1:]
            with open('{0}/MCL_I{1}_summary'.format(
                                                MCLclustdir, I), 'a') as outfl:
                if '{0}/MCL_I{1}_summary'.format(MCLclustdir, I) not in headers:
                    outfl.write('superfamily\tcluster\tsize\n')
                    headers.add('{0}/MCL_I{1}_summary'.format(MCLclustdir, I))
                with open(p[1], 'r') as clustFl:
                    clust = 0
                    for line in clustFl:
                        outfl.write('{0}\t{1}\t{2}\n'.format(superfamily, clust, 
                                                        len(line.split('\t'))))
                        clust += 1
        elif p[0].startswith('Wicker'):
            superfamily = p[0].split('_')[-1]
            Wickerclustdir = '/'.join(p[1].split('/')[:-2])
            params = '_'.join(p[0].split('_')[1:-2])
            with open('{0}/Wicker_{1}_summary'.format(
                                        Wickerclustdir, params), 'a') as outfl:
                if '{0}/Wicker_{1}_summary'.format(
                                        Wickerclustdir, params) not in headers:
                    outfl.write('superfamily\tcluster\tsize\n')
                    headers.add('{0}/Wicker_{1}_summary'.format(
                                                       Wickerclustdir, params))
                with open(p[1], 'r') as clustFl:
                    clust = 0
                    for line in clustFl:
                        outfl.write('{0}\t{1}\t{2}\n'.format(
                                    superfamily, clust, len(line.split('\t'))))
                        clust += 1


def geneconv2circoslinks(geneconvfile, 
                         ltrharvestgff, 
                         outfile, 
                         append = False, 
                         output = 'file', 
                         linksdct = None, 
                         transposeLinks = True):
    """Converts GI tract pairs from geneconvClusters() output and writes 
    a links file for Circos. The GFF3 is needed to get the scaffold 
    name. seqlengths needs to be a dictionary with the lengths of the 
    sequences whose names correspond to the sequence names in the gff 
    for the features with gene conversion tracts. Assumes 
    LTR_retrotransposon features were used. 
    
    append=True will append to the outfile if it exists.

    output    'file', or 'return'. If file, the links will be written 
                    to a file at outfile. If 'return' then the links 
                    will be returned as a dictionary like:
                links = {orf1:[linkline1, linkline2, linkline3], ... }
                     Only the first orf in each link pair needs to be 
                     considered because gene conversion was evaluated 
                     within clusters only.
    
    linksdct    If 'output'=return and a links dct is provided here it 
                will be updated and returned.

    transposeLinks=True will write links for elements in their 
                        position on the scaffolds.
    transposeLinks=False will write links for elements with a start 
                         position of 0 at one end of the element.
    """
    global paths
    seqs = {}
    starts = {}
    if append:
        mode = 'a'
    else:
        mode = 'w'
    with open(ltrharvestgff, 'r') as inFl:
        for line in inFl:
            if line.startswith('#'):
                continue
            if '\tLTR_retrotransposon\t' in line:
                gffLine = GFF3_line(line)
                element = gffLine.attributes['ID']
                scaf = gffLine.seqid
                seqs[element] = scaf
                starts[element] = int(gffLine.start)
    
    if output == 'file':
        with open(outfile, mode) as outFl:
            with open(geneconvfile, 'r') as inFl:
                for line in inFl:
                    if line.startswith('GI'):
                        rec = line.strip().split('\t')
                        el1, el2 = ['LTR_retrotransposon{0}'.format(
                                             e[1:]) for e in rec[1].split(';')]
                        if transposeLinks:
                            el1start  = int(rec[7]) + starts[el1] - 1
                            el1end  = int(rec[8]) + starts[el1] - 1
                            el2start  = int(rec[10]) + starts[el2] - 1
                            el2end  = int(rec[11]) + starts[el2] - 1
                            el1seq = copy(seqs[el1])
                            el2seq = copy(seqs[el2])
                        else:
                            el1start  = int(rec[7])
                            el1end  = int(rec[8])
                            el2start  = int(rec[10])
                            el2end  = int(rec[11])
                            el1seq = copy(el1)
                            el2seq = copy(el2)
                        if 'g0.summary' in geneconvfile:
                            outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=3,color=orange_a3\n'.format(
                                el1seq, el1start, el1end, el2seq, el2start, el2end))
                        elif 'g1.summary' in geneconvfile:
                            outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=2,color=vdpurple_a5\n'.format(
                                el1seq, el1start, el1end, el2seq, el2start, el2end))
                        elif 'g2.summary' in geneconvfile:
                            outFl.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=1,color=vlblue_a1\n'.format(
                                el1seq, el1start, el1end, el2seq, el2start, el2end))
    elif output == 'return':
        if linksdct != None:
            links = linksdct
            links_untransposed = linksdct
        else:
            links = {}
            links_untransposed = {}
        with open(geneconvfile, 'r') as inFl:
            
            el1start, el1end, el2start, el2end, el1seq, el2seq = [None]*6

            for line in inFl:
                if line.startswith('GI'):
                    rec = line.strip().split('\t')
                    el1, el2 = ['LTR_retrotransposon{0}'.format(
                                             e[1:]) for e in rec[1].split(';')]
                    if transposeLinks:
                        el1start  = int(rec[7]) + starts[el1] - 1
                        el1end  = int(rec[8]) + starts[el1] - 1
                        el2start  = int(rec[10]) + starts[el2] - 1
                        el2end  = int(rec[11]) + starts[el2] - 1
                        el1seq = copy(seqs[el1])
                        el2seq = copy(seqs[el2])
                        # Different colored links for different gscale 
                        # parameters. g values > 2 are possible but not 
                        # implemented.
                        if 'g0.summary' in geneconvfile:
                            outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=3,color=orange_a3\n'.format(
                                el1seq, el1start, el1end, el2seq, el2start, el2end)
                            g = 'g0'
                        elif 'g1.summary' in geneconvfile:
                            outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=2,color=vdpurple_a5\n'.format(
                                el1seq, el1start, el1end, el2seq, el2start, el2end)
                            g = 'g1'
                        elif 'g2.summary' in geneconvfile:
                            outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=1,color=vlblue_a1\n'.format(
                                el1seq, el1start, el1end, el2seq, el2start, el2end)
                            g = 'g2'
                        if g in links:
                            if el1 in links[g]:
                                links[g][el1].append(outline)
                            else:
                                links[g][el1] = [outline]
                        else:
                            links[g] = {el1:[outline]}
                    else:
                        el1start  = int(rec[7])
                        el1end  = int(rec[8])
                        el2start  = int(rec[10])
                        el2end  = int(rec[11])
                        # Different colored links for different gscale 
                        # parameters. g values > 2 are possible but not 
                        # implemented.
                        if 'g0.summary' in geneconvfile:
                            outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=3,color=orange_a3\n'.format(
                                el1, el1start, el1end, el2, el2start, el2end)
                            g = 'g0'
                        elif 'g1.summary' in geneconvfile:
                            outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=2,color=vdpurple_a5\n'.format(
                                el1, el1start, el1end, el2, el2start, el2end)
                            g = 'g1'
                        elif 'g2.summary' in geneconvfile:
                            outline = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tz=1,color=vlblue_a1\n'.format(
                                el1, el1start, el1end, el2, el2start, el2end)
                            g = 'g2'
                        if g in links:
                            if el1 in links[g]:
                                links[g][el1].append(outline)
                            else:
                                links[g][el1] = [outline]
                        else:
                            links[g] = {el1:[outline]}
        return links


def circosMultiprocessing(packet):
    """
    For using with multiprocessing.Pool()
    """
    global paths

    circosdir = packet[0]
    circos_call = packet[1]
    classif = packet[2]
    i = packet[3]
    outpath = packet[4]
    G = packet[5]
    current_wd = os.getcwd()
    os.chdir(circosdir)
    subprocess.call(circos_call, stdout=open('out', 'w'), stderr=open('err','w'))
    os.chdir(current_wd)
    png = '{0}/circos.png'.format(circosdir)
    svg = '{0}/circos.svg'.format(circosdir)
    newpng = '{0}/{1}.cluster_{2}.geneconv_{3}.png'.format(
                                              outpath, classif, i, '_'.join(G))
    newsvg = '{0}/{1}.cluster_{2}.geneconv_{3}.svg'.format(
                                              outpath, classif, i, '_'.join(G))
    if not os.path.isfile(newpng):
        if os.path.isfile(png):
            copyfile(png, newpng)
    if not os.path.isfile(newsvg):
        if os.path.isfile(svg):
            copyfile(svg, newsvg)
    # Clean up
    #if not KEEP_UNUSED_FILES:
    #    rmtree(circosdir)


def Circos(window = '1000000', 
           plots = 'clusters', 
           I = 6, 
           clustering_method = 'WickerFam', 
           WickerParams = {'pId':80,'percAln':80,'minLen':80}, 
           g = 'g0,g1,g2', 
           MinCircosClusterSize = 2):
    """
    Generate a Circos plot for each cluster, showing gene interelement 
    gene conversion tracts by using links. plots needs to be either 
    'clusters' or 'classifs'
    """
    global paths
    G = g.split(',')
    CLASSIFS = False
    CLUSTERS = False
    WICKERCLUST = False
    MCLCLUST = False
    if plots == 'classifs':
        CLASSIFS = True
    elif plots == 'clusters':
        CLUSTERS = True
    if clustering_method == 'WickerFam':
        WICKERCLUST = True
        ClustMethod = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                    WickerParams['pId'], 
                                                    WickerParams['percAln'], 
                                                    WickerParams['minLen'])
    elif clustering_method == 'MCL':
        MCLCLUST = True
        ClustMethod =  'MCL_I{0}'.format(I)
    if checkStatusFl('Circos_output_dir_elements_{0}'.format(
                    ClustMethod)) and checkStatusFl(
                        'Circos_output_dir_scaffolds_{0}'.format(ClustMethod)):
        append2logfile(paths['output_top_dir'], mainlogfile, 
                     'Circos() {0} has already completed.'.format(ClustMethod))
        return
    append2logfile(paths['output_top_dir'], mainlogfile, 
                                               'Beginning making Circos plots')
    MakeDir('CircosTopDir', '{0}/Circos'.format(paths['output_top_dir']))
    # Get seq lengths
    scafLengthsFlPth = '{0}/seqLengths.tab'.format(paths['CircosTopDir'])
    scafLengths = {}
    with open(paths['inputFasta'], 'r') as inFl:
        currentSeqName = None
        for line in inFl:
            if line.startswith('>'):
                currentSeqName = line.strip()[1:].split(' ')[0]
                scafLengths[currentSeqName] = 0
            else:
                scafLengths[currentSeqName] += len(line.strip())
    with open(scafLengthsFlPth, 'w') as outFl:
        for scaf in sorted(list(scafLengths.keys())):
            outFl.write('{0}\t{1}\n'.format(scaf, str(scafLengths[scaf])))
    # Create ideogram file!
    allscafs = '{0}/seqs.track'.format(paths['CircosTopDir'])
    ideogramCall = ['{0}/ideogramFromLengths.py'.format(paths['scriptsDir'])]
    makecall(ideogramCall, stdin = scafLengthsFlPth, stdout = allscafs)
    # Make track for elements!
    if CLASSIFS:
        # Separate out GFFs by classif
        allGFFoutPth = '{0}/all.gff'.format(paths['CircosTopDir'])
        with open(paths['CurrentGFF']) as gffFl:
            for line in gffFl:
                if '\tLTR_retrotransposon\t' in line:
                    gffLine = GFF3_line(line)
                    classif = classifs_by_element[gffLine.attributes['ID']]
                    MakeDir('classifDir', '{0}/{1}'.format(
                                               paths['CircosTopDir'], classif))
                    GFFoutPth = '{0}/{1}.gff'.format(
                                                  paths['classifDir'], classif)
                    with open(GFFoutPth, 'a') as GFFoutFl:
                        GFFoutFl.write(line)
                    with open(allGFFoutPth, 'a') as GFFoutFl:
                        GFFoutFl.write(line)
        scriptpath = os.path.realpath(__file__)
        lineno = getframeinfo(currentframe()).lineno + 2
        append2logfile(paths['output_top_dir'], mainlogfile, 
                     'Below log entry is from line {0} in {1}'.format(lineno, 
                                                                   scriptpath))
        append2logfile(paths['output_top_dir'], mainlogfile, ('Created GFF '
            'files for each classification for converting to Circos heatmap tracks.'))
        for classif in classifs:
            classifDir = '{0}/{1}'.format(paths['CircosTopDir'], classif)
            GFFoutPth = '{0}/{1}.gff'.format(classifDir, classif)
            gff2heatmapCall = ['{0}/gff2circos-heatmap.py'.format(
                                                  paths['scriptsDir']), 
                                                  '-gff', GFFoutPth, 
                                                  '-window', window, 
                                                  '-scafLens', scafLengthsFlPth]
            makecall(gff2heatmapCall, stdout = '{0}/{1}.heatmap.track'.format(
                                                          classifDir, classif))
            # Geneconv output to circos links track
            paths['GENECONV_{0}_dir'.format(classif)] = '{0}/{1}'.format(
                                             paths[geneconvOutputDir], classif)
            outfile = '{0}/{1}.testlinks'.format(paths['CircosTopDir'], classif)
            g0fl = '{0}/{1}_{2}.summary'.format(
                       paths['GENECONV_{0}_dir'.format(classif)],classif, 'g0')
            g1fl = '{0}/{1}_{2}.summary'.format(
                       paths['GENECONV_{0}_dir'.format(classif)],classif, 'g1')
            g2fl = '{0}/{1}_{2}.summary'.format(
                       paths['GENECONV_{0}_dir'.format(classif)],classif, 'g2')
            if os.path.isfile(g0fl):
                # Convert GENECONV output to Circos links track
                geneconv2circoslinks(g0fl, paths['CurrentGFF'], outfile)
            if os.path.isfile(g1fl):
                # Convert GENECONV output to Circos links track
                geneconv2circoslinks(g1fl, paths['CurrentGFF'], outfile, 
                                                                 append = True)
            if os.path.isfile(g2fl):
                # Convert GENECONV output to Circos links track
                geneconv2circoslinks(g2fl, paths['CurrentGFF'], outfile, 
                                                                 append = True)
            append2logfile(paths['output_top_dir'], mainlogfile, 
                            ('Created links tracks for Circos from intra-cluster'
                            ' inter-element GENECONV output'))
        gff2heatmapCall = ['{0}/gff2circos-heatmap.py'.format(
                                                paths['scriptsDir']), 
                                                '-gff', allGFFoutPth, 
                                                '-window', window, 
                                                '-scafLens', scafLengthsFlPth]
        makecall(gff2heatmapCall, stdout='{0}/all.heatmap.track'.format(
                                                        paths['CircosTopDir']))
        append2logfile(paths['output_top_dir'], mainlogfile, 
                                 'Converted GFFs to heatmap tracks for Circos')
    elif CLUSTERS:
        if WICKERCLUST:
            WickerDir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])]
            paths['Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(
                    WickerParams['pId'], 
                    WickerParams['percAln'], 
                    WickerParams['minLen'])] = '{0}/GENECONV'.format(WickerDir)
            WickerGCdirkey = 'Wicker_{0}_pId_{1}_percAln_{2}_minLen_GENECONVdir'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
            if not checkStatusFl(WickerGCdirkey):
                sys.exit('Circos() not possible: geneconvClusters() not done yet.')
            geneconvOutputDir = WickerGCdirkey
            MakeDir('CurrentTopDir', '{0}/WickerFam_{1}_pId_{2}_percAln_{3}_minLen'.format(
                                            paths['CircosTopDir'], 
                                            WickerParams['pId'], 
                                            WickerParams['percAln'], 
                                            WickerParams['minLen']))
        elif MCLCLUST:
            MCLdir = paths['MCL_I{0}'.format(I)]
            paths['MCL_I{0}_GENECONVdir'.format(I)] = '{0}/GENECONV'.format(MCLdir)
            if not checkStatusFl('MCL_I{0}_GENECONVdir'.format(I)):
                sys.exit('Circos() not possible: geneconvClusters() not done yet.')
            geneconvOutputDir = 'MCL_I{0}_GENECONVdir'.format(I)
            MakeDir('CurrentTopDir', '{0}/MCL_I{1}'.format(paths['CircosTopDir'], I))
        # Create a Circos plot for each cluster
        heatmapcalls = []
        tilecalls = []
        circoscallsScafs = []
        circoscallsElmts = []
        for classif in classifs:
            # Geneconv output to circos links track
            paths['GENECONV_{0}_dir'.format(classif)] = '{0}/{1}'.format(
                                             paths[geneconvOutputDir], classif)
            outfile = '{0}/{1}.testlinks'.format(paths['CurrentTopDir'], classif)
            outfile_untransposed = '{0}/{1}.testlinks.untransposed'.format(
                                               paths['CurrentTopDir'], classif)
            g0fl = '{0}/{1}_{2}.summary'.format(
                       paths['GENECONV_{0}_dir'.format(classif)],classif, 'g0')
            g1fl = '{0}/{1}_{2}.summary'.format(
                       paths['GENECONV_{0}_dir'.format(classif)],classif, 'g1')
            g2fl = '{0}/{1}_{2}.summary'.format(
                       paths['GENECONV_{0}_dir'.format(classif)],classif, 'g2')
            links = {}
            links_untransposed = {}
            if os.path.isfile(g0fl):
                # Convert GENECONV output to Circos links track
                links = geneconv2circoslinks(g0fl, 
                                             paths['CurrentGFF'], 
                                             outfile, 
                                             append = False, 
                                             output = 'return', 
                                             linksdct = None)
                links_untransposed  =  geneconv2circoslinks(g0fl, 
                                             paths['CurrentGFF'], 
                                             outfile, 
                                             append = False, 
                                             output = 'return', 
                                             linksdct = None, 
                                             transposeLinks = False)
            if os.path.isfile(g1fl):
                # Convert GENECONV output to Circos links track
                links = geneconv2circoslinks(g1fl, 
                                             paths['CurrentGFF'], 
                                             outfile, 
                                             append = True, 
                                             output = 'return', 
                                             linksdct = links)
                links_untransposed  =  geneconv2circoslinks(g1fl, 
                                             paths['CurrentGFF'], 
                                             outfile, 
                                             append = True, 
                                             output = 'return', 
                                             linksdct = links_untransposed, 
                                             transposeLinks = False)
            if os.path.isfile(g2fl):
                # Convert GENECONV output to Circos links track
                links  =  geneconv2circoslinks(g2fl, 
                                             paths['CurrentGFF'], 
                                             outfile, 
                                             append = True, 
                                             output = 'return', 
                                             linksdct = links)
                links_untransposed  =  geneconv2circoslinks(g2fl, 
                                             paths['CurrentGFF'], 
                                             outfile, 
                                             append = True, 
                                             output = 'return', 
                                             linksdct = links_untransposed, 
                                             transposeLinks = False)
            append2logfile(paths['output_top_dir'], mainlogfile, 
                    ('Created links tracks for Circos from intra-cluster '
                    'inter-element GENECONV output'))
            G_incl = [ ['g0'], ['g1'], ['g2'], ['g0', 'g1', 'g2'] ]
            for G in G_incl:
                if MCLCLUST:
                    clusterPath =  paths['MCL_{0}_I{1}'.format(classif, I)]
                elif WICKERCLUST:
                    clusterPath = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                             WickerParams['pId'], 
                                             WickerParams['percAln'], 
                                             WickerParams['minLen'], 
                                             classif)]
                clusters = [clust.split('\t') for clust in open(clusterPath,'r').read().strip().split('\n')]
                totallengths = {}
                element_coords = {}
                for i in range(len(clusters)):
                    if len(clusters[i]) < MinCircosClusterSize:
                        continue
                    clusterscafs = set()
                    outputlinks = []
                    outputlinks_untransposed = []
                    GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(
                                            paths['CurrentTopDir'], classif, i)
                    if os.path.isfile(GFFoutPth):
                        os.remove(GFFoutPth)
                    highlights_ltrs_fl = '{0}/{1}.cluster_{2}.LTR_highlights.track'.format(
                                            paths['CurrentTopDir'], classif, i)
                    if os.path.isfile(highlights_ltrs_fl):
                        os.remove(highlights_ltrs_fl)
                    with open(paths['CurrentGFF']) as gffFl:
                        for line in gffFl:
                            if '\tLTR_retrotransposon\t' in line:
                                gffLine = GFF3_line(line)
                                start = int(gffLine.start)
                                end = int(gffLine.end)
                                name = gffLine.attributes['ID']
                                element_coords[name] = (start, end)
                                el = gffLine.attributes['ID']
                                # Only add links from elements in i
                                if el not in clusters[i]:
                                    continue
                                # Add link to output links
                                for g in G:
                                    if g in links:
                                        if el in links[g]:
                                            outputlinks += links[g][el]
                                            outputlinks_untransposed += links_untransposed[g][el]
                                scaf = gffLine.seqid
                                if scaf not in clusterscafs:
                                    clusterscafs.add(scaf)
                                with open(GFFoutPth, 'a') as GFFoutFl:
                                    GFFoutFl.write(line)
                            # Write highlights track
                            elif '\tlong_terminal_repeat\t' in line:
                                GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(
                                            paths['CurrentTopDir'], classif, i)
                                gffLine = GFF3_line(line)
                                start = int(gffLine.start)
                                end = int(gffLine.end)
                                name = gffLine.attributes['Parent']
                                newstart = start - element_coords[name][0] + 1
                                newend = end - element_coords[name][0] + 1
                                with open(highlights_ltrs_fl, 'a') as outFl:
                                    outFl.write('{0}\t{1}\t{2}\tfill_color=black\n'.format(
                                                       name, newstart, newend))
                    gff2heatmapCallPacket = (['{0}/gff2circos-heatmap.py'.format(
                                                           paths['scriptsDir']), 
                                                '-gff', 
                                                GFFoutPth, 
                                                '-window', 
                                                window, 
                                                '-scafLens', 
                                                scafLengthsFlPth], 
                                                '{0}.heatmap.track'.format(
                                                        GFFoutPth), None, None)
                    gff2tileCallPacket = (['{0}/gff2circos-tile.py'.format(
                                                paths['scriptsDir']), 
                                                '-valueDef', 
                                                'LTR', 
                                                '-gff', 
                                                GFFoutPth], 
                                                '{0}.tile.track'.format(
                                                        GFFoutPth), None, None)
                    tilecalls.append(gff2tileCallPacket)
                    append2logfile(paths['output_top_dir'], mainlogfile, 
                        'gff2circos-heatmap.py:\n{0}'.format(' '.join(
                                                    gff2heatmapCallPacket[0])))
                    append2logfile(paths['output_top_dir'], mainlogfile, 
                        'gff2circos-tile.py:\n{0}'.format(' '.join(
                                                       gff2tileCallPacket[0])))
                    heatmapcalls.append(gff2heatmapCallPacket)
                    with open('{0}/{1}.cluster_{2}.geneconv_{3}.links.track'.format(
                                    paths['CurrentTopDir'], classif, 
                                                i, '_'.join(G)), 'w') as outFl:
                        outFl.write('\n'.join(outputlinks))
                    with open('{0}/{1}.cluster_{2}.geneconv_{3}.links_untransposed.track'.format(
                                    paths['CurrentTopDir'], classif, i, 
                                                   '_'.join(G)), 'w') as outFl:
                        outFl.write('\n'.join(outputlinks_untransposed))
                    # Write ideogram file for just scafs for this cluster
                    #chr - Sacu_v1.1_s0011    11    0    2262239    greys-6-seq-4
                    totalseq = 0
                    with open(allscafs, 'r') as inFl:
                        ideoOut  = '{0}/{1}.cluster_{2}.seq.track'.format(
                                            paths['CurrentTopDir'], classif, i)
                        with open(ideoOut, 'w') as outFl:
                            for line in inFl:
                                scaf = line.split()[2]
                                if scaf in clusterscafs:
                                    contents = line.split()
                                    totalseq += int(contents[-2]) - int(
                                                                  contents[-3])
                                    outFl.write(line)
                    totallengths[i] = totalseq
                    append2logfile(paths['output_top_dir'], mainlogfile, 
                            ('Created GFF files for each classification for '
                             'converting to Circos heatmap tracks.'))
                chunk_size = ceil(len(heatmapcalls)/procs)
                with Pool(processes = procs) as p:
                    p.map(makecallMultiprocessing, heatmapcalls, 
                                                          chunksize = chunk_size)
                p.join()
                append2logfile(paths['output_top_dir'], mainlogfile, 
                                'Converted GFFs to heatmap tracks for Circos.')
                chunk_size = ceil(len(tilecalls)/procs)
                with Pool(processes = procs) as p:
                    p.map(makecallMultiprocessing, tilecalls, 
                                                          chunksize = chunk_size)
                p.join()
                append2logfile(paths['output_top_dir'], mainlogfile, 
                                   'Converted GFFs to tile tracks for Circos.')
                # Circos plot 1: ideograms are scaffolds
                for i in range(len(clusters)):
                    if len(clusters[i]) < MinCircosClusterSize:
                        continue
                    GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(
                                            paths['CurrentTopDir'], classif, i)
                    tilefl = '{0}.tile.track'.format(GFFoutPth)
                    textfl = '{0}.tile.text.track'.format(GFFoutPth)
                    linksfl = '{0}/{1}.cluster_{2}.geneconv_{3}.links.track'.format(
                               paths['CurrentTopDir'], classif, i, '_'.join(G))
                    seqfl = '{0}/{1}.cluster_{2}.seq.track'.format(
                                            paths['CurrentTopDir'], classif, i)
                    # If no links for this cluster, don't draw a Circos 
                    # plot for the elements without scaffolds.unless 
                    # this is the composite with all gscale values, to 
                    # ensure at least one plot is drawn for each cluster.
                    if os.stat(linksfl).st_size == 0 and not G == ['g0', 'g1', 'g2']:
                        continue
                    if (os.path.isfile(tilefl) and os.path.isfile(linksfl) 
                        and os.path.isfile(seqfl)):
                        # Files exist. copy and run Circos.
                        # Plot with scaffolds
                        circosdir = '{0}/circos.{1}.cluster_{2}.geneconv_{3}'.format(
                               paths['CurrentTopDir'], classif, i, '_'.join(G))
                        if not os.path.exists(circosdir):
                            # copy circos conf files and dir structure
                            copytree('{0}/circos'.format(
                                               paths['scriptsDir']), circosdir)
                        newtilefl = '{0}/data/{1}'.format(circosdir, 
                                                         tilefl.split('/')[-1])
                        if not os.path.isfile(newtilefl):
                            copyfile(tilefl, newtilefl)
                        newlinksfl = '{0}/data/{1}'.format(circosdir, 
                                                        linksfl.split('/')[-1])
                        if not os.path.isfile(newlinksfl):
                            copyfile(linksfl, newlinksfl)
                        newseqfl = '{0}/data/{1}'.format(circosdir, 
                                                          seqfl.split('/')[-1])
                        if not os.path.isfile(newseqfl):
                            copyfile(seqfl, newseqfl)
                        newtextfl = '{0}/data/{1}'.format(circosdir, 
                                                         textfl.split('/')[-1])
                        if not os.path.isfile(newtextfl):
                            copyfile(tilefl, newtextfl)
                        conffl = '{0}/etc/circos.conf'.format(circosdir)
                        confbasename = conffl.split('/')[-1]
                        tileblock = """
<plot>
type    =    tile
thickness    =    30
file    =    data/{0}
color    =    dorange
r1    =    0.84r
r0    =    0.78r
</plot>
""".format(newtilefl.split('/')[-1])
                        glyphblock = """
<plot>
type    =    scatter
glyph    =    circle
glyph_size = 60
file    =    data/{0}
color    =    vdorange
orientation = out
r1    =    0.80r
r0    =    0.80r
</plot>
""".format(newtilefl.split('/')[-1])
                        if totallengths[i] > 5000000:
                            plotblock = glyphblock
                        else:
                            plotblock = tileblock
                        circos_conf_str = """<<include colors_fonts_patterns.conf>>
<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

karyotype   = data/{0}

chromosomes_units           = 1000000

<plots>
<plot>
type             = text
color            = black
file             = data/{1}

r0 = 0.84r
r1 = 0.99r

show_links     = no
link_dims      = 0p,10p,60p,10p,0p
link_thickness = 10p
link_color     = red

label_snuggle         = yes
max_snuggle_distance  = 1r
snuggle_tolerance     = 0.40r
snuggle_sampling      = 2
snuggle_link_overlap_test = yes
snuggle_link_overlap_tolerance = 2p
snuggle_refine        = yes
label_rotate = yes
label_size   = 100p
label_font   = condensed

padding  = 1p
rpadding = 1p

</plot>

{2}


</plots>

<links>

radius = 0.78r
crest  = 1
ribbon           = yes
flat             = yes
stroke_color     = vdgrey
stroke_thickness = 2
color            = grey_a3

bezier_radius        = 0r
bezier_radius_purity = 0.5

<link>
file       = data/{3}
</link>

</links>

<<include etc/housekeeping.conf>>
data_out_of_range* = trim""".format(newseqfl.split('/')[-1], 
                                    newtextfl.split('/')[-1], 
                                    plotblock, 
                                    newlinksfl.split('/')[-1])
                        with open(conffl, 'w') as outFl:
                            outFl.write(circos_conf_str)
                        
                        confbasename = conffl.split('/')[-1]
                        imagesize = totallengths[i]/10
                        if imagesize > 6000:
                            imagesize = 6000
                        conffl = '{0}/etc/image.generic.conf'.format(circosdir)
                        circosimageconfstr = """dir   = . 
file  = circos.png
png   = yes
svg   = yes

# radius of inscribed circle in image
radius = {0}

# by default angle=0 is at 3 o'clock position
angle_offset      = -96

#angle_orientation = counterclockwise

auto_alpha_colors = yes
auto_alpha_steps  = 5""".format(imagesize)
                        with open(conffl, 'w') as outFl:
                            outFl.write(circosimageconfstr)
                        circos_call = [executables['perl'], executables['circos']]
                        circoscallsScafs.append([circosdir, 
                                                 circos_call, 
                                                 classif, 
                                                 i, 
                                                 '{0}/plots.scaffolds'.format(
                                                   paths['CurrentTopDir']), 
                                                 G])
                # Circos plot 2: ideograms are elements
                for i in range(len(clusters)):
                    if len(clusters[i]) < MinCircosClusterSize:
                        continue
                    GFFoutPth  = '{0}/{1}.cluster_{2}.gff'.format(paths['CurrentTopDir'], classif, i)
                    tilefl = '{0}.tile.track'.format(GFFoutPth)
                    seqfl = '{0}/{1}.cluster_{2}.seq.track'.format(paths['CurrentTopDir'], classif, i)
                    highlights_ltrs_fl = '{0}/{1}.cluster_{2}.LTR_highlights.track'.format(paths['CurrentTopDir'], classif, i)
                    links_untransposedfl = '{0}/{1}.cluster_{2}.geneconv_{3}.links_untransposed.track'.format(paths['CurrentTopDir'], classif, i, '_'.join(G))
                    # If no links for this cluster, don't draw a Circos 
                    # plot for the elements without scaffolds.
                    if os.stat(links_untransposedfl).st_size == 0:
                        continue
                    if (os.path.isfile(highlights_ltrs_fl) 
                        and os.path.isfile(tilefl) 
                        and os.path.isfile(links_untransposedfl) 
                        and os.path.isfile(seqfl)):
                        # Files exist. copy and run Circos.
                        circosdir = '{0}/circos.{1}.cluster_{2}.geneconv_{3}.justelements'.format(
                               paths['CurrentTopDir'], classif, i, '_'.join(G))
                        if not os.path.exists(circosdir):
                            # copy circos conf files and dir structure
                            copytree('{0}/circos'.format(
                               paths['scriptsDir']), circosdir)
                        totallengthsLTRs = 0
                        newseqfl = '{0}/data/{1}'.format(circosdir, 
                            '{0}.seq.track'.format('.'.join(
                                       tilefl.split('/')[-1].split('.')[:-2])))
                        newhlfl = '{0}/data/{1}'.format(circosdir, 
                            '{0}.LTR_highlights.track'.format('.'.join(
                           highlights_ltrs_fl.split('/')[-1].split('.')[:-2])))
                        # Copy hl fl to circos dir
                        copyfile(highlights_ltrs_fl, newhlfl)
                        if os.path.isfile(newseqfl):
                            os.remove(newseqfl)
                        # Convert tile file to ideogram track
                        with open(tilefl, 'r') as inFl:
                            with open(newseqfl, 'w') as outFl:
                                for line in inFl:
                                    scaf, start, end, val = line.strip().split()
                                    length = int(end) - int(start) + 1
                                    totallengthsLTRs += length
                                    color = 'dorange'
                                    outline = 'chr - {0} {1} 0 {2} {3}\n'.format(
                                        'LTR_retrotransposon{0}'.format(val), 
                                        val, length, color)
                                    outFl.write(outline)
                        newlinksuntransposedfl = '{0}/data/{1}'.format(circosdir, 
                                           links_untransposedfl.split('/')[-1])
                        copyfile(links_untransposedfl, newlinksuntransposedfl)
                        conffl = '{0}/etc/circos.conf'.format(circosdir)
                        circos_conf_str = """<<include colors_fonts_patterns.conf>>
<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

karyotype   = data/{0}

chromosomes_units           = 1000000

<highlights>
 <highlight>
 file       = data/{1}
 ideogram   = yes
 color = black
 </highlight>
</highlights>


<links>
radius = 0.999r
#radius = 1r
crest  = 1
ribbon           = yes
flat             = yes
#stroke_color     = vdgrey
stroke_thickness = 2
color            = grey_a3

bezier_radius        = 0r
bezier_radius_purity = 0.5

<link>
file       = data/{2}
</link>

</links>

<<include etc/housekeeping.conf>>
data_out_of_range* = trim""".format(newseqfl.split('/')[-1], 
                                    newhlfl.split('/')[-1], 
                                    newlinksuntransposedfl.split('/')[-1])
                        with open(conffl, 'w') as outFl:
                            outFl.write(circos_conf_str)
                        confbasename = conffl.split('/')[-1]
                        imagesize = totallengthsLTRs/10
                        if imagesize > 8000:
                            imagesize = 8000
                        if imagesize < 1000:
                            imagesize = 1000
                        conffl = '{0}/etc/image.generic.conf'.format(circosdir)
                        circosimageconfstr = """dir   = . 
file  = circos.png
png   = yes
svg   = yes

# radius of inscribed circle in image
radius = {0}

# by default angle=0 is at 3 o'clock position
angle_offset      = -96

#angle_orientation = counterclockwise

auto_alpha_colors = yes
auto_alpha_steps  = 5""".format(imagesize)
                        with open(conffl, 'w') as outFl:
                            outFl.write(circosimageconfstr)
                        circos_call = [executables['perl'], executables['circos']]
                        circoscallsElmts.append([circosdir, 
                                                circos_call, 
                                                classif, 
                                                i, 
                                                '{0}/plots.elements'.format(
                                                   paths['CurrentTopDir']), 
                                                G])
                        # Clean up
                        if not KEEP_UNUSED_FILES:
                            if os.path.isfile(tilefl):
                                os.remove(tilefl)
                            if os.path.isfile(linksfl):
                                os.remove(linksfl)
                            if os.path.isfile(seqfl):
                                os.remove(seqfl)
                            if os.path.isfile(textfl):
                                os.remove(textfl)
                            if os.path.isfile(highlights_ltrs_fl):
                                os.remove(highlights_ltrs_fl)
                            if os.path.isfile(links_untransposedfl):
                                os.remove(links_untransposedfl)
                            if os.path.isfile(GFFoutPth):
                                os.remove(GFFoutPth)
        MakeDir('plotdir', '{0}/plots.scaffolds'.format(
                                                       paths['CurrentTopDir']))
        paths['Circos_output_dir_scaffolds_{0}'.format(
                                               ClustMethod)] = paths['plotdir']
        if not checkStatusFl('Circos_output_dir_scaffolds_{0}'.format(
                                                                 ClustMethod)):
            chunk_size = ceil(len(circoscallsScafs) / procs)
            if not circoscallsScafs == []:
                with Pool(processes = procs) as p:
                    p.map(circosMultiprocessing, circoscallsScafs, 
                                                          chunksize = chunk_size)
                p.join()
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format(
                    'Circos_output_dir_scaffolds_{0}'.format(ClustMethod), 
                    paths['plotdir']))
        append2logfile(paths['output_top_dir'], mainlogfile, 'Made Circos plots.')
        MakeDir('plotdir', '{0}/plots.elements'.format(paths['CurrentTopDir']))
        paths['Circos_output_dir_elements_{0}'.format(ClustMethod)] = paths['plotdir']
        if not checkStatusFl('Circos_output_dir_elements_{0}'.format(ClustMethod)):
            if not circoscallsElmts == []:
                chunk_size = ceil(len(circoscallsElmts)/procs)
                with Pool(processes = procs) as p:
                    p.map(circosMultiprocessing, circoscallsElmts, 
                                                          chunksize = chunk_size)
                p.join()
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format(
                          'Circos_output_dir_elements_{0}'.format(ClustMethod),
                           paths['plotdir']))
        append2logfile(paths['output_top_dir'], mainlogfile, 'Made Circos plots.')



def summarizeClusters(I = 6, 
                      clustering_method = 'WickerFam', 
                      WickerParams = {'pId':80,'percAln':80,'minLen':80}):
    """Creates files related to each cluster.

    1. Creates a file for either MCL or WickerFam clustering at 
       PhyLTR.output/<clustMethod>/<settings>/Clusters/<clustMethod>_<settings>.summary.tab
    2. Creates GFF files for each cluster in 
       PhyLTR.output/GFF_output/<superfamily>.<clustMethod>.<settings>/<superfamily>.<clustMethod>.<settings>.cluster_<i>.gff
    3. Creates a text file tsv, 3-col: elementName classification 
       cluster
    """
    if clustering_method == 'WickerFam':
        wicker_top_dir = paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])]
        settings = '{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        ClusterSummaryFl = 'WickerClusterSummary_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        ClusterMembershipFl = 'WickerClusterMembership_{0}_pId_{1}_percAln_{2}_minLen'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        paths[ClusterSummaryFl] = '{0}/Clusters/Wicker_{1}_pId_{2}_percAln_{3}_minLen.summary.tab'.format(
                                                        wicker_top_dir, 
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        paths[ClusterMembershipFl] = '{0}/Clusters/Wicker_{1}_pId_{2}_percAln_{3}_minLen.membership.tab'.format(
                                                        wicker_top_dir, 
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'])
        if not checkStatusFl(ClusterSummaryFl) or not checkStatusFl(
                                                          ClusterMembershipFl):
            with open(paths[ClusterSummaryFl], 'w') as outFl:
                outFl.write('classification\tcluster\tsize\n')
            with open(paths[ClusterMembershipFl], 'w') as outFl:
                outFl.write('element\tclassification\tcluster\n')
            for classif in clusters_by_classif:
                with open(paths['WickerFamDir_{0}_pId_{1}_percAln_{2}_minLen_{3}'.format(
                                                        WickerParams['pId'], 
                                                        WickerParams['percAln'], 
                                                        WickerParams['minLen'], 
                                                        classif)], 'r') as inFl:
                    c = 0
                    for line in inFl:
                        clust_members = set([el.lstrip('LTR_retrotransposon') 
                                               for el in line.strip().split()])
                        clust_size = len(clust_members)
                        # Write GFF3 for cluster
                        MakeDir('{0}_cluster_GFF3s'.format(classif), 
                                            '{0}/{1}.{2}.{3}'.format(
                                                  paths['GFFByClassification'], 
                                                  classif, 
                                                  clustering_method, 
                                                  settings))
                        with open(paths['CurrentGFF']) as gffFl:
                            with open('{0}/{1}.{2}.{3}.cluster_{4}.gff'.format(
                                                   paths['{0}_cluster_GFF3s'.format(
                                                                     classif)], 
                                                  classif, 
                                                  clustering_method, 
                                                  settings, 
                                                             c), 'w') as outFl:
                                for line in gffFl:
                                    if line.startswith('#'):
                                        continue
                                    gffLine = GFF3_line(line)
                                    if 'Parent' in gffLine.attributes:
                                        if ('LTR_retrotransposon' in 
                                                 gffLine.attributes['Parent']):
                                            elNum = gffLine.attributes[
                                                       'Parent'].lstrip(
                                                         'LTR_retrotransposon')
                                        elif ('repeat_region' in 
                                                 gffLine.attributes['Parent']):
                                            elNum = gffLine.attributes['Parent'].lstrip(
                                                               'repeat_region')
                                        else:
                                            sys.exit(('Problem with summarizeClusters() '
                                                'Parent attribute for \n{0}').format(
                                                                 str(gffLine)))
                                    elif 'ID' in gffLine.attributes:
                                        if 'LTR_retrotransposon' in gffLine.attributes['ID']:
                                            elNum = gffLine.attributes['ID'].lstrip(
                                                         'LTR_retrotransposon')
                                        elif 'repeat_region' in gffLine.attributes['ID']:
                                            elNum = gffLine.attributes['ID'].lstrip(
                                                               'repeat_region')
                                        else:
                                            sys.exit(('Problem with summarizeClusters() '
                                                'ID attribute for \n{0}').format(
                                                                 str(gffLine)))
                                    else:
                                        sys.exit(('sumarizeClusters(): No Parent '
                                            'or ID attribute for gff line: {0}').format(
                                                                 str(gffLine)))
                                    if elNum in clust_members:
                                        if gffLine.type == 'repeat_region':
                                            outFl.write('###\n')
                                        outFl.write(line)
                        # Write line to summary file
                        with open(paths[ClusterSummaryFl], 'a') as outFl:
                            outFl.write('{0}\t{1}\t{2}\n'.format(classif, 
                                                                 c, 
                                                                 clust_size))
                        with open(paths[ClusterMembershipFl], 'a') as outFl:
                            for el in clust_members:
                                outFl.write('{0}\t{1}\t{2}\n'.format(el, 
                                                                     classif, 
                                                                     c))
                        c += 1
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format(ClusterSummaryFl, 
                                                      paths[ClusterSummaryFl]))
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format(ClusterMembershipFl, 
                                                   paths[ClusterMembershipFl]))
    elif clustering_method == 'MCL':
        mcl_top_dir = paths['MCL_I{0}'.format(I)]
        settings = 'I{0}'.format(I)
        ClusterSummaryFl = 'MCL_ClusterSummary_I{0}'.format(I)
        ClusterMembershipFl = 'MCL_ClusterMembership_I{0}'.format(I)
        paths[ClusterSummaryFl] = '{0}/Clusters/MCL_I{1}.summary.tab'.format(
                                                                mcl_top_dir, I)
        paths[ClusterMembershipFl] = '{0}/Clusters/MCL_I{1}.membership.tab'.format(
                                                                mcl_top_dir, I)
        if not checkStatusFl(ClusterSummaryFl) or not checkStatusFl(
                                                          ClusterMembershipFl):
            with open(paths[ClusterSummaryFl], 'w') as outFl:
                outFl.write('superfamily\tcluster\tsize\n')
            with open(paths[ClusterMembershipFl], 'w') as outFl:
                outFl.write('element\tclassification\tcluster\n')
            for classif in clusters_by_classif:
                with open(paths['MCL_{0}_I{1}'.format(
                                                    classif, I)], 'r') as inFl:
                    c = 0
                    for line in inFl:
                        clust_members = set([el.lstrip('LTR_retrotransposon') 
                                               for el in line.strip().split()])
                        clust_size = len(clust_members)
                        # Write GFF3 for cluster
                        MakeDir('{0}_cluster_GFF3s'.format(classif), 
                                        '{0}/{1}.{2}.{3}'.format(
                                                  paths['GFFByClassification'], 
                                                  classif, 
                                                  clustering_method, 
                                                  settings))
                        with open(paths['CurrentGFF']) as gffFl:
                            with open('{0}/{1}.{2}.{3}.cluster_{4}.gff'.format(
                                                  paths['{0}_cluster_GFF3s'.format(
                                                                     classif)], 
                                                  classif, 
                                                  clustering_method, 
                                                  settings, 
                                                  c), 'w') as outFl:
                                for line in gffFl:
                                    if line.startswith('#'):
                                        continue
                                    gffLine = GFF3_line(line)
                                    if 'Parent' in gffLine.attributes:
                                        if ('LTR_retrotransposon' in 
                                                 gffLine.attributes['Parent']):
                                            elNum = gffLine.attributes['Parent'].lstrip(
                                                         'LTR_retrotransposon')
                                        elif ('repeat_region' in 
                                                 gffLine.attributes['Parent']):
                                            elNum = gffLine.attributes[
                                                'Parent'].lstrip(
                                                               'repeat_region')
                                        else:
                                            sys.exit(('Problem with summarizeClusters() '
                                                      'Parent attribute for \n{0}').format(
                                                                 str(gffLine)))
                                    elif 'ID' in gffLine.attributes:
                                        if 'LTR_retrotransposon' in gffLine.attributes['ID']:
                                            elNum = gffLine.attributes['ID'].lstrip(
                                                         'LTR_retrotransposon')
                                        elif 'repeat_region' in gffLine.attributes['ID']:
                                            elNum = gffLine.attributes['ID'].lstrip(
                                                               'repeat_region')
                                        else:
                                            sys.exit(('PROBLEM WITH summarizeClusters() '
                                                'ID attribute for \n{0}').format(
                                                                 str(gffLine)))
                                    else:
                                        sys.exit(('sumarizeClusters(): No Parent '
                                            'or ID attribute for gff line: {0}').format(
                                                                 str(gffLine)))
                                    if elNum in clust_members:
                                        if gffLine.type == 'repeat_region':
                                            outFl.write('###\n')
                                        outFl.write(line)
                        # Write line to summary file
                        with open(paths[ClusterSummaryFl], 'a') as outFl:
                            outFl.write('{0}\t{1}\t{2}\n'.format(classif, 
                                                                 c, 
                                                                 clust_size))
                        with open(paths[ClusterMembershipFl], 'a') as outFl:
                            for el in clust_members:
                                outFl.write('{0}\t{1}\t{2}\n'.format(el, 
                                                                     classif, 
                                                                     c))
                        c += 1
            with open('{0}/status'.format(
                              paths['output_top_dir']), 'a') as statusFlAppend:
                statusFlAppend.write('{0}\t{1}\n'.format(ClusterSummaryFl, 
                                                      paths[ClusterSummaryFl]))
            with open('{0}/status'.format(
                               paths['output_top_dir']), 'a') as statusflappend:
                statusflappend.write('{0}\t{1}\n'.format(ClusterMembershipFl, 
                                                   paths[ClusterMembershipFl]))


def shortHelp():

    print('''
      Usage:
      ------------
          phyltr [options] --fasta <path> [options]
      ''', file=sys.stderr)


def help2():
    print('''
                                     Default
General--------|--------------------------------------------
               | -h or --h    
               | -help or --help    
               | --logfile            log.txt
               | --procs              1
               | --output_dir         phyltr.output
               | --keep_files    
               | --min_clust_size     7
               | --nosmalls    
Input----------|--------------------------------------------
               | --fasta    
LTRharvest-----|--------------------------------------------
               | --ltrharvest    
               | --del                -3
               | --ins                -3
               | --mis                -2
               | --mat                 2
               | --xdrop               5
               | --minlenltr           100
               | --maxlenltr           1000
               | --mindistltr          1000
               | --maxdistltr          15000
               | --similar             0.0
               | --vic                 60
               | --mintsd              4
               | --maxtsd              20
LTRdigest------|--------------------------------------------
               | --ltrdigest    
               | --ltrdigest_hmms PhyLTR/RepeatDatabases/LTRdigest_HMMs/hmms
ORFs-----------|--------------------------------------------
               | --findORFs    
               | --min_orf_len         300
Classify-------|--------------------------------------------
               | --no_classification
               | --no_dfam    
               | --nhmmer_reporting_evalue    10
               | --nhmmer_inclusion_evalue    1e-5
               | --no_repbase    
               | --repbase_tblastx_evalue     1e-5
               | --keep_no_classification    
               | --keep_conflicting_classifications    
Cluster--------|--------------------------------------------
               | --wicker    
               | --wicker_pId         80
               | --wicker_pAln        80
               | --wicker_minLen      80
               | --wicker_no_ltrs    
               | --wicker_no_internals    
               | --mcl    
               | --I                  6
Solo LTR-------|--------------------------------------------
               | --soloLTRsearch    
               | --soloLTRminPid      80.0
               | --soloLTRminLen      80.0
               | --soloLTRmaxEvalue   1e-3
Alignment------|--------------------------------------------
               | --mafft_align_region           entire
               | --maxiterate_small_clusters    20
               | --maxiterate_medium_clusters   3
               | --mafft_smallAln_maxclustsize  50
               | --mafft_mediumAln_maxclustsize 500
               | --mafft_largeAln_maxclustsize  1000
LTR divergence-|--------------------------------------------
               |  (All alignment settings)
               | --modeltest                    off
               | --model                        hky85
               | --geneconvltrs    
               | --geneconvclusters    
               | --geneconv_g                   g0,g1,g2
               | --remove_GC_from_modeltest_aln    
               | --ltrdivergence    
               | --circos    
Phylo----------|--------------------------------------------
               | All alignment settings    
               | --phylo    
               | --LTT    
               | --rmhomoflank    
               | --bpflank              500
               | --flank_evalue         1e-5
               | --flank_pId            70.0
               | --flank_plencutoff     70.0
               | --auto_outgroup    
               | --bootstrap_reps       100
               | --convert_to_ultrametric    
'''.format('{0}/RepeatDatabases/LTRdigest_HMMs/hmms'.format(paths['selfDir']), 
                                                              file=sys.stderr))


def help():
    print('''

      Usage:
      ------------
      phyltr [options] -fasta <input.fasta> [options]
      
      Description:
      ------------
      The main options would be:
      phyltr --fasta fasta.fa \\
           --procs 40 \\
         --ltrharvest \\
         --ltrdigest \\
         --mcl \\
         --wicker \\
         --geneconvltrs \\
         --geneconvclusters \\
         --ltrdivergence \\
         --LTT

      ------------------------------
      Global Options:
      ------------------------------
      -f | --fasta       <path>    Sequences to analyze. Mandatory.
      -p | --procs       <int>     Number of processors (default 1)
      -o | --output_dir  <path>    Output directory name. Default is 
                                    PhyLTR.output
      --logfile          <path>    Path to where log file is written 
                                    Default <output_dir>/log.txt
      -h                           Output abbreviated help information
      -help                        Output help information
      
      -------------------------
      Program-specific Options:
      -------------------------

      LTRharvest
      ----------
      -lh | --ltrharvest          Run LTRharvest on file given by 
                                   --fasta (default ON)
        --minlenltr  <int>        minimum length allowed for LTRs for 
                                   element calling (default 100 bp)
        --maxlenltr  <int>        maximum length allowed for LTRs for 
                                   element caling (default 1000 bp)
        --mindistltr <int>        minimum distance allowed between LTRs 
                                   for element calling (default 1000 bp)
        --maxdistltr <int>        maximum distance allowed between LTRs 
                                   for element calling (default 15000 bp)
        --similar    <int|float>  minimum % similarity for LTR calling 
                                   (default 0.0)
        --vic        <int>        No. of nucleotides to left and right 
                                   to search for TSDs (default 60)
        --mintsd     <int>        minimum length allowed for TSDs (use 
                                   with --maxtsd) (default 4)
        --maxtsd     <int>        maximum length allowed for TSDs (use 
                                   with --mintsd) (default 20)
        --xdrop      <int>        xdropbelow score for extension-alignment 
                                   (default 5)
        --mat        <int>        matchscore score for extension-alignment 
                                   (default 2)
        --mis        <int>        mismatchscore score for extension-alignment 
                                   (default -2)
        --ins        <int>        insertionscore for extension-alignment 
                                   (default -3)
        --del        <int>        deletionscore for extension-alignment 
                                   (default -3)

      ORF annotation
      --------------
      --findORFs                  Turns on ORF annotation
      --min_orf_len   <int>       (default 300)

      LTRdigest
      ---------
      -ld | --ltrdigest           Run LTRdigest on file given by --fasta 
                                   and --ltrharvest results (GFF) (default ON)
      --ltrdigest_hmms   <path>   Path to a file with one or more protein 
                                   profile HMMs for LTRdigest (HMMER)
                                   (default {0})

      Classification of LTR RTs to superfamily using homology to annotated 
       sequences in Repbase and/or Dfam
      --------------------------------------------------------------------
      --no_classification                 Do not classify LTR-Rs with 
                                           Repbase or Dfam (on by default)
      --no_dfam                           Do not run hmmsearch on Dfam 
                                           database (on by default)
      --no_repbase                        Do not run tblastx on Repbase 
                                           database (on by default)
      --keep_conflicting_classifications  If an element has two 
                                           classifications that disagree 
                                            (i.e. Repbase and Dfam) and 
                                            one classification is an LTR 
                                            RT hit and one is a non-LTR 
                                            RT hit, and this flag is set, 
                                            the element will be kept.
                                            (default OFF; the element is 
                                            discarded as a false positive)
      --keep_no_classification            If an element does not have 
                                           evidence of homology to any 
                                           sequence in Repbase or Dfam, 
                                           keep it. (default OFF; elements 
                                           without homology to LTR RT in 
                                           a database are discarded as 
                                           false positives)
      --repbase_tblastx_evalue h <int|float>  Max allowed E-value for 
                                               tblastx hits (default 1e-5)
      --nhmmer_reporting_evalue  <int|float>  (default 10)
      --nhmmer_inclusion_evalue  <int|float>  (default 1e-2)

      Clustering
      -----------
      --min_clust_size      Minimum allowed cluster size. Clusters with 
                             < minclustsize elements get assembled together
                             but no model testing, gene conversion, or 
                             outgroup/LTT analysis is done.
      --mcl                 Cluster using MCL
      --I                   Inflation/granularity parameter for MCL 
                             (default 6)
      --wicker              Cluster using '80-80-80' rule, or custom values 
                             specified below.
      --wicker_pId          Minimum % ID in pairwise alignment between 
                             any two elements (default 80)
      --wicker_minLen       Minimum alignment length to considered. 
                             (default 80)
      --wicker_pAln         Minimum percentage of whole sequence (LTRs 
                             or internal region) for alignment to be 
                             considered. (default 80)
      --wicker_no_internals (default OFF)
      --wicker_no_ltrs      (default OFF)


      MAFFT (for cluster, not LTR alignment): Default for very large 
       clusters (clust_size > 200)
      --------------------------------------------------------------------

      NOTE: MAFFT uses a lot of RAM for large clusters. It repeatedly 
       failed on tests of ~2.7k seqs of length >5kb using 256Gb RAM. 
      The MAFFT algorthim FFT-NS-2 is used for small and medium clusters 
       with the user-specified --maxiterate option (see below) and 
       FFT_NS-1 for large clusters, which is very inaccurate.

      --mafft_align_region           <str> Either 'internal' or 'entire', 
                                            specifying the portions of 
                                            each element to use for 
                                            MAFFT alignments.
                                            Internal = the sequence region 
                                            between the two LTRs. 
                                            (default entire)
      --maxiterate_small_clusters    <int> Max number of iterations for 
                                            MAFFT algorithm for clusters 
                                            with < --min_clustsize_for_faster_aln 
                                            elements. 1 is fastest; greater 
                                            numbers will improve the 
                                            alignment (default 30).
      --maxiterate_medium_clusters   <int> Number of iterations for MAFFT 
                                            algorithm for large clusters 
                                            (--min_clustsize_for_faster_aln < 
                                            clust_size) (default 3)
      --mafft_smallAln_maxclustsize  <int> Clusters this size and smaller 
                                            will be aligned using the 
                                            settings from 
                                            --maxiterate_small_clusters. 
                                            (default 50)
      --mafft_mediumAln_maxclustsize <int> Clusters this size and smaller 
                                            but larger than 
                                            --mafft_smallAln_maxclustsize 
                                            will be aligned using the 
                                            settings from
      --maxiterate_medium_clusters.        (default 500)
      --mafft_largeAln_maxclustsize  <int> Clusters this size and smaller 
                                            but larger than 
                                            --mafft_mediumAln_maxclustsize 
                                            will be aligned using the 
                                            MAFFT algorithm FFT-NS-1 
                                            (--retree 1). Clusters larger 
                                            than this size will not be 
                                            aligned and the downstream 
                                            analyses that require an
                                            alignment will not be performed. 
                                            Those analyses are LTR 
                                            divergence, gene conversion,
                                            and phylogenetic analyses.
                                            (default 1000)
      --nosmalls                           Do not combine and assemble 
                                            clusters smaller than 
                                            --min_clust_size (see clustering 
                                            options)

      GENECONV
      --------------------------------------------------------------------
      --geneconv_g <comma-sep-list> A comma-separated list of the 
                                     values for the 'g' parameter 
                                     of GENECONV. Valid values are 
                                     g0, g1, and g2.
                                     (default g0,g1,g2)
      --geneconvltrs
      --geneconvclusters
      --circos                Make Circos plots for each cluster showing 
                               GENECONV results.

      LTR divergence estimation
      -------------------------
      --modeltest       Find best-supported model of nucleotide substitution and use it to estimate intra-element LTR divergences
                        Off by default.
      --model  <str>    Not implemented. Specify model to be used for 
                         estimating intra-element LTR divergences: One of:
                         jc, f81, tajnei, k2p, hky85, k3p, tamnei, gtr, 
                         logdet, upholt, neili
                        See PAUP* manual for explanation. Default hky85.
      --ltrdivergence    Estimate substitutions per site between LTRs 
                         for each element using best supported model 
                         from model test or --model (default hky85)
      --remove_GC_from_modeltest_aln  Remove elements with suspected 
                                       intra-cluster inter-element gene 
                                       conversion tracts.

      Solo LTR search
      -------------------
      --soloLTRsearch           Turn on solo LTR search.
      --soloLTRminPid    <num>  Minimum percent identity for inclusion of 
                                 a solo LTR in a cluster (80.0)
      --soloLTRminLen    <num>  Minimum percent of LTR length participating 
                                 in alignment for inclusion of LTR in a 
                                 cluster (80.0)
      --soloLTRmaxEvalue <num>  Maximum evalue allowed for blastn

      Finding pairs of elements within clusters that have homologous flanking regions
      -------------------------------------------------------------------------------
      --rmhomoflank              Remove one of each pair of elements within 
                                  each alignment (and therefore, each tree).
                                  (default OFF; Fixed ON when using --LTT)
      --bpflank      <int>       Number of bases on either side of each 
                                  element to search for homology. (default 
                                  500 bp)
      --flank_evalue <int|float> E-value ceiling for considering blastn 
                                  hits as evidence of homology. 
                                  (default 1e-5)
      --flank_pId    <int|float> Minimum percent identity in blastn 
                                  alignment to consider hit as evidence 
                                  of homology. (default 70)
      --flank_plencutoff <int|float> Minimum percentage of flanking region 
                                      required to participate in alignment 
                                      to consider blastn hit as evidence 
                                      of homology. (default 70)
      Phylogenetic analysis
      ---------------------
      --phylo                
      --nosmalls        Do not combine and perform phylogentic analyses 
                         on clusters smaller than --min_clust_size.
      --LTT             Turns on --rmhomoflank, --convert_to_ultrametric, 
                         and --auto_outgroup. Produces ultrametric
                        phylogenies which are ready for use with external 
                         R scripts for generating LTT plots and other
                        phylogenetic analyses.
      --bootstrap_reps         <int>  Number of replicates to generate 
                                       for bootstrapping (default 100)
      --convert_to_ultrametric        Convert trees to ultrametric using 
                                       PATHd8. (default OFF; ON when using 
                                       --LTT)
      --auto_outgroup                 Pick an outgroup automatically:
'''.format('{0}/RepeatDatabases/LTRdigest_HMMs/hmms'.format(
                                   paths['selfDir']), file=sys.stderr))


if __name__ == '__main__':
    args=sys.argv
    # get executable paths from CONFIG file, which should be in the 
    # same directory as this script
    executables = {}
    commentPattern = re.compile('#.*$')
    # CONFIG file specifies paths to executables
    with open('{0}/CONFIG'.format(os.path.dirname(os.path.realpath(
                                                   __file__)))) as config_file:
        paths = [re.sub(commentPattern, '', line) for line in 
                                        config_file.read().strip().split('\n')]
        for path in paths:
            if not path == '':
                p = path.split('=')
                executables[p[0]] = p[1]
    filenames = {}
    paths = {}
    paths_toClean = {}
    params = {}
    paths['selfDir'] = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    paths['scriptsDir'] = '{0}/scripts'.format(paths['selfDir'])
    if '-h' in args or '--h' in args:
        help2()
        sys.exit(0)

    if '-help' in args or '--help' in args:
        help()
        sys.exit(0)
    if len(args) < 3:
        shortHelp()
        sys.exit(0)

    if '--fasta' in args:
        paths['inputFasta'] = args[args.index('--fasta') + 1]
    elif '-f' in args:
        paths['inputFasta'] = args[args.index('-f') + 1]
    else:
        help()
        print('''

            Must specify input fasta with -f or --fasta

            ''', file=sys.stderr)
        sys.exit(0)
    # If any of these flags are specified, then only up to those 
    # processes are performed, otherwise they are all run.
    main_flags = ('--ltrharvest',
                  '--ltrdigest',
                  '--classify',
                  '--wicker',
                  '--mcl',
                  '--geneconvclusters',
                  '--circos',
                  '--sololtrsearch',
                  '--geneconvltrs',
                  '--ltrdivergence',
                  '--phylo',
                  '--LTT')
    DEFAULT = True
    for flag in main_flags:
        if flag in args:
            DEFAULT = False
            break
    # Set default flags
    if DEFAULT:
        print('Running default settings (all analyses)', file=sys.stderr)
        LTRHARVEST = True
        LTRDIGEST = True
        FINDORFS = True
        WICKER = True
        USEMCL = True
        LTRDIVERGENCE = True
        GENECONVLTRS = True
        GENECONVCLUSTERS = True
        CIRCOS = True
        PHYLO = True
        AUTO_OUTGROUP = True
        RMHOMOFLANK = True
        LTT = True
        ULTRAMETRIC = True
        SOLOLTR = True
    if '--logfile' in args:
        mainlogfile = args[args.index('--logfile')+1]
    else:
        mainlogfile = 'log.txt'
    filenames['inputFasta'] = paths['inputFasta'].split('/')[-1]
    KEEP_UNUSED_FILES = True
    # Max number of processors for parallel flows
    procs = 1
    if '--procs' in args:
        procs = int(args[args.index('--procs') + 1])
    if '-p' in args:
        procs = int(args[args.index('-p') + 1])
    if '--output_dir' in args:
        MakeDir('output_top_dir', args[args.index('--output_dir') + 1])
    if '-o' in args:
        MakeDir('output_top_dir', args[args.index('-o') + 1])
    else:
        MakeDir('output_top_dir', 'PhyLTR.output')
    MakeDir('FastaOutputDir', '{0}/FASTA_output'.format(
                                              paths['output_top_dir']))
    MakeDir('GFFOutputDir', '{0}/GFF_output'.format(
                                              paths['output_top_dir']))
    # Do not align and infer phylogenies for small clusters (size(cluster) 
    # < --min_clust_size)
    if '--nosmalls' in args:
        SMALLS = False
    else:
        SMALLS = True
    # 1. LTRharvest
    if '--ltrharvest' in args or '-lh' in args or DEFAULT:
        LTRHARVEST = True
    else:
        LTRHARVEST = False
    # LTRharvest controls
    if '--del' in args:
        ltrharvest_del = int(args[args.index('--del')+1])
    else:
        ltrharvest_del = -3
    if '--ins' in args:
        ltrharvest_ins = int(args[args.index('--ins')+1])
    else:
        ltrharvest_ins = -3
    if '--mis' in args:
        ltrharvest_mis = int(args[args.index('--mis')+1])
    else:
        ltrharvest_mis = -2
    if '--mat' in args:
        ltrharvest_mat = int(args[args.index('--mat')+1])
    else:
        ltrharvest_mat = 2
    if '--xdrop' in args:
        ltrharvest_xdrop = int(args[args.index('--xdrop')+1])
    else:
        ltrharvest_xdrop = 5 
    if '--minlenltr' in args:
        ltrharvest_minlenltr = int(args[args.index('--minlenltr')+1])
    else:
        ltrharvest_minlenltr = 100
    if '--maxlenltr' in args:
        ltrharvest_maxlenltr = int(args[args.index('--maxlenltr')+1])
    else:
        ltrharvest_maxlenltr = 1000
    if '--mindistltr' in args:
        ltrharvest_mindistltr = int(args[args.index('--mindistltr')+1])
    else:
        # Terminal-repeat retrotransposons in miniature (TRIM) are 
        # involved in restructuring plant genomes. Witte et al. 2001,
        # PNAS. TRIM elements have internal regions of only 100-300 bp
        ltrharvest_mindistltr = 100
    if '--maxdistltr' in args:
        ltrharvest_maxdistltr = int(args[args.index('--maxdistltr')+1])
    else:
        ltrharvest_maxdistltr = 15000
    if '--similar' in args:
        ltrharvest_similar = float(args[args.index('--similar')+1])
    else:
        ltrharvest_similar = 0.0
    if '--vic' in args:
        ltrharvest_vic = int(args[args.index('--vic')+1])
    else:
        ltrharvest_vic = 60
    if '--mintsd' in args:
        ltrharvest_mintsd = int(args[args.index('--mintsd')+1])
    else:
        ltrharvest_mintsd = 4
    if '--maxtsd' in args:
        ltrharvest_maxtsd = int(args[args.index('--maxtsd')+1])
    else:
        ltrharvest_maxtsd = 20


    # 2. LTRdigest
    if '--ltrdigest' in args or '-ld' in args or DEFAULT:
        LTRDIGEST = True
    else:
        LTRDIGEST = False
    # Check for user-supplied location of HMMs for LTRdigest, set 
    # default otherwise (comes with package)
    if '--ltrdigest_hmms' in args:
        paths['LTRdigestHMMs'] = args[args.index('--ltrdigest_hmms') + 1]
    else:
        paths['LTRdigestHMMs'] = '{0}/RepeatDatabases/LTRdigest_HMMs/hmms'.format(
                                                             paths['selfDir'])
    # 3. ORF annotation
    if '--findORFs' in args or DEFAULT:
        FINDORFS = True
    else:
        FINDORFS = False
    if '--min_orf_len' in args:
        min_orf_len = int(args[args.index('--min_orf_len')+1])
    else:
        min_orf_len = 300
    # 4. Classification
    # Solves can't find blastdbs problem
    os.environ['BLASTDB'] = paths['FastaOutputDir']
    if '--no_classification':
        CLASSIFYDFAM = False
        CLASSIFYREPBASE = False
    else:
        CLASSIFYREPBASE = True
        CLASSIFYDFAM = True
    if '--no_dfam' in args:
        CLASSIFYDFAM = False
    else:
        CLASSIFYDFAM = True
    if '--no_repbase' in args:
        CLASSIFYREPBASE = False
    else:
        CLASSIFYREPBASE = True
    if '--keep_conflicting_classifications' in args:
        KEEPCONFLICTS=True
    else:
        KEEPCONFLICTS=False
    if '--keep_no_classification' in args:
        KEEPNOCLASSIFICATION=True
    else:
        KEEPNOCLASSIFICATION=False
    # Dfam control
    if '--nhmmer_reporting_evalue' in args:
        nhmmer_reporting_evalue = float(
                               args[args.index('--nhmmer_reporting_evalue')+1])
    else:
        nhmmer_reporting_evalue = 10
    # Dfam control
    if '--nhmmer_inclusion_evalue' in args:
        nhmmer_inclusion_evalue = float(
                               args[args.index('--nhmmer_inclusion_evalue')+1])
    # Repbase control
    else:
        nhmmer_inclusion_evalue = 1e-5
    if '--repbase_tblastx_evalue' in args:
        repbase_tblastx_evalue = float(
                                args[args.index('--repbase_tblastx_evalue')+1])
    else:
        repbase_tblastx_evalue = 1e-5
    # 5. Clustering
    if '--wicker' in args or DEFAULT:
        WICKER = True
    else:
        WICKER = False
    if '--wicker_pId' in args:
        wicker_pId = int(args[args.index('--wicker_pId')+1])
    else:
        wicker_pId = 80
    if '--wicker_pAln' in args:
        wicker_pAln = int(args[args.index('--wicker_pAln')+1])
    else:
        wicker_pAln = 80
    if '--wicker_minLen' in args:
        wicker_minLen = int(args[args.index('--wicker_minLen')+1])
    else:
        wicker_minLen = 80
    if '--wicker_no_ltrs' in args:
        wicker_use_ltrs = False
    else:
        wicker_use_ltrs = True
    if '--wicker_no_internals' in args:
        wicker_use_internal = False
    else:
        wicker_use_internal = True
    if '--mcl' in args or DEFAULT:
        USEMCL = True
    else:
        USEMCL = False
    if '--I' in args:
        MCL_I = args[args.index('--I')+1]
    else:
        MCL_I = '6'
    if '--min_clust_size' in args:
        MinClustSize = int(args[args.index('--min_clust_size')+1])
    else:
        MinClustSize = 7
    # 6. Model testing for LTR divergence
    if '--modeltest' in args:
        MODELTEST = True
    else:
        MODELTEST = False
    if '--remove_GC_from_modeltest_aln' in args:
        remove_GC_from_modeltest_aln = True
    else:
        remove_GC_from_modeltest_aln = False
    # 7. LTR divergence and model testing
    if '--ltrdivergence' in args or DEFAULT:
        LTRDIVERGENCE = True
        MODELTEST = True
    else:
        LTRDIVERGENCE = False
    if '--model' in args:
        model = args[args.index('--model')+1]
    else:
        model = 'hky85'
    # 8. Intra-element LTR gene conversion
    if '--geneconvltrs' in args or DEFAULT:
        GENECONVLTRS = True
    else:
        GENECONVLTRS = False
    # 9. Inter-element intra-cluster LTR gene conversion
    if '--geneconvclusters' in args or DEFAULT:
        GENECONVCLUSTERS = True
    else:
        GENECONVCLUSTERS = False
    if '--geneconv_g' in args:
        gcparams = args[args.index('--geneconv_g')+1].split(',')
        if 'g0' in gcparams:
            GENECONV_G0 = True
        else:
            GENECONV_G0 = False
        if 'g1' in gcparams:
            GENECONV_G1 = True
        else:
            GENECONV_G1 = False
        if 'g2' in gcparams:
            GENECONV_G2 = True
        else:
            GENECONV_G2 = False
    else:
        GENECONV_G0 = True
        GENECONV_G1 = True
        GENECONV_G2 = True
    # 10. Circos
    if '--circos' in args or DEFAULT:
        CIRCOS = True
    else:
        CIRCOS = False
    # 10. Circos
    if '--nocircos' in args:
        CIRCOS = False
    # 11. Phylogenetic inference on clusters
    if '--phylo' in args or DEFAULT:
        PHYLO = True
    else:
        PHYLO = False
    if '--bootstrap_reps' in args:
        bootstrap_reps = int(args[args.index('--bootstrap_reps')+1])
    else:
        bootstrap_reps = 100
    if '--rmhomoflank' in args:
        RMHOMOFLANK = True
    else:
        RMHOMOFLANK = False
    if '--convert_to_ultrametric' in args:
        ULTRAMETRIC = True
    else:
        ULTRAMETRIC = False
    if '--auto_outgroup' in args:
        AUTO_OUTGROUP = True
    else:
        AUTO_OUTGROUP = False
    if '--LTT' in args or DEFAULT:
        AUTO_OUTGROUP = True
        RMHOMOFLANK = True
        LTT = True
        ULTRAMETRIC = True
    else:
        AUTO_OUTGROUP = False
        RMHOMOFLANK = False
        LTT = False
        ULTRAMETRIC = False
    if '--bpflank' in args:
        bpflank = int(args[args.index('--bpflank')+1])
    else:
        bpflank = 500
    if '--flank_evalue' in args:
        flank_evalue = int(args[args.index('--flank_evalue')+1])
    else:
        flank_evalue = 1e-5
    if '--flank_pId' in args:
        flank_pId = float(args[args.index('--flank_pId')+1])
    else:
        flank_pId = 70.0
    if '--flank_plencutoff' in args:
        flank_plencutoff = float(args[args.index('--flank_plencutoff')+1])
    else:
        flank_plencutoff = 70.0
    # 12. Solo LTR search parameters
    if '--soloLTRsearch' in args or DEFAULT:
        SOLOLTR = True
    else:
        SOLOLTR = False
    if '--soloLTRminPid' in args:
        soloLTRminPid = str(float(args[args.index('--soloLTRminPid')+1]))
    else:
        soloLTRminPid = 80.0
    if '--soloLTRminLen' in args:
        soloLTRminLen = float(args[args.index('--soloLTRminLen')+1])
    else:
        soloLTRminLen = 80.0
    if '--soloLTRmaxEvalue' in args:
        soloLTRmaxEvalue = float(args[args.index('--soloLTRmaxEvalue')+1])
    else:
        soloLTRmaxEvalue = 1e-3
    # MAFFT parameters
    if '--mafft_align_region' in args:
        mafft_align_region = args[args.index('--mafft_align_region')+1]
        if not mafft_align_region == 'entire' or mafft_align_region == 'internal':
            sys.exit(('Error: --mafft_align_region must be either \'internal\' '
                        'or \'entire\''))
    else:
        mafft_align_region = 'entire' 
    if '--maxiterate_small_clusters' in args:
        mafft_smallAln_maxiterate = int(
                             args[args.index('--maxiterate_small_clusters')+1])
    else:
        mafft_smallAln_maxiterate = 20
    if '--maxiterate_medium_clusters' in args:
        mafft_mediumAln_maxiterate = int(
                            args[args.index('--maxiterate_medium_clusters')+1])
    else:
        mafft_mediumAln_maxiterate = 3
    if '--mafft_smallAln_maxclustsize' in args:
        mafft_smallAln_maxclustsize = int(
                           args[args.index('--mafft_smallAln_maxclustsize')+1])
    else:
        mafft_smallAln_maxclustsize = 50
    if '--mafft_mediumAln_maxclustsize' in args:
        mafft_mediumAln_maxclustsize = int(
                          args[args.index('--mafft_mediumAln_maxclustsize')+1])
    else:
        mafft_mediumAln_maxclustsize = 500
    if '--mafft_largeAln_maxclustsize' in args:
        mafft_largeAln_maxclustsize = int(
                           args[args.index('--mafft_largeAln_maxclustsize')+1])
    else:
        mafft_largeAln_maxclustsize = 1000
    paths['RepbaseDB'] = '{0}/RepeatDatabases/Repbase/Repbase_ERV_LTR.fasta'.format(
                                                              paths['selfDir'])
    paths['RepbaseTruePosLTRlist'] = '{0}/RepeatDatabases/Repbase/Repbase_ERV_LTR.list'.format(
                                                              paths['selfDir'])
    paths['RepbaseShortNames'] = '{0}/RepeatDatabases/Repbase/Repbase_ERV_LTR.SF'.format(
                                                              paths['selfDir'])
    paths['DfamDB'] = '{0}/RepeatDatabases/Dfam/Dfam_ERV_LTR.hmm'.format(
                                                              paths['selfDir'])
    paths['DfamTruePosLTRlist'] = '{0}/RepeatDatabases/Dfam/Dfam_ERV_LTR.list'.format(
                                                              paths['selfDir'])
    paths['DfamShortNames'] = '{0}/RepeatDatabases/Dfam/Dfam_ERV_LTR.SF'.format(
                                                              paths['selfDir'])
    LTR_SFs = ['Copia', 'Gypsy', 'ERV', 'Pao', 'BEL', 'Tas', 'Suzu', 'Sinbad', 
                                                                     'Unknown']
    # This path will have the path to the best GFF3 to use.
    paths['CurrentGFF'] = None
    try:
        statusFlRead = open('{0}/status'.format(paths['output_top_dir']), 'r')
    except:
        pass

    def write2summary(text):
        with open('{0}/summary'.format(paths['output_top_dir']), 'a') as summaryFl:
            summaryFl.write('{0}\n'.format(text))
    # Check for status file. if exists parse it and skip the sections 
    # that are done.
    try:
        statusContents = statusFlRead.read()
        if not statusContents == '':
            for pair in statusContents.strip().split('\n'):
                pair = pair.split('\t')
                paths[pair[0]] = pair[1]
    except:
        pass
    if 'LTRdigestClassifiedNoFP' in paths:
        paths['CurrentGFF'] = paths['LTRdigestClassifiedNoFP']
    elif 'GFFwithRepbaseClassification' in paths:
        paths['CurrentGFF'] = paths['GFFwithRepbaseClassification']
    elif 'GFFwithDfamClassification' in paths:
        paths['CurrentGFF'] = paths['GFFwithDfamClassification']
    elif 'WithORFsGFF' in paths:
        paths['CurrentGFF'] = paths['WithORFsGFF']
    elif 'LTRdigestGFF' in paths:
        paths['CurrentGFF'] = paths['LTRdigestGFF']
    elif 'LTRharvestGFF' in paths:
        paths['CurrentGFF'] = paths['LTRharvestGFF']
    sys.setrecursionlimit(50000) # For WickerFam() recursive subroutine
    # Predict LTR retrotransposons using structural criteria
    ltrharvest()
    # Identify parts of element internal regions with evidence of 
    # homology to LTR RT protein coding domains
    ltrdigest()
    if FINDORFS:
        AnnotateORFs(minLen=min_orf_len)
    # Extract LTR_retrotransposon sequences for classification using homology
    classify_by_homology(KEEPCONFLICTS = KEEPCONFLICTS, 
                         KEEPNOCLASSIFICATION = KEEPNOCLASSIFICATION, 
                         repbase_tblastx_evalue = repbase_tblastx_evalue, 
                         nhmmer_reporting_evalue = nhmmer_reporting_evalue, 
                         nhmmer_inclusion_evalue = nhmmer_inclusion_evalue) 
    clusters_by_classif = shortClassif()
    classifs_by_element = shortClassif(ElNames=True)
    classifs = set(list(clusters_by_classif.keys()))
    if WICKER:
        WickerFam(pId = wicker_pId, 
                  percAln = wicker_pAln, 
                  minLen = wicker_minLen, 
                  use_ltrs = wicker_use_ltrs, 
                  use_internal = wicker_use_internal)
        summarizeClusters(I = 6, 
                          clustering_method = 'WickerFam', 
                          WickerParams = {'pId':80,'percAln':80,'minLen':80})
        if GENECONVCLUSTERS or LTRDIVERGENCE:
            if not LTRDIVERGENCE:
                AutoAlign(I = None, 
                          part = mafft_align_region, 
                          rmgeneconv = False, 
                          minClustSize = MinClustSize, 
                          align = 'clusters', 
                          rmhomologflank = False, 
                          clustering_method = 'WickerFam', 
                          WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
                          auto_outgroup = False, 
                          bpflank = bpflank, 
                          combine_and_do_small_clusters = SMALLS, 
                          flank_pId = flank_pId, 
                          flank_evalue = flank_evalue, 
                          flank_plencutoff = flank_plencutoff, 
                          LTRSONLY = True)
            else:
                AutoAlign(I = None, 
                          part = mafft_align_region, 
                          rmgeneconv = False, 
                          minClustSize = MinClustSize, 
                          align = 'clusters', 
                          rmhomologflank = False, 
                          clustering_method = 'WickerFam', 
                          WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
                          auto_outgroup = False, 
                          bpflank = bpflank, 
                          combine_and_do_small_clusters = SMALLS, 
                          flank_pId = flank_pId, 
                          flank_evalue = flank_evalue, 
                          flank_plencutoff = flank_plencutoff, 
                          LTRSONLY = False)
                if MODELTEST:
                    modeltest(iters = 1, 
                              I = None, 
                              removegeneconv = remove_GC_from_modeltest_aln, 
                              part = mafft_align_region, 
                              clustering_method = 'WickerFam', 
                              WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
                              minClustSize = MinClustSize, 
                              bpflank = bpflank, 
                              combine_and_do_small_clusters = SMALLS)
            if GENECONV_G0:
                geneconvClusters(g = '/g0', 
                                 clust = None, 
                                 I = None, 
                                 minClustSize = MinClustSize, 
                                 clustering_method = 'WickerFam', 
                                 WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
                                 combine_and_do_small_clusters = SMALLS)
            if GENECONV_G1:
                geneconvClusters(g = '/g1', 
                                 clust = None, 
                                 I = None, 
                                 minClustSize = MinClustSize, 
                                 clustering_method = 'WickerFam', 
                                 WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
                                 combine_and_do_small_clusters = SMALLS)
            if GENECONV_G2:
                geneconvClusters(g = '/g2', 
                                 clust = None, 
                                 I = None, 
                                 minClustSize = MinClustSize, 
                                 clustering_method = 'WickerFam', 
                                 WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
                                 combine_and_do_small_clusters = SMALLS)
        if SOLOLTR:
            SoloLTRsearch(I = 6, 
                          clustering_method = 'WickerFam', 
                          WickerParams = {'pId':80,'percAln':80,'minLen':80})
        if CIRCOS:
            Circos(window = '1000000', 
                   plots = 'clusters', 
                   I = None, 
                   clustering_method = 'WickerFam', 
                   WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
        # Runs if need to use geneconvLTRs or estimate divergences
        align_ltrs(I = None, 
                   clustering_method = 'WickerFam', 
                   WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}) 
        if GENECONV_G0:
            geneconvLTRs(g = '/g0', 
                         I = None, 
                         clustering_method = 'WickerFam', 
                         WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
        if GENECONV_G1:
            geneconvLTRs(g = '/g1', 
                         I = None, 
                         clustering_method = 'WickerFam', 
                         WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
        if GENECONV_G2:
            geneconvLTRs(g = '/g2', 
                         I = None, 
                         clustering_method = 'WickerFam', 
                         WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen})
        ltr_divergence(I = None, 
                         clustering_method = 'WickerFam', 
                         WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
                         model = model)
        phylo(removegeneconv = False, 
              BOOTSTRAP = True, 
              I = None, 
              align = 'cluster', 
              removehomologouspair = RMHOMOFLANK, 
              part = mafft_align_region, 
              clustering_method = 'WickerFam', 
              WickerParams = {'pId':wicker_pId,'percAln':wicker_pAln,'minLen':wicker_minLen}, 
              auto_outgroup = AUTO_OUTGROUP, 
              bootstrap_reps = bootstrap_reps, 
              minClustSize = MinClustSize, 
              convert_to_ultrametric = ULTRAMETRIC, 
              bpflank = bpflank, 
              combine_and_do_small_clusters = SMALLS)
    if USEMCL:
        MCL(I = MCL_I, 
            minClustSize = MinClustSize, 
            CombineIfTooFew = False)    
        summarizeClusters(I = MCL_I, 
                          clustering_method = 'MCL', 
                          WickerParams = {'pId':80,'percAln':80,'minLen':80})
        if GENECONVCLUSTERS or LTRDIVERGENCE:
            if not LTRDIVERGENCE:
                AutoAlign(I = MCL_I, 
                          part = mafft_align_region, 
                          rmgeneconv = False, 
                          minClustSize = MinClustSize, 
                          align = 'clusters', 
                          rmhomologflank = False, 
                          clustering_method = 'MCL', 
                          WickerParams = None, 
                          auto_outgroup = False, 
                          bpflank = bpflank, 
                          combine_and_do_small_clusters = SMALLS, 
                          flank_pId = flank_pId, 
                          flank_evalue = flank_evalue, 
                          flank_plencutoff = flank_plencutoff, 
                          LTRSONLY = True)
            else:
                AutoAlign(I = MCL_I, 
                          part = mafft_align_region, 
                          rmgeneconv = False, 
                          minClustSize = MinClustSize, 
                          align = 'clusters', 
                          rmhomologflank = False, 
                          clustering_method = 'MCL', 
                          WickerParams = None, 
                          auto_outgroup = False, 
                          bpflank = bpflank, 
                          combine_and_do_small_clusters = SMALLS, 
                          flank_pId = flank_pId, 
                          flank_evalue = flank_evalue, 
                          flank_plencutoff = flank_plencutoff, 
                          LTRSONLY = False)
                if MODELTEST:
                    modeltest(iters = 1, 
                              I = MCL_I, 
                              removegeneconv = remove_GC_from_modeltest_aln, 
                              part = mafft_align_region, 
                              clustering_method = 'MCL', 
                              WickerParams = None, 
                              minClustSize = MinClustSize, 
                              bpflank = bpflank, 
                              combine_and_do_small_clusters = SMALLS)
            if GENECONV_G0:
                geneconvClusters(g = '/g0', 
                                 clust = None, 
                                 I = MCL_I, 
                                 minClustSize = MinClustSize, 
                                 clustering_method = 'MCL', 
                                 WickerParams = None, 
                                 combine_and_do_small_clusters = SMALLS)
            if GENECONV_G1:
                geneconvClusters(g = '/g1', 
                                 clust = None, 
                                 I = MCL_I, 
                                 minClustSize = MinClustSize, 
                                 clustering_method = 'MCL', 
                                 WickerParams = None, 
                                 combine_and_do_small_clusters = SMALLS)
            if GENECONV_G2:
                geneconvClusters(g = '/g2', 
                                 clust = None, 
                                 I = MCL_I, 
                                 minClustSize = MinClustSize, 
                                 clustering_method = 'MCL', 
                                 WickerParams = None, 
                                 combine_and_do_small_clusters = SMALLS)
        if SOLOLTR:
            SoloLTRsearch(I = MCL_I, 
                          clustering_method = 'MCL', 
                          WickerParams = {'pId':80,'percAln':80,'minLen':80})
        if CIRCOS:
            Circos(window = '1000000', 
                   plots = 'clusters', 
                   I = MCL_I, 
                   clustering_method = 'MCL', 
                   WickerParams = None)
        # Runs if need to use geneconvLTRs or estimate divergences
        align_ltrs(I = MCL_I, 
                   clustering_method = 'MCL', 
                   WickerParams = None)
        if GENECONV_G0:
            geneconvLTRs(g = '/g0', 
                         I = MCL_I, 
                         clustering_method = 'MCL', 
                         WickerParams = None)
        if GENECONV_G1:
            geneconvLTRs(g = '/g1', 
                         I = MCL_I, 
                         clustering_method = 'MCL', 
                         WickerParams = None)
        if GENECONV_G2:
            geneconvLTRs(g = '/g2', 
                         I = MCL_I, 
                         clustering_method = 'MCL', 
                         WickerParams = None)
        ltr_divergence(I = MCL_I, 
                         clustering_method = 'MCL', 
                         WickerParams = None, 
                         model = model)
        phylo(removegeneconv = False, 
              BOOTSTRAP = True, 
              I = MCL_I, 
              align = 'cluster', 
              removehomologouspair = RMHOMOFLANK, 
              part = mafft_align_region, 
              clustering_method = 'MCL', 
              WickerParams = None, 
              auto_outgroup = AUTO_OUTGROUP, 
              bootstrap_reps = bootstrap_reps, 
              minClustSize = MinClustSize, 
              convert_to_ultrametric = ULTRAMETRIC, 
              bpflank = bpflank, 
              combine_and_do_small_clusters = SMALLS)
    print('fin!')
    sys.exit()
