![](https://github.com/mcsimenc/PhyLTR/blob/master/GraphicTree.png)
The main program is written in Python 3 and tested on Scientific and Ubuntu Linux. Many of the routines in PhyLTR are parallelized, but not for MPI. As the pipeline runs, paths to intermediate results like alignments are stored in the file `PhyLTR.output/status`. If the execution is interrupted, this file is used to allow PhyLTR to resume more or less where it left off.

## Input
A nucleotide FASTA file
## Output
PhyLTR populates a directory structure, keeping results for different clusterings separate. The results obtained prior to clutering are used for any post-clustering analysis. Within `PhyLTR.output` (default):
```
LTRharvest/		LTRharvest results
suffixerator/		suffix array used by LTRharvest
LTRdigest/		LTRdigest results
AnnotateORFs/		intermediate files for ORF annotation
Circos/			Circos intermediate files and plots
DfamClassification/	nhmmer search of Dfam for LTR-R homologs
RepbaseClassification/	tblastx search of Repbase for LTR-R homologs
FASTA_output/		some intermediate FASTA files
GFF_output/		various intermediate and final GFF3 files
WickerFamDir/		WickerFam clusterings and downstream analyses
MCL/			MCL clusterings and downstream analyses
```
Within each clustering directory, e.g: `MCL/I6/`
```
Alignments/		All alignments
Clusters/		Results of clustering
FASTAs/			Intermediate FASTA files
GENECONV/		Gene conversion analyses
GFFs/			Intermediate GFF3 files
LTR_divergence/		LTR divergence estimation analyis
Modeltest/		Model testing files
SoloLTRsearch/		"Solo LTR" search files
Trees/			Phylogenetic analyses
```
---
## Default settings
If phyltr is run without any flags specifying a task, all tasks are run (below). The following two calls are equivalent. The processes specified by the flags in the second call are explained below with additional optional flags. Some of the processes modify the GFF3 file that is used for downstream analyses.
```
phyltr --fasta <input> --procs <int>

phyltr --fasta <input> --procs <int> \
	--ltrharvest \
	--ltrdigest \
	--classify \
	--wicker \
	--mcl \
	--geneconvclusters \
	--circos \
	--sololtrsearch \
	--geneconvltrs \
	--ltrdivergence \
	--phylo \
	--LTT
```
---
## 1. Identify candidate long terminal repeat retrotransposon (LTR-R) loci
#### Run LTRharvest: `--ltrharvest`
###### Description
LTRharvest finds loci with the expected structure of full-length LTR-Rs, TSD-LTR-LTR-TSD (TSD = targest site duplication). LTRharvest searches a suffix array to make things fast, which PhyLTR creates from the FASTA input with the GenomeTools program suffixerator.
###### Output
* Annotations of possible LTR-Rs, including the subfeatures TSDs, LTRs, primer-binding sites (GFF3)
###### Options
```
--minlenltr (100)	Minimum LTR length (bp)
--maxlenltr (1000)	Maximum LTR length (bp)
--mindistltr (1000)	Minimum distance between LTRs (bp)
--maxdistltr (15000)	Maximum distance between LTRs (bp)
--similar (0.0)		Minimum % similarity between LTRs
--vic (60)		Distance (bp) beyond LTRs to search for TSDs
--mintsd (4)		Minimum length for each TSD
--maxtsd (20)		Maximum length for each TSD
--xdrop	(5)		xdropbelow score for extension-alignment
--mat (2)		Match score for extension-alignment
--mis (-2)		Mismatch score for extension-alignment
--insi (-3)		Insertion score for extension-alignment
--del (-3)		Deletion score for extension-alignment
```
###### External dependencies
* GenomeTools
---
## 2. Identify putatve protein-coding domains in LTR-R internal regions.
If both are run, LTRdigest runs first, then the ORF-finding routine.
#### A. Run LTRdigest: `--ltridgest`
###### Description
LTRdigest coordinates HMMER3 searches for transposable element protein coding sequence homologs in the internal regions (between the LTRs) of the putative LTR-Rs from step 1. A set of TE-related pHMMs is included with PhyLTR.
###### Output
* Modified version of step 1 GFF3 that also has protein-coding domain annotations (GFF3)
* Nucleotide sequences for each domain (FASTA)
* TSD motifs and domain orders for each element (TSV)
###### Options
```
--ltrdigest_hmms (PhyLTR/RepeatDatabases/LTRdigest_HMMs/hmms)	path to pHMMs
```
###### External dependencies
* GenomeTools
* HMMER3
* pHMMs
#### B. Run Open reading frame (ORF) annotation: `--findORFs`
###### Description
Internal regions are searched for ORFs that don't overlap any preexisting annotation from step 2A and are longer than a user-defined threshold.
###### Output
* Modified version of step 2A or 1 output with ORF annotations + translated sequences (GFF3)
###### Options
```
--min_orf_len (300)	The minimum length (bp) of ORF to find
```
###### External dependencies
* BEDtools
* EMBOSS
---
## 3. Classify elements using homology to LTR-Rs in Dfam and/or Repbase and remove false positives
#### A. Run both Repbase and Dfam classification: `--classify`
###### Description
Both methods use homology-based evidence for classifying elements as one of the classifications obtained from the database records: BEL, Copia, DIRS, Endogenous Retrovirus, ERV1, ERV2, ERV3, ERV4, Gypsy, Lentivirus (Repbase) and Copia, DIRS, ERV?, ERV1, ERV1?, ERV-Foamy, ERVK, ERVK?, ERVL, ERVL?, ERVL-MaLR, Gypsy, Gypsy?, Ngaro, Pao, Undefined, Unknown (Dfam Superfamily). As currently implemented, Dfam hits trump Repbase because they are expected to be longer, and only the highest scoring hits are considered. This worked well for our test genomes and classifications matched domain annotations (e.g. gypsy elements always and only had gypsy domains).
###### Output
* LTR-R annotations with false positives removed and Dfam and/or Repbase annotations (GFF3)
* LTR-R annotations separated by classification (GFF3s)
#### B. Run Dfam classification: `--classify_dfam`
###### Description
Finds homologs in Dfam using nhmmer
###### Output
* LTR-R annotations with best Dfam hit in attributes (GFF3)
###### Options
```
--keep_no_classifications 		Retain elements without homology to known LTR-Rs
--nhmmer_reporting_evalue (10)		See HMMER3 documentation: nhmmer -E
--nhmmer_inclusion_evalue (1e-2)	See HMMER3 documentation: nhmmer -incE
```
###### External dependencies
* BEDtools
* HMMER3
* Dfam
#### C. Run Repbase classification: `--classify_repbase`
###### Description
Finds homologs in Repbase using tblastx
###### Output
* LTR-R annotations with best Repbase hit in attributes (GFF3)
###### Options
```
--keep_no_classifications	Retain elements without homology to known LTR-Rs
--repbase_tblastx_evalue (1e-5)	Maximum E-value for tblastx hits
```
###### External dependencies
* BEDtools
* NCBI BLAST+
* Repbase
---
## 4. Cluster LTR-Rs
#### A. Run WickerFam clustering: `--wicker`
###### Description
An implementation of the method suggested in Wicker et al. (2007) for circumscribing putative LTR-R families. Elements are assigned to the same cluster if in an alignment of either their LTRs or internal regions share a sequence of some minimum % identity that is at least so many base pairs in length and over some percentage of the length of the inputs, with the defaults being 80%, 80 bp, and 80%. Separate all-by-all blastns are performed for LTRs and/or internal regions of elements in each classification and used to construct a graph, subject to constraints (see options below). Clusters are assigned the elements in the connected components in the graph which are discovered using depth first search.
###### Output
* Cluster assignments for each element (TSV)
###### Options
```
--wicker_minLen (80)	Minimum length of blastn alignment
--wicker_pAln (80)	Minimum percent of LTR or internal region required in alignment
--wicker_pId (80)	Minimum %identity in alignment
--wicker_no_internals	Turns off use of internal region alignments for clustering
--wicker_no_ltrs	Turns off use of LTR alignments for clustering
```
###### External dependencies
* BEDtools
* NCBI Blast+
#### B. Run MCL clustering: `--mcl`
###### Description
An implementation of the protocol "Clustering similarity graphs encoded in BLAST results" in the protocols section on the MCL website (https://micans.org/mcl/). Separate all-by-all blastns of entire LTR-R sequences are performed for elements in each classification and MCL programs are used to carry out MCL clustering. A single parameter set by --I controls the output. See MCL documentation for a discussion of setting this parameter.
###### Output
* Cluster assignments for each element (TSV)
###### Options
```
--I (6)		MCL inflation paramter. Larger values result in smaller clusters.
```
###### External dependencies
* BEDtools
* NCBI Blast+
* MCL
---
## 5. Estimate LTR divergences
#### A. Run gene conversion assessment on LTR pairs for each element: `--geneconvltrs`
###### Description
Finds evidence of gene conversion between the LTRs of each element using GENECONV as in Cossu et al. (2017). If this step is run, gene conversion estimates will be used to scale LTR divergence estimates as in Casola et al. (2010) in the next step (5B).
###### Output
* Aggregated GENECONV results (TSV)
###### Options
```
--geneconv_g (g1,g2,g3)	Comma-separated list, g1, g2, and/or g3. Stringency for mismatch-free gene conversion tracts: g0 > g2 > g1
```
###### External dependencies
* BEDtools
* MAFFT
* trimAl
* GENECONV
#### B. Run LTR divergence estimation: `--ltr_divergence`
###### Description
Estimates sequence divergences (substitutions per site) between LTRs for each element. Models for nucleotide subsitution are selected using jModeltest2. Input for jModeltest2 is a multiple sequence alignment of LTRs (generated by MAFFT) and tree. Trees are inferred using FastTree2 with GTR from multiple sequence alignments of entire elements from each cluster performed by MAFFT and trimmed by trimAl. PAUP\* is used to estimate divergences using the best models and corresponding ML-estimated parameters. Divergences are scaled linearly by the proportion of aligment represented by a g0 GENECONV tract.
###### Output
* Raw and scaled divergence estimates and best supported models (TSV)
###### Options
```
--remove_GC_from_modeltest_aln	Remove elements with gene conversion (--geneconvclusters)
```
###### External dependencies
* BEDtools
* MAFFT
* trimAl
* FastTree2
* jModelTest2
* PAUP\*
---
## 6. "Solo LTR" search
#### Run "Solo LTR" search: `--soloLTRsearch`
###### Description
Finds candidate solo LTRs in the input fasta. All-by-all blastn, highest scoring non-overlapping hit.
###### Output
###### Options
```
--soloLTRminPid (80.0)		Minimum %identity in blastn alignment to associate LTR with a cluster
--soloLTRminLen	(80.0)		Minimum % of length of LTR required in alignment to associate LTR with a cluster
--soloLTRmaxEvalue (1e-3)	Maximum E-value for blastn
```
###### External dependencies
* BEDtools
* NCBI BLAST+
---
## 7. Gene conversion assessment between elements in clusters
#### Run GENECONV: `--geneconvclusters`
###### Description
###### Output
###### Options
```
--geneconv_g (g1,g2,g3)	Comma-separated list, g1, g2, and/or g3. Stringency for mismatch-free gene conversion tracts: g0 > g2 > g1
```
###### External dependencies
* BEDtools
* MAFFT
* trimAl
* GENECONV
#### Make Circos plots: `--circos`
###### Description
###### Output
###### External dependencies
* BEDtools
* Circos
---
## 8. Phylogenetics
#### Run phylogenetic inference: `--phylo`
###### Description
###### Output
###### Options
```
--min_clust_size (7)		Do not align clusters smaller than this.
--nosmalls			Do not analyze clusters smaller than --min_clust_size
--rmhomoflank			Exclude elements with non-unique flanking sequences from alignments.
--bpflank (500)			Length (bp) of flank searched
--flank_evalue (1e-5)		Maximum E-value for blastn for flank search
--flank_pId (70.0) 		Minimum %identity in blastn alignment for flank search
--flank_plencutoff (70.0)	Minimum % of length of flank required in alignment
--convert_to_ultrametric	Convert tree to ultrametric (can be used for LTT plots)
--auto_outgroup			Include as an outgroup a random element from the largest available other cluster in classification (e.g. gypsy)
--bootstrap_reps (100)		Number of bootstrap replicates to perform
--LTT				Turns on --rmhomoflank, --convert_to_ultrametric, and --auto_outgroup.
```
###### External dependencies
* BEDtools
* MAFFT
* trimAl
* FastTree2
* PATHd8
* PHYLIP
---
## 9. External scripts
---
#### A. Search genes for LTR-R ORF homologs
###### Description
###### Output
###### Documentation
```
 Usage:
------------
domainSearch.py -gff <gff3> -ref <fasta> -prot <fasta> -procs <int>

Description:
------------
This script was written for comparing putative protein-coding domains
and ORFs in a GFF3 as output from PhyLTR (LTRharvest+LTRdigest+AnnotateORFs).
The required programs are bedtools getfasta, and makeblastdb and blastp from
BLAST+. Their paths should be in a text file named CONFIG file located in the
same directory as this script using the same format used for PhyLTR.

Process:
-----------
1. Putative domain sequences encoded in -gff are extracted from -ref and
   translated
2. Translated sequences are compared to protein sequences in -prot, which could
   be, for example, the protein sequences for a

Mandatory flags:
-----------
-gff	 <path>	GFF3 file with LTRdigest-format LTR retrotransposon features

-ref	 <path>	Nucleotide FASTA file that is the reference associated with -gff

-prot	 <path>	Protein FASTA file to convert to a blast database and

-procs	 <int>	Number of processors to use for blastp

-out	 <path> Output file path - you may need to give full path

Optional flags:
-----------
-evalue <num>	max E-value for blastp. Default 1e-2
```
###### External dependencies
* Python 3
* BEDtools
* NCBI BLAST+
#### B. Render graphical trees annotated with LTR-R diagrams with colored ORFs
###### Documentation
```
Usage:
	ete3_tree.py -t <newick_tree> -d <divergences> -g <gff> [options] 

Description:

	Draws trees for long terminal repeat retrotransposon phylogenies with diagrams of LTR elements' domain
	architecture.

-lflabel
	Show element IDs (integers) as leaf labels.

-classif
	Add the superfamily classification for each element (which is obtained from the LTR divergence file) 
	above each LTR RT diagram.

-geneconv
	Add the word 'Yes' or 'No' to the immediately to the right of each leaf depending on whether intra-element
	gene conversion tracts were detected between the LTRs of that element.

-reroot <int>|auto
	Two options are possible for -reroot: 'auto', or an <int> corresponding to the element name (i.e.
	LTR_retrotransposon<int>) to position as the earliest diverging lineage. Only use -reroot auto if
	the newick filename contains the outgroup in the format output by PhyLTR. 

-ultrametric
	Draw tree after applying ete3's convert_to_ultrametric() function. In my experience using this option, I
	have always seen terminal taxa drawn at different horizontal positions making the tree not look ultrametric.

-1, 2, ..., -n  <path> <str>
	Coloring ORFs. <path> is list of ORF IDs, str is ETE3 color


-transcribed <path>
	A file containing a list of element IDs (e.g. LTR_retrotransposon123) to mark with a green asterix.
	Instead of an asterix, a 'T' is shown as if -classif is used.

-round
	Draw ORFs with rounded corners.

-outfmt <str>
	One of (pdf, png, svg). Default=pdf


Colors
	protease	yellow
	gag		deepskyblue
	copia gag	mediumblue
	DUF4219		lime
	rt		red
	rnaseh		darkviolet
	int		magenta
	zf-h2c2		pink
	other		silver w/black text
	orf		dimgray
	annotated orf	custom
```
###### External dependencies
* Python 3
* Python 3 modules: ETE3
---
#### C. Visualize insertion ages
###### Description
###### Output
###### External dependencies
* R
* R packages:  hash, ggplot2
---
#### D. Lineage through time plots
###### Description
###### Output
###### External dependencies
* R
* R packages: ape, hash, ggplot2
---
#### E. Transposition rate analyses
###### Description
###### Output
###### External dependencies
* R
* R packages: ape, hash, ggplot2
---
#### F. Tree shape analyses
###### Description
###### Output
###### External dependencies
* R
* R packages: ape, apTreeshape, hash, ggplot2
---
#### G. Transposition/deletion rate modeling
###### Description
###### Output
###### External dependencies
* R
* R packages: ape, phangorn, hash, ggplot2, LASER (functions included as LASER is deprecated)
---
#### H. Other scripts
###### Description
###### Output
---
## APPENDIX A. Global MAFFT options
This step has been the limiting process in my experience. It can be sped up by reducing the number of iterations performed for each alignment and by reducing the maximum number of elements for classiying elements as medium and small clusters. MAFFT exhausted 256 Gb RAM with ~2.7k seqs of length >5kb. Depending on resources available to you, you may need to cap the size of clusters to align using `--mafft_largeAln_maxclustsize`. Default is to not align clusters with >1000 elements. The MAFFT algorthim FFT-NS-2 is used for small and medium clusters and FFT_NS-1 for large clusters, which is much more inaccurate.
###### Options
```
--maxiterate_small_clusters (30)	MAFFT iterations. More will improve alignment quality.
--maxiterate_medium_clusters (3)	MAFFT iterations. More will improve alignment quality.
--mafft_smallAln_maxclustsize (50)	Max elements to consider a cluster small.
--mafft_mediumAln_maxclustsize (500)	Max elements to consider a cluster medium.
--mafft_largeAln_maxclustsize (1000)	Max elements to consider a cluster large. Clusters larger than this will not be aligned.
```
---
## APPENDIX B. Other global options
```
--keep_files				Keeps intermediate files that are otherwise removed.
--output_dir (PhyLTR.output)		Output directory. Default is "PhyLTR.output
--logfile (<output_dir>/log.txt)	Path to where log file is written 
--min_clust_size (7)			Do not align clusters smaller than this.
--nosmalls				Do not combine and assemble all clusters smaller than --min_clust_size
```
---
## APPENDIX C. All options
---
## APPENDIX D. Example output
---
## APPENDIX E. References
* GenomeTools
* etc.
---
