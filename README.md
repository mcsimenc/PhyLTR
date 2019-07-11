![](https://github.com/mcsimenc/PhyLTR/blob/master/_web/GraphicTree.png)
PhyLTR was designed and created mainly in the [Der Lab](http://www.fullerton.edu/biology/People/faculty/derjp.php) at California State University, Fullerton. The main program is written in Python 3 and tested on Scientific and Ubuntu Linux. Many of the routines in PhyLTR are parallelized for a single computer using Python's multiprocessing module. As the pipeline runs, paths to intermediate results are stored in the file `PhyLTR.output/status`. If the execution is interrupted and restarted, this file is used to allow PhyLTR to resume approximately where it left off.

## Input
A nucleotide FASTA file (e.g. genome assembly)
## Output
* LTR-R annotations in GFF3 format
* Putative solo LTR annotations GFF3 format
* Clusterings
* Rooted, bootstrapped phylogenies as phylograms and ultrametric trees for each cluster in newick format
* LTR divergence estimates
* Gene conversion assessments
* Visualizations and modeling of branching dynamics (external scripts)

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
Within each clustering directory, e.g `MCL/I6/`, are:
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
---
## Default settings
 If none of the 11 flags below (excluding --fasta and --procs) are specified, then the default settings are used: All 11 flags are included. Other flags may be used to modify settings in either default or specific mode. The processes specified by the flags are explained below with additional optional flags. If any of these 11 flags are included in the call then it and all prerequisite processes are carried out, but not any of the others unless they are also included.
```
phyltr --fasta <input> --procs <int> \
	--ltrharvest \
	--ltrdigest \
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
---
## Global options
```
--keep_files				Keeps more intermediate files.
--output_dir (PhyLTR.output)		Output directory.
--logfile (<output_dir>/log.txt)	Desired file for log 
--min_clust_size (7)			Do not align/infer phylogenies from clusters smaller than this.
--nosmalls				Do not combine and assemble all clusters smaller than --min_clust_size
```
---
---
## MAFFT options
These apply to all steps that do alignments. Aligning can fail to complete if the input data are too big and/or not enough RAM is available, for example, MAFFT exhausted 256 Gb RAM with ~2.7k seqs of length >5kb. You can cap the size of clusters to align using `--mafft_largeAln_maxclustsize`. The MAFFT algorthim FFT-NS-2 is used for small and medium clusters and FFT_NS-1, which is much more inaccurate, for large clusters. The alignment of clusters has been the limiting process in terms of computation time in my experience, but aligning can be sped up by reducing the number of improvement iterations performed. You can use the following options to designate ranges of cluster sizes and set the number of alignment improvement iterations for each size (fewer iterations = faster).
###### Options
```
--mafft_align_region (entire)		Sequence region to align. entire or internal. Internal = between LTRs
--maxiterate_small_clusters (30)	MAFFT iterations for small clusters. (more = better alignment = slower)
--maxiterate_medium_clusters (3)	MAFFT iterations for medium clusters. (more = better alignment = slower)
--mafft_smallAln_maxclustsize (50)	Max elements to consider a cluster small.
--mafft_mediumAln_maxclustsize (500)	Max elements to consider a cluster medium.
--mafft_largeAln_maxclustsize (1000)	Max elements to consider a cluster large. Clusters larger than this will not be aligned.
```
---
---
# PROCESSES
---
---
## 1. Identify candidate long terminal repeat retrotransposon (LTR-R) loci
#### Run LTRharvest: `--ltrharvest`
![](https://github.com/mcsimenc/PhyLTR/blob/master/_web/step1output.png)
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
---
## 2. Identify putatve protein-coding domains in LTR-R internal regions.
If both are run, LTRdigest runs first, then the ORF-finding routine.
#### A. Run LTRdigest: `--ltridgest`
![](https://github.com/mcsimenc/PhyLTR/blob/master/_web/step2output.png)
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
#### B. Run open reading frame (ORF) annotation: `--findORFs`
![](https://github.com/mcsimenc/PhyLTR/blob/master/_web/step3output.png)
###### Description
Internal regions are searched for ORFs that don't overlap any preexisting annotation from step 2A and are longer than a user-defined threshold.
###### Output
* Modified version of step 2A or 1 output with ORF annotations + translated sequences (GFF3)
###### Options
```
--min_orf_len (300)	The minimum length (bp) of ORF to annotate
```
###### External dependencies
* BEDtools
* EMBOSS
---
---
## 3. Classify elements using homology to LTR-Rs in Dfam and/or Repbase and keeping only elements classified as LTR-Rs is on by default. Note, some Dfam and Repbase LTR-R entries are labelled 'Unknown'.
#### A. Turn off both Repbase and Dfam classification: `--no_classification`
![](https://github.com/mcsimenc/PhyLTR/blob/master/_web/step4output.png)
###### Description
Both methods use homology-based evidence for classifying elements as one of the classifications obtained from the database records: BEL, Copia, DIRS, Endogenous Retrovirus, ERV1, ERV2, ERV3, ERV4, Gypsy, Lentivirus (Repbase) and Copia, DIRS, ERV?, ERV1, ERV1?, ERV-Foamy, ERVK, ERVK?, ERVL, ERVL?, ERVL-MaLR, Gypsy, Gypsy?, Ngaro, Pao, Undefined, Unknown (Dfam Superfamily). As currently implemented, if both classifications are not exactly the same, the Dfam classification is used for the annotated LTR-R label. This worked well for our test genomes and classifications matched domain annotations (e.g. gypsy elements always and only had gypsy domains).
###### Output
* LTR-R annotations with false positives removed and Dfam and/or Repbase annotations (GFF3)
* LTR-R annotations separated by classification (GFF3s)
#### B. Do not run Dfam-based classification: `--no_dfam`
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
#### C. Do not run Repbas-based classification: `--no_repbase`
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
---
## 4. Cluster
#### A. Run WickerFam clustering: `--wicker`
###### Description
An implementation of the method suggested by Wicker et al. (2007) for circumscribing putative LTR-R families. Elements are assigned to the same cluster if at least **A** percentage of an alignment of either their LTRs or internal regions at least **B** bp in length has at least **C** % identity, with defaults being 80%, 80 bp, and 80%. Separate all-by-all blastns are performed for LTRs and/or internal regions of elements in each classification and used to construct a graph, subject to constraints (see options below). Clusters are assigned the elements in the connected components in the graph which are discovered using depth first search.
###### Output
* Cluster assignments for each element (TSV)
* Subsets of main annotation file for every cluster (GFF3s)
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
An implementation of the protocol _Clustering similarity graphs encoded in BLAST results_ in the protocols section on the [MCL website](https://micans.org/mcl/). Separate all-by-all blastns of entire LTR-R sequences are performed for elements in each classification and MCL programs are used to carry out MCL clustering. A single parameter set by --I controls the output. See MCL documentation for a discussion of setting this parameter.
###### Output
* Cluster assignments for each element (TSV)
* Subsets of main annotation file for every cluster (GFF3s)
###### Options
```
--I (6)		MCL inflation paramter. Larger values result in smaller clusters.
```
###### External dependencies
* BEDtools
* NCBI Blast+
* MCL
---
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
Estimates sequence divergences (substitutions per site) between LTRs for each element. Models for nucleotide subsitution are selected using jModeltest2. Input for jModeltest2 is a multiple sequence alignment of LTRs (generated by MAFFT) and a tree. Trees are inferred using FastTree2 with GTR from multiple sequence alignments of entire elements from each cluster performed by MAFFT and trimmed by trimAl. PAUP\* is used to estimate divergences using the best models and corresponding ML-estimated parameters. Divergences are scaled linearly by the proportion of aligment represented by a g0 GENECONV tract.
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
---
## 6. "Solo LTR" search
#### Run "Solo LTR" search: `--soloLTRsearch`
###### Description
Finds candidate solo LTRs in the input fasta by performs blastn of each LTR to the input FASTA and retains the highest scoring hits that do not overlap a full-length LTR-R or other blastn hits. These hits may be referred to as "solo LTRs", but they may represent the following categories of features: (1) LTRs formed by unequal recombination (2) LTRs from truncated elements formed by illegitimate recombination (3) LTRs of elements that were discarded as false positives. (4) LTRs of elements that were missed by LTRharvest. The fraction of candidate solo LTRs in category 3 and 4 is expected to be small. The constraints of this search are consistent with those in the WickerFam classification.
###### Output
* Summary file with cluster, cluster size (full-length elements), and number of solo LTRs (TSV)
* Separate solo LTR annotation files for each classification and for each cluster (GFF3s)
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
---
## 7. Gene conversion assessment between elements in clusters
#### Run GENECONV: `--geneconvclusters`
###### Description
Finds evidence of gene conversion in multiple alignments of entire elements for each cluster using GENECONV.
###### Output
Aggregated GENECONV results (TSV)
Summary file (TSV)
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
Runs Circos to make plots showing gene conversion tracts between elements as links. (g0=orange, g2=blue, g1=purple)
###### Output
* Circos plots showing the position of LTR-Rs on their parent sequence with GC tracts as links
* Circos plots showing extent of GC tracts along lengths of elements (i.e. the parent sequences are the elements in a cluster.
###### External dependencies
* BEDtools
* Circos
---
---
## 8. Trees
#### Infer phylogenies for each cluster: `--phylo`
###### Description
Infers phylogenies from alignments of entire elements for each cluster, optionally with an outgroup, using FastTree2 w/GTR. The alignment is resampled for bootstrapping using SEQBOOT from PHYLIP and optionally converted into an ultrametric tree using PATHd8.
###### Output
* Trees (Newick/parenthetical)
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
---
## 9. External scripts
---
---
#### A. Search genes for LTR-R ORF homologs
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
---
---
#### B. Render trees annotated with LTR-R diagrams with colored ORFs
![](https://github.com/mcsimenc/PhyLTR/blob/master/_web/Ete3example.png)
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
	A file containing a list of element IDs (e.g. LTR_retrotransposon123) to mark with a green asterisk.
	Instead of an asterisk, a 'T' is shown as if -classif is used.

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
---
#### C. R code
* Distributions
* Lineage through time plots
* Transposition/deletion rate modeling
* Tree shape analysis
---
---
# Appendices
---
---
## Appendix A. All options

|General|Default|
|:---:|:---:|
|-h or --h||
|-help or --help||
|--logfile|log.txt|
|--procs|1|
|--output_dir|PhyLTR.output|
|--keep_files||
|--min_clust_size|7|
|--nosmalls||

|Input||
|:---:|:---:|
|--fasta||

|LTRharvest|Default|
|:---:|:---:|
|--ltrharvest||
|--del|-3|
|--ins|-3|
|--mis|-2|
|--mat|2|
|--xdrop|5|
|--minlenltr|100|
|--maxlenltr|1000|
|--mindistltr|1000|
|--maxdistltr|15000|
|--similar|0.0|
|--vic|60|
|--mintsd|4|
|--maxtsd|20|

|LTRdigest|Default|
|:---:|:---:|
|--ltrdigest||
|--ltrdigest_hmms|PhyLTR/RepeatDatabases/LTRdigest_HMMs/hmms|

|ORFs|Default|
|:---:|:---:|
|--findORFs|
|--min_orf_len|300|

|Classification|Default|
|:---:|:---:|
|--classify||
|--classify_dfam||
|--nhmmer_reporting_evalue|10|
|--nhmmer_inclusion_evalue|1e-5|
|--classify_repbase||
|--repbase_tblastx_evalue|1e-5|
|--keep_conflicting_classifications||
|--keep_no_classification||

|Cluster|Default|
|:---:|:---:|
|--wicker||
|--wicker_pId|80|
|--wicker_pAln|80|
|--wicker_minLen|80|
|--wicker_no_ltrs||
|--wicker_no_internals||
|--mcl||
|--I|6|

|Find solo LTRs|Default|
|:---:|:---:|
|--soloLTRsearch||
|--soloLTRminPid|80.0|
|--soloLTRminLen|80.0|
|--soloLTRmaxEvalue|1e-3|

|Alignment settings|Default|
|:---:|:---:|
|--maxiterate_small_clusters|20|
|--maxiterate_medium_clusters|3|
|--mafft_smallAln_maxclustsize|50|
|--mafft_mediumAln_maxclustsize|500|
|--mafft_largeAln_maxclustsize|1000|

|LTR divergence and gene conversion|Default|
|:---:|:---:|
|All alignment settings||
|--geneconvltrs||
|--geneconvclusters||
|--geneconv_g|g0,g1,g2|
|--remove_GC_from_modeltest_aln||
|--ltrdivergence||
|--circos||

|Phylogenetics|Default|
|:---:|:---:|
|All alignment settings||
|--phylo||
|--LTT||
|--rmhomoflank||
|--bpflank|500|
|--flank_evalue|1e-5|
|--flank_pId|70.0|
|--flank_plencutoff|70.0|
|--auto_outgroup||
|--bootstrap_reps|100|
|--convert_to_ultrametric||

## Appendix B. References
Bao, W., Kojima, K. K., & Kohany, O. (2015). Repbase Update, a database of repetitive elements in eukaryotic genomes. Mobile DNA, 6(1), 11. http://doi.org/10.1186/s13100-015-0041-9

Britton, T., Anderson, C. L., Jacquet, D., Lundqvist, S., & Bremer, K. (2007). Estimating divergence times in large phylogenetic trees. Systematic Biology, 56(5), 741-752. http://doi.org/10.1080/10635150701613783

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10(1), 421. http://doi.org/10.1186/1471-2105-10-421

Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), 1972-1973. http://doi.org/10.1093/bioinformatics/btp348

Casola, C., Ganote, C. L. & Hahn, M. W. Nonallelic Gene Conversion in the Genus Drosophila. Genetics 185, 95-103 (2010).

Cossu, R. M., Casola, C., Giacomello, S., Vidalis, A., Scofield, D. G., & Zuccolo, A. (2017). LTR Retrotransposons Show Low Levels of Unequal Recombination and High Rates of Intraelement Gene Conversion in Large Plant Genomes. Genome Biology and Evolution, 9(12), 3449-3462. http://doi.org/10.1093/gbe/evx260

Darriba, D., Taboada, G. L., Doallo, R., & Posada, D. (2012). jModelTest 2: more models, new heuristics and parallel computing. Nature Methods, 9(8), 772-772. http://doi.org/10.1038/nmeth.2109

van Dongen, Stijn. Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht, May 2000.

Eddy, S. R. A new generation of homology search tools based on probabilistic inference. Genome Inform 23, 205–211 (2009).

Ellinghaus, D., Kurtz, S., & Willhoeft, U. (2008). LTRharvest, an efficient and flexible software for de novo detection of LTR retrotransposons. BMC Bioinformatics, 9(1), 18-14. http://doi.org/10.1186/1471-2105-9-18

Felsenstein, J. 2005. PHYLIP (Phylogeny Inference Package) version 3.6. Distributed by the author. Department of Genome Sciences, University of Washington, Seattle.

Finn, R. D., Coggill, P., Eberhardt, R. Y., Eddy, S. R., Mistry, J., Mitchell, A. L., et al. (2016). The Pfam protein families database: towards a more sustainable future. Nucleic Acids Research, 44(D1), D279-85. http://doi.org/10.1093/nar/gkv1344

Gremme, G., Steinbiss, S., & Kurtz, S. (2013). GenomeTools: a comprehensive software library for efficient processing of structured genome annotations. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 10(3), 645-656. http://doi.org/10.1109/TCBB.2013.68

Hubley, R., Finn, R. D., Clements, J., Eddy, S. R., Jones, T. A., Bao, W., et al. (2016). The Dfam database of repetitive DNA families. Nucleic Acids Research, 44(D1), D81-9. http://doi.org/10.1093/nar/gkv1272

Huerta-Cepas, J., Serra, F., & Bork, P. (2016). ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Molecular Biology and Evolution, 33(6), 1635-1638. http://doi.org/10.1093/molbev/msw046

Katoh, K., & Standley, D. M. (2013b). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution, 30(4), 772-780. http://doi.org/10.1093/molbev/mst010

Krzywinski, M., Schein, J., Birol, I., Connors, J., Gascoyne, R., Horsman, D., et al. (2009). Circos: an information aesthetic for comparative genomics. Genome Research, 19(9), 1639-1645. http://doi.org/10.1101/gr.092759.109

Llorens, C., Futami, R., Covelli, L., Domínguez-Escribá, L., Viu, J. M., Tamarit, D., et al. (2011). The Gypsy Database (GyDB) of mobile genetic elements: release 2.0. Nucleic Acids Research, 39(Database issue), D70-4. http://doi.org/10.1093/nar/gkq1061

Price, M. N., Dehal, P. S. & Arkin, A. P. FastTree 2--approximately maximum-likelihood trees for large alignments. PLoS ONE 5, e9490 (2010).

Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842. http://doi.org/10.1093/bioinformatics/btq033

Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: the European Molecular Biology Open Software Suite. Trends in Genetics, 16(6), 276-277.

Sawyer, S.A. (1999) GENECONV: A computer package for the statistical detection of gene conversion. Distributed by the author, Department of Mathematics, Washington University in St. Louis, available at http://www.math.wustl.edu/~sawyer.

Steinbiss, S., Willhoeft, U., Gremme, G., & Kurtz, S. (2009). Fine-grained annotation and classification of de novo predicted LTR retrotransposons. Nucleic Acids Research, 37(21), 70027013. http://doi.org/10.1093/nar/gkp759

Swofford, D. L. 2003. PAUP\*. Phylogenetic Analysis Using Parsimony (\* = and Other Methods). Version 4. Sinauer Associates, Sunderland, Massachusetts.

Wicker, T. et al. A unified classification system for eukaryotic transposable elements. Nat Rev Genet 8, 973–982 (2007).
