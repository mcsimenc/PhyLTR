# PhyLTR: de novo LTR retrotransposon annotation, classification, and phylogenetic analysis

PhyLTR is a software pipeline built from open source software. The main program is written in Python 3. Many of the routines in PhyLTR are parallelized, but it is not written for MPI so parallel components can only run on CPUs that share memory, i.e. on a single node.

As the pipeline runs, paths to intermediate results like alignments are stored in the file PhyLTR/status. If the execution is interrupted, this file is used to allow PhyLTR to resume more or less where it left off.

## Default settings
If phyltr is run without any flags specifying a task, all tasks are run. The following two calls are equivalent. The processes specified by the flags in the second call are explained below with additional optional flags. Some of the processes modify the GFF3 file that is used for downstream analyses.
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
	--DTT
```
## 1. Identifying candidate LTR-R loci with LTRharvest
#### Turn on `--ltrharvest`
###### External dependencies
* GenomeTools
###### Available options (explained in the LTRharvest documentation).
```
--minlenltr (100)
--maxlenltr (1000)
--mindistltr (1000)
--maxdistltr (15000)
--similar (85.0)
--vic (60)
--mintsd (4)
--maxtsd (20)
--xdrop	(5)
--mat (2)
--mis (-2)
--insi (-3)
--del (-3)
```
## 2. Identifying putatve protein-coding domains in LTR-R internal regions.
#### A. Turn on `--ltridgest`
###### External dependencies
* GenomeTools
* HMMER3
* pHMMs (database)
###### Required PhyLTR output
* --ltrharvest
###### Available options
--ltrdigest_hmms (/home/joshd/scripts/PhyLTR/LTRdigest_HMMs/hmms)	path to pHMMs
#### B. Turn on --findORFs
###### External dependencies
* GenomeTools
* EMBOSS
###### Required PhyLTR output
* --ltrharvest
###### Available options
```
--min_orf_len (300)	The minimum length (bp) of ORF to find
```
## 3. Classify elements using homology to LTR-Rs in Dfam and/or Repbase
#### A. Turn on both Repbase and Dfam classification `--classify`
###### Possible external dependencies
* GenomeTools
* BEDtools
* NCBI BLAST+
* HMMER3
* Repbase (database)
* Dfam (database)
###### Required PhyLTR output
* --ltrharvest
###### Available global options
```
--keep_no_classifications Retain elements without homology to known LTR-Rs
```
#### B. Turn on Dfam classification `--classify_dfam`
###### Available options (explained in the HMMER3 documentation)
```
--nhmmer_reporting_evalue (10)
--nhmmer_inclusion_evalue (1e-2)
```
#### C. Turn on Repbase classification `--classify_repbase`
###### Available options (explained in the tblastx documentation)
```
--repbase_tblastx_evalue (1e-5)
```
## 4. Cluster LTR-Rs
###### Possible external dependencies
* NCBI Blast+
* BEDtools
* MCL
###### Required PhyLTR output
* --ltrharvest
#### Turn on WickerFam clustering `--wicker`
###### Available options
--wicker_minLen (80)	Minimum length of blastn alignment
--wicker_pAln (80)	Minimum percent of LTR or internal region required in alignment
--wicker_pId (80)	Minimum %identity in alignment
--wicker_no_internals	Turns off use of internal region alignments for clustering
--wicker_no_ltrs	Turns off use of LTR alignments for clustering
#### B. Turn on MCL clustering `--mcl`
###### Available options
--I (6)
## 5. LTR divergence estimation
###### Possible external dependencies
* MAFFT
* trimAl
* BEDtools
* PAUP\*
* jModelTest2
* GENECONV
#### A. Turn on GENECONV for intra-element LTR assessment `--geneconvltrs`
###### See options for MAFFT below
###### Available options (explained in GENECONV documentation)
--geneconv_g (g1,g2,g3)	Comma-separated list, g1, g2, and/or g3
#### B. Estimate LTR divergences `--ltr_divergence`
###### Available options
## 6. "Solo LTR" search
###### External dependencies
* NCBI BLAST+
#### Turn on "solo LTR" search `--soloLTRsearch`
## 7. Gene conversion between LTR-Rs in a cluster
###### Possible external dependencies
* GENECONV
* BLAST?
* Circos
#### Turn on GENECONV `--geneconvclusters`
#### Turn on Circos `--circos`

## 8. Phylogenetics
###### Possible external dependencies
* MAFFT
* trimAl
* FastTree2 + scripts
* PATHd8
* PHYLIP



LTR divergence estimation
-------------------------
--ltrdivergence			Find statistially best supported (BIC) substitution model for each cluster (default ON)
					and estimate substitutions per site between LTRs for each element. 
--modeltest_criterion		<str>	AIC, AICc, or BIC. (default BIC)
      				MULTIPLE CHOICES NOT YET IMPLEMENTED
--gc_ltr_div_scaling		<int>	For reporting scaled divergence estimates to account for effects of gene conversion, if observed. (default 1)
      				MULTIPLE CHOICES NOT YET IMPLEMENTED
					1. divC = div * len(aln)/(len(aln)-len(gc_trac))
--default_model		<str>	When model testing is not possible (i.e. cluster is too small) HKY85 or JC,
					if HKY85 is not possible due to dinucleotide LTRs.
					MULTIPLE CHOICES NOT YET IMPLEMENTED
--remove_GC_from_modeltest_aln	Remove elements with suspected intra-cluster inter-element gene conversion tracts.

Solo LTR search
-------------------
soloLTRsearch				Turn on solo LTR search.
soloLTRminPid			<num>	Minimum percent identity for inclusion of a solo LTR in a cluster (80.0)
soloLTRminLen			<num>	Minimum percent of LTR length participating in alignment for inclusion of LTR in a cluster (80.0)
soloLTRmaxEvalue		<num>	Maximum evalue allowed for blastn



Finding pairs of elements within clusters that have homologous flanking regions
-------------------------------------------------------------------------------
--rmhomoflank				Remove one of each pair of elements within each alignment (and therefore, each tree).
					(default OFF; Fixed ON when using --DTT)
--bpflank		<int>		Number of bases on either side of each element to search for homology. (default 500 bp)
--flank_evalue	<int|float>	E-value ceiling for considering blastn hits as evidence of homology. (default 1e-5)
--flank_pId		<int|float>	Minimum percent identity in blastn alignment to consider hit as evidence of homology. (default 70)
--flank_plencutoff	<int|float>	Minimum percentage of flanking region required to participate in alignment to consider
					blastn hit as evidence of homology. (default 70)
Phylogenetic analysis
---------------------
--phylo				##### not implemented yet
--nosmalls				Do not combine and perform phylogentic analyses on clusters smaller than --min_clust_size.
--DTT					Turns on --rmhomoflank, --convert_to_ultrametric, and --auto_outgroup. Generates and attempts to run
					Rscript that generates a DTT (LTT) plot for each cluster for which rooting is possible.
--bootstrap_reps		<int>	Number of replicates to generate for bootstrapping (default 100)
--convert_to_ultrametric		Convert trees to ultrametric using PATHd8. (default OFF; ON when using --DTT)
--auto_outgroup			Pick an outgroup automatically:
      				 -The outgroup shall be a random element from cluster k where cluster k is the largest of the clusters
      				  that is not j if j is the first cluster then the next smallest cluster is k if there is no other
      				  cluster, no outgroup is used.

## APPENDIX A. Global MAFFT options
This step has been the limiting process in my experience. It can be sped up by reducing the number of iterations performed for each alignment and by reducing the maximum number of elements for classiying elements as medium and small clusters. MAFFT exhausted 256 Gb RAM with ~2.7k seqs of length >5kb. Depending on resources available to you, you may need to cap the size of clusters to align using `--mafft_largeAln_maxclustsize`. Default is to not align clusters with >1000 elements. The MAFFT algorthim FFT-NS-2 is used for small and medium clusters and FFT_NS-1 for large clusters, which is much more inaccurate.

###### Available options
--maxiterate_small_clusters (30)	MAFFT iterations. More will improve alignment quality.
--maxiterate_medium_clusters (3)	MAFFT iterations. More will improve alignment quality.
--mafft_smallAln_maxclustsize (50)	Max elements to consider a cluster small.
--mafft_mediumAln_maxclustsize (500)	Max elements to consider a cluster medium.
--mafft_largeAln_maxclustsize (1000)	Max elements to consider a cluster large. Clusters larger than this will not be aligned.
--min_clust_size (7)			Do not align clusters smaller than this.
--nosmalls				Do not combine and assemble all clusters smaller than --min_clust_size

## APPENDIX B. Other global options
--keep_files				Default=no. Removes some large intermediate files, including raw
--output_dir		<path>	Output directory. Default is "PhyLTR.output
--logfile			<path>  Path to where log file is written (default <output_dir>/log.txt)
