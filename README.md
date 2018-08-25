### For installation instructions read INSTALL.md

# PhyLTR: de novo LTR retrotransposon annotation, classification, and analysis using phylogenetic tools.

PhyLTR is a software pipeline built from open source software. The main program is written in Python 3. Many of the routines in PhyLTR are parallelized, but it is not written for MPI so parallel components can only run on CPUs that share memory, i.e. on a single node.

As the pipeline runs, paths to intermediate results like alignments are stored in the file PhyLTR/status. If the execution is interrupted, this file is used to allow PhyLTR to resume more or less where it left off.

## Default settings

If phyltr is run without any flags specifying a task, all tasks are run. The following two calls are equivalent. The processes specified by the flags in the second call are explained below with additional optional flags.

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

## All settings. Defaults are in parentheses

#### 1. Identifying candidate LTR-R loci with LTRharvest

###### Turn on using `--ltrharvest`

###### The following options are available and explained in the LTRharvest documentation.
* --minlenltr (100)
* --maxlenltr (1000)
* --mindistltr (1000)
* --maxdistltr (15000)
* --similar (85.0)
* --vic (60)
* --mintsd (4)
* --maxtsd (20)
* --xdrop	(5)
* --mat (2)
* --mis (-2)
* --insi (-3)
* --del (-3)

	  Option			    ArgType	       Default
	----------------------------------------------------------------
	--fasta				    <path>		NONE
	--procs				    <int>		1
	--keep_files			    BINARY		OFF
	--output_dir			    <path>		PhyLTR.output
	--ltrharvest			    BINARY		OFF
	--minlenltr			    <int>		100
	--maxlenltr			    <int>		1000
	--mindistltr			    <int>		1000
	--maxdistltr			    <int>		15000
	--similar			    <num>		85.0
	--vic				    <int>		60
	--mintsd			    <int>		4
	--maxtsd			    <int>		20
	--xdrop				    <int>		5
	--mat				    <int>		2
	--mis				    <int>		-2
	--ins				    <int>		-3
	--del				    <int>		-3
	--ltrdigest			    BINARY		OFF
	--ltrdigest_hmms		    <path>		/home/joshd/scripts/PhyLTR/LTRdigest_HMMs/hmms
	--classify			    BINARY		OFF
	--classify_dfam			    BINARY		OFF
	--classify_repbase		    BINARY		OFF
	--nhmmer_reporting_evalue	    <num>		10
	--nhmmer_inclusion_evalue	    <num>		1e-2
	--repbase_tblastx_evalue	    <num>		1e-5
	--keep_conflicting_classificaitons  BINARY		OFF
	--keep_no_classifications	    BINARY		OFF
	--min_clust_size		    <int>		7
	--wicker			    BINARY		OFF
	--wicker_pId			    <num>		80
	--wicker_minLen			    <int>		80
	--wicker_pAln			    <num>		80
	--wicker_no_internals		    BINARY		OFF
	--wicker_no_ltrs		    BINARY		OFF
	--mcl				    BINARY		MCL
	--I				    <num>		6
	--nosmalls			    BINARY		OFF
	--geneconvltrs			    BINARY		OFF
	--geneconv_g			    <str>		g0,g1,g2
	--ltrdivergence			    BINARY		OFF
	--remove_GC_from_modeltest_aln	    BINARY		OFF
	--modeltest_criterion		    <str>		BIC
	--gc_ltr_div_scaling		    <int>		1
	--maxiterate_small_clusters	    <int>		20
	--maxiterate_medium_clusters	    <int>		3
	--mafft_smallAln_maxclustsize	    <int>	 	50
	--mafft_mediumAln_maxclustsize	    <int>		500
	--mafft_largeAln_maxclustsize  	    <int>		1000
	--geneconvclusters		    BINARY		OFF
	--DTT				    BINARY		OFF
	--phylo				    BINARY		OFF
	--bootstrap_reps		    <int>		100
	--bpflank			    <int>		500
	--flank_evalue			    <num>		1e-5
	--flank_pId			    <num>		70
	--flank_plencutoff		    <num>		70
	--min_orf_len			    <int>		300
	--soloLTRsearch			    BINARY		OFF
	--soloLTRminPid			    <num>		80.0
	--soloLTRminLen			    <num>		80.0
	--soloLTRmaxEvalue		    <num>		1e-3

	
