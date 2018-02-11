# Phylogenetic analysis of long terminal repeat retrotransposons
## PhALTR
  
### Dependencies
Program			Version tested with PhALTR	Process using program  
----------------------------------------------------------------------------------  
BEDTools			2.26.0			many  
MAFFT				7.310			many  
SAMtools			1.3			making FASTA index  
FastTree			2.1.10			inferring phylogenies  
trimAl				1.2			many  
jModeltest2			2.1.9			model selection for LTR divergence  
GenomeTools			1.5.9			LTRharvest, LTRdigest, suffixerator, gff3  
HMMER3				3.1b2			classifying to superfamily. must have >=3.1  
BLAST+				2.2.31+			many  
MCL				14-137			Clustering using MCL  
GENECONV			1.81a			finding gene conversion  
PAUP*				4a159			LTR divergence estimation  
PHYLIP				3.697			bootstrapping (seqboot)  
Rscript				3.3.3			running R scripts  
Perl				5.24.1			running bootstrapping scripts  
PATHd8							converting to ultrametric trees  
EMBOSS (getorf)			6.6.0			finding long ORFs in elements  
  
  
R packages needed (tested with R 3.3.3 on a Mac)  
---------------------------------------------------------------------------------  
For divergence plots					ggplot2  
For LTT plots						ape  
							phangorn  
							phytools  
