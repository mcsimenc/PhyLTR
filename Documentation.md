PhyLTR is a software pipeline built from open source software. The main program is written in Python 3. Many of the routines in PhyLTR are parallelized, but it is not written for MPI so parallel components can only run on CPUs that share memory, i.e. on a single node.

Installation

`PhyLTR` = The location of the main PhyLTR repository on your computer

PhyLTR's Python 3 code runs external programs which it finds on your computer according to the paths in the text file
```
PhyLTR/CONFIG
```
The `CONFIG` file has format `key=path` where `key` needs to be exactly as shown below and `path` is either `file` or `directory` containing the program or the directory containing the program shown after the `#`

```
bedtools=file # bedtools executable
mafft=file # mafft executable
fasttree=file # fasttree executable
trimal=file # trimal executable
jmodeltest2=file # jModelTest.jar
genometools=file # gt executable
geneconv=file # geneconv executable
paup=file # paup executable
rscript=file # Rscript executable
perl=file # perl executable
circos=file # circos executable
pathd8=file #PATHd8 executable
getorf=file # EMBOSS getorf executable
phylip=directory # the bin/ directory in the PHYLIP installation
hmmer=directory # the binaries/ directory in the HMMER3 installation
blast=directory # the bin/ directory in the BLAST+ installation
mcl=directory # the bin/ directory in the MCL installation
```

PhyLTR uses two databases for classying LTR RTs by homology using HMMER3 (`nhmmer`)and BLAST+ (`tblastx`), Dfam and Repbase. Dfam used to be available for download, but the page doesn't seem to be up as of June 24, 2018. The files required by PhyLTR are

```
PhyLTR/RepeatDatabases/Dfam/Dfam_ERV_LTR.hmm
PhyLTR/RepeatDatabases/Dfam/Dfam_ERV_LTR.SF
PhyLTR/RepeatDatabases/Dfam/Dfam_ERV_LTR.list
PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.fasta
PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.SF
PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.list
```

`Dfam_ERV_LTR.hmm`  needs to contain the elements from the full Dfam DB, Dfam.hmm, that are annotated as Class: LTR.
`Dfam_ERV_LTR.SF` needs be a two-column file with the ID of a given element in the first column and its Superfamily annotation in the Dfam database. This file can be made using the following PhyLTR Python 3 script like this: `PhyLTR/scripts/DfamHMM2SuperFamTable.py < Dfam_ERV_LTR.hmm > Dfam_ERV_LTR.SF`
`Dfam_ERV_LTR.list` needs to be just the first column of the `Dfam_ERV_LTR.SF` file
`Repbase_ERV_LTR.fasta` needs to be the concatenation of all LTR retrotransposon and Enogenous Retrovirus features from Repbase. An account with GIRI is required. This is how I did it:

1. Get an account with GIRI
2. Go to http://www.girinst.org/repbase/update/browse.php
3. Select LTR Retrotransposon from the Repeat class dropdown list.
4. Select FASTA from the Output format drop down list.
5. Click the Download button, sign in, and download the text page that opens.
6. Repeat steps 3-5 but select Endogenous Retrovirus from the Repeat class dropdown list.
7. Concatenate the files using `cat` in BASH, name the concatenation `Repbase_ERV_LTR.fasta`, and put it in `PhyLTR/RepeatDatabases/Repbase`

`Repbase_ERV_LTR.SF` needs to be a two-column file with the ID of a given element in the first column and its Superfamily annotation in the Repbase database. To make this file, download the LTR retrotransposon and Endogenous Retrovirus features from Repbase just like shown above for `Repbase_ERV_LTR.fasta` but select IG format instead of FASTA format. Then concatenate them and use the PhyLTR Python 3 script to create the two-colum file like this: `PhyLTR/scripts/RepbaseIG2superfamilies.py < Repbase.LTR-ERV-concatenated.IG > Repbase_ERV_LTR.SF`
`Repbase_ERV_LTR.list` needs to be just the first column of the `Repbase_ERV_LTR.SF` file
