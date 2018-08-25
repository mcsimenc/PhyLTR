# How to install PhyLTR: Three steps,

1. Clone repository
2. Install dependencies
3. Download databases
4. Optionally add pHMMs for domain annotation

## Clone Repository

## Install dependencies

#### 1. Install [Python 3](https://www.python.org/) and [Biopython](https://biopython.org/)

#### 2. Install these however you can, then add the program paths to the CONFIG file in the PhyLTR root directory.
Parts of PhyLTR require only certain dependencies. See README.md for an explanation of dependency requirements for each process.

* [BEDtools](https://bedtools.readthedocs.io/en/latest/)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [FastTree2](http://www.microbesonline.org/fasttree/)
* [trimAl](http://trimal.cgenomics.org/)
* [jModelTest2](https://github.com/ddarriba/jmodeltest2/)
* [GenomeTools](http://genometools.org/)
* [GENECONV](https://www.math.wustl.edu/~sawyer/geneconv/)
* [PAUP*](https://paup.phylosolutions.com/)
* [Circos](http://circos.ca/)
* [PATHd8](https://www2.math.su.se/PATHd8/)
* [EMBOSS](http://emboss.sourceforge.net/)
* [PHYLIP](http://evolution.genetics.washington.edu/phylip.html)
* [HMMER3](http://www.hmmer.org/)
* [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
* [MCL](https://micans.org/mcl/)

#### 3. Edit CONFIG file and add dependency paths

The CONFIG file has format: `key=path` where `key` needs to be exactly as shown below and `path` is expected to point to either the `file` of the program itself or the `directory` containing the program, depending on the dependency.

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

## Download databases: Dfam and Repbase

### Dfam

##### 1.Download http://dfam.org/web_download/Release/Dfam_2.0/Dfam.hmm.gz and unpack it

##### 2. Run: `PhyLTR/scripts/DfamExtractLTRelements.py < Dfam.hmm > Dfam_ERV_LTR.hmm`

##### 3. Run: `PhyLTR/scripts/DfamHMM2SuperFamTable.py < Dfam_ERV_LTR.hmm > Dfam_ERV_LTR.SF`

##### 4. Run: `cut -f1  < Dfam_ERV_LTR.SF > Dfam_ERV_LTR.list`

##### 5. Move the files from B,C,D to the following locations:
```
PhyLTR/RepeatDatabases/Dfam/Dfam_ERV_LTR.hmm
PhyLTR/RepeatDatabases/Dfam/Dfam_ERV_LTR.SF
PhyLTR/RepeatDatabases/Dfam/Dfam_ERV_LTR.list
```

### Repbase

##### 1. Get an account with GIRI
1. Go to http://www.girinst.org/repbase/update/browse.php
2. Select LTR Retrotransposon from the Repeat class dropdown list.
3. Select FASTA from the Output format drop down list.
4. Click the Download button, sign in, and download the text page that opens.
5. Repeat steps 2-4 but select Endogenous Retrovirus from the Repeat class dropdown list.
6. Run: `cat <LTR.fa> <ERV.fa> >> Repbase_ERV_LTR.fasta`
7. Move the new file from 6 to: `PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.fasta`

##### 2. Run: `PhyLTR/scripts/RepbaseIG2superfamilies.py < Repbase.LTR-ERV-concatenated.IG > Repbase_ERV_LTR.SF`

##### 3. Run: `cut -f1  < Repbase_ERV_LTR.SF > Repbase_ERV_LTR.list`

##### 4. Move the files from B,C to the following locations:
```
PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.fasta
PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.SF
PhyLTR/RepeatDatabases/Repbase/Repbase_ERV_LTR.list
```

## Add pHMMs for domain annotation (optional)
##### Append any HMMs you want to include to `PhyLTR/LTRdigest_HMMs/hmm`
The version included in repository contains pHMMs for TE-related domains from Pfam and from gydb.org, downloaded Summer 2018.
