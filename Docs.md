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
