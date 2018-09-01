# PhyLTR scripts
Some of these scripts are called during PhyLTR's execution and some are not, but all are standalone and contain helps.
## External scripts
#### Visualizing trees and diagrammatic elements
```
ete3_tree.py                  Requires a graphics device. Renders phylogenetic trees of LTR-Rs with
                              diagrams showing element ORF and domain content.
```
---
## For preparing the LTR-R databases
```
DfamExtractLTRelements.py     Script for preparing the Dfam databse for PhyLTR (see INSTALL.md)
DfamHMM2SuperFamTable.py      Script for preparing the Dfam databse for PhyLTR (see INSTALL.md)
RepbaseIG2superfamilies.py    Script for preparing the Repbase databse for PhyLTR (see INSTALL.md)

```
## Used in PhyLTR
```
best_blast_hit.py             Returns the highest scoring hit for each query in tabular BLAST output
CompareToBootstrap.pl         FastTree2 script for calculating bootstrap support values
gffAddAttr.py                 Adds new key-value pairs to the attributes field of a GFF3
hmmer_table2columns.py        Parses hmmer output (deprecated in favor of nhmmer)
nhmmer_table2columns.py       Parses nhmmer output, used in Dfam classification step
gff2circos-heatmap.py         Creates Circos-format file for heatmap features from GFF3 input
gff2circos-tile.py            Creates Circos-format file for tile features from GFF3 input
ideogramFromLengths.py        Creates Circos-format file for scaffolds
geneconv2links.py             Converts GENECONV output to Circos links format
```
---
### Circos-related
```
domainSearch.py               Searchers for homologs to LTR-R ORFs in a protein FASTA
```
