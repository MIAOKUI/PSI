## Introduction
Alternative splicing generates thousands of isoforms to multiply the diversity of the proteome and non-coding RNA transcripts. Current high-throughput sequencing technologies are capable of capturing this post-transcriptional processes with great sensitvity but only yield rather short read lengths. It thus becomes very challenging to predict the entire exon composition of RNA molecules. To qualitatively and quantitatively assess alternative splicing on a genome-wide level without prior knowledge, we developed this pipeline that quantifies the the usage of unknown or known splice sites as well as the expression level of exons to calculate a splicing index called Percent Spliced In (PSI) from RNA-seq and Ribo-seq data. Ranging from 0% for non-included exons to 100% for constitutive exons, it indicates how efficiently an exon sequence is spliced into the pool of transcripts transcribed from a genetic locus. This technology is capable of identifying and quantifying novel and complex exon skipping, alternative 5 and 3 splice sites and mutually exclusive splicing events. 

For detail information please refer to our publication:
* [Schafer S, Miao K, Benson CC, Heinig M, Cook SA & Hubner N (2015) Alternative Splicing Signatures in RNA-seq Data: Percent Spliced in (PSI).Curr Protoc Hum Genet 87, 11.16.1-11.16.14.](http://onlinelibrary.wiley.com/doi/10.1002/0471142905.hg1116s87/abstract)

## Quick Start
### Dependency and preparation
1.  Runing environment 
* Linux based system(linux, mac, unix)
 
2.  For python psi script 
* [Python 2.7](https://www.python.org)
* [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html)

3.  For shell psi script
* [SAMtools](http://samtools.sourceforge.net/)
* [BEDtools](http://bedtools.readthedocs.org/en/latest/) 

4. # For Visulization
* R
* R package dependency: reshape2, ggplot2, rtracklayer, gridExtra, grid, plyr

5.  About the input file. 
* junction.bed file and bam file can be generated directly  from tophat/tophat2
* flattened_gff_file.gtf can be generated from dexseq_prepare_annotation.py script included in [DEXseq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html)

## General Usage
### For Shell PSI script:
```bash
./psi_shell/PSI.sh StartPSIStrictFilter flattened_gff_file.gtf readlen alignment_file.bam junctions.bed baseName
```
#### Output
```bash
ENSRNOG00000033734:001	38	46548	0	1
ENSRNOG00000033734:002	11	47993	0	1
ENSRNOG00000033734:003	12	7772	32484	0.175863
ENSRNOG00000033734:004	30	54112	2843	0.935926
ENSRNOG00000033734:005	87	113977	0	1
``` 
*  Exonic part ID
*  Length of exnic part
*  Inclusion count
*  Exclusion count 
*  PSI value 


### For python PSI script: 
1. Counting inclusion count

 ```bash 
 ./psi_python/dexseq_count.py -p no -s no -r pos -f bam flattened_gff_file.gtf  alignment_file.bam basename.inclusion
 ```
2. Counting exclusion count

 ```bash 
 ./psi_python/exclusion_count.py  flattened_gff_file.gtf junctions.bed basename.exclusion
 ```
3. Calculate the psi 

 ```bash 
 ./psi_python/psi_calculation.py -l 100  flattened_gff_file.gtf basename.inclusion basename.exclusion basename
 ```
**NOTICE**: step 1 and step two can run parallel. 

#### Output
```shell 
ENSRNOG00000033734:001	1.0	1	44487	0	13	48872971	48873009	+	38
ENSRNOG00000033734:002	1.0	1	45933	0	13	48873964	48873975	+	11
ENSRNOG00000033734:003	0.170570528295	1	7490	32484	13	48874085	48874097	+	12
ENSRNOG00000033734:004	0.933449251441	1	51960	2843	13	48876463	48876493	+	30
```
*  Exonic part ID
*  PSI value
*  Low inclusion filter tag, 1 pass, 0.2 failed 
*  Inclusion raw count from dexseq_count.py
*  Exclusion raw count from exlcusion_count.py 
*  Chromosome ID
*  Start coordinate of exonic part 
*  End coordinate of exonic part 
*  Strand tag
*  Length of exonic part

If the set **--detail** **yes**, then only showing first three fields like following.
```shell
ENSRNOG00000033734:001	1.0	1
ENSRNOG00000033734:002	1.0	1
ENSRNOG00000033734:003	0.170570528295	1
ENSRNOG00000033734:004	0.933449251441	1
``` 

### Visulization
We also provide R and ggplot2 base visulization script for psi result visuliztion, located on psi_visulization directory. 
*  psi_plot_vis.r: including all the function used for ploting PSI result
*  psi.plot.r: psi plot wrapper, which import functions from psi_plot_vis.r, user can do some configuration on it to fit local environment. 

#### Configuration and running. 
For simple visulization, user only need to run psi.plot.r. Before running, some options need to be changed to fit local environment. 
  R package dependency: reshape2, ggplot2, rtracklayer, gridExtra, grid, plyr

1. Source psi_plot_vis.r script
```R
source('./psi_plot_vis.r') ## Here should be the correct path to the psi_plot_vis.r script
```
2. Provide sample experiment table like DESeq style
```R
sampleTablePath = './sampleTable.csv' ## Here should be the path to sampleTable, user can refer to the testing data to create own one. 
```
3. Provide psi result folder
```R
psiFolder = './psi_files/'
```
4. Provide the path to exonic part matrix which can be created by dexseq_prepare_annotation.py script from DEXSeq packages. 
```R
exonicPartMatrixPath ="./homo_sapiens.GRCh38.78_exonic.gtf"
```
5. gene id 
```R
gene_id = "ENSG00000099810+ENSG00000264545+ENSG00000274055+ENSG00000264801+ENSG00000240498"
```
6. Optional if aggregate the PSI by experiment group. 
```R
as.group = TRUE
```

8 For ploting a single genes, refering following code. 
```R
psiTable <- importPSI(sampleTablePath, psiFolder)
exonicPartMatrix <- import(exonicPartMatrixPath)
geneAnnot <- getGenesAnot(gene_id, exonicPartMatrix)
genePsi <- getGenesPsi(gene_id = gene_id,
                       psiTable = psiTable, 
                       sampleTablePath = sampleTablePath,
                       as.group = TRUE)
## ploting the psi 
plot(psiPlot(geneAnnot, gene_id, genePsi))
```

7. For batch ploting a gene list, user can refering following code. 
```R
## Gene list
gene_list =readLines("./gene_list.txt")
psiTable <- importPSI(sampleTablePath, psiFolder)
for(sample in  gene_list ){
  geneAnnot <- getGenesAnot(sample, exonicPartMatrix)
  genePsi <- getGenesPsi(gene_id = sample,
                         psiTable = psiTable, 
                         sampleTablePath = sampleTablePath,
                         as.group = TRUE)
  path.out <- file.path("./plots/")
  if(!file.exists(path.out)){dir.create(path.out)}
  pdf(file=paste("./plots/",sample,'.pdf',sep = ''),width=40,height=30 ) 
  plot(psiPlot(geneAnnot, sample, genePsi))
  dev.off()
```


