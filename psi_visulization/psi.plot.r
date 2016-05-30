## library require to pre installed 
## reshape2 ggplot2 rtracklayer gridExtra grid plyr


## Import plotting main functions 
source("./psi_plot_vis.r")

##  the location of sampleTable, which provide sample name and grouping factor information
##  Which can be mordified based on need
sampleTablePath = "./sampleTable.csv"   

## the path to the psi result fold, the result file should be output directortly from psi.sh script 
## noted that the  end '/' shoubd be included
psiFolder ="./troponint_psi/"

## the path to the exonic part matrix
exonicPartMatrixPath ="./genes.flat.gtf"

## the gene_id which you want to plot 
gene_id = "ENSRNOG00000033734"

## if only plot the mean of each group 
as.group = TRUE


## calling the following function sequentilly 
psiTable <- importPSI(sampleTablePath, psiFolder)
geneAnnot <- getGenesAnot(gene_id, exonicPartMatrixPath)
genePsi <- getGenesPsi(gene_id = gene_id,
                       psiTable = psiTable, 
                       sampleTablePath = sampleTablePath,
                       as.group = TRUE)





## ploting the psi 
plot(psiPlot(geneAnnot, gene_id, genePsi))
