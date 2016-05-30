require(reshape2)
require(ggplot2)
require(rtracklayer)
require(gridExtra)
require(grid)
require(plyr)


importPSI <- function(sampleTableName,psiFolder){
        sampleTable <- read.csv(sampleTableName, stringsAsFactors = FALSE)
        psiTable <- list() 
        for( i in 1:nrow(sampleTable)){
                tmpPsiPath <- paste0(psiFolder,sampleTable[i,"FileName"])
                print(tmpPsiPath)
                tmpPsiFile <-  read.table(tmpPsiPath,header = TRUE,stringsAsFactors = FALSE)
                psiTable[[sampleTable[i,"SampleName"]]] <- tmpPsiFile[,5]
        }
        psiTable <- do.call("cbind", psiTable)
        rownames(psiTable) <- tmpPsiFile$exon_ID 
        return(psiTable)
}


getGenesAnot <- function(gene_id, exonicPartMatrixPath){
        exonicPartMatrix <- import(exonicPartMatrixPath)
        geneAnnot <- exonicPartMatrix[ exonicPartMatrix$gene_id == gene_id & !is.na(exonicPartMatrix$transcripts),]
        
        return(geneAnnot)
        
}


getGenesPsi <- function(gene_id, psiTable, sampleTablePath, as.group){
        sampleTable <- read.csv(sampleTablePath, stringsAsFactors = FALSE)
        vGroup <- sampleTable$Group
        names(vGroup) <- sampleTable$SampleName
        genesPsi <- psiTable[ grep(gene_id, rownames(psiTable)),]
        if(as.group == TRUE){ 
                genePsiM <- melt(genesPsi)
                names(genePsiM) <- c("exonID", "SampleName", "PSI")
                genePsiM$group <- vGroup[ as.character(genePsiM$SampleName) ]
                genePsiSum <- ddply(genePsiM, .(exonID, group), 
                                    .fun = function(x){
                                            psi_mean = mean(x$PSI)
                                    })
                names(genePsiSum)[3] <- "PSI"
                genePsiSumLong <- reshape(genePsiSum, 
                                          idvar ="exonID",
                                          direction = 'wide', 
                                          timevar = 'group')
                rownames(genePsiSumLong) <- genePsiSumLong$exonID 
                genePsiSumLong$exonID <- NULL
                return(genePsiSumLong)
        }else{
                return(genesPsi)    
        }    
}




reConGrange <- function(inputGr,intronSize){
        reConIntron <- function(inputGr,intronSize){
                intronList <- list()
                for( i in 1:(length(inputGr)-1)){
                        if(start(inputGr)[i+1] - end(inputGr)[i] == 1){
                                intronList[i] <-  0
                        }else intronList[i] <- intronSize
                }
                return(do.call(c,intronList))
        }
        reConInterval <- function(exonWid,intronList){
                start <- list(1)  
                end <- list(start[[1]] + exonWid[1]-1) 
                for( i in 2:length(exonWid)){
                        start[i] <- end[[i-1]] + intronList[i-1] 
                        end[i]  <- start[[i]] + exonWid[i]-1
                }
                start <- do.call(c,start)
                end <- do.call(c,end)
                return(IRanges(start,end))
        }
        
        strandTag <- unique(strand(inputGr))
        if(strandTag == "-"){
                ord <- length(inputGr):1
                ordIntron <- (length(inputGr)-1) : 1
        }else{
                ord <- 1:length(inputGr)
                ordIntron <- 1:(length(inputGr)-1)}
        intronList <- reConIntron(inputGr,intronSize)
        exonWid <- width(inputGr)[ord]
        seqnames <- seqnames(inputGr)
        ranges <- reConInterval(exonWid,intronList[ordIntron])
        newGrange <- GRanges(seqnames,ranges,strand(inputGr)) 
        coordinate <- paste(seqnames(inputGr),paste(start(inputGr),end(inputGr),sep="-"),sep=":")
        exonic <- paste("E",c(1:length(inputGr)),sep=":")
        return(list(newGrange=newGrange,
                    intronList=intronList[ordIntron],
                    coordinate=coordinate,
                    exonic=exonic,
                    exonicOrd=exonic[ord]))
}

addValue <- function(newGrange,valueDF){
        strandTag <- unique(strand(newGrange))
        if(strandTag == "-"){
                mcols(newGrange) <- valueDF[length(newGrange):1,]
        }else 
                mcols(newGrange) <- valueDF 
        return(newGrange)
}



## geneAnot is the annotation part for certain genes 
## gene is gene ID 
## genePsi is Pis table from Psi Calculation result for certain genes 
psiPlot <- function(geneAnot,gene_id,genePsi){
        title<- gene_id 
        if(length(geneAnot) <2) return() 
        newGr <- reConGrange(geneAnot, 100)
        newInList <- newGr[["intronList"]]
        newGrange <- addValue(newGr$newGrange,genePsi)
        df <- as.data.frame(newGrange, stringsAsFactors = FALSE )

        begin <- seq(1,width(range(newGrange)),by=width(range(newGrange))/length(newGrange))
        wid <- width(range(newGrange))/length(newGrange) * 0.8
        gap <- width(range(newGrange))/length(newGrange) * 0.2
        
        df.segment <- data.frame(df[,c(-1:-5)], begin=begin,
                                 endFragment = begin + wid, 
                                 endGap =begin + wid + gap)

        df.shift <- df[c(2:nrow(df),NA),c(-1:-5)]
      
        
        df.shift.melt <- melt(df.shift,value.name="PSI_next")
        df.segment.melt <- melt(data=df.segment,
                                id.var = c("begin","endFragment","endGap"),
                                value.name="PSI",
                                variable.name = "Group")

        names(df.segment.melt)[c(4,5)] <- c("Group","PSI") 

        names(df.shift.melt)[2] <- "PSI_next"

        df.merge <- cbind(df.segment.melt, PSI_next = df.shift.melt$PSI_next)

        
        df.link <- data.frame(beginPoint = df$start + (df$width)/2-1 ,
                              endPoint=df.segment$begin + wid/2 -1 ) 

        upPlot <- ggplot(df.merge) + 
                geom_segment(aes(x=begin,y=PSI,xend=endFragment,yend=PSI,colour = Group),size = 1.5) + 
                geom_segment(aes(x=endFragment,xend=endGap,y=PSI,yend=PSI_next, colour = Group),linetype=2,alpha=0.6) + 
                xlim(1,max(df.merge$endGap)) +
                scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1),label = seq(0,100,by=10))+
                ggtitle(title)+ 
                theme(plot.margin=unit(c(0,-0.6,0,0), "cm"),
                      panel.background = element_blank(),
                      axis.title.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(),
                      panel.grid.major = element_line(linetype= "dashed",colour = "gray"),
                      axis.line=element_line(colour = "black",lineend="square"),
                      axis.line.x = element_blank(),
                      axis.text.y = element_text(colour = "black",size = 20,face = "bold"),
                      axis.title = element_text(size = 25, face = "bold"),
                      title = element_text(size = 30, face = "bold"),
                      legend.background = element_rect(colour= "black"),
                      legend.title = element_blank(),
                      legend.text = element_text(size = 11, face = "bold"),
                      legend.key.size = unit(0.8,"cm"))


        
        medPlot <- ggplot(df.link) + 
                geom_segment(aes(x=beginPoint,xend=endPoint, y=0, yend=1),arrow = arrow(ends="first",length = unit(0.1, "inches"),type="closed")) + 
                xlim(1,max(df.merge$endGap)) +
                theme(panel.background=element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin=unit(c(-0.6,-0.6,0,0), "cm"),
                      axis.text = element_blank(),
                      axis.title = element_blank())
      

        textPlot <- ggplot(df.link,environment = environment()) + 
                geom_text(aes(x=endPoint,y=0,label = newGr$exonicOrd,fontface = "bold"),angle = 90)  +  
                xlim(1,max(df.merge$endGap)) +
                theme(panel.background=element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin=unit(c(-1,-1,0,0), "cm"),
                      axis.text = element_blank(),
                      axis.title = element_blank())
        
        bottonPlot <- ggplot(df,environment = environment()) + 
                geom_rect(aes(xmin=start,xmax=end,ymin=0.2,ymax=1.2),fill ="gray",colour="black") +  
                geom_segment(aes(x=end,xend=end + c(newInList,0), y=0.7,yend=0.7)) +
                geom_segment(aes(x=1,xend=1000,y=0.15,yend=0.15)) + 
                geom_text(aes(x=500,y=0.1, label = "1Kb")) +
                xlim(1,max(df.merge$endGap)) +
                theme(panel.background=element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin=unit(c(-0.6,0,0,0), "cm"),
                      axis.text = element_blank(),
                      axis.title= element_blank(),
                      axis.line = element_blank())
        
        cor_dat <- data.frame(x = rep(1.5,length(newGr$newGrange)),
                              y=seq(from=3,to=1,length.out=length(newGr$newGrange)),
                              label = paste(newGr$exonic,newGr$coordinate,sep = "  "))
    
        
        corPlot <- ggplot(cor_dat) + 
                geom_text(aes(x=x,y=y,label=label,fontface="bold")) + 
                theme(panel.background=element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin=unit(c(2,1,1,1), "cm"),
                      axis.text = element_blank(),
                      axis.title= element_blank(),
                      axis.line = element_blank())

        g.upPlot <- ggplotGrob(upPlot)
        g.textPlot <- ggplotGrob(textPlot)
        g.medPlot <- ggplotGrob(medPlot)
        g.bottonPlot <- ggplotGrob(bottonPlot)
        
        cleg <- g.upPlot$grobs[[8]]
        aleg <- ggplotGrob(corPlot)
        
        
        upPlot <- upPlot  + theme(legend.position = "none")
        g.upPlot <- ggplotGrob(upPlot)
        
        
        maxWidth <- unit.pmax(g.upPlot$widths,g.textPlot$widths,g.medPlot$widths,g.bottonPlot$width)
        g.upPlot$widths <- maxWidth
        g.textPlot$widths <- maxWidth
        g.medPlot$widths <- maxWidth
        g.bottonPlot$widths <- maxWidth
        
        total.notLangend <- arrangeGrob(g.upPlot, 
                                        g.textPlot,
                                        g.medPlot,
                                        g.bottonPlot,
                                        heights = c(5.8/10,0.1/10,1/10,2/10),
                                        nrow = 4)

        total.legend <- arrangeGrob(cleg,aleg,nrow=2)
        
        total.grob <- arrangeGrob(total.notLangend,
                                  total.legend,nrow=1,  widths = c(9/10,2/10))

        
        return(total.grob)  

}
