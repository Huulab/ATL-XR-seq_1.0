library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(stringi)
library(grid)
library(reshape2)


##CORRELATION================
library(ggplot2)
library(ggcorrplot)
library(corrplot)
library(dplyr)
library(RColorBrewer)
library(stringr)
i=0
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/correlation/")
sample <- list.files()[grep("OVER-end_3000_SpearmanCorr_readCounts.tab",list.files())] 


col0 = colorRampPalette(brewer.pal(11,"RdYlBu"))
col0 = colorRampPalette(c("#FC4E2A","#FEB24C","#9ECAE1","#4292C6","#2171B5"))
col0 = colorRampPalette(c("#A50026","#D73027","#F46D43","#FDAE61","#6BAED6","#4292C6","#2171B5","#08519C"))
col0 = colorRampPalette(brewer.pal(9,"YlOrRd"))
#col1 = colorRampPalette(c("white","#08306B"))

cgname <- function(n){
  ee <- c()
  for (x in n){
    
    dd <- str_c( c(unlist(str_split(x,"-"))[c(1,2,3)]),collapse = " ")
    if(unlist(str_split(dd,"XR"))[1]=="Lo"){
      dd <- str_c(c("ATL-XR",unlist(str_split(dd,"XR"))[2]),collapse = "-")
    }else{dd <- dd }
    ee <- c(ee,dd)
  }
  return(ee)
}


cor_matrix <- read.table(sprintf("%s",sample),sep ="\t",row.names = 1,
                         header = T,skip = 1,stringsAsFactors = F) %>% as.matrix()



colnames(cor_matrix) <- cgname(rownames(cor_matrix))
rownames(cor_matrix) <- cgname(rownames(cor_matrix))

##good
# p <- corrplot(cor_matrix,method="color",col.lim=c(0.5,1),col=col0(100),
#               outline=F,addCoef.col="black",number.cex=rel(0.7),
#               tl.col="black",cl.cex=rel(0.8),tl.offset=rel(0.3),cl.align.text = "c",
#               order = "hclust", addrect = 2 ,is.corr = T)

#9:9
#ggsave(plot = p,filename = "correlation_all.pdf",width=6,height = 6,device = "pdf")


col0 = colorRampPalette(brewer.pal(11,"RdYlBu"))
col0 = colorRampPalette(c("#FC4E2A","#FEB24C","#9ECAE1","#4292C6","#2171B5"))
col0 = colorRampPalette(c("#A50026","#D73027","#F46D43","#FDAE61","#6BAED6","#4292C6","#2171B5","#08519C"))
col0 = colorRampPalette(brewer.pal(9,"YlOrRd"))
col0 = colorRampPalette(c("#A50026","#D73027","#F46D43","#FDAE61","#6BAED6","#4292C6","#2171B5","#08519C"))

##red to purples
col0 = colorRampPalette(c("#F7FBFF", "#DEEBF7" ,"#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C"))
col0 = colorRampPalette(brewer.pal(9,"Blues"))
col0 = colorRampPalette(brewer.pal(11,"PiYG"))
##orange to blues
col0 = colorRampPalette(c("#E31A1C", "#FC4E2A" ,"#FD8D3C", "#FEB24C", "#FED976", "#FFEDA0", "#FFFFCC", 
                          "#F7FBFF" ,"#DEEBF7" ,"#C6DBEF", "#9ECAE1", "#6BAED6" ,"#4292C6", "#2171B5"))

cordata <- melt(cor_matrix)
p <- ggplot(data = cordata, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  geom_text(aes(Var2, Var1,label=sprintf("%.2f", value)),color="black",size=rel(6),parse = TRUE)+
  labs(x="" ,y="") +
  coord_fixed()+
  scale_fill_gradientn(colours = col0(100)
                       ,name="Spearman\nCorrelation"
                       ,limit=c(0.5,1))+
  theme_minimal()+ 
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.2),face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        strip.text = element_text(size=rel(1),face="bold"),
        strip.background = element_blank(),
        axis.text.x =element_text(size=rel(1.4),angle=45,hjust = 1,vjust=1),
        axis.text.y =element_text(size=rel(1.4),angle=0),
        axis.title=element_text(size=rel(1.3)),
        legend.title = element_text(size=rel(1.2),face = "bold"),
        legend.text = element_text(size=rel(0.8),angle=0),
        axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))
  
#ggsave(plot = p,filename = "correlation_all.pdf",width=8,height = 6,device = "pdf")
i=i+1
assign(paste("p", i, sep = ""),p)

cgname <- function(n){
  ee <- c()
  for (x in n){
    bb <- unlist(str_split(x,"-"))[3]
    cc <- unlist(str_split(bb,"_"))[2]
    if (cc=="plus"){
      cc <- "+"
    } else {
      cc <- "-"
    }
    dd <- str_c( c(unlist(str_split(x,"-"))[c(1,2)],unlist(str_split(bb,"_"))[1],cc),collapse = " ")
    ee <- c(ee,dd)
  }
  return(ee)
}



##oligo length===============
library(ggplot2)
library(stringr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/length/")
temp <- data.frame(lengths =seq(21,31),
                   distribu=rep(0,31-21+1))

samples <- list.files()[grep("_length_distribution.txt",list.files())] 
samples <- samples[grep("^XRseq",samples)]

for (sample in samples) {
  aa <- read.table(sprintf("%s",sample),sep ="\t")
  colnames(aa) <- c("lengths","amount")
  aa <- inner_join(aa,temp,by="lengths")
  
  sumcounts <- sum(aa$amount)
  aa$distribu <- aa$amount/sumcounts*100
  rownames(aa) <- aa$lengths
  aa <- subset(aa,select = -amount)
  aa$state <- rep("XR-seq",nrow(aa))
  
  bb <- read.table(sprintf("ATL-%s-rep1-markdup_length_distribution.txt",unlist(str_split(sample,"-rep1",simplify=T))[,1]),sep ="\t")
  colnames(bb) <- c("lengths","amount")
  bb <- inner_join(bb,temp,by="lengths")
  
  sumcounts <- sum(bb$amount)
  bb$distribu <- bb$amount/sumcounts*100
  rownames(bb) <- bb$lengths
  bb <- subset(bb,select = -amount)
  bb$state <- rep("ATL-XR-seq rep1",nrow(bb))
  
  cc <- read.table(sprintf("ATL-%s-rep2-markdup_length_distribution.txt",unlist(str_split(sample,"-rep1",simplify=T))[,1]),sep ="\t")
  colnames(cc) <- c("lengths","amount")
  cc <- inner_join(cc,temp,by="lengths")
  
  sumcounts <- sum(cc$amount)
  cc$distribu <- cc$amount/sumcounts*100
  rownames(cc) <- cc$lengths
  cc <- subset(cc,select = -amount)
  cc$state <- rep("ATL-XR-seq rep2",nrow(cc))
  
  
  length_distri <- bind_rows(bb,aa,cc)
  
  p <- ggplot(data=length_distri,aes(x=length_distri$lengths,y=length_distri$distribu)) +
    facet_wrap(.~ state)+
    geom_bar(stat="identity",fill="#08519C",width=0.85) +
    labs(x="excised oligo length (nt)", y="Read count (%)") +
    scale_x_continuous(breaks=c(21,26,31))+
    expand_limits(y=seq(0,25,by = 5))+
    ggtitle(sprintf("%s",unlist(str_split(sample,"-",simplify=T))[,2]))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.3),face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          strip.text = element_text(size=rel(1.5),face="bold"),
          strip.background = element_blank(),
          axis.text.x =element_text(size=rel(1.4),angle=0),
          axis.text.y =element_text(size=rel(1.2),angle=0),
          axis.title=element_text(size=rel(1.3),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))
  i=i+1
  assign(paste("p", i, sep = ""),p)
  ggsave(filename = sprintf("rep2_%s.pdf",unlist(str_split(sample,"-",simplify=T))[,2]),width=6,height = 4)
}



##dipyrimidine============
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/dipyrimidine/")

samples <- list.files()[grep("26nt_dinucleotideTable",list.files())]    ##
samples <- samples[grep("^XRseq",samples)]

#plotlist <- list()
for (sample in samples) {

  aa <- read.table(sprintf("%s",sample),sep ="",header = T)
  aa[c(6,8,14,16),] %>% .$X19 %>% sum() %>% print()
  rownames(aa) <- aa$kmer
  aa <- aa %>% filter(aa$kmer %in% c("TC","TT","CT","CC")) %>%
    subset(,select = -kmer) %>% as.matrix()  ##ctrl+shift+m
  aa <- as.data.frame(as.table(aa*100))
  aa$state <- rep("XR-seq",nrow(aa))
  
  bb <- read.table(sprintf("ATL-%s-rep1-markdup_26nt_dinucleotideTable.txt",unlist(str_split(sample,"-rep1",simplify=T))[,1]),sep ="\t",header = T)
  bb[c(6,8,14,16),] %>% .$X19 %>% sum() %>% print()
  rownames(bb) <- bb$kmer
  bb <- bb %>% filter(bb$kmer %in% c("TC","TT","CT","CC")) %>%
    subset(,select = -kmer) %>% as.matrix()  ##ctrl+shift+m
  bb <- as.data.frame(as.table(bb*100))
  bb$state <- rep("ATL-XR-seq rep1",nrow(bb))
  
  cc <- read.table(sprintf("ATL-%s-rep2-markdup_26nt_dinucleotideTable.txt",unlist(str_split(sample,"-rep1",simplify=T))[,1]),sep ="\t",header = T)
  cc[c(6,8,14,16),] %>% .$X19 %>% sum() %>% print()
  rownames(cc) <- cc$kmer
  cc <- cc %>% filter(cc$kmer %in% c("TC","TT","CT","CC")) %>%
    subset(,select = -kmer) %>% as.matrix()  ##ctrl+shift+m
  cc <- as.data.frame(as.table(cc*100))
  cc$state <- rep("ATL-XR-seq rep2",nrow(cc))
  

  dipyrimidine <- bind_rows(bb,aa,cc)
  dipyrimidine$Var2 <-  str_c(unlist(str_extract_all(dipyrimidine$Var2, "\\d+")),
                              as.numeric(unlist(str_extract_all(dipyrimidine$Var2, "\\d+")))+1,sep = "-")
  colnames(dipyrimidine) <- c("dipyrimidine","pos","Freq","state")
  dipyrimidine$pos <- factor(dipyrimidine$pos, levels=str_c(seq(1,25),seq(2,26),sep = "-"), ordered=TRUE)
  
  p <- ggplot(dipyrimidine,mapping = aes(pos,Freq,fill=dipyrimidine))+
    facet_wrap(.~ state)+
    geom_bar(stat='identity',position='stack') +
    labs(x = 'Position',y = 'Dinucleotide rate (%)') +
    coord_cartesian(ylim=c(0,100))+ 
    scale_fill_manual(values=c("TC" = "#E41A1C", "TT" = "#984EA3", "CT" = "#377EB8" , "CC" = "#A65628"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
          panel.border = element_blank(), 
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          
          strip.text = element_text(size=rel(1.5),face="bold"),
          strip.background = element_blank(),
          
          axis.text.x =element_text(size=rel(1.0),angle=90),
          axis.text.y =element_text(size=rel(1.4),angle=0),
          axis.title=element_text(size=rel(1.3),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
    )+
    ggtitle(unlist(str_split(sample,"-",simplify=T))[,2])
  i=i+1
  assign(paste("p", i, sep = ""),p)
  #plotlist[[i]] <- p
  ggsave(filename = sprintf("%s_rep2.pdf",unlist(str_split(sample,"-",simplify=T))[,2]),width=8,height = 5)
}

##compare dipyrimidine============
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/dipyrimidine/")

samples <- list.files()[grep("26nt_dinucleotideTable",list.files())]    ##
samples <- samples[grep("^XRseq",samples)]
nub=0
#plotlist <- list()
for (sample in samples) {
  nub=nub+1
  aa <- read.table(sprintf("%s",sample),sep ="",header = T)
  aa[c(6,8,14,16),] %>% .$X20 %>% sum() %>% print()
  rownames(aa) <- aa$kmer
  aa <- aa %>% filter(aa$kmer %in% c("TC","TT","CT","CC")) %>%
    subset(,select = -kmer) %>% as.matrix() %>% .[,seq(18,23)] ##ctrl+shift+m
  for (n in 1:ncol(aa)){aa[,n] <- aa[,n]/sum(aa[,n])}
  
  
  bb <- read.table(sprintf("ATL-%s-rep1-markdup_26nt_dinucleotideTable.txt",unlist(str_split(sample,"-rep1",simplify=T))[,1]),sep ="\t",header = T)
  bb[c(6,8,14,16),] %>% .$X20 %>% sum() %>% print()
  rownames(bb) <- bb$kmer
  bb <- bb %>% filter(bb$kmer %in% c("TC","TT","CT","CC")) %>%
    subset(,select = -kmer) %>% as.matrix() %>% .[,seq(18,23)] ##ctrl+shift+m
  for (n in 1:ncol(bb)){bb[,n] <- bb[,n]/sum(bb[,n])}
  
  dipyrimidine <- log2(bb/aa) 
  
  
  
  dipyrimidine <- as.data.frame(as.table(dipyrimidine))
  
  dipyrimidine$Var2 <-  str_c(unlist(str_extract_all(dipyrimidine$Var2, "\\d+")),
                              as.numeric(unlist(str_extract_all(dipyrimidine$Var2, "\\d+")))+1,sep = "-")
  colnames(dipyrimidine) <- c("dipyrimidine","pos","Freq")
  dipyrimidine$pos <- factor(dipyrimidine$pos, levels=str_c(seq(1,25),seq(2,26),sep = "-"), ordered=TRUE)
  damage <- unlist(str_split(sample,"-",simplify=T))[,2]
  pic <- ggplot(dipyrimidine,mapping = aes(pos,Freq,fill=dipyrimidine))+
    geom_bar(stat='identity',width = 0.8,position = position_dodge(width=0.8)) +
    labs(x = 'Position',y = 'ATL-XR-seq / XR-seq dipyrimidine rate (log2)')+
    coord_cartesian(ylim=c(-1,2))+ 
    scale_fill_manual(values=c("TC" = "#E41A1C", "TT" = "#984EA3", "CT" = "#377EB8" , "CC" = "#A65628"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
          panel.border = element_blank(), 
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          
          strip.text = element_text(size=rel(1.5),face="bold"),
          strip.background = element_blank(),
          
          axis.text.x =element_text(size=rel(1.4),angle=90),
          axis.text.y =element_text(size=rel(1.4),angle=0),
          axis.title=element_text(size=rel(1.3),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.6), linetype = "solid")
    )+
    ggtitle(unlist(str_split(sample,"-",simplify=T))[,2])
  assign(paste("pic", nub, sep = ""),pic)
  #plotlist[[i]] <- p
}
#15:6
p <- ggpubr::ggarrange(pic1,pic2,
                       ncol = 2,nrow=1)
i=i+1
assign(paste("p", i, sep = ""),p)
ggsave(plot = p,filename = "compare_ditt.pdf",width=15,height = 6,device = "pdf")

##baserate=========
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/baserate/")
samples <- list.files()[grep("baserate.txt",list.files())]
samples <- samples[grep("^XRseq",samples)]

i=0
for (sample in samples){
  
  bb <- read.table(sprintf("%s",sample),header = T,row.names =1)
  sumbase <- sum(bb[1,])
  bb <- bb*100/sumbase 
  bb <- as.matrix(bb[,c(1,4,3,2)])
  bb <- as.data.frame(as.table(bb))
  bb$seq <- rep("XR seq",nrow(bb))
  
  aa <- read.table(sprintf("ATL-%s-rep1-markdup_26nt_baserate.txt",unlist(str_split(sample,"-rep1",simplify=T))[,1]),header = T,row.names =1)
  sumbase <- sum(aa[1,])
  aa <- aa*100/sumbase 
  aa <- as.matrix(aa[,c(1,4,3,2)])
  aa <- as.data.frame(as.table(aa))
  aa$seq <- rep("ATL-XR-seq rep1",nrow(aa))
  
  cc <- read.table(sprintf("ATL-%s-rep2-markdup_26nt_baserate.txt",unlist(str_split(sample,"-rep1",simplify=T))[,1]),header = T,row.names =1)
  sumbase <- sum(cc[1,])
  cc <- cc*100/sumbase 
  cc <- as.matrix(cc[,c(1,4,3,2)])
  cc <- as.data.frame(as.table(cc))
  cc$seq <- rep("ATL-XR-seq rep2",nrow(cc))
  
  baserate <- bind_rows(bb,aa,cc)
  colnames(baserate) <- c("pos","base","Freq","seq")
  
  p <- ggplot(baserate,mapping = aes(pos,Freq,fill=base))+
    facet_wrap(.~ seq)+
    geom_bar(stat='identity',position='stack') +
    labs(x = 'Position',y = 'base rate (%)') +
    coord_cartesian(ylim=c(0,100))+ 
    scale_fill_manual(values=c("A" = "#006837", "T" = "#FC4E2A", "C" = "#542788" , "G" = "#FED976"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.7),face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          legend.title = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.text = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3),face="bold"),
          strip.background = element_blank(),
          
          axis.text.x =element_text(size=rel(1.0),angle=0),
          axis.text.y =element_text(size=rel(1.5),angle=0),
          axis.title=element_text(size=rel(1.4),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
    )+
  ggtitle(unlist(str_split(sample,"-",simplify=T))[,c(2,3)])
  ggsave(filename = sprintf("%s_comb.pdf",unlist(str_split(sample,"_26",simplify=T))[,1]),width=15,height = 5)
  ##
  #print("ok")
  #i=i+1
  #assign(paste("p", i, sep = ""),p)
}


##compare-baserate=========
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/baserate/")
samples <- list.files()[grep("baserate.txt",list.files())]
samples <- samples[grep("^XRseq",samples)]

i=0
for (sample in samples){
  
  bb <- read.table(sprintf("ATL-%s",sample),header = T,row.names =1,sep="\t",stringsAsFactors = T) %>% .[c(3,2)]
  bb <- bb*100/sum(bb[1,])

  

  samples2 <- unlist(str_split(sample,"-rep1",simplify=T))[,1]
  aa <- read.table(sprintf("%s",sample),header = T,row.names =1,sep="\t",stringsAsFactors = T) %>% .[c(3,2)]
  aa <- aa*100/sum(aa[1,])
 
   baserate <- log2(bb/aa) %>% .[seq(18,24),]
  
   baserate <- t(baserate) %>% as.table(.) %>% as.data.frame()

  

  colnames(baserate) <- c("base","pos","Freq")
  
  p <- ggplot(baserate,mapping = aes(pos,Freq,fill=base))+
    geom_bar(stat='identity',width = 0.8,position = position_dodge(width=0.8)) +
    labs(x = 'Position',y = 'ATL-XR-seq / XR-seq base rate (log2)')+
    coord_cartesian(ylim=c(-1,2))+ 
    scale_fill_manual(values=c("A" = "#006837", "T" = "#FC4E2A", "C" = "#542788" , "G" = "#FED976"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          legend.title = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.text = element_text(size=rel(1.4)),
          strip.text = element_text(size=rel(1),face="bold"),
          strip.background = element_blank(),
          
          axis.text.x =element_text(size=rel(2),angle=0),
          axis.text.y =element_text(size=rel(1.4),angle=0),
          axis.title=element_text(size=rel(1.5),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
    )+
    ggtitle(unlist(str_split(sample,"-",simplify=T))[,c(2,3)])
  #ggsave(filename = sprintf("%s_comb.pdf",unlist(str_split(sample,"_26",simplify=T))[,1]),width=15,height = 5)
  ##
  #print("ok")
  i=i+1
  assign(paste("p", i, sep = ""),p)
}
#15:6
p <- ggpubr::ggarrange(p1,p2,
                  ncol = 2,nrow=1, font.label = list(color = 'black'))
ggsave(plot = p,filename = "compare.pdf",width=15,height = 5,device = "pdf")


##19~20nt TTor TC baserate=========
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/1920/")
samples <- list.files()[grep("baserate.txt",list.files())]
samples <- samples[grep("^XRseq",samples)]
#samples <- samples[grep("CT",samples)]
i=0
for (sample in samples){
  pos <- unlist(str_split(sample,"_"))[2]
  ditt <- unlist(str_split(sample,"_"))[3]
  damage <- unlist(str_split(unlist(str_split(sample,"_"))[1],"-"))[2]
  a <- "5\'nt"
  b <- unlist(str_split(ditt,""))
  d <- "3\'nt"
  
  bb <- read.table(sprintf("%s",sample),header = T,row.names =1)
  sumbase <- sum(bb[1,])
  bb <- bb*100/sumbase 
  bb <- as.matrix(bb[,c(1,4,3,2)])
  bb <- as.data.frame(as.table(bb))
  bb$seq <- rep("oXR-seq",nrow(bb))
  
  sample2 <- sprintf("ATL-%s",sample)
  aa <- read.table(sprintf("%s",sample2),header = T,row.names =1)
  sumbase <- sum(aa[1,])
  aa <- aa*100/sumbase 
  aa <- as.matrix(aa[,c(1,4,3,2)])
  aa <- as.data.frame(as.table(aa))
  aa$seq <- rep("ATL-XR-seq rep1",nrow(aa))
  
  sample3 <- sprintf("ATL-%s-rep2%s",unlist(str_split(sample,"-rep1",simplify=T))[,1],unlist(str_split(sample,"-rep1",simplify=T))[,2])
  cc <- read.table(sprintf("%s",sample3),header = T,row.names =1)
  sumbase <- sum(cc[1,])
  cc <- cc*100/sumbase 
  cc <- as.matrix(cc[,c(1,4,3,2)])
  cc <- as.data.frame(as.table(cc))
  cc$seq <- rep("ATL-XR-seq rep2",nrow(cc))
  
  dd <- read.table(sprintf("genome-markdup_%s_%s_26nt_baserate.txt",pos,ditt),header = T,row.names =1)
  sumbase <- sum(dd[1,])
  dd <- dd*100/sumbase 
  dd <- as.matrix(dd[,c(1,4,3,2)])
  dd <- as.data.frame(as.table(dd))
  dd$seq <- rep("genome",nrow(dd))
  
  baserate <- bind_rows(bb,aa,cc,dd)
  colnames(baserate) <- c("pos","base","Freq","seq")
  
  baserate$seq <- factor(baserate$seq,levels = c("ATL-XR-seq rep1","ATL-XR-seq rep2","oXR-seq","genome"),ordered=TRUE)
  
  p <- ggplot(baserate,mapping = aes(pos,Freq,fill=base))+
    facet_wrap(.~ seq,nrow=1)+
    geom_bar(stat='identity',position='stack') +
    labs(x = 'Position',y = 'base rate (%)') +
    coord_cartesian(ylim=c(0,100))+ 
    scale_x_discrete(breaks=c(1,2,3,4),labels = c(a,b,d))+##
    scale_fill_manual(values=c("A" = "#006837", "T" = "#FC4E2A", "C" = "#542788" , "G" = "#FED976"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          legend.title = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.text = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3),face="bold"),
          strip.background = element_blank(),
          
          axis.text.x =element_text(size=rel(1.7),angle=0),
          axis.text.y =element_text(size=rel(1.4),angle=0),
          axis.title=element_text(size=rel(1.5),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
    )+
    ggtitle(sprintf("%s %s-%s %s",damage,pos,as.numeric(pos)+1,ditt))
  #ggsave(filename = sprintf("%s_comb.pdf",unlist(str_split(sample,"_26",simplify=T))[,1]),width=10,height = 5)
  ##
  #print("ok")
  i=i+1
  assign(paste("p", i, sep = ""),p)
}
#20:20
p <- ggpubr::ggarrange(p1,p2,p11,p3,p4,p12,p5,p6,p13,p7,p8,p14,p9,p10,p15,
                  nrow=5,ncol = 3) 
ggsave(plot = p,filename = "TTANDTC_ALL2.pdf",width=33,height = 20,device = "pdf")
#20 and 21
p <- ggpubr::ggarrange(p7,p8,p14,
                       nrow=1,ncol = 3, font.label = list(color = 'black')) 

ggsave(plot = p,filename = "TTANDTC_20-21.pdf",width=25,height = 4,device = "pdf")







##chromstate1=====================
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/chromstate/")
name_exchange <-  function(x){
  unlist(str_split(x[1],"_")) %>% str_c(.,collapse=" ") 
}

samples <- list.files()[grep("_ChromstateData.txt",list.files())]
samples <- samples[grep("^XRseq",samples)]
#samples <- samples[grep("hotpeak",samples)]
#samples <- samples[grep("64PP",samples)] ##64PP;CPD
i=0
for (sample in samples){
  x <- sample
  n <- unlist(str_split(x,"-",simplify=T))
  damage <- n[length(n)-2]
  dup <- unlist(str_split(n[length(n)],"_",simplify=T))[1]
  rp <- n[length(n)-1]
  seq <- str_sub(unlist(str_split(x,damage,simplify=T)[1]),1,-2)
  seq2 <- str_c(c(unlist(str_split(seq,"XRseq")[1]),"XR-seq"),collapse = "")
  
  LoXR <- read.table(sprintf("ATL-XRseq-%s-%s-markdup_ChromstateData.txt",damage,rp),sep = "\t",header = T) %>% .[c(1,9,6,2,3,4,10,7,8,5),]
  LoXR$chromstate <- apply(LoXR,1,name_exchange)
  LoXR$chromstate <- str_c(seq(1,10),LoXR$chromstate,sep=".")
  LoXR$state <- rep("ATL-XR-seq rep1",nrow(LoXR))
  LoXR$chromstate = factor(LoXR$chromstate,levels = LoXR$chromstate)
  
  LoXR2 <- read.table(sprintf("ATL-XRseq-%s-rep2-markdup_ChromstateData.txt",damage),sep = "\t",header = T) %>% .[c(1,9,6,2,3,4,10,7,8,5),]
  LoXR2$chromstate <- apply(LoXR2,1,name_exchange)
  LoXR2$chromstate <- str_c(seq(1,10),LoXR2$chromstate,sep=".")
  LoXR2$state <- rep("ATL-XR-seq rep2",nrow(LoXR2))
  LoXR2$chromstate = factor(LoXR2$chromstate,levels = LoXR2$chromstate)
  
  XR <- read.table(sprintf("%s-%s-%s-markdup_ChromstateData.txt",seq,damage,rp),sep = "\t",header = T)%>% .[c(1,9,6,2,3,4,10,7,8,5),]
  XR$chromstate <- apply(XR,1,name_exchange)
  XR$state <- rep("oXR-seq",nrow(XR))
  XR$chromstate <- str_c(seq(1,10),XR$chromstate,sep=".")
  XR$chromstate = factor(XR$chromstate,levels = XR$chromstate)
  
  aa <- bind_rows(LoXR[,c(4,5,6)],LoXR2[,c(4,5,6)],XR[,c(4,5,6)])
  #aa <- bind_rows(LoXR[,c(4,5,6)],XR[,c(4,5,6)])
  
  p <- ggplot(data=aa,mapping=aes(x=chromstate,y=RPKM,fill=chromstate))+
    facet_wrap(.~ state)+
    geom_bar(stat="identity")+
    #expand_limits(y=c(0,20,40,60,80,100))+
    expand_limits(y=c(0,10,20,30))+
    ggtitle(sprintf("comparsion %s %s",damage,rp))+
    labs(x="chromatin state", y=("Relative value"))+
    scale_x_discrete(labels= as.character(seq(1,nrow(XR))))+
    # scale_fill_manual(values =c("#E31A1C","#FF7F00" , "#FED976", "#08306B","#2171B5",
    #                             "#006837", "#41AB5D","#B15928","#525252","#969696"))+
    scale_fill_manual(values =c("#A50F15","#EF3B2C" , "#54278F", "#FD8D3C","#FED976",
                                "#2171B5","#006837", "#41AB5D","#525252","#969696"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.2),face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.title = element_blank(),
          strip.text = element_text(size=rel(1),face="bold"),
          strip.background = element_blank(),
          strip.text.x = element_text(size = rel(1.3), color = "black", face = "bold"),
          panel.spacing.x = unit(0.3, "inches"),
          axis.text.x =element_text(size=rel(1.4),angle=0),
          axis.text.y =element_text(size=rel(1.2),angle=0),
          axis.title=element_text(size=rel(1.3),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))
  ggsave(filename = sprintf("%s_%s_all_comb.pdf",damage,rp),width=18,height = 5)
  #i=i+1
  #assign(paste("p", i, sep = ""),p)
}


##chromstate2==============
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/chromstate2/")
name_exchange <-  function(x){
  unlist(str_split(x[1],"_")) %>% str_c(.,collapse=" ") 
}

samples <- list.files()[grep("_ChromstateData.txt",list.files())]
samples <- samples[grep("^XRseq",samples)]
i=0
for (sample in samples){
  LoXR <- read.table(sprintf("ATL-%s",sample),sep = "\t",header = T) %>% .[c(1,6,9,2,3,10,7,4,5,8),]
  LoXR$chromstate <- apply(LoXR,1,name_exchange)
  LoXR$chromstate <- str_c(seq(1,10),LoXR$chromstate,sep=".")
  LoXR$state <- rep("ATL-XR-seq rep1",nrow(LoXR))
  LoXR$chromstate = factor(LoXR$chromstate,levels = LoXR$chromstate)
  
  LoXR2 <- read.table(sprintf("ATL-%s-rep2%s",unlist(str_split(sample,"-rep1",simplify=T))[,1],unlist(str_split(sample,"-rep1",simplify=T))[,2]),sep = "\t",header = T) %>% .[c(1,6,9,2,3,10,7,4,5,8),]
  LoXR2$chromstate <- apply(LoXR2,1,name_exchange)
  LoXR2$chromstate <- str_c(seq(1,10),LoXR2$chromstate,sep=".")
  LoXR2$state <- rep("ATL-XR-seq rep2",nrow(LoXR2))
  LoXR2$chromstate = factor(LoXR2$chromstate,levels = LoXR2$chromstate)
  
  XR <- read.table(sprintf("%s",sample),sep = "\t",header = T)%>% .[c(1,6,9,2,3,10,7,4,5,8),]
  XR$chromstate <- apply(XR,1,name_exchange)
  XR$state <- rep("XR-seq",nrow(XR))
  XR$chromstate <- str_c(seq(1,10),XR$chromstate,sep=".")
  XR$chromstate = factor(XR$chromstate,levels = XR$chromstate)
  
  aa <- bind_rows(LoXR[,c(4,5,6)],LoXR2[,c(4,5,6)],XR[,c(4,5,6)])
  
  p <- ggplot(data=aa,mapping=aes(x=chromstate,y=foldcg,fill=chromstate))+
    facet_wrap(.~ state)+
    geom_bar(stat="identity")+
    expand_limits(y=c(-1.0,-0.5,0,0.5,1.0,1.5))+
    labs(x="chromatin state", y=sprintf("log2 %s / input ",unlist(str_split(sample,"-",simplify=T)[,2] )))+
    scale_x_discrete(labels= as.character(seq(1,nrow(XR))))+
    scale_fill_manual(values =c("#E31A1C","#FF7F00" , "#FED976", "#08306B","#2171B5",
                                "#006837", "#41AB5D","#B15928","#525252","#969696"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.2),face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.title = element_blank(),
          strip.text = element_text(size=rel(1),face="bold"),
          strip.background = element_blank(),
          strip.text.x = element_text(size = rel(1.3), color = "black", face = "bold"),
          panel.spacing.x = unit(0.3, "inches"),
          axis.text.x =element_text(size=rel(1.4),angle=0),
          axis.text.y =element_text(size=rel(1.2),angle=0),
          axis.title=element_text(size=rel(1.3),face="bold"),
          axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))
  ggsave(filename = sprintf("%s_chromstate_comb.pdf",unlist(str_split(sample,"-markdup",simplify=T))[,1]),width=13,height = 5)
  #i=i+1
  #assign(paste("p", i, sep = ""),p)
}



##tss and tes -3kb and 10kb=================
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/tssANDtes/")
library(ggplot2)
library(reshape2)
library(readr)
library(ggpubr)
library(cowplot)
library(stringr)
library(dplyr)
library(RColorBrewer)

#RPKM
rpkm <- read.table("data_reads.count",sep="",stringsAsFactors = F)  ##
rpkm$name <- apply(rpkm,1,function(x)unlist(str_split(x[2],"[/]",simplify = F))[3])
rpkm$name <- apply(rpkm,1,function(x)unlist(str_split(x[3],".bed",simplify = F))[1]) 

minus <- 3000
plus <- 10000
tile <- 1

##============-----64PP--------==========##
  samples <- list.files()[grep("-markdup_NTS_TSS_cov",list.files())]
  samples <- samples[grep("64PP",samples)]
  i=0
  #plotlist <- list()
  for (sample in samples){
    x <- sample
    n <- unlist(str_split(x,"-",simplify=T))
    damage <- n[length(n)-2]
    strand <- unlist(str_split(n[length(n)],"_",simplify=T))[2]
    pos <- unlist(str_split(n[length(n)],"_",simplify=T))[3]
    rp <- n[length(n)-1]
    seq <- str_sub(unlist(str_split(x,damage,simplify=T)[1]),1,-2)
    seq2 <- str_c(c(unlist(str_split(seq,"XRseq")[1]),"XR-seq"),collapse = "")
    
    
    #tss
    tss <- data.frame(seq(-minus,plus)) 
    NTS <- read_table2(sprintf("%s",sample),
                       col_names = F,
                       col_types = cols(
                         x1="c",
                         x2="i",
                         x3="i"
                       )) %>% .[,-1]
    TS <- read_table2(sprintf("%s_TS_TSS_cov",unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1])),
                      col_names = F,
                      col_types = cols(
                        x1="c",
                        x2="i",
                        x3="i"
                      )) %>% .[,-1]
    #print(filter(rpkm,name==unlist(str_split(sample,"_NTS_cov",simplify=T)[1]))[,1]*tile*0.001)
    NTS$X3 <- NTS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    NTS [,1] <- NTS [,1]-minus-200
    NTS <- NTS[seq(200,200+minus+plus),]
    TS$X3 <- TS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    colnames(NTS)[2] <- "NTS"  #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    TS[,1] <- TS[,1]-minus-200
    TS <- TS[seq(200,200+minus+plus),]
    colnames(TS)[2] <- "TS" #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    
    tss <- cbind(tss,NTS[,2])
    tss <- cbind(tss,TS[,2])
    colnames(tss)[1] <- "X2"
    
    tss <- melt(tss,id="X2",measure=colnames(tss[,seq(2,ncol(tss))]))
    tss$state <- rep("TSS",nrow(tss))
    
    #tes
    
    tes <- data.frame(seq(-plus,minus)) 
    NTS <- read_table2(sprintf("%s-markdup_NTS_TES_cov",unlist(str_split(sample,"-markdup_NTS_TSS",simplify=T))[,1]),
                       col_names = F,
                       col_types = cols(
                         x1="c",
                         x2="i",
                         x3="i"
                       )) %>% .[,-1]
    TS <- read_table2(sprintf("%s-markdup_TS_TES_cov",unlist(str_split(sample,"-markdup_NTS_TSS",simplify=T))[,1]),
                      col_names = F,
                      col_types = cols(
                        x1="c",
                        x2="i",
                        x3="i"
                      )) %>% .[,-1]
    #print(filter(rpkm,name==unlist(str_split(sample,"_NTS_cov",simplify=T)[1]))[,1]*tile*0.001)
    NTS$X3 <- NTS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    NTS [,1] <- NTS [,1]-plus-200
    NTS <- NTS[seq(200,200+plus+minus),]
    TS$X3 <- TS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    colnames(NTS)[2] <- "NTS"  #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    TS[,1] <- TS[,1]-plus-200
    TS <- TS[seq(200,200+plus+minus),]
    colnames(TS)[2] <- "TS" #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    
    tes <- cbind(tes,NTS[,2])
    tes <- cbind(tes,TS[,2])
    colnames(tes)[1] <- "X2"
    
    tes <- melt(tes,id="X2",measure=colnames(tes[,seq(2,ncol(tes))]))
    tes$state <- rep("TES",nrow(tes))
    
    genome <- bind_rows(tss,tes)
    genome$state <- factor(genome$state,levels = c("TSS","TES"))
    colnames(genome) <- c( "position","strand","value","state")
    p <- ggplot(genome, aes(x=position, y=value,group=strand,colour=strand)) + 
      geom_line(size=0.6) +
      facet_wrap(.~ state,scales = "free_x")+
      scale_x_continuous(breaks=c(-10000,-3000,0,3000,10000),labels = c("-10k","-3k","TSS","3k","10k")) +##
      expand_limits(y=c(0.025,0.05,0.075,0.1))+
      scale_color_manual(values =c("#4393C3","#EF6548"))+
      labs(title=sprintf("%s %s",seq2,rp),x="64PP" ,y="RPKM") +
      theme(plot.title = element_text(hjust = 0.5,size=rel(1.2),face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            strip.text = element_text(size=rel(1),face="bold"),
            strip.background = element_blank(),
            axis.text.x =element_text(size=rel(1.4),angle=0),
            axis.text.y =element_text(size=rel(1.2),angle=0),
            axis.title=element_text(size=rel(1.3),face="bold"),
            axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))
    
    i=i+1
    assign(paste("p", i, sep = ""),p)
  }
  
##============-----CPD--------==========##
  samples <- list.files()[grep("-markdup_NTS_TSS_cov",list.files())]
  samples <- samples[grep("CPD",samples)]
  #plotlist <- list()
  for (sample in samples){
    x <- sample
    n <- unlist(str_split(x,"-",simplify=T))
    damage <- n[length(n)-2]
    strand <- unlist(str_split(n[length(n)],"_",simplify=T))[2]
    pos <- unlist(str_split(n[length(n)],"_",simplify=T))[3]
    rp <- n[length(n)-1]
    seq <- str_sub(unlist(str_split(x,damage,simplify=T)[1]),1,-2)
    seq2 <- str_c(c(unlist(str_split(seq,"XRseq")[1]),"XR-seq"),collapse = "")
    
    #tss
    tss <- data.frame(seq(-minus,plus)) 
    NTS <- read_table2(sprintf("%s",sample),
                       col_names = F,
                       col_types = cols(
                         x1="c",
                         x2="i",
                         x3="i"
                       )) %>% .[,-1]
    TS <- read_table2(sprintf("%s_TS_TSS_cov",unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1])),
                      col_names = F,
                      col_types = cols(
                        x1="c",
                        x2="i",
                        x3="i"
                      )) %>% .[,-1]
    #print(filter(rpkm,name==unlist(str_split(sample,"_NTS_cov",simplify=T)[1]))[,1]*tile*0.001)
    NTS$X3 <- NTS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    NTS [,1] <- NTS [,1]-minus-200
    NTS <- NTS[seq(200,200+minus+plus),]
    TS$X3 <- TS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    colnames(NTS)[2] <- "NTS"  #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    TS[,1] <- TS[,1]-minus-200
    TS <- TS[seq(200,200+minus+plus),]
    colnames(TS)[2] <- "TS" #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    
    tss <- cbind(tss,NTS[,2])
    tss <- cbind(tss,TS[,2])
    colnames(tss)[1] <- "X2"
    
    tss <- melt(tss,id="X2",measure=colnames(tss[,seq(2,ncol(tss))]))
    tss$state <- rep("TSS",nrow(tss))
    
    #tes
    
    tes <- data.frame(seq(-plus,minus)) 
    NTS <- read_table2(sprintf("%s-markdup_NTS_TES_cov",unlist(str_split(sample,"-markdup_NTS_TSS",simplify=T))[,1]),
                       col_names = F,
                       col_types = cols(
                         x1="c",
                         x2="i",
                         x3="i"
                       )) %>% .[,-1]
    TS <- read_table2(sprintf("%s-markdup_TS_TES_cov",unlist(str_split(sample,"-markdup_NTS_TSS",simplify=T))[,1]),
                      col_names = F,
                      col_types = cols(
                        x1="c",
                        x2="i",
                        x3="i"
                      )) %>% .[,-1]
    #print(filter(rpkm,name==unlist(str_split(sample,"_NTS_cov",simplify=T)[1]))[,1]*tile*0.001)
    NTS$X3 <- NTS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    NTS [,1] <- NTS [,1]-plus-200
    NTS <- NTS[seq(200,200+plus+minus),]
    TS$X3 <- TS$X3*10^6/(filter(rpkm,name==unlist(str_split(sample,"_NTS_TSS_cov",simplify=T)[1]))[,1]*tile*1000)
    colnames(NTS)[2] <- "NTS"  #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    TS[,1] <- TS[,1]-plus-200
    TS <- TS[seq(200,200+plus+minus),]
    colnames(TS)[2] <- "TS" #sprintf("%s",unlist(str_split(sample,"[_]",simplify=T))[length(str_split(sample,"[_]",simplify=T))-1])
    
    tes <- cbind(tes,NTS[,2])
    tes <- cbind(tes,TS[,2])
    colnames(tes)[1] <- "X2"
    
    tes <- melt(tes,id="X2",measure=colnames(tes[,seq(2,ncol(tes))]))
    tes$state <- rep("TES",nrow(tes))
    
    genome <- bind_rows(tss,tes)
    genome$state <- factor(genome$state,levels = c("TSS","TES"))
    colnames(genome) <- c( "position","strand","value","state")
    p <- ggplot(genome, aes(x=position, y=value,group=strand,colour=strand)) + 
      geom_line(size=0.6) +
      ggtitle(sprintf("%s %s",seq2,rp))+
      facet_wrap(.~ state,scales = "free_x")+
      scale_x_continuous(breaks=c(-10000,-3000,0,3000,10000),labels = c("-10k","-3k","TSS","3k","10k")) +##
      expand_limits(y=c(0.05,0.1,0.15))+
      scale_color_manual(values =c("#4393C3","#EF6548"))+
      labs(title=sprintf("%s %s",seq2,rp),x="CPD", y="RPKM") +
      theme(plot.title = element_text(hjust = 0.5,size=rel(1.2),face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            strip.text = element_text(size=rel(1),face="bold"),
            strip.background = element_blank(),
            axis.text.x =element_text(size=rel(1.4),angle=0),
            axis.text.y =element_text(size=rel(1.2),angle=0),
            axis.title=element_text(size=rel(1.3),face="bold"),
            axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))
    
    i=i+1
    assign(paste("p", i, sep = ""),p)
  }
  

  
  
  
  
  #3x2  21:6
  p <- ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,
                    ncol = 3,nrow=2,labels = LETTERS[1:i], font.label = list(color = 'black'))
  ggsave(plot = p,filename = "ALLGENE.pdf",width=21,height = 6,device = "pdf")
  
  

  
##TSandNTS correlation=====
  setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/cor_bin/eachgene/")
  samples <- list.files()[grep("_cor.bed",list.files())]
  samples <- samples[grep("NTS",samples)]
  rpkm <- read.table("data_reads.count",sep="",stringsAsFactors = F)  ##
  rpkm$name <- apply(rpkm,1,function(x)unlist(str_split(x[2],"[/]",simplify = F))[3])
  rpkm$name <- apply(rpkm,1,function(x)unlist(str_split(x[3],"-markdup",simplify = F))[1]) 
  
  i=0
  for (sample in samples){
    x <- sample
    n <- unlist(str_split(x,"-",simplify=T))
    damage <- n[length(n)-2]
    strand <- unlist(str_split(n[length(n)],"_",simplify=T))[2]
    rp <- n[length(n)-1]
    seq <- str_sub(unlist(str_split(x,damage,simplify=T)[1]),1,-2)
    seq2 <- str_c(c(unlist(str_split(seq,"XRseq")[1]),"XR-seq"),collapse = "")
    
    aa <- read_table2(sprintf("%s",sample),
                      col_names = F,
                      col_types=cols(
                        x1="c",
                        x2="i",
                        x3="i",
                        x4="c",
                        x5="c",
                        x6="c",
                        x7="i"
                      ))
    aa$X7 <- aa$X7*10^6/(filter(rpkm,name==unlist(str_split(sample,"-markdup",simplify=T)[1]))[,1])
    aa <- data.frame(genes=aa$X4,num1=aa$X7)
    
    sample2 <- sprintf("%sTS%s",unlist(str_split(sample,"NTS",simplify=T))[,1],unlist(str_split(sample,"NTS",simplify=T))[,2])
    bb <- read_table2(sprintf("%s",sample2),
                      col_names = F,
                      col_types=cols(
                        x1="c",
                        x2="i",
                        x3="i",
                        x4="c",
                        x5="c",
                        x6="c",
                        x7="i"
                      ))
    bb$X7 <- bb$X7*10^6/(filter(rpkm,name==unlist(str_split(sample,"-markdup",simplify=T)[1]))[,1])
    bb <- data.frame(genes=bb$X4,num2=bb$X7)
    
    cor <- cbind(aa,bb) %>% .[,c(1,2,4)]
    
    cor <- within(cor,{fold <- NA
    fold[cor$num1==0 & cor$num2==0] <- 0
    fold[cor$num1==0 & cor$num2!=0] <- log2(1.5)
    fold[cor$num1!=0 & cor$num2==0] <- log2(1/1.5)
    fold[cor$num1!=0 & cor$num2!=0] <- NA})
    cor$fold <- apply(cor,1,function(x){
      if(is.na(x[4]))
      {return(log2(as.numeric(x[3])/as.numeric(x[2])))}else{return(x[4])}
    })
    
    cor <- within(cor,{color <- NA
    color[as.numeric(cor$fold)>=log2(1/1.5) & as.numeric(cor$fold) <=log2(1.5)] <- "black"
    color[as.numeric(cor$fold) >=log2(1.5)] <- "red"
    color[as.numeric(cor$fold) <=log2(1/1.5)] <- "blue"})
    
    if (seq2=="XR-seq"){
      title <- sprintf("%s %s",seq2,damage)
    }else{title <- sprintf("%s %s %s",seq2,rp,damage)}
    
    p <- ggscatter(cor, x = "num1", y = "num2",
                   color = cor$color, size = 1.2, # 点的颜色与大???
                   add = "reg.line",  # 添加回归???  
                   add.params = list(color = "red", fill = "lightgray",linetype=5), # 回归线的颜色设置为红色，区间颜色设置为灰???
                   conf.int = F, # 添加回归线的置信区间
                   cor.coef = TRUE, # 添加相关系数
                   cor.coeff.args = list(method = "spearman", label.x = 50,label.y=340,label.sep = "\n",r.digits=4))+#选择Pearson相关+
      labs(x = "NTS",
           y = "TS")+
      ggtitle(sprintf("%s",title))+
      coord_cartesian(ylim=c(0,400),xlim=c(0,400))+ 
      theme(plot.title = element_text(hjust = 0.5,size=rel(1.2),face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            legend.title = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            
            strip.text = element_text(size=rel(1),face="bold"),
            strip.background = element_blank(),
            
            axis.text.x =element_text(size=rel(1.2),angle=0),
            axis.text.y =element_text(size=rel(1.2),angle=0),
            axis.title=element_text(size=rel(1.3),face="bold"),
            axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
      )
    i=i+1
    assign(paste("p", i, sep = ""),p)
  }
  
  p <- ggpubr::ggarrange(p1,p2,p5,p3,p4,p6,
                         ncol = 3,nrow=2)
ggsave(plot = p,filename = "eachgeneTSandNTS.png",width=12,height = 8,device = "png")



##qPCR PICTURES=============================
library(ggpmisc)
ntNN <- data.frame(fmol=c(0,1,4,8,32,64),q=c(0,0.004507604,0.028833331,0.061309175,0.325394924,0.761666717
))


p <- ggscatter(ntNN, x = "fmol",y = "q",
                size = 2, 
               add = "reg.line",  # 添加回归???  
               add.params = list(color = "red", fill = "lightgray",linetype=5), 
               conf.int = F, 
               cor.coef = TRUE, 
               cor.coeff.args = list(method = "pearson", label.x = 8,label.y=0.75,label.sep = "\n",r.digits=5,size=9))+
  labs(x = "5P-26ntNN (fmol)",
       y = "")+
  
  scale_x_continuous(breaks=ntNN$fmol,labels = as.character(ntNN$fmol))+
  ggtitle("ATL-XR-qPCR")+
  #geom_text(x=40,y=0.75,label = paste0("y=",format(lm(q ~ fmol, ntNN)$coef[1],digits = 2),"+",format(lm(q ~ fmol, ntNN)$coef[2],digits = 2),"x"))+
  coord_cartesian(ylim=c(0,1),xlim=c(0,70))+ 
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.2),face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        
        strip.text = element_text(size=rel(1),face="bold"),
        strip.background = element_blank(),
        
        axis.text.x =element_text(size=rel(0.9),angle=0),
        axis.text.y =element_text(size=rel(1.2),angle=0),
        axis.title=element_text(size=rel(1.3),face="bold"),
        axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
  )


##polyA_lengthdistribution=====================
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/polyA_length/")
samples <- list.files()[grep("Alength_distribution3.txt",list.files())]
samples <- samples[grep("rep1",samples)]



for (sample in samples){
damage <- unlist(str_split(unlist(str_split(sample,"_"))[1],"-"))[3]
rep <- unlist(str_split(unlist(str_split(sample,"_"))[1],"-"))[4]
aa <- read.table(sprintf("%s",sample),header = T,sep = "\t",stringsAsFactors = F)
aa$count <- aa$count*100/sum(aa$count)
aa$rep.n <- "rep1"

bb <- read.table(sprintf("ATL-XRseq-%s-rep2_Alength_distribution3.txt",damage),header = T,sep = "\t",stringsAsFactors = F)
bb$count <- bb$count*100/sum(bb$count)
bb$rep.n <- "rep2"

pic <- rbind(aa,bb)

p <- ggplot(pic,mapping = aes(A.length,count))+
  facet_wrap(.~rep.n)+
  geom_bar(stat='identity',position='stack',width = 1) +
  labs(x = sprintf('%s A-tail length',damage),y = 'A-tail rate (%)') +
  coord_cartesian(ylim=c(0,10))+ 
  scale_x_continuous(breaks=c(10,30,50,80,100),labels=c(10,30,50,80,100))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        
        strip.text = element_text(size=rel(1.5),face="bold"),
        strip.background = element_blank(),
        
        axis.text.x =element_text(size=rel(1.6),angle=0),
        axis.text.y =element_text(size=rel(1.6),angle=0),
        axis.title=element_text(size=rel(1.3),face="bold"),
        axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
  )+
  ggtitle("ATL-XR-seq A-tail length distribution")


ggsave(filename = sprintf("ATL-XR-seq_%s_ALength3.pdf",damage),width=10,height = 5)

}


##ATL-XR-qPCR=============
library(readr)
library(tidyverse)
setwd("E:/LAB/IBS/experiment/LoXR-seq/end_data/figure-6/")

aa <- read_table2("qPCR-3.txt",
                  col_names = T,
                  col_types=cols(
                    x1="c",
                    x2="c",
                    x3="i",
                    x4="i",
                    x5="i"
                  )) 
aa$state <- factor(aa$state,levels=c("UV+DMSO","UV+1000nMTPL"))
#aa$percent <- aa$count/c(rep(0.087412686,3),rep(0.061948391,3),rep(0.029365369,3),rep(0.032291896,3))
aa$rp <- str_split_fixed(aa$group,"-",2)[,2]
aa$damage <- str_split_fixed(aa$group,"-",2)[,1]

aa$sd <- as.numeric(apply(aa[,3:5], 1, sd,na.rm=T))
aa$ebtop <- aa$count+aa$sd
aa$ebbottom <- aa$count-aa$sd

for (rps in c("rep1","rep2")){
  temp <- filter(aa,rp==rps)
  temp$relative <- temp$count/temp$count[1]
p <- ggplot(data=temp,mapping=aes(x=damage,y=relative,fill=state))+
  # geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=0.3, 
  #               color="black", ## 误差线颜色
  #               position=position_dodge(0.6))+
  geom_bar(stat="identity",width = 0.5,position = "dodge")+
  expand_limits(y=c(0,0.25,0.5,0.75,1))+
  ggtitle(sprintf("ATL-XR-qPCR %s",rps))+
  labs(x="", y="Relative repair signal")+
  #scale_x_discrete(labels=c("UV\n+DMSO","UV\n+1000 nM TPL"))+
# scale_fill_manual(values =c("#6BAED6", "#08306B"),
#                     labels=c("UV\n+0 nM TPL","UV\n+1000 nM TPL"))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.4),face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.1),angle=0),
        strip.text = element_text(size=rel(1),face="bold"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = rel(1.3), color = "black", face = "bold"),
        panel.spacing.x = unit(0.3, "inches"),
        axis.text.x =element_text(size=rel(1.5),angle=0,face = "bold"),
        axis.text.y =element_text(size=rel(1.9),angle=0),
        axis.title=element_text(size=rel(1.3),face="bold"),
        axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))

ggsave(filename = sprintf("%s_ATL-XR-qpcr.pdf",rps),width=6,height = 5)


}












# create viewports

library(grid)
#24:16
grid.newpage()
pushViewport(viewport(layout = grid.layout(16, 24)))
print(p1, vp=viewport(layout.pos.row=1:6, layout.pos.col=1:10))
print(p2, vp=viewport(layout.pos.row=1:3, layout.pos.col=13:24))
print(p3, vp=viewport(layout.pos.row=4:6, layout.pos.col=13:24))
print(p4, vp=viewport(layout.pos.row=7:10, layout.pos.col=1:14))
print(p5, vp=viewport(layout.pos.row=12:16, layout.pos.col=1:14))
print(p6, vp=viewport(layout.pos.row=10:13, layout.pos.col=15:24))





# direct the charts into the specified viewport
print(p.hist.len, vp=vp.len)
print(p.hist.wid, vp=vp.wid)
print(p.scatter, vp=vp.scatter)
##barcode pictures of RPM=======================
library(readr)
library(stringr)
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(plyr)

setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/barcode/forth8over1/")

chrom_state <- data.frame(abbr=c("Tss","TssF",
                                 "PromF",
                                 "PromP",
                                 "Enh","EnhF",
                                 "EnhWF", "EnhW", "DnaseU", "DnaseD", "FaireW",
                                 "CtcfO", "Ctcf",
                                 "Gen5'", "Elon", "ElonW", "Gen3'", "Pol2", "H4K20",
                                 "Low",
                                 "ReprD", "Repr", "ReprW",
                                 "Quies", "Art"),
                          details=c(rep("Active_Promoter",2),
                                    "Promoter_Flanking",
                                    "Inactive_Promoter",
                                    rep("Candidate_Strong_enhancer",2),
                                    rep("Candidate_Weak_enhancer/DNase",5),
                                    rep("Distal_CTCF/Candidate_Insulator",2),
                                    rep("Transcription_associated",6),
                                    "Low_activity_proximal_to_active_states",
                                    rep("Polycomb_repressed",3),
                                    rep("Heterochromatin/Repetitive/Copy_Number_Variation",2)))

name_exchange <-  function(x){
  dd <- c()
  for (i in c(x)){
    cc <- unlist(str_split(i,"_")) %>% str_c(.,collapse=" ") 
    dd <- c(dd,cc)
  }
  return(dd)
}
samples <- list.files()[grep(".txt",list.files())]
samples <- samples[grep("barcode",samples)]
samples <- samples[grep("hotpeakcount",samples)]
samples <- samples[grep("rep1",samples)] 
#samples <- samples[grep("CPD",samples)]##CPD;64PP
i=0
for (sample in samples){
  
  x <- sample
  n <- unlist(str_split(x,"-",simplify=T))
  damage <- n[length(n)-2]
  rp <- n[length(n)-1]
  PEAK <- unlist(str_split(n[length(n)],"_",simplify=T)[4])
  PEAK <- unlist(str_split(PEAK,".txt",simplify=T)[1])
  seq <- str_sub(unlist(str_split(x,damage,simplify=T)[1]),1,-2)
  seq2 <- str_c(c(unlist(str_split(seq,"XRseq")[1]),"XR-seq"),collapse = "")
  
  
  aa <- read_table2(sprintf("ATL-XRseq-%s-%s-barcode_bin_count_%s.txt",damage,rp,PEAK),
                    col_names = F,
                    col_types=cols(
                      x1="c",
                      x2="i"
                    )) 
  #aa <- aa %>% group_by(X8) %>% summarise(X10=n())
  #aa$pos <- apply(aa,1,function(x){temp <- str_c(x[1],x[2],x[3],x[6],sep = "-") return(temp)})
  
  aa$dup <- "A barcode"
  
  bb <- read_table2(sprintf("ATL-XRseq-%s-%s-markdup_bin_count_%s.txt",damage,rp,PEAK),
                    col_names = F,
                    col_types=cols(
                      x1="c",
                      x2="i"
                    )) 
  #bb$pos <- apply(bb,1,function(x){temp <- str_c(x[1],x[2],x[3],x[6],sep = "-")return(temp)})
  
  bb$dup <- "markdup"
  #bb <- bb %>% group_by(X8) %>% summarise(X10=n())
  temp <- bind_rows(aa,bb)
  
  temp$details <- apply(temp,1,function(x){return(filter(chrom_state,abbr==x[1])$details)})
  temp$details <- name_exchange(as.character(temp$details))
  
 
  
  for (dupways in c("A barcode","markdup")){
    circle.pic <- filter(temp,dup==dupways)
    
    circle.pic <- ddply(circle.pic,.(details),nrow)
    circle.pic[is.na(circle.pic)] <- 0
    circle.pic <- data.frame(details=c("others","Heterochromatin/Repetitive/Copy Number Variation")
                             ,V1=c(sum(filter(circle.pic,details!="Heterochromatin/Repetitive/Copy Number Variation")$V1),filter(circle.pic,details=="Heterochromatin/Repetitive/Copy Number Variation")$V1))
    total <- sum(circle.pic$V1)
    circle.pic$V1 <- circle.pic$V1*100/sum(circle.pic$V1)
    circle.pic$num <- seq(1,nrow(circle.pic))
    mylabel <- paste(circle.pic$details,"(",round(circle.pic$V1,2),"%)   ",sep="")
    p <- ggplot(circle.pic,aes(x="",y=V1,fill=details))+geom_bar(stat="identity")+coord_polar(theta = "y")+
      ggtitle(sprintf("%s %s %s",damage,rp,dupways))+
      labs(x="",y=sprintf("n=(%s)",total))+
      scale_fill_manual(name = "chromatin state",
                        values = c("#08306B","#969696"),
                        breaks = c("others","Heterochromatin/Repetitive/Copy Number Variation"),
                        labels=mylabel) +
      theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            #legend.title = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            legend.position="right",
            legend.text = element_text(size=rel(1)),
            strip.text = element_text(size=rel(1),face="bold"),
            strip.background = element_blank(),
            
            axis.text.x =element_blank(),
            axis.title=element_text(size=rel(1.3),face="bold"),
            axis.line = element_blank())
    
    i=i+1
    assign(paste("p", i, sep = ""),p)
  }
  for (dupways in c("A barcode","markdup")){
    circle.pic <- filter(temp,dup==dupways)
    circle.pic <- filter(circle.pic,details!="Heterochromatin/Repetitive/Copy Number Variation")
    chrom10 <- data.frame(details=factor(unique(chrom_state$details),levels = unique(chrom_state$details))
                          ,count=rep(0,10)) %>% .[-10,]
    chrom10$details <- name_exchange(as.character(chrom10$details))
    circle.pic <- ddply(circle.pic,.(details),nrow)
    circle.pic <- left_join(chrom10,circle.pic,by="details")
    circle.pic[is.na(circle.pic)] <- 0
    total <- sum(circle.pic$V1)
    circle.pic$V1 <- circle.pic$V1*100/sum(circle.pic$V1)
    circle.pic$details <- factor(circle.pic$details,levels=name_exchange(unique(as.character(chrom_state$details))))
    circle.pic$num <- seq(1,nrow(circle.pic))
    mylabel <- paste(circle.pic$details,"(",round(circle.pic$V1,2),"%)   ",sep="")
    # mylabel <- paste(circle.pic$num,circle.pic$details,sep=".")
    # mylabel2 <- paste(circle.pic$num, ". (",round(circle.pic$V1,2),"%)", sep = '')
    
    p <- ggplot(circle.pic,aes(x="",y=V1,fill=details))+geom_bar(stat="identity")+coord_polar(theta = "y")+
      ggtitle(sprintf("%s %s %s others",damage,rp,dupways))+
      labs(x="",y=sprintf("n=(%s)",total))+
      scale_fill_manual(name = "chromatin state",
                        values = c("#A50F15","#EF3B2C" , "#54278F", "#FD8D3C","#FED976",
                                   "#2171B5","#006837", "#41AB5D","#525252"),
                        breaks = c("Active Promoter","Promoter Flanking" , "Inactive Promoter", "Candidate Strong enhancer","Candidate Weak enhancer/DNase",
                                   "Distal CTCF/Candidate Insulator","Transcription associated", "Low activity proximal to active states",
                                   "Polycomb repressed"),
                        labels=mylabel) +
      theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            #legend.title = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            legend.position="right",
            legend.text = element_text(size=rel(1)),
            strip.text = element_blank(),
            strip.background = element_blank(),
            
            axis.text.x =element_blank(),
            axis.title=element_text(size=rel(1.3),face="bold"),
            axis.line = element_blank())
    
    i=i+1
    assign(paste("p", i, sep = ""),p)
    
  }
}

p <- ggpubr::ggarrange(p1,p3,p2,p4,p5,p7,p6,p8,
                       ncol = 2,nrow=4)
#ggsave(plot = p,filename = sprintf("dian_comp_barandmarkdup.png",damage),width=16,height = 12,device = "png")
#ggsave(plot = p,filename = sprintf("hotpeak_all",damage),width=8,height = 12,device = "pdf")

ggsave(plot = p,filename = "hotpeak_all.pdf",width=16,height = 12,device = "pdf")


##genome_chromatin state==========
name_exchange <-  function(x){
  dd <- c()
  for (i in c(x)){
    cc <- unlist(str_split(i,"_")) %>% str_c(.,collapse=" ") 
    dd <- c(dd,cc)
  }
  return(dd)
}
setwd("E:/LAB/IBS/experiment/infomatics/XRseq_try/OVER/barcode/third/")
aa <- read.table("hg19chromstate_DistributionFreq.txt",sep = "\t",header = T,stringsAsFactors = F)
aa$details <- name_exchange(aa$details)
aa$details <- factor(aa$details,levels=name_exchange(as.character(unique(chrom_state$details))))



circle.pic <- data.frame(details=c("others","Heterochromatin/Repetitive/Copy Number Variation")
                         ,V1=c(sum(filter(aa,details!="Heterochromatin/Repetitive/Copy Number Variation")$freqsum),filter(aa,details=="Heterochromatin/Repetitive/Copy Number Variation")$freqsum))

i=0

mylabel <- paste(circle.pic$details,"(",round(circle.pic$V1,2),"%)   ",sep="")

p <- ggplot(circle.pic,aes(x="",y=V1,fill=details))+geom_bar(stat="identity")+coord_polar(theta = "y")+
  ggtitle("The whole genome")+
  labs(x="",y="")+
  scale_fill_manual(name = "chromatin state",
                    values = c("#08306B","#969696"),
                    breaks = c("others","Heterochromatin/Repetitive/Copy Number Variation"),
                    labels=mylabel) +
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        #legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right",
        legend.text = element_text(size=rel(1)),
        strip.text = element_text(size=rel(1),face="bold"),
        strip.background = element_blank(),
        
        axis.text.x =element_blank(),
        axis.title=element_text(size=rel(1.3),face="bold"),
        axis.line = element_blank())

i=i+1
assign(paste("p", i, sep = ""),p)

circle.pic <- filter(aa,details!="Heterochromatin/Repetitive/Copy Number Variation")
chrom10 <- data.frame(details=factor(unique(chrom_state$details),levels = unique(chrom_state$details))
                      ,count=rep(0,10)) %>% .[-10,]
chrom10$details <- name_exchange(as.character(chrom10$details))
chrom10 <- left_join(chrom10,circle.pic,by="details")
chrom10$V1 <- chrom10$freqsum*100/sum(chrom10$freqsum)

chrom10$details <- factor(chrom10$details,levels=name_exchange(unique(as.character(chrom_state$details))))
mylabel <- paste(chrom10$details,"(",round(chrom10$V1,2),"%)   ",sep="")
# mylabel <- paste(circle.pic$num,circle.pic$details,sep=".")
# mylabel2 <- paste(circle.pic$num, ". (",round(circle.pic$V1,2),"%)", sep = '')

p <- ggplot(chrom10,aes(x="",y=V1,fill=details))+geom_bar(stat="identity")+coord_polar(theta = "y")+
  ggtitle("The whole genome")+
  labs(x="",y="")+
  scale_fill_manual(name = "chromatin state",
                    values = c("#A50F15","#EF3B2C" , "#54278F", "#FD8D3C","#FED976",
                               "#2171B5","#006837", "#41AB5D","#525252"),
                    breaks = c("Active Promoter","Promoter Flanking" , "Inactive Promoter", "Candidate Strong enhancer","Candidate Weak enhancer/DNase",
                               "Distal CTCF/Candidate Insulator","Transcription associated", "Low activity proximal to active states",
                               "Polycomb repressed"),
                    labels=mylabel) +
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        #legend.title = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right",
        legend.text = element_text(size=rel(1)),
        strip.text = element_blank(),
        strip.background = element_blank(),
        
        axis.text.x =element_blank(),
        axis.title=element_text(size=rel(1.3),face="bold"),
        axis.line = element_blank())

i=i+1
assign(paste("p", i, sep = ""),p)
p <- ggpubr::ggarrange(p1,p2,
                       ncol = 2,nrow=1)
ggsave(plot = p,filename = "genome.pdf",width=16,height = 3,device = "pdf")

##A-barcode and markdup barpictures==========
setwd("E:/LAB/IBS/experiment/LoXR-seq/end_data/figure-5/")
aa <- data.frame(dedup=c("64PP rep1","64PP rep2","CPD rep1","CPD rep2"),M=c(44.38,39.09,35.93,22.53)
                 ,A=c(55.53,51.97,46.56,28.17))
aa$A_barcode <- aa$A/aa$M
aa$markdup <- aa$M/aa$M
temp <- melt(aa[,c(1,4,5)])

p <- ggplot(data=temp,mapping=aes(x=dedup,y=value,fill=variable))+
  geom_bar(stat="identity",width = 0.5,position = "dodge")+
  expand_limits(y=c(0,0.25,0.5,0.75,1,1.5))+
  ggtitle("The number of reads after dedup in ATL-XR-seq")+
  labs(x="", y="Relative value")+
  #scale_x_discrete(labels=c("A barcode","UV\n+1000 nM TPL"))+
  scale_fill_manual(values =c("#35978F", "#BF812D"),
                       labels=c("A barcode","Markdup"))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.4),face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.1),angle=0),
        strip.text = element_text(size=rel(1),face="bold"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = rel(1.3), color = "black", face = "bold"),
        panel.spacing.x = unit(0.3, "inches"),
        axis.text.x =element_text(size=rel(1.7),angle=0,face = "bold"),
        axis.text.y =element_text(size=rel(1.7),angle=0),
        axis.title=element_text(size=rel(1.6),face="bold"),
        axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid"))

ggsave(filename = "readsnumberAandM.pdf",width=7,height = 5)

