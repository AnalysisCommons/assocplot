# Script generalized from 'make.plot_P.value_forManning.R' script written by Jaeyoung Hong

## Load libraries and source files
setwd("Commons Apps/assocplot/")
library(data.table)
source("modified.assocplot_P.value.R")

## Setup paths
#data.path<-""
# index.path<-"/restricted/project/aagile/Jae_MANTRA/AAGILE_T2D/regional_plot/"
# ld.path<-"/restricted/projectnb/aagile/Jae_MANTRA/AAGILE_T2D/regional_plot/LD/"
# plot.path<-"/restricted/projectnb/aagile/Jae_MANTRA/AAGILE_T2D/regional_plot/plot/"

index.path
ld.path
plot.path <- getwd()

## Setup parameters
snp.name <- "MarkerName"
gene.name <- ""
chr.name <- ""
p.name <- "P.value"
freq.name <- "Freq1"
df.name <- "HetDF"

wh.snp <- "16:56997233"
pop <- "CEU"

## Load data -- current datatype is METAL file
## read-in columns - MarkerName/snp, chr, pos, P.value/logBF ##
# fg.aagile.metal0<-fread(paste0("cut -f1,4,14,16-18,10 /restricted/project/aagile/meta/fg/fg.aagile.metal_chr.pos.tbl"),header=T,data.table=F)
# fib.aagile.metal0<-fread(paste0("cut -f1,4,14,16-18,10 /restricted/project/aagile/meta/fib/fib.aagile.metal_chr.pos.tbl"),header=T,data.table=F)
# fg.magic.metal0<-fread(paste0("cut -f1,4,14,16-18,10 /restricted/project/aagile/meta/fg/fg.magic.metal_chr.pos.tbl"),header=T,data.table=F)
# fib.magic.metal0<-fread(paste0("cut -f1,4,14,16-18,10 /restricted/project/aagile/meta/fib/fib.magic.metal_chr.pos.tbl"),header=T,data.table=F)


## read-in datasets ##
# read-in index snps - snp, chr, P.value.index, pos, Locus  #
print("Reading in data")
# aagile.index<-read.csv(paste(index.path,"index_top_snplist_add.indexPval.csv",sep=""),as.is=T)
# data0<-aagile.index


## inclusion criteria ##
# fg.aagile.metal<-subset(fg.aagile.metal0,Freq1>0.01 & Freq1<0.99 & HetDf>=4 | MarkerName%in%data0$indexsnp)
# fib.aagile.metal<-subset(fib.aagile.metal0,Freq1>0.01 & Freq1<0.99 & HetDf>=4 | MarkerName%in%data0$indexsnp)
# fg.magic.metal<-subset(fg.magic.metal0,(Freq1>0.01 & Freq1<0.99 & HetDf>=4) | MarkerName%in%data0$indexsnp ) 
# fib.magic.metal<-subset(fib.magic.metal0,Freq1>0.01 & Freq1<0.99 & HetDf>=4 | MarkerName%in%data0$indexsnp)

data0 <- fread("C:/Users/mbrown16/Desktop/Temp/EA.TestData.txt",data.table = FALSE)



# trait<-c("fg","fi")
# group<-c("MAGIC","AAGILE")
# pop<-c("CEU","YRI")

# t<-1	# 1:fg,	2:fib
# g<-2	# 1:MAGIC + CEU,	2:AAGILE + YRI




# for(g in 1:2){
  
  #if(g==2) data0$neg.log10.pval <- -log10(data0$P.value.index.aa)
  #if(g==1) data0$neg.log10.pval <- -log10(data0$Pval.index.ma)
  
  
  #for(i in 1:nrow(data0)){
  # for(i in c(26,52)){
    
    i <- which(data0$MarkerName==wh.snp)
    # print(i)
    
    snp0 <- data0[i,]$indexsnp
    pos0 <- data0[i,]$pos
    chr0 <- data0[i,]$chr
    gene0 <- data0[i,]$Locus
    #if(g==2) pval0<-data0[i,]$P.value.index.aa
    #if(g==1) pval0<-data0[i,]$Pval.index.ma
    pval0<-min(data0[i,]$P.value.index.aa, data0[i,]$Pval.index.ma)
    # qt0<-data0[i,]$QT
    # source0<-data0[i,]$Source
    
    tmp10<-fg.aagile.metal; tmp20<-fg.magic.metal;
    tmp30<-fib.aagile.metal; tmp40<-fib.magic.metal;
    
    tmp1 <- tmp10[tmp10$pos>=pos0-250000 & tmp10$pos<=pos0+250000 & !is.na(tmp10$pos) & tmp10$chr==chr0,]
    tmp2 <- tmp20[tmp20$pos>=pos0-250000 & tmp20$pos<=pos0+250000 & !is.na(tmp20$pos) & tmp20$chr==chr0,]
    tmp3 <- tmp30[tmp30$pos>=pos0-250000 & tmp30$pos<=pos0+250000 & !is.na(tmp30$pos) & tmp30$chr==chr0,]
    tmp4 <- tmp40[tmp40$pos>=pos0-250000 & tmp40$pos<=pos0+250000 & !is.na(tmp40$pos) & tmp40$chr==chr0,]
    
    
    # if(qt0=="fg" & g==2){ data1<-fg.aagile.metal; t<-1; print("fg aagile")}
    # if(qt0=="fib"& g==2){ data1<-fib.aagile.metal; t<-2; print("fib aagile") }
    # if(qt0=="fg" & g==1){ data1<-fg.magic.metal; t<-1; print("fg magic") }
    # if(qt0=="fib"& g==1){ data1<-fib.magic.metal; t<-2; print("fib magic") }
    
    # data1<-data1[,c("MarkerName","chr","pos","P.value")]
    data1 <- data0[,c("MarkerName","Chr","Pos","P.value")]
    data1$neg.log10.pval <- -log10(data1$P.value)
    
    data <- data1[data1$pos>=pos0-250000 & data1$pos<=pos0+250000 & !is.na(data1$pos) & data1$chr==chr0,]
    names(data) <- c("SNP","CHR","POS","P.value","neg.log10.pval")
    data <- data[order(data$P.value),]
    
    
    # print(paste(snp0,gene0,pop[g]))
    ld <- read.table(paste(ld.path,"ld",snp0,"_",pop,".txt",sep=""), header=T, sep="\t")
    
    print("Preparing data files")
    locus <- merge(data,ld,by.x="SNP",by.y="SNP",all=TRUE)
    if(is.na(locus[locus$SNP==snp0,"CHR"])) locus[locus$SNP==snp0,"CHR"] <- chr0
    if(is.na(locus[locus$SNP==snp0,"POS"])) locus[locus$SNP==snp0,"POS"] <- pos0
    locus <- locus[!is.na(locus$POS),]
    locus <- locus[order(-locus$neg.log10.pval),]
    for(r in 1:nrow(locus)){
      if(is.na(locus$RSQR[r])) {locus$RSQR[r] <- -1}
    } # for r
    
    
    # traitname<-toupper(trait[t])
    # groupname<-group[g]
    # ldpop<-pop[g]
    
    traitname <- ""
    groupname <- ""
    ldpop <- pop
    
    # plot y-axix #
    #topPval<-max(data$neg.log10.pval)
    # topPval <- -log10( min(min(tmp1$P.value), min(tmp2$P.value), min(tmp3$P.value), min(tmp4$P.value)) )
    topPval <- -log10( min( data$P.value, na.rm=T ))
    yrange<-topPval+topPval/10
    yaxis_scale<-ceiling(yrange/3)
    
    # if(g==2) plot.path<-"/restricted/projectnb/aagile/Jae_MANTRA/AAGILE_T2D/regional_plot/plot/AAGILE_YRI/"
    # if(g==1) plot.path<-"/restricted/projectnb/aagile/Jae_MANTRA/AAGILE_T2D/regional_plot/plot/MAGIC_CEU/"
    
    print("Making plot")
    #png(paste0(plot.path,toupper(trait[t]),"_",source0,"_chr",chr0,"_",gene0,"_",snp0,"_",groupname,"_",pop[g],".png"),width=6,height=4,res=300)
    # png(paste0(plot.path,toupper(trait[t]),"_",source0,"_chr",chr0,"_",gene0,"_",snp0,"_",groupname,"_",pop[g],".png"),width=3600,height=2400,res=300)
    png("test_assoc.png",width=720,height=480)
    make.plot(snp0,gene0,chr0,locus,yrange,yaxis_scale,traitname,groupname,ldpop)
    dev.off()
    
#  }  # for i
  
# } # for g




