# Script generalized from 'make.plot_P.value_forManning.R' script written by Jaeyoung Hong

LOCALTEST <- FALSE
if (LOCALTEST) {
    ## Load libraries and source files
    library(data.table)
    setwd("/home/mbrown/DNAnexus/assocplot/resources/")
    source("modified.assocplot_P.value.R")
    
    ## Setup paths
    datafile <- "/home/mbrown/DNAnexus/data/mb_test_data_chr8.Rda"
    # index.path
    ld.path <- "/home/mbrown/DNAnexus/"
    # ld.path <- "/home/mbrown/DNAnexus/test.ped.LD"
    # ld.path <- "/home/mbrown/DNAnexus/ppp1r3b_locus_EU.ped.LD.txt"
    # plot.path <- getwd()

    ## Setup parameters
    snp.name <- "Name"
    gene.name <- "gene"
    chr.name <- "chr"
    pos.name <- "pos"
    p.name <- "p"
    freq.name <- "maf"
    df.name <- "HetDF"
    wh.snp <- "8:9044912"
    pop <- "CEU"
    traitname <- ""
    groupname <- ""

} else {
    
    ## Parse the arguments
    args <- commandArgs(trailingOnly=TRUE)
    args <- paste(args,collapse=" ")
    print(args)
    args <- unlist(strsplit(args,";"))
    print(args)
    
    ## Mandatory parameters
    datafile <- args[1]
    ld.path <- args[2]

    ## Optional parameters
    ld.type <- args[3]
    snp.name <- args[4]
    gene.name <- args[5]
    chr.name <- args[6]
    pos.name <- args[7]
    p.name <- args[8]
    freq.name <- args[9]
    df.name <- args[10]
    wh.snp <- args[11]
    region.width <- as.numeric(args[12])
    pop <- args[13]
    chr.num <- args[14]
    wh.gene <- args[15]
    wh.region <- args[16]
    traitname <- args[17]
    groupname <- args[18]
    output.type <- args[19]
    output.args <- args[20]

    ## Load libraries and source files
    install.packages("/chron_2.3-47.tar.gz", repos=NULL, type="source")
    install.packages("/data.table_1.9.6.tar.gz", repos=NULL, type="source")
    library(data.table)
    source("/modified.assocplot_P.value.R")
}

# Functions
get.delim <- function(f) {
        d <- readLines(f,n=1)
        delim <- ifelse(grepl(" ", d)," ",ifelse(grepl(",",d),",","\t"))
        return(delim)
}

## Load data -- current datatype is seqMeta file
print("Reading in data...")

# Check filetype
if (grepl("Rda$",datafile)) {
        load(datafile)
} else {
	delim <- get.delim(datafile)
	sing <- read.table(datafile, header=T, as.is=T, sep=delim)
}

# data0 <- fread(datafile,data.table = FALSE)

    # i <- which(data0$MarkerName==wh.snp)
    i <- which(sing[,snp.name]==wh.snp)

    # snp0 <- data0[i,]$indexsnp
    # pos0 <- data0[i,]$pos
    # chr0 <- data0[i,]$chr
    # gene0 <- data0[i,]$Locus
    # pval0<-min(data0[i,]$P.value.index.aa, data0[i,]$Pval.index.ma)
    
    snp0 <- sing[i,snp.name]
    pos0 <- sing[i,pos.name]
    chr0 <- sing[i,chr.name]
    gene0 <- sing[i,gene.name]
    pval0 <- sing[i,p.name]
    
    # data1 <- data0[,c("MarkerName","Chr","Pos","P.value")]
    # data1$neg.log10.pval <- -log10(data1$P.value)
    data1 <- sing[,c(snp.name,chr.name,pos.name,p.name)]    
    data1$logp <- -log10(sing[,p.name])

    # data <- data1[data1$pos>=pos0-250000 & data1$pos<=pos0+250000 & !is.na(data1$pos) & data1$chr==chr0,]
    # names(data) <- c("SNP","CHR","POS","P.value","neg.log10.pval")
    # data <- data[order(data$P.value),]
    data <- data1[data1[,pos.name]>=pos0-250000 & data1[,pos.name]<=pos0+250000 & !is.na(data1[,pos.name]) & data1[,chr.name]==chr0,]
    names(data) <- c("SNP","CHR","POS","P.value","neg.log10.pval")
    data <- data[order(data$P.value),]

    # ld <- read.table(paste(ld.path,"ld",snp0,"_",pop,".txt",sep=""), header=T, sep="\t")    
    # Read in the LD data
    ld <- read.table(ld.path, header=T, as.is=T) 
    if (ld.type=="PLINK") {
      ld1 <- ld[,c("SNP_A","SNP_B","R2")]
      names(ld1) <- c("L1","L2","r.2")
    } else {
      names(ld)[3] <- "D"   
      ld1 <- ld[ld$L1==snp0,c("L1","L2","r.2")]
      ld1 <- rbind(c(ld1$L1[1],ld1$L1[1],1),ld1)
    }

    print("Preparing data files")
    locus <- merge(data,ld1,by.x="SNP",by.y="L2")
    names(locus)[which(names(locus)=="r.2")] <- "RSQR"

    if(is.na(locus[locus$SNP==snp0,"CHR"])) locus[locus$SNP==snp0,"CHR"] <- chr0
    if(is.na(locus[locus$SNP==snp0,"POS"])) locus[locus$SNP==snp0,"POS"] <- pos0
    locus <- locus[!is.na(locus$POS),]
    locus <- locus[order(-locus$neg.log10.pval),]
    for(r in 1:nrow(locus)){
      if(is.na(locus$RSQR[r])) {locus$RSQR[r] <- -1}
    } # for r
    
    ldpop <- pop
    
    # plot y-axix #
    topPval <- -log10( min( data$P.value, na.rm=T ))
    yrange<-topPval+topPval/10
    yaxis_scale<-ceiling(yrange/3)
    
    if (output.args!="" & !is.na(output.args) & !is.null(output.args) & output.args!="NA" & output.args!="NULL") {
       output.args <- paste(",",output.args)
    } else {
       if (output.type=="pdf") {
          output.args <- ", width=7,height=4.67"
       }
       else {
          output.args <- ", width=720,height=480,units=\"px\""
       }
    }

    print("Making plot")
    # png("assoc.png",width=720,height=480)
    eval(parse(text=paste0(output.type,"(\"assoc.",output.type,"\"",output.args,")")))
    make.plot(snp0,gene0,chr0,locus,region.width,yrange,yaxis_scale,traitname=traitname,groupname=groupname,ldpop,ld.type,LOCALTEST)
    dev.off()
    



