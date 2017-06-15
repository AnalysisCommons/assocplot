make.plot<-function(snp,gene,chr,locus,region.width,range,p.scale,traitname,groupname,ldpop,ld.type,LOCALTEST) {

hit<-locus[locus$SNP==snp,]

#			#
# size of the region	#
#			#
w <- region.width/2
# works out to +/-250k, +/-270k for default region.width of 500kb
min.pos<-hit$POS-(w+w/12.5)
max.pos<-hit$POS+(w+w/12.5)
min.snp.pos<-hit$POS-w
max.snp.pos<-hit$POS+w
size.pos<-max.pos-min.pos
#center.pos<-min.pos+(size.pos/2)
center.pos<-hit$POS
center.100kb.pos<-round(center.pos/1000)*1000
offset.100kb.pos<-round((size.pos/3)/1000)*1000

#			#
# range of y-axis	#
#			#
# this dedicates 33% of the yaxis to the genes, labels, recomb rate
offset<-(range*4/3)-range
big.range<-range+offset 
ystart.gene<- -offset
ystart.recomb<- -offset+(big.range/8)

#			#
# recombination rate #
#			#
recomb.dir <- ifelse(LOCALTEST, "hapmap_recomb_rates_2008-03_rel22_B36/","/hapmap_recomb_rates_2008-03_rel22_B36/")
if(chr>=1 & chr<=22){
	# recomb <- read.table(paste("/restricted/projectnb/jae-thesis/hapmap/phasing/2007-08_rel22/recomb_rates/latest/genetic_map_chr",chr,"_b36.txt",sep=""),header=T)
	recomb <- read.table(paste0(recomb.dir,"genetic_map_chr",chr,"_b36.txt"),header=T,as.is=T)
	keep.recomb <- subset(recomb,recomb[,1]>min.pos & recomb[,1]<max.pos)
}
if(chr==23){
	stop("Chr==23 not implemented at this time.")
	# recomb <- read.table(paste("/restricted/projectnb/aagile/hapmap_genotypes_2008-10_phaseII_fwd_strand_non-redundant/plink_format/recomb_rates/2006-10_rel21_phaseI+II/genetic_map_chrX_rel21_phaseI+II_all.txt",sep=""),header=T)
	recomb <- read.table(paste0(recomb.dir,"genetic_map_chr",chr,"_b36.txt"),header=T,as.is=T)
	keep.recomb <- subset(recomb,recomb[,1]>min.pos & recomb[,1]<max.pos)
}

# And varLD rates #
# varlds <- read.table(paste0("/restricted/projectnb/aagile/varld/chr/Chr",chr,".txt"),header=T, sep="\t")
# keep.varlds <- subset(varlds,varlds[,3]>min.pos & varlds[,3]<max.pos)
# keep.varlds1<-keep.varlds[,c(3,4)]
# keep.varlds1[,2]<-ifelse(keep.varlds1$raw_score>36.0929,keep.varlds1$raw_score,NA)



#			 #
# genes in the region #
#			 #
genelist.file <- ifelse(LOCALTEST,"hg19_refGene_from_annovar_but_reformatted.csv","/hg19_refGene_from_annovar_but_reformatted.csv")
#genelist.file <- ifelse(LOCALTEST,"hg19_to_18_annotation.csv","/hg19_to_18_annotation.csv")
#genelist<-read.table("/restricted/projectnb/glycemic/known_genes/known_genes_b36.csv",header=T,sep=",")
genelist<-read.table(genelist.file,header=T,sep=",")
genelist$CHR <- as.character(genelist$CHR)
genelist[genelist$CHR=="X","CHR"] <- "23"
genelist<-genelist[genelist$CHR==chr,]
genelist <- subset(genelist, !is.na(genelist$GENE))

genelist<-subset(genelist,(genelist$START>min.snp.pos & genelist$START<max.snp.pos) | (genelist$STOP>min.snp.pos & genelist$STOP<max.snp.pos))
genes.in.locus<-genelist[genelist$CHR==chr,]
genes.in.locus<-genes.in.locus[order(genes.in.locus$START),]
y0<-c(rep(c(0,1,2,3,4)-0.5,100))
y<-y0[1:nrow(genes.in.locus)]
if(nrow(genes.in.locus)>0){
	genes.in.locus$y<-y
}
#print(head(genes.in.locus))



#		#
# markers	#
#		#
markers.in.strong.ld<-subset(locus,(SNP != snp & locus$RSQR >= 0.8))
markers.in.moderate.ld<-subset(locus,(SNP != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8))
markers.in.weak.ld<-subset(locus,(SNP != snp & locus$RSQR >= 0.3 & locus$RSQR < 0.5))
markers.not.in.ld<-subset(locus,(SNP != snp & locus$RSQR < 0.3 & locus$RSQR >= 0))
markers.in.no.ld<-subset(locus,(SNP != snp & locus$RSQR == -1))

print(paste("total markers to print =",nrow(markers.in.strong.ld)+nrow(markers.in.moderate.ld)+nrow(markers.in.weak.ld)+nrow(markers.in.no.ld)+nrow(markers.not.in.ld)+1))
print(paste("total markers in the datasets =",nrow(locus)))

par(mar=c(5,4,4,4))

# Assign main title for plot
main.title <- ifelse(traitname=="" & groupname=="", paste0(gene,": ",snp,"  (LD: @)"),
			paste0(gene,": ",snp,"  (",groupname," ",traitname,", LD: @)"))
ld.sub <- ifelse(ld.type=="Haploview", paste0("HapMap2 ",ldpop), "PLINK")
main.title <- sub("@", ld.sub, main.title)

#								#
# start plot with recombination rate (in background)	#
#								#
par(oma=c(1,1,1,1))

plot(keep.recomb[,1],ystart.recomb+((keep.recomb[,2]/60)*(6*big.range/8)),type="l",col="lightblue",lwd=1,xlim=c(min.pos, max.pos),ylim=c(-offset,range),
xlab="",ylab="",main=main.title,font.main=3,cex=1.5,axes=F)

# And varLD rates plot #
# points(keep.varlds[,3],ystart.recomb+(abs(keep.varlds[,4])/60)*(6*big.range/8),type="l",col="palegreen",pch=46)
# points(keep.varlds1[,1],ystart.recomb+(abs(keep.varlds1[,2])/60)*(6*big.range/8),type="l",col="brown",pch=46)





#				#
# axes, titles and legends	#
#				#
mtext(paste("Chromosome",chr,"position (kb)",sep=" "),side=1,line=2.5,cex=1.5)
axis(1,at=c(center.100kb.pos-offset.100kb.pos,center.100kb.pos,center.100kb.pos+offset.100kb.pos),labels=c((center.100kb.pos-offset.100kb.pos)/1000,
center.100kb.pos/1000,(center.100kb.pos+offset.100kb.pos)/1000),las=1,mgp=c(3,.7,0)) 

axis(2,at=seq(0,range,p.scale),labels=seq(0,range,p.scale),las=1,mgp=c(3,.7,0)) 
mtext("-log10(P.value)",side=2,at=(range/2),line=2,cex=1.5)

axis(4,at=c(ystart.recomb,ystart.recomb+(big.range/4),ystart.recomb+(2*big.range/4),ystart.recomb+(3*big.range/4)),labels=c("0","20","40","60"),las=1,
mgp=c(3,.7,0))
mtext("Recombination rate (cM/Mb) | varLD",side=4,at=(-offset+big.range/2),line=2,cex=1.5)

box()
lines(c(min.pos, max.pos),c(0,0),lty="dotted",lwd=1,col="black")





#			#
# plot the markers	#
#			#
points(markers.in.no.ld$POS,markers.in.no.ld$neg.log10.pval,pch=21,cex=1.0,bg="white")
points(markers.not.in.ld$POS,markers.not.in.ld$neg.log10.pval,pch=21,cex=1.0,bg="skyblue")
points(markers.in.weak.ld$POS,markers.in.weak.ld$neg.log10.pval,pch=21,cex=1.5,bg="brown")
points(markers.in.moderate.ld$POS,markers.in.moderate.ld$neg.log10.pval,pch=21,cex=1.7,bg="orange")
points(markers.in.strong.ld$POS,markers.in.strong.ld$neg.log10.pval,pch=21,cex=1.9,bg="red")

#			#
# this is the hit	#
#			#
if(!is.na(hit$neg.log10.pval)) points(hit$POS,hit$neg.log10.pval,pch=23,cex=2.1,bg="red")

#								#
# top snp in the region unless it's the index snp	#
#								#

topsnp <- locus[max(locus$neg.log10.pval)==locus$neg.log10.pval,"SNP"]


maxlogpval<-locus[locus$SNP==topsnp,]$neg.log10.pval

if(!(snp==topsnp) & !is.na(topsnp)){
	pos1<-locus[locus$SNP==topsnp & !is.na(locus$neg.log10.pval),]$POS; 
	pval1<-locus[locus$SNP==topsnp & !is.na(locus$neg.log10.pval),]$neg.log10.pval; 
	plab1<-round(pval1,4); lpos1=1; loft1=1.5; spos1=3; soft1=1; col1="red";

#	points(pos1,pval1,pch=21,cex=3,bg=col1)
	text(pos1,pval1,paste("topsnp:",locus[locus$POS==pos1 & !is.na(locus$neg.log10.pval),1]),pos=spos1,offset=soft1)
#	text(pos1,pval1,plab1,pos=lpos1,offset=loft1)
 }

#			#
# plot the genes	#
#			#
if(nrow(genes.in.locus)>0){

for ( i in 1:nrow(genes.in.locus) ) { 
	if ( genes.in.locus[i,]$STRAND=="+" ) {
		arrows(max(genes.in.locus[i,]$START,min.snp.pos),-offset+genes.in.locus[i,]$y*(offset/6),min(genes.in.locus[i,]$STOP,max.snp.pos),-offset+genes.in.locus[i,]$y*(offset/6),length=0.03,lwd=1,code=2,lty="solid",col="darkgreen")
		} else {		
		arrows(max(genes.in.locus[i,]$START,min.snp.pos),-offset+genes.in.locus[i,]$y*(offset/6),min(genes.in.locus[i,]$STOP,max.snp.pos),-offset+genes.in.locus[i,]$y*(offset/6),length=0.03,lwd=1,code=1,lty="solid",col="darkgreen")
		}
	text.start <- max( genes.in.locus[i,]$START, min.snp.pos)
	new.gene.size <- min(genes.in.locus[i,]$STOP,max.snp.pos) - text.start
	text( text.start+new.gene.size/2,-offset+genes.in.locus[i,]$y*(offset/6)+(offset/12),labels=genes.in.locus[i,]$GENE,font=3,cex=0.8)
	}
}

#				#
# legend and reference	#
#				#
legend("topright",c(expression(R^2*" 0.8-1.0"),expression(R^2*" 0.5-0.8"),expression(R^2*" 0.3-0.5"),expression(R^2*" 0.0-0.3"),"Not available"),
fill=c("red","orange","brown","skyblue","white"),cex=1.1,bty="n",title="LD",inset=0.02)

}
