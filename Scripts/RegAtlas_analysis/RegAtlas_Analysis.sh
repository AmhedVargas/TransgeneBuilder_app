##Analysis of promoters
## get current minimal data of RegAtlas

wget https://github.com/js2264/RegAtlas/raw/master/dashboard.Ahringer/releases/dashboard.Ahringer_v0.5.3/data/minimal-data.RData

R
##R analysis
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)



load("minimal-data.RData")
table(all$regulatory_class)
table(all$is.prom)
table(all$which.tissues)

promoters = all.deconv[which(all.deconv$is.prom),]

table(promoters$which.tissues)

#get names of tissues in the strangest way
tissues=names(table(promoters$which.tissues)[order(table(promoters$which.tissues))])

PromByT=list()


for(tissue in tissues){
	temp=promoters[which(promoters$which.tissues == tissue),c("locus","TSS.fwd","TSS.rev","which.tissues","uniqueWormBaseID")]
	temp=cbind(temp,LCAP[as.character(temp$uniqueWormBaseID),])
	PromByT=append(PromByT,list(temp))	
	names(PromByT)[which(tissue == tissues)]= tissue
}

###Example
dat=PromByT[["Germline"]]
plot(dat$Germline,(dat$Germline/(dat$Germline + dat$Neurons + dat$Muscle + dat$Hypod. + dat$Intest.)))

plot(log(dat$Germline),(dat$Germline/(dat$Germline + dat$Neurons + dat$Muscle + dat$Hypod. + dat$Intest.)))

tt=cbind(dat,(dat$Germline/(dat$Germline + dat$Neurons + dat$Muscle + dat$Hypod. + dat$Intest.)))

plot(tt$Germline,tt[,ncol(tt)], main="Germline promoters RNA-seq reads", ylab="Germline mRNA / All mRNA", xlab="mRNA reads in germline (RPKM)")

abline(h = median(tt[,ncol(tt)],na.rm=TRUE))
abline(v = median(tt[,"Germline"],na.rm=TRUE))

points(c(3223.355,3963.466,9543.230,3639.432),c(0.5254071,0.1388887,0.2042827,0.5031886), col=c("blue","green","gold","red"))

legend("topright",c("cpg-1","eft-3","rpl-7A","cpg-2"),col=c("blue","green","gold","red"),bty="n",pch="o")
dev.off()


###Do for other tissues
#Germline, Neurons, Muscle, Hypod., Intest.
pdf("Promoter_specificity-vs-expressivity_median.pdf")
par(mfrow=c(2,3))
for (tissue in c("Germline", "Neurons", "Muscle", "Hypod.", "Intest.")){
	dat = PromByT[[tissue]]
	tt=cbind(dat,(dat[,tissue]/(dat$Germline + dat$Neurons + dat$Muscle + dat$Hypod. + dat$Intest.)))
	plot(tt[,tissue],tt[,ncol(tt)], main=paste(tissue,"promoters RNA-seq reads"), ylab=paste(tissue,"mRNA / All mRNA"), xlab=paste("mRNA reads in", tissue, "(RPKM)"))
	abline(h = median(tt[,ncol(tt)],na.rm=TRUE))
	abline(v = median(tt[,tissue],na.rm=TRUE))
}
dev.off()

#########GGplot figures

plots=list()

for (tissue in c("Germline", "Neurons", "Muscle", "Hypod.", "Intest.")){
	dat = PromByT[[tissue]]
	tt=cbind(dat,(dat[,tissue]/(dat$Germline + dat$Neurons + dat$Muscle + dat$Hypod. + dat$Intest.)))
	dt= data.frame(Gene=tt$uniqueWormBaseID,Expressivity=tt[,tissue], Sensitivity=tt[,ncol(tt)])

	dt=unique(dt)

		filt=dt[which((dt[,2] > median(dt$Expressivity, na.rm=TRUE))&(dt[,3]> median(dt$Sensitivity, na.rm=TRUE))),]
	
	filt=filt[rev(order(filt[,2]))[1:5],]

p=	ggplot(dt, aes(x=Expressivity,y=Sensitivity)) + 
	geom_point(colour = "black", fill = "white", alpha = .5, shape=21, size = 1.5, stroke = 1.5) + 
	geom_hline(yintercept = median(dt$Sensitivity, na.rm=TRUE), color="grey", linetype = "dashed") +
	geom_vline(xintercept = median(dt$Expressivity, na.rm=TRUE), color="grey", linetype = "dashed")	+
	geom_point(data=filt, colour=brewer.pal(n = 5, name = 'Dark2')) + 
	annotate(geom="text", x=rep(max(dt$Expressivity, na.rm=TRUE) - .2*mean(dt$Expressivity, na.rm=TRUE),5), y=c(1,.95,.9,.85,.8), label=convID[as.character(filt$Gene),2],color=brewer.pal(n = 5, name = 'Dark2')) +
	ggtitle(paste(tissue,"promoters RNA-seq reads")) +
  xlab(paste("mRNA in", tissue, "(TPM)")) + ylab(paste(tissue,"mRNA / All mRNA"))

#ggsave(paste("Top5PromsIn",tissue,".png",sep=""),p, width = 16, height = 12, units = "cm")
plots=append(plots,list(p))
names(plots)[which(tissue == c("Germline", "Neurons", "Muscle", "Hypod.", "Intest."))]= tissue
}

##Ubiquitous
dat = PromByT[["Ubiq."]]
uno=apply(dat,1,function(x){as.numeric(x[6])/sum(as.numeric(x[6:10]))})
dos=apply(dat,1,function(x){as.numeric(x[7])/sum(as.numeric(x[6:10]))})
tres=apply(dat,1,function(x){as.numeric(x[8])/sum(as.numeric(x[6:10]))})
cuatro=apply(dat,1,function(x){as.numeric(x[9])/sum(as.numeric(x[6:10]))})
cinco=apply(dat,1,function(x){as.numeric(x[10])/sum(as.numeric(x[6:10]))})
med=apply(cbind(uno,dos,tres,cuatro,cinco), 1, function(x){median(x,na.rm=TRUE)})
susu=apply(dat, 1, function(x){sum(as.numeric(x[6:10]),na.rm=TRUE)})

	tt=cbind(dat,med)
	dt= data.frame(Gene=tt$uniqueWormBaseID,Expressivity=susu, Sensitivity=med)

	dt=unique(dt)

		filt=dt[which((dt[,3] > .15 )&(dt[,3] < .25 )),]
	
	filt=filt[rev(order(filt[,2]))[1:5],]
tissue="Ubiquitous"
p=	ggplot(dt, aes(x=Expressivity,y=Sensitivity)) + 
	geom_point(colour = "black", fill = "white", alpha = .5, shape=21, size = 1.5, stroke = 1.5) + 
	geom_hline(yintercept = .15, color="grey", linetype = "dashed") +
	geom_hline(yintercept = .25, color="grey", linetype = "dashed") +
	geom_vline(xintercept = median(dt$Expressivity, na.rm=TRUE), color="grey", linetype = "dashed")	+
	geom_point(data=filt, colour=brewer.pal(n = 5, name = 'Dark2')) + 
	annotate(geom="text", x=rep(max(dt$Expressivity, na.rm=TRUE) - .2*mean(dt$Expressivity, na.rm=TRUE),5), y=c(1,.95,.9,.85,.8), label=convID[as.character(filt$Gene),2],color=brewer.pal(n = 5, name = 'Dark2')) +
	ggtitle(paste(tissue,"promoters RNA-seq reads")) +
  xlab(paste("mRNA in", tissue, "(TPM)")) + ylab(paste(tissue,"mRNA / All mRNA"))
plots=append(plots,list(p))
names(plots)[length(names(plots))]= tissue

ggsave(paste("Top5PromsInAll.png",sep=""),ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],ncol=3, nrow=2), width = 36, height = 18, units = "cm")
#ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],ncol=3, nrow=2)

##Plot into single figure

###Now with normalized data


PromByT=list()


for(tissue in tissues){
	temp=promoters[which(promoters$which.tissues == tissue),c("locus","TSS.fwd","TSS.rev","which.tissues","uniqueWormBaseID")]
	temp=cbind(temp,LCAP_normalized[as.character(temp$uniqueWormBaseID),])
	PromByT=append(PromByT,list(temp))	
	names(PromByT)[which(tissue == tissues)]= tissue
}

plots=list()

for (tissue in c("Germline", "Neurons", "Muscle", "Hypod.", "Intest.")){
	dat = PromByT[[tissue]]
	tt=cbind(dat,(dat[,tissue]/(dat$Germline + dat$Neurons + dat$Muscle + dat$Hypod. + dat$Intest.)))
	dt= data.frame(Gene=tt$uniqueWormBaseID,Expressivity=tt[,tissue], Sensitivity=tt[,ncol(tt)])

	dt=unique(dt)

		filt=dt[which((dt[,2] > median(dt$Expressivity, na.rm=TRUE))&(dt[,3]> median(dt$Sensitivity, na.rm=TRUE))),]
	
	filt=filt[rev(order(filt[,2]))[1:5],]

p=	ggplot(dt, aes(x=Expressivity,y=Sensitivity)) + 
	geom_point(colour = "black", fill = "white", alpha = .5, shape=21, size = 1.5, stroke = 1.5) + 
	geom_hline(yintercept = median(dt$Sensitivity, na.rm=TRUE), color="grey", linetype = "dashed") +
	geom_vline(xintercept = median(dt$Expressivity, na.rm=TRUE), color="grey", linetype = "dashed")	+
	geom_point(data=filt, colour=brewer.pal(n = 5, name = 'Dark2')) + 
	annotate(geom="text", x=rep(max(dt$Expressivity, na.rm=TRUE) - .2*mean(dt$Expressivity, na.rm=TRUE),5), y=c(1,.95,.9,.85,.8), label=convID[as.character(filt$Gene),2],color=brewer.pal(n = 5, name = 'Dark2')) +
	ggtitle(paste(tissue,"promoters RNA-seq reads (normalized)")) +
  xlab(paste("mRNA in", tissue, "(TPM)")) + ylab(paste(tissue,"mRNA / All mRNA"))

#ggsave(paste("Top5PromsIn",tissue,".png",sep=""),p, width = 16, height = 12, units = "cm")
plots=append(plots,list(p))
names(plots)[which(tissue == c("Germline", "Neurons", "Muscle", "Hypod.", "Intest."))]= tissue
}


dat = PromByT[["Ubiq."]]
uno=apply(dat,1,function(x){as.numeric(x[6])/sum(as.numeric(x[6:10]))})
dos=apply(dat,1,function(x){as.numeric(x[7])/sum(as.numeric(x[6:10]))})
tres=apply(dat,1,function(x){as.numeric(x[8])/sum(as.numeric(x[6:10]))})
cuatro=apply(dat,1,function(x){as.numeric(x[9])/sum(as.numeric(x[6:10]))})
cinco=apply(dat,1,function(x){as.numeric(x[10])/sum(as.numeric(x[6:10]))})
med=apply(cbind(uno,dos,tres,cuatro,cinco), 1, function(x){median(x,na.rm=TRUE)})
susu=apply(dat, 1, function(x){median(as.numeric(x[6:10]),na.rm=TRUE)})

	tt=cbind(dat,med)
	dt= data.frame(Gene=tt$uniqueWormBaseID,Expressivity=susu, Sensitivity=med)

	dt=unique(dt)

		filt=dt[which((dt[,3] > .15 )&(dt[,3] < .25 )),]
	
	filt=filt[rev(order(filt[,2]))[1:5],]
tissue="Ubiquitous"
p=	ggplot(dt, aes(x=Expressivity,y=Sensitivity)) + 
	geom_point(colour = "black", fill = "white", alpha = .5, shape=21, size = 1.5, stroke = 1.5) + 
	geom_hline(yintercept = .15, color="grey", linetype = "dashed") +
	geom_hline(yintercept = .25, color="grey", linetype = "dashed") +
	geom_vline(xintercept = median(dt$Expressivity, na.rm=TRUE), color="grey", linetype = "dashed")	+
	geom_point(data=filt, colour=brewer.pal(n = 5, name = 'Dark2')) + 
	annotate(geom="text", x=rep(max(dt$Expressivity, na.rm=TRUE) - .2*mean(dt$Expressivity, na.rm=TRUE),5), y=c(1,.95,.9,.85,.8), label=convID[as.character(filt$Gene),2],color=brewer.pal(n = 5, name = 'Dark2')) +
	ggtitle(paste(tissue,"promoters RNA-seq reads (normalized)")) +
  xlab(paste("mRNA in", tissue, "(TPM)")) + ylab(paste(tissue,"mRNA / All mRNA"))
plots=append(plots,list(p))
names(plots)[length(names(plots))]= tissue


ggsave(paste("Top5PromsInAll_Normalized.png",sep=""),ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],ncol=3, nrow=2), width = 36, height = 18, units = "cm")



###Gene-Tissue-Operonic-associated_promoter-Read_in_tissue

tissues=names(table(promoters$which.tissues)[order(table(promoters$which.tissues))])

PromByT=list()


for(tissue in tissues){
	temp=promoters[which(promoters$which.tissues == tissue),c("locus","TSS.fwd","TSS.rev","which.tissues","uniqueWormBaseID")]
	temp=cbind(temp,LCAP[as.character(temp$uniqueWormBaseID),])
	PromByT=append(PromByT,list(temp))	
	names(PromByT)[which(tissue == tissues)]= tissue
}


for (tissue in c("Germline", "Neurons", "Muscle", "Hypod.", "Intest.")){
	dat = PromByT[[tissue]]
	tt=cbind(dat,(dat[,tissue]/(dat$Germline + dat$Neurons + dat$Muscle + dat$Hypod. + dat$Intest.)))
	dt= data.frame(Gene=tt$uniqueWormBaseID,Expressivity=tt[,tissue], Sensitivity=tt[,ncol(tt)])

	dt=unique(dt)

		filt=dt[which((dt[,2] > median(dt$Expressivity, na.rm=TRUE))&(dt[,3]> median(dt$Sensitivity, na.rm=TRUE))),]
	
	filt=filt[rev(order(filt[,2]))[1:5],]

idx=c()
for(i in 1:nrow(filt)){
  idx=append(idx,which(as.character(filt[i,1])==as.character(tt$uniqueWormBaseID)))
  }

ntt = tt[idx,c("locus","which.tissues","uniqueWormBaseID")]
ntt = cbind(ntt,convID[as.character(ntt$uniqueWormBaseID),2])

ntt = cbind(ntt,ATAC[as.character(ntt$locus),])
colnames(ntt)[4]="Gene"

write.table(ntt,paste(tissue,"-Top5_info.tsv",sep=""),sep="\t",row.names=FALSE)
}

##Ubiquitous
dat = PromByT[["Ubiq."]]
uno=apply(dat,1,function(x){as.numeric(x[6])/sum(as.numeric(x[6:10]))})
dos=apply(dat,1,function(x){as.numeric(x[7])/sum(as.numeric(x[6:10]))})
tres=apply(dat,1,function(x){as.numeric(x[8])/sum(as.numeric(x[6:10]))})
cuatro=apply(dat,1,function(x){as.numeric(x[9])/sum(as.numeric(x[6:10]))})
cinco=apply(dat,1,function(x){as.numeric(x[10])/sum(as.numeric(x[6:10]))})
med=apply(cbind(uno,dos,tres,cuatro,cinco), 1, function(x){median(x,na.rm=TRUE)})
susu=apply(dat, 1, function(x){sum(as.numeric(x[6:10]),na.rm=TRUE)})

	tt=cbind(dat,med)
	dt= data.frame(Gene=tt$uniqueWormBaseID,Expressivity=susu, Sensitivity=med)

	dt=unique(dt)

		filt=dt[which((dt[,3] > .15 )&(dt[,3] < .25 )),]
	
	filt=filt[rev(order(filt[,2]))[1:5],]
tissue="Ubiquitous"

idx=c()
for(i in 1:nrow(filt)){
  idx=append(idx,which(as.character(filt[i,1])==as.character(tt$uniqueWormBaseID)))
  }

ntt = tt[idx,c("locus","which.tissues","uniqueWormBaseID")]
ntt = cbind(ntt,convID[as.character(ntt$uniqueWormBaseID),2])

ntt = cbind(ntt,ATAC[as.character(ntt$locus),])
colnames(ntt)[4]="Gene"

write.table(ntt,paste(tissue,"-Top5_info.tsv",sep=""),sep="\t",row.names=FALSE)



quit()

##Sequences
mkdir Fastas

for file in `ls *.tsv`; do echo ${file%.tsv}; perl -pe 's/\"//g' ${file} | grep -v "locus" | awk -F"\t" '{split($1,coso,"_chr"); print coso[2]"_"$3";"$4}' | awk -F"_" '{OFS="\t"; print $1,$2,$3,$4}' | sort -k1,1 -k2,2n > Fastas/${file%.tsv}.bed ;done

cd Fastas


cp /home/velazqam/Documents/Projects/Fragmenst/c_elegans.PRJNA13758.WS275.genomic.fa.gz .

##Prepare genome for next steps
gzip -d c_elegans.PRJNA13758.WS275.genomic.fa.gz

####Be aware that directoniality has been lost
for file in `ls *.bed`; do bedtools getfasta -fi c_elegans.PRJNA13758.WS275.genomic.fa -bed ${file} -name > ${file%.bed}.fasta; done

cat *.fasta > mutli_promoter.fas

cd ..

mkdir Prom-UTR
cd Prom-UTR

##Download from reg atlas
#tissue-specific.ATAC-seq.dataset.txt
#tissue-specific.RNA-seq.dataset.txt
egrep "promoter|enhancer" tissue-specific.ATAC-seq.dataset.txt | perl -pe 's/^chr//g' | awk -F"\t" '{OFS="\t"; print $1,$2,$3,$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11}' > ATAC_prom-enh.bed

#DOwnload data
 wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.annotations.gff3.gz
 wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.genomic.fa.gz

zcat c_elegans.PRJNA13758.WS275.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=" '{print $1"\t"$4"\t"$5"\t"$12"\t.\t"$7}' | perl -pe 's/Transcript://g' | awk -F"\t" '{OFS="\t"; split($4,trans,","); for(i=1;i<=length(trans);i++){print $1,$2,$3,trans[i],$5,$6}}' | sort -k1,1 -k2,2n > CDS_WS275.bed

bedtools intersect -v -a ATAC_prom-enh.bed -b CDS_WS275.bed | awk -F"\t" '{OFS="\t"; split($4,info,";"); print info[1]"_chr"$1"_"$2"_"$3,info[1],info[2],info[3],info[4],info[5],info[6],info[7],info[8]}' > ATAC_prom-enh_noG.txt

cd ..

R
##R analysis
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

load("minimal-data.RData")

ftab=read.table("Prom-UTR/ATAC_prom-enh_noG.txt",sep="\t",header=F, stringsAsFactors=F)

tab=cbind(all[as.character(ftab[,1]),],ATAC[as.character(ftab[,1]),])

tab=cbind(tab,apply(tab[,c("Germline","Neurons","Muscle","Hypod.","Intest.")],1,sum))

colnames(tab)[ncol(tab)]="Sum.ATAC"

divis=t(apply(tab[,c("Germline","Neurons","Muscle","Hypod.","Intest.","Sum.ATAC")],1,function(x){return(c(x/x[length(x)]))}))

colnames(divis)=paste(colnames(divis),"-fraction",sep="")

tab=cbind(tab,divis[,-ncol(divis)])

boxplot(tab[,(ncol(tab)-4):ncol(tab)], ylab="Tissue/Sum ATAC-reads", main="ATAC-seq signal per tissue per not coding region")

###Use .85 as proxy
abline(h=.85, col = "red")

genelists=list()
i=1
for (cn in colnames(divis)[-ncol(divis)]){
	gname=tab[which(tab[,cn] > .85),"WormBaseID"]
	gname=unique(unlist(strsplit(as.character(gname),",")))
	genelists=append(genelists,list(gname))
	names(genelists)[i]=cn
	i=i+1
}


###Prep for shiny visualization
##Filt for only ATCs with associated genes

tab=tab[-which(as.character(tab$WormBaseID) == "."),]
genes=lcap.dt[,c(1:12)]
genes=cbind(genes,apply(genes[,c("Germline","Neurons","Muscle","Hypod.","Intest.")],1,sum))
colnames(genes)[ncol(genes)]="SumRNA"
divis=t(apply(genes[,c("Germline","Neurons","Muscle","Hypod.","Intest.","SumRNA")],1,function(x){return(c(x/x[length(x)]))}))
colnames(divis)=paste(colnames(divis),"-fraction",sep="")
genes=cbind(genes,divis[,-ncol(divis)])

##Obtain fastas to make app
#write.table(tab,"temp.txt",sep="\t", col.names=F, row.names=F, quote=F)
#quit()
#cd Fastas; cp ../temp.txt .
#perl -pe 's/^chr//g' temp.txt | awk -F"\t" '{OFS="\t"; if($6 == "."){$6="+"}; print $1,$2,$3,$4,$5,$6}' > temp.bed
#bedtools getfasta -fi c_elegans.PRJNA13758.WS275.genomic.fa -bed temp.bed -name -s -tab | perl -pe 's/\(/\t/g' | perl -pe 's/\)//g' > temp.fastab
fastab=read.table("Fastas/temp.fastab",sep="\t",header=F, stringsAsFactors=F)
rownames(fastab) = as.character(fastab[,1])

###Median ATAC
