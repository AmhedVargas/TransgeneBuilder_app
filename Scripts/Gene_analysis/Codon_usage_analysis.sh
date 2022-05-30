###Redo analysis for Gene Builder

###First Process data (already but place code)

##Main figures
#Codon usage in different tissues
#Location of first intron (size distribution of exons)
#Identification of minimal but strongest promoters based on Ahringer data sets
##Specifity vs strenght plot


##TO CONSIDER:
###Identification of 5's and 3's so you can add to transgene builder
###Identification of minimal but strongest promoters based on Ahringer data sets
###Identification of 3' UTRS from ribosomal genes
###Identification of short promoters corresponding to neurons, muscle, hypoderm, intestine from the Ahringerlab datasets.
###Classifier algorithm where you can expand
###C. elegans, are theb synthesis of proteins equal in all tissues?
####http://104.131.81.59/data/sequence_lib_scores.db

###Working directory####
cd /home/velazqam/Documents/Projects/CodonOptimization/Analysis

mkdir codon_usage
cd codon_usage

############Codon analysis###############


###First get sequences again from wormbase###
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.WS276.CDS_transcripts.fa.gz
zcat c_elegans.PRJNA13758.WS276.CDS_transcripts.fa.gz | awk '{if(/>/){print "\n"$0}else{printf $0}}END{print ""}' | tail -n+2 > Celeg_CDS.fasta
##Get CDS of largest isoform
awk '{if(/^>/ ){split($0,gene,"=")}else{print gene[2]"\t"$0}}' Celeg_CDS.fasta | awk -F"\t" '{if(array[$1] == ""){array[$1] = $2; size[$1]=length($2)}else{if(size[$1] < length($2)){array[$1]=$2; size[$1]=length($2)} }} END{for(gene in array){print gene"\t"array[gene]}}' > LongestCDS_genetab_eleg.tsv
##Copy refined lists coming from Ahringers paper
cp /home/velazqam/Documents/Projects/CodonOptimization/combined_fasta/*.tsv .
##Obtain their fasta
for file in `ls Top500_*.tsv`; do echo ${file%.tsv}; awk '{if(array[$1] eq == ""){array[$1]=$2}else{print $0"\t"array[$1]}}' LongestCDS_genetab_eleg.tsv ${file} | awk -F"\t" '{print ">"$1"\n"$5}' > ${file%.tsv}.fasta; done
#Do it again but now combine all fastas into a single sequence for further analysis
for file in `ls Top500_*.tsv`; do echo ${file%.tsv}; awk '{if(array[$1] eq == ""){array[$1]=$2}else{print $0"\t"array[$1]}}' LongestCDS_genetab_eleg.tsv ${file} | awk -F"\t" 'BEGIN{print ">combined"}{printf $5}END{printf "\n"}' > comb_${file%.tsv}.fasta; done

##Now R to plot
R

library(GeneGA)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

aa=read.fasta("comb_Top500_ubiq_lcapdev.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

ubiq_amino=amino
ubiqW=w_value

aa=read.fasta("comb_Top500_soma_lcapdev.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

soma_amino=amino
somaW=w_value

aa=read.fasta("comb_Top500_germline_lcapdev.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

germ_amino=amino
germW=w_value

aa=read.fasta("comb_Top500_neurons_lcapdev.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

neur_amino=amino
neurW=w_value

data=data.frame(Gemline=germW,Neurons=neurW,Soma=somaW,Ubiq=ubiqW)
data[,5]=rownames(data)
colnames(data)[5]="Amino"

mm=melt(data,id="Amino")
colnames(mm)[2]="Expression"

p=ggplot(mm,aes(Amino,value,colour=Expression)) + geom_point() + ggtitle("Codon usage by expression group") + labs(y="Frequency", x = "Codon")
p
ggsave("Codon_Frequencies_top500.png", p, width = 40, height = 12, units = "cm")

df1 <- lapply(data.frame((data[,-5])), function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
})


df1=data.frame(df1)
mt1=as.matrix(df1)
heatmap(mt1)

codons= data[,5]
aminos= sapply(data[,5], function(x) translate(s2c(x)))

stab=cbind(1:length(codons),codons,aminos)
##Now sort by AA and codon
stab=stab[order(stab[,3]),]

mt1=mt1[as.numeric(stab[,1]),]

colnames(mt1)=c("Germline","Neurons","Somatic","Ubiquitous")

heatmap(mt1, Rowv =NA, labRow= paste(stab[,2]," (",stab[,3],")",sep=""), main="Clustering of highly expressed gene sets by codon usage",
	xlab="Gene set", ylab="Codon Frequencies", col= c(colorRampPalette(brewer.pal(8, "Oranges"))(11)) )

legend(x="bottomright", legend=c(as.character(seq(0,1,.1))), 
     fill=c(colorRampPalette(brewer.pal(8, "Oranges"))(11)) )

png("Codon_heatmap_top500.png", width = 1200, height = 800,res=100)
heatmap(t(mt1), Colv =NA, labCol= paste(stab[,2]," (",stab[,3],")",sep=""), main="Clustering of highly expressed gene sets by codon usage",
	ylab="", xlab="Codon Frequencies", col= c(colorRampPalette(brewer.pal(3, "RdBu"))(6)) )

legend(x="bottomleft", ncol=3, bty="n", legend=rev(as.character(seq(0,1,.2))), 
     fill=rev(colorRampPalette(brewer.pal(3, "RdBu"))(6)) )
dev.off()



tiff("Codon_heatmap_top500.tiff", width = 1200, height = 800)
heatmap(t(mt1), Colv =NA, labCol= paste(stab[,2]," (",stab[,3],")",sep=""), main="Clustering of highly expressed gene sets by codon usage",
	ylab="", xlab="Codon Frequencies", col= c(colorRampPalette(brewer.pal(3, "RdBu"))(6)) )

legend(x="bottomleft", ncol=3, bty="n", legend=rev(as.character(seq(0,1,.2))), 
     fill=rev(colorRampPalette(brewer.pal(3, "RdBu"))(6)) )
dev.off()

quit()

###Remove protein sequences with similarity usind CD-HIT*requires installation
##Convert first into protein sequences
mkdir Remove_outliers

for file in `ls Top500*.fasta`; do echo ${file%.fasta}; transeq ${file} Remove_outliers/${file%.fasta}_peptide.fasta; done

cd Remove_outliers/
##Run CD-HIT
for file in `ls Top500*.fasta`; do echo ${file%.fasta}; cd-hit -i ${file} -o ${file%.fasta}_representatives.fasta -c .8 -n 5 -M 32000 -T 20; done
#Germline
#495  finished        485  clusters
#Neurons
#495  finished        492  clusters
#Somatic
#496  finished        489  clusters
#Ubiquitous
#500  finished        494  clusters

##Get lists of IDS
for file in `ls Top500*_peptide_representatives.fasta`; do echo ${file%.fasta}; grep ">" ${file} | perl -pe 's/^>//g' | awk -F"_" '{print $1}' > ${file%_peptide_representatives.fasta}_representative_ids.txt; done

##Perform new analyses by getting again sequences by ID and combining them into single fasta
cd ..
mkdir filtered
cd -
for file in `ls Top500_*_representative_ids.txt`; do echo ${file%.txt}; awk '{if(array[$1] eq == ""){array[$1]=$2}else{print $0"\t"array[$1]}}' ../LongestCDS_genetab_eleg.tsv ${file} | awk -F"\t" 'BEGIN{print ">combined"}{printf $2}END{printf "\n"}' > ../filtered/${file%_representative_ids.txt}_filt_largestGeneBody.fasta; done
cd ../filtered/
##Perform R again
R

library(GeneGA)
library(ggplot2)
library(reshape2)


aa=read.fasta("Top500_ubiq_lcapdev_filt_largestGeneBody.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

ubiq_amino=amino
ubiqW=w_value

aa=read.fasta("Top500_soma_lcapdev_filt_largestGeneBody.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

soma_amino=amino
somaW=w_value

aa=read.fasta("Top500_germline_lcapdev_filt_largestGeneBody.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

germ_amino=amino
germW=w_value

aa=read.fasta("Top500_neurons_lcapdev_filt_largestGeneBody.fasta")
rscu=uco(unlist(aa), index="rscu")
w_value=rscu # w_value is the w table we need computing
names(w_value)=names(rscu)
names=sapply(names(rscu), function(x) translate(s2c(x)))
amino=hash()
for(i in unique(names)){
amino[[i]]=length(rscu[which(names==i)])
}
for(i in 1:length(names)){
w_value[i]=rscu[[names(rscu)[i]]]/amino[[translate(s2c(names(rscu)[i]))]]
}

neur_amino=amino
neurW=w_value

data=data.frame(Gemline=germW,Neurons=neurW,Soma=somaW,Ubiq=ubiqW)
data[,5]=rownames(data)
colnames(data)[5]="Amino"

mm=melt(data,id="Amino")
colnames(mm)[2]="Expression"

p=ggplot(mm,aes(Amino,value,colour=Expression)) + geom_point() + ggtitle("Codon usage by expression group") + labs(y="Frequency", x = "Codon")
p
ggsave("Codon_Frequencies_filtered.png",p, width = 40, height = 12, units = "cm")
quit()

cd ..

mkdir intron_distances
cd intron_distances
############Intron analysis###############
###get data
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/gff/c_elegans.PRJNA13758.WS276.annotations.gff3.gz
##Copy refined lists coming from Ahringers paper
cp /home/velazqam/Documents/Projects/CodonOptimization/combined_fasta/*.tsv .




###################################### Dans analysis#######################################
wget http://104.131.81.59/data/sequence_lib_scores.db

#Dump database of 4**12 sequences
db_dump -p sequence_lib_scores.db > sequence_lib_scores.txt

##Get only sequences
tail -n+5 sequence_lib_scores.txt | grep -v "END" | awk '{if((NR%2)==0){print $0}else{printf $0"\t"}}' > Dans_seqs_scores.txt

##Get only sequences with values larger than 0 and remove space at the beggining of seqs
grep -v "0" Dans_seqs_scores.txt | perl -pe 's/^\s+//g' > Dans_seqs_scores_noZeroes.txt 

##Translate them
awk -F"\t" '{print ">"NF"\n"$1}' Dans_seqs_scores_noZeroes.txt > test

#use commandline
transeq -sequence test -outseq test.pep

#Fuse
grep -v ">" test.pep | paste - Dans_seqs_scores_noZeroes.txt > Dans_tab.txt

#Sort and obtain top 5
sort -k3,3n Dans_tab.txt| awk -F"\t" '{count[$1]++; if(count[$1]==1){first[$1]=$0}; if(count[$1]==2){second[$1]=first[$1]; first[$1]=$0}; if(count[$1]==3){third[$1]=second[$1]; second[$1]=first[$1]; first[$1]=$0}; if(count[$1]==4){fourth[$1]=third[$1]; third[$1]=second[$1]; second[$1]=first[$1]; first[$1]=$0}; if(count[$1]>=5){fifth[$1]=fourth[$1]; fourth[$1]=third[$1]; third[$1]=second[$1]; second[$1]=first[$1]; first[$1]=$0};} END{for (keys in count){print first[keys]"\n"second[keys]"\n"third[keys]"\n"fourth[keys]"\n"fifth[keys]}}' - > Top5_Dans.txt

##################################################################



