#Wget 276
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS276/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS276.annotations.gff3.gz

zcat c_elegans.PRJNA13758.WS276.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=" '{print $1"\t"$4"\t"$5"\t"$12"\t.\t"$7}' | perl -pe 's/Transcript://g' | awk -F"\t" '{OFS="\t"; split($4,trans,","); for(i=1;i<=length(trans);i++){print $1,$2,$3,trans[i],$5,$6}}' | awk -F"\t" '{if(array[$4] != 0){if(start[$4] > $2){start[$4]=$2};if(end[$4] < $3){end[$4]=$3};array[$4]=$1"\t"start[$4]"\t"end[$4]"\t"$4"\t"$5"\t"$6}else{start[$4]=$2;end[$4]=$3;array[$4]=$0}}END{for(key in array){print array[key]}}' - | sort > Cel_GeneBodies.bed

zcat c_elegans.PRJNA13758.WS276.annotations.gff3.gz | grep -v "#" | grep "WormBase" | awk -F"\t" '{if ($3=="mRNA"){print $9}}'| awk -F";" '{split($1,ts,":"); split($2,gs,":"); split($5,loc,"="); print ts[2]"\t"gs[2]";"ts[2]";"loc[2]}' | awk -F"\t" '{OFS="\t"; print $1,$2,$2,$1}' > Cel_mRNAID.txt

awk -F"\t" '{OFS="\t"; if(array[$4]==0){array[$4]=$2}else{print $1,($2-1),$3,array[$4],$5,$6}}' Cel_mRNAID.txt Cel_GeneBodies.bed > Cel_Genes.bed

awk -F"\t" '{OFS="\t"; split($4,info,";"); wbid=info[1]; if(data[wbid] ==""){data[wbid]=$0; siz[wbid]=($3-$2);}else{if(($3-$2) > siz[wbid]){data[wbid]=$0; siz[wbid]=($3-$2);}}} END{for(key in data){print data[key]}}' Cel_Genes.bed | sort -k1,1 -k2,2n > Cel_largest_Genebodies.bed



zcat c_elegans.PRJNA13758.WS278.annotations.gff3.gz | awk -F"\t" '{OFS="\t"; if($2=="WormBase"){if($3=="intron"){print $1,$4,$5,$9,".",$7}}}' | awk -F"\t" '{array[$1"-"$2"-"$3"-"$6]=$0} END{for(key in array){print array[key]}}' | sort -k1,1 -k2,2n | awk -F"\t" '{OFS="\t"; split($4,info,";"); split(info[1],data,":"); print $1,$2,$3,data[2],$5,$6}' > Cel_Introns.bed

zcat c_elegans.PRJNA13758.WS276.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=" '{print $1"\t"$4"\t"$5"\t"$12"\t.\t"$7}' | perl -pe 's/Transcript://g' | awk -F"\t" '{OFS="\t"; split($4,trans,","); for(i=1;i<=length(trans);i++){print $1,$2,$3,trans[i],$5,$6}}' | awk -F"\t" '{if(array[$4] != 0){if(strand[$4] == "-"){if(start[$4]<$2){start[$4]=$2;end[$4]=$3;array[$4]=$0;}}else{if(start[$4]>$2){start[$4]=$2;end[$4]=$3;array[$4]=$0;}}}else{start[$4]=$2;end[$4]=$3;array[$4]=$0; strand[$4]=$6}}END{for(key in array){print array[key]}}' - | sort -k1,1 -k2,2n  > Cel_firstCDS.bed

awk -F"\t" '{OFS="\t"; if(array[$4]==0){array[$4]=$2}else{print $1,($2-1),$3,array[$4],$5,$6}}' Cel_mRNAID.txt Cel_firstCDS.bed > Cel_1stCDS.bed


for file in `ls Top*.list`; do echo ${file%.list}; awk -F"\t" '{OFS="\t"; split($4,info,";"); print info[1],$0}' Cel_largest_Genebodies.bed | awk -F"\t" '{if(data[$1]==""){data[$1]=$0}else{print data[$1]}}' - ${file} | awk -F"\t" '{if(NF == 7){print data[$5]}else{data[$4]=$0}}' Cel_1stCDS.bed - > ${file%.list}_firstCDSinlargestgb.bed; done

##R
R

files=c("ubiq","germline","soma", "neurons")
dists=list()

for(file in files){
	temp=read.table(paste("Top500_",file,"_lcapdev_firstCDSinlargestgb.bed",sep=""),sep="\t",header=F)
	nn=temp[,3]-temp[,2]
	dists=c(dists,list(nn))
	names(dists)[length(dists)]=file
}

df = data.frame(Tissue=c(rep("Ubiquitous",length(dists$ubiq)),rep("Germline",length(dists$germline)),rep("Somatic",length(dists$soma)),rep("Neuronal",length(dists$neurons))), Distance=c(dists$ubiq,dists$germline,dists$soma,dists$neurons))

library(ggplot2)
library(plyr)
mu <- ddply(df, "Tissue", summarise, grp.mean=mean(Distance))

 p<-ggplot(df, aes(x=Distance))+
  geom_histogram(color="black", fill="white", breaks=seq(0,max(df$Distance)+50,50))+
  scale_x_continuous(breaks = seq(0,max(df$Distance)+50,50), labels = seq(0,max(df$Distance)+50,50)) +
  facet_grid(Tissue ~ .)
p
# Add mean lines
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")



