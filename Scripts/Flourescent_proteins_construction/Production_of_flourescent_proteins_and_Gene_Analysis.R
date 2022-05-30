##At all vals directory
setwd("/home/velazqam/Documents/Projects/CodonOptimization/Analysis/All_vals")

data=read.table(paste("Top500_","ubiq","_lcapdev.tsv",sep=""),sep="\t",header=F, stringsAsFactors= F)
data[,3]=rep("Ubiquitous",nrow(data))


files=c("germline","soma", "neurons")


for(file in files){
  temp=read.table(paste("Top500_",file,"_lcapdev.tsv",sep=""),sep="\t",header=F, stringsAsFactors= F)
  temp[,3]=rep(file,nrow(temp))
  data=rbind(data,temp)
}

colnames(data)=c("WBID","Sequence","Tissue")

###Load functions / data
##Codon Adapt
CAIS=read.table("DATA/Ultimate_aminos2.txt",sep="\t",header=T)
rownames(CAIS)= toupper(as.character(CAIS$Codon))
codons=unique(CAIS$Amino)

AAtoCodF=list()

for(i in 1:length(codons)){
  AAtoCodF=append(AAtoCodF, list(CAIS[which(CAIS$Amino == codons[i]),]))
  names(AAtoCodF)[i]=as.character(codons[i])
}

IntronSeqs=read.table("DATA/Introns.csv",sep=",",header=F,row.names=1)

##piRNAs

Pies=readLines("DATA/HengPies.txt")

PiesNA=readLines("DATA/HengNames.txt")

PiesFin=cbind(Pies,PiesNA)
rownames(PiesFin)=as.character(Pies)

##Enzymes
enzy=read.table("DATA/Enzymes.txt", sep="\t", colClasses = "character")
colnames(enzy)=c("Enzyme", "Site")
###TO have proper limited columns. In case that is visually aberrant, remove this as well as in server.R
enzy$Enzyme=gsub("\\s", ".", format(enzy$Enzyme, width=max(nchar(enzy$Enzyme))))
rownames(enzy)=as.character(enzy$Enzyme)

####Functions#######
##Codon
sampcod=function(aa,list,cai){
  newcod=sample((list[[aa]])[,6],1,prob=(list[[aa]])[,cai])
  return(toupper(as.character(newcod)))
}

repcds=function(x,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),list,cai)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

sampnewcod=function(aa,oldcodon,list,cai){
  oldcodon=toupper(oldcodon)
  if(nrow(list[[aa]]) < 2 ){return(oldcodon)}
  oldcodon =which(as.character(rownames(list[[aa]])) == oldcodon)
  newcod=sample((list[[aa]])[-c(oldcodon),6],1,prob=(list[[aa]])[-c(oldcodon),cai])
  return(toupper(as.character(newcod)))
}

repnewcds=function(x,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[i:(i+2)]),sep="",collapse=""),list,cai)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

modbyposiz=function(x,starts,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  if(length(starts)>0){starts=unique(starts - (starts %% 3) +1 )}
  for(pos in c(starts)){
    nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),list,cai)
    vecseq[pos:(pos+2)]=unlist(strsplit(nncod,""))
  }
  return(paste(vecseq,sep="",collapse=""))
}

###Create random probs for GFP sequence
rfrac= function(n){
  val=runif(n)
  return(val/sum(val))
}

##Codon
sampcodtwo=function(aa,list,cai){
  newcod=sample((list[[aa]])[,6],1,prob=rfrac(length((list[[aa]])[,cai])))
  return(toupper(as.character(newcod)))
}

repcdstwo=function(x,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampcodtwo(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),list,cai)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

sampnewcodtwo=function(aa,oldcodon,list,cai){
  oldcodon=toupper(oldcodon)
  if(nrow(list[[aa]]) < 2 ){return(oldcodon)}
  oldcodon =which(as.character(rownames(list[[aa]])) == oldcodon)
  newcod=sample((list[[aa]])[-c(oldcodon),6],1,prob=rfrac(length((list[[aa]])[-c(oldcodon),cai])))
  return(toupper(as.character(newcod)))
}

repnewcdstwo=function(x,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampnewcodtwo(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[i:(i+2)]),sep="",collapse=""),list,cai)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

modbyposiztwo=function(x,starts,tabibi,list,cai){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  if(length(starts)>0){starts=unique(starts - (starts %% 3) +1 )}
  for(pos in c(starts)){
    nncod=sampnewcodtwo(as.character(tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),list,cai)
    vecseq[pos:(pos+2)]=unlist(strsplit(nncod,""))
  }
  return(paste(vecseq,sep="",collapse=""))
}

###Do by bias
###Create bias probs for CAI-GFP sequence

##Codon
sampcodbias=function(aa,list,cai,bias){
  pps=as.numeric((list[[aa]])[,cai])
  if((bias > runif(1))&(length(pps)>1)){
    pps=1-pps
    pps=pps/sum(pps)
  }
  newcod=sample((list[[aa]])[,6],1,prob=pps)
  return(toupper(as.character(newcod)))
}

repcdsbias=function(x,tabibi,list,cai,bias){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampcodbias(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),list,cai,bias)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

sampnewcodbias=function(aa,oldcodon,list,cai,bias){
  oldcodon=toupper(oldcodon)
  if(nrow(list[[aa]]) < 2 ){return(oldcodon)}
  oldcodon =which(as.character(rownames(list[[aa]])) == oldcodon)
  pps=as.numeric((list[[aa]])[-c(oldcodon),cai])
  if((bias > runif(1))&(length(pps)>1)){
    pps=1-pps
    pps=pps/sum(pps)
  }
  newcod=sample((list[[aa]])[-c(oldcodon),6],1,prob=pps)
  
  return(toupper(as.character(newcod)))
}

repnewcdsbias=function(x,tabibi,list,cai,bias){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  nnseq=c()
  for(i in seq(1,length(vecseq),by=3)){
    nncod=sampnewcodbias(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[i:(i+2)]),sep="",collapse=""),list,cai,bias)
    nnseq=append(nnseq,nncod)
  }
  return(paste(nnseq,sep="",collapse=""))
}

modbyposizbias=function(x,starts,tabibi,list,cai,bias){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  if(length(starts)>0){starts=unique(starts - (starts %% 3) +1 )}
  for(pos in c(starts)){
    nncod=sampnewcodbias(as.character(tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),list,cai,bias)
    vecseq[pos:(pos+2)]=unlist(strsplit(nncod,""))
  }
  return(paste(vecseq,sep="",collapse=""))
}


##Pis
countpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 21){return(c())}
  return((stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")))
}

Strcountpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 21){return(c())}
  con=0
  for(i in 1:(length(vecseq)-20)){
    con=con+countpies(paste(vecseq[i:(i+20)],collapse=""),y)
  }
  return(con)
}

countmatpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 21){return(c())}
  return(length(which(stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")<6)))
}

Strcountmatpies=function(x,y){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 21){return(c())}
  con=0
  for(i in 1:(length(vecseq)-20)){
    con=con+countmatpies(paste(vecseq[i:(i+20)],collapse=""),y)
  }
  return(con)
}

condepies=function(x,pies,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 20){return(c())}
  return(sum(stringdist(x,pies,method="hamming") <= mm))
}

Strcondepies=function(x,y,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 20){return(c())}
  con=0
  for(i in 1:(length(vecseq)-19)){
    con=con+condepies(paste(vecseq[i:(i+19)],collapse=""),y,mm)
  }
  return(con)
}

findpies=function(x,pies,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) != 20){return(c())}
  idx=c(stringdist(x,pies,method="hamming") <= mm)
  if(sum(idx)>0){return(pies[idx])}else{return()}
}

Strfindpies=function(x,y,mm){
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if(length(vecseq) < 20){return(c())}
  con=c()
  for(i in 1:(length(vecseq)-19)){
    con=append(con,findpies(paste(vecseq[i:(i+19)],collapse=""),y,mm))
  }
  return(con)
}

###HTML Gene info
CalculateGC = function(x){
  if(!is.character(x)){return(c())}
  x=toupper(x)
  vecseq=unlist(strsplit(x,""))
  return((countPattern("C",x)+countPattern("G",x))/length(vecseq))
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


###Replicate Bringmans values
CalculateCAI = function(x,tabibi,cai,list){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  if(cai == 5 ){cai=11}
  CAIvalues=c()
  for(pos in seq(1,length(vecseq),by=3)){
    am=tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]
    tabs=list[[as.character(am)]]
    CAIvalues=append(CAIvalues,tabs[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),cai]/max(tabs[,cai]))
  }
  return(gm_mean(CAIvalues))
}

###Re-analyze everything
##Redo tables
nSeqs=data$Sequence

caiss <- data.frame(ca1 = rep(NA,length(nSeqs)),ca2 = rep(NA,length(nSeqs)),ca3 = rep(NA,length(nSeqs)),ca4 = rep(NA,length(nSeqs)),ca5 = rep(NA,length(nSeqs)))

for(i in 1:length(nSeqs)){
  caiss[i,]=c(CalculateCAI(nSeqs[i],CAIS,1,AAtoCodF),CalculateCAI(nSeqs[i],CAIS,2,AAtoCodF),CalculateCAI(nSeqs[i],CAIS,3,AAtoCodF),CalculateCAI(nSeqs[i],CAIS,4,AAtoCodF),CalculateCAI(nSeqs[i],CAIS,5,AAtoCodF)) 
}


nbpi <- data.frame(pi3 = rep(NA,length(nSeqs)),pi4 = rep(NA,length(nSeqs)),pi5 = rep(NA,length(nSeqs)))

for(i in 1:length(nSeqs)){
  nbpi[i,]=c(length(unique(Strfindpies(nSeqs[i],Pies,3))), length(unique(Strfindpies(nSeqs[i],Pies,4))), length(unique(Strfindpies(nSeqs[i],Pies,5)))) 
}

Values=cbind(caiss,nbpi)

writeLines(nSeqs,"TopSets_4GLO.fasta")

write.table(Values, "Values_TopSets.tsv",sep="\t", row.names=F, col.names=F)

##System call
###perl Score_Dans.pl TopSets_4GLO.fasta | paste - Values_TopSets.tsv > TopSets_vals.tsv

### Read new values
Tops = read.table("TopSets_vals.tsv", sep="\t", header=F,stringsAsFactors=F)

colnames(Tops) =c("Sequence","GLO_Score","GLO_LowestScore","GLO_LowestWordScore","CAI_Ubiquitous","CAI_Germline","CAI_Neurons","CAI_Somatic","CAI_1","TargetPies_3MM","TargetPies_4MM","TargetPies_5MM")

##Do ramp
write.table(data,"TopSets_seqs.tsv",sep="\t",col.names=F,row.names=F,quote=F)


TopSets=cbind(data,Tops)

load("minimal-data.RData")

head(lcap.dt)

Hepipi=read.table("piRNA-scan_all_transcripts-data.tsv",sep="\t", header=T, stringsAsFactors=F)

GeNames=read.table("LargestGB-names-relation.tsv",sep="\t",header=F,stringsAsFactors=F)

rownames(GeNames)=as.character(GeNames[,2])
rownames(Hepipi) = as.character(Hepipi[,1])

Hepis=cbind(GeNames,Hepipi[rownames(GeNames),])

#####Octuber 15 graphs
library(ggplot2)
library(plyr)

##Names of graphs are with lower score; fix that when plotting real figures
##C.A.I.
df=data.frame(Tissue=TopSets$Tissue, Value=TopSets$CAI_1)

p<-ggplot(df, aes(x=Value))+
  geom_histogram(color="black", fill="white", breaks=seq(0,1,.1))+
  scale_x_continuous(breaks = seq(0,1,.1), labels = seq(0,1,.1)) +
  facet_grid(Tissue ~ .) + xlab("C.A.I.") 
p

##piRNA 3MM matches
df=data.frame(Tissue=TopSets$Tissue, Value=TopSets$TargetPies_3MM)

p<-ggplot(df, aes(x=Value))+
  geom_histogram(color="black", fill="white", breaks=seq(0,max(df$Value)+5,5))+
  scale_x_continuous(breaks = seq(0,max(df$Value)+5,5), labels = seq(0,max(df$Value)+5,5)) +
  facet_grid(Tissue ~ .) + xlab("piRNA 3MM") 

p

##piRNA 4MM matches
df=data.frame(Tissue=TopSets$Tissue, Value=TopSets$TargetPies_4MM)

p<-ggplot(df, aes(x=Value))+
  geom_histogram(color="black", fill="white", breaks=seq(0,max(df$Value)+5,5))+
  scale_x_continuous(breaks = seq(0,max(df$Value)+5,25), labels = seq(0,max(df$Value)+5,25)) +
  facet_grid(Tissue ~ .) + xlab("piRNA 4MM") 

p

##piRNA 5MM matches
df=data.frame(Tissue=TopSets$Tissue, Value=TopSets$TargetPies_5MM)

p<-ggplot(df, aes(x=Value))+
  geom_histogram(color="black", fill="white", breaks=seq(0,max(df$Value)+5,25))+
  scale_x_continuous(breaks = seq(0,max(df$Value)+5,125), labels = seq(0,max(df$Value)+5,125)) +
  facet_grid(Tissue ~ .) + xlab("piRNA 5MM") 

p

##piRNA sites in heng chis site
#Note that HengCHis table does not seem to have data for all the transcripts analyzed previouslyand so numbers differ
testo=Hepis[which(!(is.na(Hepis$No..Predidcted.piRNA.target.sites))),]
rownames(data)=as.character(data[,1])
table(data[as.character(testo[,1]),"Tissue"])


testo=cbind(testo,data[as.character(testo[,1]),"Tissue"])
colnames(testo)[ncol(testo)]="Tissue"

df=data.frame(Tissue=testo$Tissue, Value=testo$No..Predidcted.piRNA.target.sites)

p<-ggplot(df, aes(x=Value))+
  geom_histogram(color="black", fill="white", breaks=seq(0,max(df$Value)+5,5))+
  scale_x_continuous(breaks = seq(0,max(df$Value)+5,25), labels = seq(0,max(df$Value)+5,25)) +
  facet_grid(Tissue ~ .) + xlab("Heng Chis predicted piRNAs sites") 

p

df=data.frame(Tissue=testo$Tissue, Value=testo$No..CLASH.Identified.piRNA.target.sites)

p<-ggplot(df, aes(x=Value))+
  geom_histogram(color="black", fill="white", breaks=seq(0,max(df$Value)+5,5))+
  scale_x_continuous(breaks = seq(0,max(df$Value)+5,25), labels = seq(0,max(df$Value)+5,25)) +
  facet_grid(Tissue ~ .) + xlab("Heng Chis identified CLASH piRNAs") 

p

## GLO score
df=data.frame(Tissue=TopSets$Tissue, Value=TopSets$GLO_Score)

p<-ggplot(df, aes(x=Value))+
  geom_histogram(color="black", fill="white", breaks=seq(0,max(df$Value)+5,500))+
  scale_x_continuous(breaks = seq(0,max(df$Value)+50,1000), labels = seq(0,max(df$Value)+50,1000)) +
  facet_grid(Tissue ~ .) + xlab("GLO score") 

p

#We used an RNAseq dataset from one-cell embryos (Gerstein et al., 2014; Hillier et al., 2009) as a proxy for germline-expressed mRNAs, because one-cell embryos are not transcriptionally active and thus their entire mRNA content should be derived from the maternal germline.
## GLO score is not necerarilly germline dependent but Ubiquitous aswell the two highest genes are eff-1 and one atp-ase

######################RAMP###
###All but is better if group by tissue expression
##Function to generate SD values
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

topo=read.table("TopSets_seqs.ramp",sep=",",header=T, stringsAsFactors=F)
topo=cbind(topo,data[gsub(pattern=">",replacement="",x=as.character(topo$seq)),"Tissue"])

colnames(topo)[ncol(topo)]="Tissue"

df= data.frame(topo)

df2 <- data_summary(df,varname="speed_value",groupnames=c("position", "Tissue"))

ggplot(df2, aes(x=position, y=speed_value)) +  geom_line() + geom_point() + geom_errorbar(aes(ymin=speed_value-sd, ymax=speed_value+sd), width=.2, position=position_dodge(0.05)) +
  facet_grid(Tissue ~ .) 

df3=head(df2,200)

ggplot(df3, aes(x=position, y=speed_value)) +  geom_line() + geom_point() + geom_errorbar(aes(ymin=speed_value-sd, ymax=speed_value+sd), width=.2, position=position_dodge(0.05)) +
  facet_grid(Tissue ~ .) 

dev.off()

#############################################################################
###Save intermediary working environment
save.image(file = "WorkingSpaceImages_Ranalysis_on_Genes.RData")

###################################################
###Free energy










