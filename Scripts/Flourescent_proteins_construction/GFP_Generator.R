########GeneBuilder####
library(ggvis)
library(ggplot2)
library(DT)
library(Biostrings)
library(Cairo)
library(stringdist)


#####Data##########################
##set working directory
setwd("/home/velazqam/Documents/Projects/CodonOptimization/GFPs_iterator")

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

CalculateCAI = function(x,tabibi,cai,list){
  x=toupper(x)
  if(!is.character(x)){return(c())}
  vecseq=unlist(strsplit(x,""))
  if((length(vecseq) %% 3) != 0){return(c())}
  if(cai == 5 ){cai=10}
  CAIvalues=c()
  for(pos in seq(1,length(vecseq),by=3)){
    am=tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]
    tabs=list[[as.character(am)]]
    CAIvalues=append(CAIvalues,tabs[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),cai]/max(tabs[,cai]))
  }
  return(gm_mean(CAIvalues))
}

#############Testing
CodonA=1

##Sequence
GFPseq="atgAGAtccAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTgTTATGGTGTTCAATGCTTcTCgAGATACCCAGATCATATGAAACgGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAAggaTCT"
#Convert to uppercase
GFPseq=toupper(GFPseq)

####Codon Adapt once GFP for Ubiqituos expression
nGFP=repcds(GFPseq,CAIS,AAtoCodF,CodonA)

####Count number of piRNA matches given 5 mismatches
Strcondepies(nGFP,Pies,5)

###Calculate CAI a la Bringmann
CalculateCAI(nGFP,CAIS,5,AAtoCodF)

###Calculate CAI based on codon algorithm
CalculateCAI(nGFP,CAIS,CodonA,AAtoCodF)

###Check consistency
###Calculate CAI based on codon algorithm
CalculateCAI(repcds(GFPseq,CAIS,AAtoCodF,CodonA),CAIS,5,AAtoCodF)

CalculateCAI(repcdsbias(GFPseq,CAIS,AAtoCodF,CodonA,0),CAIS,5,AAtoCodF)

###All good now lets generate sequences
oriseq="CTTTTTACTGGAGTTGTCCCAATTCTTGTCGAACTTGATGGAGATGTCAACGGACACAAATTCTCGGTTTCCGGAGAAGGAGAGGGTGATGCTACTTACGGAAAGCTTACCCTCAAATTTATTTGTACTACGGGTAAGCTCCCTGTACCATGGCCAACCCTCGTCACCACTTTCTGCTATGGAGTTCAATGTTTCTCTCGGTATCCAGATCACATGAAGAGACATGACTTTTTTAAATCGGCCATGCCAGAGGGATACGTACAGGAGAGGACTATCTTTTTTAAGGACGACGGAAACTACAAGACTCGAGCCGAAGTGAAGTTCGAGGGCGATACACTTGTTAACCGGATTGAGCTAAAGGGAATTGATTTCAAGGAGGATGGAAATATTCTGGGACACAAACTTGAATATAATTACAACTCGCACAATGTCTACATTATGGCTGATAAGCAAAAGAACGGAATTAAGGTCAATTTCAAGATTCGACATAACATTGAGGATGGAAGTGTCCAATTGGCTGATCATTATCAACAGAATACTCCAATAGGGGACGGTCCTGTTTTGTTGCCTGATAACCACTATTTGTCCACTCAAAGCGCTTTGTCAAAGGATCCAAACGAAAAACGCGATCATATGGTTCTCCTTGAGTTCGTAACGGCCGCTGGAATCACTCACGGTATGGATGAACTCTATAAATGA"
seqs=c()
#10 per optimization per bias
for(al in 1:5){
for(i in c(0:10)/10){
for(j in 1:10){
  nnseq=repcdsbias(oriseq,CAIS,AAtoCodF,al,i)
  seqs=append(seqs,nnseq)
}
}}

caiss =c()

for(seq in seqs){
caiss=rbind(caiss,cbind(CalculateCAI(seq,CAIS,1,AAtoCodF),CalculateCAI(seq,CAIS,2,AAtoCodF),CalculateCAI(seq,CAIS,3,AAtoCodF),CalculateCAI(seq,CAIS,4,AAtoCodF),CalculateCAI(seq,CAIS,5,AAtoCodF)))  
}

nbpi = c()

for(seq in seqs){
nbpi = append(nbpi,length(Strfindpies(seq,Pies,3)))
  }

data=cbind(caiss,nbpi)
##Left there for further analysis
#load("~/Documents/Projects/CodonOptimization/GFPs_iterator/.RData")


#writeLines(seqs,"GFPs.fasta")

#write.table(data, "Values.tsv",sep="\t", row.names=F, col.names=F)

###Lets do it  once again but with higher number

###All good now lets generate sequences
oriseq="CTTTTTACTGGAGTTGTCCCAATTCTTGTCGAACTTGATGGAGATGTCAACGGACACAAATTCTCGGTTTCCGGAGAAGGAGAGGGTGATGCTACTTACGGAAAGCTTACCCTCAAATTTATTTGTACTACGGGTAAGCTCCCTGTACCATGGCCAACCCTCGTCACCACTTTCTGCTATGGAGTTCAATGTTTCTCTCGGTATCCAGATCACATGAAGAGACATGACTTTTTTAAATCGGCCATGCCAGAGGGATACGTACAGGAGAGGACTATCTTTTTTAAGGACGACGGAAACTACAAGACTCGAGCCGAAGTGAAGTTCGAGGGCGATACACTTGTTAACCGGATTGAGCTAAAGGGAATTGATTTCAAGGAGGATGGAAATATTCTGGGACACAAACTTGAATATAATTACAACTCGCACAATGTCTACATTATGGCTGATAAGCAAAAGAACGGAATTAAGGTCAATTTCAAGATTCGACATAACATTGAGGATGGAAGTGTCCAATTGGCTGATCATTATCAACAGAATACTCCAATAGGGGACGGTCCTGTTTTGTTGCCTGATAACCACTATTTGTCCACTCAAAGCGCTTTGTCAAAGGATCCAAACGAAAAACGCGATCATATGGTTCTCCTTGAGTTCGTAACGGCCGCTGGAATCACTCACGGTATGGATGAACTCTATAAATCTCGGCGGCGTAAGGCCAACCCTACTAAGTTGAGCGAGAATGCAAAAAAGCTTGCGAAGGAAGTGGAGAATTGA"
seqs2=c()
#100 per optimization per bias; it seems it breaks around bias equal to 1... Not enough replacement?
for(al in 1:5){
  for(i in c(0:10)/10){
    for(j in 1:100){
      nnseq=repcdsbias(oriseq,CAIS,AAtoCodF,al,i)
      seqs2=append(seqs2,nnseq)
    }
  }}

rema="ATGCCCAAAAAGAAACGGAAAGTTTCTAAGGGCGAAGAA"
seqs2=paste(rema,seqs2,sep="")

caiss2 <- data.frame(ca1 = rep(NA,length(seqs2)),ca2 = rep(NA,length(seqs2)),ca3 = rep(NA,length(seqs2)),ca4 = rep(NA,length(seqs2)),ca5 = rep(NA,length(seqs2)))


for(i in 1:length(seqs2)){
  caiss2[i,]=c(CalculateCAI(seqs2[i],CAIS,1,AAtoCodF),CalculateCAI(seqs2[i],CAIS,2,AAtoCodF),CalculateCAI(seqs2[i],CAIS,3,AAtoCodF),CalculateCAI(seqs2[i],CAIS,4,AAtoCodF),CalculateCAI(seqs2[i],CAIS,5,AAtoCodF)) 
}


nbpi2 <- data.frame(pi3 = rep(NA,length(seqs2)),pi4 = rep(NA,length(seqs2)),pi5 = rep(NA,length(seqs2)))

for(i in 1:length(seqs2)){
  nbpi2[i,]=c(length(unique(Strfindpies(seqs2[i],Pies,3))), length(unique(Strfindpies(seqs2[i],Pies,4))), length(unique(Strfindpies(seqs2[i],Pies,5)))) 
}



data2=cbind(caiss2,nbpi2)

writeLines(seqs2,"GFPs_3.fasta")

write.table(data2, "Values_3.tsv",sep="\t", row.names=F, col.names=F)


