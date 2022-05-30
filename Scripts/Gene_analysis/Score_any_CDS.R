####Libraries you need to install
#install.packages("stringdist")
#install.packages("Biostrings")
library(stringdist)
library(Biostrings)


#################### Internals #############################
###Load Amino tables and PiRNAseqs and format them for analysis
#Codon Adapt
        CAIS=read.table("Ultimate_aminos.txt",sep="\t",header=T)
        rownames(CAIS)= toupper(as.character(CAIS$Codon))
        codons=unique(CAIS$Amino)
        
        AAtoCodF=list()
        
        for(i in 1:length(codons)){
          AAtoCodF=append(AAtoCodF, list(CAIS[which(CAIS$Amino == codons[i]),]))
          names(AAtoCodF)[i]=as.character(codons[i])
        }
        
        ##piRNA
        pi2UGs=read.table("Ultimate_seed-remainder_iterations_2GUs.txt",sep="\t",header=F)
        
        Piseqs=list()
        for(i in 1:nrow(pi2UGs)){
          Piseqs=append(Piseqs,strsplit(as.character(pi2UGs[i,2]),";"))
          names(Piseqs)[i]=as.character(pi2UGs[i,1])
        }

        ##New approach
        Pies=readLines("Pies20.txt")
###In house functions
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

########Initial parameters###
##Codon algorithm: 1 for mimic Ubiq freqs, 2 for Germ freqs, 3 for neuron, 4 soma, and 5 for CAI equal to 1
CodonA=1

##Sequence
GFPseq="atgAGAtccAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTgTTATGGTGTTCAATGCTTcTCgAGATACCCAGATCATATGAAACgGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAAggaTCT"
#Convert to uppercase
GFPseq=toupper(GFPseq)

########EXAMPLES###########

####Codon Adapt once GFP for Ubiqituos expression
nGFP=repcds(GFPseq,CAIS,AAtoCodF,CodonA)

####Count number of piRNA matches given 5 mismatches
Strcondepies(nGFP,Pies,5)

############Iterations#########
###Given GFP sequence, iterate to generate 10, 100, 1000, 10K seqs, obtain their different matches with piRNA and save the time it took

##Storage variables
time=c()
data=list()

##For i having the different numbers of seqs
for(i in c(10,100,1000,10000)){

##Initialize table that contains results per group of sequences
dattab=data.frame(matrix(ncol = 5, nrow = i))

##Encapsulate everything into a system time function, from which only the elapsed time will be saved
ntime=system.time({

###Produce seqs iteratively
seqs=c()
for(j in 1:i){
nnseq=repcds(GFPseq,CAIS,AAtoCodF,1)
seqs=append(seqs,nnseq)
}

##From sets of sequences calculate number of piRNA matches given mismatch
for(j in 1:length(seqs)){
dattab[j,1]=Strcondepies(seqs[j],Pies,1)
dattab[j,2]=Strcondepies(seqs[j],Pies,2)
dattab[j,3]=Strcondepies(seqs[j],Pies,3)
dattab[j,4]=Strcondepies(seqs[j],Pies,4)
dattab[j,5]=Strcondepies(seqs[j],Pies,5)
}
})[3]

##Apend matrix generated
data=append(data,list(dattab))
names(data)=paste(c("Iter-",i),sep="",collapse="")

#Save execution time
time=append(time,ntime)
}

###Done
data[1]

