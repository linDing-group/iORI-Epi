setwd("iORI-Epi")

library(pROC)
library(glmnet)
library(ranger)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(JASPAR2018)
library(TFBSTools)


readGFBed<- function(GFBedFile,seqInfoChr){

 if(class(seqInfoChr)!="Seqinfo"){print("seqInfoChr is not a seqinfo object!");return(0)}
 if(length(grep(pattern='.bed',x=GFBedFile))==0){print("GFBedFile is not a bed file!");return(0)}

 dataGF=read.table(GFBedFile,sep='\t',header=F)[,c(1:3)]
 dataGF2=dataGF[dataGF[,1]%in%seqnames(seqInfoChr),]
 chr=as.character(dataGF2[,1])
 posLeft=as.numeric(dataGF2[,2])
 posRight=as.numeric(dataGF2[,3])

 GenomicFeature.GR=NULL
 for(i in 1:length(seqnames(seqInfoChr))){
  chri=seqnames(seqInfoChr)[i]
  if(sum(chr==chri)){
   GenomicFeature.GRi <- GRanges(seqnames=chri,IRanges(start=posLeft[chr==chri],end=posRight[chr==chri]))
   if(is.null(GenomicFeature.GR)){
    GenomicFeature.GR=GenomicFeature.GRi
   }else{
    GenomicFeature.GR=c(GenomicFeature.GR,GenomicFeature.GRi)
   }
  }
 }

 seqlevels(GenomicFeature.GR) = seqlevels(seqInfoChr)
 seqinfo(GenomicFeature.GR)=seqInfoChr

 return(GenomicFeature.GR)
}

# Function to find protein binding sites predicted from DNA motifs in a sequence

findMotifs<-function(seq.GR,PFML,genome,motifname){
  
  seq=as.list(getSeq(genome, seqnames(seq.GR), start=start(seq.GR), end=end(seq.GR), as.character=TRUE))
  matMotifCount=NULL
  name_id=NULL

  for(i in motifname){ 
    
    print(i)
    PFMi = PFML@listData[[i]]
    PWMi=as.matrix(toPWM(PFMi))
    #PWMci=reverseComplement(PWMi)
    idi=ID(PFMi)
    namei=name(PFMi)
    name_id=c(name_id,paste0(namei,"_",idi))
    
    posMotifPos=unlist(lapply(lapply(seq,function(x){matchPWM(PWMi,x)}),length))
    #posMotifNeg=unlist(lapply(lapply(seq,function(x){matchPWM(PWMci,x)}),length))
    #posMotifAll=c(posMotifPos,posMotifNeg)
    posMotifAll=posMotifPos
    
    matMotifCount=cbind(matMotifCount,posMotifAll)
    
    print(paste0(namei," annotated"))
  }
  colnames(matMotifCount)=name_id
  
  return(matMotifCount)
}



#Extracting motif features

Genome=BSgenome.Hsapiens.UCSC.hg19
SeqinfoGenome=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=c(paste0("chr",1:22),"chrX")
SeqInfo=SeqinfoGenome[Chr.V]


i = 'test'
fileBedOriPos = paste(i,'pos.bed',sep="")
dataOriPos.GR=sort(readGFBed(fileBedOriPos,SeqInfo))
fileBedOriNeg=paste(i,'neg.bed',sep="")
dataOriNeg.GR=sort(readGFBed(fileBedOriNeg,SeqInfo))
dataOris.GR=c(dataOriPos.GR,dataOriNeg.GR)
  
  
opts <- list()
opts[["species"]] <- 9606 # 9609 represent human
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
  

  
bin.Mat=c(rep(1,length(dataOriPos.GR)),rep(0,length(dataOriNeg.GR)))



MotifName=c('MA0002.1','MA0028.1','MA0025.1')#More motifs can be included, here we only use three motifs as an example

matMotifs=findMotifs(dataOris.GR,PFMatrixList,BSgenome.Hsapiens.UCSC.hg19.masked,MotifName)
  
#save(matMotifs,file=paste(i,'_promoter.RData',sep=""))
#load(file="matMotifs.RData")
bin.Mat=cbind(bin.Mat,matMotifs)
rownames(bin.Mat)=1:nrow(bin.Mat)
colnames(bin.Mat)[1]="class"
  
  
# Features save as csv file
file_varimp=paste(i,'_Motif_Fea.csv',sep="")
write.table(bin.Mat,file=file_varimp,row.names=F,sep=',',quote=F)




