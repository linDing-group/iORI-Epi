
# setwd("iORI-Epi")

library(pROC)
library(glmnet)
library(ranger)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg19)
# library(ROCR)

# Read bed file containing ChIP-seq peak data
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


# Map a genomic feature to loci. 
annotateLoci <- function(loci.GR, GenomicFeature.GR)
{
 if(!is(loci.GR, "GenomicRanges")){stop("'loci.GR' must be a GenomicRanges object")}
 if(!is(GenomicFeature.GR, "GenomicRanges")){stop("'GenomicFeature.GR' must be a GenomicRanges object")}
 if(any(is.na(seqlengths(GenomicFeature.GR)))){stop("'seqlengths(GenomicFeature.GR)' contains NAs")}

 if(is.null(GenomicFeature.GR$score)){
  cvg <- coverage(GenomicFeature.GR)
 }else{
  cvg <- coverage(GenomicFeature.GR, weight="score")
 }

 chr=seqnames(seqinfo(loci.GR))
 loci.list=list()
 for(i in 1:length(chr)){
  loci.list[[chr[i]]]=loci.GR[seqnames(loci.GR)==chr[i]]
 }
 loci.GRL=GenomicRangesList(loci.list)

 if(length(chr)==1){
  averageBin=viewMeans(Views(cvg[[chr]],ranges(loci.GRL[[chr]])))
 }else{
  views_list <- RleViewsList(lapply(names(cvg),function(seqname)Views(cvg[[seqname]], ranges(loci.GRL[[seqname]]))))
  averageBin=NULL
  for(i in 1:length(views_list)){
   averageBin=c(averageBin,viewMeans(views_list)[[i]])
  }
 }

 return(averageBin)
}



##########
Genome=BSgenome.Hsapiens.UCSC.hg19
model="human"
SeqinfoGenome=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=c(paste0("chr",1:22),"chrX")
SeqInfo=SeqinfoGenome[Chr.V]


fileAnnot=list.files("../data/test_loop")
AnnotNames=as.vector(sapply(fileAnnot,function(x){strsplit(x,'_')[[1]][1]}))
print(length(AnnotNames))

# Import ORIs dataset
fileBedOriPos="testpos.bed"
dataOriPos.GR=sort(readGFBed(fileBedOriPos,SeqInfo))
fileBedOriNeg="testneg.bed"
dataOriNeg.GR=sort(readGFBed(fileBedOriNeg,SeqInfo))
dataOris.GR=c(dataOriPos.GR,dataOriNeg.GR)

# Import Chip-seq data
GenomicFeatureList.GR=list()
for(i in 1:length(AnnotNames)){
    GenomicFeatureList.GR[[i]] <- sort(unique(readGFBed(paste0("../data/test_loop/",fileAnnot[i]),SeqInfo)))
    print(paste0(AnnotNames[i]," : ",length(GenomicFeatureList.GR[[i]])))
}
names(GenomicFeatureList.GR)=AnnotNames
GenomicFeatureList.GR=GenomicFeatureList.GR[names(GenomicFeatureList.GR)%in%c("chiapetpol2")]
AnnotNames=names(GenomicFeatureList.GR)


# Map features to ORI and non-ORI regions 
bin.Mat=c(rep(1,length(dataOriPos.GR)),rep(0,length(dataOriNeg.GR)))

for(i in 1:length(GenomicFeatureList.GR)){
  GRi=GenomicFeatureList.GR[[i]]
  annotPosi=annotateLoci(dataOriPos.GR,GRi)
  annotNegi=annotateLoci(dataOriNeg.GR,GRi)
  annoti=c(annotPosi,annotNegi)
  annoti[annoti>1]=1
  bin.Mat=cbind(bin.Mat,annoti)
  rm(annoti)
  print(paste0(AnnotNames[i]," annotated"))
 }
colnames(bin.Mat)=c("class",AnnotNames)

# Features save as csv file
filename = "test_3D"
file_varimp= paste(filename,"_Fea1.csv",sep='')
write.table(bin.Mat,file=file_varimp,row.names=F,sep=',',quote=F)


