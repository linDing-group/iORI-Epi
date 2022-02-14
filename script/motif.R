setwd("iORI-Epi")

library(pROC)
library(glmnet)
library(ranger)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(JASPAR2018)
library(TFBSTools)


# Function to find protein binding sites predicted from DNA motifs in a sequence
findMotifs<-function(seq.GR,PFML,genome){

 seq=as.list(getSeq(genome, seqnames(seq.GR), start=start(seq.GR), end=end(seq.GR), as.character=TRUE))
 matMotifCount=NULL
 name_id=NULL

 for(i in 1:length(PFML)){ 

  PFMi=PFML[[i]]
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

fileBedOriPos="data/ORI/MCF-7/pos.bed"
dataOriPos.GR=sort(readGFBed(fileBedOriPos,SeqInfo))
fileBedOriNeg="data/ORI/MCF-7/neg.bed"
dataOriNeg.GR=sort(readGFBed(fileBedOriNeg,SeqInfo))
dataOris.GR=c(dataOriPos.GR,dataOriNeg.GR)

opts <- list()
opts[["species"]] <- 9606 # 9609 represent human
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)

bin.Mat=c(rep(1,length(dataOriPos.GR)),rep(0,length(dataOriNeg.GR)))
findMotifs(dataOris.GR,PFMatrixList,BSgenome.Hsapiens.UCSC.hg19.masked)

matMotifs=findMotifs(dataBreaks.GR,PFMatrixList,BSgenome.Hsapiens.UCSC.hg19.masked)
#save(matMotifs,file="matMotifs.RData")
#load(file="matMotifs.RData")
bin.Mat=cbind(bin.Mat,matMotifs)
rownames(bin.Mat)=1:nrow(bin.Mat)
colnames(bin.Mat)[1]="class"


# Features save as csv file
filename = "Motif"
file_varimp= paste(filename,"_Fea.csv")
write.table(bin.Mat,file=file_varimp,row.names=F,sep=',',quote=F)


#Predcting
#Training and testing sets
dataORI=data.frame(bin.Mat)
rownames(dataORI)=1:nrow(dataORI)
idxs=sample(1:nrow(dataORI),71088)#the ratio of train number 80%
dataORIlearn=dataORI[sort(idxs),]
dataORItest=dataORI[-idxs,]
RFall=ranger("class~.",data=dataORIlearn,importance="permutation")

###variable importance values

varimp=data.frame(Feature=names(RFall$variable.importance),VariableImportance=RFall$variable.importance)
varimp=varimp[order(varimp[,2],decreasing=T),]
file_varimp=paste0(filename,"_VarimpRF.csv")
write.table(varimp,file=file_varimp,row.names=F,sep=',',quote=F)

#####AUROC and AUPRC
pred <- prediction(predict(RFall,dataORItest)$predictions,as.factor(dataORItest[,1]))

auc_perf <- performance(pred, measure = "tpr", x.measure = "fpr")
auc_data <- data.frame(auc_perf@x.values,auc_perf@y.values)
aucfile = paste(filename,"_auc.csv")
write.table(auc_data,file=aucfile,row.names=F,sep=',',quote=F)
auc.tmp <- performance(pred,"auc")
auc <- as.numeric(auc.tmp@y.values)
print('auc')
print(auc)

prc_perf <- performance(pred, measure = "prec", x.measure = "rec")
prc_data <- data.frame(prc_perf@x.values,prc_perf@y.values)
prcfile = paste(filename,"_pr.csv")
write.table(prc_data,file=prcfile,row.names=F,sep=',',quote=F)
aucpr.tmp <- performance(pred,"aucpr")
pr <- as.numeric(aucpr.tmp@y.values)
print('prc')
print(pr)




