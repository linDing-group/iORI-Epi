
library(pROC)
library(glmnet)
library(ranger)
library(Matrix)
library(ROCR)


#Predcting
#Training and testing sets

filename = 'test_Motif_Fea.csv'#csv file of features 
bin.Mat <-read.csv(filename)
dataORI=data.frame(bin.Mat)
rownames(dataORI)=1:nrow(dataORI)
idxs=sample(1:nrow(dataORI),floor(nrow(dataORI)*0.8))# train:test = 8:2
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



