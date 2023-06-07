library(GEOquery)
library(survival)
library(ggplot2)
library(survminer)
library(readxl)
library(tidyr)
library(CMSclassifier)
library(limma)
#Get series matrix downloaded from GEO
clinical08=getGEO(filename="c:/users/user/desktop/gse data/GSE14333_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin08_df<-clinical08@phenoData@data

#Re formatting the data
clin08_df<-clin08_df[,c( "characteristics_ch1","Location:ch1")]
clin08_df<-separate(clin08_df,"characteristics_ch1",into=c("Location","DukesStage","Age_diag","Gender","DFS_time","DFS_cens","adjradio","adjchemo"),sep=";")
clin08_df$Location<-substr(clin08_df$Location,11,nchar(clin08_df$Location))
clin08_df$DukesStage<-substr(clin08_df$DukesStage,13,nchar(clin08_df$DukesStage))
clin08_df$Age_diag<-as.numeric(substr(clin08_df$Age_diag,11,nchar(clin08_df$Age_diag)))
clin08_df$Gender<-substr(clin08_df$Gender,9,nchar(clin08_df$Gender))
clin08_df$DFS_time<-as.numeric(substr(clin08_df$DFS_time,11,nchar(clin08_df$DFS_time)))
clin08_df$DFS_cens<-as.numeric(substr(clin08_df$DFS_cens,11,nchar(clin08_df$DFS_cens)))


#Convert death status to 1 or 0
clin08_df$deceased = clin08_df$"DFS_cens" <1


#Build survival object
stages08<-Surv(clin08_df$"DFS_time", clin08_df$deceased) 
fit08<- survfit(stages08 ~ DukesStage, data=clin08_df)
#Plot KM curve
ggsurvplot(fit08, data=clin08_df,pval=T,conf.int=T,risk.table=T, risk.table.col="strata")

cox<-coxph(stages08~DukesStage,data=clin08_df)
summary(cox)

##DE and CMS

##DG and CMS
## Get Gene expression
ge08<-clinical08@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table08<-clinical08@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge08<-subset(ge08,ifelse(!grepl("///",table08$ENTREZ_GENE_ID) & nchar(table08$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge08<-cbind(rownames(ge08),ge08)
table08<-subset(table08,ifelse(!grepl("///",table08$ENTREZ_GENE_ID) & nchar(table08$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge08<-pid2gid(as.data.frame(ge08),as.data.frame(table08))
saveRDS(ge08,"c:/users/user/desktop/GSE Data/ge08_processed.rds")

#Run CMS Classifier
cms08<-classifyCMS.RF(as.data.frame(ge08),center=T,minPosterior = .5)
saveRDS(cms08,"c:/users/user/desktop/GSE Data/cms08_processed.rds")

table(cms08$RF.predictedCMS)



##DE
#Call ge08_processed.rds next time if needed
#Remove metastasis= Na sample
ge08<-readRDS("c:/users/user/desktop/GSE Data/ge08_processed.RDS")


#Divide samples into stage first
stage08<-factor(clin08_df$`DukesStage`)

design08<-model.matrix(~stage08)

#Limma for DGE

fit08<-lmFit(ge08,design08)
fit08<-eBayes(fit08)
res08<-topTable(fit08,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit08)) 

saveRDS(fit08,"c:/users/user/desktop/GSE Data/dge08.rds")

saveRDS(res08,"c:/users/user/desktop/GSE Data/res08_processed.rds")


#Gene ontology/ pathway analysis KEGG
go08 <- limma::goana(fit08, coef=4, species="Hs")
topGO(go08,n=20,truncate="50")

kegg08<-limma::kegga(fit08,coef=4,species="Hs")
topKEGG(kegg08,n=20,truncate="50")

