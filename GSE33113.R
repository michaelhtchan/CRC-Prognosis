library(GEOquery)
library(survival)
library(ggplot2)
library(survminer)
library(readxl)
library(CMSclassifier)
library(limma)

###????? still no?

#Get series matrix downloaded from GEO
clinical06=getGEO(filename="c:/users/user/desktop/gse data/GSE33113_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin06_df<-clinical06@phenoData@data

#Becoz no OS, we will use DFS for survival analysis
clin06_df<-clin06_df[,c("time to meta or recurrence:ch1","meta or recurrence within 3 years:ch1","Sex:ch1"    )]

clin06_df$time<-as.numeric(clin06_df$`time to meta or recurrence:ch1`)
clin06_df$sex<-clin06_df$`Sex:ch1`
clin06_df$meta<-clin06_df$`meta or recurrence within 3 years:ch1`=="yes"

#Build survival object
stages06<-Surv(clin06_df$"time", clin06_df$meta) 
fit06<- survfit(stages06 ~ sex, data=clin06_df)
#Plot KM curve
ggsurvplot(fit06, data=clin06_df,pval=T,conf.int=T,risk.table=T, risk.table.col="strata")

cox<-coxph(stages06~sex,data=clin06_df)
summary(cox)


## Get Gene expression
ge06<-clinical06@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table06<-clinical06@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge06<-subset(ge06,ifelse(!grepl("///",table06$ENTREZ_GENE_ID) & nchar(table06$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge06<-cbind(rownames(ge06),ge06)
table06<-subset(table06,ifelse(!grepl("///",table06$ENTREZ_GENE_ID) & nchar(table06$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge06<-pid2gid(as.data.frame(ge06),as.data.frame(table06))
saveRDS(ge06,"c:/users/user/desktop/GSE Data/ge06_processed.rds")

#Run CMS Classifier
cms06<-classifyCMS.RF(as.data.frame(ge06),center=T,minPosterior = .5)
saveRDS(cms06,"c:/users/user/desktop/GSE Data/cms06_processed.rds")

table(cms06$RF.predictedCMS)



##DE
#Call ge06_processed.rds next time if needed
#Remove metastasis= Na sample
ge06<-readRDS("c:/users/user/desktop/GSE Data/ge06_processed.RDS")
ge06<-ge06[,1:90]


#Divide samples into stage first
stage06<-factor(clin06_df$`meta or recurrence within 3 years:ch1`)

design06<-model.matrix(~stage06)

#Limma for DGE

fit06<-lmFit(ge06,design06)
fit06<-eBayes(fit06)
res06<-topTable(fit06,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit06)) #Stage = metastasis or not as actual stage not provided


saveRDS(res06,"c:/users/user/desktop/GSE Data/res06_processed.rds")


#Gene ontology/ pathway analysis KEGG
go06 <- limma::goana(fit06, species="Hs")
topGO(go06,n=20,truncate="50")

kegg06<-limma::kegga(fit06,species="Hs")
topKEGG(kegg06,n=20,truncate="50")

