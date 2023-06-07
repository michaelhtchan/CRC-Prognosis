library(GEOquery)
library(survival)
library(ggplot2)
library(survminer)
library(readxl)
#Get series matrix downloaded from GEO
clinical07=getGEO(filename="c:/users/user/desktop/gse data/GSE37892_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin07_df<-clinical07@phenoData@data

#Becoz no OS, we will use DFS for survival analysis
clin07_df<-clin07_df[,c("age at diagnosis:ch1" ,          "center:ch1"  ,                  
"date at distant metastasis:ch1", "date at last contact:ch1",       "date at surgery:ch1" ,          
 "gender:ch1"            ,         "localisation:ch1"    ,           "Stage:ch1" )]
#Convert death status to 1 or 0
clin07_df$deceased = clin07_df$"date at distant metastasis:ch1" != 0

#time into numeric by subtracting two different time frame
clin07_df$time<-as.numeric(difftime(as.Date(clin07_df$"date at distant metastasis:ch1",format="%Y-%m-%d"),as.Date(clin07_df$"date at surgery:ch1",format="%Y-%m-%d"),units="days"))

names(clin07_df)[names(clin07_df)=="localisation:ch1"]<-"localisation"

#Build survival object
stages07<-Surv(clin07_df$"time", clin07_df$deceased) 
fit07<- survfit(stages07 ~ localisation, data=clin07_df)
#Plot KM curve
ggsurvplot(fit07, data=clin07_df,pval=T,conf.int=T,risk.table=T, risk.table.col="strata")

cox<-coxph(stages07~deceased,data=clin07_df)
summary(cox)

##DG and CMS
## Get Gene expression
ge07<-clinical07@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table07<-clinical07@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge07<-subset(ge07,ifelse(!grepl("///",table07$ENTREZ_GENE_ID) & nchar(table07$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge07<-cbind(rownames(ge07),ge07)
table07<-subset(table07,ifelse(!grepl("///",table07$ENTREZ_GENE_ID) & nchar(table07$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge07<-pid2gid(as.data.frame(ge07),as.data.frame(table07))
saveRDS(ge07,"c:/users/user/desktop/GSE Data/ge07_processed.rds")

#Run CMS Classifier
cms07<-classifyCMS.RF(as.data.frame(ge07),center=T,minPosterior = .5)
saveRDS(cms07,"c:/users/user/desktop/GSE Data/cms07_processed.rds")

table(cms07$RF.predictedCMS)



##DE
#Call ge07_processed.rds next time if needed
#Remove metastasis= Na sample
ge07<-readRDS("c:/users/user/desktop/GSE Data/ge07_processed.RDS")


#Divide samples into stage first
stage07<-factor(clin07_df$`Stage:ch1`)

design07<-model.matrix(~stage07)

#Limma for DGE

fit07<-lmFit(ge07,design07)
fit07<-eBayes(fit07)
res07<-topTable(fit07,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit07)) 

saveRDS(fit07,"c:/users/user/desktop/GSE Data/dge07.rds")

saveRDS(res07,"c:/users/user/desktop/GSE Data/res07_processed.rds")


#Gene ontology/ pathway analysis KEGG
go07 <- limma::goana(fit07, species="Hs")
topGO(go07,n=20,truncate="50")

kegg07<-limma::kegga(fit07,species="Hs")
topKEGG(kegg07,n=20,truncate="50")


