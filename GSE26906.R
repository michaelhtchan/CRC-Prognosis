library(GEOquery)
library(survival)
library(ggplot2)
library(survminer)
library(readxl)
library(CMSclassifier)
library(limma)
#Get series matrix downloaded from GEO
clinical05=getGEO(filename="c:/users/user/desktop/gse data/GSE26906_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin05_df<-clinical05@phenoData@data

#Get data from suppl. table then merge clin05_df to it
supp05<-as.data.frame(read_excel("c:/users/user/downloads/tlo0502_0072SD1.xls",sheet="Suppl table S2"))
names(supp05)[names(supp05)=="Tumor sample_ID"]<-"title"
#typo of original data
names(supp05)[names(supp05)=="Data at diagnosis"]<-"Date at diagnosis"

clin05_df<-merge(clin05_df,supp05,by="title")

clin05_df<-clin05_df[,c("tissue:ch1","Date at diagnosis","Date at last contact/death", "age at diagnosis:ch1" ,"Cause of death" )]
#Convert death status to 1 or 0
clin05_df$deceased = clin05_df$"Cause of death" != 0

#time into numeric by subtracting two different time frame
clin05_df$time<-as.numeric(-difftime(as.Date(clin05_df$"Date at diagnosis"),as.Date(clin05_df$"Date at last contact/death"),units="days"))

#Age at diagnosis
clin05_df$older68.5<-as.numeric(clin05_df$`age at diagnosis:ch1`)>median(as.numeric(clin05_df$`age at diagnosis:ch1`))

#Build survival object
stages05<-Surv(clin05_df$"time", clin05_df$deceased) 
fit05<- survfit(stages05 ~ older68.5, data=clin05_df)
#Plot KM curve
ggsurvplot(fit05, data=clin05_df,pval=T,conf.int=T,risk.table=T, risk.table.col="strata")

cox<-coxph(stages05~tiss,data=clin05_df)
summary(cox)


##DG and CMS
## Get Gene expression
ge05<-clinical05@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table05<-clinical05@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge05<-subset(ge05,ifelse(!grepl("///",table05$ENTREZ_GENE_ID) & nchar(table05$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge05<-cbind(rownames(ge05),ge05)
table05<-subset(table05,ifelse(!grepl("///",table05$ENTREZ_GENE_ID) & nchar(table05$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge05<-pid2gid(as.data.frame(ge05),as.data.frame(table05))
saveRDS(ge05,"c:/users/user/desktop/GSE Data/ge05_processed.rds")

#Run CMS Classifier
cms05<-classifyCMS.RF(as.data.frame(ge05),center=T,minPosterior = .5)
saveRDS(cms05,"c:/users/user/desktop/GSE Data/cms05_processed.rds")

table(cms05$RF.predictedCMS)



##DE
#Call ge05_processed.rds next time if needed
#Remove metastasis= Na sample
ge05<-readRDS("c:/users/user/desktop/GSE Data/ge05_processed.RDS")


#Divide samples into stage first
stage05<-factor(clin05_df$`tissue:ch1`)

design05<-model.matrix(~stage05)

#Limma for DGE

fit05<-lmFit(ge05,design05)
fit05<-eBayes(fit05)
res05<-topTable(fit05,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit05)) 

saveRDS(fit05,"c:/users/user/desktop/GSE Data/dge05.rds")

saveRDS(res05,"c:/users/user/desktop/GSE Data/res05_processed.rds")


#Gene ontology/ pathway analysis KEGG
go05 <- limma::goana(fit05, coef=4, species="Hs")
topGO(go05,n=20,truncate="50")

kegg05<-limma::kegga(fit05,coef=4,species="Hs")
topKEGG(kegg05,n=20,truncate="50")
