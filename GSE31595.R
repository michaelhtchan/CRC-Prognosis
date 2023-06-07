library(GEOquery)
library(survival)
library(ggplot2)
#Get series matrix downloaded from GEO
clinical02=getGEO(filename="c:/users/user/desktop/gse data/GSE31595_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin02_df<-clinical02@phenoData@data
clin02_df<-clin02_df[,c("tnm-stage:ch1" ,
                        "age:ch1","death:ch1", "recurrence:ch1" ,"gender:ch1" ,"ethnicity:ch1","relapse free survival time:ch1"   
                        )]
#Convert death status to 1 or 0
clin02_df$deceased = clin02_df$"death:ch1" == "yes"
#Convert stage into new column for easier visualization
clin02_df$stage<-ifelse(clin02_df$`tnm-stage:ch1`=="II",as.numeric(2),ifelse(clin02_df$`tnm-stage:ch1`=="III",as.numeric(3),ifelse(clin02_df$`tnm-stage:ch1`=="IV",as.numeric(4),"")))

#new column for male or female
clin02_df$gender<-ifelse(substr(clin02_df$`gender:ch1`,2,2)=="a","M",ifelse(substr(clin02_df$`gender:ch1`,2,2)=="e","F",""))

#new column for >=73 or <73 where 73 = median 
clin02_df$older73<-(as.numeric(clin02_df$`age:ch1`)>72)

#RFS formatting to remove unit in original data and convert into double
clin02_df$"relapse free survival time:ch1"<-substr(clin02_df$"relapse free survival time:ch1",1,nchar(clin02_df$"relapse free survival time:ch1")-5)
clin02_df$"relapse free survival time:ch1"<-as.numeric(clin02_df$"relapse free survival time:ch1")

#Build survival object
stages02<-Surv(clin02_df$"relapse free survival time:ch1"  , clin02_df$deceased) 
fit02<- survfit(stages02 ~ older73, data=clin02_df)

#Plot KM curve
ggsurvplot(fit02, data=clin02_df,pval=T,conf.int=T,risk.table=T, risk.table.col="strata")

cox02<-coxph(stages02~stage,data=clin02_df)
summary(cox02)

##DG and CMS
## Get Gene expression
ge02<-clinical02@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table02<-clinical02@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge02<-subset(ge02,ifelse(!grepl("///",table02$ENTREZ_GENE_ID) & nchar(table02$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge02<-cbind(rownames(ge02),ge02)
table02<-subset(table02,ifelse(!grepl("///",table02$ENTREZ_GENE_ID) & nchar(table02$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge02<-pid2gid(as.data.frame(ge02),as.data.frame(table02))
saveRDS(ge02,"c:/users/user/desktop/GSE Data/ge02_processed.rds")

#Run CMS Classifier
cms02<-classifyCMS.RF(as.data.frame(ge02),center=T,minPosterior = .5)
saveRDS(cms02,"c:/users/user/desktop/GSE Data/cms02_processed.rds")

table(cms02$RF.predictedCMS)



##DE
#Call ge02_processed.rds next time if needed
#Remove metastasis= Na sample
ge02<-readRDS("c:/users/user/desktop/GSE Data/ge02_processed.RDS")


#Divide samples into stage first
stage02<-factor(clin02_df$stage)

design02<-model.matrix(~stage02)

#Limma for DGE

fit02<-lmFit(ge02,design02)
fit02<-eBayes(fit02)
res02<-topTable(fit02,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit02)) 

saveRDS(fit02,"c:/users/user/desktop/GSE Data/dge02.rds")
saveRDS(res02,"c:/users/user/desktop/GSE Data/res02_processed.rds")


#Gene ontology/ pathway analysis KEGG
go02 <- limma::goana(fit02, coef=4, species="Hs")
topGO(go02,n=20,truncate="50")

kegg02<-limma::kegga(fit02,coef=4,species="Hs")
topKEGG(kegg02,n=20,truncate="50")


