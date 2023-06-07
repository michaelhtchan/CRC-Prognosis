library(GEOquery)
library(survival)
library(ggplot2)
library(survminer)
#Get series matrix downloaded from GEO
clinical01=getGEO(filename="c:/users/user/desktop/gse data/GSE17536_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin01_df<-clinical01@phenoData@data
clin01_df<-clin01_df[,c("ajcc_stage:ch1"                                              
                        ,"dfs_event (disease free survival; cancer recurrence):ch1"    
                        ,"dfs_time:ch1"                                                
                        , "dss_event (disease specific survival; death from cancer):ch1"
                        , "dss_time:ch1"                                                
                        , "ethnicity:ch1"                                               
                        , "gender:ch1"                                                  
                        , "grade:ch1"                                                   
                        , "overall survival follow-up time:ch1"                         
                        , "overall_event (death from any cause):ch1",       "age:ch1")]
#Convert death status to 1 or 0
clin01_df$deceased = clin01_df$"overall_event (death from any cause):ch1" == "death"
#Convert stage into new column for easier visualization
clin01_df$stage=as.numeric(clin01_df$"ajcc_stage:ch1")
#time into numeric
clin01_df$"overall survival follow-up time:ch1"<-as.numeric(clin01_df$"overall survival follow-up time:ch1")

clin01_df$stage<-as.numeric(clin01_df$`ajcc_stage:ch1`)
#Build survival object
stages01<-Surv(clin01_df$"overall survival follow-up time:ch1", clin01_df$deceased) 
fit01<- survfit(stages01 ~ stage, data=clin01_df)
#Plot KM curve
ggsurvplot(fit01, data=clin01_df,pval=T,conf.int=T,risk.table=T, risk.table.col="strata")

cox<-coxph(stages01~stage,data=clin01_df)
summary(cox)


##DG and CMS
## Get Gene expression
ge01<-clinical01@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table01<-clinical01@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge01<-subset(ge01,ifelse(!grepl("///",table01$ENTREZ_GENE_ID) & nchar(table01$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge01<-cbind(rownames(ge01),ge01)
table01<-subset(table01,ifelse(!grepl("///",table01$ENTREZ_GENE_ID) & nchar(table01$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge01<-pid2gid(as.data.frame(ge01),as.data.frame(table01))
saveRDS(ge01,"c:/users/user/desktop/GSE Data/ge01_processed.rds")

#Run CMS Classifier
cms01<-classifyCMS.RF(as.data.frame(ge01),center=T,minPosterior = .5)
saveRDS(cms01,"c:/users/user/desktop/GSE Data/cms01_processed.rds")

table(cms01$RF.predictedCMS)



##DE
#Call ge01_processed.rds next time if needed
#Remove metastasis= Na sample
ge01<-readRDS("c:/users/user/desktop/GSE Data/ge01_processed.RDS")


#Divide samples into stage first
stage01<-factor(clin01_df$stage)

design01<-model.matrix(~stage01)

#Limma for DGE

fit01<-lmFit(ge01,design01)
fit01<-eBayes(fit01)
res01<-topTable(fit01,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit01)) 

saveRDS(fit01,"c:/users/user/desktop/GSE Data/dge01.rds")
saveRDS(res01,"c:/users/user/desktop/GSE Data/res01_processed.rds")


#Gene ontology/ pathway analysis KEGG
go01 <- limma::goana(fit01, coef=4, species="Hs")
topGO(go01,n=20,truncate="50")

kegg01<-limma::kegga(fit01,coef=4,species="Hs")
topKEGG(kegg01,n=20,truncate="50")

