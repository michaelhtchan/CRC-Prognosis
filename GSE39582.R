library(GEOquery)
library(survival)
library(ggplot2)
library(survminer)
library(edgeR)
library(affy)
library(CMSclassifier)
#Get series matrix downloaded from GEO
clinical04=getGEO(filename="c:/users/user/desktop/gse data/GSE39582_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin04_df<-clinical04@phenoData@data
clin04_df<-clin04_df[,c("tnm.stage:ch1"  ,"age.at.diagnosis (year):ch1","Sex:ch1", "os.delay (months):ch1", "os.event:ch1"  
)]
#Convert death status to 1 or 0
clin04_df$`os.event:ch1` = clin04_df$`os.event:ch1` == "1"

clin04_df$`os.delay (months):ch1`<-as.numeric(clin04_df$`os.delay (months):ch1`)

clin04_df$stage<-as.numeric(clin04_df$`tnm.stage:ch1`)

#Build survival object
stages04<-Surv(clin04_df$`os.delay (months):ch1`  , clin04_df$`os.event:ch1`) 
fit04<- survfit(stages04 ~ stage, data=clin04_df)

#Plot KM curve
ggsurvplot(fit04, data=clin04_df,pval=T,risk.table=T, risk.table.col="strata")

cox04<-coxph(stages04~stage,data=clin04_df)
summary(cox04)

##CMS classification
pid2gid = function(exp, probeid_geneid) 
{
  exp=as.matrix(exp)
  index = which(rowSums(is.na(exp)) == 0)
  subexp = exp[index, ]   
  probeid_geneid = unique(probeid_geneid)
  index1 = which(rowSums(is.na(probeid_geneid)) == 0)
  index2 = which(!is.na(as.numeric(as.character(probeid_geneid[, 2]))))
  id_index = intersect(index1, index2)
  pid_gid = probeid_geneid[id_index, ]
  pid = intersect(subexp[, 1], pid_gid[, 1])
  index2 = match(pid, subexp[, 1])
  subexp1 = subexp[index2, -1]
  mode(subexp1) = "numeric"
  index3 = match(pid, pid_gid[, 1])
  ### calculate mean value for genes annotated with multiple probes
  subgeneid = factor(pid_gid[index3, 2])
  geneexp = apply(subexp1, 2, function(x) tapply(x, subgeneid, mean))
  return(geneexp)
} 




clinical04=getGEO(filename="c:/users/user/desktop/gse data/GSE39582_series_matrix.txt.gz")

## Get Gene expression
ge04<-clinical04@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table04<-clinical04@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge04<-subset(ge04,ifelse(!grepl("///",table04$ENTREZ_GENE_ID) & nchar(table04$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge04<-cbind(rownames(ge04),ge04)
table04<-subset(table04,ifelse(!grepl("///",table04$ENTREZ_GENE_ID) & nchar(table04$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge04<-pid2gid(as.data.frame(ge04),as.data.frame(table04))
saveRDS(ge04,"c:/users/user/desktop/GSE Data/ge04_processed.rds")

#Run CMS Classifier
cms04<-classifyCMS.RF(as.data.frame(ge04),center=T,minPosterior = .5)
saveRDS(cms04,"c:/users/user/desktop/GSE Data/cms04_processed.rds")

table(cms04$RF.predictedCMS)



##DE
#Call ge04_processed.rds next time if needed
#Remove stage= Na sample
ge04<-readRDS("c:/users/user/desktop/GSE Data/ge04_processed.RDS")
ge04<-ge04[,-which(!((clin04_df$`tnm.stage:ch1`)=="1" |(clin04_df$`tnm.stage:ch1`)=="4") )]
clin04_df<-clin04_df[-which(!((clin04_df$`tnm.stage:ch1`)=="1" |(clin04_df$`tnm.stage:ch1`)=="4") ),]


#Divide samples into stage first
stage04<-factor(clin04_df$`tnm.stage:ch1`)

design04<-model.matrix(~stage04)

#Limma for DGE

fit04<-lmFit(ge04,design04)
fit04<-eBayes(fit04)
res04<-topTable(fit04,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit04))

saveRDS(fit04,"c:/users/user/desktop/GSE Data/dge04.rds")

saveRDS(res04,"c:/users/user/desktop/GSE Data/res04_processed.rds")


#Gene ontology/ pathway analysis KEGG
go04 <- limma::goana(fit04, species="Hs")
topGO(go04,truncate="50")

kegg04<-limma::kegga(fit04,pecies="Hs")
topKEGG(kegg04,n=50,truncate="50")



