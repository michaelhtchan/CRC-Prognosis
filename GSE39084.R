library(GEOquery)
library(survival)
library(ggplot2)
#Get series matrix downloaded from GEO
clinical03=getGEO(filename="c:/users/user/desktop/gse data/GSE39084_series_matrix.txt.gz")
#Create datafraome with clinical data only for easier processing 
clin03_df<-clinical03@phenoData@data
clin03_df<-clin03_df[,c("tnm.stage:ch1"  ,"death:ch1","age.year:ch1","gender:ch1" ,"surgery.date:ch1","lastnews or death.date:ch1","tp53.gene.mutation.status:ch1","kras.gene.mutation.status:ch1","braf.gene.mutation.status:ch1","pik3ca.gene.mutation.status:ch1","lynch.syndrom:ch1","cimp.status:ch1" )]

#Find the survival time given the surgery date and death date
clin03_df$time<-as.numeric(difftime(as.Date(clin03_df$`lastnews or death.date:ch1`),as.Date(clin03_df$`surgery.date:ch1`),units="days"))

#Turn death into boolean
clin03_df$`death:ch1`<-clin03_df$`death:ch1`=="1"

#Stage into numeric 
clin03_df$stage<-as.numeric(clin03_df$`tnm.stage:ch1`)

# gene mutation
clin03_df$tp53<-clin03_df$`tp53.gene.mutation.status:ch1`=="M"
clin03_df$kras<-clin03_df$`kras.gene.mutation.status:ch1`=="M"
clin03_df$braf<-clin03_df$`braf.gene.mutation.status:ch1`=="M"
clin03_df$pik3ca<-clin03_df$`pik3ca.gene.mutation.status:ch1`=="M"
clin03_df$lynch<-clin03_df$"lynch.syndrom:ch1"  =="1"

#age >66.3 ; statisticlaly signifiant- patients older than 66.3 have poorer survival
clin03_df$older66.3<-as.numeric(clin03_df$`age.year:ch1`)>66.3

#Build survival object
stages03<-Surv(clin03_df$"time"  , clin03_df$`death:ch1`) 
fit03<- survfit(stages03 ~ older66.3, data=clin03_df)
#Plot KM curve
ggsurvplot(fit03, data=clin03_df,pval=T,risk.table=T, risk.table.col="strata")

cox03<-coxph(stages03~older66.3,data=clin03_df)
summary(cox03)

#DE and CMS
##DG and CMS
## Get Gene expression
ge03<-clinical03@assayData$exprs

##Get probe id and entrez gene id and delete duplicating gene and empty gene
table03<-clinical03@featureData@data[c("ID","ENTREZ_GENE_ID")]
ge03<-subset(ge03,ifelse(!grepl("///",table03$ENTREZ_GENE_ID) & nchar(table03$ENTREZ_GENE_ID)>0,TRUE,FALSE))
ge03<-cbind(rownames(ge03),ge03)
table03<-subset(table03,ifelse(!grepl("///",table03$ENTREZ_GENE_ID) & nchar(table03$ENTREZ_GENE_ID)>0,TRUE,FALSE))

#probe id to gene id function
ge03<-pid2gid(as.data.frame(ge03),as.data.frame(table03))
saveRDS(ge03,"c:/users/user/desktop/GSE Data/ge03_processed.rds")

#Run CMS Classifier
cms03<-classifyCMS.RF(as.data.frame(ge03),center=T,minPosterior = .5)
saveRDS(cms03,"c:/users/user/desktop/GSE Data/cms03_processed.rds")

table(cms03$RF.predictedCMS)



##DE
#Call ge03_processed.rds next time if needed
#Remove metastasis= Na sample
ge03<-readRDS("c:/users/user/desktop/GSE Data/ge03_processed.RDS")


#Divide samples into stage first
stage03<-factor(clin03_df$stage)

design03<-model.matrix(~stage03)

#Limma for DGE

fit03<-lmFit(ge03,design03)
fit03<-eBayes(fit03)
res03<-topTable(fit03,number=Inf,adjust.method="none",coef=1)
summary(decideTests(fit03)) 

saveRDS(fit03,"c:/users/user/desktop/GSE Data/dge03.rds")
saveRDS(res03,"c:/users/user/desktop/GSE Data/res03_processed.rds")


#Gene ontology/ pathway analysis KEGG
go03 <- limma::goana(fit03, coef=4, species="Hs")
topGO(go03,n=20,truncate="50")

kegg03<-limma::kegga(fit03,coef=4,species="Hs")
topKEGG(kegg03,n=20,truncate="50")


