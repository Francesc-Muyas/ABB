suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(gamlss.data))
suppressPackageStartupMessages(library(gamlss.dist))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tools))

options(warn=-1)

parser <- ArgumentParser()

# setting parameters
parser$add_argument("--tsv_file", type="character", help="Tabulated file extracted from vcf", metavar="file", nargs=1, required=TRUE)
parser$add_argument("--cases", type="character", help="Cases file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("--controls", type="character", help="ABB score file", metavar="Controls file", nargs=1, required=TRUE)
parser$add_argument("--genes", type="character", help="Gene-variant association file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("--abb", type="character", help="ABB score file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("--out_prefix", type="character", help="Prefix of the output", metavar="file", nargs=1, required=TRUE)

# reading parameters
args <- parser$parse_args()

# ASSOCIATED VARIANTS 
DATA <- as.data.frame(fread(args$tsv_file))
#DATA <- as.data.frame(fread('/users/so/fmuyas/programs/ABB/rscripts/test/test.txt'))
DATA <- DATA[!grepl("\\.",DATA$GT),]

DATA$ALT_COUNT <- as.numeric(DATA$ALT_COUNT)
DATA$DP <- as.numeric(DATA$DP)
DATA$AB <- DATA$ALT_COUNT/DATA$DP

# IDs of the cases
CASES <- as.data.frame(fread(args$cases))
#CASES <- as.data.frame(fread("/users/so/fmuyas/programs/ABB/rscripts/test/CASES.txt",header = F))
CASES$CASE <- 1
colnames(CASES) <- c("SAMPLE","CASE")

# IDs of the controls
CONTROLS <- as.data.frame(fread(args$controls))
#CONTROLS <- as.data.frame(fread("/users/so/fmuyas/programs/ABB/rscripts/test/CONTROLS.txt",header = F))
CONTROLS$CASE <- 0
colnames(CONTROLS) <- c("SAMPLE","CASE")

CASE_CONTROL <- rbind(CASES,CONTROLS)

# Gene annotation of the variants 
GENES <- as.data.frame(fread(args$genes))
#GENES <- as.data.frame(fread("/users/so/fmuyas/programs/ABB/rscripts/test/VARIANTS_GENES.txt",header = F))
colnames(GENES) <- c("CHROM","POS","GENE")

# Gene annotation of the variants 
ABB <- as.data.frame(fread(args$abb))
#ABB <- as.data.frame(fread("/users/so/fmuyas/programs/ABB/source/ABB_SCORE.txt",header = T))
colnames(ABB) <- c("CHROM","POS","ABB")

# Prefix of the outfiles
OUTp <- args$out_prefix
#OUTp <- "/users/so/fmuyas/programs/ABB/rscripts/test/RESULTS/test"

######################################
##############  LOADING TOP SELECTED VARIANTS CLL
######################################

########## FUNCTIONS
#### FOR SNPS
GET_SNP_MATRIX <- function(TABLE,CHROM,POS){
  tab <- TABLE[TABLE$CHROM == CHROM & TABLE$POS == POS,]
  CALL_CONTROL <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 0 & tab$GT != 0,])
  CALL_CASES <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 1 & tab$GT != 0,])
  GAMLSS_CONTROL <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 0 & tab$GT != 0,]) - CALL_CONTROL
  GAMLSS_CASES <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 1 & tab$GT != 0,]) - CALL_CASES
  
  MATRIX <- matrix(c(CALL_CASES, CALL_CONTROL, GAMLSS_CASES, GAMLSS_CONTROL),
                   nrow=2,
                   dimnames= list(CASE_CONTROL = c("Case","Control"),
                                  MISS_CALL = c("Call","Miss")))
  
  return (MATRIX)
}

PLOT_SNP_PANEL <- function(ANALYSIS,CHROM,POSITION){
  POS <- ANALYSIS[ANALYSIS$CHROM == CHROM & ANALYSIS$POS == POSITION & ANALYSIS$NEW_GT != "0/0",]
  
  POS[POS$CASE == 1,]$CASE <- "CASE"
  POS[POS$CASE == 0,]$CASE <- "CONTROL"
  
  POS$MISS_CALL <- (POS$GT == "0/0" & POS$NEW_GT != POS$GT)*1
  POS[POS$MISS_CALL == 1,]$MISS_CALL <- "MISSED"
  POS[POS$MISS_CALL == 0,]$MISS_CALL <- "CALLED"
  
  POS_CASES <- POS[POS$CASE == "CASE",]
  POS_CONTROLS <- POS[POS$CASE == "CONTROL",]
  
  
  GT_AB_CASES <- ggplot(POS_CASES,aes(x=as.numeric(AB),fill=as.factor(MISS_CALL))) +
    geom_histogram(alpha=0.5, position="identity",binwidth = 0.025, color="darkgreen", aes(y = (..count..)/sum(..count..)*100)) +
    xlim(0,1) +
    ggtitle(paste("CASES:",paste(CHROM,POSITION,sep=";"),sep=" ")) +
    guides(fill = guide_legend(title = "")) +
    xlab('Allele balance') + 
    ylab('Proportion (%)') + 
    theme_bw()
  
  GT_AB_CONTROL <- ggplot(POS_CONTROLS,aes(x=as.numeric(AB),fill=as.factor(MISS_CALL))) +
    geom_histogram(alpha=0.5, position="identity",binwidth = 0.025, color="darkgreen",aes(y = (..count..)/sum(..count..)*100)) +
    xlim(0,1) +
    ggtitle(paste("CONTROLS:",paste(CHROM,POSITION,sep=";"),sep=" ")) +
    guides(fill = guide_legend(title = "")) +
    xlab('Allele balance') + 
    ylab('Proportion (%)') + 
    theme_bw()  
  
  # PANEL PLOT
  
  # Get the gtables
  GA <- ggplotGrob(GT_AB_CONTROL)
  GB <- ggplotGrob(GT_AB_CASES)
  
  # Set the widths
  GA$widths <- GB$widths
  
  # Arrange the two charts.
  # The legend boxes are centered
  GT_AB <- grid.arrange(GA, GB, nrow = 2)
  
  ggsave(file=paste("PANEL",CHROM,POSITION,"pdf",sep="."),plot=GT_AB, useDingbats=FALSE)
  
  #return(GT_AB)
}

#### FOR GENES
GET_GENE_MATRIX <- function(TABLE,gene){
  tab <- TABLE[TABLE$GENE == gene,]
  CALL_CONTROL <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 0 & tab$GT != 0,])
  CALL_CASES <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 1 & tab$GT != 0,])
  GAMLSS_CONTROL <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 0 & tab$GT != 0,]) - CALL_CONTROL
  GAMLSS_CASES <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 1 & tab$GT != 0,]) - CALL_CASES
  
  MATRIX <- matrix(c(CALL_CASES, CALL_CONTROL, GAMLSS_CASES, GAMLSS_CONTROL),
                   nrow=2,
                   dimnames= list(CASE_CONTROL = c("Case","Control"),
                                  MISS_CALL = c("Call","Gamlss")))
  
  return (MATRIX)
}

GET_GENE_MATRIX_ASSOC_BROAD <- function(TABLE,gene){
  tab <- TABLE[TABLE$GENE == gene,]
  CALL_CONTROL <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 0 & tab$GT != 0,])
  CALL_CASES <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 1 & tab$GT != 0,])
  NO_CALL_CONTROL <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 0 & tab$GT == 0,])
  NO_CALL_CASES <- nrow(tab[tab$TYPE == "BROAD" & tab$CASE == 1 & tab$GT == 0,])
  
  MATRIX <- matrix(c(CALL_CASES, CALL_CONTROL, NO_CALL_CASES, NO_CALL_CONTROL),
                   nrow=2,
                   dimnames= list(CASE_CONTROL = c("Case","Control"),
                                  NO_CALL = c("Call","No_call")))
  
  return (MATRIX)
}

GET_GENE_MATRIX_ASSOC_GAMLSS <- function(TABLE,gene){
  tab <- TABLE[TABLE$GENE == gene,]
  CALL_CONTROL <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 0 & tab$GT != 0,])
  CALL_CASES <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 1 & tab$GT != 0,])
  NO_CALL_CONTROL <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 0 & tab$GT == 0,])
  NO_CALL_CASES <- nrow(tab[tab$TYPE == "GAMLSS" & tab$CASE == 1 & tab$GT == 0,])
  
  MATRIX <- matrix(c(CALL_CASES, CALL_CONTROL, NO_CALL_CASES, NO_CALL_CONTROL),
                   nrow=2,
                   dimnames= list(CASE_CONTROL = c("Case","Control"),
                                  NO_CALL = c("Call","No_call")))
  
  return (MATRIX)
}

PLOT_GENE_PANEL <- function(ANALYSIS,GENE_ID){
  POS <- ANALYSIS[ANALYSIS$GENE == GENE_ID & ANALYSIS$NEW_GT != "0/0",]
  
  POS[POS$CASE == 1,]$CASE <- "CASE"
  POS[POS$CASE == 0,]$CASE <- "CONTROL"
  
  POS$MISS_CALL <- (POS$GT == "0/0" & POS$NEW_GT != POS$GT)*1
  POS[POS$MISS_CALL == 1,]$MISS_CALL <- "MISSED"
  POS[POS$MISS_CALL == 0,]$MISS_CALL <- "CALLED"
  
  POS_CASES <- POS[POS$CASE == "CASE",]
  POS_CONTROLS <- POS[POS$CASE == "CONTROL",]
  
  
  GT_AB_CASES <- ggplot(POS_CASES,aes(x=as.numeric(AB),fill=as.factor(MISS_CALL))) +
    geom_histogram(alpha=0.5, position="identity",binwidth = 0.025, color="darkgreen", aes(y = (..count..)/sum(..count..)*100)) +
    xlim(0,1) +
    ggtitle(paste("CASES:",GENE_ID,sep=" ")) +
    guides(fill = guide_legend(title = "")) +
    xlab('Allele balance') + 
    ylab('Proportion (%)') + 
    theme_bw()
  
  GT_AB_CONTROL <- ggplot(POS_CONTROLS,aes(x=as.numeric(AB),fill=as.factor(MISS_CALL))) +
    geom_histogram(alpha=0.5, position="identity",binwidth = 0.025, color="darkgreen",aes(y = (..count..)/sum(..count..)*100)) +
    xlim(0,1) +
    ggtitle(paste("CONTROLS:",GENE_ID,sep=" ")) +
    guides(fill = guide_legend(title = "")) +
    xlab('Allele balance') + 
    ylab('Proportion (%)') + 
    theme_bw()  
  
  # PANEL PLOT
  
  # Get the gtables
  GA <- ggplotGrob(GT_AB_CONTROL)
  GB <- ggplotGrob(GT_AB_CASES)
  
  # Set the widths
  GA$widths <- GB$widths
  
  # Arrange the two charts.
  # The legend boxes are centered
  GT_AB <- grid.arrange(GA, GB, nrow = 2)
  
  ggsave(file=paste("PANEL",GENE_ID,"pdf",sep="."),plot=GT_AB, useDingbats=FALSE)
  
  #return(GT_AB)
}

######################################
##############  MERGING FILES
######################################

MERGE_0 <- merge(DATA,CASE_CONTROL,by="SAMPLE")

MERGE <- merge(MERGE_0,GENES,by=c("CHROM","POS"))

ANALYSIS <- MERGE

####### PARAMETERS 0/0 GT

MU_0 <- 0.03270537
SIGMA_0 <- 0.1452151
NU_0 <- 1.689303e+01
TAU_0 <- 1.011730e-14

## RE-GENOTYPING
ANALYSIS_0 <- ANALYSIS[ANALYSIS$AB == 0,]
ANALYSIS_0$P_HOM_REF <- 1
ANALYSIS_0$NEW_GT <- "0/0"

ANALYSIS_1 <- ANALYSIS[ANALYSIS$AB > 0,]
ANALYSIS_1$P_HOM_REF <- pBEINF((ANALYSIS_1$ALT_COUNT - 1)/ANALYSIS_1$DP, MU_0, SIGMA_0, NU_0, TAU_0, lower.tail=FALSE)
ANALYSIS_1$NEW_GT <- apply(ANALYSIS_1, 1, function(x) ifelse(x["GT"] != "0/0" , x["GT"], if (as.numeric(x["P_HOM_REF"]) > 0.05){ "0/0"}else{"0/1"}))

ANALYSIS <- rbind(ANALYSIS_0,ANALYSIS_1)

ANALYSIS$GT_BROAD <- (ANALYSIS$GT != "0/0")*1
ANALYSIS$GT_GAMLSS <- (ANALYSIS$NEW_GT != "0/0")*1

######################
## SNP ABB-ASSOCIATION ANALYSIS
######################

SNP_BROAD <- aggregate(GT_BROAD ~ CHROM+POS+SAMPLE+CASE, data = ANALYSIS, sum)
colnames(SNP_BROAD)[5] <- "GT"

SNP_GAMLSS <- aggregate(GT_GAMLSS ~ CHROM+POS+SAMPLE+CASE, data = ANALYSIS, sum)
colnames(SNP_GAMLSS)[5] <- "GT"

SNP_BROAD$TYPE <- "BROAD"
SNP_GAMLSS$TYPE <- "GAMLSS"

SNP <- rbind(SNP_BROAD,SNP_GAMLSS)

## Positions to be analysed
POSITIONS <- unique(SNP[c("CHROM", "POS")])

## ABB intersection
temp_ABB <- ABB[ABB$POS %in% POSITIONS$POS,]
POSITIONS <- merge(POSITIONS,temp_ABB,by=c("CHROM","POS"),all.x=T)

## Running the test for each position
POSITIONS$P_VAL1 <- apply(POSITIONS,1,function(x) fisher.test(GET_SNP_MATRIX(SNP,as.numeric(x["CHROM"]),as.numeric(x["POS"])))[[1]]) 

# FDR correction
POSITIONS$P_VAL1_adj <- p.adjust(POSITIONS$P_VAL1,"fdr")
POSITIONS <- merge(POSITIONS,GENES,by=c("CHROM","POS"))
POSITIONS_SIG <- POSITIONS[POSITIONS$P_VAL1_adj < 0.1,]

######################
## GENE ABB-ASSOCIATION ANALYSIS
######################

GENE_BROAD <- aggregate(GT_BROAD ~ GENE+SAMPLE+CASE, data = ANALYSIS, sum)
colnames(GENE_BROAD)[4] <- "GT"
GENE_GAMLSS <- aggregate(GT_GAMLSS ~ GENE+SAMPLE+CASE, data = ANALYSIS, sum)
colnames(GENE_GAMLSS)[4] <- "GT"

GENE_BROAD$TYPE <- "BROAD"
GENE_BROAD$GT <- (GENE_BROAD$GT > 0)*1
GENE_GAMLSS$TYPE <- "GAMLSS"
GENE_GAMLSS$GT <- (GENE_GAMLSS$GT > 0)*1

GENE <- rbind(GENE_BROAD,GENE_GAMLSS)

## GENES to be analysed
genes <- unique(GENE["GENE"])

##### Ratio test
## Running the test for each gene
genes$P_VAL1 <- apply(genes,1,function(x) fisher.test(GET_GENE_MATRIX(GENE,x["GENE"]))[[1]]) 

# FDR correction
genes$P_VAL1_adj <- p.adjust(genes$P_VAL1,"fdr")

##### Normal chisq test (Provided results)
## Running the test for each gene
genes$P_VALn <- apply(genes,1,function(x) chisq.test(GET_GENE_MATRIX_ASSOC_BROAD(GENE,x["GENE"]))[[3]]) 

# FDR correction
genes$P_VALn_adj <- p.adjust(genes$P_VALn,"fdr")

##### Chisq test with the regenotyped samples
## Running the test for each gene
genes$P_VAL2 <- apply(genes,1,function(x) chisq.test(GET_GENE_MATRIX_ASSOC_GAMLSS(GENE,x["GENE"]))[[3]]) 

# FDR correction
genes$P_VAL2_adj <- p.adjust(genes$P_VAL2,"fdr")

##### Chisq test ignoring significant positions obtained before
ANALYSIS_FILTERED <- ANALYSIS[!(paste(ANALYSIS$CHROM,ANALYSIS$POS,sep=';') %in% paste(POSITIONS_SIG$CHROM,POSITIONS_SIG$POS,sep=';')),]
GENE_FILTERED <- aggregate(GT_BROAD ~ GENE+SAMPLE+CASE, data = ANALYSIS_FILTERED, sum)
colnames(GENE_FILTERED)[4] <- "GT"

GENE_FILTERED$TYPE <- "BROAD"
GENE_FILTERED$GT <- (GENE_FILTERED$GT > 0)*1

## Running the test for each gene
genes$P_VAL3 <- apply(genes,1,function(x) chisq.test(GET_GENE_MATRIX_ASSOC_BROAD(GENE_FILTERED,x["GENE"]))[[3]]) 

# FDR correction
genes$P_VAL3_adj <- p.adjust(genes$P_VAL3,"fdr")

# Genes prone to be errors
genes_SIG <- genes[genes$P_VAL1_adj < 0.1 | (genes$P_VALn_adj < 0.1 & (genes$P_VAL2_adj > 0.1 | genes$P_VAL3_adj > 0.1)),]

# Selecting the main columns
genes <- genes[,c("GENE","P_VAL1_adj","P_VAL2_adj","P_VAL3_adj")]
genes_SIG <- genes_SIG[,c("GENE","P_VAL1_adj","P_VAL2_adj","P_VAL3_adj")]

colnames(genes) <- c("GENE","Missed-Called_ratio(FDR)","Association_regenotyped(FDR)","Association_ABB(FDR)")
colnames(genes_SIG) <- c("GENE","Missed-Called_ratio(FDR)","Association_regenotyped(FDR)","Association_ABB(FDR)")


POSITIONS <- POSITIONS[,c("CHROM","POS","ABB","P_VAL1_adj","GENE")]
POSITIONS_SIG <- POSITIONS_SIG[,c("CHROM","POS","ABB","P_VAL1_adj","GENE")]

colnames(POSITIONS) <- c("CHROM","POS","ABB","Missed-Called_ratio(FDR)","GENE")
colnames(POSITIONS_SIG) <- c("CHROM","POS","ABB","Missed-Called_ratio(FDR)","GENE")

#######################
# WRITING RESULTS
#######################

write.table(POSITIONS,file = paste(OUTp,"SNPs.txt",sep='.'),quote=F,sep='\t',col.names = T,row.names = F)
write.table(genes,file = paste(OUTp,"GENEs.txt",sep='.'),quote=F,sep='\t',col.names = T,row.names = F)

write.table(POSITIONS_SIG,file = paste(OUTp,"SNPs.significant.txt",sep='.'),quote=F,sep='\t',col.names = T,row.names = F)
write.table(genes_SIG,file = paste(OUTp,"GENEs.significant.txt",sep='.'),quote=F,sep='\t',col.names = T,row.names = F)


### Create folders for plots
OUTP <- file_path_as_absolute(dirname(OUTp))
FOLDER <- paste(OUTP,"PLOTS",sep='/')
dir.create(FOLDER, showWarnings = FALSE)
print (FOLDER)

SUBFOLDER_SNP <- paste(FOLDER,"SNP",sep='/')
print (SUBFOLDER_SNP)
dir.create(SUBFOLDER_SNP, showWarnings = FALSE)

SUBFOLDER_GENE <- paste(FOLDER,"GENE",sep='/')
print (SUBFOLDER_GENE)
dir.create(SUBFOLDER_GENE, showWarnings = FALSE)

##########
# PLOTING SIGNIFICANT POSITIONS
##########

## SNPs
setwd(SUBFOLDER_SNP)
apply(POSITIONS_SIG,1,function(x) PLOT_SNP_PANEL(ANALYSIS,as.numeric(x["CHROM"]),as.numeric(x["POS"])))

## GENES
setwd(SUBFOLDER_GENE)
apply(genes_SIG,1,function(x) PLOT_GENE_PANEL(ANALYSIS,x["GENE"]))

