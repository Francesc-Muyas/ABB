suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser()

# setting parameters
parser$add_argument("--positions_file", type="character", help="Tabulated coordinates of the variants to be annotated", metavar="file", nargs=1, required=TRUE)
parser$add_argument("--abb_file", type="character", help="ABB list", metavar="file", nargs=1, required=TRUE)
parser$add_argument("--out_file", type="character", help="Out file", metavar="file", nargs=1, required=TRUE)

# reading parameters
args <- parser$parse_args()

# POSITIONS 
DATA <- as.data.frame(fread(args$positions_file, header=F, colClasses = c('character','numeric')))
#DATA <- as.data.frame(fread('/users/so/fmuyas/programs/ABB/rscripts/test2/POSITIONS.temp', colClasses = 'character'))
colnames(DATA) <- c("CHROM","POS")
DATA$ORDER <- seq(1,nrow(DATA),1)

# ABB score list
ABB <- as.data.frame(fread(args$abb_file, header=F, colClasses = c('character','numeric','numeric')))
#ABB <- as.data.frame(fread("/users/so/fmuyas/programs/ABB/source/ABB_SCORE.txt",header = F, colClasses = c('character','numeric','numeric')))
colnames(ABB) <- c("CHROM","POS","ABB")

# Reducing the size of the list
ABB_temp <- ABB[ABB$POS %in% DATA$POS,]

# Merging positions with ABB scores
MERGE <- merge(DATA,ABB_temp,by=c("CHROM","POS"),all.x=T, sort = F)
MERGE <- MERGE[sort(MERGE$ORDER,decreasing = F),]
MERGE$ORDER <- NULL

# writing table 
OUTFILE <- args$out_file
write.table(MERGE, OUTFILE, sep='\t', col.names = T, row.names=F, quote = F)
