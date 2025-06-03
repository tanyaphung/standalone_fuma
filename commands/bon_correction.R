library(data.table)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = 'Bonferroni correction')

# Add arguments
parser$add_argument("--filedir", type = "character", required = TRUE,
                    help = "Path")

# Parse arguments
args <- parser$parse_args()

filedir=args$filedir

if (!grepl("/$", filedir)) {
  filedir <- paste0(filedir, "/")
}

if(file.exists(paste0(filedir, "magma.gsa.out"))){
  magmaset <- fread(paste0(filedir, "magma.gsa.out"), skip=3, data.table=F)
}else{
  magmaset <- fread(paste0(filedir, "magma.sets.out"), skip=3, data.table=F)
}

magmaset$Pbon <- p.adjust(magmaset$P)
magmaset <- magmaset[order(magmaset$P),]
magmaset <- subset(magmaset, select=c("FULL_NAME", "NGENES", "BETA", "BETA_STD", "SE", "P", "Pbon"))
if(length(which(magmaset$Pbon<0.05))<10){
  magmaset <- magmaset[1:10,]
}else{
  magmaset <- magmaset[magmaset$Pbon<0.05,]
}
write.table(magmaset, paste(filedir, "magma.sets.top", sep=""), sep="\t", quote=F, row.names=F)
