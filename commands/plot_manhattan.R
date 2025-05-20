library(qqman)
library(data.table)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-d", "--job_dir", 
    help="Path to the job directory")

args <- parser$parse_args()

# print(args$job_dir)

input = fread(paste0(args$job_dir, "input.snps"))
png(paste0(args$job_dir, "manhattan.png"), width=1200, height = 900, units = "px")
manhattan(input, chr="chr", bp="bp", snp="rsID", p="p")
dev.off()

if (file.exists(paste0(args$job_dir, "magma.genes.out"))) {
    input = fread(paste0(args$job_dir, "magma.genes.out"))
    png(paste0(args$job_dir, "gene_Manhattan.png"), width=1200, height = 900, units = "px")
    manhattan(input, chr="CHR", bp="START", snp="SYMBOL", p="P", annotatePval = 3.005e-6)
    dev.off()
}