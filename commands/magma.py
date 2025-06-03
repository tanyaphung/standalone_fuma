import os
import pandas as pd
import subprocess
import sys
import logging


class MagmaParams:
    def __init__(self, cfg, param): 
        self.cfg = cfg
        self.param = param
        self.sample_size = self.get_sample_size()
        self.magmafiles = self.cfg.get('magma', 'magmafiles')
        self.magmadir = self.cfg.get('magma', 'magmadir')
        self.mhcopt = self.param.get('params', 'MHCopt')
        self.mhc = self.get_mhc()
        self.mhc_region = self.param.get('params', 'extMHC')
        self.mhc_start = self.get_mhc_region()[0]
        self.mhc_end = self.get_mhc_region()[1]
        self.refpanel = self.param.get('params', 'refpanel')
        self.pop = self.param.get('params', 'pop')
        self.ensg_v = self.param.get('params', 'ensembl')
        self.magma_window = self.param.get('magma', 'magma_window').replace(" ","")
        self.magma_exp = self.param.get('magma', 'magma_exp').split(":")
        
    def get_sample_size(self):
        sample_size = self.param.get('params', 'N')
        if sample_size == "NA":
            sample_size = None
        return sample_size
    
    def get_mhc(self):
        mhc = int(self.param.get('params', 'exMHC'))
        if mhc == 1 and self.mhcopt == "annot":
            mhc = 0
        return mhc
    
    def get_mhc_region(self):
        mhc_start = 29614758
        mhc_end = 33170276
        if self.mhc_region != "NA":
            mhc_start, mhc_end = self.mhc_region.split("-")
        return int(mhc_start), int(mhc_end)
    
def magma_gene_annotation(magma_params):
    """ run magma gene annotation if necessary"""
    if "," in magma_params.magma_window:
        w = magma_params.magma_window
    else:
        # setting up some paths
        magma_exec = os.path.join(magma_params.magmadir, "magma")
        bim = os.path.join(magma_params.magmafiles, magma_params.refpanel, magma_params.pop + ".bim")
        gene = os.path.join(magma_params.magmafiles, "ENSG" + magma_params.ensg_v + ".coding.genes.txt")
        out = os.path.join(magma_params.magmafiles, magma_params.refpanel, magma_params.pop + "_ENSG" + magma_params.ensg_v + "_w" + magma_params.magma_window.replace(",", "_"))
        
        w = magma_params.magma_window + "," + magma_params.magma_window
        
        command = magma_exec + " --annotate window=" + w + " --snp-loc " + bim + " --gene-loc " + gene + " --out " + out
        
        os.system(command)
        
def prepare_magma_file(filedir, magma_params, col_class):
    """ prepare magma.input file from input.snps file """
    
    magma_file = open(os.path.join(filedir, "magma.input"), "w")
    
    with open(os.path.join(filedir, "input.snps"), "r") as f:
        next(f)
        for line in f:
            items = line.rstrip("\n").split("\t")
            chr = items[0]
            pos = items[1]
            if magma_params.mhc == 1:
                if chr == str(6):
                    if int(pos) >= magma_params.mhc_start and int(pos) <= magma_params.mhc_end:
                        continue
            alleles = sorted([items[2], items[3]])
            if not magma_params.sample_size is None: 
                out = [chr + ":" + pos + ":" + alleles[0] + ":" + alleles[1], str(items[5]), str(magma_params.sample_size)]
            else:
                out = [chr + ":" + pos + ":" + alleles[0] + ":" + alleles[1], str(items[5]), str(items[col_class.Ncol_idx])]
            print("\t".join(out), file=magma_file)
    magma_file.close()
            
def magma_gene_analysis(filedir, magma_params):
    """ run magma gene analysis"""
    # setting up some paths
    magma_exec = os.path.join(magma_params.magmadir, "magma")
    bfile = os.path.join(magma_params.magmafiles, magma_params.refpanel, magma_params.pop )
    pval = os.path.join(filedir, "magma.input")
    gene_annot = os.path.join(magma_params.magmafiles, magma_params.refpanel, magma_params.pop + "_ENSG" + magma_params.ensg_v + "_w" + magma_params.magma_window.replace(",", "_")) + ".genes.annot"
    out = os.path.join(filedir, "magma")
    
    command = magma_exec + " --bfile " + bfile + " --pval " + pval + " N=20000" + " --gene-annot " + gene_annot + " --out " + out
    
    os.system(command)
    
def magma_gene_set(filedir, magma_params):
    """ run magma gene set analysis """
    # setting up some paths
    magma_exec = os.path.join(magma_params.magmadir, "magma")
    gene_results = os.path.join(filedir, "magma.genes.raw")
    gene_sets = os.path.join(magma_params.magmafiles, "magma_GS.txt")
    out = os.path.join(filedir, "magma")
    
    command = magma_exec + " --gene-results " + gene_results + " --set-annot " + gene_sets + " --out " + out
    os.system(command)
    
def magma_gene_property(filedir, magma_params):
    """ run magma gene property for tissues """
    
    # setting up some paths
    magma_exec = os.path.join(magma_params.magmadir, "magma")
    gene_results = os.path.join(filedir, "magma.genes.raw")
    
    for gene_exp in magma_params.magma_exp:
        gene_exp_name = gene_exp.split("/")[-1]
        gene_covar = os.path.join(magma_params.magmafiles, gene_exp + ".txt")
        out = os.path.join(filedir, "magma_exp_" + gene_exp_name)
        
        command = magma_exec + " --gene-results " + gene_results + " --gene-covar " + gene_covar + " --model direction=greater condition-hide=Average" + " --out " + out
        os.system(command)
        
def look_up_gene_symbol(filedir, magma_params):
    """ add a SYMBOL column to the file magma.genes.out"""
    genes_out = pd.read_csv(os.path.join(filedir, "magma.genes.out.tmp"), sep="\t")
    ensg = pd.read_csv("/mnt/d/standalone_fuma/ENSG/v102/ENSG.genes.txt", sep="\t", usecols=["ensembl_gene_id", "external_gene_name"]) #TODO: remove hard coding
    ensg.rename(columns={"external_gene_name": "SYMBOL", "ensembl_gene_id": "GENE"}, inplace=True)
    print(ensg.head())
    print(genes_out.head())
    genes_out_merged = pd.merge(genes_out, ensg, on="GENE")
    genes_out_merged.to_csv(os.path.join(filedir, "magma.genes.out"), sep="\t", index=False)
    
    
def magma(filedir, config_class, col_class, logger):
    
    # create a class to store parameters specific to magma analyses
    magma_params = MagmaParams(config_class.cfg, config_class.param)
    
    # check if the gene annotation exists. If not, create it using magma
    gene_annot = os.path.join(magma_params.magmafiles, magma_params.refpanel, magma_params.pop + "_ENSG" + magma_params.ensg_v + "_w" + magma_params.magma_window.replace(",", "_")) + ".genes.annot"
    if not os.path.exists(gene_annot):
        logger.info("Creating annotation file.")
        magma_gene_annotation(magma_params)
    
    # prepare magma.input file where:
    # column 1 is chr:pos:allele1:allele2
    # column 2 is p value
    # column 3 is N
    prepare_magma_file(filedir, magma_params, col_class)
    
    # creates magma.genes.raw and magma.genes.out
    magma_gene_analysis(filedir, magma_params)
    if not os.path.exists(os.path.join(filedir, "magma.genes.out")):
        sys.exit("MAGMA ERROR")
    
    #TODO: do this better
    os.system("sed 's/ \\+/\\t/g' "+filedir+"/magma.genes.out > "+filedir+"/magma.genes.out.tmp")
    
    # perform magma gene set analyses
    magma_gene_set(filedir=filedir, magma_params=magma_params)
    
    # perform magma gene property analyses
    magma_gene_property(filedir=filedir, magma_params=magma_params)
    
    # add gene symbol to the file magma.genes.out
    look_up_gene_symbol(filedir=filedir, magma_params=magma_params)
    
    # Bonferroni correction
    r_cmd = ["Rscript", os.path.join(os.path.dirname(os.path.realpath(__file__)), "bon_correction.R"), "--filedir", filedir]
    subprocess.run(r_cmd, capture_output=True, text=True)
    
    
    