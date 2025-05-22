import os


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
    magma_file = open(os.path.join(filedir, "magma.input"), "w")
    
    with open(os.path.join(filedir, "input.snps"), "r") as f:
        next(f)
        for line in f:
            items = line.rstrip("\n").split("\t")
            chr = items[0]
            pos = items[1]
            alleles = sorted([items[2], items[3]])
            if not magma_params.sample_size is None: 
                out = [chr + ":" + pos + ":" + alleles[0] + ":" + alleles[1], str(items[5]), str(magma_params.sample_size)]
            else:
                out = [chr + ":" + pos + ":" + alleles[0] + ":" + alleles[1], str(items[5]), str(items[col_class.Ncol_idx])]
            print("\t".join(out), file=magma_file)
            
def magma(filedir, config_class, col_class):
    magma_params = MagmaParams(config_class.cfg, config_class.param)
    
    #TODO: add checking if the annotation file exists
    # magma_gene_annotation(magma_params)
    
    #TODO: check if the gene annotation is successful
    
    prepare_magma_file(filedir, magma_params, col_class)