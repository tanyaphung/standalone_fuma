import configparser
import click
import os
import sys
from commands.log import log

@click.command()
@click.option('--filedir', '-d', required=True, help='Path to the directory folder.'
              'The file input.gwas is expected in this folder.')

# options for [inputfiles]
@click.option('--chr_colname', default="NA",
              help='Name of chromosome column. Default is NA.')
@click.option('--pos_colname', default="NA",
              help='Name of position column. Default is NA.')
@click.option('--rsid_colname', default="NA",
              help='Name of rsID column. Default is NA.')
@click.option('--p_colname', default="NA",
              help='Name of p value column. Default is NA.')
@click.option('--ea_colname', default="NA",
              help='Name of effect allele column. Default is NA.')
@click.option('--nea_colname', default="NA",
              help='Name of non effect allele column. Default is NA.')
@click.option('--or_colname', default="NA",
              help='Name of or column. Default is NA.')
@click.option('--be_colname', default="NA",
              help='Name of be column. Default is NA.')
@click.option('--se_colname', default="NA",
              help='Name of se column. Default is NA.')
@click.option('--lead_snps_file', default="NA",
              help='File path to the lead snps file. Default is NA.')
@click.option('--add_lead_snps', default=1,
              help='Identify additional independent lead SNPs. Default is 1.')
@click.option('--regions_file', default="NA",
              help='File path to the region file. Default is NA.')

# options for [params]
@click.option('--grch38', default=0,
              help='Genome build. Default is 0 for GRCh37.' 
              'Put in 1 for GRCh38.')
@click.option('--sample_size', default="NA",
              help='Sample size. Default is NA.')
@click.option('--n_colname', default="NA",
              help='Name of the sample size column if exists.'
              'Default is NA.')
@click.option('--exmhc', default=1,
              help='Whether or not to exclude MHC. Default is 1 (exclude).')
@click.option('--mhcopt', default="annot",
              help='Which step to exclude MHC. Default is annot.')
@click.option('--extmhc', default="NA",
              help='Extended MHC region. Default is NA.')
@click.option('--ensembl', default="v102",
              help='Ensemble version. Default is v102.')
@click.option('--genetype', default="protein_coding",
              help='Gene type. Default is protein_coding.')
@click.option('--leadp', default=5e-8,
              help='P value threshold. Default is 5e-8.')
@click.option('--gwasp', default=0.05,
              help='GWAS p value threshold. Default is 0.05.')
@click.option('--r2', default=0.6,
              help='Value for r2. Default is 0.6.')
@click.option('--r2_2', default=0.1,
              help='Value for r2_2. Default is 0.1.')
@click.option('--refpanel', default="1KG/Phase3",
              help='Specify ref panel. Default is 1KG/Phase3.')
@click.option('--pop', default="EUR",
              help='Specify population. Default is EUR.')
@click.option('--maf', default=0,
              help='MAF threshold. Default is 0.')
@click.option('--refsnps', default=1,
              help='refSNPs. Default is 1.')
@click.option('--mergedist', default=250,
              help='mergeDist. Default is 250.')

# options for [magma]
@click.option('--magma', default=0,
              help='Option to run MAGMA' 
              'Default is 0. Put in 1 to run MAGMA')
@click.option('--magma_window', default="NA",
              help='Specify magma window. Default is NA.')
@click.option('--magma_exp', default="NA",
              help='Specify files to use for MAGMA gene expression analyses.'
              'Default is NA.')




def params(filedir, chr_colname, pos_colname, rsid_colname, p_colname, ea_colname, nea_colname, or_colname, be_colname, se_colname, lead_snps_file, add_lead_snps, regions_file, grch38, sample_size, n_colname, exmhc, mhcopt, extmhc, ensembl, genetype, leadp, gwasp, r2, r2_2, refpanel, pop, maf, refsnps, mergedist, magma, magma_window, magma_exp):
    """Create the required params.config file."""
    
    logger = log(filedir=filedir) #define logger
    
    if sample_size=="NA" and n_colname == "NA":
        sys.exit("Error: Either sample size (--sample_size / -n) or N col (--n_colname) needs to be provided.")
    if magma == 0:
        logger.info("You did not select to run MAGMA. If this is incorrect, please resubmit specifying --magma 1")
    if not os.path.exists(os.path.join(filedir, "input.gwas")):
        sys.exit("Error: File input.gwas is not found in the file directory specified.")
    
    config = configparser.RawConfigParser()
    config.optionxform = str

    # Add sections and key-value pairs
    config['inputfiles'] = {'gwasfile': 'input.gwas', 
                            'chrcol': chr_colname,
                            'poscol': pos_colname,
                            'rsIDcol': rsid_colname,
                            'pcol': p_colname,
                            'eacol': ea_colname,
                            'neacol': nea_colname,
                            'orcol': or_colname,
                            'becol': be_colname,
                            'secol': se_colname,
                            'leadSNPsfile': lead_snps_file,
                            'addleadSNPs': add_lead_snps,
                            'regionsfile': regions_file}
    config['params'] = {'GRCh38': grch38,
                        'N': sample_size,
                        'Ncol': n_colname,
                        'exMHC': exmhc,
                        'MHCopt': mhcopt,
                        'extMHC': extmhc,
                        'ensembl': ensembl,
                        'genetype': genetype,
                        'leadP': leadp,
                        'gwasP': gwasp,
                        'r2': r2,
                        'r2_2': r2_2,
                        'refpanel': refpanel,
                        'pop': pop,
                        'MAF': maf,
                        'refSNPs': refsnps,
                        'mergeDist': mergedist}
    
    if magma == 1 and magma_window == "NA":
        magma_window = 0
    
    if not magma_exp == "NA":
        magma_exp_all = ":".join("GTEx/v8/" + i for i in magma_exp.split(","))
    elif magma == 1 and magma_exp == "NA":
        magma_exp_all = "GTEx/v8/gtex_v8_ts_avg_log2TPM:GTEx/v8/gtex_v8_ts_general_avg_log2TPM"
    else:
        magma_exp_all = "NA"    
    config['magma'] = {'magma': magma,
                        'magma_window': magma_window,
                        'magma_exp': magma_exp_all}
    

    # Write the configuration to a file
    with open(os.path.join(filedir, "params.config"), 'w') as configfile:
        config.write(configfile)