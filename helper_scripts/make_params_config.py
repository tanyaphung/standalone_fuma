import argparse

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Generate pseudonymonized ID.")
    
    parser.add_argument("--job_folder", default=None, help="REQUIRED. Path to the job folder.")
    
    #optional arguments
    ## inputfiles section
    parser.add_argument("--snp2genegrch38", default=None, help="OPTIONAL. Enter 1 if the input GWAS summary statistics are in GRCh38; otherwise enter 0.")
    parser.add_argument("--chrcol", default="NA", help="OPTIONAL. Column name for chromosome in the input gwas sumstat. Default: NA")
    parser.add_argument("--poscol", default="NA", help="OPTIONAL. Column name for position in the input gwas sumstat. Default: NA")
    parser.add_argument("--rsIDcol", default="NA", help="OPTIONAL. Column name for rsID in the input gwas sumstat. Default: NA")
    parser.add_argument("--pcol", default="NA", help="OPTIONAL. Column name for p-value in the input gwas sumstat. Default: NA")
    parser.add_argument("--eacol", default="NA", help="OPTIONAL. Column name for effect allele in the input gwas sumstat. Default: NA")
    parser.add_argument("--neacol", default="NA", help="OPTIONAL. Column name for non-effect allele in the input gwas sumstat. Default: NA")
    parser.add_argument("--orcol", default="NA", help="OPTIONAL. Column name for odds ratio in the input gwas sumstat. Default: NA")
    parser.add_argument("--becol", default="NA", help="OPTIONAL. Column name for beta in the input gwas sumstat. Default: NA")
    parser.add_argument("--secol", default="NA", help="OPTIONAL. Column name for standard error in the input gwas sumstat. Default: NA")
    parser.add_argument("--leadSNPsfile", default="NA", help="OPTIONAL. Put in input.lead if exists. Default: NA")
    parser.add_argument("--addleadSNPs", default="1", help="OPTIONAL. Indicate 0 if you do NOT WANT TO IDENTIFY additional lead SNPs than the provided ones. Default: 0")
    parser.add_argument("--regionsfile", default="NA", help="OPTIONAL. Put in input.regions if exists. Default: NA")
    
    ### params section
    parser.add_argument("--N", default=None, help="OPTIONAL. Sample size of the GWAS sumstat. Either N or Ncol must be provided, but not both.")
    parser.add_argument("--Ncol", default=None, help="OPTIONAL. Column name for sample size in the input gwas sumstat. Either N or Ncol must be provided, but not both.")
    parser.add_argument("--exMHC", default="1", help="OPTIONAL. Indicate 0 if you want to include the MHC region in the analysis. Default: 1")
    parser.add_argument("--MHCopt", default="annot", help="OPTIONAL. Indicate either 'annot', 'all', or 'magma' to specify how to handle the MHC region. Default: annot")
    parser.add_argument("--extMHC", default="NA", help="OPTIONAL. User defined MHC region (e.g. 25000000-33000000). When this option is not given, the default MHC region will be used. Default: NA")
    parser.add_argument("--ensembl", default="v102", help="OPTIONAL. Specify the Ensembl version to use for gene annotation. Options are v102 and v110. Default: v102")
    parser.add_argument("--genetype", default="protein_coding", help="OPTIONAL. Specify the gene type to use for gene annotation. Options are all, protein_coding, lincRNA:antisense:retained_intronic:sense_intronic:sense_overlapping:macro_lncRNA, miRNA:piRNA:rRNA:siRNA:snRNA:snoRNA:tRNA:vaultRNA, pseudogene:processed_pseudogene:unprocessed_pseudogene:polymorphic_pseudogene:IG_C_pseudogene:IG_D_pseudogene:ID_V_pseudogene:IG_J_pseudogene:TR_C_pseudogene:TR_D_pseudogene:TR_V_pseudogene:TR_J_pseudogene, IG_C_gene:TG_D_gene:TG_V_gene:IG_J_gene, or TR_C_gene:TR_D_gene:TR_V_gene:TR_J_gene . Default: protein_coding")
    parser.add_argument("--leadP", default="5e-8", help="OPTIONAL. Lead P-value threshold. Default: 5e-8")
    parser.add_argument("--gwasP", default="0.05", help="OPTIONAL. GWAS P-value threshold. Default: 0.05")
    parser.add_argument("--r2", default="0.6", help="OPTIONAL. r2 threshold to define independent significant SNPs. Default: 0.6")
    parser.add_argument("--r2_2", default="0.1", help="OPTIONAL. r2 threshold for second clumping to define lead SNPs from independent significant SNPs. Default: 0.1")
    parser.add_argument("--pop", default="EUR", help="OPTIONAL. Population to use for LD reference panel. Options are ALL, AFR, AMR, EAS, EUR, SAS. Default: EUR")
    parser.add_argument("--MAF", default="0", help="OPTIONAL. Minimum Minor Allele Frequency. This threshold defines the minimum MAF of the SNPs to be included in the annotation. MAFs are based on the selected reference panel.Default: 0")
    parser.add_argument("--refSNPs", default="1", help="OPTIONAL. Indicate 0 if you do not want to include variants in reference panel (non-GWAS tagged SNPs in LD). Default: 1")
    parser.add_argument("--mergeDist", default="250", help="OPTIONAL. Maximum distance between LD blocks to merge into a locus. Default: 250kb")
    
    ## magma section
    parser.add_argument("--magma", default="0", help="OPTIONAL. Indicate if you want to run MAGMA gene analysis. Default: 0")
    parser.add_argument("--magma_window", default="0", help="OPTIONAL. Specify the window size for MAGMA gene analysis. Default: 0. One value will set same window size both sides, two values separated by comma will set different window sizes for up- and downstream. e.g. 2,1 will set window sizes 2kb upstream and 1kb downstream of the genes.")
    parser.add_argument("--magma_exp", default="GTEx/v8/gtex_v8_ts_avg_log2TPM:GTEx/v8/gtex_v8_ts_general_avg_log2TPM", help="OPTIONAL. Specify the expression data to use for MAGMA gene analysis. Default: GTEx/v8/gtex_v8_ts_avg_log2TPM:GTEx/v8/gtex_v8_ts_general_avg_log2TPM. Choose from these options: GTEx/v8/gtex_v8_ts_avg_log2TPM, GTEx/v8/gtex_v8_ts_general_avg_log2TPM, BrainSpan/bs_age_avg_log2RPKM, or BrainSpan/bs_dev_avg_log2RPKM. If more than 1 datasets are selected, separate them with a colon (:).")
    
    ## posMap section
    parser.add_argument("--posMap", default="1", help="OPTIONAL. Indicate if you want to run positional mapping. Default: 1")
    parser.add_argument("--posMapWindowSize", default=None, help="OPTIONAL. Distance to genes to map. Default: NA.")
    parser.add_argument("--posMapAnnot", default=None, help="OPTIONAL. Functional consequences of SNPs on genes to map. Default: NA. Choose from these options: exonic, splicing, intronic, UTR3, UTR5, upstream, and downstream")
    parser.add_argument("--posMapCADDth", default="0", help="OPTIONAL. CADD score threshold for positional mapping. Default: 0")
    parser.add_argument("--posMapRDBth", default="NA", help="OPTIONAL. RegulomeDB score threshold for positional mapping. Default: NA. Choose from these options: 1a, 1b, 1c, 1d, 1e, 1f, 2a, 2b, 2c, 3a, 3b, 4, 5, 6")
    parser.add_argument("--posMapChr15", default="NA") #TODO: add later
    parser.add_argument("--posMapChr15Max", default="NA") #TODO: add later
    parser.add_argument("--posMapChr15Meth", default="NA") #TODO: add later
    parser.add_argument("--posMapAnnoDs", default="NA") #TODO: add later
    parser.add_argument("--posMapAnnoMeth", default="NA") #TODO: add later
    
    

    parsed_args = parser.parse_args(args)

    if parsed_args.job_folder is None:
        parsed_args.job_folder = input("Enter the job folder path: ").strip()
    if not parsed_args.job_folder:
        parser.error("--job_folder is required.")

    if parsed_args.snp2genegrch38 is None:
        response = input("Is the input GWAS summary statistics in GRCh38? Enter 1 for yes or 0 for no [0]: ").strip()
        parsed_args.snp2genegrch38 = response if response else "0"

    parsed_args.snp2genegrch38 = str(parsed_args.snp2genegrch38).strip()
    if parsed_args.snp2genegrch38 not in {"0", "1"}:
        parser.error("--snp2genegrch38 must be either 0 or 1.")

    if parsed_args.N is None and parsed_args.Ncol is None:
        response = input("Provide sample size using --N or a sample-size column name using --Ncol. Enter either a numeric N or a column name: ").strip()
        if response:
            try:
                int(response)
                parsed_args.N = response
                parsed_args.Ncol = "NA"
            except ValueError:
                parsed_args.N = "NA"
                parsed_args.Ncol = response

    parsed_args.N = str(parsed_args.N).strip() if parsed_args.N is not None else "NA"
    parsed_args.Ncol = str(parsed_args.Ncol).strip() if parsed_args.Ncol is not None else "NA"

    if parsed_args.posMapWindowSize is None and parsed_args.posMapAnnot is None:
        response = input("Provide either a positional map window size or a positional annotation. Enter a number for window size or a value such as exonic/upstream: ").strip()
        if response:
            try:
                float(response)
                parsed_args.posMapWindowSize = response
                parsed_args.posMapAnnot = "NA"
            except ValueError:
                parsed_args.posMapWindowSize = "NA"
                parsed_args.posMapAnnot = response

    parsed_args.posMapWindowSize = str(parsed_args.posMapWindowSize).strip() if parsed_args.posMapWindowSize is not None else "NA"
    parsed_args.posMapAnnot = str(parsed_args.posMapAnnot).strip() if parsed_args.posMapAnnot is not None else "NA"

    if parsed_args.N.upper() == "NA" and parsed_args.Ncol.upper() == "NA":
        parser.error("Either --N or --Ncol must be supplied, but not both.")
    if parsed_args.N.upper() != "NA" and parsed_args.Ncol.upper() != "NA":
        parser.error("Specify either --N or --Ncol, not both.")
    if parsed_args.Ncol.upper() != "NA":
        if not parsed_args.Ncol or parsed_args.Ncol.upper() == "NA":
            parser.error("--Ncol cannot be empty.")
        parsed_args.Ncol = parsed_args.Ncol.upper()
    else:
        parsed_args.Ncol = "NA"

    if parsed_args.posMapWindowSize.upper() == "NA" and parsed_args.posMapAnnot.upper() == "NA":
        parser.error("Either --posMapWindowSize or --posMapAnnot must be supplied, but not both else one must be NA.")
    if parsed_args.posMapWindowSize.upper() != "NA" and parsed_args.posMapAnnot.upper() != "NA":
        parser.error("Specify either --posMapWindowSize or --posMapAnnot, not both.")
    if parsed_args.posMapWindowSize.upper() != "NA":
        parsed_args.posMapWindowSize = parsed_args.posMapWindowSize
    else:
        parsed_args.posMapWindowSize = "NA"
    if parsed_args.posMapAnnot.upper() != "NA":
        parsed_args.posMapAnnot = parsed_args.posMapAnnot
    else:
        parsed_args.posMapAnnot = "NA"

    return parsed_args

def main():
    args = parse_args()
    job_folder = args.job_folder
    
    with open(f"{job_folder}/params.config", "w") as f:
        ## inputfiles section
        f.write("[inputfiles]\n")
        f.write(f"snp2genegrch38 = {args.snp2genegrch38}\n")
        f.write(f"chrcol = {args.chrcol}\n")
        f.write(f"poscol = {args.poscol}\n")
        f.write(f"rsIDcol = {args.rsIDcol}\n")
        f.write(f"pcol = {args.pcol}\n")
        f.write(f"eacol = {args.eacol}\n")
        f.write(f"neacol = {args.neacol}\n")
        f.write(f"orcol = {args.orcol}\n")
        f.write(f"becol = {args.becol}\n")
        f.write(f"secol = {args.secol}\n")
        f.write(f"leadSNPsfile = {args.leadSNPsfile}\n")
        f.write(f"addleadSNPs = {args.addleadSNPs}\n")
        f.write(f"regionsfile = {args.regionsfile}\n")
        
        ## params section
        f.write("\n[params]\n")
        f.write(f"N = {args.N}\n")
        f.write(f"Ncol = {args.Ncol}\n")
        f.write(f"exMHC = {args.exMHC}\n")
        f.write(f"MHCopt = {args.MHCopt}\n")
        f.write(f"extMHC = {args.extMHC}\n")
        f.write(f"ensembl = {args.ensembl}\n")
        f.write(f"genetype = {args.genetype}\n")
        f.write(f"leadP = {args.leadP}\n")
        f.write(f"gwasP = {args.gwasP}\n")
        f.write(f"r2 = {args.r2}\n")
        f.write(f"r2_2 = {args.r2_2}\n")
        f.write(f"refpanel = 1KG/Phase3\n")
        f.write(f"pop = {args.pop}\n")
        f.write(f"MAF = {args.MAF}\n")
        f.write(f"refSNPs = {args.refSNPs}\n")
        f.write(f"mergeDist = {args.mergeDist}\n")
        
        ## magma section
        f.write("\n[magma]\n")
        f.write(f"magma = {args.magma}\n")
        f.write(f"magma_window = {args.magma_window}\n")
        f.write(f"magma_exp = {args.magma_exp}\n")
        
        ## posMap section
        f.write("\n[posMap]\n")
        f.write(f"posMap = {args.posMap}\n")
        f.write(f"posMapWindowSize = {args.posMapWindowSize}\n")
        f.write(f"posMapAnnot = {args.posMapAnnot}\n")
        f.write(f"posMapCADDth = {args.posMapCADDth}\n")
        f.write(f"posMapRDBth = {args.posMapRDBth}\n")
        f.write(f"posMapChr15 = {args.posMapChr15}\n")
        f.write(f"posMapChr15Max = {args.posMapChr15Max}\n")
        f.write(f"posMapChr15Meth = {args.posMapChr15Meth}\n")
        f.write(f"posMapAnnoDs = {args.posMapAnnoDs}\n")
        f.write(f"posMapAnnoMeth = {args.posMapAnnoMeth}\n")

if __name__ == "__main__":
    main()