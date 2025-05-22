import re
import csv
import sys
import pandas as pd
import logging
import os
import numpy as np
from collections import defaultdict
import configparser
import tabix
from bisect import bisect_left
import subprocess
# from profilehooks import profile

logger = logging.getLogger()

def check_inputs(arguments):
    if len(arguments) < 2:
        sys.exit('ERROR: not enough arguments\nUSAGE ./gwas_file.py <filedir>')
    if not os.path.isdir(arguments[1]):
        sys.exit(f"The job directory {arguments[1]} does not exist on the FUMA server.")

def log(filedir): 
    # create and configure logger
    logging.basicConfig(filename=os.path.join(filedir, "snp2gene_logs.txt"),
                        format='%(asctime)s %(levelname)s %(message)s',
                        filemode='w')
    
    logger = logging.getLogger() #creating an object
    logger.setLevel(logging.DEBUG) #setting the threshold of logger to DEBUG
    return logger

def initialize_outfile(filedir, col_class): 
    out = open(os.path.join(filedir, "input.snps"), "w")
    out.write('\t'.join(["chr", "bp", "non_effect_allele", "effect_allele", "rsID", "p"] + write_col_header(col_class.orcol_idx, col_class.becol_idx, col_class.secol_idx, col_class.Ncol_idx)) + "\n") #format header
    return out

class Configuration: 
    def __init__(self, filedir): 
        self.filedir = filedir
        self.ctg = None
        self.param = None
        self.get_ctg()
        self.get_param()
        
    def get_ctg(self): 
        app_config_fp = os.path.join(os.path.dirname(os.path.realpath(__file__)), "app.config")
        if os.path.exists(app_config_fp):
            self.cfg = configparser.ConfigParser()
            self.cfg.read(app_config_fp)
        else:
            logger.error("app.config file is not found in the script folder")
            sys.exit("app.config file is not found in the script folder")
    
    def get_param(self):
        params_config_fp = os.path.join(self.filedir, 'params.config')
        if os.path.exists(params_config_fp):
            self.param = configparser.ConfigParser()
            self.param.read(os.path.join(self.filedir, 'params.config'))
        else: 
            logger.error("params.config file is not found in the job folder")
            sys.exit("params.config file is not found in the job folder")

class HeaderInfo: 
    def __init__(self, filedir):
        self.filedir = filedir
        self._header_line = None
        self._delim = None
        self._header = None
        self._nheader = None
        self.get_header()
        self.get_delimiter()
        self.format_header()
    
    def get_header(self):
        """ Assign the first line that does not start with # to be the header"""
        with (open(os.path.join(self.filedir, "input.gwas"), 'r')) as fin:
            self._header_line = fin.readline()
            while re.match("^#", self._header_line):
                self._header_line = fin.readline()
    
    def get_delimiter(self):
        if self._header_line:
            if re.match(r'.*\s\s.*', self._header_line) is not None:
                self._delim = r'\s+'
            sniffer = csv.Sniffer()
            dialect = sniffer.sniff(self._header_line)
            self._delim = dialect.delimiter
        
    def format_header(self):
        if self._header_line:
            header = self._header_line.rstrip("\n").split(self._delim)
            self._header = [i.upper() for i in header]
            self._nheader = len(self._header)
        
        
class ColIndex: 
    def __init__(self, chrcol_idx, poscol_idx, neacol_idx, eacol_idx, rsIDcol_idx, pcol_idx, orcol_idx, becol_idx, secol_idx, Ncol_idx):
        self.chrcol_idx = chrcol_idx
        self.poscol_idx = poscol_idx
        self.neacol_idx = neacol_idx
        self.eacol_idx = eacol_idx
        self.rsIDcol_idx = rsIDcol_idx
        self.pcol_idx = pcol_idx
        self.orcol_idx = orcol_idx
        self.becol_idx = becol_idx
        self.secol_idx = secol_idx
        self.Ncol_idx = Ncol_idx


def float_safe(s):
    """
    Return the float of a string if true, otherwise return empty
    """
    try:
        if "e" in s and float(s) < 1e-308:
            return float(1e-308)
        return float(s)
    except ValueError:
        return


def check_format(file, filetype): 
    """This function checks format of pre-defined files that users can submit.
    Currently it can be the lead SNPs file and the genomic regions files."""
    if os.stat(file).st_size == 0:
        logger.error(f"You uploaded a {filetype} file but this file is empty.")
        return "Error_empty"
    data = pd.read_csv(file, delim_whitespace=True)
    if data.shape[0] == 0: 
        logger.error(f"You uploaded a {filetype} file but this file only has column name and no data.")
        return "Error_nodata"
    if data.shape[1] < 3:
        logger.error(f"You uploaded a {filetype} file but this file has less than 3 columns.")
        return "Error_notenoughcol"


def write_col_header(orcol_idx, becol_idx, secol_idx, Ncol_idx):
    """
    Return a list of or, beta, se, and N if any of these are present
    This function is used for formatting the header
    """
    out_list = []
    if str(orcol_idx):
        out_list.append("or")
    if str(becol_idx):
        out_list.append("beta")
    if str(secol_idx):
        out_list.append("se")
    if str(Ncol_idx):
        out_list.append("N")
    return out_list
		
  
def write_col(l, orcol_idx, becol_idx, secol_idx, Ncol_idx):
    """
    Return the list of values for orcol, becol, secol, and Ncol if any of these exists
    """
    out_list = []
    if str(orcol_idx):
        out_list.append(l[orcol_idx])
    if str(becol_idx):
        out_list.append(l[becol_idx])
    if str(secol_idx):
        out_list.append(l[secol_idx])
    if str(Ncol_idx):
        out_list.append(l[Ncol_idx])
    return out_list


def get_user_params(param):
    """returns a dictionary where:
    key: user defined value
    value: name of the column such as chrcol, rsIDcol, poscol, eacol, neacol, etc..."""
    out_dict = {}
    columns = ["chrcol", "poscol", "neacol", "eacol", "rsIDcol", "pcol", "orcol", "becol",
               "secol"]  # these are the columns under [inputfiles] while Ncol is under [params]
    for i in columns:
        out_dict[i] = param.get('inputfiles', i).upper()
    out_dict["Ncol"] = param.get('params', 'Ncol').upper()
    out_dict_inverse = {v: k for (k, v) in out_dict.items() if v != "NA"}
    return out_dict_inverse


def resolve_default_names(col_name):
    defaults = {
        "CHR": "chrcol",
        "CHROMOSOME": "chrcol",
        "CHROM": "chrcol",
        "SNP": "rsIDcol",
        "MARKERNAME": "rsIDcol",
        "RSID": "rsIDcol",
        "SNPID": "rsIDcol",
        "BP": "poscol",
        "POS": "poscol",
        "POSITION": "poscol",
        "A1": "eacol",
        "EFFECT_ALLELE": "eacol",
        "ALLELE1": "eacol",
        "ALLELEB": "eacol",
        "A2": "neacol",
        "NON_EFFECT_ALLELE": "neacol",
        "ALLELE2": "neacol",
        "ALLELEA": "neacol",
        "P": "pcol",
        "PVAL": "pcol",
        "PVALUE": "pcol",
        "P-VALUE": "pcol",
        "P_VALUE": "pcol",
        "P.VALUE": "pcol",
        "OR": "orcol",
        "BETA": "becol",
        "BE": "becol",
        "SE": "secol",
        "N": "Ncol"
    }
    return defaults.get(col_name, "NA")


def input_sanitation_helper(line, header_class, col_class, rmSNPs):
    """
    Performs input sanitation per line. The line is from the input gwas.
    rmSNPS is the file storing the lines that are removed.
    """
    if line.startswith("#"):
        rmSNPs.write(line.rstrip("\n") + " DROPPED BY FUMA: because line starts with #.\n")
        return None
    
    if line==header_class._header_line:
        return None
    
    items = line.strip().split(header_class._delim)
    if len(items) != header_class._nheader:
        rmSNPs.write(line.rstrip("\n") + " DROPPED BY FUMA: because the number of columns in this line is not equal to the number of columns in the header.\n")
        return None

    # sanitize p 
    line_pval = float_safe(items[col_class.pcol_idx])
    if not line_pval:
        rmSNPs.write(line.rstrip("\n") + " DROPPED BY FUMA: because of string in the p value column. Number is expected.\n")
        return None
    if line_pval < 0 or line_pval > 1:
        rmSNPs.write(line.rstrip("\n") + " DROPPED BY FUMA: because p value out of range (<0 or >1).\n")
        return None
    if line_pval == 0 and re.match("^0", items[col_class.pcol_idx]):
        rmSNPs.write(line.rstrip("\n") + " DROPPED BY FUMA: because P value is 0.\n")
        return None
    
    items[col_class.pcol_idx] = str(line_pval)
    
    # sanitize chromosome 
    if str(col_class.chrcol_idx):
        line_chr = items[col_class.chrcol_idx].replace("chr", "").replace("CHR", "") #strip chr
        if re.match("x", line_chr, re.IGNORECASE): #convert x to 23
            items[col_class.chrcol_idx] = '23'
        if not line_chr.isdigit(): #check if chrcol is a digit
            rmSNPs.write(
            line.rstrip("\n") + " DROPPED BY FUMA: because value for the chromosome column is not a digit.\n")
            return None
        if int(line_chr) not in range(1, 24): #checking for 1-23
            rmSNPs.write(
            line.rstrip("\n") + " DROPPED BY FUMA: because value for the chromosome column is not between 1 and 23.\n")
            return None
    
    # sanitize position
    if str(col_class.poscol_idx):
        if not items[col_class.poscol_idx].isdigit():
            rmSNPs.write(
                    line.rstrip("\n") + " DROPPED BY FUMA: because value for the position column is not a digit.\n")
            return None
    
    # sanitize rsID
    if str(col_class.rsIDcol_idx):
        if not items[col_class.rsIDcol_idx].startswith("rs"):
            rmSNPs.write(line.rstrip("\n") + " DROPPED BY FUMA: because in the rsID column, expect to start with rs for rsID format.\n")
            return None
    return items

def database_lookup(items, snpd, col_class):
    """This function replaces Tabix function. It looks up the rsID based on uid from the database. If it's not found then it returns uid"""
    uid = ":".join(
        [items[col_class.chrcol_idx], items[col_class.poscol_idx]] + sorted([items[col_class.neacol_idx].upper(), items[col_class.eacol_idx].upper()]))
    id_val = snpd.get(uid)
    if not snpd.get(uid): 
        logger.info(f"chr:pos:nea:ea for this row does not have the corresponding rsID from the database used by FUMA. Use UID instead: {items}")
        id_val = uid
    updated_out = [items[col_class.chrcol_idx], items[col_class.poscol_idx], items[col_class.neacol_idx].upper(), items[col_class.eacol_idx].upper(), id_val, items[col_class.pcol_idx]] + write_col(items, col_class.orcol_idx, col_class.becol_idx, col_class.secol_idx, col_class.Ncol_idx)
    return updated_out


def run_database_lookup(tb, chrom, start, end, snps, col_class, out):
    """ This function runs database look up """
    refSNP = np.array(list(tb.querys(str(chrom) + ":" + str(start) + "-" + str(end))))
    snps = np.array(snps)
    snpd = {refSNP[i, 2]: refSNP[i, 3] for i in range(refSNP.shape[0])}
    rsid_lookup = list(map(lambda l: database_lookup(items=l, snpd=snpd, col_class=col_class), snps))
    out.write('\n'.join(['\t'.join(s) for s in rsid_lookup]) + '\n')
    
def allele_lookup(row, snpd, col_class, rmSNPs):
    """
    This function looks up the allele
    """
    out = None
    chrom_pos = ":".join([row[col_class.chrcol_idx], row[col_class.poscol_idx]])
    rsID_data = {}
    if str(col_class.rsIDcol_idx): 
        rsID_data = {"rsID": str(row[col_class.rsIDcol_idx])}
    
    if not chrom_pos in snpd: 
        to_print = np.append(row, "DROPPED BY FUMA: because the user's chr:pos is not found in the database used by FUMA.\n")
        to_print.tofile(rmSNPs, sep='\t', format='%s')
        return None
    ref_snps = snpd[chrom_pos]
    if str(col_class.eacol_idx):
        found_nea = False
        for rsID_A1_A2 in ref_snps: #in the case that a chrom_pos is multiallelic
            A1_A2 = [rsID_A1_A2[1].upper(), rsID_A1_A2[2].upper()]
            if row[col_class.eacol_idx] in A1_A2:
                out = [str(row[col_class.chrcol_idx]), str(row[col_class.poscol_idx]), str(list(set(A1_A2).difference([row[col_class.eacol_idx]]))[0].upper()), row[col_class.eacol_idx], rsID_data.get("rsID", str(rsID_A1_A2[0])), str(row[col_class.pcol_idx])]
                found_nea = True
                break
        if not found_nea:
            to_print = np.append(row, "DROPPED BY FUMA: because the user's effect allele is not the same as in the database used by FUMA.\n")
            to_print.tofile(rmSNPs, sep='\t', format='%s')
    else: 
        if len(ref_snps) == 1:
            out = [str(row[col_class.chrcol_idx]), str(row[col_class.poscol_idx]), str(ref_snps[0][1].upper()), str(ref_snps[0][2].upper()), rsID_data.get("rsID", str(ref_snps[0][0])), str(row[col_class.pcol_idx])] 
        else:
            if rsID_data.get("rsID"):
                found_rsID_match = False
                for ref_snp in ref_snps:
                    if rsID_data.get("rsID") == ref_snp[0]:
                        out = [str(row[col_class.chrcol_idx]), str(row[col_class.poscol_idx]), str(ref_snp[1].upper()), str(ref_snp[2].upper()), rsID_data.get("rsID"), str(row[col_class.pcol_idx])]
                        found_rsID_match = True
                        break
                if not found_rsID_match:
                    to_print = np.append(row, "DROPPED BY FUMA: this variant is multiallelic but none of the rsIDs for this variant in the database matches with the rsID in the user's input summary statistics.\n")
                    to_print.tofile(rmSNPs, sep='\t', format='%s')
            else:
                to_print = np.append(row, "DROPPED BY FUMA: this variant is multiallelic but rsID does not exist in the user's input summary statistics.\n")
                to_print.tofile(rmSNPs, sep='\t', format='%s')
                
    return out

# @profile
def run_allele_lookup(tb, chrom, start, end, snps, col_class, out, rmSNPs):
    """ This function runs allele look up """
    refSNP = np.array(list(tb.querys(str(chrom) + ":" + str(start) + "-" + str(end))))
    snps = np.array(snps)
    snpd = defaultdict(list) #key is chr:pos; value is a list; to account for the same chr:pos having multiple entries
    for i in refSNP: 
        key = i[0] + ":" + i[1] #chr:pos
        value = (i[2], i[3], i[4]) #i[2] is rsID, i[3] is A1, i[4] is A2
        snpd[key].append(value)
    results = list(map(lambda l: allele_lookup(row=l, snpd=snpd, col_class=col_class, rmSNPs=rmSNPs), snps))
    if results != [None]:
        out.write('\n'.join(['\t'.join(s) for s in results if s]) + '\n')
    
    
def check_input_gwas_size(filedir):
    """ check to make sure that the input gwas summary statistics exits and is not empty """
    try: 
        input_gwas_size = os.path.getsize(os.path.join(filedir, "input.gwas"))
        if input_gwas_size == 0: 
            logger.error("The input gwas summary statistics file is empty.")
            sys.exit("The input gwas summary statistics file is empty.")
    except FileNotFoundError as e:
        logger.error("The input gwas summary statistics file is not found. Something went wrong.")
        sys.exit("The input gwas summary statistics file is not found. Something went wrong.")
        
        
def check_header(header_class):
    """ check to make sure that the header is not none """
    if header_class._header == None: 
        logger.error("No header in the input summary statistics is found.")
        sys.exit("No header in the input summary statistics is found.")
    else: 
        logger.info(f"Header detected: {header_class._header}")
    if header_class._delim == None:
        logger.error("Delimiter could not be resolved.")
        sys.exit("Delimiter could not be resolved.")
    if header_class._nheader == None:
        logger.error("Error because number of columns in the header is None.")
        sys.exit("Error because number of columns in the header is None.")            
        
def detect_bad_user_defined_col(user_defined_col, header_class):
    """ detect when a user supplies a column name but it does not exist in the header of the input file """
    bad_user_defined_col = list(set(user_defined_col.keys()).difference(header_class._header))
    if bad_user_defined_col:
        bad_user_defined_col_str = ", ".join(bad_user_defined_col)
        logger.error(f"You specified that the following are the column names but these names do not exist in the header of the file: {bad_user_defined_col_str}")
        sys.exit("ERROR: The following header(s) was not detected in your input file: " + bad_user_defined_col_str)
        
        
def check_allele(col_idx_dict):
    """ If the file contains 1 column for the allele, force this column to be the effect allele column"""
    eacol_idx = col_idx_dict.get("eacol", "")
    neacol_idx = col_idx_dict.get("neacol", "")
    if neacol_idx and not eacol_idx:
        col_idx_dict["eacol"] = neacol_idx
        col_idx_dict["neacol"] = ""
        logger.info("The input file contains only one column for allele, FUMA will proceed with assumption that this is the effect allele.")
    return col_idx_dict


def check_mandatory_cols(col_class): 
    # mandatory header check
    if not str(col_class.pcol_idx): #convert to str because python interprets bool(0) as false
        logger.error("P-value column was not found. A column for P-value is mandatory.")
        sys.exit("ERROR: P-value column was not found. A column for P-value is mandatory.")
    if not str(col_class.chrcol_idx) and not str(col_class.rsIDcol_idx):
        logger.error("Both chromosome and rsID column was not found. FUMA requires chr and pos OR rsID.")
        sys.exit("ERROR: Both chromosome and rsID column was not found. FUMA requires chr and pos OR rsID.")
    if not str(col_class.poscol_idx) and not str(col_class.rsIDcol_idx):
        logger.error("Both position and rsID column was not found. FUMA requires chr and pos OR rsID.")
        sys.exit("ERROR: Both position and rsID column was not found. FUMA requires chr and pos OR rsID.")
    
def log_user_specified_columns(user_defined_col):
    logger.info(f"User specified the following as columns:")
    
    user_friendly_colnames = {"chrcol": "Chromosome",
                              "poscol": "Position",
                              "rsIDcol": "rsID",
                              "pcol": "P-value",
                              "eacol": "Effect allele",
                              "neacol": "Non effect allele",
                              "orcol": "OR",
                              "becol": "Beta",
                              "secol": "SE",
                              "Ncol": "Ncol"
                              }
    for key in user_defined_col:
        logger.info(f"For {user_friendly_colnames[user_defined_col[key]]}, user specified column name: {key}")
    
def log_automatic_columns_detection(user_defined_col, col_idx_dict):
    logger.info(f"FUMA automatically detects the following as columns:")
    
    user_friendly_colnames = {"chrcol": "Chromosome",
                              "poscol": "Position",
                              "rsIDcol": "rsID",
                              "pcol": "P-value",
                              "eacol": "Effect allele",
                              "neacol": "Non effect allele",
                              "orcol": "OR",
                              "becol": "Beta",
                              "secol": "SE"
                              }
    
    auto_col = list(set(col_idx_dict.keys()).difference(user_defined_col.values()))
    for i in auto_col:
        logger.info(f"For {user_friendly_colnames[i]}, FUMA automatically detects column index: {col_idx_dict[i]}")
    
    
def validate(filedir):
    # create a config class
    config_class = Configuration(filedir=filedir)
    check_input_gwas_size(filedir=filedir)

    # create a header class
    header_class = HeaderInfo(filedir=filedir)
    check_header(header_class=header_class)

    # get user defined columns
    user_defined_col = get_user_params(
        config_class.param)  #this returns a dict where key is the user defined value and value is the name of the col such as chrcol, rsIDcol, poscol, eacol, neacol, etc...
    
    log_user_specified_columns(user_defined_col=user_defined_col)

    detect_bad_user_defined_col(user_defined_col=user_defined_col, header_class=header_class)

    # create a col class
    col_idx_dict_tmp = {user_defined_col.get(col, resolve_default_names(col)): idx for (idx, col) in enumerate(header_class._header)} #resolve column names: col is the value in the header. If it's in defined by the user, it will be obtained by the command user_defined_col.get(col). Otherwise, we will try to resolve it by the default names. This is a dictionary where key is the type of column (chrcol, poscol, etc...) and value is the index in the header
    col_idx_dict_tmp2 = {k: col_idx_dict_tmp[k] for k in col_idx_dict_tmp if k != "NA"} #remove NA. 
    # example: {'rsIDcol': 0, 'chrcol': 1, 'poscol': 2, 'eacol': 3, 'orcol': 6, 'secol': 7, 'pcol': 11, 'neacol': 12}
    
    # allele column check: if only one allele is defined, this has to be alt (effect) allele
    col_idx_dict = check_allele(col_idx_dict=col_idx_dict_tmp2)
    
    # print to log the columns that were inferred
    log_automatic_columns_detection(user_defined_col=user_defined_col, col_idx_dict=col_idx_dict)
    
    # create a class to store the index
    col_class = ColIndex(chrcol_idx=col_idx_dict.get("chrcol", ""), poscol_idx=col_idx_dict.get("poscol", ""), neacol_idx=col_idx_dict.get("neacol", ""), eacol_idx=col_idx_dict.get("eacol", ""), rsIDcol_idx=col_idx_dict.get("rsIDcol", ""), pcol_idx=col_idx_dict.get("pcol", ""), orcol_idx=col_idx_dict.get("orcol", ""), becol_idx=col_idx_dict.get("becol", ""), secol_idx=col_idx_dict.get("secol", ""), Ncol_idx=col_idx_dict.get("Ncol", ""))

    check_mandatory_cols(col_class=col_class)
        
    return config_class, header_class, col_class


def sanitize(filedir, header_class, col_class):
    """ sanitize inputs """
    gwas_sanitized = open(os.path.join(filedir, "input.gwas.sanitized"), "w")
    rmSNPs = open(os.path.join(filedir, "input.snps.rm"), "w")
    with open(os.path.join(filedir, "input.gwas"), "r") as gwasIn:
        new_line = list(map(lambda l: input_sanitation_helper(line=l, header_class = header_class, col_class=col_class, rmSNPs=rmSNPs), gwasIn))
        for i in new_line:
            if i:
                print("\t".join(i), file=gwas_sanitized)
    gwas_sanitized.close()
    rmSNPs.close()
    

def check_sanitized_gwas_size(filedir):
    """ check if there are any snps in the sanitized file """
    try: 
        file_size = os.path.getsize(os.path.join(filedir, "input.gwas.sanitized"))
        if file_size == 0: 
            logger.error("After sanitizing the input gwas summary statistics, there are no variants left.")
            sys.exit("After sanitizing the input gwas summary statistics, there are no variants left.")
    except FileNotFoundError as e:
        logger.error("The sanitized gwas summary statitstics file is not found. Something went wrong.")
        sys.exit("The sanitized gwas summary statitstics file is not found. Something went wrong.")
        
def sort(filedir, col_class):
    """ If the user's input summary statistics has chr and pos, sort by chr and pos"""
    if str(col_class.chrcol_idx) and str(col_class.poscol_idx): 
        gwas_sanitized_fp = os.path.join(filedir, "input.gwas.sanitized")
        tempfile = os.path.join(filedir, "temp.txt")
        sort_command = "sort -k " + str(col_class.chrcol_idx + 1) + "n -k " + str(col_class.poscol_idx + 1) + "n " + gwas_sanitized_fp +" > " + tempfile
        os.system(sort_command)
        os.system("mv "+tempfile+" "+gwas_sanitized_fp)
                    
def process_all_columns_present(filedir, col_class, out): 
    """ Process input gwas summary statistics where all columns are present """
    logger.info("Chromosome, position, rsID, effect allele, and noneffect allele columns are all present in the input gwas summary statistics.")
    with open(os.path.join(filedir, "input.gwas.sanitized"), "r") as f:
        for line in f:
            s = line.rstrip("\n").split("\t")
            updated_line =  [s[col_class.chrcol_idx], s[col_class.poscol_idx], s[col_class.neacol_idx].upper(), s[col_class.eacol_idx].upper(), s[col_class.rsIDcol_idx], s[col_class.pcol_idx]] + write_col(s, col_class.orcol_idx, col_class.becol_idx, col_class.secol_idx, col_class.Ncol_idx)
            out.write('\t'.join(updated_line) + "\n")
    out.close()
    
    
def process_missing_rsID(filedir, col_class, refpanel, pop, out):
    logger.info("Both column for effect allele and non-effect allele are present. Need to look up rsID from database.")
    with open(os.path.join(filedir, "input.gwas.sanitized"), "r") as gwasIn:
        cur_chr = 1
        start = 0
        end = 0
        cur_snps = []
        for line in gwasIn:
            s = line.rstrip("\n").split("\t")
            if int(s[col_class.chrcol_idx]) == cur_chr:  # if it's the current chromosome
                if start == 0:
                    start = int(s[col_class.poscol_idx])
                if int(s[col_class.poscol_idx]) - start <= 1000000 and int(s[col_class.poscol_idx]) - start >= 0:
                    end = int(s[col_class.poscol_idx])
                    cur_snps.append(s)
                else:
                    tb = tabix.open(os.path.join(refpanel, pop, pop + ".chr" + str(cur_chr) + ".rsID.gz"))
                    run_database_lookup(tb=tb, chrom=cur_chr, start=start, end=end, snps=cur_snps, col_class=col_class, out=out)
                    start = int(s[col_class.poscol_idx])
                    end = int(s[col_class.poscol_idx])
                    cur_snps = [s]
            else:
                if start != 0 and end != 0:
                    tb = tabix.open(os.path.join(refpanel, pop, pop + ".chr" + str(cur_chr) + ".rsID.gz"))
                    run_database_lookup(tb=tb, chrom=cur_chr, start=start, end=end, snps=cur_snps, col_class=col_class, out=out)
                    cur_chr = int(s[col_class.chrcol_idx])
                    start = int(s[col_class.poscol_idx])
                    end = int(s[col_class.poscol_idx])
                    cur_snps = [s]
        if len(cur_snps) > 0:
            tb = tabix.open(os.path.join(refpanel, pop, pop + ".chr" + str(cur_chr) + ".rsID.gz"))
            run_database_lookup(tb=tb, chrom=cur_chr, start=start, end=end, snps=cur_snps, col_class=col_class, out=out)
    out.close()
            
def process_missing_allele(filedir, col_class, refpanel, pop, out):
    rmSNPs = open(os.path.join(filedir, "input.snps.rm"), "a+")
    with open(os.path.join(filedir, "input.gwas.sanitized"), "r") as gwasIn:
        cur_chr = 1
        start = 0
        end = 0
        cur_snps = []
        for line in gwasIn:
            row = line.rstrip("\n").split("\t")
            if int(row[col_class.chrcol_idx]) == cur_chr:  # if it's the current chromosome
                # print("INFO: Processing: ", cur_chr)
                if start == 0: #if it is the first row
                    # print("INFO: This is the first line, updating the start position")
                    start = int(row[col_class.poscol_idx])
                if int(row[col_class.poscol_idx]) - start <= 1000000 and int(row[col_class.poscol_idx]) - start >= 0:
                    end = int(row[col_class.poscol_idx])
                    # print("INFO: Updating maxpos. Go to the next line.")
                    cur_snps.append(row)
                else:
                    # print("INFO: Finishing getting the boundaries for 1MB. Extract alleles now.", cur_chr, start, end)
                    tb = tabix.open(os.path.join(refpanel, pop, pop + ".chr" + str(cur_chr) + ".frq.gz"))
                    run_allele_lookup(tb=tb, chrom=cur_chr, start=start, end=end, snps=cur_snps, col_class=col_class, out=out, rmSNPs=rmSNPs)
                    # initialize the next set
                    print("Initilize the next 1MB.")
                    start = int(row[col_class.poscol_idx])
                    end = int(row[col_class.poscol_idx])
                    cur_snps = [row]
            else:
                # print("INFO: Chromosome is changing. Extract allele & reset.", cur_chr, start, end)
                if start != 0 and end != 0:
                    tb = tabix.open(os.path.join(refpanel, pop, pop + ".chr" + str(cur_chr) + ".frq.gz"))
                    run_allele_lookup(tb=tb, chrom=cur_chr, start=start, end=end, snps=cur_snps, col_class=col_class, out=out, rmSNPs=rmSNPs)
                    cur_chr = int(row[col_class.chrcol_idx])
                    start = int(row[col_class.poscol_idx])
                    end = int(row[col_class.poscol_idx])
                    cur_snps = [row]
        if len(cur_snps) > 0:
            # print("INFO: Run the left over.", cur_chr, start, end)
            tb = tabix.open(os.path.join(refpanel, pop, pop + ".chr" + str(cur_chr) + ".frq.gz"))
            run_allele_lookup(tb=tb, chrom=cur_chr, start=start, end=end, snps=cur_snps, col_class=col_class, out=out, rmSNPs=rmSNPs)
    out.close()
    
def read_gwas_df(filedir, col_class, config_class):
    """ Read in the input gwas summary statistics as a pandas dataframe. Update arhived rsID if needed"""
    gwas = pd.read_csv(os.path.join(filedir, "input.gwas.sanitized"), sep="\t", header=None)
    gwas = gwas.to_numpy()
    gwas = gwas[gwas[:, col_class.rsIDcol_idx].argsort()]

    # check if the rsID from the input summary statistics has the archived ID. If so, it convert to the new ID. 
    rsID = list(gwas[:, col_class.rsIDcol_idx])
    rsIDs = set(rsID)
    
    with open(os.path.join(config_class.cfg.get('data', 'dbSNP'), "RsMerge146.txt"), 'r') as rsID_dbSNP: #TODO: make the file name RsMerge146.txt into a variable
        for l in rsID_dbSNP:
            l = l.rstrip("\n").split()
            if l[0] in rsIDs: #if the archived rsID is found in the input summary statistics
                logger.info(f"Archived rsID {l[0]} found. Updating the input summary statistics with the rsID from dbSNP.")
                gwas[bisect_left(rsID, l[0]), col_class.rsIDcol_idx] = l[1] #update the archived ID with the updated IDs

    # sort input snps by rsID for index look up
    gwas = gwas[gwas[:, col_class.rsIDcol_idx].argsort()]
    return gwas

    
def process_missing_chr_pos(filedir, col_class, config_class, out):
    gwas = read_gwas_df(filedir=filedir, col_class=col_class, config_class=config_class)
    rsID = list(gwas[:, col_class.rsIDcol_idx])
    rsIDs = set(rsID)
    checked = [] #list contains the rsIDs from the user's input GWAS summary statistics that are found in the database
    skip = []
    
    rmSNPs = open(os.path.join(filedir, "input.snps.rm"), "a")

    # process per chromosome 
    for chrom in range(1, 24):
        for chunk in pd.read_csv(os.path.join(config_class.cfg.get('data', 'dbSNP'), "dbSNP146.chr" + str(chrom) + ".vcf.gz"), header=None, sep="\t", dtype=str, chunksize=10000):
            chunk = np.array(chunk)
            for l in chunk:
                if not l[2] in rsIDs: 
                    continue
                alt = l[4].split(",") #for non-monoallelic
                checked.append(l[2])
                j = bisect_left(rsID, l[2]) #bisect_left functionally is basically gettting the index but is much faster than index
                if col_class.eacol_idx and col_class.neacol_idx:
                    if (gwas[j, col_class.eacol_idx].upper() == l[3] and gwas[j, col_class.neacol_idx].upper() in alt) or gwas[j, col_class.eacol_idx].upper() in alt and gwas[j, col_class.neacol_idx].upper() == l[3]:
                        new_line = [str(chrom), str(l[1]), gwas[j, col_class.neacol_idx].upper(), gwas[j, col_class.eacol_idx].upper(), l[2], str(gwas[j, col_class.pcol_idx])] + write_col(gwas[j,], col_class.orcol_idx, col_class.becol_idx, col_class.secol_idx, col_class.Ncol_idx)
                        out.write('\t'.join(new_line) + "\n")
                    else:
                        logger.info(f"rsID {l[2]} is found in dbSNP but is dropped because (effect allele, non effect allele) from the user's input GWAS summary statistics does not match with the database.")
                        skip.append(l[2])
                        print("\t".join([l[2], " DROPPED BY FUMA: because rsID is found but (effect allele, non effect allele) from the user's input GWAS summary statistics does not match with the database."]), file=rmSNPs)                      
                else: 
                    if len(alt) > 1:
                        logger.info(f"rsID {l[2]} is found in dbSNP but is dropped because this variant has many options for alternative allele and only effect allele is provided.")
                        skip.append(l[2])
                        print("\t".join([l[2], " DROPPED BY FUMA: because rsID is found but is dropped because this variant has many options for alternative allele and only effect allele is provided."]), file=rmSNPs)
                        continue
                    else: 
                        if col_class.eacol_idx:
                            database_alleles = set([l[3], l[4]])
                            if gwas[j, col_class.eacol_idx].upper() in database_alleles:
                                nea = list(database_alleles.difference(gwas[j, col_class.eacol_idx].upper()))[0]
                                new_line = [str(chrom), str(l[1]), nea, gwas[j, col_class.eacol_idx].upper(), l[2], str(gwas[j, col_class.pcol_idx])] + write_col(gwas[j,], col_class.orcol_idx, col_class.becol_idx, col_class.secol_idx, col_class.Ncol_idx)
                                out.write('\t'.join(new_line) + "\n")
                            else:
                                skip.append(l[2])
                                print("\t".join([l[2], " DROPPED BY FUMA: because rsID is found but effect allele from the user's input GWAS summary statistics does not match with the database."]), file=rmSNPs)
                            
                        else:
                            new_line = [str(chrom), str(l[1]), l[3], l[4], l[2], str(gwas[j, col_class.pcol_idx])] + write_col(gwas[j,], col_class.orcol_idx, col_class.becol_idx, col_class.secol_idx, col_class.Ncol_idx)
                            out.write('\t'.join(new_line) + "\n")
                
        if len(gwas) == len(checked):
            break
    
    unresolved_rsIDs = list(set(rsID).difference(checked))
    logger.info(f"These variants were not resolved: {unresolved_rsIDs}")
    logger.info(f"Number of variants that were not resolved: {len(unresolved_rsIDs)}")
    
    for i in unresolved_rsIDs:
        print("\t".join([i, " DROPPED BY FUMA: because rsID is not found in the database."]), file=rmSNPs)

    
    out.close()
    
def check_number_of_lines(filedir):
    """ This function is to check the number of lines in the input.snps file and input.snps.rm file adding up together is equal to the number of lines in the file input.gwas """
    gwas_length = int(subprocess.check_output("cat " + os.path.join(filedir, "input.gwas") + " | wc -l", shell=True).split()[0])
    snps_length = int(subprocess.check_output("cat " + os.path.join(filedir, "input.snps") + " | wc -l", shell=True).split()[0])
    snps_rm_length = int(subprocess.check_output("cat " + os.path.join(filedir, "input.snps.rm") + " | wc -l", shell=True).split()[0])
    
    if snps_length + snps_rm_length != gwas_length:
        logger.warning("The number of variants in the processed gwas file is not equal to the input file.")