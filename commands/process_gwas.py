#!/usr/bin/python
import sys
import os
import re
import pandas as pd
import numpy as np
import configparser
import time
from bisect import bisect_left
import tabix
import csv
import subprocess
import logging
from commands.process_input_helper import validate, sanitize, check_sanitized_gwas_size, write_col_header, process_all_columns_present, process_missing_rsID, process_missing_allele, process_missing_chr_pos, log, check_inputs, initialize_outfile, sort, check_number_of_lines


def process_gwas(filedir, config_class, header_class, col_class, logger):
    start = time.time()
    
    # logger = log(filedir=filedir) #define logger
    
    # config_class, header_class, col_class = validate(filedir=filedir)
    
    sanitize(filedir=filedir, header_class=header_class, col_class=col_class)
    check_sanitized_gwas_size(filedir=filedir)
    sort(filedir=filedir, col_class=col_class)
    
    out = initialize_outfile(filedir=filedir, col_class=col_class) #initialize output file input.snps

	# process input file
    if str(col_class.chrcol_idx) and str(col_class.poscol_idx) and str(col_class.rsIDcol_idx) and str(col_class.eacol_idx) and str(col_class.neacol_idx):
        process_all_columns_present(filedir=filedir, col_class=col_class, out=out)
        
    elif str(col_class.chrcol_idx) and str(col_class.poscol_idx):
        logger.info("Chromosome column and position columns are detected.")

        refpanel = config_class.cfg.get('data', 'refgenome') + "/" + config_class.param.get('params', 'refpanel')
        pop = config_class.param.get('params', 'pop')
        if str(col_class.neacol_idx) and str(col_class.eacol_idx):
            process_missing_rsID(filedir=filedir, col_class=col_class, refpanel=refpanel, pop=pop, out=out)
        else:
            process_missing_allele(filedir=filedir, col_class=col_class, refpanel=refpanel, pop=pop, out=out)
    elif not col_class.chrcol_idx or not col_class.poscol_idx:
        logger.info("Either chr or pos is not provided. Extract chr and pos from rsID using the database")
        process_missing_chr_pos(filedir=filedir, col_class=col_class, config_class=config_class, out=out)
        
    check_number_of_lines(filedir=filedir)
        
    print(f"Processing time: {time.time()-start}")