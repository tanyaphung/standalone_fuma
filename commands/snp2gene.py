import click
import os
import subprocess
import configparser
from commands.process_gwas import process_gwas
from commands.magma import magma
from commands.process_input_helper import log, validate
# from getLD import get_ld

@click.command()
@click.option('--filedir', '-d', required=True, help='Path to the directory folder.'
              'The file input.gwas is expected in this folder.')

def snp2gene(filedir):
    logger = log(filedir=filedir) #define logger
    config_class, header_class, col_class = validate(filedir=filedir)
    # process_gwas(filedir, config_class, header_class, col_class, logger)
    if str(config_class.param.get('magma', 'magma')) == str(1):
        magma(filedir)
    # subprocess.call("Rscript " + os.path.dirname(os.path.realpath(__file__)) + "/plot_manhattan.R --job_dir " + filedir, shell=True)