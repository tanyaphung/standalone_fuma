import click
import subprocess
import os
from commands import snp2gene
from commands import params
# import dash_bio
# import pandas as pd


@click.group()
def fuma():
    pass

# Register the subcommand
fuma.add_command(snp2gene.snp2gene)
fuma.add_command(params.params)

if __name__ == '__main__':
    fuma()