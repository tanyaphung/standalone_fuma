#!/usr/bin/python
import shutil

import pandas as pd
import numpy as np
import sys
import os
import re
import subprocess

def is_float(s):
	try:
		return float(s)
	except ValueError:
		return False

def main():
	##### check argument #####
	if len(sys.argv)<2:
		sys.exit('ERROR: not enough arguments\nUSAGE ./allSNPs.py <filedir>')

	##### add '/' to the filedir #####
	filedir = sys.argv[1]

	for snps in pd.read_csv(
		os.path.join(filedir, "input.snps"),
		dtype={'chr': 'Int64', 'bp': 'Int64', 'p': 'float64'},
		sep="\t",
		usecols=['chr', 'bp', 'p'],
		chunksize=50000
	):
		snps = snps.dropna()
		snps.to_csv(os.path.join(filedir, "all.txt"), header=False, index=False, sep="\t", mode="a")

	subprocess.run(
    ["bgzip", "-f", os.path.join(filedir, "all.txt")],
    check=True
)

	subprocess.run(
    [
        "tabix",
        "-s", "1",
        "-b", "2",
        "-e", "2",
        os.path.join(filedir, "all.txt.gz")
    ],
    check=True
)
	# path = os.path.abspath(os.path.join(filedir, "all.txt.gz"))

	# print("FILE:", path)
	# print("EXISTS:", os.path.exists(path))
	# print("SIZE:", os.path.getsize(path))
	# print("TABIX:", shutil.which("tabix"))

	# print(
	# 	subprocess.run(
	# 		["tabix", "--version"],
	# 		capture_output=True,
	# 		text=True
	# 	).stdout
	# )

	# result = subprocess.run(
	# 	["tabix", "-p", "vcf", path],
	# 	capture_output=True,
	# 	text=True
	# )

	# print("RETURN CODE:", result.returncode)
	# print("STDERR:", result.stderr)

if __name__ == "__main__": main()
