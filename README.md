# standalone_fuma
This repo contains a standalone version of FUMA that is intended to run on the command line. 
Usage: for organizations/groups where GWAS summary statistics can't be uploaded to FUMA for security version
**Important notes:**
- Currently support SNP2GENE module (as of September 2026)
- Mirror version 2.2.0
    - Newer updates of SNP2GENE might not be incorporated. 

## Instructions: 
1. Email Tanya Phung (email can be found at https://cncr.nl/ctg/people/) to request sharing of database
- Please include affiliation and reasons why you need to run SNP2GENE on your server and not on FUMA

2. Clone this repo

3. Generate the params.config file
- Create a virtual environment
```
python3 -m venv .
```
- Activate the virtual environment
```
source bin/activate
```
- Install dependency
```
python -m pip install -r venv_requirements.txt
```

- Launch the flask app
```
flask --app flask_app.py --debug run
```

- Go to: http://127.0.0.1:5000/
- Follow the steps to generate the parameter files. Save as `params.config` in a job folder directory such as `snp2gene_jobs/1/params.config`.

4. Install docker if you do not have it installed already

5. Build docker images
```
cd standalone_fuma/docker_images
docker compose -f compose.yaml build gwas_file_cli allsnps_cli magma_cli manhattan_filt_cli qqsnps_filt_cli getld_cli snpannot_cli getgwascatalog_cli geteqtl_cli qtls_map_cli getci_cli genemap_cli create_circos_plot_cli
```

6. Place the input GWAS summary statistics in the folder `snp2gene_jobs/1/` and rename it to `input.gwas`

7. Run: 
```
python fuma_cli.py --filedir `snp2gene_jobs/1/`
```