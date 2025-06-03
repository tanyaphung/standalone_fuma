#!/usr/bin/python
import time
import os
import subprocess
import sys
import re
import pandas as pd
import numpy as np
import tabix
import glob
import configparser
from bisect import bisect_left
from threading import Thread
from collections import defaultdict

class LdParams:
    def __init__(self, filedir, cfg, param):
        self.filedir = filedir
        self.cfg = cfg
        self.param = param
        self.lead_snps = self.check_lead_snps()
        self.add_lead_snps = int(self.param.get('inputfiles', 'addleadSNPs'))
        self.regions = self.get_regions()
        self.ref_panel = self.param.get('params', 'refpanel')
        self.pop = self.param.get('params', 'pop')
        self.lead_p = float(self.param.get('params', 'leadP'))
        self.ref_snps = int(self.param.get('params', 'refSNPs')) #1 to add, 0 to not add
        self.gwas_p = float(self.param.get('params', 'gwasP'))
        self.maf = float(self.param.get('params', 'MAF'))
        self.r2 = float(self.param.get('params', 'r2'))
        self.r2_2 = float(self.param.get('params', 'r2_2'))
        self.merge_dist = int(self.param.get('params', 'mergeDist'))*1000
        self.mhcopt = self.param.get('params', 'MHCopt')
        self.mhc = self.get_mhc()
        self.mhc_region = self.param.get('params', 'extMHC')
        self.mhc_start = self.get_mhc_region()[0]
        self.mhc_end = self.get_mhc_region()[1]
        self.refgenome_dir = self.cfg.get('data', 'refgenome')
        self.annot_dir = os.path.join(self.refgenome_dir, self.ref_panel, "annot")
        
    def check_lead_snps(self):
        lead_snps = self.param.get('inputfiles', 'leadSNPsfile')
        if lead_snps == "NA":
            print("predefined lead SNPs are not provided")
            return None
        else:
            print("predefined lead SNPs are provided")
            return os.path.join(self.filedir, "input.lead")
        
    def get_regions(self):   
        regions = self.param.get('inputfiles', 'regionsfile')
        if regions == "NA":
            print("predefined genomic regions are not provided")
            return None
        else:
            print("predefined genomic regions are provided")
            return os.path.join(self.filedir, "input.regions")
    
    def get_mhc(self):
        mhc = int(self.param.get('params', 'exMHC'))
        if mhc == 1 and self.mhcopt == "magma":
            mhc = 0
        return mhc

    def get_mhc_region(self):
        mhc_start = 29614758
        mhc_end = 33170276
        if self.mhc_region != "NA":
            mhc_start, mhc_end = self.mhc_region.split("-")
        return int(mhc_start), int(mhc_end)

class ColIndexSanitized:
    """ return index of the input.snps file """
    def __init__(self, filedir):
        self.filedir = filedir
        self.chr_idx = 0
        self.pos_idx = 1
        self.nea_idx = 2
        self.ea_idx = 3
        self.rsid_idx = 4
        self.p_idx = 5
        self.or_idx = None
        self.be_idx = None
        self.se_idx = None
        self.get_idx()
        
    def get_idx(self):
        with open(os.path.join(self.filedir, "input.snps"), "r") as f:
            header = f.readline().strip().split("\t")
            if "or" in header:
                self.or_idx = header.index("or")
            if "beta" in header:
                self.be_idx = header.index("beta")
            if "se" in header: 
                self.se_idx = header.index("se")
            

##### Return index of a1 which exists in a2 #####
def array_in(a1, a2):
    return np.flatnonzero(np.isin(a1, a2))

##### return unique element in list #####
def unique(a):
    unique = []
    [unique.append(s) for s in a if s not in unique]
    return unique


##### get annotations for candidate SNPs #####
def get_annot(snps, annot_dir):
    chroms = unique(snps[:,2].astype(int))
    out = []
    ### process per chromosome
    for chrom in chroms:
        annotfile = annot_dir+"/chr"+str(chrom)+".annot.gz"

        tmp = snps[snps[:,2].astype(int)==chrom]
        if len(tmp)==0:
            continue

        ## split snps into chunks
        ranges = []
        start = min(tmp[:,3].astype(int))
        end = min(tmp[:,3].astype(int))
        cur_start = start
        cur_end = end
        for l in tmp:
            if int(l[3])-cur_start < 1000000:
                cur_end = int(l[3])
            else:
                ranges.append([cur_start, cur_end])
                cur_start = int(l[3])
                cur_end = int(l[3])
        ranges.append([cur_start, cur_end])

        tmp = tmp[tmp[:,0].argsort()]
        suid = set(tmp[:,0])

        ## get annotations
        tmp_out = []
        for i in range(0, len(ranges)):
            tb = tabix.open(annotfile)
            annot_tb = tb.querys(str(chrom)+":"+str(ranges[i][0])+"-"+str(ranges[i][1]))
            for l in annot_tb:
                uid = ":".join([l[0], l[1]]+sorted([l[2], l[3]]))
                if uid in suid:
                    j = bisect_left(tmp[:,0], uid)
                    tmp_out.append([tmp[j,2], tmp[j,3], uid]+l[4:])
        tmp_out = np.array(tmp_out)
        tmp_out = tmp_out[np.lexsort((tmp_out[:,0], tmp_out[:,1])), 2:]

        if len(out)==0:
            out = tmp_out
        else:
            out = np.r_[out, tmp_out]
    return out


class LeadSnp:
    """ a class to store information about a lead snp"""
    def __init__(self, lead_snps, iterator, snps, ld):
        self.lead_snps = np.array(lead_snps)
        self.iterator = iterator
        self.snps = np.array(snps)
        self.ld = np.array(ld)
        self.uid = lead_snps[iterator, 0]
        self.rsid = lead_snps[iterator, 1]
        self.chrom = lead_snps[iterator, 2]
        self.position = lead_snps[iterator, 3]
        self.pval = lead_snps[iterator, 4]
        self.rsids, self.uids, self.snps_tmp, self.start, self.end = self.parse_lead_snp()
        
    def parse_lead_snp(self):
        rsids = self.lead_snps[self.iterator,6].split(";") #split rsID

        rsids_loc = array_in(self.snps[:, 1], rsids)
        uids = self.snps[rsids_loc, 0]
        uids_loc = array_in(self.ld[:, 0], uids)
        linked_uids = self.ld[uids_loc, 1]
        
        linked_uids_loc = array_in(self.snps[:, 0], linked_uids)
        snps_tmp = self.snps[linked_uids_loc]

        # Convert position column to int and get range
        positions = snps_tmp[:, 3].astype(int)
        start = positions.min()
        end = positions.max()
        
        return rsids, uids, snps_tmp, start, end
        

class Locus:
    def __init__(self, merge_dist, gidx, lead_snp):
        """ gidx is the index for all chromosomes
        while locus_idx is the index of the locus for each chromosome so it restarts at 0 for the new chromosome"""
        self.merge_dist = merge_dist
        self.gidx = gidx + 1 # index of genomic locus (1-based)
        self.uid = lead_snp.uid # uniqueID
        self.rsid = lead_snp.rsid # rsID
        self.chrom = lead_snp.chrom
        self.position = lead_snp.position
        self.pval = lead_snp.pval
        self._start = lead_snp.start
        self._end = lead_snp.end
        self.in_ind = set(list(lead_snp.rsids)) 
        self.in_lead = set([lead_snp.rsid])
        self.non_gwas_snps = set(list(lead_snp.snps_tmp[lead_snp.snps_tmp[:,7]=="NA", 0]))
        self.gwas_snps = set(list(lead_snp.snps_tmp[lead_snp.snps_tmp[:,7]!="NA", 0]))
    
    def overlaps(self, other) -> bool:
        return self._start - other._end < self.merge_dist
     
    def set_start(self, new_start):
        self._start = min(self._start, new_start)
    
    def start(self):
        return str(self._start)
    
    def set_end(self, new_end):
        self._end = max(self._end, new_end)
        
    def update_ind_sig_snps(self, new_ind_sig_snps):
        self.in_ind.update(list(new_ind_sig_snps))
        
    def update_lead_list(self, new_lead):
        self.in_lead.update([new_lead])
        
    def update_non_gwas_snps(self, new):
        self.non_gwas_snps.update(list(new))
    
    def update_gwas_snps(self, new):
        self.gwas_snps.update(list(new))
        
    def update_lead_snp(self, new_lead_snp):
        self.uid = new_lead_snp.uid
        self.rsid = new_lead_snp.rsid
        self.position = new_lead_snp.position
        self.pval = new_lead_snp.pval
    
    def update_locus(self, new_locus):
        """ update the locus with the next locus if overlapping"""
        self._start = min(self._start, new_locus._start)
        self._end = max(self._end, new_locus._end)
        self.in_ind.update(new_locus.inInd)
        self.in_lead.update([new_locus.inLead])
        self.non_gwas_snps.update(new_locus.non_gwas_snps)
        self.gwas_snps.update(new_locus.gwas_snps)
        
        if float(new_locus.pval) < float(self.pval):
            self.update_lead_snp(new_locus)
    
    # def write_to_file(out_file):
    #     np.savetext(out_file, )
        
class VariantInfo:
    def __init__(self, col_sanitized, v_row):
        self.col_sanitized = col_sanitized
        self.v_row = v_row
        self.v_p = v_row[col_sanitized.p_idx]
        self.v_chr = v_row[col_sanitized.chr_idx]
        self.v_pos = v_row[col_sanitized.pos_idx]
        self.v_rsID = v_row[col_sanitized.rsid_idx]
        self.v_nea = v_row[col_sanitized.nea_idx]
        self.v_ea = v_row[col_sanitized.ea_idx]
        self.v_uid = self.v_uid()

    def v_uid(self):
        chr = str(self.v_chr)
        pos = str(self.v_pos)
        alleles = sorted([self.v_nea, self.v_ea])
        return ":".join([chr, pos, alleles[0], alleles[1]])

class ChrProcessor(Thread):
    def __init__(self, filedir, cfg, param, ld_params, col_sanitized):
        self.filedir = filedir
        self.cfg = cfg
        self.param = param
        self.ld_params = ld_params
        self.col_sanitized = col_sanitized
        self.ld_out_file = os.path.join(self.filedir, 'ld.txt')
        self.snps_out_file = os.path.join(self.filedir, 'snps.txt')
        self.gwas_df = self.get_gwas_df()
        self.gidx = 0
        self.ind_sig_idx = 0
        self.lead_idx = 0

    def get_gwas_df(self):
        return pd.read_csv(os.path.join(self.filedir, "input.snps"), sep='\t')

    def get_ld_tb(self, chr_curr):
        return tabix.open(os.path.join(self.ld_params.refgenome_dir, self.ld_params.ref_panel, self.ld_params.pop, self.ld_params.pop + ".chr" + str(chr_curr) + ".ld.gz"))

    def get_maf_tb(self, chr_curr):
        return tabix.open(os.path.join(self.ld_params.refgenome_dir, self.ld_params.ref_panel, self.ld_params.pop, self.ld_params.pop + ".chr" + str(chr_curr) + ".frq.gz"))
    
    def write_ld(self, chr_curr, ld):
        if chr_curr == 1:
            outfile = open(os.path.join(self.filedir, "ld.txt"), "w")
            print("\t".join(["SNP1", "SNP2", "r2"]), file=outfile)
        else:
            outfile = open(os.path.join(self.filedir, "ld.txt"), "a")
        for i in ld:
            print("\t".join(i), file=outfile)
            
    def write_ind_sig_snps(self, chr_curr, ind_sig_idx, ind_sig_snps, uid2gl):
        addcol = []
        for i in range(0,len(ind_sig_snps)):
            addcol.append([str(ind_sig_idx+i+1), str(uid2gl[ind_sig_snps[i,0]])])
        ind_sig_snps = np.c_[addcol, ind_sig_snps]
        
        ind_sig_idx += len(ind_sig_snps)
        
        if chr_curr == 1:
            outfile = open(os.path.join(self.filedir, "IndSigSNPs.txt"), "w")
            print("\t".join(["No", "GenomicLocus", "uniqID", "rsID", "chr", "pos", "p", "nSNPs", "nGWASSNPs"]), file=outfile)
        else:
            outfile = open(os.path.join(self.filedir, "IndSigSNPs.txt"), "a")
        for i in ind_sig_snps:
            print("\t".join(i), file=outfile)
        
        return ind_sig_idx, ind_sig_snps
    
    def write_lead_snps(self, chr_curr, lead_idx, lead_snps, uid2gl):
        addcol = []
        for i in range(0,len(lead_snps)):
            addcol.append([str(lead_idx+i+1), str(uid2gl[lead_snps[i,0]])])
        lead_snps = np.c_[addcol, lead_snps]
        
        lead_idx += len(lead_snps)
        
        if chr_curr == 1:
            outfile = open(os.path.join(self.filedir, "leadSNPs.txt"), "w")
            print("\t".join(["No", "GenomicLocus", "uniqID", "rsID", "chr", "pos", "p", "nIndSigSNPs", "IndSigSNPs"]), file=outfile)
        else:
            outfile = open(os.path.join(self.filedir, "leadSNPs.txt"), "a")
        for i in lead_snps:
            print("\t".join(i), file=outfile)
        
        return lead_idx
    
    def write_can_snps(self, chr_curr, ld, can_snps, ind_sig_snps):
        pd_ld = pd.DataFrame(ld)
        pd_ld[[2]] = pd_ld[[2]].astype(float)
        idx = pd_ld.groupby(1)[2].transform(max) == pd_ld[2]
        uid1 = np.array(pd_ld[0][idx].tolist())
        uid2 = np.array(pd_ld[1][idx].tolist())
        r2 = np.array(pd_ld[2][idx].tolist())
        tmp = list(can_snps[:,0])
        uid2 = list(uid2)
        idx = [uid2.index(x) for x in tmp]
        uid1 = uid1[idx]
        r2 = r2[idx]
        rsids = can_snps[[tmp.index(x) for x in uid1],1]
        tmp = list(ind_sig_snps[:,2])
        gl = ind_sig_snps[[tmp.index(x) for x in uid1],1]
        can_snps = np.c_[can_snps, r2, rsids, gl]
        
        if chr_curr == 1:
            outfile = open(os.path.join(self.filedir, "snps.txt"), "w")
            print("\t".join(["uniqID", "rsID", "chr", "pos", "non_effect_allele", "effect_allele", "MAF", "gwasP", "r2", "IndSigSNP", "GenomicLocus"]), file=outfile)
        else:
            outfile = open(os.path.join(self.filedir, "snps.txt"), "a")
        
        for i in can_snps:
            print("\t".join(i), file=outfile)
            
            
    def write_loci(self, loci, chr_curr):
        if chr_curr == 1:
            outfile = open(os.path.join(self.filedir, "GenomicRiskLoci.txt"), "w")
            print("\t".join(["GenomicLocus", "uniqID", "rsID", "chr", "pos", "p", "start", "end", "nSNPs", "nGWASSNPs", "nIndSigSNPs", "IndSigSNPs", "nLeadSNPs", "LeadSNPs"]), file=outfile)
        else:
            outfile = open(os.path.join(self.filedir, "GenomicRiskLoci.txt"), "a")
        
        for locus in loci:
            out = [
                str(locus.gidx),
                str(locus.uid),
                locus.rsid,
                str(locus.chrom),
                str(locus.position),
                str(locus.pval),
                str(locus._start),
                str(locus._end),
                str(len(locus.non_gwas_snps) + len(locus.gwas_snps)), #nSNPs
                str(len(locus.gwas_snps)), #nGWASSNPs
                str(len(locus.in_ind)), #nIndSigSNPs
                str(", ".join(map(str, locus.in_ind))), 
                str(len(locus.in_lead)), #nLeadSNPs
                str(", ".join(map(str, locus.in_lead)))
            ]
            
            print("\t".join(out), file=outfile)
            
    def process_chromosomes(self):
        """ calls process_chromosome function to process per chromosome"""
        for chr_curr in self.gwas_df['chr'].unique():
            self.process_chromosome(chr_curr)
        # with open("whateeer", "w") as ouf_file:
        #     [l.write_to_file(ouf_file) for l in loci]

    def process_chromosome(self, chr_curr):
        """ process for each chromosome"""
        self.chr_process(chr_curr)

    def extract_gwas_by_chr(self, chr_curr):
        return self.gwas_df[self.gwas_df['chr'] == chr_curr]

    def filter_by_p(self, gwasfile_chr, p_threshold):
        """
        filter the gwas file by p value threshold
        """
        return gwasfile_chr[
            gwasfile_chr["p"] < p_threshold]

    def check_uid(self, variant_uid, maf_threshold):
        """ check if variant uid exists in the file {pop}.chr*.frq.gz and maf passes"""
        lead_id = False
        lead_maf = False
        chr_curr = str(variant_uid.split(":")[0])
        pos_curr = str(variant_uid.split(":")[1])
        check_id = self.get_maf_tb(chr_curr).querys(str(chr_curr + ":" + pos_curr + "-" + pos_curr))
        for m in check_id:
            tmp_uid = ":".join([m[0], m[1]] + sorted([m[3], m[4]]))
            if tmp_uid == variant_uid:
                lead_id = True
                if float(m[5]) >= maf_threshold:
                    lead_maf = True
                break
        if not lead_id or not lead_maf:
            return "Variant not passed"

    def ld_query(self, v_item: VariantInfo):
        snps_inLD = self.get_ld_tb(v_item.v_chr).querys(str(v_item.v_chr)
                                                              + ":" + str(v_item.v_pos) + "-" + str(v_item.v_pos))
        return [snp for snp in snps_inLD if int(snp[1]) == v_item.v_pos and float(snp[6]) >= self.ld_params.r2]  # get only the variants in LD defined by
        # threshold >= r2

    def maf_query(self, chr, min_pos, max_pos):
        maf_tb_snps_inLD = self.get_maf_tb(chr).querys(
            str(chr) + ":" + str(min_pos) + "-" + str(max_pos))
        return [snp for snp in maf_tb_snps_inLD]
    
    def make_maf_frq_dict(self, maf_tb_ls):
        maf_tb_snps_alleles_dict = defaultdict(list)  # key is the position, value is (allele1, allele2)
        maf_tb_snps_frq_dict = {}
        snpid_rsid_dict = {}
        for snp in list(maf_tb_ls):
            if chr == 6 and self.ld_params.mhc == 1:
                if int(snp[1]) >= self.ld_params.mhc_start and int(snp[1]) <= self.ld_params.mhc_end:
                    continue
            if sorted([snp[3], snp[4]]) not in maf_tb_snps_alleles_dict[int(snp[1])]:
                maf_tb_snps_alleles_dict[int(snp[1])].append(sorted([snp[3], snp[4]]))
            maf_tb_snps_frq_dict[int(snp[1])] = float(snp[5])
            snpid_rsid_dict[snp[0]+":"+snp[1] +":"+ sorted([snp[3], snp[4]])[0] +":"+ sorted([snp[3], snp[4]])[1]] = snp[2]
            
        return maf_tb_snps_alleles_dict, maf_tb_snps_frq_dict, snpid_rsid_dict
    
    def get_lead_snps(self, ind_sig_snps):
        """ get lead snps from independent significant snps"""
        ind_sig_snps = np.array(ind_sig_snps)
        lead_snps = []
        checked = []
        for indsigsnp in ind_sig_snps:
            if indsigsnp[1] in checked:
                continue
            indsigsnp_all = [indsigsnp[1]] #list of snps that are in ld with the snp being examined, including itself
            indsigsnp_inLD = self.get_ld_tb(indsigsnp[2]).querys(str(indsigsnp[2])
                                                              + ":" + str(indsigsnp[3]) + "-" + str(indsigsnp[3]))
            indsigsnp_ld = [snp for snp in indsigsnp_inLD if str(snp[1]) == str(indsigsnp[3]) and float(snp[6]) >= self.ld_params.r2_2] #outputs from the tabix command 
            for line in indsigsnp_ld:
                if int(line[4]) in ind_sig_snps[:,3].astype(int):
                    rsID = ind_sig_snps[ind_sig_snps[:,3].astype(int)==int(line[4]),1][0]
                    checked.append(rsID)
                    indsigsnp_all.append(rsID)
            lead_snps.append([indsigsnp[0], indsigsnp[1], indsigsnp[2], indsigsnp[3], indsigsnp[4], str(len(indsigsnp_all)), ";".join(indsigsnp_all)])
            
        lead_snps = np.array(lead_snps)
        lead_snps = lead_snps[lead_snps[:,3].astype(int).argsort()]
        return lead_snps
    
    def get_genomic_risk_loci(self, gidx, snps, ld, leadSNPs, params):
        loci = []
        iloci = 0

        uid2gl = {}
        for i in range(0, len(leadSNPs)):
            
            lead_snp = LeadSnp(leadSNPs, i, snps, ld)

            new_locus = (i > 0) and not (lead_snp.start <= int(loci[iloci]._end) or lead_snp.start-int(loci[iloci]._end)<params.merge_dist)

            if new_locus:
                gidx += 1
                iloci += 1
            for s in lead_snp.uids:
                uid2gl[s] = gidx+1
            if i==0 or new_locus: #initialize a locus class
                locus = Locus(params.merge_dist, gidx, lead_snp)
                loci.append(locus)
            else:
                locus = loci[iloci]
                locus.set_start(lead_snp.start)
                locus.set_end(lead_snp.end)
                locus.update_ind_sig_snps(lead_snp.rsids)
                locus.update_lead_list(lead_snp.rsid)
                locus.update_non_gwas_snps(lead_snp.snps_tmp[lead_snp.snps_tmp[:,7]=="NA", 0])
                locus.update_gwas_snps(lead_snp.snps_tmp[lead_snp.snps_tmp[:,7]!="NA", 0])
                if float(lead_snp.pval) < float(loci[iloci].pval):
                    locus.update_lead_snp(lead_snp)
                if iloci > 0 and locus.overlaps(loci[iloci-1]):
                    loci[iloci-1].update_locus(locus)
                    iloci -= 1
                    
        loci = np.array(loci)
        gidx += 1
        
        return gidx, loci, uid2gl
    
    def filter_gwas_mhc(self, gwas_currchr, chr_curr):
        """ remove mhc region """
        if chr_curr == 6 and self.ld_params.mhc == 1: 
            pos_col = gwas_currchr.iloc[:, self.col_sanitized.pos_idx].astype(int)
            gwas_currchr = gwas_currchr[(pos_col < self.ld_params.mhc_start) | (pos_col > self.ld_params.mhc_end)]
            return gwas_currchr
        else:
            return gwas_currchr
        
    def filter_gwas(self, chr_curr):
        """ extract by chromosome """
        """ filter by genome-wide p, default is 5e-8 """
        # Extract GWAS data for the current chromosome
        gwas_currchr_tmp = self.extract_gwas_by_chr(chr_curr)
        
        # Filter MHC region if applicable
        gwas_currchr = self.filter_gwas_mhc(gwas_currchr_tmp, chr_curr)
          
        # Remove variants where p values < genome-wide threshold (default of leadP is 5e-8) & then sorted by p values and then by position if p values are the same
        gwas_pfiltered = self.filter_by_p(gwas_currchr, self.ld_params.lead_p)
        gwas_pfiltered.sort_values(by=['p','bp'],ascending=[True,True],inplace=True)
        
        gwas_currchr_allpos = set(gwas_currchr['bp'])
        
        return gwas_currchr, gwas_pfiltered, gwas_currchr_allpos
    
    def get_snps_in_ld(self, v_info): 
        """ get snps in ld """
        # get SNPs in LD, divide into SNPs that are present in the input gwas file and SNPs that are NOT present in the
        ld_records = self.ld_query(v_info) ###chr1  bp1     SNP1    chr2    bp2     SNP2    r2
        
        # get [position, rsID, r2]
        snps_in_ld = [[v_info.v_pos, v_info.v_rsID, 1.0]]
        
        for snp in ld_records:
            if int(snp[1]) != v_info.v_pos:
                continue
            if v_info.v_chr==6 and self.ld_params.mhc==1:
                if int(snp[1]) >= self.ld_params.mhc_start and int(snp[1]) <= self.ld_params.mhc_end:
                    continue
            snps_in_ld.append([int(snp[4]), snp[5], snp[6]])
        
        snps_in_ld = np.array(snps_in_ld, dtype=object)
        snps_in_ld_pos = snps_in_ld[:, 0].astype(int)  # this is a list of all positions
            
        return snps_in_ld, snps_in_ld_pos


    def chr_process(self, chr_curr): 
        
        gwas_currchr, gwas_pfiltered, gwas_currchr_allpos = self.filter_gwas(chr_curr)

        # initialize variables
        ld, can_snps, ind_sig_snps = [], [], []
        checkeduid = set()  # contains a set of checked uid defined as chr:pos:a1:a2

        # process each genome-wide significant snp, starting with the most significant
        for _, row in gwas_pfiltered.iterrows():  
            v_info = VariantInfo(v_row=row, 
                                 col_sanitized=self.col_sanitized)

            if v_info.v_uid in checkeduid:  # if this variant has already been examined, skip it
                continue
            
            if self.check_uid(variant_uid=v_info.v_uid, 
                              maf_threshold=self.ld_params.maf) == "Variant not passed":
                continue

            n_gwas_snps = 0
            n_nongwas_snps = 0

            snps_in_ld, snps_in_ld_pos = self.get_snps_in_ld(v_info)

            # get MAF
            maf_records = self.maf_query(v_info.v_chr, min(snps_in_ld_pos), max(snps_in_ld_pos))
            alleles_dict, frq_dict, snpid_rsid_dict = self.make_maf_frq_dict(maf_records)

            for idx in range(len(snps_in_ld_pos)):
                snp_in_ld_pos = snps_in_ld_pos[idx]
                # snp_in_ld_idx = np.where(snps_in_ld_pos == snp_in_ld_pos)[0][0]
                if not snp_in_ld_pos in alleles_dict:
                    continue
                if snp_in_ld_pos in gwas_currchr_allpos:
                    # get the alleles from input gwas file
                    gwas_alleles = sorted(
                        [gwas_currchr.loc[gwas_currchr['bp'] == snp_in_ld_pos, 'non_effect_allele'].iloc[0],
                         gwas_currchr.loc[gwas_currchr['bp'] == snp_in_ld_pos, 'effect_allele'].iloc[0]])
                    snp_in_ld_uid = ":".join(
                        [str(v_info.v_chr), str(snp_in_ld_pos), gwas_alleles[0], gwas_alleles[1]])

                    if not gwas_alleles in alleles_dict[snp_in_ld_pos]:
                        continue
                    if not frq_dict[snp_in_ld_pos] > self.ld_params.maf:
                        continue
                    if gwas_currchr.loc[gwas_currchr['bp'] == snp_in_ld_pos, 'p'].iloc[0] >= self.ld_params.gwas_p:
                        continue
                    
                    ld.append([v_info.v_uid, snp_in_ld_uid, str(snps_in_ld[idx][2])])
                    
                    if snp_in_ld_uid in checkeduid:
                        continue
                    checkeduid.add(snp_in_ld_uid)
                    n_gwas_snps += 1
                    can_snps.append([snp_in_ld_uid,
                                    gwas_currchr.loc[gwas_currchr['bp'] == snp_in_ld_pos, 'rsID'].iloc[0],
                                    str(chr_curr), 
                                    str(snp_in_ld_pos), 
                                    gwas_currchr.loc[gwas_currchr['bp'] == snp_in_ld_pos, 'non_effect_allele'].iloc[0], 
                                    gwas_currchr.loc[gwas_currchr['bp'] == snp_in_ld_pos, 'effect_allele'].iloc[0],
                                    str(frq_dict[snp_in_ld_pos]),
                                    str(gwas_currchr.loc[gwas_currchr['bp'] == snp_in_ld_pos, 'p'].iloc[0])])
                    
                elif self.ld_params.ref_snps == 1:  # Include variants in reference panel (non-GWAS tagged SNPs in LD)
                    for i in alleles_dict[snp_in_ld_pos]:
                        snp_in_ld_nongwas_uid = ':'.join([str(v_info.v_chr), str(snp_in_ld_pos), i[0], i[1]])
                        
                        ld_toadd = [v_info.v_uid, snp_in_ld_nongwas_uid, str(snps_in_ld[idx][2])]
                        if ld_toadd not in ld: #to avoid duplicates
                            ld.append(ld_toadd)
                            
                        if snp_in_ld_nongwas_uid in checkeduid:
                            continue
                        n_nongwas_snps += 1
                        checkeduid.add(snp_in_ld_nongwas_uid)
                        
                        can_snps.append([snp_in_ld_nongwas_uid, 
                                        snpid_rsid_dict.get(snp_in_ld_nongwas_uid, "NA"), 
                                        str(chr_curr), 
                                        str(snp_in_ld_pos), 
                                        i[0], 
                                        i[1],
                                        str(frq_dict[snp_in_ld_pos]), 
                                        "NA"])
            ind_sig_snps.append(
                [v_info.v_uid, 
                 v_info.v_rsID, 
                 v_info.v_chr, 
                 int(v_info.v_pos), 
                 v_info.v_p, 
                 n_gwas_snps + n_nongwas_snps,
                 n_gwas_snps
                 ])
            
        if len(can_snps) > 0:
            can_snps = np.array(can_snps)
            can_snps = can_snps[can_snps[:, 3].astype(int).argsort()]
            
            
            
            # get leadSNPs from IndSigSNPs
            lead_snps = self.get_lead_snps(ind_sig_snps=ind_sig_snps)
            lead_snps = lead_snps[lead_snps[:, 3].astype(int).argsort()]
            
            self.gidx, loci, uid2gl = self.get_genomic_risk_loci(self.gidx, can_snps, ld, lead_snps, self.ld_params)
            
            # write ld.txt
            self.write_ld(chr_curr, ld)
            
            # write IndSigSNPs.txt
            ind_sig_snps = np.array(ind_sig_snps)
            ind_sig_snps = ind_sig_snps[ind_sig_snps[:, 3].astype(int).argsort()]
            self.ind_sig_idx, ind_sig_snps = self.write_ind_sig_snps(chr_curr, self.ind_sig_idx, ind_sig_snps, uid2gl)
            
            # write leadSNPs.txt
            self.lead_idx = self.write_lead_snps(chr_curr, self.lead_idx, lead_snps, uid2gl)
            
            # write snps.txt
            self.write_can_snps(chr_curr, ld, can_snps, ind_sig_snps)
            
            # write GenomicRiskLoci.txt
            self.write_loci(loci, chr_curr)
            
            print("tmp")
        
        
        
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
            sys.exit("app.config file is not found in the script folder")
    
    def get_param(self):
        params_config_fp = os.path.join(self.filedir, 'params.config')
        if os.path.exists(params_config_fp):
            self.param = configparser.ConfigParser()
            self.param.read(os.path.join(self.filedir, 'params.config'))
        else: 
            sys.exit("params.config file is not found in the job folder")


def main():
    
    # temp for testing
    filedir = "input_dir"
    config_class = Configuration(filedir=filedir)
    
    ld_params = LdParams(filedir, config_class.cfg, config_class.param)
    col_class_sanitized = ColIndexSanitized(filedir=filedir)
    processor = ChrProcessor(filedir=filedir, 
                             cfg=config_class.cfg, 
                             param=config_class.param,
                             ld_params=ld_params,
                             col_sanitized=col_class_sanitized)
    processor.process_chromosomes()
    
    #TODO: add logics for creating the output files
if __name__ == "__main__":
    main()

    
    