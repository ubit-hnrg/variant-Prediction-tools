import pandas as pd
import numpy as np
import os
import argparse 
import time
import gc
# utility functions
def parseArgs():
    parser = argparse.ArgumentParser(description='This script extracts rare variants from control populations in ExAC, GenomAD, EPS6500 project, 1000Genomas and UK cohort')
    parser.add_argument('-i','--input_file',required=True,help='input gzip file without index (dbNSFP4.0b1c dic/2018 or greater)')
    parser.add_argument('-o','--output_file',required=False,default = '.')
    parser.add_argument('-c','--chunksize',dest = 'chunksize',default = 10**6,required = False,type=int)
    parser.add_argument('-m','--max_frequency_threshold',default = 0.01,required = False,type=float)
    args = parser.parse_args()
    return args
    

# Criteria used for rare variants:
# those variants 
def extract_rareVariants_from_dbNSFP(in_f, out_f, size, freqcols , max_frequency_threshold):
    reader = pd.read_table(in_f, chunksize=size,na_values = '.')
    
    t0 = time.time()
    for k,chunk in enumerate(reader):
        minutes = (time.time()-t0)/float(60)
        if k&(k-1)==0:
            print 'processed %s*%s rows. Elapsed, %s minutes'%(k+1,size,minutes)
        ii = (chunk[freqcols].isna()).sum(axis=1)>0
        chunkf = chunk[ii] 
        max_freq_observed = chunkf[freqcols].max(axis=1)
        chunkf['max_freq_observed'] =  max_freq_observed
        
        # keep variants with max_freq_observed under threshold  
        # ALSO exclude those with maximum freq = 0.0000 (only for sanity)
        flag_rare_variant =  (max_freq_observed < max_frequency_threshold)&(max_freq_observed >0 )
        rares = chunkf[flag_rare_variant]

        if k==0:
            rares.to_csv(out_f, index=False,sep='\t') 
        else:
            rares.to_csv(out_f, index=False, header=False, mode='a',sep='\t')
        gc.collect()
    return(None)

# 1) get variants with at least 1 frequency field not nul
# 2) filter out those with the maximum allele frequency across all population greather than the threshold, 0.01 by default. 
def define_freq_control_fields():

    KGP = [u'1000Gp3_AF', u'1000Gp3_AFR_AF', u'1000Gp3_EUR_AF', \
    u'1000Gp3_AMR_AF',u'1000Gp3_EAS_AF', u'1000Gp3_SAS_AF']     # 1000 Genomes Project

    UK = [u'UK10K_AF']                                          #Alternative allele frequency in combined genotypes in UK10K cohort (TWINSUK+ALSPAC).
    ESP6500 = ['ESP6500_AA_AF',u'ESP6500_EA_AF']                #Alternative allele frequency in the samples of the NHLBI GO Exome Sequencing Project (ESP6500 data set).

    ExAC=[u'ExAC_AF', u'ExAC_Adj_AF', u'ExAC_AFR_AF',u'ExAC_AMR_AF',\
    u'ExAC_EAS_AF', u'ExAC_FIN_AF', u'ExAC_NFE_AF',u'ExAC_SAS_AF']      # ExAC populations (60,706 samples)
                                                                        # "_Adj_" fields correspond to Adjusted Alt allele frequency (DP >= 10 & GQ >= 20)
    GENOMAD = [u'gnomAD_exomes_controls_POPMAX_AF',u'gnomAD_genomes_controls_POPMAX_AF'] #Maximum allele frequency across GenomAD CONTROL populations 
                                                                                        # (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry)
                                                                                        #  exome (54,704 samples) and 
                                                                                        # genome (~5000 samples)
    fields = KGP + UK + ESP6500 + ExAC + GENOMAD
    return fields


def main():
    #read args
    args = parseArgs()
    input_file, output_file, size, max_frequency_threshold = args.input_file, args.output_file, args.chunksize, args.max_frequency_threshold
    if output_file == '.':
        output_file =  './rare_control_variants_'+os.path.basename(input_file).split('.gz')[0]+'.tsv'

    # frequency fields of considered control populations. (1000GProj, UK cohort, EPS6500, ExAC, GenomAD)
    fields = define_freq_control_fields()


    #process dbNSFP data:
    rares = extract_rareVariants_from_dbNSFP(in_f=input_file,out_f=output_file,size=size,freqcols=fields,max_frequency_threshold = max_frequency_threshold)




if __name__ == "__main__":
    main()
