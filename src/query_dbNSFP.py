#!/usr/bin/python
import pandas as pd
import numpy as np
import os
import argparse 
import time
import gc

in_f = '/home/vpt/data/dbNSFP/dbNSFP4.0a/dbNSFP4.0a.txt.gz'
out_f = '/home/vpt/data/benchmark/Jinchen2018/jinchen_dbNSFP40a.tsv'
query_f = '/home/vpt/data/benchmark/Jinchen2018/supplementary-Table-S2-The-testing-missense-variants-of-benchmark-datasets.tsv'
size = 1000000


def extract_rareVariants_from_dbNSFP(in_f, out_f, size, query):
    reader = pd.read_csv(in_f, chunksize=size,na_values = '.',sep = '\t',dtype={'hg19_pos(1-based)':str,'hg19_chr':str})
    
    t0 = time.time()
    for k,chunk in enumerate(reader):
        minutes = (time.time()-t0)/float(60)
        if k&(k-1)==0:
            print 'processed %s*%s rows. Elapsed, %s minutes'%(k+1,size,minutes)
#        if k==5:
#            break

        chunk.index = chunk.hg19_chr.map(str) + '-' + chunk['hg19_pos(1-based)'].map(str) + '-' + chunk.ref.map(str) + '-' + chunk.alt.map(str)
        query_chunk = pd.merge(chunk,query,left_index = True, right_index= True, how = 'inner') 

        # keep variants with max_freq_observed under threshold  
        # ALSO exclude those with maximum freq = 0.0000 (only for sanity)

        if k==0:
            query_chunk.to_csv(out_f, index=True,sep='\t') 
	    dflist = []
	    dflist.append(query_chunk)
        else:
            query_chunk.to_csv(out_f, index=True, header=False, mode='a',sep='\t')
	    dflist.append(query_chunk)
        gc.collect()
    df = pd.concat(dflist)
    return(df)


query = pd.read_csv(query_f,sep = '\t')
query['Chr'] = query.Chr.str.strip('chr')
query.index = query.Chr.map(str) +'-' + query.Start.map(str) +'-' + query.Ref.map(str) +'-' + query.Alt.map(str)


df = extract_rareVariants_from_dbNSFP(in_f = in_f,out_f = out_f,size = size,query = query)
df.to_csv(out_f+'.backup', index=True,sep='\t')
