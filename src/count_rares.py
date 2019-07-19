# este script lo use para contar variantes raras y sin frequencia en exomas clinicos del hospital TSO anotados con gnomad exome and genome via annovar.  
import pandas as pd
import sys
import numpy as np

thr = 0.01
ffile = sys.argv[1]

colnames = ['gnomAD_exome_ALL','gnomAD_genome_ALL']

freq = pd.read_csv(ffile,usecols=colnames,sep = '\t',na_values='.')
rares = (freq[colnames].max(axis = 1) < thr).sum()
unknown = (freq.max(axis = 1).isnull().values).sum()

res = '%s,%s,%s'%(freq.shape[0],rares,unknown)
sys.stdout.write(str(res) + '\n')
