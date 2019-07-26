#import dbfsnp_ADAR_tools.py
import re
import vcf
import gzip
import subprocess
import pandas as pd
import numpy as np

def condition_set(rareness,clin,pato,rev,columns,cured_ADAR_file="df_cln_info_pat2",cured_dbfsnp="~/rare_control_variants_pops_dbNSFP4.0b2a.txt.tsv",homo=None):
    '''
    rareness= binario si tomar el subset con AF segun definimos anteriormente
    homo=vector de homozigosidad minima y maxima para filtrar
    pato= vector con tipos de patogeneidad a seleccionar
    rev= status de review permitidos de la clasificacion
    clin= si usar clinvar para definir categoria
                 '''
    if clin:
        #originalmente
        #df_cln_info_pat=pd.read_csv(cured_file,low_memory=False)
        #el subscript2 solo considera terminos de omim y ademas tiene el Last reviewed 
    
        df_cln_info_pat=pd.read_csv(cured_ADAR_file,low_memory=False)
        
        df_cln_info_pat['patho_rev'] = df_cln_info_pat['patho_rev'].apply(eval)
        

        df=df_cln_info_pat.loc[df_cln_info_pat['patho_status'].isin(pato)]
        
        #print(df.head)


        #df=df.loc[df['patho_rev'].isin(rev)]
        if rev is not None:
            df=df.loc[df.patho_rev.apply(lambda x:x in rev).values]

        #print(df.head)
    if rareness:
        reader=pd.read_csv(cured_dbfsnp,low_memory=False, chunksize=100000,sep="\t",usecols=["clinvar_id","Ensembl_geneid"]+columns+["gnomAD_exomes_controls_nhomalt"])
        
        #variants_dbfnsp=variants_dbfnsp[["clinvar_id","Ensembl_geneid"]+columns+["gnomAD_exomes_controls_nhomalt"]]
        lchunk=[] 
        for k,chunk in enumerate(reader):
            
             if homo != None:
        
                chunkf=chunk.loc[(chunk["gnomAD_exomes_controls_nhomalt"]>homo[0])&(chunk["gnomAD_exomes_controls_nhomalt"]<homo[1])]
                lchunk.append(chunkf)
             else:
                 lchunk.append(chunk)

        variants_dbfnsp = pd.concat(lchunk)
       
    else:
       
        reader = pd.read_csv('/home/felipe/dbNSFP4.0b2a.txt.gz',sep="\t", chunksize=100000,low_memory=False,na_values = '.',usecols=["clinvar_id","Ensembl_geneid"]+columns)
    
        lchunk=[] 
        for k,chunk in enumerate(reader):
            
            ii = chunk["clinvar_id"].isin(df["clinvar_id"])
            if sum(ii)==0:
                continue
            
            chunkf = chunk.loc[ii] 
           # print(chunkf.head)

            #chunkf = chunkf[["clinvar_id"]+columns]

            lchunk.append(chunkf)

        variants_dbfnsp = pd.concat(lchunk)

   
    if clin:
            
            df=df[['clinvar_id','auto_dominant','auto_recessive', 'x_recessive', 'x_dominant']]
        
            variants_dbfnsp=variants_dbfnsp.merge(df,how="right",on="clinvar_id",)
    
    return(variants_dbfnsp)
        
def impute_classes(dataset,patogenas,benignas=None):
    '''
    imputar las clases
        '''
    dataset["target"]=None
    if len(patogenas)==1:
         dataset["target"].loc[dataset[patogenas[0]]]=1
        
    else:
        dataset["target"].loc[dataset[patogenas].sum(axis=1)>0]=1
    if benignas==None:  
        dataset=dataset.loc[~dataset["target"].isna()].copy()
        return(dataset)
    if len(benignas)==1:
         dataset["target"].loc[dataset[benignas[0]]]=0

    else:
        dataset["target"].loc[dataset[benignas].sum(axis=1)>0]=0
    dataset=dataset.loc[~dataset["target"].isna()].copy()
    
    return(dataset)
    
    

def process_proteinwise_scores(word,fun=max):
    ''' para variantes que dan lugar a varias proteinas me quedo con 
             el scoremas alto'''

    if pd.isna(word):
        return(None)
    if not (type(word) is str):
        return(word)
        
    averga=re.split(";",word)
    averga=set(averga)
    averga.discard(";")
    averga.discard(".")
    averga.discard(",")
    averga.discard("")
    
    if averga==set():
        return(None)
    
    mara=map(float,averga)

    av=list(mara)

    return(fun(av))

    
def get_ens(obj):
    '''me quedo con el primer ensemble id
               '''
    #aver=re.split(r';',obj)
    #aver=aver[0]
    obj=str(obj)
    
    if not (';' in obj):
        return(obj)
    aver=re.search(r'(.*?);', obj).group(1)
    #k=""
    #for i in aver:
     #   k=k+i
    #return(k)
    return(aver)        



