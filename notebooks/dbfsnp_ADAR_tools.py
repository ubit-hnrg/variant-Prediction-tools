#import dbfsnp_ADAR_tools.py
import re
import vcf
import gzip
import subprocess
import pandas as pd
import numpy as np

def condition_set(rareness,clin,pato,rev,columns,cured_ADAR_file="../src/df_cln_info_pat2",cured_dbfsnp="~/rare_control_variants_pops_dbNSFP4.0b2a.txt.tsv",homo=None):
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
            
            df=df[['clinvar_id','auto_dominant','auto_recessive', 'x_recessive', 'x_dominant',"Last reviewed"]]
        
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



def forge_graph_features(df_rare_pat,target="target",biogrid_dir="~/BIOGRID/BIOGRID-ORGANISM-3.5.168.tab2/BIOGRID-ORGANISM-Homo_sapiens-3.5.168.tab2.txt"):
    import igraph



    biogrid=pd.read_csv(biogrid_dir,low_memory=False,sep="\t")



    df_biogrid=biogrid[['Entrez Gene Interactor A','Entrez Gene Interactor B']]

    tuples = [tuple(x) for x in df_biogrid.values]
    
    Gm = igraph.Graph.TupleList(tuples, directed = True, edge_attrs = ['weight'])


    dbfsnp_gen= pd.read_csv('~/dbNSFP4.0b2_gene.gz', compression='gzip',sep="\t",low_memory=False)


    dbfsnp_gen.rename(columns={"Ensembl_gene":"Ensembl_geneid"},inplace=True)



    df_ensem_to_entres=dbfsnp_gen[["Ensembl_geneid","Entrez_gene_id"]]
################################################################################
    df_rare_pat=df_rare_pat.merge(df_ensem_to_entres,how="left",on=["Ensembl_geneid"])




    Gm.vs["pathogenicity"]=[id in df_rare_pat.loc[ df_rare_pat[target]>0]["Entrez_gene_id"] for id in Gm.vs["name"]   ]


    #v_patho=[vertex["name"] for vertex in Gm.vs if vertex["pathogenicity"]]
    v_patho=[vertex for vertex in Gm.vs if vertex["pathogenicity"]]


    def get_min_dist_to_patho(Vtx):
        path_to_patho=[len(path)for path in Vtx.get_shortest_paths(v_patho)]
        len_to_patho=list(set(path_to_patho))

        len_to_patho.sort()
        if(len(len_to_patho)<=1):
            return(-1000)

        return(len_to_patho[1])



    Gm.vs["dist_to_patho"]=[ get_min_dist_to_patho(Vtx)  for Vtx in Gm.vs]

   
    def get_number_patho_first_neigh(Vtx):
        if(len(Gm.neighbors(Vtx))):
           return(0)
        sum(Gm.vs[list(set(Gm.neighbors(Vtx)))]["pathogenicity"])


    Gm.vs["patho_first_neigh"]=[ get_number_patho_first_neigh(Vtx)  for Vtx in Gm.vs]


    Gm.save("bio_grid_graph",format="picklez")
    
    
    
    susu=pd.DataFrame( {"Entrez_gene_id":Gm.vs["name"] , "net_dis":Gm.vs["dist_to_patho"],"net_nn":Gm.vs["patho_first_neigh"]})
    
   # df_rare_pat.Entrez_gene_id.loc[df_rare_pat.Entrez_gene_id.isna()]=-2
    susu["Entrez_gene_id"]=susu.Entrez_gene_id.astype(str)
    
    susu=df_rare_pat.merge(susu,how="left",on=["Entrez_gene_id"])
    
    susu.net_dis.fillna(1000,inplace=True)
    
    return(susu)

    
    

    

