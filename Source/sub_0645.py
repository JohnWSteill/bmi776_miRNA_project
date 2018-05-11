import unittest
import os
import sys
import tempfile
import filecmp
import glob
import subprocess
import time
import mock
import pdb
import shutil
import time
from gup import *
import Analysis.usageScripts.utilities as ut
import pandas as pd
from itertools import groupby

pd.DataFrame.tabOut = lambda df, f: pd.DataFrame.to_csv(
    df, f, sep='\t', quoting=3, index=True, float_format='%.2f')

def do_miRNA(db, fcid, subNum):
    sampList = g.db.getAllSamples(subNum=subNum)
    ut.doBclThruFfqFromSamps(db=g.db, fcid=fcid, sampList=sampList)
    miClc = mi_rna_ss.MiRnaSs(
        db=g.db,
        fcid=fcid,
        sampList=sampList)
    ct = mi_rna_ss.CollateMiRNACalcTuple(
        db=g.db,
        fcid=fcid,
        subNum=subNum)
    ColMiClc = mi_rna_ss.CollateMiRna(db=g.db, fcid=fcid, subNum=subNum)
    return (ct, ct.getMetadata())

def tamFilter(df):
    ''' use only miRNAs which are above 5 in at least 25% of samples '''
    return (df[df.apply(lambda x: len([el for el in x if el >5])
                >len(df.columns)/4, axis=1)])
    
def get_tam_UQ_norm(df):
    return df/tamFilter(df).quantile(q=.75)

def get_cpm_norm(df):
    return 1e6*df/df.sum()

if __name__ == "__main__":
    '''
    645:           miRNA TS
    '''
    g = run_gup.Gup()
    run = config.FC190
    g.db.importFlowcell(run)
    fcid = g.db.getFcidFromIlluminaDir(run)
    ct, metaData = do_miRNA(db=g.db, fcid=fcid, subNum='645')
    odir = metaData['outdir']
    raw_counts = pd.read_csv(
        glob.glob(os.path.join(odir, '*.tab'))[0], sep='\t', index_col=0) 
    tam_UQ = get_tam_UQ_norm(raw_counts)
    tam_UQ.tabOut(os.path.join(odir,'tam_UQ.tsv'))
    cpm_norm = get_cpm_norm(raw_counts)
    cpm_norm.tabOut(os.path.join(odir,'cpm_norm.tsv'))
    print odir
    print "leaving out copy step:")
    print "cp /mnt/isi4gup/.GUP_HOME/RUNS/HNJ7TBBXX/Collate/Sub_0645_microRNA_Tests_hg19__7b52440fbe646f6a/* Analysis/Sub_0645_miRNA_Norm/data/"






