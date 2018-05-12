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


def tamFilter(df):
    ''' use only miRNAs which are above 5 in at least 25% of samples '''
    return (df[df.apply(lambda x: len([el for el in x if el >5])
                >len(df.columns)/4, axis=1)])
    
def get_tam_UQ_norm(df):
    return df/tamFilter(df).quantile(q=.75)

def get_cpm_norm(df):
    return 1e6*df/df.sum()

def do_670_with_100gn_goldStandard(db):
    samps = g.db.getAllSamples(subNum='670')
    samps += [el for el in db.getAllSamples(sub='645') 
            if '100ng' in db.getAttrFromSamp('sample_name',el) ]
    ct, _ = ut.do_miRNA(db=db,
            fcid = db.getFcidFromIlluminaDir(config.FC192),
            sampList = samps,
            subNum = 's670_p_645_100ng')
    do_miRNA_postProc(ct)

def do_miRNA_postProc(col_ct):
    odir = col_ct.getMetadata('outdir')
    raw_counts = pd.read_csv(
        glob.glob(os.path.join(odir, '*.tab'))[0], sep='\t', index_col=0) 
    tam_UQ = get_tam_UQ_norm(raw_counts)
    tam_UQ.tabOut(os.path.join(odir,'tam_UQ_sub_{}.tsv'.format(sub)))
    cpm_norm = get_cpm_norm(raw_counts)
    cpm_norm.tabOut(os.path.join(odir,'cpm_norm_sub_{}.tsv'.format(sub)))
    print odir



if __name__ == "__main__":
    '''
    645:           miRNA TS
    '''
    g = run_gup.Gup()
    cts = g.db.getAllCalcTuples(node='CollateMiRna')
    runs = [('645', config.FC190), ('670',config.FC192)]
    for sub,run in runs:
        g.db.importFlowcell(run)
    do_670_with_100gn_goldStandard(g.db)
    for sub,run in runs:
        fcid = g.db.getFcidFromIlluminaDir(run)
        ct, metaData = ut.do_miRNA(db=g.db, fcid=fcid, subNum=sub)
        do_miRNA_postProc(ct)

