import pipeline_calc as pClc
import calc_tuple
import excel_report
from pool_before_align import getSubNumFromPoolDict
import time
import copy
import subprocess
import os
p_join = os.path.join
import config
import bash_mgr
import errno
import pickle
import pdb
import multiprocessing
import math
import pandas as pd

pd.DataFrame.pos_msk = lambda df, g: df[df.gene_id.isin(g)]
pd.DataFrame.neg_msk = lambda df, g: df[~df.gene_id.isin(g)]
# 3 = magic number for csv.QUOTE_NONE
pd.DataFrame.tabOut = lambda df, f: pd.DataFrame.to_csv(
    df, f, sep='\t', quoting=3, index=False, float_format='%.2f')

def sampSum(df): return df.groupby(by=['Sample'])['TPM'].sum()


class CollateExpressionMeasures(pClc.PipelineCalc):
    """ Manages the collation of Expression Measures."""

    def getCalcMetrics(self, ct):
        """ Best indications of a successful collation:
        1) Sufficient Size of directory
        2) Existence of isoforms.no_mt.ec.tab"""

        od = ct.getMetadata('outdir')
        metrics = {}
        metrics['size'] = sum(os.path.getsize(os.path.join(od, f))
                              for f in os.listdir(od))
        metrics['isoOutExists'] = 'isoforms.no_mt.ec.tab' in os.listdir(od)
        ct.putMetadata(calcMetrics=metrics)

    def makeDirPref(self, calcName=None):
        return super(CollateExpressionMeasures, self).makeDirPref(
            calcName='Collate')

    def buildCalcTupleForArbitraryPooling(self, **kwargs):
        pool = kwargs.pop('pool')
        subNum = getSubNumFromPoolDict(pool)
        if 'subNumPrefix' in kwargs and kwargs['subNumPrefix']:
            subNum = kwargs['subNumPrefix'] + '_' + subNum
        ct = calc_tuple.CalcTuple(
            db=self.db,
            node='Collate',
            fcid=self.fcid,
            subNum=subNum,
            poolId=sorted(pool.keys()),
            pool=pool,
            **kwargs)
        repSamp = ct.dct['pool'].values()[0]
        # TODO a) make this hack a method of collate ct, or b) singleton pools
        # get encased in list
        if isinstance(repSamp, list):
            repSamp = repSamp[0]
        self.calcTuples.append(ct)
        _outdir = p_join(
            self.dirPref,
            "Sub_{0}_{1}_{2}__{3}".format(
                subNum,
                ct.getRefGenome(),
                self.db.getAttrFromSamp(
                    'project_name',
                    repSamp),
                ct.hsh[
                    :config.ODIR_HSH_LEN]))
        if self.tryToLoadFinishedCalc(ct):
            return
        self.buildCalcInfoWithErrAndWrn(ct, _outdir)
        ct.putMetadata(Sample=None)
        subprocess.call(['mkdir', '-p', '-m', '777',
                         ct.getMetadata('outdir')])

    def buildCalcTupleForSampList(self, **kwargs):
        assert kwargs['subNumPrefix'], "need it for sampList collate."
        ct = calc_tuple.CalcTuple(
            db=self.db,
            node='Collate',
            fcid=self.fcid,
            **kwargs)
        self.calcTuples.append(ct)
        if self.tryToLoadFinishedCalc(ct):
            return
        _outdir = p_join(
            self.dirPref,
            "Sub_{0}_{1}_{2}".format(
                ct.dct['subNum'],
                ct.getRefGenome(),
                ct.hsh[:config.ODIR_HSH_LEN])
        )
        self.buildCalcInfoWithErrAndWrn(ct, _outdir)
        ct.putMetadata(Sample=ct.dct['Sample'])
        subprocess.call(['mkdir', '-p', '-m', '777',
                         ct.getMetadata('outdir')])

    def buildCalcTuples(self, **kwargs):
        """ We get the submission numbers, and build a dictionary mapping
        each subNum's sample list. Finally we check if there are multiple
        species, making a calcTuple for each unique genome.

        Planning on three use cases at the moment:
        1) FCID - group all samples in this flowcell by their submission number
        2) By SubNum - Ideally this would span multiple fcids, now just a
                filter. TODO JWS
        3) Custom Pool - This generates its own specific hexidecimal subNum

        There's also three cases for algorithms, commented as Case_I, ...
        1) Sample - Each alignment entity is single sample.
        2) DefaultPooling - Each Alignment entity is a pool based on sample Ids
                I believe all samples in pool have same sample_name in db
        3) ArbPooling - Each Alignment entity is a arbitrary pooling. For
                sample_name we will again use pooId. This is confusing,
                because throughout this code Sample often means 'alignment
                entity'

        """

        # TODO Refactor this, by getting ct's first.
        # A little tricky to test throughly:
        # is always used, even for pools of 1 with the same name as the sample.
        # refGenome (matching all, none, or some of default vaules)
        # Sample Sets: SubNum 1 Fcid, SubNum multFC, custom sample set,
        #    default pool, custom pool
        # Homogenous/Heterogeneous submissions
        # 30 cases
        self.calcTuples = []
        doPooling = False
        if 'pool' in kwargs:
            if isinstance(kwargs['pool'], dict):
                self.buildCalcTupleForArbitraryPooling(**kwargs)  # Case III
                return
            elif kwargs['pool'] is False:
                kwargs.pop('pool')
            else:
                doPooling = True
                assert kwargs['pool'] is True, \
                    "allowed vals for pool: True (False)  or dct"
        elif 'sampList' in kwargs:
            self.buildCalcTupleForSampList(**kwargs)
            return
        if 'refGenome' in kwargs:
            assert 'subNum' in kwargs or 'sampList' in kwargs, \
                "must specify single subNum with refG"
            #refGenome = kwargs.pop('refGenome')
        # else:
            #refGenome = None
        # TODO fcidOnly or MultFC? two descriptions of same thing
        ''' if not ('multFC' in kwargs and kwargs['multFC']):
            kwargs['fcid'] = self.fcid
            kwargs['multFC'] = False '''

        subNumToSamps = self.getSubNumToSampsDict(**kwargs)
        # TODO this shouldn't be here. find When samples are first imported
        # from stemcell and lower them there.
        for samp in self.db.tables['Samp']:
            self.db.tables['Samp'][samp]['genome'] = self.db.tables[
                'Samp'][samp]['genome'].lower()

        for subNum in subNumToSamps:
            if 'subNum' in kwargs:
                kwargs.pop('subNum')
            for genome in set([self.db.getAttrFromSamp('genome', el)
                               for el in subNumToSamps[subNum]]):
                subSampList = sorted([el for el in subNumToSamps[subNum] if
                                      self.db.getAttrFromSamp('genome', el)
                                      == genome])
                prj = self.db.getAttrFromSamp('project_name', subSampList[0])
                pool = {}
                if doPooling:
                    pool = self.getPoolFromSubSamps(subSampList)

                ct = calc_tuple.CalcTuple(
                    db=self.db,
                    node='Collate',
                    subNum=subNum,
                    Sample=subSampList,
                    pool=pool,
                    poolId=sorted(pool.keys()),
                    fcid=self.fcid,
                    **kwargs)

                self.calcTuples.append(ct)
                try:
                    subNumStr = "{0:04d}".format(int(subNum))
                except BaseException:
                    subNumStr = subNum

                _outdir = p_join(self.dirPref,
                                 "Sub_{0}_{1}_{2}__{3}".format(
                                     subNumStr,
                                     prj,
                                     ct.getRefGenome(),
                                     ct.hsh[:config.ODIR_HSH_LEN]))
                if self.tryToLoadFinishedCalc(ct):
                    continue
                self.buildCalcInfoWithErrAndWrn(ct, _outdir)
                ct.putMetadata(Sample=subSampList)
                subprocess.call(['mkdir', '-p', '-m', '777',
                                 ct.getMetadata('outdir')])

    def getPoolFromSubSamps(self, subSampList):
        """ Get default pooling, based on 1st half of sample name.
        this function could be replaced with a call to the imported
        getPoolDictFromSamps, but I'd like to verify I'm not mixing
        and matching at the moment. """
        pool = calc_tuple.getPoolDictFromSamps(subSampList)  # imported
        for poolId, samps in pool.items():
            for samp in samps:
                if samp in subSampList:
                    subSampList.remove(samp)
        assert not subSampList, "haven't got pooling for mix of "\
            "pooled and unpooled samps."
        return pool

    def getSubNumToSampsDict(self, **kwargs):
        """ get samples by subNum, sorted for uniqueness. Default is entire
        flowcell, but will do multiple flowcells for all samples whwn given
        a subNum argument. """
        if 'multFC' in kwargs and kwargs['multFC']:
            if 'fcid' in kwargs:
                kwargs.pop('fcid')
        else:
            kwargs['fcid'] = self.fcid
        samps = self.db.getAllSamples(**kwargs)
        subNumToSamps = {subN: list()
                         for subN in set([self.db.getAttrFromSamp('seq_num', el)
                                          for el in samps])}
        for samp in samps:
            subNum = self.db.getAttrFromSamp('seq_num', samp)
            subNumToSamps[subNum].append(samp)
        for subNum in subNumToSamps:
            subNumToSamps[subNum] = sorted(subNumToSamps[subNum])
        return subNumToSamps

    def configureCalc(self):
        """ We perform this calculation by building dataframes, holding
        various values for all samples. The samples are columns headers. Row
        headers are genes and isoforms with and without mitochondrial labels.
        values include expected counts, tpm, and fpkm, the later rescaled
        after mito removal."""
        def mitoProt(fi, fo, mp):
            pd.read_csv(fi, sep='\t', quoting=3).pos_msk(mp).tabOut(fo)
        for ct in self.filter_calclist_by_finished_and_errors():
            _od = ct.getMetadata('outdir')
            refGenome = ct.getRefGenome()
            refDir = p_join(config.GUP_HOME, 'Ref')
            self.writeHeadersAndSampleNames(ct)
            sampToName = ct.getSampToName()
            try:
                (ensToSym, mitoProtList, genesDscr, isoDscr) =\
                    getRefs(refDir=refDir, genome=refGenome)
                self.processGenes(
                    ct, _od, ensToSym, genesDscr, isoDscr, sampToName)
                self.processIso(ct, _od, ensToSym, isoDscr, sampToName)
                files_in = [
                    'genes.no_mt.tpm.rescale.tab',
                    'genes.no_mt.ec.tab',
                    'isoforms.no_mt.tpm.rescale.tab',
                    'isoforms.no_mt.ec.tab']
                files_out = [el.replace('no_mt', 'mitoprot')
                             for el in files_in]
                try:
                    for f in files_in:
                        assert(os.path.exists(p_join(_od, f)))
                except BaseException:
                    assert ct.getMetadata('errors'), \
                        "no no_mt files, no Err!"
                    continue
                for (f_in, f_out) in zip(files_in, files_out):
                    mitoProt(p_join(_od, f_in),
                             p_join(_od, f_out), mitoProtList)
            except IOError as e:
                self.processNMOGenes(ct, _od, sampToName)
                self.processNMOIso(ct, _od, sampToName)
                wrn = ct.getMetadata('warnings')
                wrn.append('No Annotaion: {}'.format(refGenome))
                ct.putMetadata(warnings=wrn)
            self.move_other_expr_files(ct)
            if not ct.getMetadata('errors'):
                excel_report.CollateReport(ct)

    def move_other_expr_files(self, ct):
        odir = ct.getMetadata('outdir')
        other_exp_files_dir = os.path.join(
            ct.getMetadata('outdir'),
            'OtherExpressionFiles')
        subprocess.call(['mkdir', '-p', '-m', '777', other_exp_files_dir])

        expr_files = [el for el in os.listdir(odir) if el.endswith('.tab')]
        isNMO = not any(['no_mt' in el for el in expr_files])
        if isNMO:
            top_level_suffixes = ['ec.tab', 'tpm.tab']
        else:
            top_level_suffixes = ['no_mt.ec.tab', 'no_mt.tpm.rescale.tab']
        for f in expr_files:
            keep = any([el in f for el in top_level_suffixes])
            if not keep:
                os.rename(os.path.join(odir, f),
                          os.path.join(other_exp_files_dir, f))

    def processNMOGenes(self, ct, _od, sampToName):
        genesDf = self.buildDataFrame(ct, resType='genes')
        if ct.getRefGenome().lower() == 'galgal4':
            geneIdToSym = getGeneIdToSymForGalgal()
        else:
            geneIdToSym = None
        for (tabFilename, valOfInterest) in zip(
                ['genes.ec.tab', 'genes.tpm.tab', 'genes.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            makeGeneIdPivotTabs(df_in=genesDf,
                                fn=p_join(_od, tabFilename),
                                valOfInterest=valOfInterest,
                                sampToName=sampToName,
                                geneIdToSym=geneIdToSym)

    def processGenes(self, ct, _od, ensToSym, genesDscr, isoDscr, sampToName):
        genesDf = self.buildDataFrame(ct, resType='genes')
        genes_mito = genesDf.pos_msk(ensToSym.keys())
        for (tabFilename, valOfInterest) in zip(
                ['genes.mito.ec.tab',
                 'genes.mito.tpm.tab',
                 'genes.mito.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            makeMitoPivotTabs(df_in=genes_mito,
                              fn=p_join(_od, tabFilename),
                              valOfInterest=valOfInterest,
                              isoDscr=isoDscr,
                              ensToSym=ensToSym,
                              sampToName=sampToName)
        genes_no_mt = genesDf.neg_msk(ensToSym.keys())
        for (tabFilename, valOfInterest) in zip(
                ['genes.no_mt.ec.tab',
                 'genes.no_mt.tpm.tab',
                 'genes.no_mt.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            makeGeneIdPivotTabs(df_in=genes_no_mt,
                                fn=p_join(_od, tabFilename),
                                valOfInterest=valOfInterest,
                                genesDscr=genesDscr,
                                sampToName=sampToName)
        rescaleFactors = pd.concat([
            sampSum(genesDf), sampSum(genes_mito)], axis=1)
        rescaleFactors.columns = ['genes_tpm', 'mito_tpm']
        rescaleFactors['no_mtF'] = rescaleFactors['genes_tpm'] / (
            rescaleFactors['genes_tpm'] - rescaleFactors['mito_tpm'])
        df = pd.read_csv(p_join(_od, 'genes.no_mt.tpm.tab.tmpWSampId'),
                         sep='\t', quoting=3)
        os.remove(p_join(_od, 'genes.no_mt.tpm.tab.tmpWSampId'))
        for col in df.columns[1:-1]:
            try:
                df[col] = (df[col] * rescaleFactors['no_mtF'][col])
            except Exception as e:
                print ("temp try to find {}".format(e))
                raise Exception("Un-narrowed Excpt. in collate.py ln 127")
        df = df.rename(columns=sampToName)
        df.tabOut(p_join(_od, 'genes.no_mt.tpm.rescale.tab'))

    def processNMOIso(self, ct, _od, sampToName):
        isoDf = self.buildDataFrame(ct, resType='isoforms')
        if ct.getRefGenome().lower() == 'galgal4':
            geneIdToSym = getGeneIdToSymForGalgal()
        else:
            geneIdToSym = None
        for (tabFilename, valOfInterest) in zip(
                ['isoforms.ec.tab', 'isoforms.tpm.tab', 'isoforms.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            makeIsoPivotTabs(df_in=isoDf,
                             fn=p_join(_od, tabFilename),
                             valOfInterest=valOfInterest,
                             sampToName=sampToName,
                             geneIdToSym=geneIdToSym)

    def processIso(self, ct, _od, ensToSym, isoDscr, sampToName):
        isoDf = self.buildDataFrame(ct, resType='isoforms')
        iso_no_mt = isoDf.neg_msk(ensToSym.keys())
        for (tabFilename, valOfInterest) in zip(
                ['isoforms.no_mt.ec.tab',
                 'isoforms.no_mt.tpm.tab',
                 'isoforms.no_mt.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            makeIsoPivotTabs(df_in=iso_no_mt,
                             fn=p_join(_od, tabFilename),
                             valOfInterest=valOfInterest,
                             isoDscr=isoDscr, sampToName=sampToName)
        iso_mt = isoDf.pos_msk(ensToSym.keys())
        for (tabFilename, valOfInterest) in zip(
                ['isoforms.mito.ec.tab',
                 'isoforms.mito.tpm.tab',
                 'isoforms.mito.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            makeIsoPivotTabs(df_in=iso_mt,
                             fn=p_join(_od, tabFilename),
                             valOfInterest=valOfInterest,
                             isoDscr=isoDscr, sampToName=sampToName)
        rescaleFactors = pd.concat([
            sampSum(isoDf), sampSum(iso_mt)], axis=1)
        rescaleFactors.columns = ['iso_tpm', 'mito_tpm']
        rescaleFactors['no_mtF'] = rescaleFactors['iso_tpm'] / (
            rescaleFactors['iso_tpm'] - rescaleFactors['mito_tpm'])
        df = pd.read_csv(p_join(_od, 'isoforms.no_mt.tpm.tab.tmpWSampId'),
                         sep='\t', quoting=3)
        for col in df.columns[2:-1]:
            df[col] = (df[col] * rescaleFactors['no_mtF'][col])
        df = df.rename(columns=sampToName)
        df.tabOut(p_join(_od, 'isoforms.no_mt.tpm.rescale.tab'))
        os.remove(p_join(_od, 'isoforms.no_mt.tpm.tab.tmpWSampId'))

    def buildDataFrame(self, ct, resType):
        assert resType in ['genes', 'isoforms']
        genesDf = pd.DataFrame()
        dfs = []
        for upCst in ct.getUpstream():
            upStrErr = upCst.getErrors()
            if upStrErr:
                ct.appendToErrors(upStrErr)
                continue
            od = upCst.getMetadata('outdir')
            if 'poolId' in upCst.dct:
                samp = upCst.dct['poolId']
            else:
                samp = upCst.dct['Sample']
            fn = p_join(od, '{}.{}.results'.format(samp, resType))
            try:
                df = pd.read_csv(fn, sep='\t')
                df['Sample'] = samp
                #genesDf = genesDf.append(df)
                dfs.append(df)
            except IOError as e:
                ct.appendToErrors(
                    "{} is not in collation. Err: {}".format(samp, e)
                )
        genesDf = pd.concat(dfs)
        if genesDf.empty:
            assert ct.getMetadat('errors'), "empty genesDf, no Err!"
            raise ValueError("Empty GenesDf")
        return genesDf

    def writeHeadersAndSampleNames(self, ct):
        # TODO not sure why Sample needs to be in info, shouldn't it just be
        # ct.dct['Sample']?
        # info[Sample] == None ? Then use poolIds for Samp Names
        # must be sorted, so don't use ct.dct['pool'] keys
        if 'poolId' in ct.dct:
            sampNames = sorted(ct.dct['poolId'])
            sampIds = sampNames
        else:
            sampNamesAndIds = sorted([(
                self.db.getAttrFromSamp('sample_name', el), el)
                for el in ct.dct['Sample']])
            (sampNames, sampIds) = zip(*sampNamesAndIds)
        _od = ct.getMetadata('outdir')
        with config.openGupOutfile(p_join(_od, 'isoforms.header.tab')) as f, \
                config.openGupOutfile(p_join(_od, 'genes.header.tab')) as g, \
                config.openGupOutfile(p_join(_od, 'sample_order.tab')) as h:
            f.write('#transcript_id\tsymbol\t' + '\t'.join(sampNames))
            g.write('#symbol\t' + '\t'.join(sampNames))
            h.write('\t'.join(sampIds) + '\n')


class CollateExpressionMeasuresNMO(CollateExpressionMeasures):
    ''' NMO: Assumes Non-Model Organisms - no mito symbols, possibly gene or
    isoform names. '''

    def configureCalc(self):
        assert len(self.calcTuples) == 1, 'NMO assumes just 1 collation'
        ct = self.calcTuples[0]
        assert ct.getRefGenome() not in ('mm10', 'mm09', 'hg19'), "model org?"
        ref = ct.getRefFromGenome()
        self.writeHeadersAndSampleNames(ct)
        sampToName = ct.getSampToName()
        od = ct.getMetadata('outdir')
        self.processGenes(ct, od, sampToName)
        self.processIso(ct, od, sampToName)

    def processGenes(self, ct, _od, sampToName):
        genesDf = self.buildDataFrame(ct, resType='genes')
        for (tabFilename, valOfInterest) in zip(
                ['genes.ec.tab', 'genes.tpm.tab', 'genes.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            df = pd.pivot_table(genesDf, values=valOfInterest,
                                index=['gene_id'], columns='Sample')
            df.rename(columns=sampToName).tabOut(fn)

    def processIso(self, ct, _od, sampToName):
        isoDf = self.buildDataFrame(ct, resType='iso')
        for (tabFilename, valOfInterest) in zip(
                ['iso.ec.tab', 'iso.tpm.tab', 'iso.fpkm.tab'],
                ['expected_count', 'TPM', 'FPKM']):
            df = pd.pivot_table(isoDf, values=valOfInterest,
                                index=['iso_id'], columns='Sample')
            df.rename(columns=sampToName).tabOut(fn)


def makeMitoPivotTabs(df_in, fn, valOfInterest,
                      isoDscr, ensToSym, sampToName):
    """ Pivot on sample, bring trans Id's 1st as multiIndex, then reset
    as Column for description merge """
    df = pd.pivot_table(df_in, values=valOfInterest,
                        index=['gene_id', 'transcript_id(s)'],
                        columns=['Sample']).reset_index(level=1)
    pd.merge(df, isoDscr, how='left',
             left_on='transcript_id(s)', right_index=True).\
        sort_index().\
        rename(ensToSym, columns=sampToName).\
        drop('transcript_id(s)', axis=1).\
        reset_index().\
        tabOut(fn)


def makeGeneIdPivotTabs(
        df_in,
        fn,
        valOfInterest,
        sampToName,
        geneIdToSym=None,
        genesDscr=None):
    df = pd.pivot_table(df_in, values=valOfInterest,
                        index=['gene_id'], columns=['Sample'])
    if genesDscr is not None:
        df = pd.merge(df, genesDscr,
                      how='left', left_index=True, right_index=True)
    # Add the 0, 1, 2, 3 back as indx, tabOut ignores it.
    if geneIdToSym:
        df = df.rename(geneIdToSym)
    df.reset_index().\
        rename(columns=sampToName).\
        tabOut(fn)
    if fn.endswith("no_mt.tpm.tab") and genesDscr is not None:
        """ TODO Delete me
        pd.merge(
            df,
            genesDscr,
            how='left',
            left_index=True,
            right_index=True). reset_index(). tabOut(
            fn +
            '.tmpWSampId')
        """
        # rename(ensToSym, columns=sampToName).\
        df.reset_index().tabOut(fn + '.tmpWSampId')


def makeIsoPivotTabs(
        df_in,
        fn,
        valOfInterest,
        sampToName,
        geneIdToSym=None,
        isoDscr=None):
    df = pd.pivot_table(df_in, values=valOfInterest,
                        index=['transcript_id', 'gene_id'],
                        columns=['Sample']).reset_index()
    if isoDscr is not None:
        df = pd.merge(df, isoDscr,
                      how='left', left_on='transcript_id', right_index=True)
    if geneIdToSym:
        df = df.replace({"gene_id": geneIdToSym})
    df.sort_values('gene_id').\
        rename(columns=sampToName).\
        tabOut(fn)
    if fn.endswith("no_mt.tpm.tab") and isoDscr is not None:
        """ TODO Delete me
        pd.merge(df, isoDscr, how='left', left_on='transcript_id',
                 right_index=True).\
            sort_values('gene_id').\
            tabOut(fn + '.tmpWSampId')
        """
        df.sort_values('gene_id').tabOut(fn + '.tmpWSampId')


def getRefs(refDir, genome):
    """
    Retrieve metadata:
        ensToSym:  ensamble:name dict
        mitoProtList: which (ensemble id | name ?) are mito
        genesDesc: multi-word description of genes
        isoDesc: " " " " of isoform function
    """
    ensToSym = getEnsToSym(refDir, genome)
    with open(p_join(refDir, genome + '.mitoprot.symbols.txt')) as f:
        mitoProtList = [el.rstrip() for el in f.readlines()]
    genesDscr = pd.read_csv(p_join(
        refDir, genome + '.geneDescriptions.tab'),
        sep='\t', header=None, index_col=0, names=['gene_id', 'description'])
    genesDscr.drop_duplicates(inplace=True)
    isoDscr = pd.read_csv(p_join(
        refDir, genome + '.isoformDescriptions.tab'), sep='\t', header=None,
        index_col=0, names=['transcript_id', 'description'])
    isoDscr.drop_duplicates(inplace=True)
    return ((ensToSym, mitoProtList, genesDscr, isoDscr))


def getEnsToSym(refDir, genome):
    ensToSymFile = p_join(refDir, genome + '.ensToSym.tab')
    try:
        with open(ensToSymFile) as f:
            return dict([el.rstrip().split('\t') for el in f])
    except IOError as e:
        return None
        return getEnsToSymFromIdxFaFile(genome)  # TODO orphaned


def getEnsToSymFromIdxFaFile(genome):
    """ build an idempotent dictionary. not efficient, but avoids rewrite
    and condition testing in the collation algorithm. """
    # TODO orphaned for now
    # protein_coding_Galgal4_79.idx.fa: >ENSGALT00000000003\nATC....
    def getId(x): return x.strip()[1:]
    idxFile = config.REF_GENOMES[genome] + '.idx.fa'
    #idxFile = p_join(refDir, genome + 'idx.fa')
    ensToSym = {}
    with open(idxFile) as f:
        for line in f:
            ensId = getId(line)
            ensToSym[ensId] = ensId
    return ensToSym


def copyRawRiles(ct, dest=None):
    sampList = ct.dct['Sample']
    if not dest:
        proj = ct.db.getAttrFromSamp('project_name', sampList[0])
        dest = '/isiseqruns/GUP_Deliveries/Sub_0{}_RawFQ_{}'.format(
            subNum, proj)
    fcid_ln = list(set([(
        ct.db.getAttrFromSamp('fcid', s),
        ct.db.getAttrFromSamp('flowcell_lane', s))
        for s in sampList]))
    for (fcid, ln) in fcid_ln:
        bcl_ct = calc_tuple.CalcTuple(
            db=ct.db, node='Bcl2fastq', fcid=fcid, laneNum=ln)
        subSampList = [s for s in sampList if
                       s in ct.db.getAllSamples(**bcl_ct.dct)]
        for samp in subSampList:
            odir = bcl_ct.getSampOutdir(samp)
            intermedDir = os.path.relpath(
                odir, os.path.join(config.GUP_HOME, 'RUNS'))
            toDir = os.path.join(dest, intermedDir)
            subprocess.call(
                ['mkdir', '-p', '-m', '777', os.path.split(toDir)[0]])
            try:
                shutil.copytree(odir, toDir)
            except BaseException:
                print "{} already there?".format(sampSrc)


def getGeneIdToSymForGalgal():
    fn = '/isitools/references/Genomes/gg4/refseq/REFSEQtranscript_gene_id_to_name_LUT.tab'
    with open(fn, 'r') as f:
        nestLst = [el.split() for el in f.readlines()]
    return {el[1]: el[2] for el in nestLst}
