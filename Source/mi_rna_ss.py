'''
we have a spike-in c.elegans 67
1) let the reference grow
2) create a new reference and add parameter to calctuple
3) Add spike-in argument

for now let reference grow.

plan:
    * delete the 2) choice
    * delete CTs
    * rerun

TODO: remove this before checkin
'''


import config
import bash_mgr
import pipeline_calc as pClc
import calc_tuple
import errno
import os
import glob
import pdb
import subprocess
import time
import multiprocessing
import bt_rsem
import calc_tuple
import fastq_filter
import collate
from collate import pd as pd  # keeps my tabOut


class miRNACalcTuple(calc_tuple.BtRsemCalcTuple):

    def buildCalcTup(self, **kwargs):
        """ The miRNA-specific logic for building CalcTuple. """

        self.node = 'miRNA'
        kwargUniverse = ['Sample', 'refGenome']
        self.dct = {
            k: v for k,
            v in self.dct.items() if k in kwargUniverse and v}
        assert 'Sample' in kwargs, "Need Sample for miRNA ct"
        if 'refGenome' in kwargs:
            self.popRefGenomeIfDefault(**kwargs)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def getUpstream(self, node=None):
        'similar to FilterCalcTuple but no ReadRescue'
        dct = self.dct.copy()
        samp = dct.pop('Sample')
        dct['fcid'] = self.db.getAttrFromSamp('fcid', samp)
        dct['laneNum'] = self.db.getAttrFromSamp('flowcell_lane', samp)
        return calc_tuple.CalcTuple(db=self.db, node='Bcl2fastq', **dct)


class CollateMiRNACalcTuple(calc_tuple.CollateCalcTupleSampleList):

    def buildCalcTup(self, **kwargs):
        self.node = 'CollateMiRna'
        kwargUniverse = ['fcid', 'subNum', 'multFC', 'Sample', 'refGenome',
                         'mismatch']
        assert 'fcid' in kwargs, "Need fcid for collate ct"
        assert 'subNum' in kwargs, "Need fcid for collate ct"
        assert 'sampList' not in kwargs, "no more sampList"
        assert 'subNumPrefix' not in kwargs, "no more subNumPrefix"
        if 'Sample' in kwargs and kwargs['Sample']:
            self.getSubNum(**kwargs)
            self.dct['Sample'] = sorted(kwargs['Sample'])
        else:
            self.putSamplesInDctCheckingForMultFC(**kwargs)
        self.checkForRefGenome(**kwargs)
        self.filterDct(kwargUniverse)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def getUpstream(self, node=None):
        dct = self.dct.copy()
        samps = dct.pop('Sample')
        return [miRNACalcTuple(
            db=self.db, node='miRNACalcTuple', Sample=el, **dct)
            for el in samps]


class CollateMiRna(collate.CollateExpressionMeasures):

    def buildCalcTuples(self, **kwargs):
        self.calcTuples = []
        if 'Sample' in kwargs and kwargs['Sample']:
            self.buildCalcTupleForSampList(**kwargs)
        else:
            if 'Sample' in kwargs:
                kwargs.pop('Sample')
            self.buildCalcTupleForSubNum(**kwargs)

    def buildCalcTupleForSampList(self, **kwargs):
        ct = CollateMiRNACalcTuple(db=self.db, fcid=self.fcid, **kwargs)
        self.calcTuples.append(ct)
        if self.tryToLoadFinishedCalc(ct):
            return
        _outdir = os.path.join(
            self.dirPref,
            "Sub_{0}_{1}_{2}".format(
                ct.dct['subNum'],
                ct.getRefGenome(),
                ct.hsh[:config.ODIR_HSH_LEN])
        )
        self.buildCalcInfoWithErrAndWrn(ct, _outdir)
        subprocess.call(['mkdir', '-p', '-m', '777',
                         ct.getMetadata('outdir')])

    def buildCalcTupleForSubNum(self, **kwargs):
        subNumToSamps = self.getSubNumToSampsDict(**kwargs)
        kwargs.pop('subNum')
        for subNum in subNumToSamps:
            for genome in set([self.db.getAttrFromSamp('genome', el)
                               for el in subNumToSamps[subNum]]):
                subSampList = sorted([el for el in subNumToSamps[subNum] if
                                      self.db.getAttrFromSamp('genome', el)
                                      == genome])
                ct = CollateMiRNACalcTuple(
                    db=self.db,
                    subNum=subNum,
                    Sample=subSampList,
                    fcid=self.fcid,
                    **kwargs)
                self.calcTuples.append(ct)
                if self.tryToLoadFinishedCalc(ct):
                    continue
                try:
                    subNumStr = "{0:04d}".format(int(subNum))
                except BaseException:
                    subNumStr = subNum
                prj = self.db.getAttrFromSamp('project_name', subSampList[0])
                _outdir = os.path.join(self.dirPref,
                                       "Sub_{0}_{1}_{2}__{3}".format(
                                           subNumStr,
                                           prj,
                                           ct.getRefGenome(),
                                           ct.hsh[:config.ODIR_HSH_LEN]))
                self.buildCalcInfoWithErrAndWrn(ct, _outdir)
                subprocess.call(['mkdir', '-p', '-m', '777',
                                 ct.getMetadata('outdir')])

    def configureCalc(self):
        '''
        For s in samples:
            read in fullcounts -> list of dfs
            read in stats -> list of dfs
        Add them together
        print
        '''
        for ct in self.calcTuples:
            _od = ct.getMetadata('outdir')
            sampToName = ct.getSampToName()
            miRNA_df = pd.DataFrame()
            dfs = []
            for up_ct in ct.getUpstream():
                if up_ct.getErrors():
                    ct.appendToErrors(up_ct.getErrors())
                    continue
                samp = sampToName[up_ct.dct['Sample']]
                df = pd.read_csv(
                    os.path.join(
                        up_ct.getMetadata('outdir'),
                        'miRNA.fullcounts.tab'),
                    sep='\t',
                    index_col=3,
                    header=None,
                    names=['Sequence', 'NUniq', samp, '#symbol'])
                df = df.drop(['Sequence', 'NUniq'], axis=1)
                dfs.append(df)
            miRNA_df = pd.concat(dfs, axis=1)
            miRNA_df.to_csv(os.path.join(_od, 'miRNA.ec.tab'), sep='\t',
                            quoting=3, index=True, index_label='#symbol',
                            float_format='%.0f')

    def getCalcMetrics(self):
        self.ct.putMetadata(
            calcMetrics={'nSamps': len(self.ct.getUpstream())})


class MiRnaSs(bt_rsem.BtRsem):

    def makeDirPref(self, calcName=None):
        ''' want grandparent, super will change it to BtRsem'''
        return super(bt_rsem.BtRsem, self).makeDirPref(calcName='MiRnaSs')

    def buildCalcTuples(self, sampList=None, **kwargs):
        self.calcTuples = []
        sampList = self.getSampList(sampList=sampList, **kwargs)
        for samp in sampList:
            ct = miRNACalcTuple(
                db=self.db, Sample=samp, **kwargs)
            self.calcTuples.append(ct)
            _outdir = os.path.join(self.dirPref, str(samp) + "__" +
                                   ct.hsh[:config.ODIR_HSH_LEN])
            if self.tryToLoadFinishedCalc(ct):
                continue
            self.buildCalcInfoWithErrAndWrn(ct=ct, outdir=_outdir)
            subprocess.call(['mkdir', '-p', '-m', '777', _outdir])

    def configureCalc(self):
        """ Build self.jobs bash job array """
        _calcList = self.filter_calclist_by_finished_and_errors()
        for ct in _calcList:
            ct.putMetadata(status='Running')
            _jb = self.buildMiRNAFilterBashJobs(ct, nproc=1)
            if _jb:  # Errors will return None
                self.jobs.append(_jb)
        if self.jobs:
            bash_mgr.BashMgr(self.jobs)
        self.jobs = []
        for ct in _calcList:
            _jb = self.buildMiRNAAlignBashJobs(ct, nproc=1)
            if _jb:  # Errors will return None
                self.jobs.append(_jb)
        if self.jobs:
            bash_mgr.BashMgr(self.jobs)
        self.jobs = []

    def buildMiRNAFilterBashJobs(self, ct, nproc=1):
        '''
        TS.filter.minlen15.sh 055082_0076295 TTAGGC &

        zcat ../../$1*_R1_*.gz \
                | ../../scripts/ffq_TS_mRNA.pl --minlen 15 --muxIdx $2 \
                -l ffq.out
                --ftags ffq.tags.tab
                --ftrimmed_adapters ffq.trimmedAdapters.tab \
                > $1.fq 2> ffq.log

        '''
        samp = ct.dct['Sample']
        odir = ct.getMetadata('outdir')
        bcl_ct = ct
        while bcl_ct.getUpstream():
            bcl_ct = bcl_ct.getUpstream()
        bcl_samp_dir = bcl_ct.getSampOutdir(samp)
        fqz_file = glob.glob(os.path.join(bcl_samp_dir, '*_R1*.gz'))[0]
        sampIndexType = self.db.getIndexTypeFromSamp(ct.dct['Sample'])
        assert sampIndexType in ['TS', 'LM']
        sampIndexLabel = self.db.getAttrFromSamp(
            'index_label', ct.dct['Sample'])
        cmd1 = fastq_filter.CMD_ZCAT.format(fqz_file)
        # add --minlen 15 argument to filter command:
        split_key = ' --mux'
        cmd_pieces = fastq_filter.CMD_LM_NE_TS.split(split_key)
        cmd_pieces[0] += ' --minlen 15 '
        filter_cmd = split_key.join(cmd_pieces)
        cmd2 = filter_cmd.format(
            fastq_filter.SCRIPT[sampIndexType],
            sampIndexLabel,
            odir)
        bashCommand = [cmd1, cmd2]
        return bash_mgr.BashJob(
            bashCommand=bashCommand,
            errFile=os.path.join(odir, 'ffq.err'),
            outFile=os.path.join(odir, ct.dct['Sample'] + '.fq'),
            workingDir=odir
        )

    def buildMiRNAAlignBashJobs(self, ct, nproc=1):
        ''' run_hg19.short.miRNAs.sh
        short.miRNAs.sh hg19 055083_0076296 > hg19.short.miRNAs.sh.out \
                2> hg19.short.miRNAs.sh.log &

        cat ../../$1.miRNAs.tab div.txt  ../ffq.tags.tab | perl -nle \
            'if ($f){
                ($s,$c)=split /\t /;
                $m=substr($s,0,15);
                if ($t{$m}){
                    $k{$m}++;
                    $kk{$m}+=$c
                }
            }elsif(m/^#$/){
                $f=1
            }else{
                ($m, $c,$t,@f)=split /\t/;
                $t{$m}=$t; $kk{$m}=$k{$m}=0
            } END {
                map {print join("\t", $_, $k{$_}, $kk{$_}, $t{$_})} keys %k}'\
                | sort -ke,3nr -k2,2nr -k4,4 > miRNA.fullcounts.tab
        '''

        samp = ct.dct['Sample']
        refGenome = ct.getRefGenome()
        assert ct.getRefGenome() in ['mm10', 'hg19']
        ref_file = os.path.join(config.GUP_SRC_HOME, 'Ref',
                                '{}.miRNAs.tab'.format(refGenome))
        odir = ct.getMetadata('outdir')
        infile = os.path.join(odir, 'ffq.tags.tab')
        script = os.path.join(config.GUP_SRC_HOME, 'FilterScripts',
                              'short.miRNAs.sh')
        outLog = '{}.short.miRNAs.sh.out'.format(refGenome)
        errLog = '{}.short.miRNAs.sh.log'.format(refGenome)
        bashCommand = '{0} {1} {2}'.format(script, ref_file, infile)

        return bash_mgr.BashJob(
            bashCommand=bashCommand,
            errFile=os.path.join(odir, errLog),
            outFile=os.path.join(odir, outLog),
            workingDir=odir
        )

    def getCalcMetrics(self, ct):
        ''' nReads: total aligned.
            percAligned '''

        metrics = ct.getMetadata('calcMetrics')
        samp = ct.dct['Sample']
        odir = ct.getMetadata('outdir')

        miRNAcounts = pd.read_csv(
            os.path.join(odir, 'miRNA.fullcounts.tab'),
            sep='\t',
            index_col=3,
            header=None,
            names=['Sequence', 'NUniq', 'NReads', 'miRNA'],
        )

        metrics['nReads'] = miRNAcounts['NReads'].sum()
        reads_in = ct.getUpstream().getMetadata('calcMetrics')['nReads'][samp]
        metrics['percAligned'] = metrics['nReads'] * 100. / reads_in
        ct.putMetadata(calcMetrics=metrics)
