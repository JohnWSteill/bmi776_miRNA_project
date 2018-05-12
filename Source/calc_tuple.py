import md5
import xmltodict
import glob
import os
import copy
import config
import hashlib
from itertools import groupby


def CalcTuple(db, node=None, ctStr=None, **kwargs):
    """
        Creates a tuple filled with data
        @param node - type of calculation
        @param ctStr - name of the calc tuple
        @param kwargs - calculation specific parameters
    """

    assert node or ctStr, "CalcTuple needs node or ctStr to be specified"
    if ctStr:
        node = ctStr.split('\n')[-1].split(':')[0]
    assert node in ['Bcl2fastq', 'ReadRescue', 'Filter', 'Pool',
                    'BtRsem', 'Collate'], "invalid NodeType: {}".format(node)
    if node == 'Bcl2fastq':
        return BclCalcTuple(db=db, ctStr=ctStr, **kwargs)
    if node == 'ReadRescue':
        return ReadRescCalcTuple(db=db, ctStr=ctStr, **kwargs)
    if node == 'Filter':
        return FilterCalcTuple(db=db, ctStr=ctStr, **kwargs)
    if node == 'Pool':
        return PoolCalcTuple(db=db, ctStr=ctStr,
                             **kwargs)
    if node == 'BtRsem':
        return BtRsemCalcTuple(db=db, ctStr=ctStr, **kwargs)
    if node == 'Collate':
        try:
            if 'sampList' in kwargs:
                return CollateCalcTupleSampleList(db=db, **kwargs)
            return CollateCalcTuple(db=db, ctStr=ctStr, **kwargs)
        except AssertionError as e:
            genomes = db.getGenomeFromSubNum(subNum=kwargs['subNum'])
            assert len(genomes) > 1, \
                "except is to catch mult genomes! {}".format(e)
            return [CollateCalcTuple(
                db=db, genome=el, **kwargs) for el in genomes]


class _CalcTuple(object):
    """ The CalcTuple object is an economical representation of a unique
    atomic calculation.

    It has a tuple object, a unique string representation, and a correspoding
    hash value. (All of which are immutable and potentially useful for keys)
    The tuple values are strings, format: "kw:value(s)". It also has a
    dictionary representation exposing these values.
    """

    def __init__(self, db, ctStr=None, **kwargs):
        """
            Initializes calcTuple object

            @param db - GUP database object
            @param ctStr - name of the calc_tuple
            @param kwargs - calculation specific parameters
        """
        self.db = db
        self.dct = kwargs.copy()
        if ctStr:
            kwargs.update(dict([el.split(':') for el in ctStr.split('\n')]))
        self.tup = 'Gup:{}'.format(config.VERSION),
        self.buildCalcTup(**kwargs)  # can edit self.dct
        # Omit the Gup version in dct
        for k, v in [el.split(':') for el in self.tup[1:]]:
            if k not in self.dct:
                self.dct[k] = v
        self.hsh = self.getHash()

    def buildCalcTupFromDctAndOrderedKwargUniv(self, kwargUniverse,
                                               silentKeys=list()):
        """ 2 requirements here: 1) dct keys should be a subset of
            kwargUniverse, and 2) kwargUniverse specifies CalcTup.tup order,
            which in turn gives uniq str representation, giving uniq hsh.
        """

        dct = copy.deepcopy(self.dct)
        for k in kwargUniverse:
            if k in dct:
                if isinstance(dct[k], list):
                    # in the tups, we want str representation only
                    dct[k] = ','.join(dct[k])
                self.tup = self.tup + ('{}:{}'.format(k, dct.pop(k)),)
        if not silentKeys:
            silentKeys = list()
        for k in dct:
            assert k in silentKeys, \
                "an unblessed kwarg slipped through: {}".format(dct)

    def __str__(self):
        """
            Returns the calc_tuple name
        """
        return '\n'.join(["Gup CalcTuple Object: " + self.node] +
                         ['\t' + el for el in self.tup]) + '\n\n'

    def __repr__(self):
        """
            returns calc_tuple's representation
        """
        return str(self)

    def __eq__(self, other):
        """
            How to know if 2 calc_tuple's are equal
            Needed for hashing
        """
        return str(self) == str(other)

    def __ne__(self, other):
        """
            How to know if 2 calc_tuple's are not equal
        """
        return not(self == other)

    def __hash__(self):
        """
            Hash the calc_tuple
            tuple is the key, data is stored in dict
        """
        return hash(str(self))

    def __iter__(self):
        self.iterDelivered = False
        return self

    def next(self):
        if not self.iterDelivered:
            self.iterDelivered = True
            return self
        raise StopIteration
        return

    def getCalcInfoFilename(self):
        '''
            returns /path/calcInfo_hsh.p
        '''

        return os.path.join(
            config.GUP_CALCINFO,
            'calcInfo_' + self.hsh + '.p')

    def findIndexOf(self, attr):
        """
            Returns index of the data
        """
        return [el.split(':')[0] for el in self.tup].index(attr)

    def getUpstream(self, node=None):
        """
            Returns previous data processing step
            Implemented in children
        """
        raise NotImplementedError

    def getNode(self):
        """
            retrieve current data processing step
        """
        return str(self).split('\n')[-1].split(':')[0]

    def buildCalcTup(self, **kwargs):
        """
            Creates ct object
            Implemented in children
        """
        raise NotImplementedError

    def popMismatchIfDefault(self, **kwargs):
        ''' trace back popRefGenomeIfDefault so I can imitate'''
        default = self._getDefaultMismatch(**kwargs)
        if default == kwargs['mismatch']:
            self.dct.pop('mismatch')

    def getDefaultMismatch(self):
        ''' Wrapper for _getDefaultMismatch

        In buildCalcTup I need this functionality, even though self is not
        yet fully formed, and so there I can call the private function with the
        needed fcid and lane
        '''
        return self._getDefaultMismatch(
            fcid=self.dct['fcid'], laneNum=self.dct['laneNum'])

    def _getDefaultMismatch(self, **kwargs):
        try:
            assert 'fcid' in kwargs and 'laneNum' in kwargs
        except AssertionError:
            if 'Sample' in self.dct:
                if isinstance(self.dct['Sample'], str):
                    samp = self.dct['Sample']
                else:
                    samp = self.dct['Sample'][0]
            elif 'pool' in self.dct:
                samp = self.dct['pool'].values()[0][0]
            kwargs['fcid'] = self.db.getAttrFromSamp('fcid', samp)
            kwargs['laneNum'] = self.db.getAttrFromSamp('laneNum', samp)
        index_type = self.db.getIndexTypeFromFcidLn(
            fcid=kwargs["fcid"],
            laneNum=kwargs["laneNum"])
        return config.INDEXTYPE_TO_MISMATCH[index_type]

    def getMismatch(self):
        if 'mismatch' in self.dct:
            return self.dct['mismatch']
        return self.getDefaultMismatch()

    def popRefGenomeIfDefault(self, **kwargs):
        try:
            default = self._getDefaultGenome(self.dct['Sample'])
        except BaseException:
            default = None
        if default == kwargs['refGenome']:
            self.dct.pop('refGenome')

    def getDefaultGenome(self):
        """ Wrapper for _getDefaultGenome

        in buildCalcTup I need this functionality, even though self is not
        yet fully formed, and so there I can call the private function with the
        a needed sample list.
        """
        assert 'Sample' in self.dct or 'pool' in self.dct, \
            "have to have Sample or pool to call this"
        if 'Sample' in self.dct:
            return self._getDefaultGenome(self.dct['Sample'])
        else:
            sampsTop = self.dct['pool'].values()
            # This is fine if this is a list of strings,
            # but may need to flatten:
            samps = []
            for top in sampsTop:
                # (falseVal, trueVal)[boolean_expresion] trick
                samps += ([top], top)[isinstance(top, list)]
            return self._getDefaultGenome(samps)

    def _getDefaultGenome(self, samps):
        """ samps could be a list (for a pool) or just a str for one sample.
        To avoid writing two similar functions, recast str to one-element
        list.
        """
        if isinstance(samps, str):
            samps = [samps]
        genomes = list(set([self.db.getAttrFromSamp('genome', el)
                            for el in samps]))
        assert len(genomes) == 1, \
            "All samps should be same genome, or specify refGenome"
        return genomes[0].lower()

    def getRefGenome(self):
        """ Mainly to avoid checking in code. Also looking forward to the
        future when I'll need to differentiate between species and references.
        """
        if 'refGenome' in self.dct:
            return self.dct['refGenome']
        else:
            return self.getDefaultGenome()

    def getRefFromGenome(self, genome=None):
        ''' TODO: write unittest in calctuple '''
        if not genome:
            genome = self.getRefGenome()
        try:
            ref = config.REF_GENOMES[genome]
        except KeyError:
            ref = os.path.join(
                config.GUP_HOME,
                'CustomReferences',
                genome,
                genome)
            if not os.path.exists(os.path.dirname(ref)):
                raise KeyError("{} not found ".format(genome))
        return ref

    def getHash(self):
        """
            Return hash value
        """
        return md5.md5(str(self)).hexdigest()

    def addDictToTup(self, ct):
        """
            Use calc_tuple to create dict
            @param ct - calc_tuple object
        """
        d = {}
        for el in ct[1:]:
            (k, v) = el.split(':')
            d[k] = v
        return (repr(d),) + ct[1:]

    def getAttr(self, attr):
        """
            Retrieves the attribute

            If attribute not in current calc_tuple,
            try from previous calc_tuples

            @param attr - attribute related to calculation
        """
        ct = self
        attrVal = None
        while ct and attr not in ct.dct:
            ct = ct.getUpstream()
        if ct:
            return ct.dct[attr]

    def getMetadata(self, fields=None):
        """
        Retrieve Calculation Metadata. If fields is a string
        just a value is returned. If multiple fields are given,
        a dictionary is returned.

        If a field is not defined, None is returned.
        """
        info = self.db.getCalc(self)
        if not info:
            return None
        if isinstance(fields, str):
            if fields not in info:
                return None
            else:
                return info[fields]
        elif fields:
            return {el: info[el] for el in fields if el in info}
        else:
            return info.copy()

    def putMetadata(self, **kwargs):
        """
        Write calculation metadata.
        """
        # TODO Validate field(s)
        calcInfo = self.db.getCalc(self)
        if not calcInfo:
            calcInfo = {}
        for kw in kwargs:
            calcInfo[kw] = kwargs[kw]
        self.db.writeCalc(self, calcInfo)

    def getErrors(self):
        calc = self.getMetadata()
        isInDb = bool(calc)
        if not isInDb:
            return "{} is not in db.".format(self)

        if 'errors' in calc and calc['errors']:
            return calc['errors']
        isOnDisk = os.path.exists(self.getCalcInfoFilename())
        if not isOnDisk:
            return "{} is in db but not filesystem".format(self)
        hasData = ('outdir' in calc and os.path.exists(calc['outdir'])
                   and os.listdir(calc['outdir']))
        if not hasData:
            return "{} has an empty calc['outdir']".format(self)
        if not ('calcMetrics' in calc and calc['calcMetrics']):
            return "{} has no calcMetrics".format(self)

    def appendToErrors(self, newErr):
        errors = self.getMetadata('errors')
        errors.append(newErr)
        self.putMetadata(errors=errors)

    def getAllUpstream(self):
        up_cts = self.getUpstream()
        all_up_cts = list(up_cts)  # copy
        for el in up_cts:
            up_up_cts = el.getAllUpstream()
            if up_up_cts:
                all_up_cts += list(up_up_cts)
        return list(set(all_up_cts))

    def getUpstreamErrors(self):
        """ only gets 1st, likely maybe no info about which one """
        upCts = self.getUpstream()
        if upCts:
            for upCt in upCts:
                err = upCt.getErrors()
                if err:
                    return(err)

    def delete_readOneOnly_if_false(self):
        if 'readOneOnly' not in self.dct:
            # not there to begin with
            return True
        if not self.dct['readOneOnly']:
            # I expect a key to evaluate as True
            del self.dct['readOneOnly']
            return True

    def check_if_readOneOnly_should_be_popped(self):
        if self.delete_readOneOnly_if_false():
            return True
        up_cts = self.getUpstream()
        if all([el.check_if_readOneOnly_should_be_popped() for el in up_cts]):
            del self.dct['readOneOnly']
            return True
        return False

    def getSampToName(self):
        """ Translates sample_id to sample_name
        3 cases:
            Sample
            pool
            PoolId/Sample"""
        if 'poolId' in self.dct:
            fail = False
            sampToName = {}
            pl = self.dct['pool']
            for k in pl:
                names = list(set([self.db.getAttrFromSamp('sample_name', s)
                                  for s in pl[k]]))
                if len(names) == 1:
                    sampToName[k] = names[0]
                else:
                    fail = True
            if fail:
                sampToName = dict([(el, el) for el in self.dct['poolId']])
        else:
            sampToName = dict([(el, self.db.getAttrFromSamp('sample_name', el))
                               for el in self.dct['Sample']])
        vals = sampToName.values()
        dups = [el for el in set(vals) if vals.count(el) > 1]
        if dups:
            for samp in sampToName:
                if samp in dups:
                    sampToName[samp] = sampToName[samp] + "_" + samp
        return sampToName


class BclCalcTuple(_CalcTuple):

    def buildCalcTup(self, **kwargs):
        """
            Creates the calc tuple obj for a bcl2FastQ calculation
            Current functional parameters are in kwargUniverse.

            @param kwargs - calculation specific params

        """
        kwargUniverse = ['fcid', 'laneNum', 'mismatch']
        self.dct = {
            k: v for k,
            v in self.dct.items() if k in kwargUniverse and v}
        self.node = 'Bcl2fastq'
        assert 'laneNum' in kwargs, 'Need laneNum for Bcl2fastq'
        assert 'fcid' in kwargs, 'Need fcid for Bcl2fastq'
        if 'mismatch' in kwargs:
            self.popMismatchIfDefault(**kwargs)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def getUpstream(self, node=None):
        """
            BCL2FastQ is the first calculation - no upstream
        """
        return []

    def getErrors(self):
        superErr = super(BclCalcTuple, self).getErrors()
        if superErr:
            return superErr
        if not self.checkFastqFilesStillThere():
            return "Missing Demux Dirs or Fastq.gz files"
        if self.db.getIs3kFromFcid(self.dct['fcid']):
            samp_order_err = self.check3kSampleOrderInCsvBug()
            if samp_order_err:
                return samp_order_err
            err_file_error = self.check_3k_err_file()
            if err_file_error:
                return err_file_error
        else:
            err_file_error = self.check_2p5_err_file()
            if err_file_error:
                return err_file_error
        return None

    def checkFastqFilesStillThere(self):
        samps = self.db.getAllSamples(**self.dct)
        sampDirs = [self.getSampOutdir(el) for el in samps]
        missingDirs = [el for el in sampDirs if not os.path.exists(el)]
        if missingDirs:
            return False
        for d in sampDirs:
            if not glob.glob(os.path.join(d, '*.fastq.gz')):
                return False
        return True

    def check3kSampleOrderInCsvBug(self):
        for samp in self.db.getAllSamples(**self.dct):
            gzFiles = os.listdir(self.getSampOutdir(samp))
            # '057471_0085050__S11_L006_I1_001.fastq.gz'
            # ['057471', '0085050', '', 'S11', 'L006', 'I1', '001.fastq.gz']
            sNums = set([el.split('_')[3] for el in gzFiles])  # S11
            if len(sNums) > 1:
                return samp + " has too many reads: CSV duplication bug"

    def getSampOutdir(self, samp):
        laneDir = self.getMetadata('outdir')
        return os.path.join(laneDir, self.getSampRelOutdir(samp))

    def getSampReadFiles(self, samp):
        d = self.getSampOutdir(samp)
        files = [el for el in os.listdir(d) if '_R1_' in el or '_R2_' in el]
        return([os.path.join(d, f) for f in files])

    def getSampRelOutdir(self, samp):
        """
            function because 2500 version prepends with 'Project_'
            @param sampInfo - information about Sample
        """
        if self.db.getIs3kFromFcid(self.dct['fcid']):
            projDir = self.db.getAttrFromSamp('project_name', samp)
            sampDir = samp
        else:
            projDir = 'Project_' + \
                self.db.getAttrFromSamp('project_name', samp)
            sampDir = 'Sample_' + samp
        return os.path.join(projDir, sampDir)

    def check_3k_err_file(self):
        """

        """
        _od = self.getMetadata('outdir')
        errF = _od + '_CBTF.err'
        with open(errF) as f:
            for line in f:
                pass
        if 'Processing completed with 0 errors' in line:
            return None
        else:
            return line

    def check_2p5_err_file(self):
        """
            open error file, and as long as there are no error messages,
        """
        _od = self.getMetadata('outdir')
        try:
            with open(os.path.join(_od, 'make.err'), 'r') as f:
                errLines = f.readlines()
        except IOError as e:
            errMsg = "no make.err in {}".format(_od)
        good_end = "INFO: all completed successfully.\n"
        # new possible last line
        # make: warning:  Clock skew detected.  Your build may be incomplete.
        if (good_end in errLines[-1] or good_end in errLines[-2]):
            errMsg = None
        else:
            errMsg = "bcl2fastq Failed: make.err"
        return errMsg

    def buildMaskString(self):
        """
        We here build the argument for the -use-bases-mask flag. See pg 27:
        http://support.illumina.com/content/dam/illumina-support/documents/\
                documentation/software_documentation/bcl2fastq/\
                bcl2fastq2-v2-17-software-guide-15051736-g.pdf

        This algorithm will be incorrect anytime there are single-end TS
        samples on a paired-end flowcell!
        """
        # TODO: improve by getting isPaired from seq_status db
        itype = self.db.getIndexTypeFromFcidLn(
            fcid=self.dct['fcid'], laneNum=self.dct['laneNum'])
        ilen = config.PROTOCOL_TO_INDEX_LENGTH[itype]
        runInfof = os.path.join(
            config.ILLUMINA_RUNS,
            self.db.getIlluminaDirFromFcid(self.dct['fcid']),
            'RunInfo.xml')
        with open(runInfof) as f:
            runInfoDict = xmltodict.parse(f.read())
        reads = runInfoDict['RunInfo']['Run']['Reads']['Read']
        assert len(reads) in [2, 3, 4]
        nCycles = [int(el['@NumCycles']) for el in reads]
        isInd = [el['@IsIndexedRead'] == 'Y' for el in reads]
        assert isInd in [
            [False, True],  # Single-End Read
            [False, True, True],  # Nextera ??
            [False, True, False],  # Paired End TS, AT
            [False, True, True, False]]  # Paired End NE`

        def nonIndRdStr(n, ignore=False):
            if ignore:
                return "n{}".format(n)
            return "Y{}".format(n)

        def indRdStr(ilen, n, ignore=False):
            if ignore:
                istr = "n{}".format(n)
            elif n == ilen:
                istr = "I{}".format(ilen)
            else:
                istr = "I{}n{}".format(ilen, n - ilen)
            return istr

        mstr = [nonIndRdStr(nCycles[0])]
        rd_indx = 1
        assert isInd[rd_indx]
        mstr.append(indRdStr(ilen=ilen, n=nCycles[rd_indx]))
        while len(reads) > rd_indx + 1:
            rd_indx += 1
            if isInd[rd_indx]:
                ignoreInd2 = (itype != 'NE')
                mstr.append(indRdStr(ilen=ilen, n=nCycles[rd_indx],
                                     ignore=ignoreInd2))
            else:
                # The assumption here is that if a read is capable of PE
                # reads, it must be doing paired end. This is not necessarily
                # true: one could have PE Flowcell with both PE and SE
                # TruSeq preps, for example.
                # Thus: TODO get the SE/PE protocol from seq_status.
                # TODO commented-out code:
                # ignoreRead2 = (itype == 'LM')
                ignoreRead2 = False  # Always get R2s
                mstr.append(
                    nonIndRdStr(n=nCycles[rd_indx], ignore=ignoreRead2))
        return ','.join(mstr)


class BtRsemCalcTuple(_CalcTuple):
    """ The BtRsemCalcTuple is upstream of collation, and
    downstream of Pool and Filter. If no pool, fcid isn't
    needed, as it's redundant. (This will need to be reconsidered
    if we decide to cast all cases as a pool of at least 1.)
    """

    def buildCalcTup(self, **kwargs):
        """ The bt_rsem-specific logic for building CalcTuple. """

        self.node = 'BtRsem'
        kwargUniverse = ['fcid', 'poolId', 'Sample', 'refGenome', 'mismatch',
                         'aligner', 'readOneOnly']
        self.dct = {
            k: v for k,
            v in self.dct.items() if k in kwargUniverse and v}
        if 'pool' in kwargs and kwargs['pool']:
            if kwargs['pool'] is True:
                poolmates = self.db.getPoolMates(kwargs['Sample'])
                if poolmates:
                    self.dct['Sample'] = poolmates
                    self.dct['poolId'] = self.db.getPoolId(kwargs['Sample'])
            else:
                assert kwargs['poolId'] and kwargs['Sample'], \
                    "Haven't written this yet, use poolID and "\
                    "Sample keywords to express arbitrary pooling."
        assert 'Sample' in kwargs, "Need Sample for BtRsem ct"
        if 'poolId' in self.dct:
            assert 'fcid' in self.dct, "Because Pools can span flowcells, the "\
                "fcid assignment of a pooled sample is arb, but must be "\
                "specified to tell which directory to write it to."
        else:
            if 'fcid' in self.dct:
                self.dct.pop('fcid')
        if 'refGenome' in kwargs:
            if kwargs['refGenome']:
                allowed_vals = (config.REF_GENOMES.keys() +
                                [os.path.basename(el) for el in os.listdir(
                                    os.path.join(config.GUP_HOME, 'CustomReferences'))])
                assert kwargs['refGenome'] in allowed_vals,\
                    "bad ref: {}".format(kwargs['refGenome'])
            self.popRefGenomeIfDefault(**kwargs)
        if 'mismatch' in kwargs:
            self.popMismatchIfDefault(**kwargs)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def getUpstream(self, node=None, **kwargs):
        """ Get's upstream Node, either Pool or Filter. """
        assert(node is None), "Haven't written logic for arb. upstrm node yet."
        dct = kwargs.copy()
        dct.update(self.dct)
        if 'poolId' in dct:
            return CalcTuple(db=self.db, node='Pool', **dct)
        return CalcTuple(db=self.db, node='Filter', **dct)

    def getAligner(self):
        ''' TODO could have a refGenome not in config.REF_GENOMES
        so anybody calling this should have some excpetion
        handling, or I do it here with None '''
        return _getAlignerFromRef(self.getRefFromGenome())

    def getErrors(self):
        '''changeover from fracAligned to percAligned will go away '''
        superErr = super(BtRsemCalcTuple, self).getErrors()
        if superErr:
            return superErr
        metrics = self.getMetadata('calcMetrics')
        if metrics and 'fracAligned' in metrics:
            metrics['percAligned'] = metrics.pop('fracAligned')
            self.putMetadata(calcMetrics=metrics)


def _getAlignerFromRef(ref):
    """
    The refGenome isnt a file, but a prefix to a family of files.

    For example: rsem.minke_plus_trans_plus_blue_SOX2 refers about a dozen
    files, including rsem.minke_plus_trans_plus_blue_SOX2.transcripts.fa,
    .grp, .bt2, etc.

    We get the aligner by examining these file extensions.
    """
    extensions = [el.split('.')[-1] for el in glob.glob(ref + '*')]
    if 'bt2' in extensions:
        return 'bowtie2'
    elif 'ebwt' in extensions:
        return 'bowtie'
    else:
        raise LookupError("Couldn't find appropriate aligner")


class CollateCalcTuple(_CalcTuple):

    def buildCalcTup(self, **kwargs):
        """
        The wrinkle here is that self.dct has a 'pool' dictionary that
        is not in self.tup. I don't want to go through machinations of creating
        a sorted str representation for uniqueness, and it could be too long
        to belong to classes string reprensentation.

        To maintain hashable uniqueneness, there is a method,
        pool_before_align.getSubNumFromPoolDict, which gives a hex subNum
        from sorted poolDict hash.
        """
        self.node = 'Collate'

        kwargUniverse = ['fcid', 'subNum', 'multFC', 'poolId', 'Sample',
                         'refGenome', 'aligner', 'mismatch', 'readOneOnly']
        silentKeys = ['pool']
        self.filterDct(kwargUniverse, silentKeys)
        assert 'subNum' in kwargs, "Need Submission # (subNum) for collate ct"
        assert 'fcid' in kwargs, "Need fcid for collate ct"
        if 'pool' in kwargs and kwargs['pool'] is True:
            # pool= True is code for default pooling
            self.dct['pool'] = getPoolDictFromSamps(
                self.db.getAllSamples(**kwargs))
            if self.dct['pool']:
                self.dct['poolId'] = sorted(self.dct['pool'].keys())
            else:
                # even though kwarg['pool'] was true, there weren't any
                # samples to pool.
                self.dct.pop('pool')
        assert not (('poolId' in self.dct) and ('Sample' in self.dct)), \
            "Collate ct takes PoolId or Sample, not both"
        if 'poolId' in self.dct:
            assert 'pool' in self.dct and isinstance(self.dct['pool'], dict), \
                'need pool dictionary to guarentee uniq'
            if 'Sample' in self.dct:
                self.dct.pop('Sample')  # Already assured if its empty if here.
            if 'mismatch' in kwargs:
                self.popMismatchIfDefault(**kwargs)
            self.buildCalcTupFromDctAndOrderedKwargUniv(
                kwargUniverse, silentKeys=['pool'])
            return
        elif 'Sample' in self.dct:
            if 'poolId' in self.dct:
                self.dct.pop('poolId')  # Already assured if its empty if here.
        else:
            self.putSamplesInDctCheckingForMultFC(**kwargs)
            if 'genome' in self.dct:
                self.dct.pop('genome')
        self.checkForRefGenome(**kwargs)
        if 'mismatch' in kwargs:
            self.popMismatchIfDefault(**kwargs)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def putSamplesInDctCheckingForMultFC(self, **kwargs):
        # If we want samples from multiple flowcells, we need to remove
        # fcid from filter. Still important, as ct uses it for id, output
        _kw = kwargs.copy()
        if 'multFC' in _kw and _kw['multFC']:
            _kw.pop('fcid')
        self.dct['Sample'] = sorted(self.db.getAllSamples(**_kw))

    def checkForRefGenome(self, **kwargs):
        """
        In the case where no ref given, samples have to all be same.
        If ref is given but it the default reference anyway, it's popped so
            will produce same calcTuple.
        """
        if 'refGenome' in kwargs:
            self.popRefGenomeIfDefault(**kwargs)
        else:
            assert self._getDefaultGenome(self.dct['Sample']), \
                "Multiple genomes in this collation, need to subdivide."

    def filterDct(self, kwargUniverse, silentKeys=[]):
        """ implicitly takes out everything that doesn't belong by rebuilding.
        Using pop() to removed undesired would take more lines but perhaps be
        more clear. """
        self.dct = {k: v for k, v in self.dct.items() if
                    k in kwargUniverse + silentKeys and v}

    def getUpstream(self, node=None, **kwargs):
        assert(node is None), "Haven't written logic for arb. upstrm node yet."
        dct = kwargs.copy()
        dct.update(copy.deepcopy(self.dct))
        if 'poolId' in dct:
            poolId = dct.pop('poolId')
            pool = dct.pop('pool')
            return [
                CalcTuple(
                    db=self.db,
                    node='BtRsem',
                    poolId=el,
                    Sample=pool[el],
                    **dct) for el in poolId]
        samps = dct.pop('Sample')
        dct.pop('subNum')
        if 'genome' in dct:
            dct.pop('genome')
        if 'fcid' in dct:
            dct.pop('fcid')
        return [CalcTuple(
                db=self.db, node='BtRsem', Sample=el, **dct)
                for el in samps]

    def getErrors(self):
        '''changeover from tsv to tab will go away '''
        superErr = super(CollateCalcTuple, self).getErrors()
        if superErr:
            return superErr
        if self.getMetadata('outdir'):
            tab_files = glob.glob(
                os.path.join(self.getMetadata('outdir'), '*.tsv'))
            if tab_files:
                return "Tsv files in outdir"


class CollateCalcTupleSampleList(CollateCalcTuple):
    """
    Non-pooled case when Samplist is specified.
    """

    def buildCalcTup(self, **kwargs):
        assert 'fcid' in kwargs, "Need fcid for collate ct"
        assert 'sampList' in kwargs, "Need sampList for this collate ct"
        assert 'pool' not in kwargs, "No pool"
        assert 'poolId' not in kwargs, "No pool"

        self.node = 'Collate'
        kwargUniverse = ['fcid', 'subNum', 'Sample', 'refGenome', 'aligner',
                         'readOneOnly']
        self.dct['Sample'] = sorted(self.dct.pop('sampList'))
        self.getSubNum(**kwargs)
        self.checkForRefGenome(**kwargs)
        self.filterDct(kwargUniverse)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def getSubNum(self, **kwargs):
        """ Assumes self.dct['Sample'] is already sorted. The attachment of
        the Sub_ or sub_ prefix happens in other places, might want to
        reconsider. """
        if 'subNumPrefix' in kwargs and kwargs['subNumPrefix']:
            subNum = kwargs['subNumPrefix']
        else:
            assert kwargs['subNum']
            subNum = kwargs['subNum']
        hsh = hashlib.md5(str(self.dct['Sample'])).hexdigest()[:5]
        self.dct['subNum'] = '{}_{}'.format(subNum, hsh)


class FilterCalcTuple(_CalcTuple):

    def buildCalcTup(self, **kwargs):
        kwargUniverse = ['Sample', 'readOneOnly', 'mismatch']
        self.dct = {
            k: v for k,
            v in self.dct.items() if k in kwargUniverse and v}
        self.node = 'Filter'
        assert 'Sample' in kwargs, "Need Sample for Filter ct"
        self.check_if_readOneOnly_should_be_popped()
        if 'mismatch' in kwargs:
            self.popMismatchIfDefault(**kwargs)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def check_if_readOneOnly_should_be_popped(self):
        if self.delete_readOneOnly_if_false():
            return
        samp = self.dct['Sample']
        fcid = self.db.getAttrFromSamp('fcid', samp)
        ln = self.db.getAttrFromSamp('flowcell_lane', samp)
        bcl_ct = BclCalcTuple(db=self.db, fcid=fcid, laneNum=ln)
        mask = bcl_ct.buildMaskString()
        is_single_end_lane = (len(mask.split(',')) == 2)
        if is_single_end_lane:
            del self.dct['readOneOnly']

    def getUpstream(self, node=None):
        assert(node is None), "Haven't written logic for arb. upstrm node yet."
        dct = self.dct.copy()
        samp = dct.pop('Sample')
        dct['fcid'] = self.db.getAttrFromSamp('fcid', samp)
        dct['laneNum'] = self.db.getAttrFromSamp('flowcell_lane', samp)
        doReadRescue = self.db.getIndexTypeFromSamp(self.dct['Sample']) == 'LM'
        upNode = ['Bcl2fastq', 'ReadRescue'][doReadRescue]
        return CalcTuple(db=self.db, node=upNode, **dct)


class PoolCalcTuple(_CalcTuple):

    def buildCalcTup(self, **kwargs):
        """
            Assume all samples come from same flowcell
            @param kwargs - calculation specific parameters
        """
        kwargUniverse = ['fcid', 'poolId', 'Sample', 'mismatch', 'readOneOnly']
        self.dct = {
            k: v for k,
            v in self.dct.items() if k in kwargUniverse and v}
        self.node = 'Pool'

        assert 'Sample' in self.dct, "Need Samples for Pool ct"
        if isinstance(self.dct['Sample'], str):
            samps = [self.dct['Sample']]
        else:
            samps = self.dct['Sample']

        if 'poolId' not in self.dct:
            poolId = list(set([samp.split('_')[0] for samp in samps]))
            assert len(poolId) == 1, "no poolId given, but Samples have "\
                "different prefixes: {}".format(Sample)
            self.dct['poolId'] = poolId[0]

        if 'fcid' not in self.dct:
            # Bold move here: if there are samples from different fcids, I'm
            # assuming they're sorted (good assumption) and that the largest
            # accession number is going to belong to latest fcid (iffy) and
            # we are working form latest possible fcid. (also iffy). If user
            # wants it done right, they should be inputting an fcid in kwargs.
            fcid = [self.db.getAttrFromSamp('fcid', samp)
                    for samp in samps][-1]

        if 'mismatch' in kwargs:
            self.popMismatchIfDefault(**kwargs)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def getUpstream(self, node=None):
        """
            return previous calculation
            node doesn't seem to be used, so why is it there
            come back
        """
        dct = self.dct.copy()
        samps = dct.pop('Sample')
        if isinstance(samps, str):
            samps = [samps]
        return [CalcTuple(
            db=self.db, node='Filter', Sample=samp, **dct)
            for samp in samps]

    def getErrors(self):
        superErr = super(PoolCalcTuple, self).getErrors()
        if superErr:
            return superErr
        info = self.db.getCalc(self)
        if not info['outdir'] or not os.path.exists(info['outdir']):
            return "No outdir found: {}".format(info)
        sampleFfqFiles = glob.glob(info['outdir'] + '/*')
        if len(sampleFfqFiles) < 1:
            return "No files in pool directory: {}".format(info['outdir'])
        for file in sampleFfqFiles:
            f = os.path.join(info['outdir'], file)
            if (os.path.islink(f) and not os.path.exists(os.readlink(f))):
                return "Dead Link in pool dir {}".format(f)
        return None


class ReadRescCalcTuple(_CalcTuple):

    def buildCalcTup(self, **kwargs):
        kwargUniverse = ['fcid', 'laneNum', 'mismatch']
        self.dct = {
            k: v for k,
            v in self.dct.items() if k in kwargUniverse and v}
        self.node = 'ReadRescue'
        assert "laneNum" in kwargs, "Need laneNum for RR"
        assert "fcid" in kwargs, "Need fcid for RR"
        assert self.db.getIndexTypeFromFcidLn(
            fcid=kwargs['fcid'], laneNum=kwargs['laneNum']) == 'LM', \
            'Read Rescue only implemented for LM lanes'
        if 'mismatch' in kwargs:
            self.popMismatchIfDefault(**kwargs)
        self.buildCalcTupFromDctAndOrderedKwargUniv(kwargUniverse)

    def getUpstream(self, node=None):
        if node:
            assert node == 'Bcl2fastq', 'Bcl2fastq only upstream node, {} \
                    is invalid.'.format(node)
        return CalcTuple(db=self.db, node='Bcl2fastq', **self.dct)

    def getErrors(self):
        superErr = super(ReadRescCalcTuple, self).getErrors()
        if superErr:
            return superErr
        info = self.getMetadata()
        sampleDirs = glob.glob(info['outdir'] + '/*/*')
        for sampDir in sampleDirs:
            files = os.listdir(sampDir)
            if len(files) < 2:
                return "No Links in RR dir"
            for file in files:
                f = os.path.join(sampDir, file)
                if (os.path.islink(f) and
                        not os.path.exists(os.readlink(f))):
                    return "Dead Link in RR dir"
        indx_type = self.db.getIndexTypeFromFcidLn(**self.dct)
        if indx_type != 'LM':
            return "{} is not LM, no Read Rescue avail.".format(indx_type)
        indx_type = self.db.getIndexTypeFromFcidLn(**self.dct)
        if indx_type != 'LM':
            return "{} is not LM, no Read Rescue avail.".format(indx_type)
        return None

    def getSampOutdir(self, samp):
        laneDir = self.getMetadata('outdir')
        return os.path.join(laneDir,
                            self.getUpstream().getSampRelOutdir(samp))


def getPoolDictFromSamps(samps):
    """ We generate a default Pool dictionary as follows: Recall we define the
    sample id string as an underscore-separated compound of the
    flowcell_sample_id followed by accession (seq_status terms). So all samples
    which share the prefix (flowcell_sample_id) are collected, and the new
    sample is just the prefix.

    For example, for the samples 1_1, 1_2, and 1_3, the key is '1' and the
    value is a list: ['1_1', 1_2' 1_3']

        @param samps - list of samples
    """

    ids = sorted(samps)  # this asssures each value list is sorted, so we can
    # later assure a unique representation for hashing and repeatability.
    poolDict = {_id: sorted(dups) for (_id, dups) in
                groupby(ids, lambda x: x.split('_')[0])}
    # Two possible approaches: 1) only include duplicate samples in pool, then
    # later make sure we find the solos, or 2) have solos here. I prefer 1.
    poolDict = {k: v for (k, v) in poolDict.items() if len(v) > 1}
    return poolDict
