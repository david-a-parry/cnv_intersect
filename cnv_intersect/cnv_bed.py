import gzip
import csv
from collections import defaultdict
from cnv_intersect.cnv import Cnv
from operator import attrgetter

bed_cols = {'dbVar': ['chrom', 'start', 'end', 'name', 'score', 'strand',
                      'thickStart', 'thickEnd', 'reserved', 'frequency',
                      'type', 'length', 'label', 'freq_range', 'call_list'],
            'DGV': ['chrom', 'start', 'end', 'name', 'score', 'strand',
                    'thickStart', 'thickEnd', 'itemRgb', 'type', 'reference',
                    'pubMedId', 'method', 'platform', 'mergedVariants',
                    'supportingVariants', 'sampleSize', 'observedGains',
                    'observedLosses', 'cohortDescription', 'genes', 'samples',
                    '_size']}


class CnvBed(object):
    '''
    Hold CNV information from BED files in memory. Currently supports BED
    files from UCSC's "NCBI dbVar Curated Common Structural Variants"
    generated as follows:

        rsync -a -P rsync://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/dbVar ./
        for BB in dbVar/*.bb
        do
            BED=$(dirname $BB)$(basename $BB .bb).bed
            bigBedToBed $BB $BED
        done

    '''

    def __init__(self, bed, bed_format='dbVar'):
        self.filename = bed
        self.cnvs = self.read_bed(bed, bed_format)

    def read_bed(self, f, bed_format):
        regions = defaultdict(list)
        if bed_format not in bed_cols:
            raise ValueError("Unrecognized BED format '{}'".format(bed_format))
        if f.endswith('.gz'):
            o_func = gzip.open
        else:
            o_func = open
        with o_func(f, 'rt') as fh:
            bed = csv.DictReader(fh,
                                 delimiter='\t',
                                 fieldnames=bed_cols[bed_format])
            for row in bed:
                if 'deletion' in row['type'] \
                   or row['type'] == 'copy number loss':
                    cnv_type = ['LOSS']
                elif row['type'] == 'duplication' \
                   or row['type'] == 'copy number gain':
                    cnv_type = ['GAIN']
                elif row['type'] == 'copy number variation':
                    cnv_type = ['LOSS', 'GAIN']  # TODO is this correct?!
                for ct in cnv_type:
                    cnv = Cnv(chrom=row['chrom'],
                              start=row['start'],
                              stop=row['end'],
                              cnv_type=ct,
                              records=[row])
                    regions[row['chrom']].append(cnv)
        for r in regions.values():
            r.sort(key=attrgetter('start', 'stop'))
        return regions

    def search(self, cnv, match_type=True):
        ''' Return list of CNVs overlapping Cnv object'''
        if cnv.chrom not in self.cnvs:
            return []
        i = self._bin_search_cnvs(cnv)
        matches = []
        if i > -1:
            for j in range(i - 1, -1, -1):
                other = self.cnvs[cnv.chrom][j]
                if cnv.overlaps(other):
                    if not match_type or cnv.cnv_type == other.cnv_type:
                        matches.insert(0, other)
                else:
                    break
            for j in range(i, len(self.cnvs[cnv.chrom])):
                other = self.cnvs[cnv.chrom][j]
                if cnv.overlaps(other):
                    if not match_type or cnv.cnv_type == other.cnv_type:
                        matches.append(other)
                else:
                    break
        return matches

    def _bin_search_cnvs(self, cnv):
        ''' Return array index of first overlapping cnv found in self.cnvs'''
        others = self.cnvs[cnv.chrom]
        low = 0
        high = len(others) - 1
        while low <= high:
            i = (low + high) // 2
            if cnv.overlaps(others[i]):
                return i
            elif cnv < others[i]:
                high = i - 1
            elif cnv > others[i]:
                low = i + 1
        return -1

    def walk(self, chrom, start, stop, cnv_type):
        '''
        Search for CNVs overlapping coordinates. Designed to be run
        sequentially for multiple look-ups performed in coordinate order.
        '''
        raise NotImplementedError()
