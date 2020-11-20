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

    def search(self, chrom, start, stop, cnv_type):
        ''' Search for CNVs overlapping coordinates'''
        raise NotImplementedError()

    def walk(self, chrom, start, stop, cnv_type):
        '''
        Search for CNVs overlapping coordinates. Designed to be run
        sequentially for multiple lookups performed in coordinate order.
        '''
        raise NotImplementedError()
