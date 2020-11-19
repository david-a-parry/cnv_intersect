#!/usr/bin/env python3
import sys
import argparse
import logging
import re
import pysam
from vase.ped_file import PedFile

logger = logging.getLogger("CNV Intersect")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)

valid_cnv_types = ['LOSS', 'GAIN']
autosome_re = re.compile(r"""^(chr)?(\d+)$""")
warned_chroms = set()


def warn_non_standard_chromosome(chrom):
    if chrom not in warned_chroms:
        logger.warn("Ignoring non-standard chromosome '{}'".format(chrom))
        warned_chroms.add(chrom)


def cnv_type_from_record(record, sample, ploidy=2):
    ''' Return 'LOSS', 'GAIN' or 'REF' based on sample call in VcfRecord'''
    if record.samples[sample]['CN'] < ploidy:
        return 'LOSS'
    elif record.samples[sample]['CN'] > ploidy:
        return 'GAIN'
    return 'REF'


class Cnv(object):
    ''' A contiguous CNV that may be made up of >1 VCF records.'''

    def __init__(self, cnv_type, records):
        self.cnv_type = cnv_type
        self.records = records
        self.chrom = records[0].chrom
        self.start = records[0].start
        self.stop = records[-1].stop
        assert(self.start < self.stop)

    def __str__(self):
        return "{}:{}-{}-{} ({} records)".format(self.chrom,
                                                 self.start,
                                                 self.stop,
                                                 self.cnv_type,
                                                 len(self.records))


class CnvVcf(object):
    ''' A class for iterating over CNVs of the same type in a VCF where '''

    def __init__(self, vcf, cnv_type='LOSS', ped=None):
        '''
            Args:
                  vcf: path to VCF/BCF file

                  cnv_type:
                       Type of CNV to check for. Must be one of 'LOSS'
                       or 'GAIN'.

                  ped: Optional PedFile object for assigning affected and
                       unaffected individuals from VCF and for
                       determining gender. If None all samples from VCF
                       will be considered to be affected and all calls
                       from sex chromosomes will be ignored.

        '''
        if cnv_type.upper() not in valid_cnv_types:
            raise ValueError("cnv_type must be either 'LOSS' or 'GAIN'")
        self.vcf = pysam.VariantFile(vcf)
        self.cnv_type = cnv_type.upper()
        self.ped = ped
        self.current_cnv = None
        self.next_record = None
        self.buffer = []
        if self.ped:
            self.affected = [x for x in self.header.samples if x in
                             self.ped.get_affected]
            self.unaffected = [x for x in self.header.samples if x in
                               self.ped.get_unaffected]
            if not self.affected:
                raise ValueError("None of the specified affected samples " +
                                 "from {} were found in VCF {}".format(
                                     self.ped.filename, vcf))
        else:
            self.affected = [x for x in self.vcf.header.samples]
            if not self.affected:
                raise ValueError("No samples found in VCF {}".format(vcf))
            self.unaffected = []

    def __iter__(self):
        return self

    def __next__(self):
        self.buffer = []
        try:
            # get the next variant that matches self.cnv_type
            if self.next_record is not None:
                if self.record_matches_type(self.next_record):
                    self.buffer.append(self.next_record)
            self.next_record = None
            if not self.buffer:
                for record in self.vcf:
                    if self.record_matches_type(record):
                        self.buffer.append(record)
                        break
            # add contiguous records to buffer
            if self.buffer:
                for record in self.vcf:
                    if record.chrom != self.buffer[-1].chrom:
                        break
                    if record.start < self.buffer[-1].stop and \
                            self.record_matches_type(record):
                        self.buffer.append(record)
                    else:
                        break
            self.next_record = record
        except StopIteration:
            pass
        if not self.buffer:
            raise StopIteration()
        return Cnv(cnv_type=self.cnv_type,
                   records=self.buffer)

    def record_matches_type(self, record):
        '''
            Return True if for every affected call the sample has copy
            number call matching self.cnv_type.
        '''
        if autosome_re.match(record.chrom):
            ploidies = [2] * len(self.affected)
            un_ploidies = [2] * len(self.unaffected)
        else:
            if not self.ped:
                return False  # no genders, ignore non-autosomes
            elif record.chrom == 'chrX':
                ploidies = [2 if self.ped.individuals[s].is_female() else 1
                            for s in self.affected]
                un_ploidies = [2 if self.ped.individuals[s].is_female() else 1
                               for s in self.unaffected]
            elif record.chrom == 'chrY':
                ploidies = [1 if self.ped.individuals[s].is_male() else 0
                            for s in self.affected]
                un_ploidies = [1 if self.ped.individuals[s].is_male() else 0
                               for s in self.unaffected]
            else:
                warn_non_standard_chromosome(record.chrom)
                return False  # ignore non-standard chromosomes
        if all(cnv_type_from_record(record, s, p) == self.cnv_type
               for s, p in zip(self.affected, ploidies)):
            if any(cnv_type_from_record(record, u) == self.cnv_type
                   for u, p in zip(self.unaffected, un_ploidies)):
                max_cn_aff = max(record[s]['CN'] for s in self.affected)
                min_cn_un = min(record[u]['CN'] for u in self.unaffected)
                if max_cn_aff >= min_cn_un:
                    return False
            return True
        return False


def get_argparser():
    parser = argparse.ArgumentParser(
        usage='%(prog)s [options] VCF VCF [VCF ...]',
        description='Generate VCF of overlapping CNVs from Canvas called VCF')
    parser.add_argument('vcfs', metavar='VCF', nargs='+',
                        help='2 or more input VCF files with CNV calls.')
    parser.add_argument('-p', '--ped', help='''PED file indicating which
                        individuals are affected. If not provided all
                        individuals in given VCFs are assumed to be
                        affected.''')
    parser.add_argument('-o', '--output', metavar='VCF/BCF output',
                        help='''Output file for intersected calls. Default is
                        STDOUT.''')
    parser.add_argument('-m', '--min_intersect', default=2, type=int,
                        metavar='N', help='''Minimum number of intersecting
                        VCFs required in order to output CNV. Default=2.''')
    parser.add_argument('-t', '--cnv_types', default=['LOSS', 'GAIN'],
                        nargs='+', help='''CNV types to identify. Valid values
                        are "LOSS" or "GAIN".''')
    return parser


def get_contig_order(vcf):
    if vcf.index is not None:
        return dict(vcf.index)
    contigs = []
    prev_contig = None
    prev_pos = None
    not_sorted_err = ValueError("Input file '{}' is not sorted ".format(vcf) +
                                "- exiting!")
    for record in vcf:
        if record.chrom not in contigs:
            contigs.append(record.chrom)
        elif prev_contig is not None and record.chrom != prev_contig:
            raise not_sorted_err
        elif prev_pos is not None and prev_pos > record.pos:
            raise not_sorted_err
        prev_contig = record.chrom
        prev_pos = record.pos
    vcf.reset()
    return dict((k, n) for n, k in enumereate(contigs))


def check_contig_order(cnv_readers):
    contig_order = get_contig_order(cnv_readers[0].vcf)
    for cr in cnv_readers[1:]:
        if get_contig_order(cr.vcf) != contig_order:
            raise ValueError("Input VCFs do not share the same contig order!")
    return contig_order


def main(vcfs, ped=None, output=None, cnv_types=['LOSS', 'GAIN'],
         min_intersect=2):
    invalid = set(cnv_types).difference(valid_cnv_types)
    if invalid:
        raise ValueError("ERROR: Invalid CNV types specified: {}".format(
            ", ".join(invalid)))
    pedfile = None
    if ped is not None:
        pedfile = PedFile(ped)
    for cnv_type in cnv_types:
        cnv_readers = [CnvVcf(x, cnv_type=cnv_type, ped=pedfile) for x in vcfs]
        contig_order = check_contig_order(cnv_readers)
        cnvs = []
        for cr in cnv_readers:
            cnvs.append(next(cr))
        while True:
            if None in cnvs:
                remove = [i for i in range(len(cnvs)) if cnvs[i] is None]
                cnvs = [cnvs[i] for i in range(len(cnvs)) if i not in remove]
                cnv_readers = [cnv_readers[i] for i in range(len(cnv_readers))
                               if i not in remove]
                if len(cnvs) < min_intersect:
                    logger.info("Fewer than {} VCFs ".format(min_intersect) +
                                "with remaining {} calls. ".format(cnv_type) +
                                "Done parsing {} CNVs.".format(cnv_type))
                    break
            min_chrom = cnvs[0].chrom
            min_stop = cnvs[0].stop
            min_i = 0
            for i in range(1, len(cnvs)):
                if contig_order[cnvs[i].chrom] < contig_order[min_chrom]:
                    min_chrom = cnvs[i].chrom
                    min_stop = cnvs[i].stop
                    min_i = i
                elif cnvs[i].chrom == min_chrom and cnvs[i].stop < min_stop:
                    min_stop = cnvs[i].stop
                    min_i = i
            others = (cnvs[i] for i in range(len(cnvs)) if i != min_i)
            overlaps = [x for x in others if x.chrom == min_chrom and
                        x.start < min_stop]
            if overlaps:
                overlaps.append(cnvs[min_i])
                if len(overlaps) >= min_intersect:
                    start = min(x.start for x in overlaps)
                    stop = max(x.stop for x in overlaps)
                    logger.info("Overlapping {} at {}:{}-{} ({})".format(
                        cnv_type,
                        cnvs[0].chrom,
                        start,
                        stop,
                        len(overlaps)))
                #TODO Output CNV!
            for i in (x for x in range(len(cnvs)) if cnvs[x].stop == min_stop):
                try:
                    cnvs[i] = next(cnv_readers[i])
                except StopIteration:
                    logger.info("Exhausted records in {}".format(vcfs[i]))
                    cnvs[i] = None


if __name__ == '__main__':
    parser = get_argparser()
    args = parser.parse_args()
    if len(args.vcfs) < 2:
        parser.print_usage()
        sys.exit('At least 2 VCFs must be provided.')
    main(**vars(args))
