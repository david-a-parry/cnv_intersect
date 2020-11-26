import logging
import pysam
import re
from collections import defaultdict
from cnv_intersect.cnv import Cnv

valid_cnv_types = ['LOSS', 'GAIN']
autosome_re = re.compile(r"""^(chr)?(\d+)$""")
warned_chroms = set()

logger = logging.getLogger("CNV Intersect")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)


def warn_non_standard_chromosome(chrom):
    if chrom not in warned_chroms:
        logger.warn("Ignoring non-standard chromosome '{}'".format(chrom))
        warned_chroms.add(chrom)


def cnv_type_from_record(record, sample, ploidy=2):
    ''' Return 'LOSS', 'GAIN' or 'REF' based on sample call in VcfRecord'''
    if 'SVTYPE' not in record.info:
        return 'REF'
    if record.info['SVTYPE'] == 'CNV':
        if record.samples[sample]['CN'] < ploidy:
            return 'LOSS'
        elif record.samples[sample]['CN'] > ploidy:
            return 'GAIN'
    else:
        if len(record.alts) > 1:
            raise ValueError("Variants must be biallelic but the variant " +
                             "{}:{}-{}/{}".format(record.chrom,
                                                  record.pos,
                                                  record.ref,
                                                  ",".join(record.alts)) +
                             " has {} ALT alleles".format(len(record.alts)))
        if 1 in record.samples[sample]['GT']:
            if record.info['SVTYPE'] == 'DEL':
                return 'LOSS'
            elif record.info['SVTYPE'] == 'DUP':
                return 'GAIN'
            elif record.info['SVTYPE'] == 'INS':
                return 'GAIN'
            return record.info['SVTYPE']
    return 'REF'


class CnvFromVcf(Cnv):
    ''' A contiguous CNV that may be made up of >1 VCF records.'''

    def __init__(self, cnv_type, records, var_samples=[]):
        self.cnv_type = cnv_type
        self.records = records
        super().__init__(chrom=records[0].chrom,
                         start=records[0].start,
                         stop=records[-1].stop,
                         cnv_type=cnv_type,
                         records=records)
        self.samples = self._copy_numbers_from_records()
        self.var_samples = var_samples
        self.__n_records = None

    @property
    def n_records(self):
        if self.__n_records is None:
            self.__n_records = len(self.records)
        return self.__n_records

    def _copy_numbers_from_records(self):
        smp2cn = defaultdict(dict)
        for s in self.records[0].samples:
            copy_numbers = [x.samples[s]['CN'] for x in self.records]
            smp2cn[s]['mean_copy_number'] = sum(copy_numbers)/len(copy_numbers)
            smp2cn[s]['max_copy_number'] = max(copy_numbers)
            smp2cn[s]['min_copy_number'] = min(copy_numbers)
            smp2cn[s]['total_bin_counts'] = sum(x.samples[s]['BC'] for x in
                                                self.records)
        return smp2cn


class CnvVcf(object):
    ''' A class for iterating over CNVs of the same type in a VCF.'''

    def __init__(self, vcf, cnv_type='LOSS', ped=None, pass_filters=False,
                 minimum_length=100):
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

                  pass_filters:
                       Only include variants with PASS in the FILTER field.

                  minimum_length:
                       Ignore variants shorter than this value.

        '''
        if cnv_type.upper() not in valid_cnv_types:
            raise ValueError("cnv_type must be either 'LOSS' or 'GAIN'")
        self.vcf = pysam.VariantFile(vcf)
        self.cnv_type = cnv_type.upper()
        self.ped = ped
        self.pass_filters = pass_filters
        self.minimum_length = minimum_length
        self.current_cnv = None
        self.next_record = None
        self.buffer = []
        if self.ped:
            self.affected = [x for x in self.vcf.header.samples if x in
                             self.ped.get_affected()]
            self.unaffected = [x for x in self.vcf.header.samples if x in
                               self.ped.get_unaffected()]
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
                    if not self.record_matches_type(record):
                        continue
                    if (record.start < self.buffer[-1].stop and
                            record.chrom == self.buffer[-1].chrom):
                        self.buffer.append(record)
                    else:
                        self.next_record = record
                        break
        except StopIteration:
            pass
        if not self.buffer:
            raise StopIteration()
        return CnvFromVcf(cnv_type=self.cnv_type,
                          records=self.buffer,
                          var_samples=self.affected)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.vcf.close()

    def record_matches_type(self, record):
        '''
            Return True if for every affected call the sample has copy
            number call matching self.cnv_type and passes filters if set.
        '''
        if self.pass_filters and 'PASS' not in record.filter:
            return False
        if (self.minimum_length > 0 and
                record.stop - record.start < self.minimum_length):
            return False
        if autosome_re.match(record.chrom):
            ploidies = [2] * len(self.affected)
            un_ploidies = [2] * len(self.unaffected)
        else:
            if not self.ped:
                return False  # no genders, ignore non-autosomes
            elif record.chrom == 'chrX' or record.chrom == 'Y':
                ploidies = [2 if self.ped.individuals[s].is_female() else 1
                            for s in self.affected]
                un_ploidies = [2 if self.ped.individuals[s].is_female() else 1
                               for s in self.unaffected]
            elif record.chrom == 'chrY' or record.chrom == 'Y':
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
                if self._unaffected_with_same_ploidy(record):
                    return False
            return True
        return False

    def _unaffected_with_same_ploidy(self, record):
        if record.info['SVTYPE'] == 'CNV':
            return self._unaffected_with_same_cnv(record)
        return self._unaffected_with_same_sv(record)

    def _unaffected_with_same_cnv(self, record):
        if self.cnv_type == 'LOSS':
            max_cn_aff = max(record.samples[s]['CN'] for s in
                             self.affected)
            min_cn_un = min(record.samples[u]['CN'] for u in
                            self.unaffected)
            return max_cn_aff >= min_cn_un
        elif self.cnv_type == 'GAIN':
            min_cn_aff = min(record.samples[s]['CN'] for s in
                             self.affected)
            max_cn_un = max(record.samples[u]['CN'] for u in
                            self.unaffected)
            return min_cn_aff <= max_cn_un

    def _unaffected_with_same_sv(self, record):
        # should already have checked that this is biallelic in
        # cnv_type_from_record function
        min_aff = min(record.samples[s]['GT'].count(1) for s in self.affected)
        max_un = min(record.samples[u]['GT'].count(1) for u in self.unaffected)
        return min_aff > max_un
