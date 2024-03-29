import logging
import pysam
import re
from collections import defaultdict
from cnv_intersect.cnv import Cnv
from vase.info_filter import InfoFilter
from vase.format_filter import FormatFilter

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
    return dict((k, n) for n, k in enumerate(contigs))


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
                         start=min(x.start for x in records),
                         stop=max(x.stop for x in records),
                         cnv_type=cnv_type,
                         records=records)
        self.samples = self._copy_numbers_from_records()
        self.var_samples = var_samples
        self.__n_records = None

    def __eq__(self, other):
        return (other is not None and self.cnv_type == other.cnv_type and
                self.chrom == other.chrom and self.start == other.start and
                self.stop == other.stop and
                len(self.records) == len(other.records) and
                set(self.var_samples) == set(other.var_samples))

    def __ne__(self, other):
        return (other is None or self.cnv_type != other.cnv_type or
                self.chrom != other.chrom or self.start != other.start or
                self.stop != other.stop and
                len(self.records) != len(other.records) and
                set(self.var_samples) != set(other.var_samples))

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return "{}:{}-{}-{}-{} ({} records)".format(self.chrom,
                                                    self.start,
                                                    self.stop,
                                                    self.cnv_type,
                                                    ",".join(self.var_samples),
                                                    len(self.records))

    @property
    def n_records(self):
        if self.__n_records is None:
            self.__n_records = len(self.records)
        return self.__n_records

    def _copy_numbers_from_records(self):
        smp2cn = defaultdict(dict)
        for s in self.records[0].samples:
            copy_numbers = [self._copy_number_from_record(x, s) for x in
                            self.records]
            smp2cn[s]['mean_copy_number'] = sum(copy_numbers)/len(copy_numbers)
            smp2cn[s]['max_copy_number'] = max(copy_numbers)
            smp2cn[s]['min_copy_number'] = min(copy_numbers)
            if 'BC' in self.records[0].samples[s]:
                smp2cn[s]['total_bin_counts'] = sum(x.samples[s]['BC'] for x in
                                                    self.records)
            else:
                smp2cn[s]['total_bin_counts'] = None
            for qs in ['QS', 'GQ']:
                if all(qs in self.records[i].samples[s] for i in
                       range(len(self.records))):
                    smp2cn[s]['mean_quality_score'] = sum(
                        x.samples[s][qs] for x in
                        self.records)/len(self.records)
                    break
                smp2cn[s]['mean_quality_score'] = None
        return smp2cn

    def _copy_number_from_record(self, record, sample, ploidy=2):
        if record.info['SVTYPE'] == 'CNV':
            return record.samples[sample]['CN']
        else:
            alts = sum(1 for x in record.samples[sample]['GT'] if x is not None
                       and x > 0)
            if record.info['SVTYPE'] == 'DEL':
                return (alts * -1) + ploidy
            elif record.info['SVTYPE'] == 'DUP':
                return alts + ploidy
            elif record.info['SVTYPE'] == 'INS':
                return alts + ploidy
            raise ValueError("Can not interpret SNVTYPE={}".format(
                record.info['SVTYPE']))


class CnvVcf(object):
    ''' A class for iterating over CNVs of the same type in a VCF.'''

    def __init__(self, vcf, cnv_type='LOSS', ped=None, pass_filters=False,
                 minimum_length=None, maximum_length=None, cnv_filters=[],
                 overlap_fraction=0.8, info_filters=None, sample_filters=None):
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

                  cnv_filters:
                       Collection of CnvBed objects to use to filter
                       matching variants. If any of these contain overlapping
                       CNVs of the same type as a record in the record will be
                       ignored if covered by a fraction defined by the
                       'overlap_fraction' option.

                 overlap_fraction:
                       Minimum overlap of a CNV from cnv_filters with a record
                       required in order to filter the record. Default=0.8.

                 info_filters:
                        Custom filter expressions for filtering on fields in
                        the INFO field of each record. Must be in the format
                        '<INFO_FIELD> <comparator> <value>'. Variants will be
                        retained if they meet the given criteria. For example,
                        to only keep records with a QD score greater than 4,
                        you would pass the expression "QD > 4". To only
                        keep records with the "DB" flag present you would pass
                        the expression "DB == True".

                        Standard python style operators (">", "<", ">=", "<=",
                        "==", "!=") are supported. Comparisons will be
                        performed using the types specified for the given field
                        in the VCF header (e.g. Float, Integer or String) or as
                        booleans for Flags.

                sample_filters:
                        As above for 'info_filters' parameter except values are
                        evaluated against individual sample calls instead. If a
                        ped is provided these expressions are only evaluated
                        against affected individuals.
        '''
        if cnv_type.upper() not in valid_cnv_types:
            raise ValueError("cnv_type must be either 'LOSS' or 'GAIN'")
        self.vcf = pysam.VariantFile(vcf)
        self.cnv_type = cnv_type.upper()
        self.ped = ped
        self.pass_filters = pass_filters
        self.minimum_length = minimum_length
        self.maximum_length = maximum_length
        self.cnv_filters = cnv_filters
        self.overlap_fraction = overlap_fraction
        self.current_cnv = None
        self.next_record = None
        self.records_read = 0
        self.records_filtered = 0
        self.cnvs_created = 0
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
        self.info_filter = None
        self.sample_filter = None
        if info_filters:
            self.info_filter = self._parse_filter_expressions(
                filters=info_filters,
                filter_class=InfoFilter)
        if sample_filters:
            self.sample_filter = self._parse_filter_expressions(
                filters=sample_filters,
                filter_class=FormatFilter)

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
                    self.records_read += 1
                    if self.record_matches_type(record):
                        self.buffer.append(record)
                        break
            # add contiguous records to buffer
            if self.buffer:
                for record in self.vcf:
                    self.records_read += 1
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
        self.cnvs_created += 1
        return CnvFromVcf(cnv_type=self.cnv_type,
                          records=self.buffer,
                          var_samples=self.affected)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.vcf.close()

    def _parse_filter_expressions(self, filters, filter_class):
        split_filters = []
        for expression in filters:
            exp = expression.split()
            if len(exp) < 3 or len(exp) > 4:
                raise ValueError(
                    "--info_filters/--sample_filters must consist of three " +
                    "or four quoted values separated by whitespace - for " +
                    "example: 'QS > 20' or 'DHFFC < 0.7 LOSS'. The provided " +
                    "expression '{}' is invalid.".format(expression))
            if len(exp) == 3:
                split_filters.append(exp)
            else:
                if exp[3] not in valid_cnv_types:
                    raise ValueError(
                        "CNV type '{}' in filter expression ".format(exp[3]) +
                        "'{}' not recognised. ".format(expression) +
                        "Valid CNV values are 'LOSS' or 'GAIN'.")
                if exp[3] == self.cnv_type:
                    split_filters.append(exp[:3])
        if split_filters:
            return filter_class(vcf=self.vcf, filters=split_filters)
        return None

    def _filter_record(self):
        self.records_filtered += 1
        return False

    def record_matches_type(self, record):
        '''
            Return True if for every affected call the sample has copy
            number call matching self.cnv_type and passes filters if set.
        '''
        if self.pass_filters and 'PASS' not in record.filter:
            return self._filter_record()
        if (self.minimum_length and
                record.stop - record.start < self.minimum_length):
            return self._filter_record()
        if (self.maximum_length and
                record.stop - record.start > self.maximum_length):
            return self._filter_record()
        if autosome_re.match(record.chrom):
            ploidies = [2] * len(self.affected)
            un_ploidies = [2] * len(self.unaffected)
        else:
            if not self.ped:
                return self._filter_record()  # no genders, ignore non-autosomes
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
                return self._filter_record()  # ignore non-standard chromosomes
        if all(cnv_type_from_record(record, s, p) == self.cnv_type
               for s, p in zip(self.affected, ploidies)):
            if any(cnv_type_from_record(record, u) == self.cnv_type
                   for u, p in zip(self.unaffected, un_ploidies)):
                if self._unaffected_with_same_ploidy(record):
                    return self._filter_record()
            if self.info_filter is not None:
                if all(self.info_filter.filter(record)):
                    return self._filter_record()
            if self.sample_filter is not None:
                if all(all(self.sample_filter.filter(record, x)) for x in
                       self.affected):
                    # ALL affected samples fail - should this be any?
                    return self._filter_record()
            if self._cnv_filter_matches(record):
                return self._filter_record()
            return True
        return self._filter_record()

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

    def _cnv_filter_matches(self, record):
        if not self.cnv_filters:
            return False
        cnv = Cnv(chrom=record.chrom,
                  start=record.start,
                  stop=record.stop,
                  cnv_type=self.cnv_type,
                  records=[record])
        for cf in self.cnv_filters:
            overlaps = cf.search(cnv)
            for ovr in overlaps:
                o_start = max(ovr.start, record.start)
                o_stop = min(ovr.stop, record.stop)
                covered = float(o_stop - o_start)/(record.stop - record.start)
                if covered > self.overlap_fraction:
                    return True
        return False
