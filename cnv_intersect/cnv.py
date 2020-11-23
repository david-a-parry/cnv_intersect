from operator import attrgetter


class Cnv(object):
    ''' A contiguous CNV that may be made up of >1 VCF/BED records.'''

    def __init__(self, chrom, start, stop, cnv_type, records=[]):
        self.cnv_type = cnv_type
        self.records = records
        self.chrom = chrom
        self.start = int(start)  # 0-based
        self.stop = int(stop)    # 1-based
        assert(self.start < self.stop)

    @property
    def length(self):
        return self.stop - self.start

    def __copy__(self):
        return Cnv(self.chrom, self.start, self.stop, self.cnv_type,
                   self.records)

    def __str__(self):
        return "{}:{}-{}-{} ({} records)".format(self.chrom,
                                                 self.start,
                                                 self.stop,
                                                 self.cnv_type,
                                                 len(self.records))

    def __eq__(self, other):
        return (other is not None and self.cnv_type == other.cnv_type and
                self.chrom == other.chrom and self.start == other.start and
                self.stop == other.stop)

    def __ne__(self, other):
        return (other is None or self.cnv_type != other.cnv_type or
                self.chrom != other.chrom or self.start != other.start or
                self.stop != other.stop)

    def __lt__(self, other):
        return (self.chrom < other.chrom or
                (self.chrom == other.chrom and
                 (self.start < other.start or self.start == other.start and
                  self.stop < other.stop)))

    def __le__(self, other):
        return (self.chrom < other.chrom or
                (self.chrom == other.chrom and self.stop <= other.stop))

    def __gt__(self, other):
        return (self.chrom > other.chrom or
                (self.chrom == other.chrom and
                 (self.start > other.start or self.start == other.start and
                  self.stop > other.stop)))

    def __ge__(self, other):
        return (self.chrom > other.chrom or
                (self.chrom == other.chrom and self.start >= other.start))

    def overlaps(self, other):
        if self.chrom != other.chrom:
            return False
        elif self.start < other.stop and self.stop > other.start:
            # start is 0-based so only overlaps if < stop
            return True
        return False

    def contains(self, other):
        if self.chrom != other.chrom:
            return False
        if self.start <= other.start and self.stop >= other.stop:
            return True
        return False

    def merge_cnv(self, other):
        '''
            Merge an overlapping CNV.

            Args:
                other:   Another Cnv object

        '''
        if not self.overlaps(other):
            raise NonOverlappingIntervalError("Can not merge non-overlapping" +
                                              " intervals '{}' and '{}'"
                                              .format(self, other))
        if other.start < self.start:
            self.start = other.start
        if other.stop > self.stop:
            self.stop = other.stop
        self.records.extend(other.records)

    def merge_new(self, others, assume_sorted=True):
        '''
        Return a new merged Cnv object where each record is a reference to
        the individual merged Cnvs. Additional CNVs must be provided as a
        list sorted in coordinate order unless assume_sorted is False.

        Args:
            others:
                A list of additional CNVs. Must be sorted in coordinate
                order unless assume_sorted is False.

            assume_sorted:
                Assume others list is already sorted in coordinate order.
                Avoids performing a sort on the list but if the list is
                not already sorted this will likely result in a
                NonOverlappingIntervalError.
        '''
        if not assume_sorted:
            others.sort(key=attrgetter('start', 'stop'))
        merged = Cnv(self.chrom,
                     self.start,
                     self.stop,
                     self.cnv_type,
                     [self] + others)
        for other in others:
            if not merged.overlaps(other):
                raise NonOverlappingIntervalError(
                    "Can not merge non-overlapping intervals '{}' and '{}'"
                    .format(self, other))
            if other.start < merged.start:
                merged.start = other.start
            if other.stop > merged.stop:
                merged.stop = other.stop
        return merged


class NonOverlappingIntervalError(ValueError):
    pass
