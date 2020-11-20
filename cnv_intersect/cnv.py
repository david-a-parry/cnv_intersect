class Cnv(object):
    ''' A contiguous CNV that may be made up of >1 VCF/BED records.'''

    def __init__(self, chrom, start, stop, cnv_type, records=[]):
        self.cnv_type = cnv_type
        self.records = records
        self.chrom = chrom
        self.start = int(start)
        self.stop = int(stop)
        assert(self.start < self.stop)

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

    def merge_interval(self, other):
        '''
            Merge an overlapping interval.

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


class NonOverlappingIntervalError(ValueError):
    pass
