from cnv_intersect.cnv import Cnv


class CnvIntersect(Cnv):
    '''
    Class for intersection of multiple CNVs. Assumes all provided CNVs
    overlap. The start and stop attributes refer to the intersecting
    region (i.e. the region overlapping in all given cnvs). The
    limit_start and limit_stop properties refer to the 5' and 3' most
    coordinates of all CNVs.
    '''

    def __init__(self, cnvs):
        super().__init__(chrom=cnvs[0].chrom,
                         start=max(x.start for x in cnvs),
                         stop=min(x.stop for x in cnvs),
                         cnv_type=cnvs[0].cnv_type,
                         records=cnvs)
        self.limit_start = min(x.start for x in cnvs)
        self.limit_stop = max(x.stop for x in cnvs)
        self.n_cnvs = len(cnvs)

    @property
    def limit_length(self):
        return self.limit_stop - self.limit_start

    def limits_contain(self, other):
        '''
        Return True if the 5' and 3' limits of this CnvIntersect contain
        the limits of other.
        '''
        if self.chrom != other.chrom:
            return False
        if self.limit_start <= other.limit_start:
            if self.limit_stop >= other.limit_stop:
                return True
        return False
