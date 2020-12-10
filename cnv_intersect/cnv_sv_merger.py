import re
from cnv_intersect.cnv_vcf import get_contig_order, valid_cnv_types
from vase.info_filter import InfoFilter
from vase.format_filter import FormatFilter

chrom_re = re.compile(r'^(chr)?(\\d+|[XY])$')


class CnvSvMerger(object):
    '''
    Class for merging CNVs from Canvas and Manta VCFs from the same sample
    set.
    '''
    def __init__(self, canvas, manta, minimum_overlap=0.8, minimum_length=None,
                 maximum_length=None, pass_filters=False, canvas_filters=[],
                 manta_filters=[]):
        if sorted(canvas.header.samples) != sorted(manta.header.samples):
            raise ValueError("Samples differ in VCFs - Canvas and Manta VCFs" +
                             " must be from the same set of samples.")
        self.canvas = canvas
        self.manta = manta
        self.header = self._combine_headers()
        self.minimum_overlap = minimum_overlap
        self.minimum_length = minimum_length
        self.maximum_length = maximum_length
        self.pass_filters = pass_filters
        self.read_filters = self._parse_filters(canvas_filters, manta_filters)
        self.current_canvas_record = None
        self.current_manta_record = None
        self._current_chrom = None
        self._written_variants = []
        self.contig_order = self._get_contig_order()
        self._canvas_finished = False
        self._manta_finished = False

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_canvas_record is None and not self._canvas_finished:
            try:
                self.current_canvas_record = self._next_canvas_non_ref()
            except StopIteration:
                pass
        if self.current_manta_record is None and not self._manta_finished:
            try:
                self.current_manta_record = self._next_manta_cnv()
            except StopIteration:
                pass
        if self._canvas_finished and self._manta_finished:
            raise StopIteration
        if self.current_canvas_record is None:
            if self.current_manta_record is None:
                raise StopIteration
            record = self.current_manta_record
            vcf = self.canvas
        elif self.current_manta_record is None:
            record = self.current_canvas_record
            vcf = self.manta
        else:
            record, vcf = self._get_5prime_record()
        if record is self.current_manta_record:
            try:
                self.current_manta_record = self._next_manta_cnv()
            except StopIteration:
                self.current_manta_record = None
        else:
            try:
                self.current_canvas_record = self._next_canvas_non_ref()
            except StopIteration:
                self.current_canvas_record = None
        return self._overlaps_from_record(record,
                                          vcf)

    def _next_canvas_non_ref(self):
        for record in self.canvas:
            if self._record_passes_filters(record, 'canvas'):
                if record.id not in self._written_variants:
                    return record
        self._canvas_finished = True
        raise StopIteration

    def _next_manta_cnv(self):
        for record in self.manta:
            if self._record_passes_filters(record, 'manta'):
                if record.id not in self._written_variants:
                    return record
        self._manta_finished = True
        raise StopIteration

    def _record_passes_filters(self, record, genotyper):
        if record.alts is None:
            return False
        if self.pass_filters and 'PASS' not in record.filter:
            return False
        length = record.stop - record.start
        if self.minimum_length and length < self.minimum_length:
            return False
        if self.maximum_length and length > self.maximum_length:
            return False
        cnv_type = self._cnv_type(record)
        if cnv_type not in valid_cnv_types:
            return False
        info_filter = self.read_filters[genotyper]['info'][cnv_type]
        format_filter = self.read_filters[genotyper]['format'][cnv_type]
        if info_filter is not None and all(info_filter.filter(record)):
            return False
        if format_filter is not None and all(all(
                format_filter.filter(record, x)) for x in record.samples):
            # ALL samples fail format filtering
            return False
        return True

    def _get_5prime_record(self):
        canvas_chrom = self.contig_order[self.current_canvas_record.chrom]
        manta_chrom = self.contig_order[self.current_manta_record.chrom]
        canvas_start = self.current_canvas_record.start
        manta_start = self.current_manta_record.start
        if (canvas_chrom == manta_chrom and canvas_start <= manta_start) or \
           canvas_chrom < manta_chrom:
            return (self.current_canvas_record, self.manta)
        return (self.current_manta_record, self.canvas)

    def _cnv_type(self, record):
        if record.id is not None:
            if record.id.startswith("Canvas"):
                return record.id.split(':')[1]
            elif record.id.startswith("Manta"):
                cnv_type = record.id.split(':')[0].replace('Manta', '')
                if cnv_type == 'DEL':
                    return 'LOSS'
                elif cnv_type == 'DUP' or cnv_type == 'INS':
                    return 'GAIN'
                return 'Other'
        raise ValueError("Could not parse variant ID for record:\n" +
                         "{}\n. Was VCF created by Manta ".format(record) +
                         "or Canvas?")

    def _overlaps_from_record(self, record, other_vcf):
        if record.chrom != self._current_chrom:
            self._written_variants = []
        cnv_type = self._cnv_type(record)
        var_iter = other_vcf.fetch(record.chrom, max(0, record.start),
                                   record.stop)
        overlaps = []
        for other_rec in var_iter:
            if self._cnv_type(other_rec) != cnv_type:
                continue
            overlap_start = max(record.start, other_rec.start)
            overlap_stop = min(record.stop, other_rec.stop)
            rec_overlap = float(overlap_stop - overlap_start)/(record.stop -
                                                               record.start)
            oth_overlap = float(overlap_stop - overlap_start)/(other_rec.stop -
                                                               other_rec.start)
            if rec_overlap >= self.minimum_overlap \
               and oth_overlap >= self.minimum_overlap:
                if self._compatible_genotypes(record, other_rec):
                    overlaps.append(other_rec)
        return self._merge_records(record, overlaps)

    def _merge_records(self, record, overlaps):
        self._written_variants.append(record.id)
        if record.id.startswith('Canvas'):
            rec_source = 'Canvas'
        else:
            rec_source = 'Manta'
        if not overlaps:
            record.info['MergeSource'] = rec_source
            return self._new_record_from_template(record)
        to_merge = []
        # prefer Manta records with presumably more precise breakpoints
        if rec_source == 'Canvas':
            to_merge.append(record)
            output_record = overlaps[0]
            for x in overlaps[1:]:  # pick longest manta overlap as main record
                if (x.stop - x.start) > (output_record.stop -
                                         output_record.start):
                    to_merge.append(output_record)
                    output_record = x
                else:
                    to_merge.append(x)
        else:
            output_record = record
            to_merge = overlaps
        output_record = self._new_record_from_template(output_record)
        max_i = None
        max_length = None
        for i in range(len(to_merge)):
            # get index of longest Canvas record in to_merge
            if to_merge[i].id.startswith('Canvas'):
                length = to_merge[i].stop - to_merge[i].start
                if max_length is None or length > max_length:
                    max_length = length
                    max_i = i
        output_record.info['CanvasSupport'] = True
        output_record.info['MergeSource'] = 'Manta'
        alts = []
        start_diffs = []
        end_diffs = []
        for x in to_merge:
            output_record.id += ';' + x.id
            alts.append(x.alts[0])
            start_diffs.append(output_record.pos - x.pos)
            end_diffs.append(output_record.stop - x.stop)
            self._written_variants.append(x.id)
        output_record.info['MergeAlts'] = alts
        output_record.info['MergeStartDiff'] = start_diffs
        output_record.info['MergeEndDiff'] = end_diffs
        for s in output_record.samples:
            sep = '|' if to_merge[max_i].samples[s].phased else '/'
            output_record.samples[s]['CGT'] = sep.join(
                str(x) if x is not None else '.' for x in
                to_merge[max_i].samples[s]['GT'])
            if 'FT' in to_merge[max_i].format:
                output_record.samples[s]['CFT'] = \
                    to_merge[max_i].samples[s]['FT']
            for f in (x for x in to_merge[max_i].format if x not in ('GT',
                                                                     'FT')):
                output_record.samples[s][f] = to_merge[max_i].samples[s][f]
        return output_record

    def _compatible_genotypes(self, record1, record2):
        for s in record1.samples:
            gt1 = record1.samples[s].allele_indices
            gt2 = record2.samples[s].allele_indices
            # Canvas call might be ./1 while Manta call is 0/1
            # Sometimes Canvas and Manta disagree on hom/heterozygosity
            # just check for presence/absence of ALT allele(?)
            # TODO(?) IMPROVE THIS?
            if (1 in gt1) != (1 in gt2):
                return False
        return True

    def _get_contig_order(self):
        contig_order_canvas = get_contig_order(self.canvas)
        main_contig_order_canvas = dict((k, v) for k, v in
                                        contig_order_canvas.items() if
                                        chrom_re.match(k))
        contig_order_manta = get_contig_order(self.manta)
        main_contig_order_manta = dict((k, v) for k, v in
                                       contig_order_manta.items() if
                                       chrom_re.match(k))
        if main_contig_order_canvas != main_contig_order_manta:
            raise ValueError("Input VCFs do not share the same contig order!")
        for k in (x for x in contig_order_manta if x not in
                  main_contig_order_manta and x not in contig_order_canvas):
            contig_order_canvas[k] = len(contig_order_canvas)
        return contig_order_canvas

    def _combine_headers(self):
        for k, v in ((k, v) for k, v in self.canvas.header.info.items() if k
                     not in self.manta.header.info):
            self.manta.header.info.add(k,
                                       v.number,
                                       v.type,
                                       v.description)
        for k, v in ((k, v) for k, v in self.canvas.header.formats.items() if k
                     not in self.manta.header.formats):
            self.manta.header.formats.add(k,
                                          v.number,
                                          v.type,
                                          v.description)
        for header in (self.canvas.header, self.manta.header):
            header.info.add('MergeSource',
                            1,
                            'String',
                            'Caller origin of original CNV/SV call.')
        self.manta.header.info.add(
            'MergeAlts',
            '.',
            'String',
            'ALT allele(s) from merged call(s).')
        self.manta.header.info.add(
            'MergeStartDiff',
            '.',
            'Integer',
            'Difference between record start and start coordinate(s) of ' +
            'merged call(s).')
        self.manta.header.info.add(
            'MergeEndDiff',
            '.',
            'Integer',
            'Difference between record end and end coordinate(s) of merged ' +
            'call(s).')
        self.manta.header.info.add('CanvasSupport',
                                   0,
                                   'Flag',
                                   'Call has matching call from Canvas')
        self.manta.header.formats.add('CGT',
                                      1,
                                      'String',
                                      'GT field from matching Canvas call')
        self.manta.header.formats.add('CFT',
                                      1,
                                      'String',
                                      'FT field from matching Canvas call')
        return self.manta.header.copy()

    def _parse_annotation_filters(self, filters, vcf):
        annotation_filters = {'info': {'GAIN': None, 'LOSS': None},
                              'format': {'GAIN': None, 'LOSS': None}}
        if filters is None:
            return annotation_filters
        for field, filter_class in zip(('info', 'format'),
                                   (InfoFilter, FormatFilter)):
            filter_expressions = {'GAIN': [], 'LOSS': []}
            for exp in filters[field]:
                if len(exp) < 3 or len(exp) > 4:
                    raise ValueError("Could not parse filter expression {}"
                                     .format(" ".join(exp)))
                if len(exp) == 3:
                    for filter_list in filter_expressions.values():
                        filter_list.append(exp)
                else:
                    if exp[3].upper() not in annotation_filters[field]:
                        raise ValueError(
                            "CNV type '{}' in filter expression ".format(exp[3]) +
                            "'{}' not recognised. ".format(' '.join(exp)) +
                            "Valid CNV values are 'LOSS' or 'GAIN'.")
                    filter_expressions[exp[3].upper()].append(exp[:3])
            for cnv_type, expressions in filter_expressions.items():
                if expressions:
                    annotation_filters[field][cnv_type] = filter_class(
                        vcf=vcf, filters=expressions)
        return annotation_filters

    def _parse_filters(self, canvas_filters, manta_filters):
        '''
        Returns dict of dicts in format:
            genotyper->vcf_field->cnv_type->AnnotationFilter
        '''
        filters = dict()
        filters['canvas'] = self._parse_annotation_filters(canvas_filters,
                                                           self.canvas)
        filters['manta'] = self._parse_annotation_filters(manta_filters,
                                                          self.manta)
        return filters

    def _new_record_from_template(self, rec):
        '''Necessary due to potential weirdness with headers'''
        new_rec = self.header.new_record()
        for att in ['alts', 'chrom', 'id', 'qual', 'ref', 'stop']:
            setattr(new_rec, att, getattr(rec, att))
        new_rec.pos = max(1, new_rec.pos)  # pysam errors when pos is 0
        for k, v in rec.info.items():
            new_rec.info[k] = v
        for s in rec.samples:
            for k, v in rec.samples[s].items():
                new_rec.samples[s][k] = v
        return new_rec
