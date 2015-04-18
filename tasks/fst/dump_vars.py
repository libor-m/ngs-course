"""using vcf file and a gff annotation dump the data 

format that should be enough for most of the analyses in R
- standard vcf columns: CHROM, POS, QUAL
- is indel flag (INDEL key in info)
- raw read depth (DP), DP4
- position in contig
- position in zebra finch
-- exact mapped with NAs for varints w/o mapping
-- binned by min-max exon
- contig name
"""

import sys
import itertools

import vcf
import pybedtools
import pysam

# workaround for the half broken pybedtools
def pybed_get_interval(pblist, interval):
    """get features in interval for sorted list of intervals

    useful when pybedtools do not work at all;)
    """
    for i in pblist:
        # skip precedin data
        if i.chrom != interval.chrom:
            continue
        if i.start < interval.start:
            continue
        # yield matching intervals
        if i.start >= interval.start and i.end <= interval.end:
            yield i
        # early exit
        if i.start > interval.end:
            return

def pybed_find_features(ilist, chrom, pos):
    """find features from list that contain pos
    """
    for i in ilist:
        if i.chrom != chrom:
            continue
        if i.end < pos:
            continue
        if i.start <= pos and i.end > pos:
            yield i
        if i.start > pos:
            return

def fetch_intervals(ti, interval):
    """fetch pybedtools.Intervals from tabix indexed file
    ti is pysam.Tabixfile
    interval is pybedtools.Interval
    """
    mx = [l.split("\t") for l in ti.fetch(interval.chrom, interval.start, interval.stop)]
    return [pybedtools.cbedtools.create_interval_from_list(fs) for fs in mx]

def print_header():
    print "\t".join([
        'chrom',
        'pos',
        'qual',
        'filter',
        'type',
        'subtype',
        'raw_depth',
        'ref_f',
        'ref_r',
        'alt_f',
        'alt_r',
        'contig_name',
        'mvz_name',
        'contig_size',
        'zf_min',
        'zf_max',
        'zf_pos',
        'mvz_seqpart',
        ])

def main():

    if len(sys.argv) < 3:
        sys.exit("use: %s gff variants" % sys.argv[0])

    # load all the annotations, because pybedtools is broken..
    annotations = pybedtools.BedTool(sys.argv[1])
    annotations._isbam = False
    an_tabix = pysam.Tabixfile(sys.argv[1])

    variants = vcf.Reader(filename=sys.argv[2])

    print_header()

    # for each mrna pick all exons with Target attribute
    # and all variants in the mrna
    # for each variant, if in exon, output with exact coords,
    # otherwise output with NA and mrna block

    def has_target(i):
        return 'Target' in i.attrs

    def exon_target(i):
        target = i.attrs['Target']
        chrom, start, end, _ = target.split(" ", 4)
        return pybedtools.Interval(chrom, int(start), int(end))

    def is_type(i, itype):
        return i.fields[2] == itype

    def is_source(i, src):
        return i.fields[1] == src

    for mrna in itertools.ifilter(lambda x: is_type(x, 'mRNA'), annotations):
        feats = fetch_intervals(an_tabix, mrna)
        exons0 = filter(has_target, filter(lambda x: is_type(x, 'exon'), feats))

        # pick only the exons mapping to the putative chromosome
        exons = filter(lambda x: exon_target(x).chrom == mrna.chrom, exons0)

        # no exons means no chromosome information, we don't need such variants
        if not exons:
            continue

        # use features from mvz pipeline
        mvz = filter(lambda x: is_source(x, 'mvz-annot'), feats)
        
        # pick the exon range
        exon_min = min(exon_target(x).start for x in exons)
        exon_max = max(exon_target(x).end for x in exons)

        mvars = list(variants.fetch(mrna.chrom, mrna.start, mrna.end))

        for var in mvars:
            f = list(pybed_find_features(exons, var.CHROM, var.POS))
            
            # if found a mapped exon, translate the coordinates
            # relative to the exon start
            if f:
                # pick the lowest exon
                minex = f[0]
                refpos = "%d" % (exon_target(minex).start + var.POS - minex.start)
            else:
                refpos = "NA"

            # list of feature types from mvz annotation pipeline (CDS, 3utr, 5utr..)
            mvzf = list(pybed_find_features(mvz, var.CHROM, var.POS))
            if mvzf:
                mvz_type = ",".join(i.fields[2] for i in mvzf)
            else:
                mvz_type = "NA"

            # try to pick mvz CDS feature and extract the Name(s) from it
            mvz_cds = filter(lambda x: is_type(x, 'CDS'), mvzf)
            if mvz_cds:
                try:
                    mvz_name = ",".join(i.attrs['Name'] for i in mvz_cds if 'Name' in i.attrs)
                except:
                    print >> sys.stderr, "error at var:", var.CHROM, var.POS
            else:
                mvz_name = "NA"

            print "\t".join(map(str, [
                var.CHROM, 
                var.POS,
                var.QUAL,
                var.FILTER,
                var.var_type,
                var.var_subtype,
                var.INFO['DP'],
                "\t".join(str(x) for x in var.INFO['DP4']),
                mrna.name,
                mvz_name,
                mrna.end - mrna.start,
                exon_min,
                exon_max,
                refpos,
                mvz_type,
                ]))

if __name__ == '__main__':
    main()
