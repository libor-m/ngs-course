"""Given VCF file and population assignments, calculate per site dxy
when a sample is not in pop file, do not consider it at all

http://binhe.org/2011/12/29/calculate-nucleotide-diversity-per-base-pair-summation-method/
sum ( 2j(n-j) / n(n-1) )

# pi calc and general thoughts here
http://seqanswers.com/forums/showthread.php?t=33639

Also, I don't know if anyone here is analyzing multiple populations, 
or would like to calculate Dxy, but I think this formula works using 
the variables stored for the pi calculations and will give you per-site 
Dxy between two populations.

((n1*n2) - sum_i(J1_i*J2_i)) / (n1*n2)

for population 1 and population 2 where J_i is the allele count (# of 0,1,2, or 3's)
at a given site in the vcf for either population.

Code:
   if n_seq_POP_1 > 0 and n_seq_POP_2 >0 :
       POP_1_POP_2_dxy = ((n_seq_POP_1 * n_seq_POP_2) - \
           ((j0_POP_1 * j0_POP_2) + \
            (j1_POP_1 * j1_POP_2) + \
            (j2_POP_1 * j2_POP_2) + \
            (j3_POP_1 * j3_POP_2)) ) / \
        (n_seq_POP_1 * n_seq_POP_2)

Charlesworth 1998: 
pi_b = sum_for_i<j(w_i * w_j * pi_ij) / sum_for_i<j(w_i * w_j)

w_x represents weight of population size, so either use samples sizes or weight equally
that is w_x = 1/n where n is number of populations, sum(w_i) = 1

"""
import sys
import itertools

import vcf

def report(var, values=[], count=3, sep="\t"):
    """report calculated value in some reasonable format
    report max `count` fields
    fill all empty fields with blanks
    """
    
    fields = [var.CHROM, var.POS]
    fields.extend(values[:count])
    blanks = [''] * (count - len(values))
    fields.extend(blanks)
    print sep.join(map(str, fields))

def harvest_alleles(samples, pops):
    """harvest are alleles for each population into {pop : [alleles]}
    """
    pop_alleles = {}

    for s in samples:
        if s.sample not in pops:
            continue

        alleles = pop_alleles.setdefault(pops[s.sample], [])
        alleles.extend(s.gt_alleles)

    return pop_alleles

def calc_dxy(pop_alleles):
    """calculate dxy
    by the dumbest approach - count differences for all comparisons
    hope it won't be slow/expensive
    """
    combs = itertools.product(*pop_alleles.values())
    
    diffs, total = (0, 0)
    for a, b in combs:
        if a != b:
            diffs += 1
        total += 1

    return (diffs, total)

def main():
    if len(sys.argv) < 3:
        sys.exit("use: %s vcf popfile (popfile is sample<whitespace>pop)")

    variants = vcf.Reader(filename=sys.argv[1])

    # read in the sample, population table
    # into a dictionary
    with open(sys.argv[2]) as f:
        pops = {s:p for s, p in [l.strip().split() for l in f.readlines()]}

    # for each variant
    for var in variants:
        # pick samples that have depth > 0
        # ie do not use imputed genotypes for the calculation
        samples = [s for s in var.samples if s['DP'] > 0]

        # take every haploid sequence as an independent sample (TODO: is this ok?)
        # assign all calls to the population
        pop_alleles = harvest_alleles(samples, pops)

        if len(pop_alleles) != 2:
            report(var)
            continue

        # get and report the dxy value
        ndiff, ntotal = calc_dxy(pop_alleles)
        report(var, [ndiff, ntotal, float(ndiff)/ntotal])

if __name__ == '__main__':
    main()