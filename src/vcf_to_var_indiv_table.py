#!/usr/bin/env python

import vcf
import sys
import pysam
import normalize # https://github.com/ericminikel/minimal_representation/blob/master/normalize.py

# suggested use: vcf_to_var_indiv_table.py $lv500genos $b37ref > genotypes.table
# or, bsub -q priority -J lv500 -o lv500.o -e lv500.e ". private_paths.bash; vcf_to_var_indiv_table.py $lv500genos $b37ref > genotypes.table"

def vcf_to_var_indiv_table(genos_vcf, fasta_path):
    genos = vcf.Reader(filename=genos_vcf,strict_whitespace=True)
    pysam_fasta = pysam.FastaFile(fasta_path)
    sys.stdout.write('\t'.join(map(str,['chrom','pos','ref','alt','filter','sample','zyg','gt','gq','dp','ad']))+'\n')
    for record in genos:
        if len(record.FILTER) > 0:
            filter_val = record.FILTER[0]
        else:
            filter_val = 'PASS'
        for alt_allele in record.ALT:
            allele_number = record.ALT.index(alt_allele) + 1
            alt_allele = str(alt_allele)
            chrom, pos, ref, alt = normalize.normalize(pysam_fasta, record.CHROM, record.POS, record.REF, alt_allele)
            for sample in record.samples:
                if sample['GT'] is None:
                    continue
                if allele_number in map(int,sample.gt_alleles):
                    if all(x==allele_number for x in map(int,sample.gt_alleles)):
                        zyg = 'hom'
                    else:
                        zyg = 'het'
                    sys.stdout.write('\t'.join(map(str,[chrom, pos, ref, alt, filter_val,sample.sample,zyg,sample['GT'],sample['GQ'],sample['DP'],sample['AD'][allele_number-1]]))+'\n')

if __name__ == '__main__':
    vcf_to_var_indiv_table(sys.argv[1], sys.argv[2])

