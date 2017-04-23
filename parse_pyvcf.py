import sys,argparse
import vcf
from pprint import pprint

parser = argparse.ArgumentParser(
                    description="parse input VCF",
                    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--input_vcf', help='Input VCF')

def convert_to_matrix(infile):
    """
        parse input VCF file into a 2D matrix
        of users x genotypes. Each entry will
        have value (0,1,2) for the number of 
        minor alleles. Additional info is stored
        in a separate dictionary.
    """
    genotypes = {}
    snp_info = {}

    vcf_reader = vcf.Reader(open(infile, 'r'))
    count = 0
    # for now, just take SNPs from chr 7
    for record in vcf_reader.fetch('7'):
        count += 1
        if count % 1000 == 0:
            print count
        snp_info[record.ID] = {"chrom":record.CHROM,
                               "pos":record.POS,
                               "ref":record.REF,
                               "alt":record.ALT}
        user_genos = {}
        for call in record.samples:
            if call.called:
                user_id = call.sample
                genos = call.data.GT
                geno_score = float(genos[0]) + float(genos[-1])
                user_genos[user_id] = geno_score
        genotypes[record.ID] = user_genos
    return genotypes, snp_info

def main(args):
    args = parser.parse_args(args)
    genotypes, snp_info = convert_to_matrix(args.input_vcf)
    pprint(snp_info)
    pprint(genotypes)

if __name__ == '__main__':
    main(sys.argv[1:])
