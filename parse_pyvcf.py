import sys,argparse

parser = argparse.ArgumentParser(
                    description="parse input VCF",
                    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--input_vcf', help='Input VCF')

def main(args):
    args = parser.parse_args(args)
    print args.input_vcf

if __name__ == '__main__':
    main(sys.argv[1:])
