from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Check and make sure the REF/ALT alleles in GWAS results corresponds to A1/A2 alleles in UKBB bim file

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2017/08/19
-------------------------------------------------------------------------
'''
import sys
import os
import argparse
import pandas as pd
import fileinput
import numpy as np

def snp_f(chr, snp_d = '/oak/stanford/groups/mrivas/private_data/ukbb/24983/snp/snp_download'):
    return os.path.join(snp_d, 'ukb_snp_chr{}_v2.bim'.format(chr))


def flipfix_A1A2(in_f, chr, is_qt, is_binary):
    if( int(is_qt) + int(is_binary) != 1):
        raise RuntimeError("is_qt xor is_binary should be 1")
        sys.exit(1)

    # read bim
    bim_df = pd.read_csv(
        snp_f(chr), sep='\t',
        names=['chr', 'rsid', 'dist', 'pos', 'A1', 'A2']
    )
    
    A1 = dict(zip(bim_df.pos.map(lambda x: str(x)), bim_df.A1))
    A2 = dict(zip(bim_df.pos.map(lambda x: str(x)), bim_df.A2))    

    if(is_qt):
        test_f = 5
    elif(is_binary):
        test_f = 6

    in_file = '/dev/stdin' if in_f is None else in_f

    for line in fileinput.input(in_file):
        if(line[0] == '#'):
            # print header line
            print(line.rstrip())
        else:
            l = line.split()
            POS=str(l[1])            
            REF=l[3]
            ALT=l[4]
            if(A1[POS] == REF and A2[POS] == ALT):
                # if there is no flip
                print(line.rstrip())
            elif(A1[POS] == ALT and A2[POS] == REF):
                # if there is a flip
                l[3] = ALT
                l[4] = REF
                if(is_qt and l[test_f] == 'ADD' and l[7] != 'NA'):
                    l[7] = -float(l[7])
                elif(is_binary and l[test_f] == 'ADD' and l[8] != 'NA'):
                    l[8] = np.exp(- np.log(float(l[8])))
                print('\t'.join([str(x) for x in l]))
            else:
                # multi-allelic, haven't implemented yet for this case. Dump error and exit
                raise RuntimeError('Unexpected entry: {}'.format(line.rstrip))
                sys.exit(1)

def main():  
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )    
    parser.add_argument('-i', metavar='i', default=None, help='input file [default: stdin]')
    parser.add_argument('-c', metavar='c', required=True,
                        help='chromosome')      
    parser.add_argument('-q', action='store_true', help='QT')
    parser.add_argument('-b', action='store_true', help='Binary')
    args = parser.parse_args()
    flipfix_A1A2(args.i, args.c, args.q, args.b)
   
if __name__ == "__main__":
    main()

