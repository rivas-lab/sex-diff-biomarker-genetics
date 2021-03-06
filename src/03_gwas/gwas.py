#!/bin/python
import os
from random import randint
_README_='''
A script for running GWAS with UK Biobank data (array/exome/imputed/cnv genotypes) using PLINK.

Author: Matthew Aguirre (SUNET: magu)

(Updated 10/14/2019 by E Flynn for sex-div analysis)

'''

def filterPopFile(outDir, popFile='', keepSexFile=''):
            # reformat the popFile to only contain the iids from one sex
        if not popFile or popFile=='':
            popFileSexFiltered=keepSexFile
        else:
            pop_iids = []
            keep_iids = []
            with open(popFile, 'r') as pop_f:
                for line in pop_f:
                    iid1, iid2 = line.split()
                    pop_iids.append(iid1.strip())
            with open(keepSexFile, 'r') as sex_f:
                for line in sex_f:
                    iid1 = line.strip()
                    keep_iids.append(iid1)
            filt_iids = set(pop_iids).intersection(set(keep_iids))
            # TODO: update to use temporary file (but that remains open thru the plink run?)
            popFileSexFiltered='{0}/tmp_pop_sex_{1}.keep'.format(outDir, randint(10000,99999))
            with open(popFileSexFiltered, 'w') as keep_f:
                for filt_iid in filt_iids:
                    keep_f.write("{0}\t{1}\n".format(filt_iid, filt_iid))
        return popFileSexFiltered

def make_plink_command(bpFile, pheFile, outFile, outDir, pop, cores=None, memory=None, related=False, plink1=False, 
                       arrayCovar=False, sexDiv=False, keepSex='', keepSexFile='', includeX=False):
    # paths to plink genotypes, input phenotypes, output directory are passed
    qcDir         = '/ukbb24983/sqc/'
    popFile       = os.path.join(qcDir,'population_stratification','ukb24983_{}.phe'.format(pop)) if pop != 'all' else ''
    unrelatedFile = os.path.join(qcDir,'ukb24983_v2.not_used_in_pca.phe') if not related else '' 
    arrayVarFile  = os.path.join(qcDir,'{}_array_variants.txt'.format('both' if arrayCovar else 'one')) if '/cal/' in bpFile else ''
    is_cnv_burden = os.path.basename(bpFile) == 'burden'
    if os.path.dirname(pheFile).split('/')[-1] == "binary": 
        is_biomarker_binary = True
    elif len(os.path.basename(pheFile).split('_')) < 2:
        is_biomarker_binary = False
    else:
        is_biomarker_binary = os.path.basename(pheFile).split('_')[1]  == 'binary'
    if is_cnv_burden:
        covarFile = 'ukbb24983/cnv/pgen/ukb24983_cnv_burden.covar'
    elif is_biomarker_binary:
        covarFile = 'biomarkers/'
    else:
        covarFile = os.path.join(qcDir, 'ukb24983_GWAS_covar.phe')
    
    # filter the population file to only contain IIDs of the specified sex
    if sexDiv:
        popFileSexFiltered=filterPopFile(outDir, popFile, keepSexFile)
                                                                                                                
    if not includeX:
        xChrStr=""
    else:
        xChrStr=",X,XY,Y"
    # paste together the command from constituent parts
    return " ".join(["plink" if plink1 else "plink2", 
                     "--threads {0}".format(cores) if cores is not None else "",
                     "--memory {0}".format(memory) if memory is not None else "",                     
                     "--bfile" if plink1 or '/imp/' in bpFile else "--bpfile", bpFile, "--chr 1-22{0}".format(xChrStr),
                     "--maf 0.005" if "/imp" in bpFile else "",
                     "--pheno", pheFile, "--pheno-quantile-normalize",
                     "--glm firth-fallback hide-covar omit-ref",
                     "--keep {0}".format(popFile) if (popFile and not sexDiv) else '', "--keep {0}".format(popFileSexFiltered) if sexDiv else "", 
                     "--remove {0}".format(unrelatedFile) if unrelatedFile else "",
                     "--extract {0}".format(arrayVarFile) if arrayVarFile else "",
                     "--covar", covarFile, 
                     "--covar-name age ", "sex " if not sexDiv else "", 
                     "Array " if arrayCovar else "", "PC1-PC4", "N_CNV LEN_CNV --covar-variance-standardize" if is_cnv_burden else "PC5-PC10 FastingTime --covar-variance-standardize" if is_biomarker_binary else "",
                     "--covar-variance-standardize --vif 100000000" if pop in ['non_british_white', 'african', 'e_asian', 's_asian'] else "",
                     #"--xchr-model 1" if includeX else "",
                     "--out", outFile]) 


def make_batch_file(batchFile, plinkCmd, cores, memory, time, partitions):
    with open(batchFile, 'w') as f:
       # formats options for sbatch header and pastes input command below it
       f.write("\n".join(["#!/bin/bash","",
                          "#SBATCH --job-name=RL_GWAS",
                          "#SBATCH --output={}".format(os.path.join(os.path.dirname(batchFile), "rl-gwas.%A-%a.out")),
                          "#SBATCH --cores={}".format(cores),
                          "#SBATCH --mem={}".format(memory),
                          "#SBATCH --time={}".format(time),
                          "#SBATCH -p {}".format(','.join(partitions)),
                          '#SBATCH --constraint="CPU_GEN:HSW|CPU_GEN:BDW|CPU_GEN:SKX"', # plink2 avx2 compatibility
                          "", plinkCmd]))
    return batchFile


def run_gwas(kind, pheFile, outDir='', pop='white_british', related=False, plink1=False, 
             logDir='', cores="4", memory="24000", time="1-00:00:00", partition=["normal","owners"], now=False,
             sexDiv=False, keepSex='', keepSexFile='', includeX=False):
    # ensure usage
    if not os.path.isfile(pheFile):
        raise ValueError("Error: phenotype file {0} does not exist!".format(pheFile))
    if not os.path.isdir(outDir):
        raise ValueError("output directory {0} does not exist!".format(outDir))
    if not os.path.isdir(logDir):
        print("Warning: Logging directory either does not exist or was unspecified. Defaulting to {0}...".format(os.getcwd()))
        logDir=os.getcwd()
    if sexDiv and not os.path.isfile(keepSexFile):
        raise ValueError("Error: keep sex file {0} does not exist!".format(keepSexFile))
    if sexDiv and kind != "genotyped":
        print("Warning: sex div GWAS has only been tested on array data! Please make sure to test this more or use --run-array.")
    if pop not in ['all', 'white_british', 'non_british_white', 'african', 's_asian', 'e_asian']:
        raise ValueError("population must be one of (all, white_british, non_british_white, african, s_asian, e_asian)")
    # paths for running gwas
    pgen_root='ukbb/24983/'
    imp_bfile_path=os.path.join(pgen_root,'imp','pgen','ukb24983_imp_chr${SLURM_ARRAY_TASK_ID}_v3')
    cal_bfile_path=os.path.join(pgen_root,'cal','pgen','ukb24983_cal_cALL_v2_hg19')
    exome_spb_path=os.path.join(pgen_root,'exome','pgen','spb','data','ukb_exm_spb')
    exome_fe_path=os.path.join(pgen_root,'exome','pgen','fe','data','ukb_exm_fe')
    cnv_bfile_path=os.path.join(pgen_root,'cnv','pgen','cnv') + ' --mac 15'
    cnv_burden_path=os.path.join(pgen_root,'cnv','pgen','burden')
    hla_bfile_path=os.path.join(pgen_root,'hla','pgen','ukb_hla_v3')
    pheParts=os.path.basename(pheFile).split('.')
    pheName=".".join(pheParts[0:(len(pheParts)-1)])
    outFile=os.path.join(outDir, 'ukb24983_v2_{0}.{1}{2}.{3}'.format('hg38' if 'exome' in kind else 'hg19', pheName, "_{0}".format(keepSex) if sexDiv else "", kind))
    # this is where the fun happens
    if kind == 'imputed':
        # needs one array job per chromosome (for compute time), hence the slurm variable
        outFile += '.chr${SLURM_ARRAY_TASK_ID}'
        cmd = make_plink_command(bpFile  = imp_bfile_path,
                                 pheFile = pheFile,
                                 outFile = outFile,
                                 outDir  = outDir,
                                 pop     = pop,
                                 related = related,
                                 plink1  = plink1,
                                 cores   = cores,
                                 memory  = memory,
                                 arrayCovar = True,
                                 sexDiv = sexDiv,
                                 keepSex = keepSex,
                                 keepSexFile = keepSexFile,
                                 includeX = includeX)
    elif kind == 'genotyped': 
        # needs one plink call with genotyping array a covariate, and one without
        outFile1 = outFile+'.both_arrays'
        cmd1 = make_plink_command(bpFile  = cal_bfile_path,
                                  pheFile = pheFile,
                                  outFile = outFile1,
                                  outDir = outDir,
                                  pop     = pop,
                                  related = related,
                                  plink1  = plink1,
                                  cores   = cores,
                                  memory  = memory,
                                  arrayCovar = True,
                                  sexDiv = sexDiv,
                                  keepSex = keepSex,
                                  keepSexFile = keepSexFile, 
                                  includeX = includeX)
        outFile2 = outFile+'.one_array'
        cmd2 = make_plink_command(bpFile  = cal_bfile_path,
                                  pheFile = pheFile,
                                  outFile = outFile2,
                                  outDir  = outDir,
                                  pop     = pop,
                                  related = related,
                                  plink1  = plink1,
                                  cores   = cores,
                                  memory  = memory,
                                  arrayCovar = False,
                                  sexDiv = sexDiv,
                                  keepSex = keepSex,
                                  keepSexFile = keepSexFile, 
                                  includeX = includeX)
        # join the plink calls, add some bash at the bottom to combine the output
        cmd = "\n\n".join([cmd1, cmd2] +  # this is the plink part, below joins the two files
                          ["if [ -f {0}.*.{3} ]; then cat {0}.*.{3} {1}.*.{3} | sort -k1,1n -k2,2n -u > {2}.{3}; fi".format(
                               outFile1, outFile2, outFile, suffix) for suffix in ['glm.linear', 'glm.logistic.hybrid']] + 
                          ["cat {0}.log {1}.log > {2}.log".format(outFile1, outFile2, outFile),
                           "rm {0}.* {1}.*".format(outFile1, outFile2)])
    elif kind == 'cnv' or kind == 'cnv-burden':
        cnv_path = cnv_bfile_path if kind == 'cnv' else cnv_burden_path
        cmd = make_plink_command(bpFile  = cnv_path,
                                 pheFile = pheFile,
                                 outFile = outFile,
                                 outDir  = outDir,
                                 pop     = pop,
                                 related = related,
                                 plink1  = plink1,
                                 cores   = cores,
                                 memory  = memory,
                                 arrayCovar = False,
                                 sexDiv = sexDiv,
                                 keepSex = keepSex,
                                 keepSexFile = keepSexFile, 
                                 includeX = includeX) 
    # more usage management, in case someone wants to import the function for use elsewhere
    elif kind == 'exome-spb' or kind == 'exome-fe':
        exome_bfile_path = exome_spb_path if kind == 'exome-spb' else exome_fe_path
        cmd = make_plink_command(bpFile  = exome_bfile_path,
                                 pheFile = pheFile,
                                 outFile = outFile,
                                 outDir  = outDir,
                                 pop     = pop,
                                 related = related,
                                 plink1  = plink1,
                                 cores   = cores,
                                 memory  = memory,
                                 arrayCovar = False,
                                 sexDiv = sexDiv,
                                 keepSex = keepSex,
                                 keepSexFile = keepSexFile, 
                                 includeX = includeX)
    elif kind == 'hla':
        cmd = make_plink_command(bpFile  = hla_bfile_path,
                                 pheFile = pheFile,
                                 outFile = outFile,
                                 outDir  = outDir,
                                 pop     = pop,
                                 related = related,
                                 plink1  = plink1,
                                 cores   = cores,
                                 memory  = memory,
                                 arrayCovar = False,
                                 sexDiv = sexDiv,
                                 keepSex = keepSex,
                                 keepSexFile = keepSexFile, 
                                 includeX = includeX)
    else:
        raise ValueError("argument kind must be one of (imputed, genotyped, exome-spb, exome-fe): {0} was provided".format(kind))
    # run immediately, OR make the batch job submission file, then call it with an appropriate array handler
    if now:
        print("Running the below: \n'''\n" + cmd + "\n'''\n") 
        os.system(cmd)
    else:
        sbatch = make_batch_file(batchFile = os.path.join(logDir, "gwas.{0}.{1}.sbatch.sh".format(kind,pheName)),
                                 plinkCmd  = cmd,
                                 cores     = cores,
                                 memory    = memory,
                                 time      = time,
                                 partitions = partition)
        os.system(" ".join(("sbatch", "--array=1{}".format("-22" if kind == 'imputed' else ""), sbatch))) 
    return


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    parser.add_argument('--plink1', dest='plink1', action='store_true',
                            help='Flag to run GWAS with PLINK v1.9 instead of v2.0')
    parser.add_argument('--run-array', dest="arr", action='store_true',
                            help='Run GWAS on directly genotyped (array) data')
    parser.add_argument('--run-exome', dest="ex1", action='store_true',
                            help="Run GWAS on exome data (Regeneron calls)")
    parser.add_argument('--run-exome-gatk', dest="ex2", action='store_true',
                            help="Run GWAS on exome data (GATK calls)")
    parser.add_argument('--run-imputed', dest="imp", action='store_true',
                            help='Run GWAS on imputed data') 
    parser.add_argument('--run-cnv', dest="cnva", action='store_true',
                            help='Run GWAS on array-derived CNV genotypes') 
    parser.add_argument('--run-cnv-burden', dest="cnvb", action='store_true',
                            help='Run CNV burden test (GWAS on 0/1 CNV overlaps gene, from array-derived CNV genotypes)') 
    parser.add_argument('--run-hla', dest="hla", action='store_true',
                            help='Run GWAS on imputed HLA allelotypes') 
    parser.add_argument('--pheno', dest="pheno", required=True, nargs='*',
                            help='Path to phenotype file(s)')
    parser.add_argument('--out', dest="outDir", required=True, nargs=1,
                            help='Path to desired output *directory*. Summary stats will be output according to phenotype name (derived from passed file) and Rivas Lab specification for GBE. Defaults to current working directory.')
    parser.add_argument('--population', dest="pop", required=False, default=["white_british"], nargs=1,
                            help='Flag to indicate which ethnic group to use for GWAS. Must be one of all, white_british, non_british_white, e_asian, s_asian, african')
    parser.add_argument('--keep-related', dest="relatives", action='store_true',
                            help='Flag to keep related individuals in GWAS. Default is to remove them.')
    parser.add_argument('--cores', dest="cores", required=False, default=[None], nargs=1,
                            help='For underlying plink command/batch job submission: Amount of cores to request. Default is 4 cores for batch jobs.')    
    parser.add_argument('--memory', dest="mem", required=False, default=[None], nargs=1,
                            help='For underlying plink command/batch job submission: Amount of memory (in MB) to request. Default is 24000 for batch jobs.')
    parser.add_argument('--batch-time', dest="sb_time", required=False, default=["24:00:00"], nargs=1,
                            help='For underlying batch job submission: Amount of time (DD-HH:MM:SS) to request. Default is 24 hours.')
    parser.add_argument('--batch-partitions', dest="sb_parti", required=False, default=["normal","owners"], nargs='*',
                            help='For underlying batch job submission: Compute partition to submit jobs. Default is normal,owners.')
    parser.add_argument('--log-dir', dest="log", required=False, nargs=1, default=[''],
                            help="Directory in which to place log files from this script and its related SLURM jobs. Default is current working directory.")
    parser.add_argument('--run-now', dest="local", action='store_true',
                            help='Flag to run GWAS immediately, on the local machine, rather than submitting a script with sbatch. Not available for use with --run-imputed.')
    parser.add_argument('--sex-div', dest="sex_div",required=False, default=False,
                            help='Whether to run sex-divided GWAS, this removes sex as a covariate from the analysis. Default is False. If specifying this, please also use --keep-sex to specify the sex to keep and --keep-sex-file to specify the file location.')
    parser.add_argument('--keep-sex', dest="keep_sex", required=False, default='',
                            help='Which sex to keep, will be used in constructing the output filenames, for use with --sex-div.')
    parser.add_argument('--keep-sex-file', dest="keep_sex_file", required=False, default='',
                            help='Location of the file specifying the IIDs to include related to that sex, for use with --sex-div.')
    parser.add_argument('--include-x', dest="include_x", required=False, default=False,
                            help='Whether to include the X chromosome, defaults to False')

    args = parser.parse_args()
    # TODO: feature add: genotype model  
    print(args) 
    # ensure handler-relevant usage (more insurance is in run_gwas()):
    flags = [args.imp, args.arr, args.ex1, args.ex2, args.cnva, args.cnvb, args.hla]
    kinds = ['imputed','genotyped','exome-spb','exome-fe','cnv','cnv-burden', 'hla']
    if not any(flags):
        raise ValueError("Error: no analysis specified, did you mean to add --run-array?")
    if args.local and args.imp:
        raise ValueError("--run-imputed cannot be present in conjunction with --run-now!")
    if (args.sex_div and args.keep_sex == '') or (args.sex_div and args.keep_sex_file == ''):
        raise ValueError("Sex-div analysis is indicated but either the sex to keep or the file for that is missing. Please use --keep-sex and --keep-sex-file to specify both.")
    # add default parameters to batch job submission
    if not args.local:
        if args.cores[0] is None:
            args.cores[0] = "4"
        if args.mem[0] is None:
            args.mem[0] = "24000"
    # lol i hope this works
    for flag,kind in filter(lambda x:x[0], zip(flags,kinds)):
        run_gwas(kind=kind, pheFile=args.pheno[0], outDir=args.outDir[0], 
                 pop=args.pop[0], related=args.relatives, plink1=args.plink1, 
                 logDir=args.log[0], cores=args.cores[0], memory=args.mem[0], 
                 time=args.sb_time[0], partition=args.sb_parti, now=args.local, 
                 sexDiv=args.sex_div, keepSex=args.keep_sex, keepSexFile=args.keep_sex_file,
                 includeX=args.include_x)
