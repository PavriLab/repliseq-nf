#!/usr/bin/env python

#######################################################################
#######################################################################
## Copied on March 25, 2020 from nf-core/atacseq
#######################################################################
#######################################################################

import os
import sys
import requests
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Reformat repliseq-nf design file and check its contents.'
Epilog = """Example usage: python check_design.py <DESIGN_FILE_IN> <DESIGN_FILE_OUT>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('DESIGN_FILE_IN', help="Input design file.")
argParser.add_argument('DESIGN_FILE_OUT', help="Output design file.")
args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def check_design(DesignFileIn,DesignFileOut):

    ERROR_STR = 'ERROR: Please check design file'
    HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2']

    ## CHECK HEADER
    fin = open(DesignFileIn,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print "{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER))
        sys.exit(1)

    numColList = []
    groupRepDict = {}
    while True:
        line = fin.readline()
        if line:
            lspl = [x.strip() for x in line.strip().split(',') if x]

            ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
            numCols = len(lspl)
            if numCols not in [3,4]:
                print "{}: Invalid number of columns (3 for single-end or 4 for paired-end)!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)
            numColList.append(numCols)

            ## CHECK GROUP COLUMN HAS NO SPACES
            group,replicate,fastQFiles = lspl[0],lspl[1],lspl[2:]
            if group.find(' ') != -1:
                print "{}: Group id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)

            ## CHECK REPLICATE COLUMN IS INTEGER
            if not replicate.isdigit():
                print "{}: Replicate id not an integer!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)

            for fastq in fastQFiles:
                ## CHECK FASTQ FILE EXTENSION
                if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                    print "{}: FastQ file has incorrect extension (has to be '.fastq.gz' or 'fq.gz') - {}\nLine: '{}'".format(ERROR_STR,fastq,line.strip())
                    sys.exit(1)

            ## CREATE GROUP MAPPING DICT = {GROUP_ID: {REPLICATE_ID:[[FASTQ_FILES]]}
            replicate = int(replicate)
            if not groupRepDict.has_key(group):
                groupRepDict[group] = {}
            if not groupRepDict[group].has_key(replicate):
                groupRepDict[group][replicate] = []
            groupRepDict[group][replicate].append(fastQFiles)

        else:
            fin.close()
            break

    ## CHECK IF DATA IS PAIRED-END OR SINGLE-END AND NOT A MIXTURE
    if min(numColList) != max(numColList):
        print "{}: Mixture of paired-end and single-end reads!".format(ERROR_STR)
        sys.exit(1)

    ## CHECK IF MULTIPLE GROUPS EXIST
    multiGroups = False
    if len(groupRepDict) > 1:
        multiGroups = True

    ## WRITE TO FILE
    numRepList = []
    fout = open(DesignFileOut,'w')
    fout.write(','.join(['sample_id','fastq_1','fastq_2']) + '\n')
    for group in sorted(groupRepDict.keys()):

        ## CHECK THAT REPLICATE IDS ARE IN FORMAT 1..<NUM_REPLICATES>
        uniq_rep_ids = set(groupRepDict[group].keys())
        if len(uniq_rep_ids) != max(uniq_rep_ids):
            print "{}: Replicate IDs must start with 1..<num_replicates>\nGroup: {}, Replicate IDs: {}".format(ERROR_STR,group,list(uniq_rep_ids))
            sys.exit(1)
        numRepList.append(max(uniq_rep_ids))

        for replicate in sorted(groupRepDict[group].keys()):
            for idx in range(len(groupRepDict[group][replicate])):
                fastQFiles = groupRepDict[group][replicate][idx]
                sample_id = "{}_R{}_T{}".format(group,replicate,idx+1)
                if len(fastQFiles) == 1:
                    fout.write(','.join([sample_id] + fastQFiles) + ',\n')
                else:
                    fout.write(','.join([sample_id] + fastQFiles) + '\n')
    fout.close()

    ## CHECK IF REPLICATES IN DESIGN
    repsExist = False
    if max(numRepList) != 1:
        repsExist = True

    ## CHECK FOR BALANCED DESIGN ACROSS MULTIPLE GROUPS.
    balancedDesign = False
    if len(set(numRepList)) == 1 and multiGroups and repsExist:
        balancedDesign = True

############################################
############################################
## RUN FUNCTION
############################################
############################################

check_design(DesignFileIn=args.DESIGN_FILE_IN,DesignFileOut=args.DESIGN_FILE_OUT)

############################################
############################################
############################################
############################################
