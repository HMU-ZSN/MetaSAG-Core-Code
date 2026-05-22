import time
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#%matplotlib inline
import seaborn as sns
import os
import shutil
import sys
import warnings
import subprocess

#Time decorator
def timeit(func):
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = func(*args,**kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
        return result
    wrapper.__doc__=func.__doc__
    return wrapper


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

####### W1 stent sequence recognition #######

# returns the hamming distance fo two strings
def hamdist(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

# returns the closest hamming distance of the starting location of 2 sequences, if larger than 5, return -1
def find_W1_loc(read, W1 = "GAGTGATTGCTTGTGACGCCTT", threshold = 5):
    length_W1 = len(W1)
    start_loc = 8
    tests = 4
    closest_loc = start_loc
    diff = length_W1
    for i in range(tests):
        diff_temp = hamdist(read[start_loc + i:start_loc + i + length_W1], W1)
        if diff > diff_temp:
            closest_loc = i + start_loc
            diff = diff_temp
        #if hamming diatance less then threshold, then the location is correct and returned

        if diff < threshold:
            return closest_loc
    return 0

def closest_index(str1, ref_list, mismatch_allow = 1):
    closest_index_list = []
    closest_ham_dist = len(str1)
    for i in range(len(ref_list)):
        ham_dist_temp = hamdist(str1, ref_list[i])
        if ham_dist_temp < closest_ham_dist:# if a closer element if found, initialize closest_index_list and closest_ham_dist
            closest_ham_dist = ham_dist_temp
            closest_index_list = [ref_list[i]]
        elif ham_dist_temp == closest_ham_dist:# if it's the same length, add it to length list
            closest_index_list += [ref_list[i]]
    if closest_ham_dist <= mismatch_allow:
        return closest_index_list
    else:
        return []

#First, write a function for reading to recognize barcodes correctly
fastq = 'GAGTTACTCCGAGTGATTGCTTGTGACGCCTTCGGCACATTCGTCGGCAGCGTCAGATCTGTATAAGAGACAGTCTCATAGGCTTAACGCCGACGTGATTGGCATTTTTCGCGAGGAGACCATTGCTGGCATAACGGTGCTTATGTTGCGTGAAGGCC'

PYTHONDIR = os.path.dirname(os.path.abspath(__file__))
barcode_temp = pd.read_excel(os.path.join(PYTHONDIR, 'Barcode.xlsx'))

bc1_temp = barcode_temp['bc1'].values
bc2_temp = barcode_temp['bc2'].values

bc1_number = 96
bc2_number = 384

bc2 = {}
for i in range(bc2_number):
    # bc2 += [str(bc2_temp[i][31:39])]
    bc2[str(bc2_temp[i][31:39])] = i

bc2_list = list(bc2.keys())

bc1 = {}
for i in range(bc1_number):
    loc = bc1_temp[i].find("AGATCGGAAGAGCGTCGTGTAGGGAAAGAG")
    # bc1 += [str(bc1_temp[i][23:loc - 1])]
    bc1[str(bc1_temp[i][23:loc - 1])] = i
bc1_list = list(bc1.keys())

def FindBarcode(read):

    warning1 = 'Len'
    warning2 = 'W1'
    warning3 = 'BC'


    #The first step is to determine if the length is greater than 80 (otherwise, return waring1)
    if len(read) <= 80:
        return warning1
    else:
        #The second step is to determine whether W1 is included (otherwise, return warning2)
        loc = find_W1_loc(read)
        if loc == 0:
            return warning2
        else:
            #Step 3: Whether barcode1 and barcode2 can be found (if either of them cannot be found, return warning3)
            bc1_temp = reverseComplement(read[:loc])
            bc2_temp = reverseComplement(read[loc + 22:loc + 30])

            if bc1_temp in bc1_list:
                read_bc1 = bc1[bc1_temp]
            else:
                #Use closest index fuzzy matching
                bc1_cl = closest_index(bc1_temp,bc1_list)
                if len(bc1_cl) == 1:
                    read_bc1 = bc1[bc1_cl[0]]
                else:
                    return warning3

            if bc2_temp in bc2_list:
                read_bc2 = bc2[bc2_temp]
            else:
                # Use closest index fuzzy matching
                bc2_cl = closest_index(bc2_temp, bc2_list)
                if len(bc2_cl) == 1:
                    read_bc2 = bc2[bc2_cl[0]]
                else:
                    return warning3

    index = read_bc1 * bc2_number + read_bc2
    return index


#If the sequence length is qualified and both Barcoded + W1 are recognized, return the Barcode code of this read
#If any one item is unqualified, return "warnings: 1 = Length unqualified". warning2 = W1 not recognized; warning3 = Barcode1 not recognized; warning4 = Barcode2 not recognized

#Design a SAGSplit(inputFastq,CellBarn,FindBarcode=FindBarcode,warning=['Len','W1','BC'],filterWarning=True)
#Among them, FindBarcode requires the input of a function written by the user themselves. This function only takes in one read string and is required to return a string. This string returns the Barcode #string to which the read belongs during normal matching. Or the unmatched reason strings eg,"Erro","Len","W1","BC", in short, can return a string
#filterWarning means that if the Barcode string corresponding to the read is in the warning, such reads will not be processed directly. It defaults to True unless specified by the user

@timeit
def SAGSplit(inputFastq,CellBarn,FindBarcode=FindBarcode,warning=['Len','W1','BC'],filterWarning=True):


    #create directory:
    if not os.path.exists(CellBarn):
        os.makedirs(CellBarn,exist_ok=True)

    #If the length of the warning is 0, filterWarning will be automatically set to False and no warning judgment will be performed
    if len(warning) == 0:
        filterWarning = False
    else:
        unique_warning = list(set(warning))
        warning = {element:0 for element in unique_warning}

    # Determine whether inputfastq is single-ended or double-ended
    if type(inputFastq) == str and inputFastq.endswith('.fastq'):
        ReadsEnd = 'Single'
        inputFastq = inputFastq
    elif type(inputFastq) == list:
        if len(inputFastq) == 2 and inputFastq[0].endswith('_R1.fastq') and inputFastq[1].endswith('_R2.fastq'):
            ReadsEnd = 'Pair'
            inputFastq1 = inputFastq[0]
            inputFastq2 = inputFastq[1]

        elif len(inputFastq) == 2 and inputFastq[0].endswith('_R2.fastq') and inputFastq[1].endswith('_R1.fastq'):
            ReadsEnd = 'Pair'
            inputFastq1 = inputFastq[1]
            inputFastq2 = inputFastq[0]

        elif len(inputFastq) == 1 and inputFastq[0].endswith('.fastq'):
            ReadsEnd = 'Single'
            inputFastq = inputFastq[0]
        else:
            warnings.warn("The input fastq file does not meet the format requirements", UserWarning)
    else:
        warnings.warn("The input fastq file does not meet the format requirements", UserWarning)



    #Determine which method to distribute to the cells based on the single or double ends
    if ReadsEnd == 'Single':
        with open(inputFastq, 'r') as f:
            lines = []
            while True:
                line = f.readline()
                if len(line) == 0:
                    break

                if len(lines) == 4:
                    bc = FindBarcode(lines[1])

                    if filterWarning and bc in warning.keys():
                        warning[bc] += 1
                    else:
                        lines[0] = lines[0].replace('@','@'+str(bc)+':')
                        lines[2] = '+\n'
                        fastq = "".join(lines)
                        with open(os.path.join(CellBarn,str(bc)) + '.fastq', 'a') as o:
                            o.write(fastq)

                    lines = []
                    lines.append(line)

                else:
                    lines.append(line)

    elif ReadsEnd == 'Pair':
        with open(inputFastq1, 'r') as i1, open(inputFastq2, 'r') as i2:
            lines1 = []
            lines2 = []
            while True:
                line1 = i1.readline()
                line2 = i2.readline()
                if len(line1) == 0 or len(line2) == 0:
                    break

                if len(lines1) == 4:
                    bc = FindBarcode(lines1[1])
                    if filterWarning and bc in warning.keys():
                        warning[bc] += 1
                    else:
                        lines1[0] = lines1[0].replace('@','@'+str(bc)+':')
                        lines1[2] = '+\n'
                        lines2[0] = lines2[0].replace('@','@'+str(bc)+':')
                        lines2[2] = '+\n'
                        fastq1 = "".join(lines1)
                        fastq2 = "".join(lines2)
                        with open(os.path.join(CellBarn,str(bc)) + '_R1.fastq', 'a') as o1, open(os.path.join(CellBarn,str(bc)) + '_R2.fastq', 'a') as o2:
                            o1.write(fastq1)
                            o2.write(fastq2)

                    lines1 = []
                    lines1.append(line1)
                    lines2 = []
                    lines2.append(line2)
                else:
                    lines1.append(line1)
                    lines2.append(line2)

    print(warning)



PYTHONDIR = os.path.dirname(os.path.abspath(__file__))
toolDir=os.path.join(PYTHONDIR,"trimmomatic")

@timeit
def trim(inputFastqDir,trimDir,ReadsEnd='Single',ILLUMINACLIP='TruSeq3-PE.fa:2:30:10:3:TRUE',LEADING=25,TRAILING=3,SLIDINGWINDOW='4:20',MINLEN=30,threads=12):

    #Create the result file directory
    trimFastqDir = os.path.join(trimDir,'Fastq')
    trimUnpairedDir = os.path.join(trimDir,'Unpaired')
    trimLogDir = os.path.join(trimDir,'Log')
    trimErrorDir = os.path.join(trimDir,'Error')

    DirCreate = [trimDir, trimFastqDir, trimUnpairedDir, trimLogDir, trimErrorDir]
    for dir in DirCreate:
        if not os.path.exists(dir):
            os.makedirs(dir, exist_ok=True)
            

    try:
        if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
            raise ValueError('ReadsEnd must be Single or Pair!')
    except ValueError as e:
        print(e)

    if ReadsEnd=='Single':
        cells=[]
        file_list = os.listdir(inputFastqDir)
        for file_temp in file_list:
            if '.fastq' in file_temp:
                cells += [file_temp.replace('.fastq','')]         
        cells = list(set(cells))
        
        for cell in cells:
            CellInputFile=os.path.join(inputFastqDir,cell)+'.fastq'
            CellTrimFile=os.path.join(trimFastqDir,cell)+'.fastq'
            CellUnpairedFile=os.path.join(trimUnpairedDir,cell)+'_unpaired.fastq'
            CellErrorFile=os.path.join(trimErrorDir,cell)+'.trim.stderr.txt'
            command = f"(java -jar {toolDir}/trimmomatic-0.36.jar SE -threads {threads} -phred33 \
             -trimlog {trimLogDir}/{cell}.log \
              {CellInputFile} {CellTrimFile} {CellUnpairedFile} \
               ILLUMINACLIP:{toolDir}/{ILLUMINACLIP} \
                LEADING:{LEADING} \
                 TRAILING:{TRAILING} \
                  SLIDINGWINDOW:{SLIDINGWINDOW} \
                   MINLEN:{MINLEN}) \
                    2> {CellErrorFile}"
            subprocess.call(command, shell=True)

    elif ReadsEnd=='Pair':
        cells = []
        file_list = os.listdir(inputFastqDir)
        for file_temp in file_list:
            if '_R2.fastq' in file_temp:
                cells += [file_temp.replace('_R2.fastq', '')]
        cells = list(set(cells))
        
        for cell in cells:
            Cell1InputFile=os.path.join(inputFastqDir, cell)+'_R1.fastq'
            Cell2InputFile = os.path.join(inputFastqDir, cell) + '_R2.fastq'
            Cell1TrimFile=os.path.join(trimFastqDir,cell)+'_R1.fastq'
            Cell2TrimFile = os.path.join(trimFastqDir, cell) + '_R2.fastq'
            Cell1UnpairedFile=os.path.join(trimUnpairedDir,cell)+'_unpaired_R1.fastq'
            Cell2UnpairedFile = os.path.join(trimUnpairedDir, cell) + '_unpaired_R2.fastq'
            CellErrorFile=os.path.join(trimErrorDir,cell)+'.trim.stderr.txt'
            command = f"(java -jar {toolDir}/trimmomatic-0.36.jar PE -threads {threads} -phred33 \
             -trimlog {trimLogDir}/{cell}.log \
              {Cell1InputFile} {Cell2InputFile} {Cell1TrimFile} {Cell1UnpairedFile}  {Cell2TrimFile} {Cell2UnpairedFile} \
               ILLUMINACLIP:{toolDir}/{ILLUMINACLIP} \
                LEADING:{LEADING} \
                 TRAILING:{TRAILING} \
                  SLIDINGWINDOW:{SLIDINGWINDOW} \
                   MINLEN:{MINLEN}) \
                    2> {CellErrorFile}"
            subprocess.call(command, shell=True)