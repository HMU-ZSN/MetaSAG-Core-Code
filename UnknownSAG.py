import os
import time
import json
import pandas as pd
import subprocess
import warnings
import scipy.spatial.distance as sdist
from scipy.cluster.hierarchy import linkage,fcluster

def timeit(func):
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = func(*args,**kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
        return result
    return wrapper

@timeit
def ClusterSAG(inputFastq, outputDir, ReadsEnd='Single', SpadesEnv=None, SourmashEnv=None, ClusterThreshold=0.95, ClusterCriterion='distance'):

    round1 = os.path.join(outputDir,'Round1')
    CellFasta=os.path.join(round1,'CellFasta')
    sig = os.path.join(round1,'Sig')
    result_fastq = os.path.join(round1,'ResultFastq')
    TempWorkDir = os.path.join(result_fastq,'TempWorkDir')
    result_fasta = os.path.join(round1,'ResultFasta')
    ClusterFile = os.path.join(round1,'cluster.json')

    DirCreate = [outputDir, round1, sig, result_fastq, TempWorkDir, result_fasta]
    for dir in DirCreate:
        if not os.path.exists(dir):
            os.makedirs(dir, exist_ok=True)

    Cells=[]
    if ReadsEnd == 'Single':
        for file in os.listdir(inputFastq):
            if file.endswith('.fastq'):
                Cells.append(file.replace('.fastq',''))

        #step_1 Assemble each cell
        for cell in Cells:
            CellFastqFile = os.path.join(inputFastq, cell) + '.fastq'
            CellFastaSC = os.path.join(CellFasta, cell) + '_sc'
            CellFastaFile = os.path.join(CellFasta, cell) + '.fasta'

            command = f"spades.py --sc --pe1-s {CellFastqFile} -o {CellFastaSC}"
            if SpadesEnv is not None:
                subprocess.run(['conda', 'run', '-n', SpadesEnv, 'bash', '-c', command])
                subprocess.call('mv ' + CellFastaSC + '/contigs.fasta ' + CellFastaFile, shell=True)
                subprocess.call('rm -rf ' + CellFastaSC, shell=True)
            else:
                subprocess.call(command, shell=True)
                subprocess.call('mv ' + CellFastaSC + '/contigs.fasta ' + CellFastaFile, shell=True)
                subprocess.call('rm -rf ' + CellFastaSC, shell=True)


        #step_2 Obtain the.sig file of each cell
        CellSigMatrix = os.path.join(sig, 'round1_cmp')

        for cell in Cells:
            CellFastaFile = os.path.join(CellFasta, cell) + '.fasta'
            CellSigFile = os.path.join(sig, cell) + '.sig'


            command = f"sourmash compute --track-abundance {CellFastaFile} --output {CellSigFile}"
            if SourmashEnv is not None:
                subprocess.run(['conda', 'run', '-n', SourmashEnv, 'bash', '-c', command])
            else:
                subprocess.call(command, shell=True)


        command = f"sourmash compare {sig}/*.sig -k 51 -o {CellSigMatrix}.npy --csv {CellSigMatrix}.csv"
        if SourmashEnv is not None:
            subprocess.run(['conda', 'run', '-n', SourmashEnv, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        mat_file = pd.read_csv(CellSigMatrix + '.csv')
        mat_file=1.0-mat_file
        for i in range(len(mat_file)):
            mat_file.iloc[i][i] = 0

        square_form = sdist.squareform(mat_file)
        assignments = fcluster(linkage(square_form, method='complete'), ClusterThreshold, criterion=ClusterCriterion)

        cluster_result={}
        for j in range(max(assignments)):
            a_temp_list = []
            for i in range(len(assignments)):
                if assignments[i] == j + 1:
                    temp=mat_file.columns.values[i]
                    temp=str.split(temp,'/')[-1]
                    temp=str(temp[:-len('.fasta')])
                    a_temp_list += [temp]
                    #a_temp_list += [mat_file.columns.values[i][:-len('.fasta')]]
            cluster_result[str(j)] = a_temp_list

        with open(ClusterFile,'w') as out:
            json.dump(cluster_result,out)

        for bin,cells in cluster_result.items():

            result_fastq_bin = os.path.join(result_fastq, str(bin))

            for cell in cells:
                CellFastqFile = os.path.join(inputFastq, cell) + '.fastq'
                subprocess.call(f"cp {CellFastqFile} {TempWorkDir}",shell=True)

            subprocess.call(f"cat {TempWorkDir}/*.fastq > {result_fastq_bin}.fastq",shell = True)
            subprocess.call(f"rm {TempWorkDir}/*", shell=True)

            result_fasta_bin=os.path.join(result_fasta,str(bin))
            command = f'spades.py --sc --pe1-s {result_fastq_bin}.fastq -o {result_fasta_bin}_sc'
            if SpadesEnv is not None:
                subprocess.run(['conda', 'run', '-n', SpadesEnv, 'bash', '-c', command])
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)
            else:
                subprocess.call(command, shell=True)
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)


    elif ReadsEnd == 'Pair':
        for file in os.listdir(inputFastq):
            if file.endswith('_R1.fastq'):
                Cells.append(file.replace('_R1.fastq', ''))  

        # step_1 Assemble each cell
        for cell in Cells:
            CellFastqFile1 = os.path.join(inputFastq, cell) + '_R1.fastq'
            CellFastqFile2 = os.path.join(inputFastq, cell) + '_R2.fastq'
            CellFastaSC = os.path.join(CellFasta, cell) + '_sc'
            CellFastaFile = os.path.join(CellFasta, cell) + '.fasta'

            command = f"spades.py --sc --pe1-1 {CellFastqFile1} --pe1-2 {CellFastqFile2} -o {CellFastaSC}"
            if SpadesEnv is not None:
                subprocess.run(['conda', 'run', '-n', SpadesEnv, 'bash', '-c', command])
                subprocess.call('mv ' + CellFastaSC + '/contigs.fasta ' + CellFastaFile, shell=True)
                subprocess.call('rm -rf ' + CellFastaSC, shell=True)
            else:
                subprocess.call(command, shell=True)
                subprocess.call('mv ' + CellFastaSC + '/contigs.fasta ' + CellFastaFile, shell=True)
                subprocess.call('rm -rf ' + CellFastaSC, shell=True)

        # step_2 Obtain the.sig file of each cell
        CellSigMatrix = os.path.join(sig, 'round1_cmp')

        for cell in Cells:
            CellFastaFile = os.path.join(CellFasta, cell) + '.fasta'
            CellSigFile = os.path.join(sig, cell) + '.sig'

            command = f"sourmash compute --track-abundance {CellFastaFile} --output {CellSigFile}"
            if SourmashEnv is not None:
                subprocess.run(['conda', 'run', '-n', SourmashEnv, 'bash', '-c', command])
            else:
                subprocess.call(command, shell=True)

        command = f"sourmash compare {sig}/*.sig -k 51 -o {CellSigMatrix}.npy --csv {CellSigMatrix}.csv"
        if SourmashEnv is not None:
            subprocess.run(['conda', 'run', '-n', SourmashEnv, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        mat_file = pd.read_csv(CellSigMatrix + '.csv')
        mat_file = 1.0 - mat_file
        for i in range(len(mat_file)):
            mat_file.iloc[i][i] = 0

        square_form = sdist.squareform(mat_file)
        assignments = fcluster(linkage(square_form, method='complete'), ClusterThreshold, criterion=ClusterCriterion)

        cluster_result = {}
        for j in range(max(assignments)):
            a_temp_list = []
            for i in range(len(assignments)):
                if assignments[i] == j + 1:
                    temp = mat_file.columns.values[i]
                    temp = str.split(temp, '/')[-1]
                    temp = str(temp[:-len('.fasta')])
                    a_temp_list += [temp]
                    #a_temp_list += [mat_file.columns.values[i][:-len('.fasta')]]
            cluster_result[str(j)] = a_temp_list

        with open(ClusterFile, 'w') as out:
            json.dump(cluster_result, out)

        for bin, cells in cluster_result.items():

            result_fastq_bin = os.path.join(result_fastq, str(bin))

            for cell in cells:
                CellFastqFile1 = os.path.join(inputFastq, cell) + '_R1.fastq'
                CellFastqFile2 = os.path.join(inputFastq, cell) + '_R2.fastq'

                subprocess.call(f"cp {CellFastqFile1} {TempWorkDir}", shell=True)
                subprocess.call(f"cp {CellFastqFile2} {TempWorkDir}", shell=True)

            subprocess.call(f"cat {TempWorkDir}/*_R1.fastq > {result_fastq_bin}_R1.fastq", shell=True)
            subprocess.call(f"cat {TempWorkDir}/*_R2.fastq > {result_fastq_bin}_R2.fastq", shell=True)
            subprocess.call(f"rm {TempWorkDir}/*", shell=True)

            result_fasta_bin = os.path.join(result_fasta, str(bin))

            command = f'spades.py --sc --pe1-1 {result_fastq_bin}_R1.fastq --pe1-2 {result_fastq_bin}_R2.fastq -o {result_fasta_bin}_sc'
            if SpadesEnv is not None:
                subprocess.run(['conda', 'run', '-n', SpadesEnv, 'bash', '-c', command])
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)
            else:
                subprocess.call(command, shell=True)
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)

    else:
        warnings.warn("The ReadsEnd parameter can only be Single or Pair", UserWarning)

@timeit
def ClusterBin(OldRoundFold, NewRoundFold, ReadsEnd='Single', SpadesEnv=None, SourmashEnv=None, ClusterThreshold=0.95, ClusterCriterion='distance'):

    with open(os.path.join(OldRoundFold,'cluster.json')) as input:
        OldCluster = json.load(input)

    OldFastaDir = os.path.join(OldRoundFold,'ResultFasta')
    OldFastqDir = os.path.join(OldRoundFold, 'ResultFastq')

    NewSig=os.path.join(NewRoundFold,'Sig')
    NewResultFastq=os.path.join(NewRoundFold,'ResultFastq')
    NewResultFasta=os.path.join(NewRoundFold,'ResultFasta')
    TempWorkDir = os.path.join(NewResultFastq, 'TempWorkDir')
    NewClusterFile_forBin = os.path.join(NewRoundFold,'ClusterFile_forBin.json')
    NewClusterFile_forSag = os.path.join(NewRoundFold,'ClusterFile_forSag.json')

    DirCreate = [NewRoundFold, NewSig, NewResultFastq, NewResultFasta, TempWorkDir]
    for dir in DirCreate:
        if not os.path.exists(dir):
            os.makedirs(dir, exist_ok=True)

    OldBins=[]
    #Run sourmash for each Fasta file to obtain the Sig file
    for file in os.listdir(OldFastaDir):
        if file.endswith('.fasta'):
            OldBins.append(str(file.replace('.fasta','')))


    BinSigMatrix = os.path.join(NewSig, 'NewRound_cmp')

    for bin in OldBins:
        BinFastaFile = os.path.join(OldFastaDir, bin) + '.fasta'
        BinSigFile = os.path.join(NewSig, bin) + '.sig'

        command = f"sourmash compute --track-abundance {BinFastaFile} --output {BinSigFile}"
        if SourmashEnv is not None:
            subprocess.run(['conda', 'run', '-n', SourmashEnv, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

    command = f"sourmash compare {NewSig}/*.sig -k 51 -o {BinSigMatrix}.npy --csv {BinSigMatrix}.csv"
    if SourmashEnv is not None:
        subprocess.run(['conda', 'run', '-n', SourmashEnv, 'bash', '-c', command])
    else:
        subprocess.call(command, shell=True)


    mat_file = pd.read_csv(BinSigMatrix + '.csv')
    mat_file = 1.0 - mat_file
    for i in range(len(mat_file)):
        mat_file.iloc[i][i] = 0

    square_form = sdist.squareform(mat_file)
    assignments = fcluster(linkage(square_form, method='complete'), ClusterThreshold, criterion=ClusterCriterion)

    cluster_bin = {}
    for j in range(max(assignments)):
        a_temp_list = []
        for i in range(len(assignments)):
            if assignments[i] == j + 1:
                temp = mat_file.columns.values[i]
                temp = str.split(temp, '/')[-1]
                temp = str(temp[:-len('.fasta')])
                a_temp_list += [temp]
                #a_temp_list += [mat_file.columns.values[i][:-len('.fasta')]]
        cluster_bin[str(j)] = a_temp_list

    with open(NewClusterFile_forBin, 'w') as out:
        json.dump(cluster_bin, out)


    cluster_sag = {}
    if ReadsEnd == 'Single':
        for newbin,oldbins in cluster_bin.items():
            newbin=str(newbin)

            if newbin not in cluster_sag.keys():
                cluster_sag[newbin] = []
            for oldbin in oldbins:
                cluster_sag[newbin] += list(OldCluster[oldbin])
                OldBinFastq=os.path.join(OldFastqDir, oldbin)
                subprocess.call(f"cp {OldBinFastq}.fastq {TempWorkDir}", shell=True)

            result_fastq_bin = os.path.join(NewResultFastq, newbin)
            subprocess.call(f"cat {TempWorkDir}/*.fastq > {result_fastq_bin}.fastq", shell=True)
            subprocess.call(f"rm {TempWorkDir}/*", shell=True)


            result_fasta_bin = os.path.join(NewResultFasta, newbin)
            command = f'spades.py --sc --pe1-s {result_fastq_bin}.fastq -o {result_fasta_bin}_sc'
            if SpadesEnv is not None:
                subprocess.run(['conda', 'run', '-n', SpadesEnv, 'bash', '-c', command])
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)
            else:
                subprocess.call(command, shell=True)
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)

        with open(NewClusterFile_forSag, 'w') as out:
            json.dump(cluster_sag, out)

    elif ReadsEnd == 'Pair':
        for newbin, oldbins in cluster_bin.items():
            if newbin not in cluster_sag.keys():
                cluster_sag[newbin] = []
            for oldbin in oldbins:
                cluster_sag[newbin] += list(OldCluster[oldbin])
                OldBinFastq = os.path.join(OldFastqDir, oldbin)
                subprocess.call(f"cp {OldBinFastq}_R1.fastq {TempWorkDir}", shell=True)
                subprocess.call(f"cp {OldBinFastq}_R2.fastq {TempWorkDir}", shell=True)

            result_fastq_bin = os.path.join(NewResultFastq, newbin)
            subprocess.call(f"cat {TempWorkDir}/*_R1.fastq > {result_fastq_bin}_R1.fastq", shell=True)
            subprocess.call(f"cat {TempWorkDir}/*_R2.fastq > {result_fastq_bin}_R2.fastq", shell=True)
            subprocess.call(f"rm {TempWorkDir}/*", shell=True)

            result_fasta_bin = os.path.join(NewResultFasta, newbin)
            command = f'spades.py --sc --pe1-1 {result_fastq_bin}_R1.fastq --pe1-2 {result_fastq_bin}_R2.fastq -o {result_fasta_bin}_sc'
            if SpadesEnv is not None:
                subprocess.run(['conda', 'run', '-n', SpadesEnv, 'bash', '-c', command])
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)
            else:
                subprocess.call(command, shell=True)
                subprocess.call(f'mv {result_fasta_bin}_sc/contigs.fasta {result_fasta_bin}.fasta', shell=True)
                subprocess.call(f'rm -rf {result_fasta_bin}_sc', shell=True)

        with open(NewClusterFile_forSag, 'w') as out:
            json.dump(cluster_sag, out)

    else:
        warnings.warn("The ReadsEnd parameter can only be Single or Pair", UserWarning)


@timeit
def FilterUnassignedCells(AnnotationFile, inputFastqDir, OutputDir):
    """
    Based on the filter of the comment file, link the FASTQ file corresponding to the cell with Type 'UnAsignedCell' to the output directory through a [soft link]. 
    Parameters: 
    AnnotationFile (str): The path of the annotation file (TSV format, including Cell, Type, and SGB columns) 
    inputFastqDir (str): The path of the input folder for storing FASTQ files split by cells. 
    OutputDir (str): The folder path where the output result is stored. The soft link will be automatically named 'unassigned_ original file name '.
    """
    
    # Automatically create the output directory (if it does not exist)
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir, exist_ok=True)

    # ==========================================
    # Part One: Read the comment file and create a hash table
    # ==========================================
    try:
        anno_df = pd.read_csv(AnnotationFile, sep='\t')
        target_cells = set(anno_df[anno_df['Type'] == 'UnAsignedCell']['Cell'])
        print(f"[Info] Successfully loaded the annotation file; a total of {len(target_cells)} unknown cells were identified for extraction.")
    except Exception as e:
        raise ValueError(f"Failed to read annotation file: {e}")

    # ==========================================
    # Part Two: Scan the folder and establish a soft connection
    # ==========================================
    total_files = 0
    linked_files = 0

    for file_name in os.listdir(inputFastqDir):
        if not file_name.endswith('.fastq'):
            continue  # Ignore non-FastQ files
            
        total_files += 1
        
        # Extract Cell ID
        base_name = file_name.replace('.fastq', '')
        if base_name.endswith('_R1') or base_name.endswith('_R2'):
            cell_id = base_name[:-3]
        else:
            cell_id = base_name

        # ==========================================
        # Part Three: If matched, establish an absolute path soft connection
        # ==========================================
        if cell_id in target_cells:
            # Key point: The absolute path of the source file must be obtained; otherwise, the soft connection is prone to failure
            src_path = os.path.abspath(os.path.join(inputFastqDir, file_name))
            dest_path = os.path.join(OutputDir, f"unassigned_{file_name}")
            
            # If a file or soft link already exists in the target path, delete it first (to prevent duplicate running and error FileExistsError)
            if os.path.exists(dest_path) or os.path.islink(dest_path):
                os.remove(dest_path)
            
            # Establish a soft connection (the Linux command is equivalent to: ln -s src_path dest_path)
            os.symlink(src_path, dest_path)
            linked_files += 1

    # statistical result

    print(f"[Success] Filtering and symbolic link creation completed!")
    print(f"Total number of FASTQ files scanned: {total_files}")
    print(f"Number of symbolic links successfully created: {linked_files}")
    print(f"All files have been stored in: {OutputDir}")