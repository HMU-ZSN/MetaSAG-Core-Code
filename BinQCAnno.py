import time
import math
import os
import pandas as pd
import numpy as np
import shutil
import subprocess
import logging
from pathlib import Path

'''
Script comment
This script requires the input of a fasta folder.
1. Provide the SGB number corresponding to each fasta file
2. Use the SGB comment file of MetaPhlan4 that comes with this Python package to categorize and name each fasta file
3. In the Checkm section, select whether to perform a secondary Bandage correction. If a secondary calibration is to be carried out, a fasta folder with the Bandage that originally exceeded the pollution standard but met the quality standards after calibration must be provided


'''

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

def custom_log(level, msg):
    import time
    current_time = time.strftime('%H:%M:%S', time.localtime())
    print(f"[{current_time} - {level}] {msg}")



PYTHONDIR = os.path.dirname(os.path.abspath(__file__))
def removeShortContig(inputfastafile,outputfastafile,minlen=500):
    fasta={}
    with open(inputfastafile) as input1:
        while True:
            line = input1.readline()
            if len(line) == 0:
                break

            line = line.replace('\n', '')

            if line.startswith('>'):
                flag = line
                fasta[flag] = []
            else:
                fasta[flag] += [line]
    input1.close()

    # contigs with a length of 500 or more
    lt_500 = []
    for key in fasta.keys():
        l = int(key.split('_')[3])
        if l >= minlen:
            lt_500 += [key]

    # Cover more than two standard deviations
    cov = []
    for node in lt_500:
        cov += [float(node.split('_')[5])]

    cov_log = [math.log(num) for num in cov]
    pj = np.mean(cov_log)
    bz = np.std(cov_log)
    threshold = pj - 2 * bz

    cov_save = []
    for i in range(len(lt_500)):
        if cov_log[i] > threshold:
            cov_save += [lt_500[i]]

    with open(outputfastafile, 'a') as output:
        for key, value in fasta.items():
            if key in cov_save:
                output.write(key + '\n')
                for line in value:
                    output.write(line + '\n')

    output.close()

def getContigRC(inputfastafile,outputfastafile,outputstatistic,readsLen = 300,ContigRCThreshold = 100):
    fasta = {}
    with open(inputfastafile) as input1:
        while True:
            line = input1.readline()
            if len(line) == 0:
                break

            line = line.replace('\n', '')

            if line.startswith('>'):
                flag = line
                fasta[flag] = []
            else:
                fasta[flag] += [line]
    input1.close()

    Node = []
    Length = []
    Cov = []
    ContigRC = []
    for key in fasta.keys():
        node = key
        l = int(key.split('_')[3])
        cov = float(key.split('_')[5])
        Node += [node]
        Length += [l]
        Cov += [cov]

        temp = (l/readsLen)*cov
        ContigRC += [temp]

    ContigRCdf = pd.DataFrame({'Node':Node,'Length':Length,'Cov':Cov,'ContigRC':ContigRC})
    ContigRCdf.to_csv(outputstatistic,index=False,sep='\t')
    CongtigSaver = ContigRCdf[(ContigRCdf['ContigRC'] >= ContigRCThreshold) & (ContigRCdf['Length'] >= readsLen)]
    ContigSaver = CongtigSaver['Node'].tolist()

    with open(outputfastafile, 'a') as output:
        for key, value in fasta.items():
            if key in ContigSaver:
                output.write(key + '\n')
                for line in value:
                    output.write(line + '\n')
    return ContigRCdf


@timeit
def Summary(FastaSGB,CheckmFile,GTDBFile,outputSummary=None):

    mpa = pd.read_csv(os.path.join(PYTHONDIR,'mpa_vOct22_CHOCOPhlAnSGB_202403_species.txt'),sep='\t',header=None)
    mpa.columns=['SGB','SGBName']

    FastaSGB = pd.read_csv(FastaSGB,sep='\t',header=0)
    FastaSGB.columns = ['Bin','SGB']
    Checkm = pd.read_csv(CheckmFile,header=0,sep='\t')
    GTDB = pd.read_csv(GTDBFile,header=0,sep='\t')

    GTDB = GTDB[['user_genome','closest_placement_reference','closest_placement_taxonomy','closest_placement_ani']]
    GTDB.columns = ['Bin','closest_reference','closest_taxonomy','closest_ani']
    Checkm = Checkm[['Bin Id','Completeness','Contamination']]
    Checkm.columns = ['Bin','Completeness','Contamination']

    SummaryResult = pd.merge(FastaSGB,mpa,on='SGB',how='left')
    SummaryResult = pd.merge(SummaryResult,Checkm,on='Bin',how='left')
    SummaryResult = pd.merge(SummaryResult,GTDB,on='Bin',how='left')

    #Fill all cells that do not match with 'no'
    SummaryResult.fillna('no',inplace=True)

    SummaryResult['ContamFlag'] = None
    SummaryResult['CompleteFlag'] = None
    SummaryResult['QC'] = None
    for index, row in SummaryResult.iterrows():
        if row['Contamination'] != 'no':
            if row['Contamination'] > 10:
                SummaryResult.loc[index,'ContamFlag'] = 'fail'
            else:
                SummaryResult.loc[index, 'ContamFlag'] = 'pass'
        else:
            SummaryResult.loc[index, 'ContamFlag'] = 'no'


        if row['Completeness'] != 'no':
            if row['Completeness'] < 50:
                SummaryResult.loc[index,'CompleteFlag'] = 'fail'
            else:
                SummaryResult.loc[index, 'CompleteFlag'] = 'pass'
        else:
            SummaryResult.loc[index, 'CompleteFlag'] = 'no'


        if row['Contamination'] != 'no' and row['Completeness'] != 'no':
            if (row['Contamination'] <= 5) & (row['Completeness'] >= 90):
                SummaryResult.loc[index,'QC'] = 'High'
            elif (row['Contamination'] <= 10) & (row['Completeness'] >= 50):
                SummaryResult.loc[index, 'QC'] = 'Medium'
            else:
                SummaryResult.loc[index, 'QC'] = 'Low'
        else:
            SummaryResult.loc[index, 'QC'] = 'no'



    '''
    SummaryResult['ContamFlag'] = np.where(SummaryResult['Contamination'] > 10,'fail','pass')
    SummaryResult['CompleteFlag'] = np.where(SummaryResult['Completeness'] < 50,'fail','pass')

    conditions = [
        (SummaryResult['Contamination'] <= 5) & (SummaryResult['Completeness'] >= 90),
        (SummaryResult['Contamination'] <= 10) & (SummaryResult['Completeness'] >= 50)
    ]
    choices = ['High','Medium']
    SummaryResult['QC'] = np.select(conditions,choices,default='Low')
    '''

    if outputSummary is not None:
        SummaryResult.to_csv(outputSummary,sep='\t',index=False)

    return SummaryResult


@timeit
def FastaQC1(FastaDir,FastaOut,minlen=500):
    #This function only performs the operation of removing short contigs
    #create directory
    if not os.path.exists(FastaOut):
        os.makedirs(FastaOut,exist_ok=True)

    files = os.listdir(FastaDir)
    for filename in files:
        if filename.endswith('.fasta'):
            removeShortContig(os.path.join(FastaDir, filename),os.path.join(FastaOut, filename),minlen=minlen)


@timeit
def FastaQC2(FastaDir,FastaOut,CheckmFile,FastgDir=None):
    #This function:
    #(1)For fasta files with short contigs removed, based on the Checkm result file, return which files are qualified genomes [Pass], which files are unqualified but can be Bandage corrected [Bandage], 
    # and which genome files can be directly abandoned due to their small completeness [Abandon]
    #(2)If the fastg folder can be provided, copy the corresponding fastg file to the Bandage folder; If the fastg folder is not provided, only place the corresponding fasta file in the Bandage folder
    ##(3)The distribution logic of the Bandage has been refined and multiple outputs are supported: the original Bandage classification has been split into low(highly contaminated) and 
    ## medium(moderately contaminated). High-quality sketches that meet specific conditions (such as completeness >=90% and 5%<= contamination <=10%) are allowed to be simultaneously output to Pass 
    ## and Bandage/medium to balance the automation advancement of the main process and the subsequent manual refinement requirements.

    #create directory
    Pass = os.path.join(FastaOut,'Pass')
    Abandon = os.path.join(FastaOut,'Abandon')
    Bandage_low = os.path.join(FastaOut,'Bandage','low')
    Bandage_medium = os.path.join(FastaOut,'Bandage','medium')

    DirCreate = [FastaOut,Pass,Abandon,Bandage_low,Bandage_medium]
    for dir in DirCreate:
        if not os.path.exists(dir):
            os.makedirs(dir,exist_ok=True)

    Checkm = pd.read_csv(CheckmFile, header=0, sep='\t')
    Checkm = Checkm[['Bin Id', 'Completeness', 'Contamination']]
    Checkm.columns = ['Bin', 'Completeness', 'Contamination']

    for index,row in Checkm.iterrows():
        Bin = row['Bin']
        Completeness = float(row['Completeness'])
        Contamination = float(row['Contamination'])

        source_fasta = os.path.join(FastaDir,Bin)+'.fasta'
        target_dirs = []

        if Completeness >= 50 and Contamination <= 10:
            target_dirs.append(Pass)
        if Completeness > 90 and 5 < Contamination < 10:
            target_dirs.append(Bandage_medium)
        if Completeness >= 50 and Contamination > 10:
            target_dirs.append(Bandage_low)
        if Completeness < 50:
            target_dirs.append(Abandon)

        for target_dir in target_dirs:
            target_fasta = os.path.join(target_dir, Bin) + '.fasta'
            shutil.copy(source_fasta, target_fasta)

            if target_dir in [Bandage_low, Bandage_medium] and FastgDir is not None:
                source_fastg = os.path.join(FastgDir, Bin) + '.fastg'
                target_fastg = os.path.join(target_dir, Bin) + '.fastg'
                if os.path.exists(source_fastg):
                    shutil.copy(source_fastg, target_fastg)

@timeit
def FastaCLQC(FastaDir, FastaOut, readsLen=300, ContigRCThreshold=100):
    # 1. Define the subdirectory path
    StacDir = os.path.join(FastaOut, 'Stac')
    
    # 2. Fix Bug: Create all directories correctly
    DirCreate = [FastaOut, StacDir]
    for d in DirCreate:
        os.makedirs(d, exist_ok=True) # exist_ok=True already includes the logic of "create if it doesn't exist"

    # 3. Traverse the file and process it
    files = os.listdir(FastaDir)
    for filename in files:
        if filename.endswith('.fasta'):
            # Fix the Bug: Specify the correct output directory StacDir and add the.tsv suffix
            StacFile = os.path.join(StacDir, filename[:-6] + '.tsv')
            
            # Execute the core filtering function
            getContigRC(
                os.path.join(FastaDir, filename),
                os.path.join(FastaOut, filename),
                StacFile,
                readsLen=readsLen,
                ContigRCThreshold=ContigRCThreshold
            )



def _run_metabat2_tnf(fasta_path: Path, output_prefix: Path, env: str = None, threads: int = 16, 
                      min_contig: int = 1500, max_edges: int = 500, min_cls_size: int = 100000) -> bool:
    """Internal helper function: Perform MetaBAT 2 secondary boxing with relaxed restrictions"""
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    command = (
        f"metabat2 -i {fasta_path} -o {output_prefix} -t {threads} "
        f"-m {min_contig} --maxEdges {max_edges} --minClsSize {min_cls_size} -v"
    )
    
    if env is not None:
        cmd = f"conda run -n {env} --no-capture-output {command}"
    else:
        cmd = command

    try:
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"MetaBAT 2 run failed for {fasta_path.name}, exit code: {e.returncode}")
        return False
    except FileNotFoundError:
        logging.error("Command 'metabat2' not found! Please ensure the environment is activated.")
        return False

@timeit
def RebinMetaBAT2(InputDir, OutputDir, MetaBat2_env=None, threads=16, 
                      min_contig=1500, max_edges=500, min_cls_size=100000):
    """
    Perform secondary purification on the fasta files within the specified input directory and summarize the generated bins to the new output directory
    
    parameters:
        InputDir (str/Path): The path of the fasta folder that needs to be processed
        OutputDir (str/Path): Summarize the path of the new folder output
        threads (int): Number of running threads (default 16)
        min_contig (int): Discard contigs shorter than this length (default 1500)
        max_edges (int): The upper limit of the number of edges connected to the graph can be relaxed to accommodate more fragments (default 500)
        min_cls_size (int): The minimum volume threshold of Bin (default: 100,000)
    """
    logging.basicConfig(level=logging.INFO, format='[%(asctime)s - %(levelname)s] %(message)s', datefmt='%H:%M:%S', force=True)

    input_path = Path(InputDir)
    export_path = Path(OutputDir)
    
    export_path.mkdir(parents=True, exist_ok=True)
    
    fasta_files = list(input_path.glob("*.fasta"))
    if not fasta_files:
        logging.warning(f"No .fasta files found in {InputDir}! Please check the path.")
        return

    logging.info(f"Found {len(fasta_files)} fasta files to process, starting secondary refinement...")

    for fasta_file in fasta_files:
        sample_name = fasta_file.stem  
        
        tmp_work_dir = input_path / f"{sample_name}_tmp_workspace"
        output_prefix = tmp_work_dir / "bin"
        
        # Pass the parameters passed in from outside to the internal metabat2 function
        success = _run_metabat2_tnf(fasta_file, output_prefix, 
                                    env=MetaBat2_env,
                                    threads=threads, 
                                    min_contig=min_contig, 
                                    max_edges=max_edges, 
                                    min_cls_size=min_cls_size)
        
        if success:
            bin_files = list(tmp_work_dir.glob("bin.*.fa"))
            
            if bin_files:
                for bin_file in bin_files:
                    bin_num = bin_file.stem.split('.')[-1]
                    new_filename = f"{sample_name}_rebin{bin_num}.fasta"
                    dest_file = export_path / new_filename
                    shutil.copy2(bin_file, dest_file)
            else:
                logging.warning(f"Sample {sample_name} failed to generate new Bins.")
                
        if tmp_work_dir.exists():
            shutil.rmtree(tmp_work_dir)

    logging.info("=" * 40)
    logging.info(f"Secondary refinement completed! Results are in: {export_path.resolve()}")

@timeit
def FilterMetaBAT2(InputDir, OutputDir, CheckmFile, Level, max_het=40):
    """
    Ultimate distribution script: Supports intersection extraction, lifts the To_Improve contamination limit, and provides a lenient rescue mechanism for the low group
    """
    if Level not in ['low', 'medium']:
        logging.error("The Level parameter must be either 'low' or 'medium'!")
        return

    Pass = os.path.join(OutputDir, 'Pass') 
    To_Improve = os.path.join(OutputDir, 'To_Improve') 
    Abandon = os.path.join(OutputDir, 'Abandon') 
    
    for d in [OutputDir, Pass, To_Improve, Abandon]:
        if not os.path.exists(d):
            os.makedirs(d, exist_ok=True)

    try:
        checkm_df = pd.read_csv(CheckmFile, header=0, sep='\t')
        checkm_df = checkm_df[['Bin Id', 'Completeness', 'Contamination', 'Strain heterogeneity']]
        checkm_df.columns = ['Bin', 'Completeness', 'Contamination', 'Heterogeneity']
    except Exception as e:
        logging.error(f"Unable to read the CheckM file: {e}")
        return

    checkm_df['Original_Sample'] = checkm_df['Bin'].apply(lambda x: x.split('_rebin')[0])
    checkm_df['Score'] = checkm_df['Completeness'] - 5 * checkm_df['Contamination']
    best_bins = checkm_df.sort_values(
        ['Original_Sample', 'Score'], ascending=[True, False]
    ).drop_duplicates('Original_Sample')

    logging.info(f"Successfully identified the best Bin for {len(best_bins)} samples. Starting distribution according to the latest rules...")

    counts = {'Pass': 0, 'Improve': 0, 'Intersection': 0, 'Abandon': 0}

    for index, row in best_bins.iterrows():
        bin_id = row['Bin']
        comp = float(row['Completeness'])
        cont = float(row['Contamination'])
        het = float(row['Heterogeneity'])
        orig_name = row['Original_Sample']

        source_fasta = os.path.join(InputDir, bin_id + '.fasta')
        if not os.path.exists(source_fasta):
            continue

        target_dirs = []
        is_pass = False
        is_improve = False

        # ==========================================
        # Rule 1: Pass (Standard Available at Any Time
        # As long as the integrity is greater than or equal to 50 and the contamination degree is less than or equal to 10
        # ==========================================
        if comp >= 50 and cont <= 10:
            target_dirs.append(Pass)
            is_pass = True

        # ==========================================
        # Rule 2: To_Improve (Refinement Potential Standard)
        # Basic threshold: Heterogeneity <= 40
        # Condition A (General) : Completeness >= 90 and contamination degree > 5 (no upper limit)
        # Condition B (for the low group) : Completeness >= 50 and contamination degree > 10 (no upper limit)
        # ==========================================
        if het <= max_het:
            if (comp >= 90 and cont > 5) or (Level == 'low' and comp >= 50 and cont > 10):
                target_dirs.append(To_Improve)
                is_improve = True

        # ==========================================
        # Rule 3: Abandon (Completely abandon the standard
        # All those who do not meet the Pass requirements and are not qualified to improve will be eliminated
        # ==========================================
        if not target_dirs:
            target_dirs.append(Abandon)
            counts['Abandon'] += 1

        # statistical data
        if is_pass: counts['Pass'] += 1
        if is_improve: counts['Improve'] += 1
        if is_pass and is_improve: counts['Intersection'] += 1

        # Perform a copy and restore the name
        for target_dir in target_dirs:
            target_fasta = os.path.join(target_dir, orig_name + '.fasta')
            shutil.copy(source_fasta, target_fasta)

# The final streamlined statistical output
    custom_log("INFO", "=" * 40)
    custom_log("INFO", f"[{Level} Group] Distribution Summary:")
    custom_log("INFO", f"  -> Selected for Pass group:        {counts['Pass']}")
    custom_log("INFO", f"  -> Selected for To_Improve group:  {counts['Improve']}")
    custom_log("INFO", f"     (Number in the intersection: {counts['Intersection']})")
    custom_log("INFO", f"  -> Assigned to Abandon group:      {counts['Abandon']}")

@timeit
def MergeFasta(dir1: str, dir2: str, output_dir: str):
    """
    Merge the fasta files from the two folders into a third folder (uniformly using the Secure Copy mode).
    Files in dir2 have a higher priority. When encountering files with the same name, they will automatically overwrite files with the same name in dir1.
    """
    path1 = Path(dir1)
    path2 = Path(dir2)
    out_path = Path(output_dir)
    
    # Create an output directory
    out_path.mkdir(parents=True, exist_ok=True)
    
    count1, count2, overwritten = 0, 0, 0
    
    # ==========================================
    # Step 1: Copy the file of the first path (base)
    # ==========================================
    if path1.exists():
        for fasta in path1.glob("*.fasta"):
            dest = out_path / fasta.name
            shutil.copy2(fasta, dest)
            count1 += 1
    else:
        logging.warning(f"Path1 does not exist or is misspelled: {path1}")
        
    # ==========================================
    # Step 2: Copy the file in the second path (priority overwriting)
    # ==========================================
    if path2.exists():
        for fasta in path2.glob("*.fasta"):
            dest = out_path / fasta.name
            
            # Record overwrite event
            if dest.exists():
                overwritten += 1
                logging.debug(f"Overwriting file: {fasta.name}") 
            
            shutil.copy2(fasta, dest)
            count2 += 1
    else:
        logging.warning(f"Path2 does not exist or is misspelled: {path2}")

    # ==========================================
    # Step 3: Output statistical results
    # ==========================================
    logging.info("=== FASTA file merging completed ===")
    logging.info(f" -> Copied from Path1: {count1} files")
    logging.info(f" -> Copied from Path2: {count2} files")
    if overwritten > 0:
        logging.info(f" -> [Notice] Path2 successfully overwrote {overwritten} old-version files!")
    logging.info(f" -> Final results are stored in: {out_path.resolve()}\n")