import os
import pandas as pd
import random
import subprocess
import time

'''
Script comment
1. This script requires the input of a "FastaDir" and an intermediate file of a tree, "Lujin TreeTemp", to return a tree file
2. The comment file Bin Phylum CellNum of the input box is used as the input for the itol drawing tree

'''
PHY_COLOR = {"Verrucomicrobiota":"#660000","Thermoproteota":"#783f04","Thermoplasmatota":"#7f6000","Spirochaetota":"#274e13","Proteobacteria":"#0c343d","Planctomycetota":"#073763","Patescibacteria":"#20124d","Myxococcota":"#4c1130","Marinisomatota":"#cc0000","Halobacteriota":"#e69138","Gemmatimonadota":"#f1c232","Firmicutes_C":"#6aa84f","Firmicutes_B":"#45818e","Firmicutes_A":"#3d85c6","Firmicutes":"#674ea7","Desulfobacterota_I":"#a64d79","Desulfobacterota":"#ea9999","Deinococcota":"#f9cb9c","Cyanobacteria":"#ffe599","Chloroflexota":"#b6d7a8","Campylobacterota":"#a2c4c9","Bacteroidota":"#9fc5e8","Actinobacteriota":"#b4a7d6","Acidobacteriota":"#d5a6bd","NA":"#000000","ERROR":"#ff0000"}



def timeit(func):
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = func(*args,**kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
        return result
    return wrapper





def getRandomColor():
    red =  random.randint(0,255)
    green = random.randint(0,255)
    blue = random.randint(0,255)
    color_code = "#{:02x}{:02x}{:02x}".format(red,green,blue)
    return color_code

@timeit
def BuildTree(FastaDir,TreeTemp,env=None):
   
    FastaAdjust = os.path.join(TreeTemp,'FastaAdjust/')
    DB = os.path.join(TreeTemp,'DB/')
    Anno = os.path.join(TreeTemp,'Anno/')

    
    DirCreat = [FastaAdjust, DB, Anno]
    for dir in DirCreat:
        os.makedirs(dir, exist_ok=True)

    
    files = os.listdir(FastaDir)
    fastafiles = []  #
    for filename in files:
        if filename.endswith('.fasta'):
            filename = filename[:-6]
            fastafiles.append(filename)

    
    for fastafile in fastafiles:
        inputfasta = os.path.join(FastaDir,fastafile) + '.fasta'
        outputdb = os.path.join(DB,fastafile) + '.db'
        command = 'anvi-gen-contigs-database -f ' + inputfasta + ' -o ' + outputdb
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command,shell=True)

    # Genomic files that cannot be constructed
    for fastafile in fastafiles:
        inputfasta = os.path.join(FastaDir, fastafile) + '.fasta'
        outputdb = os.path.join(DB, fastafile) + '.db'
        adjustfasta = os.path.join(FastaAdjust, fastafile) + '.fasta'
        adjustfastaname = os.path.join(FastaAdjust, fastafile) + '_name_conversion.txt'
        if not os.path.exists(outputdb):
            command = 'anvi-script-reformat-fasta ' + inputfasta + ' -o ' + adjustfasta + ' --simplify --report ' + adjustfastaname
            if env is not None:
                subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
            else:
                subprocess.call(command, shell=True)


            command = 'anvi-gen-contigs-database -f ' + adjustfasta +' -o ' + outputdb
            if env is not None:
                subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
            else:
                subprocess.call(command, shell=True)

    # build external.genomes.txt
    files = os.listdir(DB)
    DBfiles = []
    for filename in files:
        if filename.endswith('.db'):
            filename = filename[:-3]
            DBfiles.append(filename)

    externalfile = os.path.join(TreeTemp,'external_genomes.txt')
    concatenatedfile = os.path.join(TreeTemp,'concatenated-proteins_Bacteria_71_ribosomal6.fa')
    phylogenomicfile = os.path.join(TreeTemp, 'phylogenomic-tree_Bacteria_71_ribosomal6.txt')

    with open(externalfile, 'w') as input1:
        input1.write('name' + '\t' + 'contigs_db_path' + '\n')

    for db in DBfiles:
        with open(externalfile, 'a') as input1:
            input1.write(db + '\t' + os.path.join(DB,db) + '.db' + '\n')

    # Build an hmm model for each db
    for db in DBfiles:
        command = 'anvi-run-hmms -c ' + os.path.join(DB,db) + '.db'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

    command = 'anvi-get-sequences-for-hmm-hits --external-genomes ' + externalfile + ' -o ' + concatenatedfile + ' --hmm-source Bacteria_71 --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 --return-best-hit --get-aa-sequences --concatenate'
    if env is not None:
        subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
    else:
        subprocess.call(command, shell=True)

    command = 'anvi-gen-phylogenomic-tree -f ' + concatenatedfile + ' -o ' + phylogenomicfile
    if env is not None:
        subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
    else:
        subprocess.call(command, shell=True)


#BinAnno file format
# Bin Phylum CellNum Color The first type sets Color for each Phylum (this is relatively simple)
# Bin Phylum CellNum The second type, no Color is set for Phylum (randomly set the Color for uncommon Phylum, and directly select the color for common Phylum in PHY_COLOR)

@timeit
def itolPlot(BinAnno,Anno):
    
    if not os.path.exists(Anno):
        os.makedirs(Anno,exist_ok=True)

    BinAnno = pd.read_csv(BinAnno,sep='\t',header=0)
    #Create two itol Anno files：Label_Anno,bar_anno
    Label_Anno = os.path.join(Anno,'Label_Anno.txt')
    Bar_Anno = os.path.join(Anno,'Bar_Anno.txt')

    with open(Label_Anno, 'w') as input1:
        input1.write('TREE_COLORS\nSEPARATOR TAB\nDATA\n')

    with open(Bar_Anno, 'w') as input1:
        input1.write(
            'DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,simple bar testing\nDATASET_SCALE,100-1st line at 100-#5b5b5b-5-1-1,500-2nd line at 200-#5b5b5b-5-1-1,1000-3rd line at 1000-#5b5b5b-5-1-1,2000-4st line at 2000-#5b5b5b-5-1-1\nCOLOR,#999999\nWIDTH,1000\nMARGIN,0\nHEIGHT_FACTOR,1\nBAR_SHIFT,0\nBAR_ZERO,0\nDATA\n')

    #fill Label_Anno.txt
    if 'Color' in BinAnno.columns:
        for index, row in BinAnno.iterrows():
            phylum = row['Phylum']
            bin = row['Bin']
            color = row['Color']
            with open(Label_Anno,'a') as output:
                output.write(str(bin) + '\trange\t' + str(color) + '\t' + str(phylum) + '\n')
    else:
        for index, row in BinAnno.iterrows():

            phylum = row['Phylum']
            bin = row['Bin']

            if phylum in PHY_COLOR.keys():
                color = PHY_COLOR[phylum]
            else:
                #Assign a color to the uncommon Phylum and place it in the PHY_COLOR keyword
                while True:
                    color = getRandomColor()
                    if color not in PHY_COLOR.values():
                        PHY_COLOR[phylum] = color
                        break
                color = PHY_COLOR[phylum]

            with open(Label_Anno, 'a') as output:
                output.write(str(bin) + '\trange\t' + str(color) + '\t' + str(phylum) + '\n')


    #fill Bar_Anno.txt
    for index, row in BinAnno.iterrows():
        with open(Bar_Anno,'a') as output:
            output.write(str(row['Bin'])+','+str(row['CellNum'])+'\n')

@timeit
def GenerateBinAnno(fasta_dir, cell_anno_file, summary_file, output_file):
    print("Starting integration of BinAnno data...")

    # ==========================================
    # Step 1: Obtain the list of Fasta that are truly involved in the establishment (Information 1)
    # ==========================================
    valid_bins = []
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fasta') or filename.endswith('.fa'):
            # Remove the suffix and retain the core names of Bin/SGB
            bin_id = os.path.splitext(filename)[0]
            valid_bins.append(bin_id)
            
    valid_bins_set = set(valid_bins)
    print(f" -> Read {len(valid_bins_set)} valid sequences from FastaDir.")

    # ==========================================
    # Step 2: Calculate the CellNum for each SGB (Information 2)
    # ==========================================
    # Read CellAnno.txt (since the provided data has no header, we specify it manually)
    cell_df = pd.read_csv(cell_anno_file, sep='\t', header=None, names=['CellID', 'CellType', 'SGB'])
    
    # Automatically remove invalid data from DoubleCell, UnAsignedCell and NoSGB
    clean_cells = cell_df[~cell_df['CellType'].isin(['DoubleCell', 'UnAsignedCell'])]
    clean_cells = clean_cells[clean_cells['SGB'] != 'NoSGB']
    
    # Count the number of single cells contained in each SGB and convert it into a dictionary mapping
    sgb_counts = clean_cells['SGB'].value_counts().to_dict()
    print(f" -> After filtering invalid cells, successfully counted the CellNum for {len(sgb_counts)} valid SGBs.")

    # ==========================================
    # Step 3: Extract Phylum phylate-level classification information (Information 3)
    # ==========================================
    summary_df = pd.read_csv(summary_file, sep='\t', header=0)
    
    # After extracting p__ from the closest_taxonomy column using a regular expression, the semicolon; Previous content
    # r'p__([^;]+)' It means: Match p__ and capture all subsequent characters that are not semicolons
    summary_df['Phylum'] = summary_df['closest_taxonomy'].str.extract(r'p__([^;]+)')
    
    # Handle abnormal situations where regular expressions may not match (such as missing classification information), and fill in the default values
    summary_df['Phylum'] = summary_df['Phylum'].fillna('Unknown_Phylum')
    
    # Build a dictionary mapping from Bin to Phylum
    sgb_phylum = summary_df.set_index('Bin')['Phylum'].to_dict()
    print(" -> Successfully extracted Phylum information from summary.txt.")

    # ==========================================
    # Step 4: Cross-merge to generate the final BinAnno file
    # ==========================================
    output_data = []
    for bin_id in valid_bins_set:
        # If a Bin involved in the construction cannot be found in CellAnno, the default CellNum is 0
        cell_num = sgb_counts.get(bin_id, 0)
        # If a Bin that participated in the construction does not correspond to Phylum in the summary, enter "Unknown"
        phylum = sgb_phylum.get(bin_id, 'Unknown_Phylum')
        
        output_data.append({
            'Bin': bin_id,
            'Phylum': phylum,
            'CellNum': cell_num
        })
        
    final_df = pd.DataFrame(output_data)

    output_dir = os.path.dirname(output_file)
    # If there is a directory path and the directory has not been created yet, create it automatically (exist_ok=True to avoid reporting an error when it already exists).
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Save it as a tab-separated file for itolPlot to read
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f" Generation completed! BinAnno file has been saved to: {output_file}")
    print(final_df.head()) # Print the first few lines for checking