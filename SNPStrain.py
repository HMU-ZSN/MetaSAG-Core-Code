import shutil
from typing import Dict

import pandas as pd
import os
import subprocess
import json
import numpy as np
import umap
from sklearn import preprocessing
import time
import pandas as pd
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree, fcluster
import matplotlib
import matplotlib.pyplot as plt
import os
import random
import matplotlib.lines as mlines
import matplotlib as mpl
from pylab import *
import json
from ast import literal_eval
from matplotlib.cm import ScalarMappable
import scipy
import scipy.stats
import scipy.signal
from matplotlib.ticker import MaxNLocator
import subprocess
import sys
import pickle
from collections import defaultdict


'''
Script comment
The strain SNP script provides two strategies:
1. Strain only one species of cell box at a time. The cell box should contain the fastq file of each cell as well as the assembly genome file of the box
2. Strain all the cells in the boxes at once. It is required to uniformly provide the path of the CellBarn folder where the assembly genome files and the fastq files of each cell are stored. [Provide optional parameters to indicate which boxes in all cells are subjected to strain separation operations]
'''


def timeit(func):
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = func(*args,**kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
        return result
    return wrapper




class SingleBin():

    def __init__(self, BinDir, ResultDir, SpeciesName):
        self.BinDir = BinDir
        self.ResultDir = ResultDir
        self.SpeciesName = SpeciesName


    strain_color_list = ["#33bbff", "#0050e6", "#009999", "#777777", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                         "#8c564b",
                         "#e377c2"]
    umap_plot_shape = ['D', 'o', 's', 'v', 'd', 'p', '*', 'H', '+', 'x']
    default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2',
                          u'#7f7f7f',
                          u'#bcbd22', u'#17becf']

    def strain_split1(self,work_dirt_temp='', file_name='', cutoff_cell_ratio=0.05, cutoff_type_count=2,
                      cutoff_type_ratio=0.01,
                      cutoff_SNP_ratio=0.01):
        barcode_list = []
        contig_list = []
        location_on_contig_list = []
        reference_list = []
        alter_list = []
        sequence_list = []

        total_line = 0
        with open(os.path.join(work_dirt_temp, file_name)) as finput:
            while True:
                line_read = finput.readline()
                if len(line_read) == 0:
                    break
                total_line += 1
                if total_line == 4:
                    parsed_bam = line_read.split()
                    for i in range(len(parsed_bam)):
                        barcode_with_bam = parsed_bam[i]
                        if '.bam' == barcode_with_bam[-4:]:
                            barcode_list += [barcode_with_bam[:-4]]
                        cutoff_cell_count = int(cutoff_cell_ratio * len(barcode_list))
                        cutoff_type_count = max(cutoff_type_count, int(cutoff_type_ratio * len(barcode_list)))

                if line_read[0] != '#':
                    line_split_temp = line_read.split()
                    contig_list += [int(line_split_temp[0].split('_')[1])]
                    location_on_contig_list += [int(line_split_temp[1])]
                    reference_list += [line_split_temp[3]]
                    alter_list += [line_split_temp[4]]

                    mapping_count_list_temp = []
                    for i in range(9, len(line_split_temp)):
                        mapping_result_temp = line_split_temp[i].split(':')[-1]
                        reference_count_temp = int(mapping_result_temp.split(',')[0])
                        alter_count_temp = int(mapping_result_temp.split(',')[1])
                        total_count_temp = reference_count_temp + alter_count_temp
                        if total_count_temp < 2:
                            mapping_count_list_temp += [0]
                        elif max(reference_count_temp, alter_count_temp) / float(total_count_temp) < 0.99:
                            mapping_count_list_temp += [0]
                        else:
                            if max(reference_count_temp, alter_count_temp) == reference_count_temp:
                                mapping_count_list_temp += [1]
                            elif max(reference_count_temp, alter_count_temp) == alter_count_temp:
                                mapping_count_list_temp += [-1]
                            else:
                                print('someting wrong')
                    sequence_list += [mapping_count_list_temp]

        #print(len(sequence_list))
        #print(sequence_list[0])

        
        new_columns = {str(i + 1): sequence_list[i] for i in range(len(sequence_list))}
        snp_pd = pd.DataFrame(new_columns)
        
        snp_pd = snp_pd.T
        #print(snp_pd)
        #print(barcode_list)
        snp_pd.columns = barcode_list

        #print(snp_pd.shape)
        #print(len(barcode_list))
        #print(barcode_list)

        snp_pd['contig_list'] = contig_list
        snp_pd['location_on_contig_list'] = location_on_contig_list
        snp_pd['reference_list'] = reference_list
        snp_pd['alter_list'] = alter_list
        snp_pd.to_csv(os.path.join(work_dirt_temp, file_name) + '_snp.csv', index=False)

        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name) + '_snp.csv')

        to_drop_row_list = []
        for i in range(len(snp_pd)):
            if sum(snp_pd.iloc[i][:-4] != 0) < cutoff_cell_count:
                to_drop_row_list.append(i)
                continue
            if sum(snp_pd.iloc[i][:-4] == -1) < cutoff_type_count or sum(snp_pd.iloc[i][:-4] == 1) < cutoff_type_count:
                to_drop_row_list.append(i)
        print('File name is {},total SNPs before drop is {},SNPs to drop is {}.'.format(file_name, len(snp_pd),
                                                                                        len(to_drop_row_list)))

        self.DropSNPNum = len(to_drop_row_list)
        self.AllSNPNum = len(snp_pd)

        snp_pd = snp_pd.drop(to_drop_row_list, axis=0)
        snp_pd = snp_pd.reset_index(drop=True)
        snp_pd.to_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered.csv', index=False)

        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered.csv')
        cutoff_SNP_count = max(cutoff_SNP_ratio * len(snp_pd), 10)
        to_drop_column_list = []
        barcode_list = list(snp_pd.columns.values[:-4])
        for i in range(snp_pd.shape[1] - 4):
            if sum(snp_pd[barcode_list[i]] != 0) < cutoff_SNP_count:
                to_drop_column_list += [barcode_list[i]]

        print("Total cells before drop is {},cell to drop is {}.".format(snp_pd.shape[1] - 4, len(to_drop_column_list)))
        self.AllCellNum = snp_pd.shape[1] - 4
        self.DropCellNum = len(to_drop_column_list)

        snp_pd = snp_pd.drop(to_drop_column_list, axis=1)
        snp_pd.to_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered_bc_dropped.csv', index=False)

        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered_bc_dropped.csv')
        
        total_ratio_of_snps_list = []
        for i in range(len(snp_pd)):
            total_coverage_temp = sum(list(snp_pd.iloc[i][:-4] != 0))
            if total_coverage_temp == 0:
                print(i)
                total_ratio_of_snps_list += [1]
                continue
            total_ratio_of_snps_list += [sum(list(snp_pd.iloc[i][:-4] == 1)) / float(total_coverage_temp)]
            
        
        plt.figure()
        count_temp = plt.hist(np.array(total_ratio_of_snps_list), bins=np.arange(0, 1.01, 0.02))
        plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
        plt.xlabel('Ratio of cells matching strain 1')
        plt.ylabel('Number of alleles')
        plt.savefig(os.path.join(work_dirt_temp, file_name) + '_SNP.svg', format='svg', bbox_inches='tight')
        plt.clf()
        return (total_ratio_of_snps_list)

    def check_strain1(self, file_name, dirt_temp, output_dirt, tree_cut_threshold=2.7, min_cells=5):
        strain_color_list = ["#33bbff", "#0050e6", "#009999", "#777777", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                             "#8c564b", "#e377c2"]
        umap_plot_shape = ['D', 'o', 's', 'v', 'd', 'p', '*', 'H', '+', 'x']
        default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2',
                              u'#7f7f7f', u'#bcbd22', u'#17becf']

        # work_dirt_temp=dirt_base+'fastq/'+bin+'/bamResult/variants/'
        work_dirt_temp = dirt_temp
        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name))
        #print(snp_pd.shape)  
        snp_pd_temp = snp_pd
        snp_pd_temp = snp_pd_temp.drop(['contig_list', 'location_on_contig_list', 'reference_list', 'alter_list'],
                                       axis=1)
        snp_pd_temp = snp_pd_temp.T
        
        # Increase the interception for all SGB that have been filtered out or have very few remaining cells
        SpeciesName = file_name.split('_')[0]
        
        # Note: At this point, snp_pd_temp has been transposited. shape[0] represents the remaining number of cells, and shape[1] represents the number of SNPS
        current_cell_count = snp_pd_temp.shape[0] 
        
        if snp_pd_temp.empty or current_cell_count < min_cells or snp_pd_temp.shape[1] == 0:
            print(f"  [QC Intercept] {SpeciesName}: Only {current_cell_count} valid cells remain after quality control (minimum required: {min_cells}); insufficient for dimensionality reduction and clustering. Skipping downstream analysis.")
            
            # Generate empty files in the corresponding folders as markers to facilitate the summary report records
            with open(os.path.join(output_dirt, "QC_FAILED_NO_VALID_DATA.flag"), 'w') as f:
                f.write(f"{SpeciesName} QC failed: Only {current_cell_count} valid cells remain after filtering, which is insufficient for UMAP dimensionality reduction analysis.\n")
                
            return None, None, None
        

        X = preprocessing.normalize(snp_pd_temp)
        linked = linkage(X, method='ward')
        # print(linked)
        matplotlib.rcParams.update({'font.size': 50})
        fig, ax = plt.subplots(1, 1, figsize=(3, 9))

        plt.yticks([0, 50, 100, 150, 200])
        plt.tick_params(axis='both', direction="in", which='both', right=True, length=15, width=3)
        matplotlib.rcParams['lines.linewidth'] = 0.8

        # strain_color_list_BV = [strain_color_list[3], strain_color_list[2], strain_color_list[1],strain_color_list[0]]###############################################################
        strain_color_list_BV = [strain_color_list[2], strain_color_list[0]]
        hierarchy.set_link_color_palette(strain_color_list_BV)

        # tree_cut_threshold=1.8
        dendrogram_plot = dendrogram(linked,
                                     orientation='left',
                                     distance_sort='ascending',
                                     leaf_font_size=8,
                                     no_labels=True,
                                     color_threshold=tree_cut_threshold,
                                     above_threshold_color='black',
                                     show_leaf_counts=False)

        plt.xticks([])
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.savefig(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity_tree.svg',
                    format='svg',
                    bbox_inches='tight')
        matplotlib.rcParams.update({'font.size': 12})

        plt.figure()
        X = snp_pd_temp
        X = preprocessing.normalize(snp_pd_temp)
        reducer = umap.UMAP(random_state=42)
        embedding = reducer.fit_transform(X)
        #print(embedding.shape)
        # print(embedding)
        plt.scatter(embedding[:, 0], embedding[:, 1], s=3)

        plt.savefig(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_umap.svg', format='svg',
                    bbox_inches='tight')
        plt.clf()

        index_int_list = []

        for i in dendrogram_plot['ivl']:
            index_int_list.append(int(i))
        SNP_similarity_df = pd.DataFrame(-1.0, index=snp_pd.columns[index_int_list],
                                         columns=snp_pd.columns[index_int_list])
        for i in range(len(SNP_similarity_df)):
            column_name_a = snp_pd.columns[index_int_list[i]]
            for j in range(len(SNP_similarity_df)):
                if i == j:
                    SNP_similarity_df.iloc[i, j] = 1.0
                    continue
                column_name_b = snp_pd.columns[index_int_list[j]]
                check_snp_temp = snp_pd[[column_name_a, column_name_b]]
                check_snp_temp = check_snp_temp.loc[check_snp_temp[column_name_b] != 0]
                check_snp_temp = check_snp_temp.loc[check_snp_temp[column_name_a] != 0]

                if len(check_snp_temp) == 0:
                    continue
                SNP_similarity_df.iloc[i, j] = sum(
                    check_snp_temp[column_name_b] == check_snp_temp[column_name_a]) / float(
                    len(check_snp_temp))

        SNP_similarity_df.to_csv(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity_df.csv',
                                 index=False)
        fig, ax = plt.subplots()
        data_color = []
        fig.set_size_inches(6, 6, forward=True)
        # my_cmap=plt.cm.get_cmap(matplotlib.cm.cmap_d.keys()[58])
        # my_cmap = mpl.colors.ListedColormap(['w', default_color_list[4]])
        my_cmap = 'RdBu'
        maximum_similarity = 1

        # print(matplotlib.cm.cmap_d.keys()[58])
        SNP_similarity_df = pd.read_csv(
            os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity_df.csv')
        ax.imshow(SNP_similarity_df, cmap=my_cmap, vmax=maximum_similarity, vmin=0)
        ax.set_xlabel("Single cells", fontsize=15)
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0, maximum_similarity))
        sm.set_array([])
        cbar = plt.colorbar(sm, orientation='vertical', fraction=0.05, pad=0.05, aspect=15,
                            ticks=[0, maximum_similarity])
        cbar.set_label('Similarity of SNPs between cells', rotation=90, labelpad=-1, fontsize=15)
        plt.savefig(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity.svg', format='svg',
                    bbox_inches='tight')
        plt.clf()
        return (snp_pd, linked, embedding)

    @timeit
    def SetupRunDir(self, fastaDir, fastqDir, cellAnno, ReadsEnd='Single'):
        
        import sys

        
        try:
            if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
                raise ValueError('ReadsEnd must be Single or Pair!')
        except ValueError as e:
            print(f" [Error] {e}")
            sys.exit()

        
        print(f"\n>> Setting up an independent runtime directory for {self.SpeciesName} (SetupRunDir)...")

       
        if not os.path.exists(self.BinDir):
            os.makedirs(self.BinDir, exist_ok=True)
        if not os.path.exists(self.ResultDir):
            os.makedirs(self.ResultDir, exist_ok=True)

        
        try:
            CellAnno_df = pd.read_csv(cellAnno, sep='\t', header=0)
        except Exception as e:
            print(f" [Error] Failed to read the cell annotation table: {e}")
            sys.exit()

        if 'Cluster' not in CellAnno_df.columns:
            if 'SGB' in CellAnno_df.columns:
                CellAnno_df.rename(columns={'SGB': 'Cluster'}, inplace=True)
            else:
                print(" [Error] The CellAnno file must contain either a 'Cluster' or 'SGB' column!")
                sys.exit()

        
        target_cells = CellAnno_df[CellAnno_df['Cluster'].astype(str) == str(self.SpeciesName)]['Cell'].tolist()
        
        if len(target_cells) == 0:
            print(f" [Error] Terminated: No cell records belonging to {self.SpeciesName} were found in the annotation table {cellAnno}!")
            sys.exit()

       
        sourceFasta = os.path.join(fastaDir, self.SpeciesName) + '.fasta'
        targetFasta = os.path.join(self.BinDir, self.SpeciesName) + '.fasta'
        
        if not os.path.exists(sourceFasta):
           
            print(f" [Error] Terminated: Required reference genome file not found in {fastaDir}: {self.SpeciesName}.fasta")
            sys.exit()
            
        if not os.path.exists(targetFasta):
            os.symlink(sourceFasta, targetFasta)

       
        print(f">> A total of {len(target_cells)} cells belonging to {self.SpeciesName} were found; starting to create symbolic links for sequences...")
        
        linked_count = 0
        missing_cells = []
        
        if ReadsEnd == 'Single':
            for cell in target_cells:
                cell = str(cell)
                src_fq = os.path.join(fastqDir, cell) + '.fastq'
                tgt_fq = os.path.join(self.BinDir, cell) + '.fastq'
                
                if os.path.exists(src_fq):
                    if not os.path.exists(tgt_fq):
                        os.symlink(src_fq, tgt_fq)
                    linked_count += 1
                else:
                    missing_cells.append(cell)
                    
        elif ReadsEnd == 'Pair':
            for cell in target_cells:
                cell = str(cell)
                src_fq1 = os.path.join(fastqDir, cell) + '_R1.fastq'
                src_fq2 = os.path.join(fastqDir, cell) + '_R2.fastq'
                tgt_fq1 = os.path.join(self.BinDir, cell) + '_R1.fastq'
                tgt_fq2 = os.path.join(self.BinDir, cell) + '_R2.fastq'
                
                if os.path.exists(src_fq1) and os.path.exists(src_fq2):
                    if not os.path.exists(tgt_fq1):
                        os.symlink(src_fq1, tgt_fq1)
                    if not os.path.exists(tgt_fq2):
                        os.symlink(src_fq2, tgt_fq2)
                    linked_count += 1
                else:
                    missing_cells.append(cell)

       
        if len(missing_cells) > 0:
            
            print(f" [Warning] Fastq files for {len(missing_cells)} cells were not found in {fastqDir} and have been skipped.")
            
        if linked_count == 0:
          
            print(f" [Error] No corresponding Fastq files were found for any cells; failed to set up the working directory!")
            sys.exit()
            
        
        print(f" Runtime directory setup successful! Data for {linked_count} cells has been prepared. Current directory: {self.BinDir}")
        return True

    
    @timeit
    def SingleBinPrepare(self, bcftools=None, snap_aligner=None, samtools=None, env=None, ReadsEnd='Single', min_cells=5, sort_tool='samtools'):

        BinDir = self.BinDir
        ResultDir = self.ResultDir
        SpeciesName = self.SpeciesName

        try:
            if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
                raise ValueError('ReadsEnd must be Single or Pair!')
            if sort_tool not in ['samtools', 'snap']:
                raise ValueError("Currently, only 'samtools' or 'snap' are supported as sorting software.")
        except ValueError as e:
            print(e)
            sys.exit()

       
        ResultBam = os.path.join(ResultDir, 'ResultBam')
        ResultVCF = os.path.join(ResultDir, 'ResultVCF')
        ResultPic = os.path.join(ResultDir, 'ResultPic')

        self.PrepareBam = ResultBam
        self.PrepareVCF = ResultVCF
        self.PreparePic = ResultPic

        DirCreate = [ResultDir, ResultBam, ResultVCF, ResultPic]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        if bcftools is None:
            bcftools = 'bcftools'

        if snap_aligner is None:
            snap_aligner = 'snap-aligner'
            
        if samtools is None:
            samtools = 'samtools'

        FastaFile = os.path.join(BinDir, SpeciesName) + '.fasta'
        FastaIndex = os.path.join(ResultBam, SpeciesName) + '_index'
        
        command = snap_aligner + ' index ' + FastaFile + ' ' + FastaIndex
        subprocess.run(['conda', 'run', '-n', 'base', 'bash', '-c', command])

        
        FastqList = []
        if ReadsEnd == 'Single':
            for filename in os.listdir(BinDir):
                if filename.endswith('.fastq'):
                    FastqList.append(filename[:-6])

            for cell in FastqList:
                CellFastq = os.path.join(BinDir, cell) + '.fastq'
                CellBam = os.path.join(ResultBam, cell) + '.bam'

                if sort_tool == 'samtools':
                    
                    UnsortedBam = CellBam[:-4] + '_unsorted.bam'
                    RG_info = 'ID:FASTQ\tPL:Illumina\tPU:pu\tLB:lb\tSM:' + cell
                    
                    cmd_snap = [snap_aligner, 'single', FastaIndex, CellFastq, '-I', '-R', RG_info, '-t', '8', '-o', UnsortedBam]
                    cmd_sort = [samtools, 'sort', '-@', '8', '-o', CellBam, UnsortedBam]
                    cmd_rm = ['rm', '-f', UnsortedBam]
                    
                    if env is not None:
                        subprocess.run(['conda', 'run', '-n', env] + cmd_snap, check=True)
                        subprocess.run(['conda', 'run', '-n', env] + cmd_sort, check=True)
                        subprocess.run(['conda', 'run', '-n', env] + cmd_rm, check=True)
                    else:
                        subprocess.run(cmd_snap, check=True)
                        subprocess.run(cmd_sort, check=True)
                        subprocess.run(cmd_rm, check=True)
                elif sort_tool == 'snap':
                    
                    command = f"{snap_aligner} single {FastaIndex} {CellFastq} -I -R ID:FASTQ\\tPL:Illumina\\t:PU:pu\\tLB:lb\\tSM:{cell} -t 8 -so -o {CellBam}"
                    if env is not None:
                        subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                    else:
                        subprocess.call(command, shell=True)
    

        elif ReadsEnd == 'Pair':
            for filename in os.listdir(BinDir):
                if filename.endswith('_R1.fastq'):
                    FastqList.append(filename[:-9])

            for cell in FastqList:
                CellFastq_R1 = os.path.join(BinDir, cell) + '_R1.fastq'
                CellFastq_R2 = os.path.join(BinDir, cell) + '_R2.fastq'
                CellBam = os.path.join(ResultBam, cell) + '.bam'

                if sort_tool == 'samtools':
                   
                    UnsortedBam = CellBam[:-4] + '_unsorted.bam'
                    RG_info = 'ID:FASTQ\tPL:Illumina\tPU:pu\tLB:lb\tSM:' + cell
                    
                    cmd_snap = [snap_aligner, 'paired', FastaIndex, CellFastq_R1, CellFastq_R2, '-I', '-R', RG_info, '-t', '8', '-o', UnsortedBam]
                    cmd_sort = [samtools, 'sort', '-@', '8', '-o', CellBam, UnsortedBam]
                    cmd_rm = ['rm', '-f', UnsortedBam]
                    
                    if env is not None:
                        subprocess.run(['conda', 'run', '-n', env] + cmd_snap, check=True)
                        subprocess.run(['conda', 'run', '-n', env] + cmd_sort, check=True)
                        subprocess.run(['conda', 'run', '-n', env] + cmd_rm, check=True)
                    else:
                        subprocess.run(cmd_snap, check=True)
                        subprocess.run(cmd_sort, check=True)
                        subprocess.run(cmd_rm, check=True)
                elif sort_tool == 'snap':
                    
                    command = f"{snap_aligner} paired {FastaIndex} {CellFastq_R1} {CellFastq_R2} -I -R ID:FASTQ\\tPL:Illumina\\t:PU:pu\\tLB:lb\\tSM:{cell} -t 8 -so -o {CellBam}"
                    if env is not None:
                        subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                    else:
                        subprocess.call(command, shell=True)
                

        ###### Use bcftools to obtain the vcf file [bcftools] ######
        command = bcftools + ' mpileup -f ' + FastaFile + ' -a FORMAT/AD ' + ResultBam + '/*.bam > ' + ResultVCF + '/raw_calls.bcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        command = bcftools + ' call --ploidy 1 -v -m ' + ResultVCF + '/raw_calls.bcf > ' + ResultVCF + '/calls.vcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        command = bcftools + ' view -v snps ' + ResultVCF + '/calls.vcf > ' + ResultVCF + '/calls_snps.vcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        command = bcftools + " filter -o " + ResultVCF + '/' + SpeciesName + "_calls_snps_qual30.vcf -i 'QUAL>30' " + ResultVCF + '/calls_snps.vcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        command = "rm -rf " + ResultVCF + '/raw_calls.bcf ' + ResultBam + '/*.bam'
        subprocess.call(command,shell=True)

        ######
        self.strain_split1(ResultVCF, SpeciesName + '_calls_snps_qual30.vcf')
        
        #Check whether the SGB contains effective cells that have passed quality control
        result = self.check_strain1(SpeciesName + '_calls_snps_qual30.vcf_snp_filtered_bc_dropped.csv',
                                    ResultVCF, 
                                    ResultPic,
                                    min_cells=min_cells)

        if result[0] is None:
            return False  # Terminate the current function prematurely and report "Failure (False)" to AllBinPrepare.
        snp_pd, linked, embedding = result


        plt.clf()

        with open(ResultVCF + '/' + 'snp_pd.pkl', 'wb') as f:
            pickle.dump(snp_pd, f)

        with open(ResultVCF + '/' + 'linked.pkl', 'wb') as f:
            pickle.dump(linked, f)

        with open(ResultVCF + '/' + 'embedding.pkl', 'wb') as f:
            pickle.dump(embedding, f)

        self.snp_pd = snp_pd
        self.linked = linked
        self.embedding = embedding
        self.ReadsEnd = ReadsEnd

        return True

    @timeit
    def SingleBinSplit(self,ClusterNum):
        ClusterNum = int(ClusterNum)
        Species = self.SpeciesName
        PrepareDir = self.ResultDir

        self.PrepareBam = os.path.join(PrepareDir, 'ResultBam')
        self.PrepareVCF = os.path.join(PrepareDir, 'ResultVCF')
        self.PreparePic = os.path.join(PrepareDir, 'ResultPic')

        PrepareBam = self.PrepareBam
        PrepareVCF = self.PrepareVCF
        PreparePic = self.PreparePic

      
        if not hasattr(self, 'snp_pd') or self.snp_pd is None:
            print(f" Restoring preprocessed data for {Species} from disk...")
            try:
                # Note that here it must be self.snp_pd, forcing the data back into the object's mind
                with open(os.path.join(PrepareVCF, 'snp_pd.pkl'), 'rb') as f:
                    self.snp_pd = pickle.load(f)
                with open(os.path.join(PrepareVCF, 'linked.pkl'), 'rb') as f:
                    self.linked = pickle.load(f)
                with open(os.path.join(PrepareVCF, 'embedding.pkl'), 'rb') as f:
                    self.embedding = pickle.load(f)
            except FileNotFoundError:
                print(f" Error: PKL cache file not found in {PrepareVCF}!")
                print(" Please confirm whether SingleBinPrepare has been successfully run for this species.")
                return False



        snp_pd_old = self.snp_pd
        snp_pd = self.snp_pd
        linked = self.linked
        embedding = self.embedding

        cutree_contigs = cut_tree(linked, n_clusters=ClusterNum)
        umap_compare = pd.DataFrame({'Barcode': snp_pd.columns[:-4]})
        umap_compare['X'] = embedding[:, 0]
        umap_compare['Y'] = embedding[:, 1]

        cluster_hier = []
        for i in range(len(cutree_contigs)):
            cluster_hier.append(cutree_contigs[i][0])

        umap_compare['cluster_hier'] = cluster_hier

        for j in range(ClusterNum):
            umap_compare_bin = umap_compare.loc[umap_compare['cluster_hier'] == j]
            plt.scatter(umap_compare_bin['X'], umap_compare_bin['Y'], s=10, c=SingleBin.strain_color_list[j],
                        marker=SingleBin.umap_plot_shape[j])


        for j in range(ClusterNum):
            #Calculate the center position of each cluster
            cluster_points = umap_compare.loc[umap_compare['cluster_hier'] == j]
            x_center = cluster_points['X'].mean()
            y_center = cluster_points['Y'].mean()
            plt.text(x_center,y_center,str(j),fontsize=12,fontweight='bold',color='black',bbox=dict(facecolor='white',alpha=0.8,edgecolor='none'))
        plt.xticks([])
        plt.yticks([])
        plt.ylabel('UMAP dimention 2')
        plt.xlabel('UMAP dimention 1')
        plt.savefig(os.path.join(PreparePic, Species) + '_sp_call_SNP_umap_final.svg', format='svg',bbox_inches='tight')
        plt.clf()

        ###################################
        ###get strain share SNP genotype###
        ###################################

        # file = os.path.join(PrepareVCF,Species) + '_calls_snps_qual30.vcf_snp_filtered_bc_dropped.csv'

        for i in range(ClusterNum):
            index_of_this_cluster = []
            for j in range(len(cutree_contigs)):
                if cutree_contigs[j] == i:
                    index_of_this_cluster += [j]
                snp_pd['hierarchy_cluster_' + str(i) + 'pos'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == 1).sum(axis=1)
                snp_pd['hierarchy_cluster_' + str(i) + 'neg'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == -1).sum(axis=1)

        strain_X_genotype = [[] for i in range(ClusterNum)]  # Create an empty two-dimensional array of length ClusterNum
        for i in range(len(snp_pd)):
            for j in range(ClusterNum):
                # generate genotype of each strain
                # at each SNP, if both have less than 2 cells, assign it 0, means not sure
                # if the max genotype is less than 90% total, assign 0
                # otherwise assign the major genotype
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] <= 1 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][
                            i] <= 1:
                    strain_X_genotype[j].append(0)
                    continue
                elif snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] != 0 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i] != 0:
                    if max(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                           snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]) / float(
                        min(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                            snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i])) < 9:
                        strain_X_genotype[j].append(0)
                        continue
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] > snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]:
                    strain_X_genotype[j].append(1)
                else:
                    strain_X_genotype[j].append(-1)

        for j in range(ClusterNum):
            snp_pd['strain_' + str(j) + '_genotype'] = strain_X_genotype[j]

        rows_to_drop = []
        for i in range(len(snp_pd)):
            genotype_SNP_temp = []
            for j in range(ClusterNum):
                genotype_SNP_temp.append(snp_pd['strain_' + str(j) + '_genotype'][i])
            if not (1 in genotype_SNP_temp and -1 in genotype_SNP_temp):
                rows_to_drop.append(i)

        print("SNPs before drop is {}, SNPs to drop is {}.".format(len(snp_pd), len(rows_to_drop)))

        snp_pd = snp_pd.drop(rows_to_drop, axis=0)
        snp_pd = snp_pd.reset_index(drop=True)

        # this is to generate strain list
        list_of_ratio = [[] for i in range(ClusterNum)]
        barcode_list_strain = [[] for i in range(ClusterNum + 1)]

        for i in range(snp_pd.shape[1]):
            column_name_temp = snp_pd.columns[i]
            if column_name_temp.startswith('strain') or column_name_temp.startswith('hierarchy_cluster'):
                continue

            if 'list' not in column_name_temp:
                y = column_name_temp.split('/')[-1]
            else:
                continue

            # if not column_name_temp in snp_pd.columns:
            # continue
            temp_df_check = snp_pd.loc[snp_pd[column_name_temp] != 0]
            for j in range(ClusterNum):
                temp_df_check_strain = temp_df_check.loc[temp_df_check['strain_' + str(j) + '_genotype'] != 0]
                total_count_same = sum(temp_df_check_strain[column_name_temp] == \
                                       temp_df_check_strain['strain_' + str(j) + '_genotype'])
                list_of_ratio[j].append(total_count_same / float(len(temp_df_check_strain[column_name_temp])))

            list_of_ratio_check = [list_of_ratio[j][-1] for j in range(ClusterNum)]
            print(list_of_ratio_check)

            if max(list_of_ratio_check) < 0.95:
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            elif sum(np.array(list_of_ratio_check) > 0.95) > 1:
                print('Not sure,more then two strains match', column_name_temp)
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            else:
                # barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(column_name_temp)
                barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(y)


        for cluster in range(len(barcode_list_strain)):
            # print('look here is barcode_list_strain:', i)
            # print(barcode_list_strain[i])
            if cluster == (len(barcode_list_strain) - 1):
                for cell in barcode_list_strain[cluster]:
                    with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'a') as output:
                        output.write('Mist'+str(cluster) + '\t' + str(cell) + '\n')
            else:
                for cell in barcode_list_strain[cluster]:
                    with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'a') as output:
                        output.write(str(cluster) + '\t' + str(cell) + '\n')


        self.snp_pd = snp_pd_old
        self.snp_pd_strain = snp_pd
        self.barcode_list_strain = barcode_list_strain
        self.StrainCellsFile = os.path.join(PrepareDir, 'StrainCells.txt')

    def StrainAssem(self, StrainAssemDir="", env=None, StrainCellsFile="", ReadsEnd=None):

        Species = self.SpeciesName  

       
        if ReadsEnd is not None:
            self.ReadsEnd = ReadsEnd
        elif not hasattr(self, 'ReadsEnd'):
            print(" [Warning] ReadsEnd was not specified and no memory record was found; defaulting to 'Single'.")
            self.ReadsEnd = 'Single'

      
        if StrainAssemDir == "":
            SpeciesStrainAssemDir = os.path.join(self.ResultDir, 'StrainAssem')
        else:
            if not os.path.exists(StrainAssemDir):
                os.makedirs(StrainAssemDir, exist_ok=True)
            SpeciesStrainAssemDir = os.path.join(StrainAssemDir, Species + '_StrainAssem')

        if not os.path.exists(SpeciesStrainAssemDir):
            os.makedirs(SpeciesStrainAssemDir, exist_ok=True)

      
        if StrainCellsFile == "":
            StrainCellsFile = self.StrainCellsFile

        clusters = defaultdict(list)
        with open(StrainCellsFile) as f:
            for line in f:
                if line.strip():
                    parts = line.split()
                    if len(parts) >= 2:
                        c, cell = parts[0], parts[1]
                        if c.lower() == "cluster":
                            continue
                        if "Mist" not in c:
                            clusters[c].append(cell)

        if self.ReadsEnd == 'Single':
            for cl, cells in clusters.items():
                with open(f'{SpeciesStrainAssemDir}/{cl}.fastq', 'w') as out:
                    for cell in cells:
                        CellFastq = os.path.join(self.BinDir, cell) + '.fastq'
                        if os.path.exists(CellFastq):
                            with open(CellFastq) as in_f:
                                out.write(in_f.read())

          
            for cluster_name in clusters.keys():
                perfect_prefix = f"{Species}@{cluster_name}"
                command = f'spades.py --sc --pe1-s {SpeciesStrainAssemDir}/{cluster_name}.fastq -o {SpeciesStrainAssemDir}/{perfect_prefix}_sc'
                try:
                    if env is not None:
                        result = subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                        success = (result.returncode == 0)
                    else:
                        return_code = subprocess.call(command, shell=True)
                        success = (return_code == 0)

                    if success:
                        print(f"cluster {perfect_prefix} assem successfully")
                        os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}.fastq')
                        
                        import shutil
                        source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/scaffolds.fasta'
                        final_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}.fasta'
                        if not os.path.exists(source_fasta):
                            source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/contigs.fasta'
                        if os.path.exists(source_fasta):
                            shutil.copy(source_fasta, final_fasta)
                            print(f" -> Successfully extracted the final genome: {perfect_prefix}.fasta")
                    else:
                        print(f"cluster {perfect_prefix} failure")

                except Exception as e:
                    print(f"execute spades.py command failed: {e}")
                    success = False


        elif self.ReadsEnd == 'Pair':
            for cl, cells in clusters.items():
                with open(f'{SpeciesStrainAssemDir}/{cl}_R1.fastq','w') as out1, open(f'{SpeciesStrainAssemDir}/{cl}_R2.fastq','w') as out2:
                    for cell in cells:
                        CellFastq_R1 = os.path.join(self.BinDir, cell) + '_R1.fastq'
                        CellFastq_R2 = os.path.join(self.BinDir, cell) + '_R2.fastq'
                        if os.path.exists(CellFastq_R1):
                            with open(CellFastq_R1) as in_f1:
                                out1.write(in_f1.read())

                        if os.path.exists(CellFastq_R2):
                            with open(CellFastq_R2) as in_f2:
                                out2.write(in_f2.read())

          
            for cluster_name in clusters.keys():
                perfect_prefix = f"{Species}@{cluster_name}"
                command = f'spades.py --sc --pe1-1 {SpeciesStrainAssemDir}/{cluster_name}_R1.fastq --pe1-2 {SpeciesStrainAssemDir}/{cluster_name}_R2.fastq -o {SpeciesStrainAssemDir}/{perfect_prefix}_sc'
                try:
                    if env is not None:
                        result = subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                        success = (result.returncode == 0)
                    else:
                        return_code = subprocess.call(command, shell=True)
                        success = (return_code == 0)

                    if success:
                        print(f"cluster {perfect_prefix} assem successfully")
                        os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}_R1.fastq')
                        os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}_R2.fastq')
                        
                        import shutil
                        source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/scaffolds.fasta'
                        final_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}.fasta'
                        if not os.path.exists(source_fasta):
                            source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/contigs.fasta'
                        if os.path.exists(source_fasta):
                            shutil.copy(source_fasta, final_fasta)
                            print(f" -> Successfully extracted the final genome: {perfect_prefix}.fasta")
                    else:
                        print(f"cluster {perfect_prefix} failure")

                except Exception as e:
                    print(f"execute spades.py command failed: {e}")
                    success = False






class AllBin():
    def __init__(self,FastaDir, FastqDir, CellAnno, ResultDir):
        self.FastaDir = FastaDir
        self.FastqDir = FastqDir
        self.CellAnno = CellAnno
        self.ResultDir = ResultDir
        self.SpeciesStrainsInfo: Dict[str, str] = {} 
        self.BinClusterAnno = os.path.join(self.ResultDir, "BinClusterAnno.txt")

    strain_color_list = ["#33bbff", "#0050e6", "#009999", "#777777", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                         "#8c564b",
                         "#e377c2"]
    umap_plot_shape = ['D', 'o', 's', 'v', 'd', 'p', '*', 'H', '+', 'x']
    default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2',
                          u'#7f7f7f',
                          u'#bcbd22', u'#17becf']
    
    @timeit
    def AllBinPrepare(self, bcftools=None, snap_aligner=None, samtools=None, env=None, ReadsEnd='Single', exclude_clusters=None, min_cells=5, sort_tool='samtools'):

        FastaDir = self.FastaDir
        FastqDir = self.FastqDir
        CellAnno = self.CellAnno
        ResultDir = self.ResultDir

        ######  Organize all species cells into the directory format of the file required for SingleBinPrepare  ######
        # ResultDir：ClusterxxxInput ClusterxxxResult（It is required that the name of the Cluster in the CellAnno file be consistent with the fasta file name）
        try:
            if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
                raise ValueError('ReadsEnd must be Single or Pair!')
            if sort_tool not in ['samtools', 'snap']:
                raise ValueError("Currently, only 'samtools' or 'snap' are supported as sorting software.")
        except ValueError as e:
            print(e)
            sys.exit()

        CellAnno = pd.read_csv(CellAnno, sep='\t', header=0)
        
        if 'Cluster' not in CellAnno.columns:
            if 'SGB' in CellAnno.columns:
                CellAnno.rename(columns={'SGB': 'Cluster'}, inplace=True)
            else:
                # If neither of them exists, throw an exception directly to prevent the program from crashing halfway through
                raise ValueError("CellAnno Error: The file must contain either a 'Cluster' or 'SGB' column!")

        if exclude_clusters: 

            
            if isinstance(exclude_clusters, str):
                exclude_clusters = [exclude_clusters]

            exclude_list_str = [str(item) for item in exclude_clusters]
            CellAnno = CellAnno[~CellAnno['Cluster'].astype(str).isin(exclude_list_str)]
            
            print(f"\n>> [Data Preprocessing] Automatically filtered out classifications not required for analysis: {exclude_list_str}")


        if bcftools is None:
            bcftools = 'bcftools'

        if snap_aligner is None:
            snap_aligner = 'snap-aligner'
            
        if samtools is None:
            samtools = 'samtools'

        
        # Scan the Fasta file in advance to create folders and soft links only for valid SGB
        summary_data = []
        valid_clusters = []
        
        for cluster in CellAnno['Cluster'].unique():
            SpeciesName = str(cluster)
            sourceFasta = os.path.join(FastaDir, SpeciesName) + '.fasta'
            
            if not os.path.exists(sourceFasta):
                print(f"\n>> [Skipped] SGB: {SpeciesName} (Missing corresponding Fasta file)")
                summary_data.append([SpeciesName, "No Fasta", 0])
            else:
                valid_clusters.append(SpeciesName)
                
        # Filter CellAnno and remove the cell rows without Fasta to avoid doing unnecessary work later
        CellAnno = CellAnno[CellAnno['Cluster'].astype(str).isin(valid_clusters)]


        # Create a "valid" ClusterInput folder
        DirCreate = []
        for cluster in valid_clusters:
            dir1 = os.path.join(ResultDir, str(cluster)) + 'Input'
            dir2 = os.path.join(ResultDir, str(cluster)) + 'Result'
            DirCreate.append(dir1)
            DirCreate.append(dir2)

        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        # Only move the fastq files of the "effective" cells
        if ReadsEnd == 'Pair':
            for index, row in CellAnno.iterrows():
                Cluster = str(row['Cluster'])
                Cell = str(row['Cell'])

                sourceFastq1 = os.path.join(FastqDir, Cell) + '_R1.fastq'
                sourceFastq2 = os.path.join(FastqDir, Cell) + '_R2.fastq'
                targetFastq1 = os.path.join(ResultDir, str(Cluster) + 'Input', Cell) + '_R1.fastq'
                targetFastq2 = os.path.join(ResultDir, str(Cluster) + 'Input', Cell) + '_R2.fastq'

                if not os.path.exists(targetFastq1):
                    os.symlink(sourceFastq1, targetFastq1)
                if not os.path.exists(targetFastq2):
                    os.symlink(sourceFastq2, targetFastq2)

        elif ReadsEnd == 'Single':
            for index, row in CellAnno.iterrows():
                Cluster = str(row['Cluster'])
                Cell = str(row['Cell'])

                sourceFastq = os.path.join(FastqDir, Cell) + '.fastq'
                targetFastq = os.path.join(ResultDir, str(Cluster) + 'Input', Cell) + '.fastq'

                if not os.path.exists(targetFastq):
                    os.symlink(sourceFastq, targetFastq)

        
        self.BinClusterAnno = os.path.join(ResultDir,'BinClusterAnno.txt')
        with open(os.path.join(ResultDir,'BinClusterAnno.txt'),'w') as f:
            f.write('BinPrepareDir\tSpeciesName\tClusterNum\n')


        
        for SpeciesName in valid_clusters:
            sourceFasta = os.path.join(FastaDir, SpeciesName) + '.fasta'
            targetFasta = os.path.join(ResultDir, SpeciesName + 'Input', SpeciesName) + '.fasta'
            
            if not os.path.exists(targetFasta):
                os.symlink(sourceFasta, targetFasta)
                
            BinDir = os.path.join(ResultDir, SpeciesName) + 'Input'
            ResultDirSingle = os.path.join(ResultDir, SpeciesName) + 'Result'
            with open(os.path.join(ResultDir,'BinClusterAnno.txt'),'a') as f:
                f.write(ResultDirSingle + '\t' + SpeciesName +'\n')

            sbp = SingleBin(BinDir,ResultDirSingle,SpeciesName)
            status = sbp.SingleBinPrepare(bcftools=bcftools, 
                                          snap_aligner=snap_aligner, 
                                          samtools=samtools, 
                                          env=env, 
                                          ReadsEnd=ReadsEnd,
                                          min_cells=min_cells,
                                          sort_tool=sort_tool)
            
            passed_cells = max(0, getattr(sbp, 'AllCellNum', 0) - getattr(sbp, 'DropCellNum', 0))
            
            if status is False:
                summary_data.append([SpeciesName, "FAILED", passed_cells])
            else:
                summary_data.append([SpeciesName, "SUCCESS", passed_cells])

      
        summary_file = os.path.join(ResultDir, 'Strain_Analysis_Summary.csv')
        with open(summary_file, 'w') as f:
           
            f.write("SpeciesName,Status,Passed_Cells\n")
            
           
            for row in summary_data:
               
                f.write(f"{row[0]},{row[1]},{row[2]}\n")
            
        print(f"\nAll processing completed! Summary log available at: {summary_file}")


    def SingleBinSplit(self,PrepareDir, ClusterNum, Species):



        ClusterNum = int(ClusterNum)

        # PrepareBam = os.path.join(PrepareDir, 'ResultBam')
        PrepareVCF = os.path.join(PrepareDir, 'ResultVCF')
        PreparePic = os.path.join(PrepareDir, 'ResultPic')

        with open(PrepareVCF + '/linked.pkl', 'rb') as f:
            linked = pickle.load(f)

        with open(PrepareVCF + '/snp_pd.pkl', 'rb') as f:
            snp_pd = pickle.load(f)

        with open(PrepareVCF + '/embedding.pkl', 'rb') as f:
            embedding = pickle.load(f)

        cutree_contigs = cut_tree(linked, n_clusters=ClusterNum)
        umap_compare = pd.DataFrame({'Barcode': snp_pd.columns[:-4]})
        umap_compare['X'] = embedding[:, 0]
        umap_compare['Y'] = embedding[:, 1]

        cluster_hier = []
        for i in range(len(cutree_contigs)):
            cluster_hier.append(cutree_contigs[i][0])

        umap_compare['cluster_hier'] = cluster_hier

        for j in range(ClusterNum):
            umap_compare_bin = umap_compare.loc[umap_compare['cluster_hier'] == j]
            plt.scatter(umap_compare_bin['X'], umap_compare_bin['Y'], s=10, c=AllBin.strain_color_list[j],
                        marker=AllBin.umap_plot_shape[j])

        for j in range(ClusterNum):
            # Calculate the center position of each cluster
            cluster_points = umap_compare.loc[umap_compare['cluster_hier'] == j]
            x_center = cluster_points['X'].mean()
            y_center = cluster_points['Y'].mean()
            plt.text(x_center, y_center, str(j), fontsize=12, fontweight='bold', color='black',
                     bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        plt.xticks([])
        plt.yticks([])
        plt.ylabel('UMAP dimention 2')
        plt.xlabel('UMAP dimention 1')
        plt.savefig(os.path.join(PreparePic, Species) + '_sp_call_SNP_umap_final.svg', format='svg',
                    bbox_inches='tight')
        plt.clf()

        ###################################
        ###get strain share SNP genotype###
        ###################################

        # file = os.path.join(PrepareVCF,Species) + '_calls_snps_qual30.vcf_snp_filtered_bc_dropped.csv'

        for i in range(ClusterNum):
            index_of_this_cluster = []
            for j in range(len(cutree_contigs)):
                if cutree_contigs[j] == i:
                    index_of_this_cluster += [j]
                snp_pd['hierarchy_cluster_' + str(i) + 'pos'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == 1).sum(axis=1)
                snp_pd['hierarchy_cluster_' + str(i) + 'neg'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == -1).sum(axis=1)

        strain_X_genotype = [[] for i in range(ClusterNum)]  # Create an empty two-dimensional array of length ClusterNum
        for i in range(len(snp_pd)):
            for j in range(ClusterNum):
                # generate genotype of each strain
                # at each SNP, if both have less than 2 cells, assign it 0, means not sure
                # if the max genotype is less than 90% total, assign 0
                # otherwise assign the major genotype
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] <= 1 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][
                            i] <= 1:
                    strain_X_genotype[j].append(0)
                    continue
                elif snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] != 0 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i] != 0:
                    if max(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                           snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]) / float(
                        min(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                            snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i])) < 9:
                        strain_X_genotype[j].append(0)
                        continue
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] > snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]:
                    strain_X_genotype[j].append(1)
                else:
                    strain_X_genotype[j].append(-1)

        for j in range(ClusterNum):
            snp_pd['strain_' + str(j) + '_genotype'] = strain_X_genotype[j]

        rows_to_drop = []
        for i in range(len(snp_pd)):
            genotype_SNP_temp = []
            for j in range(ClusterNum):
                genotype_SNP_temp.append(snp_pd['strain_' + str(j) + '_genotype'][i])
            if not (1 in genotype_SNP_temp and -1 in genotype_SNP_temp):
                rows_to_drop.append(i)

        print("SNPs before drop is {}, SNPs to drop is {}.".format(len(snp_pd), len(rows_to_drop)))
        snp_pd = snp_pd.drop(rows_to_drop, axis=0)
        snp_pd = snp_pd.reset_index(drop=True)

        # this is to generate strain list
        list_of_ratio = [[] for i in range(ClusterNum)]
        barcode_list_strain = [[] for i in range(ClusterNum + 1)]

        for i in range(snp_pd.shape[1]):
            column_name_temp = snp_pd.columns[i]
            if column_name_temp.startswith('strain') or column_name_temp.startswith('hierarchy_cluster'):
                continue

            if 'list' not in column_name_temp:
                y = column_name_temp.split('/')[-1]
            else:
                continue

            # if not column_name_temp in snp_pd.columns:
            # continue
            temp_df_check = snp_pd.loc[snp_pd[column_name_temp] != 0]
            for j in range(ClusterNum):
                temp_df_check_strain = temp_df_check.loc[temp_df_check['strain_' + str(j) + '_genotype'] != 0]
                total_count_same = sum(temp_df_check_strain[column_name_temp] == \
                                       temp_df_check_strain['strain_' + str(j) + '_genotype'])
                list_of_ratio[j].append(total_count_same / float(len(temp_df_check_strain[column_name_temp])))

            list_of_ratio_check = [list_of_ratio[j][-1] for j in range(ClusterNum)]
            print(list_of_ratio_check)

            if max(list_of_ratio_check) < 0.95:
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            elif sum(np.array(list_of_ratio_check) > 0.95) > 1:
                print('Not sure,more then two strains match', column_name_temp)
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            else:
                # barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(column_name_temp)
                barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(y)

        with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'w') as output:
            output.write('Cluster\tCell\n')

        for cluster in range(len(barcode_list_strain)):
            # print('look here is barcode_list_strain:', i)
            # print(barcode_list_strain[i])
            if cluster == (len(barcode_list_strain) - 1):
                for cell in barcode_list_strain[cluster]:
                    with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'a') as output:
                        output.write('Mist'+str(cluster) + '\t' + str(cell) + '\n')
            else:
                for cell in barcode_list_strain[cluster]:
                    with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'a') as output:
                        output.write(str(cluster) + '\t' + str(cell) + '\n')

        self.SpeciesStrainsInfo[Species] = str(os.path.join(PrepareDir, 'StrainCells.txt'))


    @timeit
    def AllBinSplit(self):
        # need to input the ClusterxxxResult path of the strain you want to classify and determine how many clusters the box is divided into.
        # BinClusterAnno.txt (SpeciesName must be consistent with the fasta file name under the box)
        # BinPrepareDir    ClusterNum    SpeciesName
        # /data_alluser/.. /SGB4933Result 4(the number of clusters cannot exceed 10) SGB4933

        # Read the BinClusterAnno file and call SingleBinSplit() to separate strains box by box


        BinClusterAnnoDir = self.BinClusterAnno
        BinClusterAnno = pd.read_csv(BinClusterAnnoDir, sep='\t', header=0)
        for index, row in BinClusterAnno.iterrows():
            BinPrepareDir = row['BinPrepareDir']
            SpeciesName = row['SpeciesName']
            
            
            raw_cluster_num = row['ClusterNum']
            
            
            if pd.isna(raw_cluster_num) or str(raw_cluster_num).strip() == '':
                print(f"  [Warning] Skipping {SpeciesName}: ClusterNum not specified (cluster count is empty)")
                continue
                
           
            try:
                ClusterNum = int(raw_cluster_num)
            except ValueError:
                print(f"  [Error] Skipping {SpeciesName}: Invalid ClusterNum format (received '{raw_cluster_num}'); must be an integer!")
                continue
            # =======================================================

            print(f"\n>> Starting strain splitting for SGB: {SpeciesName} (splitting into {ClusterNum} strains)...")

            #SingleBinSplit(BinPrepareDir, ClusterNum, SpeciesName)
            self.SingleBinSplit(BinPrepareDir,ClusterNum, SpeciesName)



    @timeit
    def AllBinfStrainAssem(self, StrainAssemDir="", env=None, ReadsEnd='Single'):

        self.ReadsEnd = ReadsEnd

        # If StrainAssemDir is an empty string, then place the result of each assembled strain in the automatically created species /StrainAssem/ folder.

        if StrainAssemDir == "":
            for Species,StrainCellsFile in self.SpeciesStrainsInfo.items():
                SpeciesStrainAssemDir = os.path.join(os.path.dirname(StrainCellsFile), 'StrainAssem')
                if not os.path.exists(SpeciesStrainAssemDir):
                    os.makedirs(SpeciesStrainAssemDir,exist_ok=True)

                clusters = defaultdict(list)
                with open(StrainCellsFile) as f:
                    for line in f:
                        if line.strip():
                            
                            parts = line.split()
                            # Defensive programming: Ensure that there are at least two columns of data in this row
                            if len(parts) >= 2:
                                c, cell = parts[0], parts[1]
                                
                                # Kick out the header: If the first column of this row is called 'Cluster' or 'cluster', skip it directly
                                if c.lower() == "cluster":
                                    continue
                                    
                                # Blur to block noise: Do not use any names that contain 'Mist' (such as Mist, Mist2, Mist_1)!
                                if "Mist" not in c:
                                    clusters[c].append(cell)

                if self.ReadsEnd == 'Single':
                    for cl, cells in clusters.items():
                        with open(f'{SpeciesStrainAssemDir}/{cl}.fastq', 'w') as out:
                            for cell in cells:
                                CellFastq = os.path.join(self.FastqDir, cell) + '.fastq'
                                if os.path.exists(CellFastq):
                                    with open(CellFastq) as in_f:
                                        out.write(in_f.read())

                    # Assemble each bacterial type.fastq file
                    for cluster_name in clusters.keys():
                        perfect_prefix = f"{Species}@{cluster_name}"
                        command = f'spades.py --sc --pe1-s {SpeciesStrainAssemDir}/{cluster_name}.fastq -o {SpeciesStrainAssemDir}/{perfect_prefix}_sc'
                        try:
                            if env is not None:
                                result = subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                                success = (result.returncode == 0)
                            else:
                                return_code = subprocess.call(command, shell=True)
                                success = (return_code == 0)

                            if success:
                                print(f"cluster {perfect_prefix} assem successfully")
                                os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}.fastq')
                                
                                import shutil
                                source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/scaffolds.fasta'
                                final_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}.fasta'
                                
                                if not os.path.exists(source_fasta):
                                    source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/contigs.fasta'
                                
                                if os.path.exists(source_fasta):
                                    shutil.copy(source_fasta, final_fasta)
                                    print(f"  -> The final genome has been extracted: {perfect_prefix}.fasta")
                            else:
                                print(f"cluster {perfect_prefix} failure")

                        except Exception as e:
                            print(f"execute spades.py command failed: {e}")
                            success = False


                elif self.ReadsEnd == 'Pair':
                    for cl, cells in clusters.items():
                        with open(f'{SpeciesStrainAssemDir}/{cl}_R1.fastq', 'w') as out1, open(f'{SpeciesStrainAssemDir}/{cl}_R2.fastq',
                                                                                        'w') as out2:
                            for cell in cells:
                                CellFastq_R1 = os.path.join(self.FastqDir, cell) + '_R1.fastq'
                                CellFastq_R2 = os.path.join(self.FastqDir, cell) + '_R2.fastq'
                                if os.path.exists(CellFastq_R1):
                                    with open(CellFastq_R1) as in_f1:
                                        out1.write(in_f1.read())

                                if os.path.exists(CellFastq_R2):
                                    with open(CellFastq_R2) as in_f2:
                                        out2.write(in_f2.read())

                    # Assemble each strain type _R1.fastq and strain type _R2.fastq file
                    for cluster_name in clusters.keys():
                        perfect_prefix = f"{Species}@{cluster_name}"
                        command = f'spades.py --sc --pe1-1 {SpeciesStrainAssemDir}/{cluster_name}_R1.fastq --pe1-2 {SpeciesStrainAssemDir}/{cluster_name}_R2.fastq -o {SpeciesStrainAssemDir}/{perfect_prefix}_sc'
                        try:
                            if env is not None:
                                result = subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                                success = (result.returncode == 0)
                            else:
                                return_code = subprocess.call(command, shell=True)
                                success = (return_code == 0)

                            if success:
                                print(f"cluster {perfect_prefix} assem successfully")
                                os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}_R1.fastq')
                                os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}_R2.fastq')
                                
                                import shutil
                                source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/scaffolds.fasta'
                                final_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}.fasta'
                                
                                if not os.path.exists(source_fasta):
                                    source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/contigs.fasta'
                                
                                if os.path.exists(source_fasta):
                                    shutil.copy(source_fasta, final_fasta)
                                    print(f"  -> The final genome has been extracted: {perfect_prefix}.fasta")
                            else:
                                print(f"cluster {perfect_prefix} failure")

                        except Exception as e:
                            print(f"execute spades.py command failed: {e}")
                            success = False

       

        else:
            if not os.path.exists(StrainAssemDir):
                os.makedirs(StrainAssemDir, exist_ok=True)


            for Species,StrainCellsFile in self.SpeciesStrainsInfo.items():
                SpeciesStrainAssemDir = os.path.join(StrainAssemDir, Species+'_StrainAssem')
                if not os.path.exists(SpeciesStrainAssemDir):
                    os.makedirs(SpeciesStrainAssemDir,exist_ok=True)

                clusters = defaultdict(list)
                with open(StrainCellsFile) as f:
                    for line in f:
                        if line.strip():
                            
                            parts = line.split()
                            
                            if len(parts) >= 2:
                                c, cell = parts[0], parts[1]
                                
                              
                                if c.lower() == "cluster":
                                    continue
                                    
                              
                                if "Mist" not in c:
                                    clusters[c].append(cell)

                if self.ReadsEnd == 'Single':
                    for cl, cells in clusters.items():
                        with open(f'{SpeciesStrainAssemDir}/{cl}.fastq', 'w') as out:
                            for cell in cells:
                                CellFastq = os.path.join(self.FastqDir, cell) + '.fastq'
                                if os.path.exists(CellFastq):
                                    with open(CellFastq) as in_f:
                                        out.write(in_f.read())

                   
                    for cluster_name in clusters.keys():
                        perfect_prefix = f"{Species}@{cluster_name}"
                        command = f'spades.py --sc --pe1-s {SpeciesStrainAssemDir}/{cluster_name}.fastq -o {SpeciesStrainAssemDir}/{perfect_prefix}_sc'
                        try:
                            if env is not None:
                                result = subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                                success = (result.returncode == 0)
                            else:
                                return_code = subprocess.call(command, shell=True)
                                success = (return_code == 0)

                            if success:
                                print(f"cluster {perfect_prefix} assem successfully")
                                os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}.fastq')
                                
                                import shutil
                                source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/scaffolds.fasta'
                                final_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}.fasta'
                                
                                if not os.path.exists(source_fasta):
                                    source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/contigs.fasta'
                                
                                if os.path.exists(source_fasta):
                                    shutil.copy(source_fasta, final_fasta)
                                    print(f"  -> The final genome has been extracted: {perfect_prefix}.fasta")
                            else:
                                print(f"cluster {perfect_prefix} failure")

                        except Exception as e:
                            print(f"execute spades.py command failed: {e}")
                            success = False


                elif self.ReadsEnd == 'Pair':
                    for cl, cells in clusters.items():
                        with open(f'{SpeciesStrainAssemDir}/{cl}_R1.fastq', 'w') as out1, open(f'{SpeciesStrainAssemDir}/{cl}_R2.fastq','w') as out2:
                            for cell in cells:
                                CellFastq_R1 = os.path.join(self.FastqDir, cell) + '_R1.fastq'
                                CellFastq_R2 = os.path.join(self.FastqDir, cell) + '_R2.fastq'
                                if os.path.exists(CellFastq_R1):
                                    with open(CellFastq_R1) as in_f1:
                                        out1.write(in_f1.read())

                                if os.path.exists(CellFastq_R2):
                                    with open(CellFastq_R2) as in_f2:
                                        out2.write(in_f2.read())

                    
                    for cluster_name in clusters.keys():
                        perfect_prefix = f"{Species}@{cluster_name}"
                        command = f'spades.py --sc --pe1-1 {SpeciesStrainAssemDir}/{cluster_name}_R1.fastq --pe1-2 {SpeciesStrainAssemDir}/{cluster_name}_R2.fastq -o {SpeciesStrainAssemDir}/{perfect_prefix}_sc'
                        try:
                            if env is not None:
                                result = subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                                success = (result.returncode == 0)
                            else:
                                return_code = subprocess.call(command, shell=True)
                                success = (return_code == 0)

                            if success:
                                print(f"cluster {perfect_prefix} assem successfully")
                                os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}_R1.fastq')
                                os.remove(f'{SpeciesStrainAssemDir}/{cluster_name}_R2.fastq')
                                
                                import shutil
                                source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/scaffolds.fasta'
                                final_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}.fasta'
                                
                                if not os.path.exists(source_fasta):
                                    source_fasta = f'{SpeciesStrainAssemDir}/{perfect_prefix}_sc/contigs.fasta'
                                
                                if os.path.exists(source_fasta):
                                    shutil.copy(source_fasta, final_fasta)
                                    print(f"  -> The final genome has been extracted: {perfect_prefix}.fasta")
                            else:
                                print(f"cluster {perfect_prefix} failure")

                        except Exception as e:
                            print(f"execute spades.py command failed: {e}")
                            success = False















