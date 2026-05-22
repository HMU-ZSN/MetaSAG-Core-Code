import umap
import subprocess
from sklearn import preprocessing
import time
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram,linkage,cut_tree,fcluster
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
import shutil

'''
script comment
This document requires input:
(1) the fasta file
(2) fasta's annotation file (species, the gate where the species is located)
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

class CellHGT():
    def __init__(self,FastaDir,outputDir):
        self.FastaDir = FastaDir
        self.outputDir = outputDir

    def merge_intervals(intervals):
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
        merged = []

        for higher in sorted_by_lower_bound:
            if not merged:
                merged.append(higher)
            else:
                lower = merged[-1]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1] = (lower[0], upper_bound)
                else:
                    merged.append(higher)
        return merged

    # There are two scenarios: one is that the input Anno file is at the species level or above, and the other is the classification of species and their corresponding strain levels obtained through SNP
    # These two functions not only return some files but also the HGTNum data frame for the subsequent drawing of the HGT Tree
    @timeit
    def SpeciesHGT(self):
        fastaDir = self.FastaDir
        HGTTemp = os.path.join(self.outputDir,'SpeciesHGT')

        # create directory
        DB = os.path.join(HGTTemp, 'DB/')
        Blast = os.path.join(HGTTemp, 'Blast/')


        DirCreate = [DB, Blast]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        ######Read the list of Fasta files######

        file_list_of_assemblies = os.listdir(fastaDir)
        for file_temp in file_list_of_assemblies:
            if '.fasta' not in file_temp:
                file_list_of_assemblies.remove(file_temp)

        ######Create a reference library######
        for file_temp in file_list_of_assemblies:
            fastafile = os.path.join(fastaDir, file_temp[:-6]) + '.fasta'
            dbfile = os.path.join(DB, file_temp[:-6]) + '_blastdb'
            subprocess.call('makeblastdb -dbtype "nucl" -in ' + fastafile + ' -out ' + dbfile, shell=True)

        ######blast between each pair of fasta files using the reference library######
        for query_temp in file_list_of_assemblies:
            for subject_temp in file_list_of_assemblies:
                queryfasta = os.path.join(fastaDir, query_temp)
                subjectdb = os.path.join(DB, subject_temp[:-6]) + '_blastdb'
                blastresult = os.path.join(Blast, query_temp[:-6] + '_to_' + subject_temp[:-6]) + '_blast_result'
                subprocess.call(
                    'blastn -query ' + queryfasta + ' -db ' + subjectdb + ' -outfmt 6 -perc_identity 99.98 -out ' + blastresult,
                    shell=True)

        ######Put the blastn results into the data frame: fasta number of qcontig,qcontig number, fasta number of scontig,scontig number, length, mismatch######
        file_list_of_assemblies = os.listdir(Blast)
        a_temp_list = []
        blastn_HGT_result = []

        for i in range(len(file_list_of_assemblies)):
            file_temp = file_list_of_assemblies[i]
            if '_blast_result' not in file_temp:
                continue

            # elif file_temp.split('to')[0][:-1] != file_temp.split('to')[1][1:1 + len(file_temp.split('to')[0][:-1])]:
            elif file_temp.split('to')[0][:-1] != file_temp.split('to')[1][1:-13]:
                a_temp_list += [file_temp]

        file_list_of_assemblies = a_temp_list
        for file_temp in file_list_of_assemblies:
            try:
                temp_blastn_result = pd.read_csv(os.path.join(Blast, file_temp), header=None, sep='\t')
            except:
                continue
            temp_blastn_result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                                          'qend',
                                          'sstart', 'send', 'evalue', 'bitscore']
            temp_blastn_result['query_cluster'] = [file_temp.split('to')[0][:-1]] * len(temp_blastn_result)
            temp_blastn_result['subject_cluster'] = [file_temp.split('to')[1][1:-len('_blast_result')]] * len(
                temp_blastn_result)
            temp_blastn_result = temp_blastn_result.loc[temp_blastn_result['length'] > 500]
            if len(blastn_HGT_result) > 0:
                blastn_HGT_result = pd.concat([blastn_HGT_result, temp_blastn_result])
            else:
                blastn_HGT_result = temp_blastn_result

        blastn_HTG_result = blastn_HGT_result.sort_values(by=['query_cluster', 'subject_cluster', ])
        blastn_HGT_result = blastn_HGT_result.reset_index(drop=True)
        cols = list(blastn_HGT_result.columns.values)

        blastn_HGT_result = blastn_HGT_result[cols[-2:-1] + cols[:1] + cols[-1:] + cols[1:2] + cols[2:-2]]
        blastn_HGT_result.to_csv(os.path.join(Blast, 'blastn_HGT_result.csv'), index=False)

        ####Possible HGT strain pairs were selected based on length (5Kb) and the similarity percentage of sequence alignment consistency (99.98)
        length_threshold = 5000
        match_bp_ratio = 99.98
        # length_threshold=1000
        # match_bp_ratio=80
        # blastn_HGT_result = pd.read_csv(HGT_strain_Blastoutdirt + 'blastn_HGT_result.csv')
        blastn_result_sc = blastn_HGT_result.loc[blastn_HGT_result['length'] > length_threshold].loc[
            blastn_HGT_result.loc[blastn_HGT_result['length'] > length_threshold]['pident'] > match_bp_ratio]
        blastn_result_sc = blastn_result_sc.reset_index(drop=True)
        blastn_result_sc.to_csv(os.path.join(Blast, 'blastn_HGT_result_trimmed.csv'), index=False)

        ####
        query_cluster_list = []
        subject_cluster_list = []
        total_length_of_HGT = []
        # blastn_result_sc = pd.read_csv(Blast + 'blastn_HGT_result_trimmed.csv')
        assembly_all = list(
            set(list(set(blastn_result_sc['query_cluster'])) + list(set(blastn_result_sc['subject_cluster']))))

        for query_temp in assembly_all:
            df_temp_query = blastn_result_sc.loc[blastn_result_sc['query_cluster'] == query_temp]
            for subject_temp in list(set(df_temp_query['subject_cluster'])):
                df_temp_query_subject = df_temp_query.loc[df_temp_query['subject_cluster'] == subject_temp]
                df_temp_query_subject = df_temp_query_subject.reset_index(drop=True)
                query_cluster_list += [query_temp]
                subject_cluster_list += [subject_temp]
                all_contig_list_temp = list(set(df_temp_query_subject['qseqid']))
                total_HGT_count = 0
                for i in range(len(all_contig_list_temp)):
                    contig_to_check_temp = all_contig_list_temp[i]
                    df_temp_query_subject_contig = df_temp_query_subject.loc[
                        df_temp_query_subject['qseqid'] == contig_to_check_temp]
                    df_temp_query_subject_contig = df_temp_query_subject_contig.reset_index(drop=True)
                    all_intervals_temp = []
                    for j in range(len(df_temp_query_subject_contig)):
                        start_temp = min(df_temp_query_subject_contig['qstart'][j],
                                         df_temp_query_subject_contig['qend'][j])
                        end_temp = max(df_temp_query_subject_contig['qstart'][j],
                                       df_temp_query_subject_contig['qend'][j])
                        all_intervals_temp += [(start_temp, end_temp)]
                    merged_intervals_temp = CellHGT.merge_intervals(all_intervals_temp)
                    total_length_intervals = 0
                    for each_interval_temp in merged_intervals_temp:
                        total_length_intervals += each_interval_temp[1] - each_interval_temp[0] + 1
                    total_HGT_count += total_length_intervals

                total_length_of_HGT += [total_HGT_count]

        blastn_result_sc_bp_count = pd.DataFrame({'query_cluster_list': query_cluster_list})
        blastn_result_sc_bp_count['subject_cluster_list'] = subject_cluster_list
        blastn_result_sc_bp_count['total_length_of_HGT'] = total_length_of_HGT
        blastn_result_sc_bp_count.to_csv(os.path.join(Blast, 'blastn_result_sc_bp_count.csv'), index=False)

        ######getDiffSpeciesHGT.py######
        # blastn_result_sc = pd.read_csv(HGT_strain_Blastoutdirt + 'blastn_HGT_result_trimmed.csv', sep=',')
        cols = blastn_result_sc.columns.values
        blastn_result_diffbin = pd.DataFrame(columns=cols)
        repeat_pair = []
        for index, row in blastn_result_sc.iterrows():
            bin1 = row['query_cluster']
            bin2 = row['subject_cluster']
            pair_temp = [str(row['query_cluster']) + str(row['qseqid']),
                         str(row['subject_cluster']) + str(row['sseqid'])]
            pair_temp.sort()
            pair = "".join(pair_temp)

            if bin1 != bin2 and pair not in repeat_pair:
                blastn_result_diffbin = blastn_result_diffbin._append(blastn_result_sc.loc[index])
                repeat_pair += [pair]

        blastn_result_diffbin.to_csv(os.path.join(Blast, 'blastn_result_diffbin.txt'), sep='\t', index=False)

        columns = ['Species1', 'Species2']
        HGTNum = pd.DataFrame(columns=columns)

        for index, row in blastn_result_diffbin.iterrows():
            bin = row['query_cluster']
            start = row['qstart'] - 1
            end = row['qend']

            fasta = {}
            NODE = ''

            with open(os.path.join(fastaDir, str(bin)) + '.fasta') as input:
                while True:
                    line = input.readline()
                    if len(line) == 0:
                        break
                    line = line.replace('\n', '')
                    if line.startswith('>'):
                        NODE = line.replace('>', '')
                        fasta[NODE] = ''
                    else:
                        fasta[NODE] = ''.join((fasta[NODE], str(line)))

            NODE = row['qseqid']
            read = fasta[NODE][start:end]
            with open(os.path.join(HGTTemp, 'HGT.fasta'), 'a') as output:
                qc = row['query_cluster']
                qid = row['qseqid'].split('_')[1]
                sc = row['subject_cluster']
                sid = row['sseqid'].split('_')[1]

                output.write('>PAIR_' + str(qc) + '_' + str(qid) + 'TO' + str(sc) + '_' + str(sid) + '\n')
                output.write(read + '\n')

                species1 = qc
                species2 = sc
                HGTNumRow = [species1, species2]
                HGTNum.loc[len(HGTNum)] = HGTNumRow

        self.HGTfasta = os.path.join(HGTTemp, 'HGT.fasta')
        HGTNum.to_csv(os.path.join(HGTTemp, 'HGTNum.txt'), sep='\t', index=False)
        self.HGTNum = HGTNum
        return HGTNum

    @timeit
    def StrainHGT(self):

        fastaDir = self.FastaDir
        HGTTemp = os.path.join(self.outputDir, 'StrainHGT')

        # create directory
        DB = os.path.join(HGTTemp, 'DB/')
        Blast = os.path.join(HGTTemp, 'Blast/')

        DirCreate = [DB, Blast]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        ######Read the list of Fasta files######

        file_list_of_assemblies = os.listdir(fastaDir)
        for file_temp in file_list_of_assemblies:
            if '.fasta' not in file_temp:
                file_list_of_assemblies.remove(file_temp)

        ######Create a reference library######
        for file_temp in file_list_of_assemblies:
            fastafile = os.path.join(fastaDir, file_temp[:-6]) + '.fasta'
            dbfile = os.path.join(DB, file_temp[:-6]) + '_blastdb'
            subprocess.call('makeblastdb -dbtype "nucl" -in ' + fastafile + ' -out ' + dbfile, shell=True)

        ######blast between each pair of fasta files using the reference library######
        for query_temp in file_list_of_assemblies:
            for subject_temp in file_list_of_assemblies:
                queryfasta = os.path.join(fastaDir, query_temp)
                subjectdb = os.path.join(DB, subject_temp[:-6]) + '_blastdb'
                blastresult = os.path.join(Blast, query_temp[:-6] + '_to_' + subject_temp[:-6]) + '_blast_result'
                subprocess.call(
                    'blastn -query ' + queryfasta + ' -db ' + subjectdb + ' -outfmt 6 -perc_identity 99.98 -out ' + blastresult,
                    shell=True)

        ######Put the blastn results into the data frame: fasta number of qcontig,qcontig number, fasta number of scontig,scontig number, length, mismatch######
        file_list_of_assemblies = os.listdir(Blast)
        a_temp_list = []
        blastn_HGT_result = []

        for i in range(len(file_list_of_assemblies)):
            file_temp = file_list_of_assemblies[i]
            if '_blast_result' not in file_temp:
                continue

            # elif file_temp.split('to')[0][:-1] != file_temp.split('to')[1][1:1 + len(file_temp.split('to')[0][:-1])]:
            elif file_temp.split('to')[0][:-1] != file_temp.split('to')[1][1:-13]:
                a_temp_list += [file_temp]

        file_list_of_assemblies = a_temp_list
        for file_temp in file_list_of_assemblies:
            try:
                temp_blastn_result = pd.read_csv(os.path.join(Blast, file_temp), header=None, sep='\t')
            except:
                continue
            temp_blastn_result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                                          'qend',
                                          'sstart', 'send', 'evalue', 'bitscore']
            temp_blastn_result['query_cluster'] = [file_temp.split('to')[0][:-1]] * len(temp_blastn_result)
            temp_blastn_result['subject_cluster'] = [file_temp.split('to')[1][1:-len('_blast_result')]] * len(
                temp_blastn_result)
            temp_blastn_result = temp_blastn_result.loc[temp_blastn_result['length'] > 500]
            if len(blastn_HGT_result) > 0:
                blastn_HGT_result = pd.concat([blastn_HGT_result, temp_blastn_result])
            else:
                blastn_HGT_result = temp_blastn_result

        blastn_HTG_result = blastn_HGT_result.sort_values(by=['query_cluster', 'subject_cluster', ])
        blastn_HGT_result = blastn_HGT_result.reset_index(drop=True)
        cols = list(blastn_HGT_result.columns.values)

        blastn_HGT_result = blastn_HGT_result[cols[-2:-1] + cols[:1] + cols[-1:] + cols[1:2] + cols[2:-2]]
        blastn_HGT_result.to_csv(os.path.join(Blast, 'blastn_HGT_result.csv'), index=False)

        ####Possible HGT strain pairs were selected based on length (5Kb) and the similarity percentage of sequence alignment consistency (99.98)
        length_threshold = 5000
        match_bp_ratio = 99.98
        # length_threshold=1000
        # match_bp_ratio=80
        # blastn_HGT_result = pd.read_csv(HGT_strain_Blastoutdirt + 'blastn_HGT_result.csv')
        blastn_result_sc = blastn_HGT_result.loc[blastn_HGT_result['length'] > length_threshold].loc[
            blastn_HGT_result.loc[blastn_HGT_result['length'] > length_threshold]['pident'] > match_bp_ratio]
        blastn_result_sc = blastn_result_sc.reset_index(drop=True)
        blastn_result_sc.to_csv(os.path.join(Blast, 'blastn_HGT_result_trimmed.csv'), index=False)

        ####
        query_cluster_list = []
        subject_cluster_list = []
        total_length_of_HGT = []
        # blastn_result_sc = pd.read_csv(Blast + 'blastn_HGT_result_trimmed.csv')
        assembly_all = list(
            set(list(set(blastn_result_sc['query_cluster'])) + list(set(blastn_result_sc['subject_cluster']))))

        for query_temp in assembly_all:
            df_temp_query = blastn_result_sc.loc[blastn_result_sc['query_cluster'] == query_temp]
            for subject_temp in list(set(df_temp_query['subject_cluster'])):
                df_temp_query_subject = df_temp_query.loc[df_temp_query['subject_cluster'] == subject_temp]
                df_temp_query_subject = df_temp_query_subject.reset_index(drop=True)
                query_cluster_list += [query_temp]
                subject_cluster_list += [subject_temp]
                all_contig_list_temp = list(set(df_temp_query_subject['qseqid']))
                total_HGT_count = 0
                for i in range(len(all_contig_list_temp)):
                    contig_to_check_temp = all_contig_list_temp[i]
                    df_temp_query_subject_contig = df_temp_query_subject.loc[
                        df_temp_query_subject['qseqid'] == contig_to_check_temp]
                    df_temp_query_subject_contig = df_temp_query_subject_contig.reset_index(drop=True)
                    all_intervals_temp = []
                    for j in range(len(df_temp_query_subject_contig)):
                        start_temp = min(df_temp_query_subject_contig['qstart'][j],
                                         df_temp_query_subject_contig['qend'][j])
                        end_temp = max(df_temp_query_subject_contig['qstart'][j],
                                       df_temp_query_subject_contig['qend'][j])
                        all_intervals_temp += [(start_temp, end_temp)]
                    merged_intervals_temp = CellHGT.merge_intervals(all_intervals_temp)
                    total_length_intervals = 0
                    for each_interval_temp in merged_intervals_temp:
                        total_length_intervals += each_interval_temp[1] - each_interval_temp[0] + 1
                    total_HGT_count += total_length_intervals

                total_length_of_HGT += [total_HGT_count]

        blastn_result_sc_bp_count = pd.DataFrame({'query_cluster_list': query_cluster_list})
        blastn_result_sc_bp_count['subject_cluster_list'] = subject_cluster_list
        blastn_result_sc_bp_count['total_length_of_HGT'] = total_length_of_HGT
        blastn_result_sc_bp_count.to_csv(os.path.join(Blast, 'blastn_result_sc_bp_count.csv'), index=False)

        ######getDiffSpeciesHGT.py######
        # Special note: In StrainHGT.py, since the genome is required to be at the strain level, the input fasta genome file should be in the format of species name @ strain number
        # blastn_result_sc = pd.read_csv(HGT_strain_Blastoutdirt + 'blastn_HGT_result_trimmed.csv', sep=',')
        cols = blastn_result_sc.columns.values
        blastn_result_diffbin = pd.DataFrame(columns=cols)
        repeat_pair = []
        for index, row in blastn_result_sc.iterrows():
            bin1 = row['query_cluster']
            bin2 = row['subject_cluster']
            pair_temp = [str(row['query_cluster']) + str(row['qseqid']),
                         str(row['subject_cluster']) + str(row['sseqid'])]
            pair_temp.sort()
            pair = "".join(pair_temp)

            if bin1 != bin2 and pair not in repeat_pair:
                blastn_result_diffbin = blastn_result_diffbin._append(blastn_result_sc.loc[index])
                repeat_pair += [pair]

        blastn_result_diffbin.to_csv(os.path.join(Blast, 'blastn_result_diffbin.txt'), sep='\t', index=False)

        columns = ['Species1', 'Strain1', 'Species2', 'Strain2']
        HGTNum = pd.DataFrame(columns=columns)

        for index, row in blastn_result_diffbin.iterrows():
            bin = row['query_cluster']
            start = row['qstart'] - 1
            end = row['qend']

            fasta = {}
            NODE = ''

            with open(os.path.join(fastaDir, str(bin)) + '.fasta') as input:
                while True:
                    line = input.readline()
                    if len(line) == 0:
                        break
                    line = line.replace('\n', '')
                    if line.startswith('>'):
                        NODE = line.replace('>', '')
                        fasta[NODE] = ''
                    else:
                        fasta[NODE] = ''.join((fasta[NODE], str(line)))

            NODE = row['qseqid']
            read = fasta[NODE][start:end]
            with open(os.path.join(HGTTemp, 'HGT.fasta'), 'a') as output:
                qc = row['query_cluster']
                qid = row['qseqid'].split('_')[1]
                sc = row['subject_cluster']
                sid = row['sseqid'].split('_')[1]

                output.write('>PAIR_' + str(qc) + '_' + str(qid) + 'TO' + str(sc) + '_' + str(sid) + '\n')
                output.write(read + '\n')

                species1 = qc.split('@')[0]
                strain1 = qc.split('@')[1]
                species2 = sc.split('@')[0]
                strain2 = sc.split('@')[1]

                HGTNumRow = [species1, strain1, species2, strain2]
                HGTNum.loc[len(HGTNum)] = HGTNumRow

        self.HGTfasta = os.path.join(HGTTemp, 'HGT.fasta')
        HGTNum.to_csv(os.path.join(HGTTemp, 'HGTNum.txt'), sep='\t', index=False)
        self.HGTTemp = HGTTemp
        return HGTNum

    

    PHY_COLOR = {"Verrucomicrobiota": "#660000", "Thermoproteota": "#783f04", "Thermoplasmatota": "#7f6000",
                 "Spirochaetota": "#274e13", "Proteobacteria": "#0c343d", "Planctomycetota": "#073763",
                 "Patescibacteria": "#20124d", "Myxococcota": "#4c1130", "Marinisomatota": "#cc0000",
                 "Halobacteriota": "#e69138", "Gemmatimonadota": "#f1c232", "Firmicutes_C": "#6aa84f",
                 "Firmicutes_B": "#45818e", "Firmicutes_A": "#3d85c6", "Firmicutes": "#674ea7",
                 "Desulfobacterota_I": "#a64d79", "Desulfobacterota": "#ea9999", "Deinococcota": "#f9cb9c",
                 "Cyanobacteria": "#ffe599", "Chloroflexota": "#b6d7a8", "Campylobacterota": "#a2c4c9",
                 "Bacteroidota": "#9fc5e8", "Actinobacteriota": "#b4a7d6", "Acidobacteriota": "#d5a6bd",
                 "NA": "#000000", "ERROR": "#ff0000"}

    def getRandomColor(self):
        red = random.randint(0, 255)
        green = random.randint(0, 255)
        blue = random.randint(0, 255)
        color_code = "#{:02x}{:02x}{:02x}".format(red, green, blue)
        return color_code

    # HGT--All fasta evolutionary Tree branch drawing requires calling the BuildTree() function in Tree.py, which will not be duplicated in CellHGT.py here
    @timeit
    def HGTSpeciesPlot(self, TreeAnno):

        HGTTemp = os.path.join(self.outputDir, 'HGTSpeciesPlot')
        # create directory
        if not os.path.exists(HGTTemp):
            os.makedirs(HGTTemp, exist_ok=True)

        #HGTNum = pd.read_csv(HGTNum, sep="\t", header=0)
        HGTNum = self.HGTNum
        HGTNumSpecies = HGTNum.groupby(['Species1', 'Species2']).size().reset_index(name='HGTNum')
        HGTNumSpecies.to_csv(os.path.join(HGTTemp, 'HGTNumSpecies.txt'), sep='\t', index=False)

        SpeciesAnno = pd.read_csv(TreeAnno, header=0, sep='\t')
        
        # Dynamically identify the ID column (prioritize SpeciesID, followed by Bin)
        if 'SpeciesID' in SpeciesAnno.columns:
            id_col = 'SpeciesID'
        elif 'Bin' in SpeciesAnno.columns:
            id_col = 'Bin'
        else:
            raise ValueError("TreeAnno file format error: missing 'SpeciesID' or 'Bin' column.")
            
        # Force conversion to a string and remove duplicates, and establish a dictionary
        SpeciesAnno[id_col] = SpeciesAnno[id_col].astype(str)
        SpeciesAnno = SpeciesAnno.drop_duplicates(subset=id_col)
        SpeciesAnno_dict = SpeciesAnno.set_index(id_col)['Phylum'].to_dict()

        SpeciesList = HGTNumSpecies['Species1'].tolist() + HGTNumSpecies['Species2'].tolist()
        SpeciesList = list(set(SpeciesList))
        SpeciesPhylum = SpeciesAnno[SpeciesAnno[id_col].isin(SpeciesList)].copy()

        Phylum = SpeciesPhylum['Phylum'].tolist()
        Phylum = list(set(Phylum))

        if 'Phylum_Color' in SpeciesPhylum.columns:
            PhylumColor = SpeciesPhylum.drop_duplicates(subset='Phylum')
            PhylumColor_dict = SpeciesPhylum.set_index('Phylum')['Phylum_Color'].to_dict()
        else:
            for p in Phylum:
                if p not in CellHGT.PHY_COLOR.keys():
                    # Assign a color to the uncommon Phylum and place it in the PHY_COLOR keyword
                    while True:
                        color = self.getRandomColor()
                        if color not in CellHGT.PHY_COLOR.values():
                            CellHGT.PHY_COLOR[p] = color
                            break
            PhylumColor_dict = CellHGT.PHY_COLOR
            SpeciesPhylum['Phylum_Color'] = SpeciesPhylum['Phylum'].apply(lambda s: PhylumColor_dict[s])

        HGTNodeAnno = os.path.join(HGTTemp, 'HGTNodeAnno.txt')
        HGTLineAnno = os.path.join(HGTTemp, 'HGTLineAnno.txt')

        # Fill in the HGTNodeAnno.txt file
        with open(HGTNodeAnno, 'w') as input1:  # The first time you open it, use 'w' to ensure that all old files are cleared
            input1.write('TREE_COLORS\nSEPARATOR TAB\nDATA\n')
            for index, row in SpeciesPhylum.iterrows():
                row_str = str(row[id_col]) + '\trange\t' + str(row['Phylum_Color']) + '\t' + str(row['Phylum']) + '\n'
                input1.write(row_str)

        # Fill in the HGTLineAnno.txt file
        with open(HGTLineAnno, 'w') as input1:  # Fill in the HGTLineAnno.txt file
            input1.write('DATASET_CONNECTION\nSEPARATOR COMMA\nDATASET_LABEL,example connections\nCOLOR,#ff0ff0\nDRAW_ARROWS,0\nLOOP_SIZE,100\nMAXIMUM_LINE_WIDTH,10\nCURVE_ANGLE,-50\nCENTER_CURVES,0\nALIGN_TO_LABELS,1\nDATA\n')

            # Determine whether the gates of HGTNumSpecies Species1 and Species2 are the same. If they are the same, the color of the line is the gate color; otherwise, it is gray
            # At this point, the file remains open and can be efficiently written directly into the loop
            for index, row in HGTNumSpecies.iterrows():
                Species1Phylum = SpeciesAnno_dict[row['Species1']]
                Species2Phylum = SpeciesAnno_dict[row['Species2']]

                if Species1Phylum == Species2Phylum:
                    HGTPhylum = Species1Phylum
                    Color = PhylumColor_dict[HGTPhylum]
                else:
                    HGTPhylum = 'Diff'
                    Color = 'grey'

                row_str = str(row['Species1']) + ',' + str(row['Species2']) + ',' + '1,' + Color + ',dashed,xxxx\n'
                input1.write(row_str)

    '''
    The HGTAnno function performs clustering [CDHit] and annotation 
    of contig files from the HGT.fasta files obtained by SpeciesHGT or StrainHGT.
    '''

    # !!!!!If you want to change the input contig length set by prokka, you should modify $MAXCONTIGIDLEN in the bin/prokka script

    # HGTStrainAnno requires that the contig naming format for the input HGTFasta file be: >Pair_ species 1@ Strain 1 to species 2@ strain 2
    @timeit
    def HGTStrainAnno(self, TreeAnno, prokka_env=None, cdhit_env=None, emapper_env=None, emapper_DB=None):

        HGTTemp = os.path.join(self.outputDir,'HGTStrainAnno')
        HGTFasta = self.HGTfasta
        AnnoProkka = os.path.join(HGTTemp, 'AnnoProkka')
        AnnoEggnog = os.path.join(HGTTemp, 'AnnoEggnog')
        AnnoCDHit = os.path.join(HGTTemp, 'AnnoCDHit')

        # create directory
        DirCreate = [AnnoProkka, AnnoEggnog, AnnoCDHit]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        SpeciesAnno = pd.read_csv(TreeAnno, header=0, sep='\t')
        
        # Dynamically identify the ID column (prioritize SpeciesID, followed by Bin)
        if 'SpeciesID' in SpeciesAnno.columns:
            id_col = 'SpeciesID'
        elif 'Bin' in SpeciesAnno.columns:
            id_col = 'Bin'
        else:
            raise ValueError("TreeAnno file format error: missing 'SpeciesID' or 'Bin' column.")
            
        # Force conversion to a string and remove duplicates, and establish a dictionary
        SpeciesAnno[id_col] = SpeciesAnno[id_col].astype(str)
        SpeciesAnno = SpeciesAnno.drop_duplicates(subset=id_col)
        SpeciesAnno_dict = SpeciesAnno.set_index(id_col)['Phylum'].to_dict()

        ######Prokka######
        command = 'prokka ' + HGTFasta + ' --outdir ' + AnnoProkka + '/HGT_prokka --prefix HGT --kingdom Bacteria'
        if prokka_env is not None:
            subprocess.run(['conda', 'run', '-n', prokka_env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        ######Eggnog-Mapper######
        if emapper_DB is not None:
            command = 'emapper.py -m diamond -i ' + AnnoProkka + '/HGT_prokka/HGT.faa --output ' + AnnoEggnog + '/HGTEggnog --data_dir ' + emapper_DB + ' --cpu 10'
            #command = 'emapper.py -m diamond -i ' + AnnoProkka + '/HGT_prokka/HGT.faa --output ' + AnnoEggnog + '/HGTEggnog --data_dir /data_alluser/public/database/eggnogDB/ --cpu 10'
        else:
            command = 'emapper.py -m diamond -i ' + AnnoProkka + '/HGT_prokka/HGT.faa --output ' + AnnoEggnog + '/HGTEggnog --cpu 10'

        if emapper_env is not None:
            subprocess.run(['conda', 'run', '-n', emapper_env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        ######CDHit######
        command = 'cdhit -i ' + AnnoProkka + '/HGT_prokka/HGT.faa' + ' -o ' + AnnoCDHit + '/HGTCDHit -c 1 -aS 0 -d 0'
        if cdhit_env is not None:
            subprocess.run(['conda', 'run', '-n', cdhit_env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        ######Annotation result organization#######
        # 1. The annotation table obtained from Prokka + emapper.py results, organized into a document ORFAnno.txt: Contig Bin1 Species1 Bin2 Species2 Gene PhylumToPhylum
        HGTgff = pd.read_csv(AnnoProkka + '/HGT_prokka/HGT.gff', sep='\t', header=None, comment='#')
        HGTEggnog = pd.read_csv(AnnoEggnog + '/HGTEggnog.emapper.annotations', sep='\t', header=None, comment='#',
                                index_col=0)

        ORFAnno = pd.DataFrame(columns=['Contig', 'Bin1', 'Phylum1', 'Bin2', 'Phylum2', 'ProkkaID', 'EggnogGene'])
        EggID2ContigID = {}

        for index, row in HGTgff.iterrows():

            Contig = row.iloc[0]

            if Contig.startswith('>'):
                break

            Bin1 = Contig.split('@')[0]
            Bin1 = Bin1.replace('PAIR_', '')
            Bin1 = Bin1.split('_')[0]
            Bin2 = Contig.split('TO')[1]
            Bin2 = Bin2.split('@')[0]
            Bin2 = Bin2.split('_')[0]

            Phylum1 = SpeciesAnno_dict[str(Bin1)]
            Phylum2 = SpeciesAnno_dict[str(Bin2)]

            ProkkaID = row.iloc[8]
            ProkkaID = ProkkaID.split(';')[0]
            ProkkaID = ProkkaID.replace('ID=', '')

            EggID2ContigID[ProkkaID] = Contig
            if ProkkaID in HGTEggnog.index:
                EggnogGene = HGTEggnog.loc[ProkkaID, 7]
                temp = [Contig, Bin1, Phylum1, Bin2, Phylum2, ProkkaID, EggnogGene]
                ORFAnno.loc[len(ORFAnno)] = temp

        ORFAnno.to_csv(AnnoEggnog + '/ORFAnno.txt', sep='\t', index=False)
        # 2 The obtained comment box by prokka (HGT.gff) + cdhit (HGTCDHit.clstr) result,Organize it into a document
        # CDHitCluster.txt : Cluster Gene Contig Bin1 [Species1] Bin2 [Species2]
        # ClusterBin.txt Cluster [Bin1 Bin2 Bin3 ... ...]
        HGTclstr = AnnoCDHit + '/HGTCDHit.clstr'
        CDHitCluster = pd.DataFrame(columns=['Cluster', 'Gene', 'Contig', 'Bin1', 'Bin2'])
        ClusterBin = {}
        with open(HGTclstr, 'r') as input1:
            while True:
                line = input1.readline()
                if len(line) == 0:
                    break

                line = line.replace('\n', '')
                if line.startswith('>Cluster'):
                    cluster = line.replace('>', '')
                else:
                    Gene = line.split(', ')
                    Gene = Gene[1].split('...')[0]
                    Gene = Gene.replace('>', '')

                    Contig = EggID2ContigID[Gene]

                    Bin1 = Contig.split('@')[0]
                    Bin1 = Bin1.replace('PAIR_', '')
                    Bin1 = Bin1.split('_')[0]
                    Bin2 = Contig.split('TO')[1]
                    Bin2 = Bin2.split('@')[0]
                    Bin2 = Bin2.split('_')[0]

                    temp = [cluster, Gene, Contig, Bin1, Bin2]
                    CDHitCluster.loc[len(CDHitCluster)] = temp

                    if cluster not in ClusterBin.keys():
                        ClusterBin[cluster] = [Bin1, Bin2]
                    else:
                        ClusterBin[cluster] += [Bin1, Bin2]
                        ClusterBin[cluster] = list(set(ClusterBin[cluster]))

        CDHitCluster.to_csv(AnnoCDHit + '/CDHitCluster.txt', sep='\t', index=False)
        for key, value in ClusterBin.items():
            cluster = key
            Bins = ','.join(value)
            with open(AnnoCDHit + '/ClusterBin.txt', 'a') as output:
                output.write(cluster + '\t' + Bins + '\n')

    # HGTSpeciesAnno requires that the contig naming format for the input HGTFasta file be: > pair_species1_contig1to species2_contig2
    @timeit
    def HGTSpeciesAnno(self, TreeAnno, prokka_env=None, cdhit_env=None, emapper_env=None, emapper_DB=None):

        HGTTemp = os.path.join(self.outputDir, 'HGTSpeciesAnno')
        HGTFasta = self.HGTfasta

        AnnoProkka = os.path.join(HGTTemp, 'AnnoProkka')
        AnnoEggnog = os.path.join(HGTTemp, 'AnnoEggnog')
        AnnoCDHit = os.path.join(HGTTemp, 'AnnoCDHit')

        # create directory
        DirCreate = [AnnoProkka, AnnoEggnog, AnnoCDHit]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        SpeciesAnno = pd.read_csv(TreeAnno, header=0, sep='\t')
        
        # Dynamically identify the ID column (prioritize SpeciesID, followed by Bin)
        if 'SpeciesID' in SpeciesAnno.columns:
            id_col = 'SpeciesID'
        elif 'Bin' in SpeciesAnno.columns:
            id_col = 'Bin'
        else:
            raise ValueError("TreeAnno file format error: missing 'SpeciesID' or 'Bin' column.")
            
        # Force conversion to a string and remove duplicates, and establish a dictionary
        SpeciesAnno[id_col] = SpeciesAnno[id_col].astype(str)
        SpeciesAnno = SpeciesAnno.drop_duplicates(subset=id_col)
        SpeciesAnno_dict = SpeciesAnno.set_index(id_col)['Phylum'].to_dict()

        ######Prokka######
        command = 'prokka ' + HGTFasta + ' --outdir ' + AnnoProkka + '/HGT_prokka --prefix HGT --kingdom Bacteria'
        if prokka_env is not None:
            subprocess.run(['conda', 'run', '-n', prokka_env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        ######Eggnog-Mapper######
        if emapper_DB is not None:
            command = 'emapper.py -m diamond -i ' + AnnoProkka + '/HGT_prokka/HGT.faa --output ' + AnnoEggnog + '/HGTEggnog --data_dir ' + emapper_DB + ' --cpu 10'
            #command = 'emapper.py -m diamond -i ' + AnnoProkka + '/HGT_prokka/HGT.faa --output ' + AnnoEggnog + '/HGTEggnog --data_dir /data_alluser/public/database/eggnogDB/ --cpu 10'
        else:
            command = 'emapper.py -m diamond -i ' + AnnoProkka + '/HGT_prokka/HGT.faa --output ' + AnnoEggnog + '/HGTEggnog --cpu 10'

        if emapper_env is not None:
            subprocess.run(['conda', 'run', '-n', emapper_env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        ######CDHit######
        command = 'cdhit -i ' + AnnoProkka + '/HGT_prokka/HGT.faa' + ' -o ' + AnnoCDHit + '/HGTCDHit -c 1 -aS 0 -d 0'
        if cdhit_env is not None:
            subprocess.run(['conda', 'run', '-n', cdhit_env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        ######Annotation result organization#######
        # 1 The comment box obtained by prokka and the result of emapper.py are organized into a file named ORFAnno.txt : Contig Bin1 Species 1 Bin2 Species 2 Gene phylum To phylum
        HGTgff = pd.read_csv(AnnoProkka + '/HGT_prokka/HGT.gff', sep='\t', header=None, comment='#')
        HGTEggnog = pd.read_csv(AnnoEggnog + '/HGTEggnog.emapper.annotations', sep='\t', header=None, comment='#',
                                index_col=0)

        ORFAnno = pd.DataFrame(columns=['Contig', 'Bin1', 'Phylum1', 'Bin2', 'Phylum2', 'ProkkaID', 'EggnogGene'])
        EggID2ContigID = {}

        for index, row in HGTgff.iterrows():

            Contig = row.iloc[0]

            if Contig.startswith('>'):
                break

            Bin1 = Contig.split('TO')[0]
            Bin1 = Bin1.replace('PAIR_', '')
            Bin1 = Bin1.split('_')[0]
            Bin2 = Contig.split('TO')[1]
            Bin2 = Bin2.split('_')[0]

            Phylum1 = SpeciesAnno_dict[str(Bin1)]
            Phylum2 = SpeciesAnno_dict[str(Bin2)]

            ProkkaID = row.iloc[8]
            ProkkaID = ProkkaID.split(';')[0]
            ProkkaID = ProkkaID.replace('ID=', '')

            EggID2ContigID[ProkkaID] = Contig
            if ProkkaID in HGTEggnog.index:
                EggnogGene = HGTEggnog.loc[ProkkaID, 7]
                temp = [Contig, Bin1, Phylum1, Bin2, Phylum2, ProkkaID, EggnogGene]
                ORFAnno.loc[len(ORFAnno)] = temp

        ORFAnno.to_csv(AnnoEggnog + '/ORFAnno.txt', sep='\t', index=False)
        # 2 The annotation box obtained by prokka (HGT.gff) + cdhit (HGTCDHit.clstr) result is organized
        # CDHitCluster.txt : Cluster Gene Contig Bin1 [Species1] Bin2 [Species2]
        # ClusterBin.txt Cluster [Bin1 Bin2 Bin3 ... ...]
        HGTclstr = AnnoCDHit + '/HGTCDHit.clstr'
        CDHitCluster = pd.DataFrame(columns=['Cluster', 'Gene', 'Contig', 'Bin1', 'Bin2'])
        ClusterBin = {}
        with open(HGTclstr, 'r') as input1:
            while True:
                line = input1.readline()
                if len(line) == 0:
                    break

                line = line.replace('\n', '')
                if line.startswith('>Cluster'):
                    cluster = line.replace('>', '')
                else:
                    Gene = line.split(', ')
                    Gene = Gene[1].split('...')[0]
                    Gene = Gene.replace('>', '')

                    Contig = EggID2ContigID[Gene]

                    Bin1 = Contig.split('TO')[0]
                    Bin1 = Bin1.replace('PAIR_', '')
                    Bin1 = Bin1.split('_')[0]
                    Bin2 = Contig.split('TO')[1]
                    Bin2 = Bin2.split('_')[0]

                    temp = [cluster, Gene, Contig, Bin1, Bin2]
                    CDHitCluster.loc[len(CDHitCluster)] = temp

                    if cluster not in ClusterBin.keys():
                        ClusterBin[cluster] = [Bin1, Bin2]
                    else:
                        ClusterBin[cluster] += [Bin1, Bin2]
                        ClusterBin[cluster] = list(set(ClusterBin[cluster]))

        CDHitCluster.to_csv(AnnoCDHit + '/CDHitCluster.txt', sep='\t', index=False)
        for key, value in ClusterBin.items():
            cluster = key
            Bins = ','.join(value)
            with open(AnnoCDHit + '/ClusterBin.txt', 'a') as output:
                output.write(cluster + '\t' + Bins + '\n')
