import os
import pandas as pd
import subprocess
import time

'''

Script comment
The purpose of this script is:
For the input fastq file [different types of fastq files will be introduced later], first perform Diamond comparison to obtain the Uniref fragment annotation of the fastq file mapping. 
1. Organize the Diamond result file to obtain the Cell-Uniref Count matrix;
2. [Seurat] For the input Cell-Uniref Count matrix, use the seuratSBg.r script to Cluster the cells and perform SBG-cluster statistical plotting of the cells (the plotting part has not been implemented yet).
3. [HUMAnN] organize the Diamond annotation results of each Cluster cell, perform HUMAnN pathway annotation, and use PathClusterHeatmap.R to draw the coverage_rate heatmap of the cluster-path.

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



class HP():
    def __init__(self,FastqDir,ResultDir):
        self.FastqDir = FastqDir
        self.ResultDir = ResultDir
        if not os.path.exists(self.ResultDir):
            os.makedirs(self.ResultDir,exist_ok=True)

    PYTHONDIR = os.path.dirname(os.path.abspath(__file__))

    @timeit
    def Diamond(self, DiamondDB=None, Diamondenv=None):

        FastqDir = self.FastqDir
        DiamondDir = os.path.join(self.ResultDir,'DiamondDir')
        self.DiamondDir = DiamondDir

        # The Fastq files input by this function can be one or multiple.

        # create directory
        if not os.path.exists(DiamondDir):
            os.makedirs(DiamondDir, exist_ok=True)

        # Traverse the fastq file and perform Diamond to obtain xxx uniref
        # --- The newly added intelligent path judgment logic begins---
        fastq_files_to_process = []
        
        if os.path.isfile(FastqDir):
            if FastqDir.endswith('.fastq'):
                fastq_files_to_process.append(FastqDir)
            else:
                raise ValueError(f"The input single file is not in.fastq format: {FastqDir}")
                
        elif os.path.isdir(FastqDir):
            files = os.listdir(FastqDir)
            for filename in files:
                if filename.endswith('.fastq'):
                    fastq_files_to_process.append(os.path.join(FastqDir, filename))
        else:
            raise ValueError(f"Invalid path: {FastqDir} is neither a file nor a folder")
        # --- The newly added intelligent path judgment logic has been completed ---

        # Traverse the generated list of items to be processed for Diamond comparison
        for fastqfile in fastq_files_to_process:
            filename = os.path.basename(fastqfile) 
            DiamondFile = os.path.join(DiamondDir, filename[:-6] + '_uniref')
            if DiamondDB is not None:
                command = 'diamond blastx --db ' + DiamondDB + ' -q ' + fastqfile + ' -o ' + DiamondFile + ' -k 1'
            else:
                command = 'diamond blastx -q ' + fastqfile + ' -o ' + DiamondFile + ' -k 1'

            if Diamondenv is not None:
                subprocess.run(['conda', 'run', '-n', Diamondenv, 'bash', '-c', command])
            else:
                subprocess.call(command, shell=True)


    @timeit
    def Uniref2Matrix(self, MinUnirefNum=5, MinCellNum=5):
        DiamondDir = self.DiamondDir
        MatrixDir = os.path.join(self.ResultDir,'MatrixDir')
        self.MatrixDir = MatrixDir

        # create directory
        if not os.path.exists(MatrixDir):
            os.makedirs(MatrixDir)

        # Input all the Uniref results obtained from the comparison and organize them into count cell uniref file (the diamond_uniref file name to which the cell belongs)
        AllSGBCellUnirefCount = pd.DataFrame(columns=['count', 'cell', 'uniref', 'file'])
        for unirefFile in os.listdir(DiamondDir):
            if unirefFile.endswith('_uniref'):
                unirefFileDir = os.path.join(DiamondDir, unirefFile)

                CellUnirefCount = {}
                with open(unirefFileDir) as input1:
                    while True:
                        line = input1.readline()
                        if len(line) == 0:
                            break

                        items = line.split('\t')
                        cell = items[0].split(':')[0]
                        uniref = items[1].split('|')[0]
                        ident = float(items[2])

                        key = str(cell) + '@' + str(uniref)
                        if ident > 90:
                            if key in CellUnirefCount.keys():
                                CellUnirefCount[key] += 1
                            else:
                                CellUnirefCount[key] = 1

                CellUnirefCount = pd.DataFrame(list(CellUnirefCount.items()))
                df = CellUnirefCount[0].str.split('@', expand=True)
                CellUnirefCount.columns = ['cell_uniref', 'count']
                CellUnirefCount["cell"] = df[0]
                CellUnirefCount['uniref'] = df[1]
                CellUnirefCount = CellUnirefCount.drop(columns='cell_uniref')
                CellUnirefCount['file'] = unirefFile
                CellUnirefCount.to_csv(os.path.join(MatrixDir, unirefFile[:-7] + '_CellUnirefCount.txt'), sep='\t',
                                       index=False)
                AllSGBCellUnirefCount = pd.concat([AllSGBCellUnirefCount, CellUnirefCount], ignore_index=True)

        AllSGBCellUnirefCount.to_csv(os.path.join(MatrixDir, 'AllSGBCellUnirefCount.txt'), sep='\t', index=False)
        self.AllSGBCellUnirefCount = AllSGBCellUnirefCount

        # Organize all the cells into a cell-uniref count matrix and output the file
        # Remove Uniref90 that appears in less than 5 cells
        data = AllSGBCellUnirefCount
        Uniref10 = data['uniref'].value_counts()
        Uniref10Filter = [uniref for uniref in Uniref10.index if Uniref10[uniref] < MinCellNum]
        data = data[~data['uniref'].isin(Uniref10Filter)]
        # Remove cells with less than 5 Unirefs
        Cell10 = data['cell'].value_counts()
        Cell10Filter = [cell for cell in Cell10.index if Cell10[cell] < MinUnirefNum]
        data = data[~data['cell'].isin(Cell10Filter)]
        matrix = data.pivot(index='uniref', columns='cell', values='count').fillna(0)
        self.AllSGB_matrix = os.path.join(MatrixDir, 'AllSGB_matrix')
        matrix.to_csv(os.path.join(MatrixDir, 'AllSGB_matrix'), sep='\t')

    # Input AllSGB_matrix for cell clustering and obtain the clustering results
    @timeit
    def SeuratCluster(self,nfeatures=20000,npcs=30,FindNeighborsDim="1:10",FindClustersRes=0.5,UmapPDFWidth=6,UmapPDFHeight=4,DimPlotMethod='umap',FindAllMarkersPCT=0.25,FindAllMarkerslogFC=0.25,topN=10,HeatmapPDFWidth=12,HeatmapPDFHeight=8):

        SeuratResult = os.path.join(self.ResultDir,'SeuratResult')
        AllSGB_matrix = self.AllSGB_matrix
        # create directory
        if not os.path.exists(SeuratResult):
            os.makedirs(SeuratResult, exist_ok=True)
        SeuratSGB_R = os.path.join(self.PYTHONDIR, 'SeuratSGB.R')
        #command = 'Rscript ' + SeuratSGB_R + ' -i ' + AllSGB_matrix + ' -o ' + SeuratResult
        command =f'Rscript {SeuratSGB_R} -i {AllSGB_matrix} -o {SeuratResult} --nfeatures {nfeatures} --npcs {npcs} --FindNeighborsDim {FindNeighborsDim} --FindClustersRes {FindClustersRes} --UmapPDFWidth {UmapPDFWidth} --UmapPDFHeight {UmapPDFHeight} --DimPlotMethod {DimPlotMethod} --FindAllMarkersPCT {FindAllMarkersPCT} --FindAllMarkerslogFC {FindAllMarkerslogFC} --topN {topN} --HeatmapPDFWidth {HeatmapPDFWidth} --HeatmapPDFHeight {HeatmapPDFHeight}'
        subprocess.call(command, shell=True)

    @timeit
    def getPathCovCluster(self,HUMAnNFolderDirt):
        HUMAnNFolders = [folder for folder in os.listdir(HUMAnNFolderDirt) if folder.endswith('_humann')]

        PathCovCluster = pd.DataFrame(columns=['Path', 'Coverage', 'Cluster'])
        for folder in HUMAnNFolders:
            cluster = folder.replace('_humann', '')
            for file_name in os.listdir(os.path.join(HUMAnNFolderDirt, folder)):
                if file_name.endswith('_pathcoverage.tsv'):
                    file_path = os.path.join(HUMAnNFolderDirt, folder, file_name)
                    temp = pd.read_csv(file_path, sep='\t', header=None, comment='#')
                    temp.columns = ['Path', 'Coverage']
                    temp = temp[(~temp['Path'].str.contains('\|')) & (temp['Coverage'] > 0.1)]
                    # temp = temp[(~temp['Path'].str.contains('UNMAPPED')) | (~temp['Path'].str.contains('UNINTEGRATED'))]
                    temp = temp[~temp['Path'].str.contains('UNMAPPED|UNINTEGRATED')]
                    temp['Cluster'] = cluster

                    PathCovCluster = pd.concat([PathCovCluster, temp], ignore_index=True)
        PathCovCluster.to_csv(os.path.join(HUMAnNFolderDirt, 'PathCovCluster.txt'), sep='\t', index=False)

    # The uniref file of grouped cells：
    # Enter the CellAnno file and then input the column names of the cell groups (e.g. "SGB ", "Cluster"," Sam ", "DiseaseState ", it is recommended to use" SGB "or "Cluster").
    # The required input file：CellAnno
    @timeit
    def HUMAnNPath(self, CellAnno, Group, ExcludeGroups=None, HUMAnNenv=None):

        DiamondDir = self.DiamondDir
        HUMAnNResult = os.path.join(self.ResultDir,'HUMAnNResult')

        UnirefForHUMAnN = os.path.join(HUMAnNResult, 'UnirefForHUMAnN')
        HUMAnNPathForGroup = os.path.join(HUMAnNResult, 'HUMAnNPathForGroup')
        HeatMapPlot = os.path.join(HUMAnNResult, 'HeatMapPlot')

        # create directory
        DirCreate = [HUMAnNResult, UnirefForHUMAnN, HUMAnNPathForGroup, HeatMapPlot]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        # --- New addition: Handle blacklist parameters ---
        if ExcludeGroups is None:
            ExcludeGroups = []
        elif isinstance(ExcludeGroups, str):
            ExcludeGroups = [ExcludeGroups]
        elif not isinstance(ExcludeGroups, list):
            ExcludeGroups = list(ExcludeGroups)

        # Retrieve all uniref files under DiamondDir, and assign the Unirefs to the corresponding group's UnirefForHUMAnN/group_uniref files based on the CellAnno files and groups
        CellAnno_df = pd.read_csv(CellAnno, sep='\t', header=0)
        Cell_Group = dict(zip(CellAnno_df["Cell"], CellAnno_df[Group]))
        
        # --- New addition: Initialize the collection used to collect missing cells ---
        missing_cells = set()

        for file in os.listdir(DiamondDir):
            if file.endswith('_uniref'):
                UnirefFile = os.path.join(DiamondDir, file)
                with open(UnirefFile) as input1:
                    while True:
                        line = input1.readline()
                        if len(line) == 0:
                            break

                        cell = line.split(':')[0]
                        
                        # --- New addition: Secure extraction and interception logic ---
                        group = Cell_Group.get(cell, None)
                        
                        if group is None:
                            missing_cells.add(cell) # Record and print them collectively below
                            continue
                            
                        if group in ExcludeGroups:
                            continue                # If you hit the blacklist, skip it directly
                        # --------------------------------

                        groupUnirefFile = os.path.join(UnirefForHUMAnN, str(group) + '_uniref')

                        with open(groupUnirefFile, 'a') as output:
                            output.write(line)

        # --- New addition: Centralized printing of statistical information on unmatched cells ---
        if missing_cells:
            print(f"\n Notice: Detected {len(missing_cells)} cells without grouping information in the annotation file; they have been safely skipped.")
            print(f"Unmatched cell list (first 50): {list(missing_cells)[:50]} ...\n")

        # Annotate the HUMAnN path for each uniref file in UnirefForHUMAnN, and place the resulting files in the HUMAnNPathForGroup folder
        files = os.listdir(UnirefForHUMAnN)
        for filename in files:
            if filename.endswith('_uniref'):
                unirefFileDir = os.path.join(UnirefForHUMAnN, filename)
                humannFileDir = os.path.join(HUMAnNPathForGroup, filename.replace('_uniref', '')) + '_humann'
                command = 'humann --input ' + unirefFileDir + ' --output ' + humannFileDir
                if HUMAnNenv is not None:
                    subprocess.run(['conda', 'run', '-n', HUMAnNenv, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

        # Integrate all the group_humann files under the HUMAnNPathForGroup folder
        self.getPathCovCluster(HUMAnNPathForGroup)
        # Draw a Heatmap and place the related results of the heatmap in the HeatMapPlot folder
        PathCovCluster = os.path.join(HUMAnNPathForGroup, 'PathCovCluster.txt')
        PathClusterHeatmap_R = os.path.join(self.PYTHONDIR, 'PathClusterHeatmap.R')
        command = 'Rscript ' + PathClusterHeatmap_R + ' -i ' + PathCovCluster + ' -o ' + HeatMapPlot + '/'
        subprocess.call(command, shell=True)

