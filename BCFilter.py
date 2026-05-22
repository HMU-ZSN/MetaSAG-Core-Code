import os
import sys
import time
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statsmodels.distributions.empirical_distribution import ECDF
import warnings






class BCFilter():
    def __init__(self,inputfastq,outputdir):
        self.inputfastq = inputfastq
        self.outputdir = outputdir

        # Determine whether inputfastq is single-ended or double-ended
        if type(inputfastq) == str and inputfastq.endswith('.fastq'):
            self.ReadsEnd = 'Single'
            self.inputfastq = inputfastq
        elif type(inputfastq) == list:
            if len(inputfastq) == 2 and inputfastq[0].endswith('_R1.fastq') and inputfastq[1].endswith('_R2.fastq'):
                self.ReadsEnd = 'Pair'
                self.inputfastq1 = inputfastq[0]
                self.inputfastq2 = inputfastq[1]

            elif len(inputfastq) == 2 and inputfastq[0].endswith('_R2.fastq') and inputfastq[1].endswith('_R1.fastq'):
                self.ReadsEnd = 'Pair'
                self.inputfastq1 = inputfastq[1]
                self.inputfastq2 = inputfastq[0]

            elif len(inputfastq) == 1 and inputfastq[0].endswith('.fastq'):
                self.ReadsEnd = 'Single'
                self.inputfastq = inputfastq
            else:
                warnings.warn("The input fastq file does not meet the format requirements", UserWarning)


        else:
            warnings.warn("The input fastq file does not meet the format requirements",UserWarning)


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
    def CellCountStatistic(self):
        ######Step1,prefix_bcread statistic######

        #Create the result folder

        outputDir = self.outputdir
        ReadsEnd = self.ReadsEnd

        if ReadsEnd == 'Single':
            inputFastq = self.inputfastq
        else:
            inputFastq = self.inputfastq1


        if not os.path.exists(outputDir):
            os.makedirs(outputDir,exist_ok=True)


        BC_Count={}
        #read inputFastq
        with open(inputFastq,'r') as f:
            while True:
                line = f.readline()
                if len(line) == 0:
                    break

                if line.startswith('@'):
                    line = line.replace('@','')
                    cell = line.split(':')[0]
                    if cell not in BC_Count.keys():
                        BC_Count[cell]=1
                    else:
                        BC_Count[cell]+=1

        BC_Count=pd.DataFrame(list(BC_Count.items()))
        BC_Count.columns=['cell','num']
        BC_Count.to_csv(os.path.join(outputDir,'bcread.txt'),sep='\t',index=False)
        self.BC_Count=BC_Count

    @timeit
    def getMinReads(self,FirstDrv_xlim=[-1000,5000],SecondDrv_xlim=[-1000,5000]):
        outputDir = self.outputdir
        ######Step2,Cumulative distribution derivative plot######
        #BC_Count = pd.read_csv(BC_Count,sep='\t',header=0)
        BC_Count = self.BC_Count
        BC_Count_Sorted = BC_Count.sort_values('num')
        total_sum=BC_Count_Sorted['num'].sum()
        BC_Count_Sorted['consum'] = np.cumsum(BC_Count_Sorted['num'])
        BC_Count_Sorted['conprop'] = BC_Count_Sorted['consum']/total_sum

        x=[i for i in BC_Count_Sorted['num'] if i !=0]
        ecdf = ECDF(x)
        cell_readsnum=ecdf.x
        cumprop=ecdf.y

        plt.figure(figsize=(8, 8))
        plt.step(cell_readsnum, cumprop, where='post', label='ECDF')
        plt.xlabel('Cell_ReadsNum')
        plt.ylabel('CumProp')
        plt.title('Empirical Cumulative Distribution Function (ECDF)')
        plt.legend(loc='best')
        plt.grid(True)
        plt.savefig(os.path.join(outputDir,'CDF_BCReads.pdf'))
        plt.close()
        #plt.show()

        ######Step3,The maximum optimization algorithm for min reads######
        reads_prop_data = []  # 1. Create an empty list for collecting data

        for num in BC_Count_Sorted['num'].unique():
            temp1=num
            temp2=float(ecdf(num))
            # 2. Build the data from each loop into a dictionary and append it to the list
            reads_prop_data.append({'reads_num': temp1, 'prop': temp2})

        # 3. After the loop ends, convert the list to a DataFrame at one time
        reads_prop = pd.DataFrame(reads_prop_data)


        reads_prop['n']=reads_prop['prop']-(reads_prop['reads_num'])/(len(reads_prop))
        flag = reads_prop.loc[reads_prop['n'] != np.inf, 'n'].max()
        min_reads=reads_prop.loc[(reads_prop['n'] == flag) & (reads_prop['n'] != np.inf),'reads_num'].values[0]



        ######Step4,Calculate the first-order derivative and the second-order derivative and plot them######
        #First-order derivative of statistics
        d1=[]
        for i in range(len(reads_prop)):
            if i==0:
                d=reads_prop.iloc[0]['prop']/1
            else:
                d=(reads_prop.iloc[i]['prop']-reads_prop.iloc[i-1]['prop'])/(reads_prop.iloc[i]['reads_num']-reads_prop.iloc[i-1]['reads_num'])

            d1+=[d]

        #d1PlotData=pd.DataFrame({'reads_num':reads_prop['reads_num'],'d1':d1})
        x=reads_prop['reads_num']
        y=d1
        plt.scatter(x,y,s=1)
        xlim1=FirstDrv_xlim[0]
        xlim2=FirstDrv_xlim[1]
        plt.xlim(xlim1,xlim2)
        plt.xlabel('reads_num')
        plt.ylabel('D1')
        plt.title('Scatter Plot of reads number and D1')
        plt.grid(True)
        plt.savefig(os.path.join(outputDir,'1stDerivate.pdf'))
        plt.close()
        #plt.show()
        #Statistical second-order derivative

        d2=[]
        for i in range(len(d1)):
            if i==0:
                d=d1[0]
            else:
                d=(d1[i]-d1[i-1])/1
            d2+=[d]


        x=reads_prop['reads_num']
        y=d2
        plt.scatter(x,y,s=1)
        xlim1 = SecondDrv_xlim[0]
        xlim2 = SecondDrv_xlim[1]
        plt.xlim(xlim1,xlim2)
        plt.xlabel('reads_num')
        plt.ylabel('D2')
        plt.title('Scatter Plot of reads number and D2')
        plt.grid(True)
        plt.savefig(os.path.join(outputDir,'2ndDerivate.pdf'))
        plt.close()
        #plt.show()

        self.min_reads = min_reads

    def BCMinReads(self, filter_outdir=None, filter_outname=None):
        
        BC_Count = self.BC_Count
        min_reads = self.min_reads
        
        # 1. The output location of the statistical file (keeping the original location unchanged)
        outputDir = self.outputdir

        ######Step5, Count the cell numbers obtained through filtration######
        BC_Count_Filter = BC_Count[BC_Count['num'] >= min_reads]
        BC_Count_Filter.to_csv(os.path.join(outputDir, 'BC_Count_Filter.txt'), sep='\t', index=False)
        self.BC_Count_Filter = BC_Count_Filter

        ######Step6, Generate the filtered FASTQ file######
        # 2. Decide on the storage directory for the filtered files
        if filter_outdir is None:
            filter_outdir = outputDir  # default location
        else:
            # Adjust according to relative path and absolute path (relative path is based on outputDir by default)
            if not os.path.isabs(filter_outdir):
                filter_outdir = os.path.join(outputDir, filter_outdir)
                
            if not os.path.exists(filter_outdir):
                os.makedirs(filter_outdir, exist_ok=True)

        # 3. Determine the input files that need to be processed
        if self.ReadsEnd == 'Single':
            input_files = [self.inputfastq]
        else:
            input_files = [self.inputfastq1, self.inputfastq2]

        # 4. Extract the set of valid cells to ensure ultra-fast search of O(1)
        valid_cells = set(BC_Count_Filter['cell'].values)

        for ifq in input_files:
            # 5. Determine the output name of the current FASTQ file
            # Determine the base prefix. The default is BCFilterFile; otherwise, use the input filter_outname (without the suffix).
            base_name = 'BCFilterFile' if filter_outname is None else filter_outname
            
            # Add suffixes uniformly at the end according to the situation
            if self.ReadsEnd == 'Single':
                filename = f"{base_name}.fastq"
            else:
                if ifq == self.inputfastq1:
                    filename = f"{base_name}_R1.fastq"
                else:
                    filename = f"{base_name}_R2.fastq"

            # Combine the final complete output path
            BCFilterFile = os.path.join(filter_outdir, filename)

            # 6. Safe and efficient file block reading and writing
            with open(ifq, 'r') as f_in, open(BCFilterFile, 'w') as f_out:
                while True:
                    line1 = f_in.readline()
                    if len(line1) == 0:  
                        break
                    line2 = f_in.readline()
                    line3 = f_in.readline()
                    line4 = f_in.readline()

                    line1_clean = line1.strip('\n').replace('@', '')
                    cell = line1_clean.split(':')[0]

                    if cell in valid_cells:
                        f_out.write(line1)
                        f_out.write(line2)
                        f_out.write(line3)
                        f_out.write(line4)

    @timeit
    def Sam2Cell(self):
        inputfastq=self.inputfastq
        CellBarnDir=os.path.join(self.outputdir,'CellBarnDir')
        #Create the result directory of CellBarnDir
        if not os.path.exists(CellBarnDir):
            os.makedirs(CellBarnDir,exist_ok=True)


        for ifq in self.inputfastq:
            if ifq.endswith('_R1.fastq'):
                suffix = '_R1.fastq'
            elif ifq.endswith('_R2.fastq'):
                suffix = '_R2.fastq'
            elif ifq.endswith('.fastq'):
                suffix = '.fastq'

            with open(ifq,'r') as f:
                while True:
                    line = f.readline()
                    if len(line) == 0:
                        break
                    line = line.strip('\n')

                    if line.startswith('@'):
                        line1 = line.replace('@','')
                        cell = line1.split(':')[0]

                    with open(os.path.join(CellBarnDir,cell)+suffix,'a') as output:
                        output.write(line+'\n')






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
def SamsBCFilterStack(Sam,SavedCell,AllCell,outputDir):
    #AllCell = [3132, 4557, 2321, 3433, 4500]
    #SavedCell = [2000,3000,1734,2222,2000]
    #Sam=['S5','S6','S7','S8','S10']
    fig, ax = plt.subplots()
    ax.bar(Sam,AllCell,label='Filtered Cells')
    ax.bar(Sam, SavedCell, label='Saved Cells')
    ax.set_title('Stacked Bar Chart')
    ax.set_xlabel('Sam')
    ax.set_ylabel('Cell Number')
    ax.legend()

    # Display graphics
    plt.savefig(os.path.join(outputDir,'CellStack.pdf'))
    plt.close()
    #plt.show()





