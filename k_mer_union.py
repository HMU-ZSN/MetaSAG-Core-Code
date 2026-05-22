import numpy as np
import os
import re

def merge_kmer_npy_files(folder_path, output_filename='merged_kmer_features.npy'):
        
    folder_path = os.path.normpath(folder_path)
    
    file_pattern = re.compile(r'k-mer_M(\d+)_Contig20\.npy')
    
    files_info = []
    for file_name in os.listdir(folder_path):
        match = file_pattern.match(file_name)
        if match:
            file_num = int(match.group(1))
            file_path = os.path.join(folder_path, file_name)
            files_info.append((file_num, file_path, file_name))
    
    if not files_info:
        raise FileNotFoundError(f"No k-mer_M*_Contig20.npy files were found in the folder {folder_path}")
    
   
    files_info.sort(key=lambda x: x[0])
    
    print(f"Find {len(files_info)} matching files:")
    for file_num, _, file_name in files_info:
        print(f"  {file_name}")
    
    
    first_num, first_path, first_name = files_info[0]
    print(f"\n is loading the initial file: {first_name}")
    merged_array = np.load(first_path)
    
    
    if len(merged_array.shape) == 1:
        merged_array = merged_array.reshape(1, -1)
    
    print(f"Initial array shape: {merged_array.shape}")
    
   
    for i in range(1, len(files_info)):
        file_num, file_path, file_name = files_info[i]
        print(f"Loading and merging files: {file_name}", end=' | ')
        
      
        current_array = np.load(file_path)
        
        
        if len(current_array.shape) == 1:
            current_array = current_array.reshape(1, -1)
        
        print(f"shape: {current_array.shape}", end=' | ')
        
       
        if merged_array.shape[1] != current_array.shape[1]:
            raise ValueError(f"The number of columns {current_array.shape[1]} of the file {file_name} does not match the number of columns {merged_array.shape[1]} of the merged array.")
        
      
        merged_array = np.vstack((merged_array, current_array))
        
       
        del current_array
        
        print(f"The overall shape after merging: {merged_array.shape}")
    
   
    output_path = os.path.join(folder_path, output_filename)
    np.save(output_path, merged_array)
    print(f"\nAll files have been merged successfully!")
    print(f"The final array shape: {merged_array.shape}")
    print(f"The result has been saved to: {output_path}")
    
    return merged_array

# Optional: Versions that need to be merged by number range
def merge_kmer_npy_files_by_range(folder_path, start_idx=1, end_idx=4718, output_filename='merged_kmer_features.npy'):
    
   
    folder_path = os.path.normpath(folder_path)
    
   
    file_pattern = re.compile(r'k-mer_M(\d+)_Contig20\.npy')
    
  
    files_info = []
    for file_name in os.listdir(folder_path):
        match = file_pattern.match(file_name)
        if match:
            file_num = int(match.group(1))
         
            if start_idx <= file_num <= end_idx:
                file_path = os.path.join(folder_path, file_name)
                files_info.append((file_num, file_path, file_name))
    
    if not files_info:
        raise FileNotFoundError(f"The k-mer_M*_Contig20.npy file with numbers ranging from {start_idx} to {end_idx} was not found in the folder {folder_path}")
    
   
    files_info.sort(key=lambda x: x[0])
    
    print(f"Find {len(files_info)} files within the number range from {start_idx} to {end_idx} :")
    for file_num, _, file_name in files_info:
        print(f"  {file_name}")
    
  
    file_numbers = [num for num, _, _ in files_info]
    missing_numbers = []
    for num in range(start_idx, end_idx + 1):
        if num not in file_numbers:
            missing_numbers.append(num)
    
    if missing_numbers:
        print(f"\nWarning: Files with the following numbers are missing: {missing_numbers}")
        print(f"Actually merge {len(files_info)} files instead of the expected {end_idx - start_idx + 1}")
    
  
    first_num, first_path, first_name = files_info[0]
    print(f"\nThe initial file is being loaded: {first_name}")
    merged_array = np.load(first_path)
    
    
    if len(merged_array.shape) == 1:
        merged_array = merged_array.reshape(1, -1)
    
    print(f"Initial array shape: {merged_array.shape}")
    
   
    for i in range(1, len(files_info)):
        file_num, file_path, file_name = files_info[i]
        print(f"Files are being loaded and merged: {file_name}", end=' | ')
        
       
        current_array = np.load(file_path)
        
       
        if len(current_array.shape) == 1:
            current_array = current_array.reshape(1, -1)
        
        print(f"shape: {current_array.shape}", end=' | ')
        
        # Check if the number of columns matches
        if merged_array.shape[1] != current_array.shape[1]:
            raise ValueError(f"The number of columns {current_array.shape[1]} of the file {file_name} does not match the number of columns {merged_array.shape[1]} of the merged array.")
        
        # Stack vertically along the first axis (row)
        merged_array = np.vstack((merged_array, current_array))
        
       
        del current_array
        
        print(f"The overall shape after merging: {merged_array.shape}")
    
    
    output_path = os.path.join(folder_path, output_filename)
    np.save(output_path, merged_array)
    print(f"\nAll files have been merged successfully！")
    print(f"The final array shape: {merged_array.shape}")
    print(f"The result has been saved to: {output_path}")
    
    return merged_array


if __name__ == "__main__":
    
    folder_to_search = "your_file.npy" 
    
   
    try:
        merged_data = merge_kmer_npy_files(folder_to_search)
    except FileNotFoundError as e:
        print(e)
    
   