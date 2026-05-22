
from Bio import SeqIO
import os
from collections import defaultdict
import numpy as np
import itertools
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA



record = next(SeqIO.parse("phage.fasta", "fasta"))
sequence = str(record.seq).upper()  
print(f"sequence length: {len(sequence)}, first 20 bases: {sequence[:20]}")



def process_fasta_files(folder_path, k=3):
    
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith(('.fasta', '.fa', '.fna'))]
    print(f"find{len(fasta_files)} FASTA files")

    
    all_vectors = []
    all_sequence_info = []

    
    bases = ['A', 'T', 'C', 'G']
    all_possible_kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]

    
    for file_name in fasta_files:
        file_path = os.path.join(folder_path, file_name)

        
        for record in SeqIO.parse(file_path, "fasta"):
            sequence = str(record.seq).upper()
            seq_id = record.id
            print(f'process sequence：{seq_id}, length：{len(sequence)}, The first 20 bases：{sequence[:20]}')

            
            kmers = get_kmers(sequence, k=k)
            freq = kmer_frequency(kmers)

           
            vector = [freq.get(kmer, 0) for kmer in all_possible_kmers]
            all_vectors.append(vector)
            all_sequence_info.append(seq_id)

    return np.array(all_vectors), all_sequence_info, all_possible_kmers
def get_kmers(seq, k=3):
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    return kmers

def kmer_frequency(kmers):
    counts = defaultdict(int)
    total = len(kmers)
    for kmer in kmers:
        counts[kmer] += 1
    
    freq = {k: v/total for k, v in counts.items()}
    return freq


kmers = get_kmers(sequence, k=3)
freq = kmer_frequency(kmers)
print(f"There are a total of {len(kmers)} 3-mers, among which the 'ATG' frequency: {freq.get('ATG', 0):.4f}")



def build_kmer_vector(freq, k=3):
    
    bases = ['A', 'T', 'C', 'G']
    all_kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    
    vector = [freq.get(kmer, 0) for kmer in all_kmers]
    return np.array(vector)


kmer_vector = build_kmer_vector(freq, k=3)
print(f"Feature vector dimension: {kmer_vector.shape}, sample value :\n{kmer_vector[:10]}")  


top_kmers = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:10]
x, y = zip(*top_kmers)
plt.bar(x, y)
plt.title("Top 10 3-mers in Phage Sequence")
plt.ylabel("Frequency")
plt.xticks(rotation=45)
plt.show()


X = np.array([kmer_vector])  
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)
print(f"The dimension after dimensionality reduction: {X_pca.shape}")