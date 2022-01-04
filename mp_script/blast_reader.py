import os
import pandas as pd
import numpy as np


def divide_check(chunks):
    for i in range(chunks):
        if not os.path.exists('blastn/blastn_' + str(i) + '.csv'):
            return False
    return True


def read_blastn_and_divide(blastn_file_name: str, sep: str):

    ONE_CHUNK_SIZE = 1e7

    chunks = int(os.path.getsize(blastn_file_name) // ONE_CHUNK_SIZE + 1)
    
    chunks = chunks if chunks > os.cpu_count() else os.cpu_count()
    
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'evalue', 'bitscore', 'qlen', 'slen', 'qseq', 'sseq']
    df = pd.read_csv(blastn_file_name, sep=sep, names=columns, header=None)
    # drop seqs
    df = df.drop(df.columns[-2:], axis=1)
    df = df[df['evalue'] <= 0.02]
    df = df.sort_values(['qseqid', 'qstart'], ascending=[True, True]).reset_index(drop=True)
    
    if len(np.unique(df['qseqid'])) < chunks:
        chunks = len(np.unique(df['qseqid']))

    if not divide_check(chunks):
        print('Data not divided yet, dividing...')
        divide_data(df, chunks)
        print('Finish Dividing.')
    
    return chunks

    
def divide_data(df, chunks: int):
    groups = df['qseqid'].unique()

    for i in range(chunks):
        df_chunk = df[df['qseqid'].isin(groups[int(i/chunks*len(groups)):int((i+1)/chunks*len(groups))])]
        df_chunk.to_csv('blastn/blastn_' + str(i) + '.csv')