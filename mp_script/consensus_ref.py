import pandas as pd
import numpy as np
import os
import threading
from tqdm import tqdm

def consensus_finding_sparc_lseq(group_id: int, ref_file_name: str):
    try:
        df_seq_extract = pd.read_csv('temp_df/lseqs_df_' + str(group_id) + '.csv')
    except FileNotFoundError:
        print('No LSEQ file found')
        return

    unique_id = np.unique(df_seq_extract['qseqid'])
    result_df = pd.DataFrame()
    dummy_counter = 0

    for seqid in tqdm(unique_id):
        temp = df_seq_extract[df_seq_extract['qseqid'] == seqid]

        if temp.shape[0] > 0:
            file_name = str(seqid) + '_l.fasta'
            f = open(file_name, mode='w+')
            depth = temp.shape[0]
            
            for index, row in temp.iterrows():
                f.write(">" + str(row['lname']) + '\n')
                f.write(str(row['lseq']) + '\n')
            f.close()
            
            thread_id = str(threading.get_ident())
            thread_id = thread_id + file_name + '_L'

            if depth == 1:
                result_df.at[dummy_counter, 'qseqid'] = seqid
                result_df.at[dummy_counter, 'depth'] = depth
                result_df.at[dummy_counter, 'seq'] = temp.iloc[0]['lseq']

                dummy_counter += 1
            elif depth > 1:
                
                os.system('blasr ' + file_name + ' ' + ref_file_name + ' --bestn 1 --minMatch 5 --placeGapConsistently -m 5 --out mapped' + str(thread_id) + '.m5 --nproc ' + str(os.cpu_count()))
                
                # Start of sparc consensus generating
                
                ret_val = os.system('Sparc b ' + ref_file_name + ' m mapped' + str(thread_id) + '.m5 c 2 k 2 g 1 o ' + str(thread_id))
                
                if os.path.exists(str(thread_id) + '.consensus.fasta'):
                    consensus_file = open(str(thread_id) + '.consensus.fasta', mode='r')
                    seq = consensus_file.readlines()
                    if len(seq) > 1 and ret_val == 0:
                        seq = seq[1][:-1]
                    else:
                        seq = ''
                    consensus_file.close()

                    result_df.at[dummy_counter, 'qseqid'] = seqid
                    result_df.at[dummy_counter, 'depth'] = depth
                    result_df.at[dummy_counter, 'seq'] = seq

                    dummy_counter += 1
                    
                    os.remove(str(thread_id) + '.consensus.fasta')
  
                os.remove('mapped' + str(thread_id) + '.m5')

            os.remove(file_name)
        
    result_df.to_csv('results/ref/Result_Lseq_sparc_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' Sparc ref lseq Consensus Finding Finished.')

def consensus_finding_sparc(group_id: int, ref_file_name: str):
    try:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_' + str(group_id) + '.csv')
    except FileNotFoundError:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_no_blat_' + str(group_id) + '.csv')

    unique_id = np.unique(df_seq_extract['qseqid'])
    result_df = pd.DataFrame()
    dummy_counter = 0

    for seqid in tqdm(unique_id):
        temp = df_seq_extract[df_seq_extract['qseqid'] == seqid]

        if temp.shape[0] > 0:
            file_name = str(seqid) + '.fasta'
            f = open(file_name, mode='w+')
            depth = temp.shape[0]
            
            for index, row in temp.iterrows():
                f.write(">" + str(row['name']) + '\n')
                f.write(str(row['corr_seq']) + '\n')
            f.close()
            
            thread_id = str(threading.get_ident())
            thread_id += file_name

            if depth == 1:
                result_df.at[dummy_counter, 'qseqid'] = seqid
                result_df.at[dummy_counter, 'depth'] = depth
                result_df.at[dummy_counter, 'seq'] = temp.iloc[0]['corr_seq']

                dummy_counter += 1
            elif depth > 1:
                
                os.system('blasr ' + file_name + ' ' + ref_file_name + ' --bestn 1 --minMatch 5 --placeGapConsistently -m 5 --out mapped' + str(thread_id) + '.m5 --nproc ' + str(os.cpu_count()))
                
                # Start of sparc consensus generating
                
                ret_val = os.system('Sparc b ' + ref_file_name + ' m mapped' + str(thread_id) + '.m5 c 2 k 2 g 2 o ' + str(thread_id))
                
                if os.path.exists(str(thread_id) + '.consensus.fasta'):
                    consensus_file = open(str(thread_id) + '.consensus.fasta', mode='r')
                    seq = consensus_file.readlines()
                    if len(seq) > 1 and ret_val == 0:
                        seq = seq[1][:-1]
                    else:
                        seq = ''
                    consensus_file.close()

                    result_df.at[dummy_counter, 'qseqid'] = seqid
                    result_df.at[dummy_counter, 'depth'] = depth
                    result_df.at[dummy_counter, 'seq'] = seq

                    dummy_counter += 1
                    
                    os.remove(str(thread_id) + '.consensus.fasta')
  
                os.remove('mapped' + str(thread_id) + '.m5')

            os.remove(file_name)
        
    result_df.to_csv('results/ref/Result_sparc_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' Sparc ref Consensus Finding Finished.')


def consensus_finding_pbdagcon(group_id: int, ref_file_name: str):
    try:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_' + str(group_id) + '.csv')
    except FileNotFoundError:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_no_blat_' + str(group_id) + '.csv')
    unique_id = np.unique(df_seq_extract['qseqid'])
    result_df = pd.DataFrame()
    dummy_counter = 0

    for seqid in tqdm(unique_id):
        temp = df_seq_extract[df_seq_extract['qseqid'] == seqid]

        if temp.shape[0] > 0:
            file_name = str(seqid) + '.fasta'
            f = open(file_name, mode='w+')
            depth = temp.shape[0]
            
            for index, row in temp.iterrows():
                f.write(">" + str(row['name']) + '\n')
                f.write(str(row['corr_seq']) + '\n')
            f.close()
            
            thread_id = str(threading.get_ident())
            thread_id += file_name

            if depth == 1:
                result_df.at[dummy_counter, 'qseqid'] = seqid
                result_df.at[dummy_counter, 'depth'] = depth
                result_df.at[dummy_counter, 'seq'] = temp.iloc[0]['corr_seq']

                dummy_counter += 1
            elif depth > 1:
                
                os.system('blasr ' + file_name + ' ' + ref_file_name + ' --bestn 1 --minMatch 5 --placeGapConsistently -m 5 --out mapped' + str(thread_id) + '.m5 --nproc ' + str(os.cpu_count()))
                
                os.system('pbdagcon --min-coverage 4 --min-length 100 --threads 20 mapped' + str(thread_id) + '.m5 > consensus' + str(thread_id) + '.fasta')
                
                consensus_file = open('consensus' + str(thread_id) + '.fasta', mode='r')
                seq = consensus_file.readlines()
                if len(seq) > 1:
                    seq = seq[1][:-1]
                else:
                    seq = ''
                consensus_file.close()

                result_df.at[dummy_counter, 'qseqid'] = seqid
                result_df.at[dummy_counter, 'depth'] = depth
                result_df.at[dummy_counter, 'seq'] = seq

                dummy_counter += 1
                
                os.remove('mapped' + str(thread_id) + '.m5')
                os.remove('consensus' + str(thread_id) + '.fasta')
                
            os.remove(file_name)
                
    result_df.to_csv('results/ref/Result_sparc_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' Sparc ref Consensus Finding Finished.')

