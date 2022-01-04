import pandas as pd
import numpy as np
import os
import threading
from tqdm import tqdm
import sys
from Threshold_Consensus import muscle_generation


def no_ref_consensus_finding_sparc(group_id: int):
    try:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_' + str(group_id) + '.csv')
    except FileNotFoundError:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_no_blat_' + str(group_id) + '.csv')
    unique_id = np.unique(df_seq_extract['qseqid'])
    sp_result_df = pd.DataFrame()
    sp_dummy_counter = 0
    
    for seqid in tqdm(unique_id):
        temp = df_seq_extract[df_seq_extract['qseqid'] == seqid]

        if temp.shape[0] > 0:
            file_name = str(seqid) + '.fasta'
            ref_file_name = 'ref_seq_' + str(seqid) + '.fasta'
            f = open(file_name, mode='w+')
            ref_f = open(ref_file_name, mode='w+')
            
            longest = temp.iloc[np.argmax(temp['corr_seq'].apply(len))]
            depth = temp.shape[0]
            
            ref_f.write('>' + str(longest['name']) + '\n')
            ref_f.write(str(longest['corr_seq']) + '\n')
            ref_f.close()
            
            for index, row in temp.iterrows():
                f.write('>' + str(row['name']) + '\n')
                f.write(str(row['corr_seq']) + '\n')
            f.close()
            
            thread_id = str(threading.get_ident())
            thread_id += file_name

            if depth == 1:
                sp_result_df.at[sp_dummy_counter, 'qseqid'] = seqid
                sp_result_df.at[sp_dummy_counter, 'depth'] = depth
                sp_result_df.at[sp_dummy_counter, 'seq'] = temp.iloc[0]['corr_seq']

                sp_dummy_counter += 1
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

                    sp_result_df.at[sp_dummy_counter, 'qseqid'] = seqid
                    sp_result_df.at[sp_dummy_counter, 'depth'] = depth
                    sp_result_df.at[sp_dummy_counter, 'seq'] = seq

                    sp_dummy_counter += 1
                    
                    os.remove(str(thread_id) + '.consensus.fasta')
                
                os.remove('mapped' + str(thread_id) + '.m5')
                
            os.remove(file_name)
            os.remove(ref_file_name)
            
    sp_result_df.to_csv('results/no_ref/Result_sparc_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' no ref Sparc Consensus Finding Finished.')
    
    
def no_ref_consensus_finding_sparc_lseq(group_id: int):
    try:
        df_seq_extract = pd.read_csv('temp_df/lseqs_df_' + str(group_id) + '.csv')
    except FileNotFoundError:
        print('No LSEQ file found')
        return
    unique_id = np.unique(df_seq_extract['qseqid'])
    sp_result_df = pd.DataFrame()
    sp_dummy_counter = 0
    
    for seqid in tqdm(unique_id):
        temp = df_seq_extract[df_seq_extract['qseqid'] == seqid]

        if temp.shape[0] > 0:
            file_name = str(seqid) + '_l.fasta'
            ref_file_name = 'ref_seq_' + str(seqid) + '_l.fasta'
            f = open(file_name, mode='w+')
            ref_f = open(ref_file_name, mode='w+')
            
            longest = temp.iloc[np.argmax(temp['lseq'].apply(len))]
            depth = temp.shape[0]
            
            ref_f.write('>' + str(longest['lname']) + '\n')
            ref_f.write(str(longest['lseq']) + '\n')
            ref_f.close()
            
            for index, row in temp.iterrows():
                f.write('>' + str(row['lname']) + '\n')
                f.write(str(row['lseq']) + '\n')
            f.close()
            
            thread_id = str(threading.get_ident())
            thread_id += file_name

            if depth == 1:
                sp_result_df.at[sp_dummy_counter, 'qseqid'] = seqid
                sp_result_df.at[sp_dummy_counter, 'depth'] = depth
                sp_result_df.at[sp_dummy_counter, 'seq'] = temp.iloc[0]['lseq']

                sp_dummy_counter += 1
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

                    sp_result_df.at[sp_dummy_counter, 'qseqid'] = seqid
                    sp_result_df.at[sp_dummy_counter, 'depth'] = depth
                    sp_result_df.at[sp_dummy_counter, 'seq'] = seq

                    sp_dummy_counter += 1
                    
                    os.remove(str(thread_id) + '.consensus.fasta')
                
                os.remove('mapped' + str(thread_id) + '.m5')
                
            os.remove(file_name)
            os.remove(ref_file_name)
            
    sp_result_df.to_csv('results/no_ref/Result_Lseq_sparc_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' no ref lseq Sparc Consensus Finding Finished.')
    
    
    
def no_ref_consensus_finding_pbdagcon(group_id: int):
    try:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_' + str(group_id) + '.csv')
    except FileNotFoundError:
        df_seq_extract = pd.read_csv('temp_df/df_w_seqs_no_blat_' + str(group_id) + '.csv')
    unique_id = np.unique(df_seq_extract['qseqid'])
    pb_result_df = pd.DataFrame()
    pb_dummy_counter = 0
    
    for seqid in tqdm(unique_id):
        temp = df_seq_extract[df_seq_extract['qseqid'] == seqid]

        if temp.shape[0] > 0:
            file_name = str(seqid) + '.fasta'
            ref_file_name = 'ref_seq_' + str(seqid) + '.fasta'
            f = open(file_name, mode='w+')
            ref_f = open(ref_file_name, mode='w+')
            
            longest = temp.iloc[np.argmax(temp['corr_seq'].apply(len))]
            depth = temp.shape[0]
            
            ref_f.write('>' + str(longest['name']) + '\n')
            ref_f.write(str(longest['corr_seq']) + '\n')
            ref_f.close()
            
            for index, row in temp.iterrows():
                f.write('>' + str(row['name']) + '\n')
                f.write(str(row['corr_seq']) + '\n')
            f.close()
            
            thread_id = str(threading.get_ident())
            thread_id += file_name

            if depth == 1:
                pb_result_df.at[pb_dummy_counter, 'qseqid'] = seqid
                pb_result_df.at[pb_dummy_counter, 'depth'] = depth
                pb_result_df.at[pb_dummy_counter, 'seq'] = temp.iloc[0]['corr_seq']

                pb_dummy_counter += 1
            elif depth > 1:
                
                os.system('blasr ' + file_name + ' ' + ref_file_name + ' --bestn 1 --minMatch 5 --placeGapConsistently -m 5 --out mapped' + str(thread_id) + '.m5 --nproc ' + str(os.cpu_count()))
                
                # Start of pbdagcon consensus generating
                
                os.system('pbdagcon --min-coverage 4 --min-length 100 --threads 20 mapped' + str(thread_id) + '.m5 > consensus' + str(thread_id) + '.fasta')
                
                consensus_file = open('consensus' + str(thread_id) + '.fasta', mode='r')
                seq = consensus_file.readlines()
                if len(seq) > 1:
                    seq = seq[1][:-1]
                else:
                    seq = ''
                consensus_file.close()

                pb_result_df.at[pb_dummy_counter, 'qseqid'] = seqid
                pb_result_df.at[pb_dummy_counter, 'depth'] = depth
                pb_result_df.at[pb_dummy_counter, 'seq'] = seq

                pb_dummy_counter += 1
                
                os.remove('consensus' + str(thread_id) + '.fasta')
                
                os.remove('mapped' + str(thread_id) + '.m5')
                
            os.remove(file_name)
            os.remove(ref_file_name)
            
    pb_result_df.to_csv('results/no_ref/Result_pbdagcon_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' no ref pbdagcon Consensus Finding Finished.')


def no_ref_consensus_finding_chris(group_id: int):
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

            """
            
            if temp.shape[0] > 70:
                temp = temp.loc[list(temp['corr_seq'].apply(len).sort_values(ascending=True).iloc[:70].index)]
                
            """

            for index, row in temp.iterrows():
                seq = row['corr_seq']
                f.write(">" + str(seqid) + '\n')
                f.write(str(seq) + '\n')
            f.close()

            if depth == 1:
                result_df.at[dummy_counter, 'qseqid'] = seqid
                result_df.at[dummy_counter, 'depth'] = depth
                result_df.at[dummy_counter, 'seq'] = seq

                dummy_counter += 1

                os.remove(file_name)
            elif depth == 0:
                os.remove(file_name)
            else:
                muscle_generation(file_name)
                consensus_file = open("Consensus_No_Dashes_" + file_name, mode='r')
                lines = consensus_file.readlines()
                if len(lines) == 1:
                    seq = lines[0]
                else:
                    seq = lines[1]
                consensus_file.close()
                result_df.at[dummy_counter, 'qseqid'] = seqid
                result_df.at[dummy_counter, 'depth'] = depth
                result_df.at[dummy_counter, 'seq'] = seq

                dummy_counter += 1

                os.remove(file_name)
                os.remove("Consensus_No_Dashes_" + file_name)
                os.remove("Consensus_" + file_name)
                os.remove("MSA_Consensus_" + file_name)
                # os.rename("MSA_Consensus_" + file_name, "MSA/MSA_Consensus_" + file_name)

    result_df.to_csv('results/no_ref/Result_chris_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' no ref Chris Consensus Finding Finished.')