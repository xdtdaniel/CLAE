import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import threading
from Bio import Seq


def quality_trimming(Llen: int, seqs: dict, group_id: int):
    if os.path.exists('temp_df/Q_trimming_' + str(group_id) + '.csv'):
        print('temp_df/Q_trimming_' + str(group_id) + '.csv found. Now loading...')
        return
    df = pd.read_csv('blastn/blastn_' + str(group_id) + '.csv').drop('Unnamed: 0', axis=1)
    groups = df['qseqid'].unique()

    for qseqid in tqdm(groups):
        group_num = 1
        curr = df[df['qseqid'] == qseqid]
        qstart_prev = None
        for index, row in curr.iterrows():
            if qstart_prev is None:
                df.at[index, 'subgroup'] = 'Grp_' + str(group_num)
                df.at[index, 'qs-pqs'] = 0
                qstart_prev = row['qstart']
            else:
                qs_pqs = row['qstart'] - qstart_prev
                df.at[index, 'qs-pqs'] = qs_pqs
                qstart_prev = row['qstart']
                if qs_pqs >= 1.2 * Llen:
                    group_num += 1
                    df.at[index, 'subgroup'] = 'Grp_' + str(group_num)

    groups = df['qseqid'].unique()

    query_error_seqs = []

    # get each qseqid group
    for qseqid in tqdm(groups):

        if qseqid not in seqs:
            print(seqs.keys())
            print(qseqid)
            continue

        curr = df[df['qseqid'] == qseqid]
        error_count = 0
        group_start = list(curr.index[~curr['subgroup'].isna()]) + [curr.index[-1] + 1]
        group_count = len(group_start) - 1
        

        # process subgroups with the same qseqid
        for i in range(group_count):
            curr_group = group_start[i]
            next_group = group_start[i+1]
            row_count = next_group - curr_group
            grp_dist = df['qend'][next_group-1] - df['qstart'][curr_group]
            df.at[curr_group, '# rows'] = row_count
            df.at[curr_group, 'Grp_dist'] = grp_dist

            if row_count == 1:
                df.at[curr_group, 'Lst'] = df['qstart'][curr_group] - 120
                df.at[curr_group, 'Lend'] = df['qend'][next_group-1] + 120
                df.at[curr_group, 'Tst'] = df['qend'][next_group-1] + 20
                if group_start[i+1] == curr.index[-1] + 1:
                    df.at[curr_group, 'Tend'] = -np.inf
                else:
                    df.at[curr_group, 'Tend'] = df['qstart'][next_group] - 20
            elif row_count == 2:
                df.at[curr_group, 'Lst'] = df['qstart'][curr_group] - 100
                df.at[curr_group, 'Lend'] = df['qend'][next_group-1] + 100
                df.at[curr_group, 'Tst'] = df['qend'][next_group-1] + 20
                if group_start[i+1] == curr.index[-1] + 1:
                    df.at[curr_group, 'Tend'] = -np.inf
                else:
                    df.at[curr_group, 'Tend'] = df['qstart'][next_group] - 20
            elif 3 <= row_count <= 8:
                if grp_dist > 1.3 * Llen:
                    df.at[curr_group, 'Lst'] = -np.inf
                    df.at[curr_group, 'Lend'] = -np.inf
                    df.at[curr_group, 'Tst'] = -np.inf
                    df.at[curr_group, 'Tend'] = -np.inf

                    # mark an error
                    error_count += 1
                else:
                    df.at[curr_group, 'Lst'] = df['qstart'][curr_group] - 50
                    df.at[curr_group, 'Lend'] = df['qend'][next_group-1] + 50
                    df.at[curr_group, 'Tst'] = df['qend'][next_group-1] + 20
                    if group_start[i+1] == curr.index[-1] + 1:
                        df.at[curr_group, 'Tend'] = -np.inf
                    else:
                        df.at[curr_group, 'Tend'] = df['qstart'][next_group] - 20
            elif row_count > 8:
                df.at[curr_group, 'Lst'] = -np.inf
                df.at[curr_group, 'Lend'] = -np.inf
                df.at[curr_group, 'Tst'] = -np.inf
                df.at[curr_group, 'Tend'] = -np.inf

                # mark an error
                error_count += 1

            if df['Lst'][curr_group] != -np.inf and df['Lend'][curr_group] != -np.inf:
                curr_seq = seqs[qseqid][int(df['Lst'][curr_group]):int(df['Lend'][curr_group])+1]
                df.at[curr_group, 'Lseq'] = curr_seq
                df.at[curr_group, 'Lseq_len'] = len(curr_seq)
            else:
                df.at[curr_group, 'Lseq_len'] = np.nan

            if df['Tst'][curr_group] != -np.inf and df['Tend'][curr_group] != -np.inf:
                curr_seq = seqs[qseqid][int(df['Tst'][curr_group]):int(df['Tend'][curr_group])+1]
                df.at[curr_group, 'Tseq'] = curr_seq
                df.at[curr_group, 'Tseq_len'] = len(curr_seq)
            else:
                df.at[curr_group, 'Tseq_len'] = np.nan

        tseq_len_med = np.median(df[df['qseqid'] == qseqid]['Tseq_len'].dropna())

        if error_count > 0.5 * group_count or (error_count > 0.35 * group_count and tseq_len_med < 300):
            query_error_seqs.append(qseqid)

    df = df[~df['qseqid'].isin(query_error_seqs)]
    df.reset_index(drop=True)
    
    groups = df['qseqid'].unique()
    query_error_seqs = []

    for qseqid in tqdm(groups):
        curr = df[df['qseqid'] == qseqid]
        group_start = list(curr.index[~curr['subgroup'].isna()]) + [curr.index[-1] + 1]
        group_count = len(group_start) - 1

        lseq_len_med = np.median(df[df['qseqid'] == qseqid]['Lseq_len'].dropna())
        tseq_len_med = np.median(df[df['qseqid'] == qseqid]['Tseq_len'].dropna())

        for i in range(group_count):
            curr_group = group_start[i]

            df.at[curr_group, 'Median LseqLen'] = lseq_len_med
            df.at[curr_group, 'Median TseqLen'] = tseq_len_med

            lseq_pieces = df['Lseq_len'][curr_group] / lseq_len_med
            tseq_pieces = df['Tseq_len'][curr_group] / tseq_len_med
            df.at[curr_group, 'Lseq:medi'] = lseq_pieces
            df.at[curr_group, 'Tseq:medi'] = tseq_pieces

            # trimming
            
            # Tseq processing

            # 2 pieces
            if 1.5 <= tseq_pieces <= 2.5:
                index_to_trim = int(df['Tseq_len'][curr_group] / 2 - 0.5 * Llen)
                seq_to_trim = df['Tseq'][curr_group]
                df.at[curr_group, '1st_seqid'] = df['subgroup'][curr_group] + '_1'
                df.at[curr_group, '1st_seq'] = seq_to_trim[:index_to_trim]
                df.at[curr_group, '2nd_seqid'] = df['subgroup'][curr_group] + '_2'
                df.at[curr_group, '2nd_seq'] = seq_to_trim[-index_to_trim:]

            elif tseq_pieces > 2.5:
                index_to_trim = int(tseq_len_med - 20)
                seq_to_trim = df['Tseq'][curr_group]
                df.at[curr_group, '1st_seqid'] = df['subgroup'][curr_group] + '_1'
                df.at[curr_group, '1st_seq'] = seq_to_trim[:index_to_trim]
                df.at[curr_group, '2nd_seqid'] = df['subgroup'][curr_group] + '_2'
                df.at[curr_group, '2nd_seq'] = seq_to_trim[-index_to_trim:]
                df.at[curr_group, '3rd_seqid'] = df['subgroup'][curr_group] + '_3'
                df.at[curr_group, '3rd_seq'] = seq_to_trim[index_to_trim:-index_to_trim]
                

        curr = df[df['qseqid'] == qseqid]
        total_tseq = curr[~curr['Tseq:medi'].isna()]
        qual_tseq = total_tseq[((0.75 <= total_tseq['Tseq:medi']) & (total_tseq['Tseq:medi'] <= 1.15)) | ((1.7 <= total_tseq['Tseq:medi']) & (total_tseq['Tseq:medi'] <= 2.2)) | ((2.7 <= total_tseq['Tseq:medi']) & (total_tseq['Tseq:medi'] <= 3.4))]
        if total_tseq.shape[0] == 0 or qual_tseq.shape[0] / total_tseq.shape[0] < 0.6:
            query_error_seqs.append(qseqid)

    df_correct = df[~df['qseqid'].isin(query_error_seqs)]
    df_correct.reset_index(drop=True)
    df_correct.to_csv('temp_df/Q_trimming_' + str(group_id) + '.csv')
    print('Group ' + str(group_id) + ' Quality Trimming Finished.')


def seq_extraction(group_id: int, enable_phred_threshold: bool, fastq_filename: str, phred_threshold: int, mugio_path: str, plus_or_minus: str):
    if os.path.exists('temp_df/df_w_seqs_' + str(group_id) + '.csv'):
        print('temp_df/df_w_seqs_' + str(group_id) + '.csv found. Now loading...')
        return

    df = pd.read_csv('temp_df/Q_trimming_' + str(group_id) + '.csv').drop('Unnamed: 0', axis=1)
    
    result_columns = ['match', 'mismatch', 'rep_match', 'N_s', 'Q_gap_count', 'Q_gap_bases', 'T_gap_count', 'T_gap_bases', 'strand', 'Q_name', 'Q_size', 'Q_start', 'Q_end', 'T_name', 'T_size', 'T_start', 'T_end', 'block_count', 'block_size', 'qStarts', 'tStarts']
    
    # step 1
    groups = df['qseqid'].unique()
    seq_list = []
    lseq_list = []

    for qseqid in tqdm(groups):
        curr = df[df['qseqid'] == qseqid]
        curr = curr[~curr['subgroup'].isna()]

        for index, row in curr.iterrows():
        
            # Lseq processing
            
            lname = row['qseqid'] + '_' + row['subgroup'] + '_L'
            lseq = row['Lseq']
            
            if not pd.isna(lseq) and lseq != '':
                lseq_list.append((qseqid, lname, lseq))

        
            # Tseq processing
            tseq_len_med = row['Median TseqLen']

            # two pieces
            if 1.5 <= row['Tseq:medi'] <= 2.5:
                name_1 = row['qseqid'] + '_' + row['1st_seqid']
                seq_1 = row['1st_seq']
                try:
                    seq1_start = int(row['Tst'])
                    seq1_end = seq1_start + len(seq) - 1
                except:
                    seq1_start = -np.inf
                    seq1_end = -np.inf
                
                name_2 = row['qseqid'] + '_' + row['2nd_seqid']
                seq_2 = row['2nd_seq']
                try:
                    seq2_start = seq1_end + 1
                    seq2_end = seq2_start + len(seq_2) - 1
                except:
                    seq2_start = -np.inf
                    seq2_end = -np.inf
                
                seq_list.append((qseqid, name_1, seq_1, tseq_len_med, seq1_start, seq1_end))
                seq_list.append((qseqid, name_2, seq_2, tseq_len_med, seq2_start, seq2_end))
            # three pieces
            elif row['Tseq:medi'] > 2.5:
                name_1 = row['qseqid'] + '_' + row['1st_seqid']
                seq_1 = row['1st_seq']
                try:
                    seq1_start = int(row['Tst'])
                    seq1_end = seq1_start + len(seq) - 1
                except:
                    seq1_start = -np.inf
                    seq1_end = -np.inf
                
                name_2 = row['qseqid'] + '_' + row['2nd_seqid']
                seq_2 = row['2nd_seq']
                try:
                    seq2_start = seq1_end + 1
                    seq2_end = seq2_start + len(seq_2) - 1
                except:
                    seq2_start = -np.inf
                    seq2_end = -np.inf
                
                name_3 = row['qseqid'] + '_' + row['3rd_seqid']
                seq_3 = row['3rd_seq']
                try:
                    seq3_start = seq2_end + 1
                    seq3_end = seq3_start + len(seq_3) - 1
                except:
                    seq3_start = -np.inf
                    seq3_end = -np.inf
                
                
                seq_list.append((qseqid, name_1, seq_1, tseq_len_med, seq1_start, seq1_end))
                seq_list.append((qseqid, name_2, seq_2, tseq_len_med, seq2_start, seq2_end))
                if 0.7 * tseq_len_med <= len(seq_3) <= 1.3 * tseq_len_med:
                    seq_list.append((qseqid, name_3, seq_3, tseq_len_med, seq3_start, seq3_end))
            # one piece
            else:
                name = row['qseqid'] + '_' + row['subgroup']
                seq = row['Tseq']
                try:
                    seq1_start = int(row['Tst'])
                    seq1_end = seq1_start + len(seq) - 1
                except:
                    seq1_start = -np.inf
                    seq1_end = -np.inf
                seq_list.append((qseqid, name, seq, tseq_len_med, seq1_start, seq1_end))

    seq_df = pd.DataFrame(seq_list, columns=['qseqid', 'name', 'seq', 'Tlenmed', 'TStart', 'TEnd']).dropna().reset_index(drop=True)
    lseq_df = pd.DataFrame(lseq_list, columns=['qseqid', 'lname', 'lseq']).dropna().reset_index(drop=True)

    # get each group of qseqid seperately

    groups = seq_df['qseqid'].unique()

    for qseqid in tqdm(groups):
        file_path = os.getcwd()

        curr = seq_df[seq_df['qseqid'] == qseqid]
        tseq_len_med = list(curr['Tlenmed'])[0]

        # get all good seqs
        good_seqs = list(curr[(0.9 * tseq_len_med <= curr['seq'].apply(len).to_numpy()) & (curr['seq'].apply(len).to_numpy() <= 1.1 * tseq_len_med)]['seq'])
        good_seqs.sort(key=len, reverse=True)

        ref_seq_usable = False
        
        thread_id = threading.get_ident() + np.random.rand()

        # get queries into files
        query_file = open(file_path + '/query' + str(thread_id) + str(qseqid) + '.fasta', mode='w+')
        for index, row in curr.iterrows():
            query_file.write('>' + row['name'] + '\n')
            query_file.write(row['seq'] + '\n')
        query_file.close()

        for j in range(len(good_seqs)):
            # if there are multiple good seqs, we iterate through them until we find a one that is useable
            if len(good_seqs) == 0:
                ref_seq = list(curr['seq'])[0]
            else:
                ref_seq = good_seqs[j]

            # write ref_seq
            
            ref_file = open(file_path + '/ref' + str(thread_id) + str(qseqid) + '.fasta', mode='w+')
            ref_file.write('>' + qseqid + '\n')
            ref_file.write(ref_seq)
            ref_file.close()

            os.system('blat ' + file_path + '/ref' + str(thread_id) + str(qseqid) + '.fasta ' + file_path + '/query' + str(thread_id) + str(qseqid) + '.fasta ' + file_path + '/output' + str(thread_id) + str(qseqid) + '.psl ' + 'minIdentity=0 minScore=0 tileSize=6')

            # read blat information
            if os.stat(file_path + '/output' + str(thread_id) + str(qseqid) + '.psl').st_size == 0:
                os.remove(file_path + '/ref' + str(thread_id) + str(qseqid) + '.fasta')
                os.remove(file_path + '/output' + str(thread_id) + str(qseqid) + '.psl')
                continue
            blat_result = pd.read_csv(file_path + '/output' + str(thread_id) + str(qseqid) + '.psl', sep='\t', index_col=None).reset_index(drop=False)
            blat_result.columns = result_columns
            blat_result = blat_result[3:].drop_duplicates(subset='Q_name', keep='first').reset_index(drop=True).astype({'match': 'float64', 'Q_size': 'float64', 'Q_name': 'str'})

            curr_avg_match = np.average(blat_result['match'])
            if curr_avg_match < 250:
                os.remove(file_path + '/ref' + str(thread_id) + str(qseqid) + '.fasta')
                continue

            ref_seq_usable = True
            for index, row in blat_result.iterrows():
                curr_name = row['Q_name']
                seq_df.at[seq_df[seq_df['name'] == curr_name].index, 'strand'] = row['strand']

            # deleting files
            os.remove(file_path + '/ref' + str(thread_id) + str(qseqid) + '.fasta')
            os.remove(file_path + '/query' + str(thread_id) + str(qseqid) + '.fasta')
            os.remove(file_path + '/output' + str(thread_id) + str(qseqid) + '.psl')
            break

        if not ref_seq_usable:
            os.remove(file_path + '/query' + str(thread_id) + str(qseqid) + '.fasta')

    seq_df = seq_df.dropna(subset=['strand']).reset_index(drop=True)

    for index, row in seq_df.iterrows():
        if plus_or_minus == 'both':
            if row['strand'] == '+':
                seq_df.at[index, 'corr_seq'] = row['seq']
            else:
                seq = Seq.Seq(row['seq'])
                rev = seq.reverse_complement()
                seq_df.at[index, 'corr_seq'] = str(rev)
                
        elif plus_or_minus == 'minus':
            if row['strand'] == '-':
                seq_df.at[index, 'corr_seq'] = row['seq']
                
        elif plus_or_minus == 'plus':
            if row['strand'] == '+':
                seq_df.at[index, 'corr_seq'] = row['seq']
                
    seq_df = seq_df[(seq_df['corr_seq'] != '') & ~(seq_df['corr_seq'].isna())]
     
    # If phred score generation enabled
    
    if enable_phred_threshold and fastq_filename != '' and phred_threshold != 0:
        seq_df.loc[seq_df['TStart'] < 1, 'TStart'] = 0
        seq_df.loc[seq_df['TEnd'] < 1, 'TEnd'] = 0
        groups = seq_df['qseqid'].unique()
    
    
        for qseqid in tqdm(groups):
            phred_score_file_name = str(qseqid) + '_phred.jpg'
            os.system(f'python {mugio_path} -pp -f ' + fastq_filename + ' -uid ' + str(qseqid) + ' -o ' + phred_score_file_name)
            
            # remove the phred image as we do not need it now
            if os.path.exists(phred_score_file_name):
                os.remove(phred_score_file_name)
                
            if os.path.exists(phred_score_file_name + '.csv'):
                
                phred_df = pd.read_csv(phred_score_file_name + '.csv')
                os.remove(phred_score_file_name + '.csv')
                curr_subgroup = seq_df[seq_df['qseqid'] == qseqid]
                
                all_mean = phred_df['Score'].mean()
                all_std = phred_df['Score'].std()
                
                
                for index, row in curr_subgroup.iterrows():
                    curr_mean = phred_df.iloc[int(np.min([row['TStart'], row['TEnd']])):int(np.max([row['TStart'], row['TEnd']]))]['Score'].mean()
                    
                    seq_df.at[index, 'avg_score'] = curr_mean
                    seq_df.at[index, 'normalized_avg_score'] = (curr_mean - all_mean) / all_std
                        
        seq_df = seq_df[seq_df['avg_score'] >= phred_threshold]

    seq_df.to_csv('temp_df/df_w_seqs_' + str(group_id) + '.csv')
    lseq_df.to_csv('temp_df/lseqs_df_' + str(group_id) + '.csv')


def seq_extraction_no_blat(group_id: int):
    
    if os.path.exists('temp_df/df_w_seqs_no_blat_' + str(group_id) + '.csv'):
        print('temp_df/df_w_seqs_no_blat_' + str(group_id) + '.csv found. Now loading...')
        return

    df = pd.read_csv('temp_df/Q_trimming_' + str(group_id) + '.csv').drop('Unnamed: 0', axis=1)
    
    # step 1
    groups = df['qseqid'].unique()
    seq_list = []
    lseq_list = []

    for qseqid in tqdm(groups):
        curr = df[df['qseqid'] == qseqid]
        curr = curr[~curr['subgroup'].isna()]

        for index, row in curr.iterrows():
        
            # Lseq processing
            
            lname = row['qseqid'] + '_' + row['subgroup'] + '_L'
            lseq = row['Lseq']
            
            if not pd.isna(lseq) and lseq != '':
                lseq_list.append((qseqid, lname, lseq))
        
            tseq_len_med = row['Median TseqLen']

            # two pieces
            if 1.5 <= row['Tseq:medi'] <= 2.5:
                name_1 = row['qseqid'] + '_' + row['1st_seqid']
                seq_1 = row['1st_seq']
                try:
                    seq1_start = int(row['Tst'])
                    seq1_end = seq1_start + len(seq) - 1
                except:
                    seq1_start = -np.inf
                    seq1_end = -np.inf
                
                name_2 = row['qseqid'] + '_' + row['2nd_seqid']
                seq_2 = row['2nd_seq']
                try:
                    seq2_start = seq1_end + 1
                    seq2_end = seq2_start + len(seq_2) - 1
                except:
                    seq2_start = -np.inf
                    seq2_end = -np.inf
                
                seq_list.append((qseqid, name_1, seq_1, tseq_len_med, seq1_start, seq1_end))
                seq_list.append((qseqid, name_2, seq_2, tseq_len_med, seq2_start, seq2_end))
            # three pieces
            elif row['Tseq:medi'] > 2.5:
                name_1 = row['qseqid'] + '_' + row['1st_seqid']
                seq_1 = row['1st_seq']
                try:
                    seq1_start = int(row['Tst'])
                    seq1_end = seq1_start + len(seq) - 1
                except:
                    seq1_start = -np.inf
                    seq1_end = -np.inf
                
                name_2 = row['qseqid'] + '_' + row['2nd_seqid']
                seq_2 = row['2nd_seq']
                try:
                    seq2_start = seq1_end + 1
                    seq2_end = seq2_start + len(seq_2) - 1
                except:
                    seq2_start = -np.inf
                    seq2_end = -np.inf
                
                name_3 = row['qseqid'] + '_' + row['3rd_seqid']
                seq_3 = row['3rd_seq']
                try:
                    seq3_start = seq2_end + 1
                    seq3_end = seq3_start + len(seq_3) - 1
                except:
                    seq3_start = -np.inf
                    seq3_end = -np.inf
                
                
                seq_list.append((qseqid, name_1, seq_1, tseq_len_med, seq1_start, seq1_end))
                seq_list.append((qseqid, name_2, seq_2, tseq_len_med, seq2_start, seq2_end))
                if 0.7 * tseq_len_med <= len(seq_3) <= 1.3 * tseq_len_med:
                    seq_list.append((qseqid, name_3, seq_3, tseq_len_med, seq3_start, seq3_end))
            # one piece
            else:
                name = row['qseqid'] + '_' + row['subgroup']
                seq = row['Tseq']
                try:
                    seq1_start = int(row['Tst'])
                    seq1_end = seq1_start + len(seq) - 1
                except:
                    seq1_start = -np.inf
                    seq1_end = -np.inf
                seq_list.append((qseqid, name, seq, tseq_len_med, seq1_start, seq1_end))

    seq_df = pd.DataFrame(seq_list, columns=['qseqid', 'name', 'seq', 'Tlenmed', 'TStart', 'TEnd']).dropna().reset_index(drop=True)
    lseq_df = pd.DataFrame(lseq_list, columns=['qseqid', 'lname', 'lseq']).dropna().reset_index(drop=True)

    for index, row in seq_df.iterrows():
        seq_df.at[index, 'corr_seq'] = row['seq']
        
    seq_df.to_csv('temp_df/df_w_seqs_no_blat_' + str(group_id) + '.csv')
    lseq_df.to_csv('temp_df/lseqs_df_' + str(group_id) + '.csv')