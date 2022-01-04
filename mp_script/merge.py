import pandas as pd
from datetime import date


def merge(start: int, end: int, ref: bool, consensus_mode: str, fasta_src_name: str, lseq_enabled: bool):
    today = str(date.today())
    fasta_name = ''
    if ref:
        fasta_name = 'results/ref/fasta/consensus_ref_' + consensus_mode + '_' + today + '_' + fasta_src_name
        fasta = open(fasta_name, mode='w+')
    else:
        fasta_name = 'results/no_ref/fasta/consensus_no_ref_' + consensus_mode + '_' + today + '_' + fasta_src_name
        fasta = open(fasta_name, mode='w+')
    
    file_prior = ''
    if ref:
        if consensus_mode == 's':
            file_prior = 'results/ref/Result_sparc_'
        elif consensus_mode == 'p':
            file_prior = 'results/ref/Result_pbdagcon_'
        else:
            print('merging error.')
    else:
        if consensus_mode == 's':
            file_prior = 'results/no_ref/Result_sparc_'
        elif consensus_mode == 'p':
            file_prior = 'results/no_ref/Result_pbdagcon_'
        elif consensus_mode == 'c':
            file_prior = 'results/no_ref/Result_chris_'
        else:
            print('merging error.')
            
    total = pd.DataFrame()
    for i in range(start, end):
        total = pd.concat([total, pd.read_csv(file_prior + str(i) + '.csv').drop('Unnamed: 0', axis=1)])
    total.to_csv(file_prior + str(start) + '_' + str(end) + '_' + fasta_src_name + '.csv')
    
    for index, row in total.iterrows():
        if row['seq'] != '' and (not pd.isna(row['seq'])):
            fasta.write('>' + row['qseqid'] + '_' + str(int(row['depth'])) + '\n' + str(row['seq']) + '\n')
            
    fasta.close()
    
    if lseq_enabled:
        if ref:
            file_prior = 'results/ref/Result_Lseq_sparc_'
            fasta_name = 'results/ref/fasta/lseq_consensus_' + consensus_mode + '_' + today + '_' + fasta_src_name
        else:
            file_prior = 'results/no_ref/Result_Lseq_sparc_'
            fasta_name = 'results/no_ref/fasta/lseq_consensus_' + consensus_mode + '_' + today + '_' + fasta_src_name
        total = pd.DataFrame()
        for i in range(start, end):
            total = pd.concat([total, pd.read_csv(file_prior + str(i) + '.csv').drop('Unnamed: 0', axis=1)])
        total.to_csv(file_prior + str(start) + '_' + str(end) + '_' + fasta_src_name + '.csv')
        
        with open(fasta_name, mode='w+') as fasta:
            for index, row in total.iterrows():
                if row['seq'] != '' and (not pd.isna(row['seq'])):
                    fasta.write('>' + row['qseqid'] + '_' + str(int(row['depth'])) + '\n' + str(row['seq']) + '\n')
        
        file_prior = 'temp_df/Q_trimming_'
        total = pd.DataFrame()
        for i in range(start, end):
            total = pd.concat([total, pd.read_csv(file_prior + str(i) + '.csv').drop('Unnamed: 0', axis=1)])
        total.to_csv(file_prior + str(start) + '_' + str(end) + '_' + fasta_src_name + '.csv')
    
    print('Merge completed. Result directory: ' + fasta_name + '\n Exiting now...')