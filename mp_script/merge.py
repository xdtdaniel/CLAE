import pandas as pd
from datetime import datetime
import os

def merge(start: int, end: int, ref: bool, consensus_mode: str, fasta_src_name: str, lseq_enabled: bool):
    today = str(datetime.now().strftime("%m%d%Y_%H%M%S"))
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
    
    print('Merge completed. Result directory: ' + fasta_name + '\n')
    
    for i in range(start, end):
        if (os.path.exists(file_prior + str(i) + '.csv')):
            os.remove(file_prior + str(i) + '.csv')