import os
from shutil import which


def env_check():
    
    program_list = ['blasr', 'pbdagcon', 'Sparc', 'blat', 'blastn']
    for p in program_list:
        if which(p) is None:
            print(p + ' not found.')
            exit(1)
    
    dir_list = ['temp_df', 'blastn', 'results', 'results/ref', 'results/no_ref', 'results/ref/fasta', 'results/no_ref/fasta']
    for d in dir_list:
        if not os.path.exists(d):
            os.makedirs(d)
    
    print('Environment Check Finished. All Clear.')
    