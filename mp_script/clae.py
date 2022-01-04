import shutil
from env_check import env_check
from blast_reader import read_blastn_and_divide
from seq_reader import read_seqs
from subread import quality_trimming, seq_extraction_no_blat
from consensus_ref import consensus_finding_sparc, consensus_finding_pbdagcon
from consensus_no_ref import no_ref_consensus_finding_pbdagcon, no_ref_consensus_finding_sparc
from merge import merge
import logging, time, os
from datetime import datetime


import sys, os, argparse, warnings
from multiprocessing import Pool

def main():
    warnings.filterwarnings("ignore")
    
    
    
    parser = argparse.ArgumentParser('Subread Extractor and Consensus Sequence Generator')
    parser.add_argument('--blast', type=str, help='Blast File Name .csv: ', required=True)
    parser.add_argument('--comma', help='Include to use comma (,) seperated blast result instead of tab seperated blast', action='store_true')
    parser.add_argument('--seq', help='RCA reads input file name', type=str, required=True)
    parser.add_argument('--llen', help='Specify Llen', type=int, default=160)
    parser.add_argument('--ref', help='Include to use reference mode', action='store_true')
    parser.add_argument('--start', help='Start Chunk ID, inclusive, for debugging purposes', type=int)
    parser.add_argument('--end', help='End Chunk ID, inclusive, for debugging purposes', type=int)
    parser.add_argument('--refseq', help='Ref Sequence Filename', type=str, required='--ref' in sys.argv)
    parser.add_argument('--algo', help='Consensus Tool Selection, s for Sparc, p for pbdagcon', type=str, required=True)
    parser.add_argument('--merge', help='Include to merge autodivided chunks. RECOMMEND TO INCLUDE', action='store_true')

    

    args = parser.parse_args()
    
    blast_path = os.path.abspath(args.blast)
    seq_path = os.path.abspath(args.seq)
    ref_seq_path = os.path.abspath(args.refseq) if args.refseq != None else ''
    
    
    root_dir = os.getcwd()
    
    exec_dirname = f'{os.path.basename(args.seq)}_{datetime.now().strftime("%m%d%Y_%H%M%S")}_{args.algo}_files'
    if not os.path.exists(exec_dirname):
        os.makedirs(exec_dirname)
    os.chdir(exec_dirname)
    
    logging.basicConfig(filename='consensus.log', level=logging.DEBUG)
    
    logging.info("==================================================")
    logging.info(">>> Running Log <<<")
    logging.info("==================================================")
    t1_start = time.perf_counter()
    t2_start = time.process_time()

    sep = '\t'
    if args.comma:
        sep = ','

    env_check()
    chunks = read_blastn_and_divide(blast_path, sep)

    start = 0
    end = chunks

    if args.start is not None:
        start = args.start
    if args.end is not None:
        end = args.end + 1

    seqs = read_seqs(seq_path)
    
    pool = Pool()

    # Quality Trimming

    for i in range(start, end):
        pool.apply_async(quality_trimming, args=(args.llen, seqs, i))
        
    pool.close()
    pool.join()
    
    pool = Pool()

    for i in range(start, end):
        pool.apply_async(seq_extraction_no_blat, args=(i,))
            
    pool.close()
    pool.join()

    # Consensus Finding
    
    pool = Pool()

    for i in range(start, end):

        # Ref mode
        if args.ref:
            if 's' in args.algo:
                pool.apply_async(consensus_finding_sparc, args=(i, ref_seq_path))
            if 'p' in args.algo:
                pool.apply_async(consensus_finding_pbdagcon, args=(i, ref_seq_path))
        # No ref mode
        else:
            if 's' in args.algo:
                pool.apply_async(no_ref_consensus_finding_sparc, args=(i,))
            if 'p' in args.algo:
                pool.apply_async(no_ref_consensus_finding_pbdagcon, args=(i,))

    pool.close()
    pool.join()
    
    # clean files generated by blasr
    for f in [os.path.join(dp, fi) for dp, dn, fn in os.walk(os.path.expanduser(root_dir)) for fi in fn]:
        if 'core.' in f or '.psl' in f:
            os.remove(f)
    
    
    if args.merge:
        merge(start, end, args.ref, args.algo, os.path.basename(args.seq), False)
        
    t1_stop = time.perf_counter()
    t2_stop = time.process_time()
    logging.info("--------------------------------------------------")
    logging.info("Elapsed time: %.1f [min]" % ((t1_stop-t1_start)/60))
    logging.info("CPU process time: %.1f [min]" % ((t2_stop-t2_start)/60))
    logging.info("--------------------------------------------------")

    logging.info("Temp file cleanup...")

    # Clean-up temp files
    if os.path.isdir('temp_df'):
        shutil.rmtree('temp_df')
        
    if os.path.isdir('blastn'):
        shutil.rmtree('blastn')

    logging.info('Done.')
    logging.info("--------------------------------------------------")

if __name__ == "__main__":
    main()




