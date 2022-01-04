import os
import sys
import argparse


def generate_db(in_name: str, out_name: str):
    os.system('makeblastdb -in ' + in_name + ' -parse_seqids -blastdb_version 5 -title '+ '"' + out_name + '" -dbtype nucl -out ' + out_name)
    
    
def generate_blast(db: str, query: str, out_name: str):
    os.system('blastn -task blastn -db ' + db + ' -query ' + query + ' -out ' + out_name + '.csv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qlen slen qseq sseq" -num_threads ' + str(os.cpu_count()) + ' -word_size 7 -evalue 0.02')
    
    
def main():
    parser = argparse.ArgumentParser('One Step Blast Generator')
    parser.add_argument('--dbseqin', help='Input File Name of Seqs to Generate Blast DB', type=str, required='--db' not in sys.argv)
    parser.add_argument('--db', help='User specify DB', type=str, required='--dbseqin' not in sys.argv)
    parser.add_argument('--out', help='Output name for Blast result.', type=str, required=True)
    parser.add_argument('--query', help='Query Seq File Name', type=str, required=True)
    
    args = parser.parse_args()
    
    db_name = args.db
    
    if args.db is None and args.dbseqin is not None:
        db_name = args.dbseqin + '_DB'
        generate_db(args.dbseqin, db_name)
    
    generate_blast(db_name, args.query, args.out)
    
    
if __name__ == "__main__":
    main()