def read_seqs(seq_file_name: str):
    seq_file = open(seq_file_name, mode='r')
    lines = seq_file.readlines()

    i = 0
    seqs = dict()
    while i <= len(lines) - 2:
        name_line = lines[i]
        if ' ' in name_line:
            seq_id = name_line[1:name_line.index(' ')]
        else:
            seq_id = name_line[1:].strip()
        seq = lines[i+1][:-1]
        seqs[seq_id] = seq
        i += 2

    return seqs


