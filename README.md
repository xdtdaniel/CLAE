# SequenceExtraction

This is a complete toolkit generating consensus sequence.

## Requirements
- numpy==1.17.2
- Bio==0.4.1
- pandas==1.2.4
- tqdm==4.60.0
- halo

## Dependencies


- Sparc: https://github.com/yechengxi/Sparc
- pbdagcon: https://github.com/PacificBiosciences/pbdagcon
- blasr: https://github.com/PacificBiosciences/blasr
- BLAT: http://genome.ucsc.edu/FAQ/FAQblat.html
- BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi
- muscle: https://www.ebi.ac.uk/Tools/msa/muscle/
- mugio: https://github.com/pspealman/mugio (needed if you want to use phred score to filter sequences)

|          | Reference Mode | | Non-reference Mode |||
| -------- | -------------- |--| ------------------ |--|--|
|          | pbdagcon       | Sparc              | pbdagcon | Sparc | muscle |
| pbdagcon | √              |                    | √ |  |  |
| Sparc    |                | √                  |  | √ |  |
| muscle   |                |                    |  |  | √ |
| blasr    | √              | √                  | √ | √ |  |
| BLAT     |                |                    |  |  | √ |
