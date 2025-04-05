# Geptop3
A newly essential gene predict tool combining homology mapping and machine learning. Geptop3 is available online at (http://gepa.org.cn/Geptop3).

## Environment
1. Linux or Windows environment.
2. Python(>=3.6) and compatible Biopython, numpy, pandas, Scikit-learn, imbanlanced-learn.
3. Standalone NCBI BLAST.

## Usage
1. Whole-genomic __NUCLEOTIDE__ sequences in __FASTA__ format is the only usable type.
2. Change the root_dir variable in both file_decompress.py and Geptop3.py to the path of the Geptop3 directory on your device. Then run file_decompress.py to decompress the files.
3. Save the FASTA files in the uploadFile folder, and then run the terminal command. *python geptop3.py -i uploadFile/.fna*. Optional parameters: *–s essentiality score cutoff*, range:0~1, default: 0.28.
