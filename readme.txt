Just run on the command line with exactly 5 arguments:
sRNAMapperC sRNA_read_filename reference_seqs_filename reference_width number_of_hashtable_bins output_filename
I found 10,000,000 to be a good number of bins for the reference files given - bonus points if it is a multiple of 2
eg command:
sRNAMapperC sRNA_sample_1.txt multi_sequence_reference.fa 24 10000000 test_out.csv
