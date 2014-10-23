Just run on the command line with exactly 5 arguments:
sRNAMapperC sRNA_read_filename reference_seqs_filename reference_width number_of_hashtable_bins output_filename
I found 16,777,216 to be a good number of bins for the reference files given - bonus points if it is a multiple of 2
eg command:
sRNAMapper_suffix.exe sRNA_sample_1.txt multi_sequence_reference.fa 24 16777216 test_out.csv
