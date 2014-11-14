cimport sRNAMapperC

cdef test1():
    cdef sRNAMapperC.StringHashTable * read_counts
    num_bins = 1000000000
    read_counts = sRNAMapperC.create_hash_table_string(num_bins)
    