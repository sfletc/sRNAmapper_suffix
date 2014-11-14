cimport stdintp as si

#Maybe this can be improved change the name to import just the required types
#but this will do for now
cdef extern from "mapper_string.h":
    ctypedef struct StringHashTable:
        pass
    ctypedef struct StringLinkedList:
        pass
    StringHashTable *create_hash_table_string(si.uint32_t num_bins);
    void put_string(StringHashTable *hash_table, char *string, short string_length,
        si.uint32_t num_reads);
    si.uint32_t get_counts_from_string(StringHashTable *hash_table, char *sequence,
        short window_size);
    void iterate_and_add_from_string(StringLinkedList *existing_list, char *string,
        short string_length, si.uint32_t num_reads);
    si.uint64_t hash_string(char *string);
    void revcmp_seq(char *input, char **output, int num_chars);