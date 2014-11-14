//A basic linked list format with the list stored as a perfect hash
typedef struct _LinkedList_ {
	uint32_t num_reads; //The number of reads from RNA-seq of this exact sequence
	uint64_t string; //The sequence in two-bit format
	struct _LinkedList_ *next;
} LinkedList;
//A basic linked list format with the sequence stored as a character array

typedef struct _HashTable_ {
	uint32_t num_bins; /* the size of the table */
	LinkedList **bins; //This strictly doesn't need to
} HashTable;
