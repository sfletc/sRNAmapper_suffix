/* Define a generic linked list and hash table structs for the hash table
 */
typedef struct _StringLinkedList_ {
	uint32_t num_reads;
	char* string;
	struct _StringLinkedList_ *next;
} StringLinkedList;

typedef struct _StringHashTable_ {
	uint32_t num_bins; /* the size of the table */
	StringLinkedList **bins; //This strictly doesn't need to
} StringHashTable;

/*Declare functions here */
StringHashTable *create_hash_table_string(uint32_t num_bins);
void put_string(StringHashTable *hash_table, char *string, short string_length,
		uint32_t num_reads);
uint32_t get_counts_from_string(StringHashTable *hash_table, char *sequence,
		short window_size);
void iterate_and_add_from_string(StringLinkedList *existing_list, char *string,
		short string_length, uint32_t num_reads);
uint64_t hash_string(char *string);
void revcmp_seq(char *input, char **output, int num_chars);
