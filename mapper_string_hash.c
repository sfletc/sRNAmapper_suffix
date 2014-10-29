/*
 * A script for reading a series of sRNAs read counts and
 * mapping exact matches to positions on a reference sequence.
 * This implementation utilises a standard string hashing function
 * and string hash table to map the reads.
 */

#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fasta.h" //Header file from fasta library
#include <time.h>

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

//Analysis variables
clock_t read_time = 0;
clock_t write_time = 0;
clock_t scan_time = 0;
clock_t build_hashtable_time = 0;
clock_t read_sRNAs_time = 0;
int num_collisions = 0;
long num_sRNA_reads = 0;
/**
 * Initialise a hash table that maps sequence strings to integers.
 *
 * The function allocates space for the struct and the linked lists
 * in the bins
 */
StringHashTable *create_hash_table_string(uint32_t num_bins) {
	if (num_bins < 1) {
		printf("%d is an invalid number of bins\n", num_bins);
		return NULL ;
	}
	StringHashTable *new_table;
	StringLinkedList **new_bins;
	//Allocate memory for the table structure
	new_table = malloc(sizeof(StringHashTable));
	//Allocate memory for the bins
	new_bins = calloc(num_bins, sizeof(StringLinkedList*));
	//There were errors allocating memory
	assert(new_table != NULL && new_bins != NULL);
	new_table->bins = new_bins;
	new_table->num_bins = num_bins;
	return new_table;
}

/**
 * This function adds an sRNA sequence to the hash table. If there is nothing in the bin it creates
 * a new linked list element and adds it to the bin. If there is something in the bin it
 * continues along the list and if there is a match the read count is added to the match,
 * if there is no match the new sequence is added to the end of the list.
 */
void put_string(StringHashTable *hash_table, char *string, short string_length,
		uint32_t num_reads) {
	StringLinkedList *existing_list;
	uint32_t bin_num;
	//Find the bin from the hash (which is the 2bit sequence itself)
	uint64_t hash = hash_string(string);
	bin_num = hash % hash_table->num_bins;
	existing_list = hash_table->bins[bin_num];
	//If there is nothing in this bin create a new linked list with the sequence and put it in the bin
	if (existing_list == NULL ) {
		StringLinkedList *new_element = malloc(sizeof(*new_element));
		//Something wrong with memory allocation
		assert(new_element != NULL);
		new_element->num_reads = num_reads;
		//Copy the string and store it in the list element
		new_element->string = malloc(string_length * sizeof(char) + 1);
		new_element->string[string_length] = '\0';
		memcpy(new_element->string, string, string_length);
		new_element->next = NULL;
		hash_table->bins[bin_num] = new_element;
	}
	//There is something in the bin. Iterate through the list. If this sequence is already present then don't put it in again.
	else {
		num_collisions += 1;
		iterate_and_add_from_string(existing_list, string, string_length,
				num_reads);
	}
}

/* Access the value of an element in the hash table.
 *
 */
uint32_t get_counts_from_string(StringHashTable *hash_table, char *sequence,
		short window_size) {
	uint32_t bin_num;
	StringLinkedList *existing_list;
	StringLinkedList *listIt;
	uint64_t hash = hash_string(sequence);
	bin_num = hash % hash_table->num_bins;
	existing_list = hash_table->bins[bin_num];
	//If there is no linked list structure then the sequence hasn't been entered
	if (existing_list != NULL ) {
		//Loop through the linked list and attempt to match the sequence
		for (listIt = existing_list; listIt != NULL ; listIt = listIt->next) {
			//If the sequence matches then return the stored value
			if (strcmp(listIt->string, sequence) == 0) {
				return listIt->num_reads;
			}
		}
	}
	//Otherwise there is no match - return -1
	return -1;
}

/**Iterate through the linked list and if a match is found add the count of
 * the new element to the match. If the match isn't found add the element to
 * the end of the list
 */
void iterate_and_add_from_string(StringLinkedList *existing_list, char *string,
		short string_length, uint32_t num_reads) {
	StringLinkedList *listIt;
	for (listIt = existing_list; listIt->next != NULL ; listIt = listIt->next) {
		//If the sequence is already present add this count to the previous count and return
		if (strcmp(listIt->string, string)) {
			listIt->num_reads += num_reads;
			return;
			//This isn't a duplicate - move to the next element
		}
	}
	//Check for the sequence in the final node
	//If the sequence is already present add this count to the previous count and return
	if (strcmp(listIt->string, string)) {
		listIt->num_reads += num_reads;
		return;
	}
	//The sequence isn't present. Create a new node and add it to the list
	StringLinkedList *new_element = malloc(sizeof(StringLinkedList));
	//Copy the string into the hash table
	new_element->string = malloc(string_length * sizeof(char) + 1);
	new_element->string[string_length] = '\0';
	memcpy(new_element->string, string, string_length);
	new_element->num_reads = num_reads;
	new_element->next = NULL;
	listIt->next = new_element;
}
/**
 * Calculate the hash value for a string.
 *
 * This is an implementation of the djb2 string hash function
 */
uint64_t hash_string(char *string) {
	{
		uint64_t hash = 5381;
		int c;
		while ((c = (unsigned char) *string++) != 0)
			hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
		return hash;
	}
}

/**
 * Get the reverse complement of a sequence and return in in output.
 *
 * This could be changed to recognise the null terminator instead of requiring
 * the nummber of characters.
 */
void revcmp_seq(char *input, char **output, int num_chars) {
	char *temp;
	int i;
	i = 0;
	temp = *output;
	//Move the input pointer to the end of the char array
	input = input + (num_chars - 1);
	//We've already reached the first position
	while (i < num_chars) {
		switch (*input) {
		case ('A'):
			*temp = 'T';
			break;
		case ('T'):
			*temp = 'A';
			break;
		case ('G'):
			*temp = 'C';
			break;
		case ('C'):
			*temp = 'G';
			break;
		}
		//Dereference both pointers and assign output to input
		input--;
		temp++;
		i++;
	}
	//Return to the original index
	*output = temp - num_chars;
}
;

/** Map the sRNA reads on to a reference sequence and return a
 * forward and reverse array with the number of reads mapped to
 * each position.
 *
 * The return arrays are of length (sequence_length - window_size + 1)
 * which is the number of available positions for the read to map. The function
 * iterates through the reference sequence, and checks in the given sRNA hash
 * table if the window is present. If it is present then the reads are mapped
 * to the smallest index of the window.
 *
 */
ReadCounts *get_sequence_counts_reverse_window_string_hash(
		StringHashTable *sRNA_count_hashtable, char* sequence,
		uint32_t sequence_length, short window_size) {
	//Slice all of the characters in the sequence
	char *window, *rev_complement;
	uint32_t *fwd_counts, *rev_counts;
	uint32_t count;
	ReadCounts *read_counts;
	read_counts = malloc(sizeof(ReadCounts));
	//Get this correct
	int num_positions = (sequence_length - window_size + 1);
	fwd_counts = calloc(num_positions, sizeof(uint32_t));
	rev_counts = calloc(num_positions, sizeof(uint32_t));
	window = malloc(sizeof(char) * window_size + 1);
	window[window_size] = '\0';
	rev_complement = malloc(sizeof(char) * window_size + 1);
	rev_complement[window_size] = '\0';
	int i;
	int reset = 0; //The number of bases that need to be pushed through to flush out the ambiguous base(s) in
	// the window
	for (i = 0; i < sequence_length - window_size + 1; ++i) {
		//If the next character is an ambiguous character make sure the hashing is skipped
		if (sequence[window_size - 1] == 'N') {
			reset = window_size;
		}
		if (reset != 0) {
			--reset;
			++sequence;
			continue;
		}
		//Get the window
		//TODO This copying may not be necessary - we are just getting a hash.
		memcpy(window, sequence, window_size);
//		printf("%s\n", window);
		++sequence;
		count = get_counts_from_string(sRNA_count_hashtable, window,
				window_size);
		if (count != -1) {
			fwd_counts[i] += count;
		}
		//Get the reverse compliment sequence
		revcmp_seq(window, &rev_complement, window_size);
		count = get_counts_from_string(sRNA_count_hashtable, rev_complement,
				window_size);
		if (count != -1) {
			rev_counts[i] += count;
		}
	}
	//Assign the
	read_counts->fwd_counts = fwd_counts;
	read_counts->rev_counts = rev_counts;
	//Free the sequence window
	free(window);
	free(rev_complement);
	return read_counts;
}
;


StringHashTable *read_sRNAs_string(char *filename, short required_read_length,
		uint32_t num_bins) {
	FILE *file;
	char *sequence; //Max read length + 1 for null terminator
	uint32_t num_reads_32;
	int num_reads;
	short read_length;
	clock_t begin;
	uint64_t max_bins_required;
	StringHashTable *seq_read_nums;
	//Assign space for the sequences to be read into (with null terminator)
	sequence = malloc(sizeof(char) * (required_read_length + 1));
	//Don't waste space - this is a perfect hash function up to 4^width-1 bins
	//The bit shift is the equivalent to 4 ^ required_read_length
	max_bins_required = 1;
	max_bins_required <<= (required_read_length * 2);
	if (max_bins_required <= num_bins) {
		printf("The perfect hash function only requires %d bins\n",
				(int) max_bins_required);
		num_bins = max_bins_required;
	}
	begin = clock();
	seq_read_nums = create_hash_table_string(num_bins);
	build_hashtable_time += clock() - begin;
	begin = clock();
	//Open the file
	file = fopen(filename, "r");
	if (file == NULL ) {
		printf("Error reading file %s", filename);
		exit(1);
	}
	while (!feof(file)) {
		++num_sRNA_reads;
		//Read in the sequence width and RNA seq read count
		if (fscanf(file, "%hd\t%d\t", &read_length, &num_reads) != 2) {
			printf(
					"Invalid format (should be short\\tint\\tstring) at line %d of %s\nUsing only prior reads\n",
					(int) num_sRNA_reads, filename);
			break;
		}
		//If the read length isn't correct then move to the next line
		if (read_length != required_read_length) {
			fscanf(file, "%*[^\n]");
			continue;
		} else {
			//Read the sequence until the next space or newline (just in case there is trailing whitespace)
			if (fscanf(file, "%[^ \n]", sequence) != 1) {
				printf(
						"Invalid format (should be short\\tint\\tstring) at line %d of %s\nUsing only prior reads\n",
						(int) num_sRNA_reads, filename);
				break;
			}
		}
		num_reads_32 = (uint32_t) num_reads;
		put_string(seq_read_nums, sequence, read_length, num_reads_32);
	}
	read_sRNAs_time = clock() - begin;
	//Deallocate memory used by sequence
	free(sequence);
	//Free the file
	fclose(file);
	return seq_read_nums;
}
;

void calc_and_write_output_string(StringHashTable *sRNA_read_counts,
		char *input_filename, char* output_filename, short window_size) {
	char *seq, *seq_window;
	char *name;
	clock_t begin;
	uint32_t count;
	int length, i;
	FILE *output_file;
	//Open the output file
	output_file = fopen(output_filename, "w");
	if (output_file == NULL ) {
		//Error reading the file
		printf("Error reading file %s: %d (%s)\n", output_filename, errno,
				strerror(errno));
	}
	FASTAFILE *fasta_file;
	ReadCounts *read_counts;
	//Assign memory for the sequence window
	seq_window = malloc(sizeof(char) * (window_size + 1));
	//Null terminate the string
	seq_window[window_size] = '\0';
	fasta_file = OpenFASTA(input_filename);
	begin = clock();
	while (ReadFASTA(fasta_file, &seq, &name, &length) != 0) {
		read_time += clock() - begin;
		//Sequences less than the window size won't be bound to
		if (length < window_size) {
			continue;
		}
		begin = clock();
//		read_counts = get_sequence_counts_direct_hash_rev_cmp(sRNA_read_counts,
//				seq, length, window_size);
		read_counts = get_sequence_counts_reverse_window_string_hash(
				sRNA_read_counts, seq, length, window_size);
//		read_counts - get_sequence_counts_reverse_window_rev_cmp()
		scan_time += (clock() - begin);
		if (read_counts == NULL ) {
			//Invalid char in sequence - complete the error message
			printf(" %s in file %s\n", name, input_filename);
			exit(1);
		}
		begin = clock();
		for (i = 0; i < (length - window_size + 1); ++i) {
			//Format: Ref_seq_ID', 'sRNA_seq', 'count', "pos_from_5'_fwd", Fwd_or_Rvs
			count = read_counts->fwd_counts[i];
			if (count != 0) {
				//Make a copy of the sequence at the current position
				memcpy(seq_window, &seq[i], window_size);
				if (fprintf(output_file, "%s,%s,%d,%d,%s\n", name, seq_window,
						count, i, "f") < 0) {
					printf("There was an error writing the output file\n");
					exit(1);
				}
			}
			count = read_counts->rev_counts[i];
			if (count != 0) {
				//Flip the sequence to match the output
				revcmp_seq(&seq[i], &seq_window, window_size);
				//TODO - the i+1 is to keep the file the same as the python version
				if (fprintf(output_file, "%s,%s,%d,%d,%s\n", name, seq_window,
						-count, i + window_size, "r") < 0) {
					printf("There was an error writing the output file\n");
					exit(1);
				}
			}
		}
		write_time += clock() - begin;
		begin = clock();
		//Free the read files
		free(seq);
		free(name);
		//Free the memory allocated in reading the sequence
		free(read_counts->fwd_counts);
		free(read_counts->rev_counts);
		free(read_counts);
		read_counts = NULL;
	}
	//Free the memory allocated for the struct
	free(seq_window);
	free(read_counts);
	//Free the input and output files
	CloseFASTA(fasta_file);
	fclose(output_file);
}
;

int main(int argc, char *argv[]) {
//##########################Test reading of sequences########################

	//TODO investigate why bin sizes of ~100,000,000 lead to a dramatic slow down
	//even though there are less hash collisions and its not allocation of space
	//for the hash table
	uint32_t num_bins;

	char *sRNA_filename, *reference_filename, *output_filename;
	short width;
	clock_t begin;
	printf("SRNA mapper: ");
	if (argc != 6) {
		printf("Five arguments required");
		return -1;
	}
//	FILE *out = fopen("Bin_size_ex_time", "w");
	sRNA_filename = argv[1];
	reference_filename = argv[2];
	output_filename = argv[5];
	num_bins = strtol(argv[4], NULL, 10);
	width = (short) strtol(argv[3], NULL, 10);
	StringHashTable *seq_read_counts;
	seq_read_counts = read_sRNAs_string(sRNA_filename, width, num_bins);
//	HashTable *seq_read_counts;
//	seq_read_counts = read_sRNAs(sRNA_filename, width, num_bins);
	printf("Time spent building hashtable: %2.2fs\n",
			((double) build_hashtable_time) / CLOCKS_PER_SEC);
	printf("Time spent building hashtable: %2.2fs\n",
			((double) build_hashtable_time) / CLOCKS_PER_SEC);
	printf("Time taken reading sRNAs: %2.2fs\n",
			((double) read_sRNAs_time) / CLOCKS_PER_SEC);
	begin = clock();
	calc_and_write_output_string(seq_read_counts, reference_filename,
			output_filename, width);
	printf("Number of hash collisions: %d from %d sequences - %0.2f%%\n",
			(int) num_collisions, (int) num_sRNA_reads,
			((float) num_collisions) / ((float) num_sRNA_reads) * 100);
	printf("Time taken scanning and matching to sequences: %fs\n",
			((double) (clock() - begin)) / CLOCKS_PER_SEC);
	printf("Time spent reading in sequences: %f2.2s\n",
			((double) read_time) / CLOCKS_PER_SEC);
	printf("Time spent scanning sequences for matches: %f2.2s\n",
			((double) scan_time) / CLOCKS_PER_SEC);
	printf("Time spent writing output file: %2.2fs\n",
			((double) write_time) / CLOCKS_PER_SEC);
	//Different bin sizes characterization
	//	int i;
//	for (i = 0; i < 500000000; i+=10000) {
//		sRNA_filename = argv[1];
//		reference_filename = argv[2];
//		output_filename = argv[5];
//		num_bins = strtol(argv[4], NULL, 10);
//		width = (short) strtol(argv[3], NULL, 10);
//		seq_read_counts = read_sRNAs(sRNA_filename, width, num_bins);
//		printf("Time spent building hashtable: %2.2fs\n",((double) build_hashtable_time)/ CLOCKS_PER_SEC);
//		fprintf("Time spent building hashtable: %2.2fs\n",((double) build_hashtable_time)/ CLOCKS_PER_SEC);
//		printf("Time taken reading sRNAs: %2.2fs\n", ((double) read_sRNAs_time)/ CLOCKS_PER_SEC);
//		begin = clock();
//		calc_and_write_output(seq_read_counts,reference_filename, output_filename, width);
//		printf("Number of hash collisions: %d from %d sequences - %0.2f%%\n", num_collisions, num_sRNA_reads, ((float) num_collisions) / ((float) num_sRNA_reads) * 100);
//		printf("Time taken scanning and matching to sequences: %fs\n", ((double)(clock()-begin))/ CLOCKS_PER_SEC);
//		printf("Time spent reading in sequences: %f2.2s\n",((double) read_time)/ CLOCKS_PER_SEC);
//		printf("Time spent scanning sequences for matches: %f2.2s\n",((double) scan_time)/ CLOCKS_PER_SEC);
//		printf("Time spent writing output file: %2.2fs\n",((double) write_time)/ CLOCKS_PER_SEC);
//	}

//
//	if (failed){
//		printf("Five arguments required: sRNA_file, reference_file, read_width, num_bins(the number of bins for the hash table), csv_output_file");
//		return -1;
//	}

//######################Test getting reverse complement sequence##############
//	char *sequence = "GGCAT";
//	char *revcmp = malloc(sizeof(char)*6);
//	revcmp[5]= '\0';
//	char *correct = "ATGCC";
//	revcmp_seq(sequence, 5, &revcmp);
//	assert(strcmp(revcmp,correct));
//	sequence = "A";
//	correct = "T";
//	revcmp[1]= '\0';
//	revcmp_seq(sequence, 1, &revcmp);
//	assert(strcmp(revcmp,correct));
//###########################Test the sRNA hashing############################
//	uint32_t num_bins;
//	HashTable *seq_read_counts;
//	short width;
//	num_bins = 100000;
//	width = 5;
//	seq_read_counts = read_sRNAs("sRNA_test_1.txt", width, num_bins);
//	//Test that the sRNAs are read properly
//	assert(get(seq_read_counts, hash_seq("ATGGT", 5)) == 20);
//	assert(get(seq_read_counts, hash_seq("ACCAT", 5)) == 27);
//	assert(get(seq_read_counts, hash_seq("CTCGT", 5)) == 28);
//	assert(get(seq_read_counts, hash_seq("TGCTC", 5)) == 25);
//	assert(get(seq_read_counts, hash_seq("ATGGA", 5)) == 22);
//	assert(get(seq_read_counts, hash_seq("AAAAA", 5)) == 27);
//	assert(get(seq_read_counts, hash_seq("TTTTT", 5)) == 200);
//	assert(get(seq_read_counts, hash_seq("TCCAT", 5)) == 111);
//	assert(get(seq_read_counts, hash_seq("ACGTA", 5)) == 345);
//###########################Test building of hash table######################
//	char *filename;
//	short width;
//	uint32_t num_bins;
//	HashTable *read_counts;
//	num_bins = 100000;
//	width = 5;
//	read_counts = read_sRNAs(filename, width, num_bins, true);
//	filename = "sRNA_test_1.txt";
	//#######################Test reverse bits################################
//	uint64_t bits, rev_bits;
//	bits = 12514235;
//	rev_bits = reverse_bits_fwd(bits);
//	assert(rev_bits == 0xBEF3BB);
//	assert(rev_bits == reverse_bits_fwd(bits));
//  ########################Test reverse char array###########################
//	char *seq, *correct;
//	char *output;
//	output = malloc(sizeof(char) * (7 + 1));
//	output[7] = '\0';
//	seq = malloc(sizeof(char) * (7 + 1));
//	seq[7] = '\0';
//	correct = malloc(sizeof(char) * (7 + 1));
//	correct[7] = '\0';
//	seq = "ABCDEFG";
//	reverse_char_array(seq, &output, 7);
//	reverse_char_array(output, &correct, 7);
//	assert(strcmp(seq,correct));
//	printf("%s", output);

	//#############################Test the hash function########################
//	char *seq;
//	uint64_t hash_val;
//	seq = "ATGGTCGTAA";
//	hash_val = hash(seq, 10);
//	assert(hash_val == 586186);

//##############################Test the hash table###########################
//	HashTable *hashTable;
//	uint32_t num_bins;
//	unsigned short sRNA_width;
//	uint64_t sRNA1, sRNA2, sRNA3, sRNA4, sRNA5;
//	uint32_t count1, count2, count3, count4, count5, output;
//	char *testWindow1, *testWindow2;
//	//Test the hash table
//	sRNA_width = 4;
//	sRNA1 = 6534563;
//	sRNA2 = 6534563;
//	sRNA3 = 653454523;
//	num_bins = 10000;
//	sRNA4 = sRNA1 % num_bins;
//	sRNA5 = sRNA1 % num_bins;
//	count1 = 1000;
//	count2 = 21222;
//	count3 = 58;
//	count4 = 6534;
//	count5 = 145243;
//	testWindow1 = malloc(sizeof(char)*100);
//	testWindow2 = malloc(sizeof(char)*100);
//	hashTable = create_hash_table(num_bins);
//	//Add a value
//	put(hashTable, sRNA1, count1);
//	//Retrieve a value
//	output = get(hashTable, sRNA1);
//	assert(output == count1);
//	//Add a duplicate value in the first position
//	put(hashTable, sRNA2, count2);
//	//Make sure that the counts were added
//	output = get(hashTable, sRNA1);
//	assert(output == count1 + count2);
//	//Add a new value
//	put(hashTable, sRNA3, count3);
//	//Add a new value with the same hash
//	put(hashTable, sRNA4, count4);
//	//Add a duplicate value with the same hash
//	output = get(hashTable, sRNA4);
//	assert(output = count4);
//	put(hashTable, sRNA5, count5);
//	//Check that the outputs have been added
//	output = get(hashTable, sRNA5);
//	assert(output = count5+count4);
	return 1;
}
;

