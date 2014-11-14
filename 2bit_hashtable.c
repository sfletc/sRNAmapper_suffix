#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "2bit_hashtable_def.h"
#include "exit_codes.h"


//######################Hash table functions########################################
/**The creation function for the hash table*/
HashTable *create_hash_table_2bit(int num_bins) {
	if (num_bins < 1) {
		printf("%d is an invalid number of bins\n", num_bins);
		return NULL ;
	}
	HashTable *new_table;
	LinkedList **new_bins;
	//Allocate memory for the table structure
	new_table = malloc(sizeof(HashTable));
	//Allocate memory for the bins
	new_bins = calloc(num_bins, sizeof(LinkedList*));
	//There were errors allocating memory
	if (new_table == NULL || new_bins == NULL){
		printf("%s", "Out of memory");
		exit(EXIT_OUT_OF_MEM);
	}
	new_table->bins = new_bins;
	new_table->num_bins = num_bins;
	return new_table;
}

/**Get the hash of the complement to the hash for a sequence.
 Because of the way we've structured the bit equivalents of
 the characters the bitwise negation gives the complementary
 Sequence with leading ones which then need to be pushed off the hash*/
void put_2bit(HashTable *hash_table, uint64_t sequence, uint32_t num_reads) {
	LinkedList *existing_list;
	uint32_t bin_num;
	//Find the bin from the hash (which is the 2bit sequence itself)
	bin_num = sequence % hash_table->num_bins;
	existing_list = hash_table->bins[bin_num];
	//If there is nothing in this bin create a new linked list with the sequence and put it in the bin
	if (existing_list == NULL ) {
		LinkedList *new_element = malloc(sizeof(*new_element));
		//Something wrong with memory allocation
		assert(new_element != NULL);
		new_element->num_reads = num_reads;
		new_element->string = sequence;
		new_element->next = NULL;
		hash_table->bins[bin_num] = new_element;
	}
	//There is something in the bin. Iterate through the list. If this sequence is already present then don't put it in again.
	else {
		num_collisions += 1;
		iterate_and_add_2bit(existing_list, sequence, num_reads);
	}
}

/**
 *
 */
uint32_t get_2bit(HashTable *hash_table, uint64_t sequence) {
	uint32_t bin_num;
	LinkedList *existing_list;
	LinkedList *listIt;
	bin_num = sequence % hash_table->num_bins;
	existing_list = hash_table->bins[bin_num];
	//If there is no linked list structure then the sequence hasn't been entered
	if (existing_list != NULL ) {
		//Loop through the linked list and attempt to match the sequence
		for (listIt = existing_list; listIt != NULL ; listIt = listIt->next) {
			//If the sequence matches then return the stored value
			if (listIt->string == sequence) {
				return listIt->num_reads;
			}
		}
	}
	return -1;
}

/**Iterate through the linked list and if a match is found add the count of
 * the new element to the match. If the match isn't found add the element to
 * the end of the list
 */
void iterate_and_add_2bit(LinkedList *existing_list, uint64_t sequence,
		uint32_t num_reads) {
	LinkedList *listIt;
	for (listIt = existing_list; listIt->next != NULL ; listIt = listIt->next) {
		//If the sequence is already present add this count to the previous count and return
		if (listIt->string == sequence) {
			listIt->num_reads += num_reads;
			return;
			//This isn't a duplicate - move to the next element
		}
	}
	//Check for the sequence in the final node
	//If the sequence is already present add this count to the previous count and return
	if (listIt->string == sequence) {
		listIt->num_reads += num_reads;
		return;
	}
	//The sequence isn't present. Create a new node and add it to the list
	LinkedList *new_element = malloc(sizeof(LinkedList));
	new_element->string = sequence;
	new_element->num_reads = num_reads;
	new_element->next = NULL;
	listIt->next = new_element;
}
