//     Copyright 2013 Jason Jooste Licensed under the
//     Educational Community License, Version 2.0 (the "License"); you may
//     not use this file except in compliance with the License. You may
//     obtain a copy of the License at
//
//     http://www.osedu.org/licenses/ECL-2.0
//
//     Unless required by applicable law or agreed to in writing,
//     software distributed under the License is distributed on an "AS IS"
//     BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
//     or implied. See the License for the specific language governing
//     permissions and limitations under the Lircense.
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fasta.h" //Header file from fasta library
#include <time.h>
//A basic linked list format
typedef struct _LinkedList_ {
	uint32_t num_reads; //The number of reads from RNA-seq of this exact sequence
	uint64_t sequence; //The sequence in two-bit format
	struct _LinkedList_ *next;
} LinkedList;

typedef struct _HashTable_ {
	uint32_t num_bins; /* the size of the table */
	LinkedList **bins; //This strictly doesn't need to
} HashTable;

typedef struct _ReadCounts_ {
	uint32_t *fwd_counts;
	uint32_t *rev_counts;
} ReadCounts;
//Define boolean operators
typedef enum {
	false, true
} bool;

HashTable *create_hash_table(uint32_t);
void iterate_and_add(LinkedList*, uint64_t, uint32_t);
void iterate_and_stop(LinkedList*, uint64_t, uint32_t);
void put(HashTable*, uint64_t, uint32_t);
uint32_t get(HashTable*, uint64_t);
uint64_t hash_seq(char*, short);
void reverse_char_array(char*, char**, int);
uint64_t reverse_bits(uint64_t);

//Timing variables
clock_t read_time = 0;
clock_t write_time = 0;
clock_t scan_time = 0;
clock_t build_hashtable_time = 0;
clock_t read_sRNAs_time = 0;
int num_collisions = 0;
long num_sRNA_reads = 0;


//#############Begin borrowed code for fasta file reading#########################
/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id: fasta.c,v 1.1 2003/10/05 18:43:39 eddy Exp $
 */
/* Function: OpenFASTA(), ReadFASTA(), CloseFASTA().
 * Date:     SRE, Sun Sep  8 06:39:26 2002 [AA2721, transatlantic]
 *
 * Purpose:  A very rudimentary FASTA file reading API. Designed
 *           for simplicity and clarity, not for robustness.
 *
 *           The API is:
 *
 *           FASTAFILE *ffp;
 *           char      *seq;
 *           char      *name;
 *           int        seqlen;
 *
 *           ffp = OpenFASTA(seqfile);
 *           while (ReadFASTA(ffp, &seq, &name, &seqlen)
 *           {
 *             do stuff with sequence;
 *             free(name);
 *             free(seq);
 *           }
 *           CloseFASTA(ffp);
 *
 * Args:
 *           seqfile   - name of a FASTA file to open.
 *           seq       - RETURN: one sequence
 *           name      - RETURN: name of the sequence
 *           seqlen    - RETURN: length of the sequence in residues
 *           ffp       - ptr to a FASTAFILE object.
 *
 * Commentary:
 *           The basic problem with reading FASTA files is that there is
 *           no end-of-record indicator. When you're reading sequence n,
 *           you don't know you're done until you've read the header line
 *           for sequence n+1, which you won't parse 'til later (when
 *           you're reading in the sequence n+1). One common trick for
 *           this is to implement a one-line "lookahead" buffer that you
 *           can peek at, before parsing later.
 *
 *           This buffer is kept in a small structure (a FASTAFILE), rather
 *           than in a static char[] in the function. This allows
 *           us to have multiple FASTA files open at once. The static approach
 *           would only allow us to have one file open at a time. ANSI C
 *           predates the widespread use of parallel programming. It was
 *           not overly concerned about the drawbacks of statics. Today,
 *           though, you should keep in mind that you may someday want to
 *           turn your program into a multithreaded, parallel program, and
 *           all functions in parallelized code must be "reentrant": able to
 *           be called a second time - with different arguments,
 *           and while the code in the first function call is still executing! -
 *           without overwriting or corrupting any static storage in the
 *           function. Statics have fewer uses now (for example, to
 *           test that some initialization code for a function is run once
 *           and only once.)
 *
 * Limitations:
 *           There is no error handling, for clarity's sake. Also,
 *           the parser is brittle. Improper FASTA files (for instance,
 *           blank lines between records) will cause unexpected
 *           behavior. Real file parsers are more complex.
 *           In real life, they have to deal with absolutely anything the user might
 *           pass as a "FASTA file"; and either parse it correctly,
 *           or detect that it's an invalid format and fail cleanly.
 *
 *           Lines are read in from the file using ANSI C's fgets(). fgets()
 *           requires a maximum buffer length (here, FASTA_MAXLINE, which is
 *           defined as 512 in bio5495.h). Some FASTA files have very long
 *           description lines, however; notably the NCBI NR database. Static
 *           limitations on things like line or sequence lengths should be
 *           avoided. An example of a replacement for fgets() that dynamically
 *           allocates its buffer size and allows any line length is
 *           SQUID's sre_fgets().
 *
 *           We use ANSI C's strtok() to parse the sequence name out of the line.
 *           strtok() is deprecated in modern programs because it is not threadsafe.
 *           (See comments above.) An example of a threadsafe version is
 *           SQUID's sre_strtok().
 *
 * Returns:
 *           OpenFASTA() returns a FASTAFILE pointer, or NULL on failure (for
 *           instance, if the file doesn't exist, or isn't readable).
 *
 *           ReadFASTA() returns 1 on success, or a 0 if there are no
 *           more sequences to read in the file.
 *
 *           CloseFASTA() "always succeeds" and returns void.
 */
FASTAFILE *
OpenFASTA(char *seqfile)
{
  FASTAFILE *ffp;

  ffp = malloc(sizeof(FASTAFILE));
  ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
  if (ffp->fp == NULL) { free(ffp); return NULL; }
  if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL)
    { free(ffp); return NULL; }
  return ffp;
}

int
ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, int *ret_L)
{
  char *s;
  char *name;
  char *seq;
//  int *seq_length;
  int   n;
  int   nalloc;

  /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
     */
    if (ffp->buffer[0] != '>') return 0;

    /* Parse out the name: the first non-whitespace token after the >
     * Jason: Altered this to allow spaces and terminate at newline
     */
    s  = strtok(ffp->buffer+1, "\t\n");
    name = malloc(sizeof(char) * (strlen(s)+1));
    strcpy(name, s);

    /* Everything else 'til the next descline is the sequence.
     * Note the idiom for dynamic reallocation of seq as we
     * read more characters, so we don't have to assume a maximum
     * sequence length.
     */
    seq = malloc(sizeof(char) * 128);     /* allocate seq in blocks of 128 residues */
    nalloc = 128;
    n = 0;
    while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp))
      {
        if (ffp->buffer[0] == '>') break;	/* a-ha, we've reached the next descline */

        for (s = ffp->buffer; *s != '\0'; s++)
  	{
      //Only accept sequence characters
      //Jason: altered this so that it explicitly deals with the expected characters in a fasta file instead of using isalpha
	  switch (*s)
	  {
	  case('A'): case('T'): case('G'): case('C'): case('N'):
		seq[n] = *s;
		break;
	  case('K'): case('M'): case('R'): case('Y'): case('S'): case('W'): case('B'): case('V'): case('H'): case('D'): case('X'):
		//Replace the character with out canonical way of representing an ambiguous base - 'N'
	    seq[n] = 'N';
	    break;
	  case('\n'):
		continue;
	  default:
		  printf("Invalid character \"%c\" at position %d in sequence %s", *s, n, name);
		  return 0;
	  }
	  n++;
	  if (nalloc == n)	        /* are we out of room in seq? if so, expand */
	    {			        /* (remember, need space for the final '\0')*/
	      nalloc += 128;
	      seq = realloc(seq, sizeof(char) * nalloc);
	    }
	}
    }
  seq[n] = '\0';

  *ret_name = name;
  *ret_seq  = seq;
  *ret_L    = n;
  return 1;
}

void
CloseFASTA(FASTAFILE *ffp)
{
  fclose(ffp->fp);
  free(ffp);
}
//######################End borrowed code for fasta file reading####################

//The creation function for the hash table
HashTable *create_hash_table(uint32_t num_bins) {
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
	assert(new_table != NULL && new_bins != NULL);
	new_table->bins = new_bins;
	new_table->num_bins = num_bins;
	return new_table;
}
//Get the hash of the complement to the hash for a sequence.
//Because of the way we've structured the bit equivalents of
//the characters the bitwise negation gives the complementary
//Sequence with leading ones which then need to be pushed off the hash
void put(HashTable *hash_table, uint64_t sequence, uint32_t num_reads) {
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
		new_element->sequence = sequence;
		new_element->next = NULL;
		hash_table->bins[bin_num] = new_element;
	}
	//There is something in the bin. Iterate through the list. If this sequence is already present then don't put it in again.
	else {
		num_collisions += 1;
		iterate_and_add(existing_list, sequence, num_reads);
	}
}

uint32_t get(HashTable *hash_table, uint64_t sequence) {
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
			if (listIt->sequence == sequence) {
				return listIt->num_reads;
			}
		}
	}
	return -1;
}

uint64_t get_cmp_hash(uint64_t hash, short num_leading_bits) {
	int64_t revcmp_hash;
	//Invert the bits - the hash is set up for just the compliment - not the reverse
	revcmp_hash = ~hash;
	//Push the leading ones off the hash
	revcmp_hash = revcmp_hash << num_leading_bits;
	revcmp_hash = revcmp_hash >> num_leading_bits;
	return revcmp_hash;
}

//Gets the reverse complement of a hash given the number of leading zeros that should be present
uint64_t get_revcmp_hash(uint64_t hash, short num_leading_bits) {
	int64_t revcmp_hash;

	//Reverse the bit string
	revcmp_hash = reverse_bits(hash);
	//Invert the bits
	revcmp_hash = ~revcmp_hash;
	//Push the trailing ones off
	revcmp_hash = revcmp_hash >> num_leading_bits;
	return revcmp_hash;
}
//Iterate through the linked list and if a match is found add the count of the new element to the match. If the match isn't found add the element to the end of the list
void iterate_and_add(LinkedList *existing_list, uint64_t sequence,
		uint32_t num_reads) {
	LinkedList *listIt;
	for (listIt = existing_list; listIt->next != NULL ; listIt = listIt->next) {
		//If the sequence is already present add this count to the previous count and return
		if (listIt->sequence == sequence) {
			listIt->num_reads += num_reads;
			return;
			//This isn't a duplicate - move to the next element
		}
	}
	//Check for the sequence in the final node
	//If the sequence is already present add this count to the previous count and return
	if (listIt->sequence == sequence) {
		listIt->num_reads += num_reads;
		return;
	}
	//The sequence isn't present. Create a new node and add it to the list
	LinkedList *new_element = malloc(sizeof(LinkedList));
	new_element->sequence = sequence;
	new_element->num_reads = num_reads;
	new_element->next = NULL;
	listIt->next = new_element;
}
//Iterate through the list - if a match is found return, otherwise add the new element to the end of the list
//void iterate_and_stop(LinkedList *existing_list, uint64_t sequence,
//		uint32_t num_reads) {
//	LinkedList *listIt;
//	for (listIt = existing_list; listIt->next != NULL ; listIt = listIt->next) {
//		//If the sequence is already present the abort
//		if (listIt->sequence == sequence) {
//			return;
//		}
//	}
//	//Check for the sequence in the final node
//	//If the sequence is already present then abort
//	if (listIt->sequence == sequence) {
//		return;
//	}
//	//The sequence isn't present. Create a node and add it to the list
//	LinkedList *new_element = malloc(sizeof(LinkedList));
//	new_element->sequence = sequence;
//	new_element->num_reads = num_reads;
//	new_element->next = NULL;
//	listIt->next = new_element;
//}
//Accepts the sequence in two-bit format (which is used as the hash) and put
//Search the hash table for the given sequence in 2bit format. If the sequence isn't present -1 is returned
void reverse_char_array(char *input, char **output, int num_chars) {
	char *temp;
	int i;
	i = 0;
	temp = *output;
	//Move the input pointer to the end of the char array
	input = input + (num_chars - 1);
	//We've already reached the first position
	while (i < num_chars) {
		//Dereference both pointers and assign output to input
		*temp = *input;
		input--;
		temp++;
		i++;
	}
	//Return to the original index
	*output = temp - num_chars;
}

ReadCounts *get_sequence_counts_reverse_window_rev_cmp(HashTable *seq_read_nums,
		char* sequence, uint32_t sequence_length, short window_size) {
	//Slice all of the characters in the sequence
	char *window, *rev_complement;
	uint32_t *fwd_counts, *rev_counts;
	uint32_t count;
	uint64_t hash;
	ReadCounts *read_counts;
	read_counts = malloc(sizeof(ReadCounts));
	//Get this correct
	fwd_counts = calloc(sequence_length - window_size, sizeof(uint32_t));
	rev_counts = calloc(sequence_length - window_size, sizeof(uint32_t));
	window = malloc(sizeof(char) * window_size + 1);
	window[window_size] = '\0';
	rev_complement = malloc(sizeof(char) * window_size + 1);
	rev_complement[window_size] = '\0';
	int i;
	for (i = 0; i < sequence_length; ++i) {
		//Get the window
		memcpy(window, sequence, window_size);
		++sequence;
		hash = hash_seq(window, window_size);
		count = get(seq_read_nums, hash);
		if (count != -1) {
			fwd_counts[i + window_size / 2] += count;
		}
		//Get the reverse compliment hash
		reverse_char_array(window, &rev_complement, window_size);
		hash = hash_seq(rev_complement, window_size);
		count = get(seq_read_nums, hash);
		if (count != -1) {
			rev_counts[i + window_size / 2] += count;
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

//Canonical method of reversing bits - push left
uint64_t reverse_bits(uint64_t input) {
	int count = sizeof(input) * 8 - 1;
	uint64_t rev_input = input;
	input >>= 1;
	while (input) {
		rev_input <<= 1;
		rev_input |= input & 1;
		input >>= 1;
		count--;
	}
	rev_input <<= count;
	return rev_input;
}

//ReadCounts *get_sequence_counts_reverse_hash_rev_cmp(HashTable *seq_read_nums,
//		char* sequence, uint32_t sequence_length, short window_size) {
//	//Slice all of the characters in the sequence
//	char *window;
//	uint32_t *fwd_counts, *rev_counts;
//	uint32_t count;
//	uint64_t hash;
//	short num_leading_bits;
//	num_leading_bits = sizeof(uint64_t) - window_size;
//	ReadCounts *read_counts;
//	read_counts = malloc(sizeof(ReadCounts));
//	//Get this correct
//	fwd_counts = calloc(sequence_length - window_size, sizeof(uint32_t));
//	rev_counts = calloc(sequence_length - window_size, sizeof(uint32_t));
//	window = malloc(sizeof(char) * window_size);
//	int i;
//	for (i = 0; i < sequence_length; ++i) {
//		//Get the window
//		memcpy(window, sequence, window_size);
//		++sequence;
//		hash = hash_seq(window, window_size);
//		count = get(seq_read_nums, hash);
//		if (count != -1) {
//			fwd_counts[i + window_size / 2] += count;
//		}
//		//Get the reverse compliment hash
//		//Invert the bits
//		hash = ~hash;
//		//Reverse the bit string
//		hash = reverse_bits(hash);
//		//Push the trailing ones off
//		hash = hash >> num_leading_bits;
//		count = get(seq_read_nums, hash);
//		if (count != -1) {
//			rev_counts[i + window_size / 2] += count;
//		}
//	}
//	//Assign the
//	read_counts->fwd_counts = fwd_counts;
//	read_counts->rev_counts = rev_counts;
//	return read_counts;
//}
//;

ReadCounts *get_sequence_counts_direct_hash_rev_cmp(HashTable *seq_read_nums,
		char* sequence, int sequence_length, short window_size) {
	uint32_t *fwd_counts, *rev_counts;
	uint32_t count;
	uint64_t hash, revcmp_hash;
	short num_leading_bits;
	//Number of bits in the 64 bit integer - number of bases * 2 bits per base
	num_leading_bits = (sizeof(hash)*8) - (window_size*2);
	ReadCounts *read_counts;
	read_counts = malloc(sizeof(ReadCounts));
	int num_positions = (sequence_length - window_size+1);
	fwd_counts = calloc(num_positions, sizeof(uint32_t));
	rev_counts = calloc(num_positions, sizeof(uint32_t));
	if (fwd_counts == NULL || rev_counts == NULL){
		printf("Error allocating memory");
		exit(1);
	}
	int i; // The number of bases from the sequence iterated through
	hash = 0;
	revcmp_hash = 0;
	i = 0;
	while(i < sequence_length){
		//Build the hash with a full window from the next window of characters
		for (i = 0; i<(window_size-1); ++i) {
			switch (*sequence) {
			case 'T':
				revcmp_hash += 3LL<<(sizeof(revcmp_hash)*8-2); //11 at the first postion 11...00
				break;
			case 'C':
				hash += 1;
				revcmp_hash += 2LL<<(sizeof(revcmp_hash)*8-2); //10 at the first position 10...00
				break;
			case 'G':
				hash += 2;
				revcmp_hash += 1LL<<(sizeof(revcmp_hash)*8-2); //01 at the first position 01...00
				break;
			case 'A':
				hash += 3;
				break;
			default:
				printf("Invalid character \"%c\" in sequence", *sequence);
				return NULL;
			}
			revcmp_hash = revcmp_hash >> 2;
			hash = hash << 2;
			sequence++;
		}
		//Test the window against stored sRNAs. If a ambiguous character is encountered then rebuild the hash and continue.
		while (i < sequence_length) {
			switch (*sequence) {
			case 'T':
				revcmp_hash += 3LL<<(sizeof(revcmp_hash)*8-2); //01 at the first position 01...00
				hash += 0;
				break;
			case 'C':
				hash += 1;
				revcmp_hash += 2LL<<(sizeof(revcmp_hash)*8-2); //01 at the first position 01...00
				break;
			case 'G':
				hash += 2;
				revcmp_hash += 1LL<<(sizeof(revcmp_hash)*8-2); //01 at the first position 01...00
				break;
			case 'A':
				hash += 3;
				revcmp_hash += 0;
				break;
			default:
				printf("Invalid character \"%c\" in sequence", sequence[i]);
				return NULL;
			}
			//Get counts from the hash
			count = get(seq_read_nums, hash);
			//Put the count in the first position of the window (5' end)
			if (count != -1) {
				fwd_counts[i] += count;
			}
			//Push the reverse complement hash back
			revcmp_hash = revcmp_hash>>num_leading_bits;
			count = get(seq_read_nums, revcmp_hash);
			if (count != -1) {
				//Put the count in the same position
				rev_counts[i] += count;
			}
			//Clear the leading bits from the hash and leave space for the next assignment
			revcmp_hash = revcmp_hash<<(num_leading_bits-2);
			hash = hash << (num_leading_bits+2);
			hash = hash >> num_leading_bits;
			sequence++;
		}
	}





	//Each hash from this point onwards will be a new window of sequence in forward and reverse

	//Assign the output counts
	read_counts->fwd_counts = fwd_counts;
	read_counts->rev_counts = rev_counts;
	return read_counts;
}
;


//The hash is a left-ordered 64 bit integer of the sequence in two-bit format
//Going to be clever with our hash function and deviate from the 2bit format. If we reorder the bases as
//T = 00
//C = 01
//G = 10
//A = 11
//Then the hash of the complementary sequence is just ~hash with the leading ones removed
uint64_t hash_seq(char *sequence, short sequence_length) {
	uint64_t hash;
	int i;
	hash = 0;
	//Don't shift the last element
	for (i = 0; i < sequence_length - 1; ++i) {
		switch (*sequence) {
		case 'T':
			hash += 0;
			break;
		case 'C':
			hash += 1;
			break;
		case 'G':
			hash += 2;
			break;
		case 'A':
			hash += 3;
			break;
		default:
			printf("Invalid character \"%c\"", *sequence);
			return -1;
		}
		++sequence;
		//Edge case of length 1
		if (sequence_length == 1){
			return hash;
		}
		hash = hash << 2;
	}
	switch (*sequence) {
	case 'T':
		break;
	case 'C':
		hash += 1;
		break;
	case 'G':
		hash += 2;
		break;
	case 'A':
		hash += 3;
		break;
	case 'N':
	default:
		printf("Invalid character \"%c\"", *sequence);
		return -1;
	}
	return hash;
}

HashTable *read_sRNAs(char *filename, short required_read_length,
		uint32_t num_bins) {
	FILE *file;
	char *sequence; //Max read length + 1 for null terminator
	uint32_t num_reads_32;
	int num_reads;
	short read_length;
	clock_t begin;
	uint64_t hash, max_bins_required;
	HashTable *seq_read_nums;
	//Assign space for the sequences (with null terminator)
	sequence = malloc(sizeof(char) * (required_read_length + 1));
	//Don't waste space - this is a perfect hash function up to 4^width-1 bins
	//The bit shift is the equivalent to 4 ^ required_read_length
	max_bins_required = 1;
	max_bins_required <<= (required_read_length * 2);
	if (max_bins_required <= num_bins) {
		printf("The perfect hash function only requires %d bins", max_bins_required);
		num_bins = max_bins_required;
	}
	begin = clock();
	seq_read_nums = create_hash_table(num_bins);
	build_hashtable_time += clock()-begin;
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
					num_sRNA_reads, filename);
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
						num_sRNA_reads, filename);
				break;
			}
		}
		num_reads_32 = (uint32_t) num_reads;
		//Add the sequence to the hash table
		hash = hash_seq(sequence, read_length);
		if (hash == -1){
			printf(" at line %d in file %s\n",num_sRNA_reads, filename);
			exit(1);
		}
		put(seq_read_nums, hash, num_reads_32);
	}
	read_sRNAs_time = clock() - begin;
	//Deallocate memory used by sequence
	free(sequence);
	//Free the file
	fclose(file);
	return seq_read_nums;
}

//void read_seqs_and_write_output_whole_file(char *seqs_filename,
//		char *output_filename, HashTable *seq_read_counts) {
//	FILE *seqs_file;
//	FILE *output_file;
//	char *file_contents;
//	long file_size; //The size of the file
//	//Open the output file
//	output_file = fopen(output_filename, "w");
//	if (output_file == NULL ) {
//		//Error reading the file
//		printf("Error reading file %s: %d (%s)\n", output_filename, errno,
//				strerror(errno));
//	}
//	seqs_file = fopen(seqs_filename, "rb");
//	if (seqs_file == NULL ) {
//		//Error reading the file
//		printf("Error reading file %s: %d (%s)\n", seqs_filename, errno,
//				strerror(errno));
//	}
//
////This was for reading the entire file at once - fscanf should be easier
////	//Find the end of the sequence file
////	fseek(seqs_file, 0, SEEK_END);
////	//Get the file size
////	file_size = ftell(seqs_file);
////	//Return to the start of the file
////	rewind(seqs_file);
////	//Assign space for the file contents
////	file_contents = malloc(file_size + 1);
////	if (file_contents == NULL){
////		//Allocation has failed
////		printf("Memory allocation has failed");
////		exit(1);
////	}
////	if (fread(file_contents, file_size, 1, seqs_file) != 1){
////		//File reading has failed
////		printf("Reading contents of file %s : %d (%s)\n",seqs_filename, errno, strerror(errno));
////	}
////	//0 terminate the string
////	file_contents[file_size] = 0;
//	//Iterate through the string to find sequences and headers
//	//Split the string on >
//
//	//Scan the file and store the contents in an array of pointers
//	const int START_SIZE = 50;
//	int seq_ind = 0;
//	char *seqs;
//	char *headers;
//
//	while (!feof(seqs_file)) {
//		if (fscanf(seqs_file, "%[^\n]\n%[^>]", seqs[seq_ind], headers[seq_ind]) != 2) {
//			//Measuring line numbers would be a pain so the sequence number will have to do
//			fprintf("Invalid format in sequence number %d in file %s", (seq_ind+1), seqs_filename);
//		}
//	}
//
//	//Free memory and unbind parameters
////	free(file_contents);
//	fclose(output_file);
//}

void revcmp_seq(char *input, int num_chars, char **output){
	static const char REVERSED_CHAR_HASHES[] = {'A', 'G', 'C', 'T'};
	char *temp;
	int i;
	i = 0;
	temp = *output;
	//Move the input pointer to the end of the char array
	input = input + (num_chars - 1);
	//We've already reached the first position
	while (i < num_chars) {
		//Dereference both pointers and assign output to input
		*temp = REVERSED_CHAR_HASHES[hash_seq(input, 1)];
		input--;
		temp++;
		i++;
	}
	//Return to the original index
	*output = temp - num_chars;
};

void calc_and_write_output(HashTable *sRNA_read_counts, char *input_filename, char* output_filename, short window_size){
	char *seq;
	char *name;
	char *seq_window;
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
	seq_window = malloc(sizeof(char) * (window_size+1));
	//Null terminate the string
	seq_window[window_size] = '\0';
	fasta_file = OpenFASTA(input_filename);
	begin = clock();
	while (ReadFASTA(fasta_file, &seq, &name, &length) != 0){
		read_time += clock()-begin;
		//Sequences less than the window size won't be bound to
		if (length < window_size){
			continue;
		}
		begin = clock();
		read_counts = get_sequence_counts_direct_hash_rev_cmp(sRNA_read_counts, seq, length, window_size);

		scan_time += (clock()-begin);
		if (read_counts == NULL){
			//Invalid char in sequence - complete the error message
			printf(" %s in file %s\n", name, input_filename);
			exit(1);
		}
			begin = clock();
			for (i = 0; i < (length-window_size+1); ++i) {
			//Format: Ref_seq_ID', 'sRNA_seq', 'count', "pos_from_5'_fwd", Fwd_or_Rvs
				count = read_counts->fwd_counts[i];
				if (count != 0){
					//Make a copy of the sequence at the current position
					memcpy(seq_window, &seq[i], window_size);
					if (fprintf(output_file, "%s,%s,%d,%d,%s\n", name, seq_window, count, i,  "f") < 0){
						printf("There was an error writing the output file\n");
						exit(1);
					}
				}
				count = read_counts->rev_counts[i];
				if (count != 0){
					//Flip the sequence to match the output
					revcmp_seq(&seq[i], window_size, &seq_window);
					//TODO - the i+1 is to keep the file the same as the python version
					if (fprintf(output_file, "%s,%s,%d,%d,%s\n", name, seq_window, -count, i+window_size, "r") < 0){
						printf("There was an error writing the output file\n");
						exit(1);
					}
				}
			}
			write_time += clock()-begin;
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
};

//Output format (ascii): pos (32bit int), reverse (bool)
//char *read2bit(char *filename){
//	FILE *file;
//	file = fopen(filename, 'rb');
//
//
//}



//Testing
int main(int argc, char *argv[]) {
//##########################Test reading of sequences########################

	//TODO investigate why bin sizes of ~100,000,000 lead to a dramatic slow down
	//even though there are less hash collisions and its not allocation of space
	//for the hash table
	uint32_t num_bins;
	HashTable *seq_read_counts;
	char *sRNA_filename, *reference_filename, *output_filename;
	short width;
	clock_t begin;
	bool failed;
	printf("SRNA mapper: ");
	if (argc != 6){
		printf("Five arguments required");
		return -1;
	}
	sRNA_filename = argv[1];
	reference_filename = argv[2];
	output_filename = argv[5];
	num_bins = strtol(argv[4], NULL, 10);
	width = (short) strtol(argv[3], NULL, 10);
	seq_read_counts = read_sRNAs(sRNA_filename, width, num_bins);
	printf("Time spent building hashtable: %2.2fs\n",((double) build_hashtable_time)/ CLOCKS_PER_SEC);
	printf("Time taken reading sRNAs: %2.2fs\n", ((double) read_sRNAs_time)/ CLOCKS_PER_SEC);
	begin = clock();
	calc_and_write_output(seq_read_counts,reference_filename, output_filename, width);
	printf("Number of hash collisions: %d from %d sequences - %0.2f%%\n", num_collisions, num_sRNA_reads, ((float) num_collisions) / ((float) num_sRNA_reads) * 100);
	printf("Time taken scanning and matching to sequences: %fs\n", ((double)(clock()-begin))/ CLOCKS_PER_SEC);
	printf("Time spent reading in sequences: %f2.2s\n",((double) read_time)/ CLOCKS_PER_SEC);
	printf("Time spent scanning sequences for matches: %f2.2s\n",((double) scan_time)/ CLOCKS_PER_SEC);
	printf("Time spent writing output file: %2.2fs\n",((double) write_time)/ CLOCKS_PER_SEC);
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

