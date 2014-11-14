/*
 * util.h
 *
 *  Created on: 23/10/2014
 *      Author: Jason
 */
#include <inttypes.h>



extern void reverse_char_array(char*, char**, int);

typedef struct _ReadCounts_ {
	uint32_t *fwd_counts;
	uint32_t *rev_counts;
} ReadCounts;
//Define boolean operators
typedef enum {
	false, true
} bool;
