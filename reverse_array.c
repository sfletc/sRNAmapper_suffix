/*
 * reverse_array.c
 *
 *  Created on: 23/10/2014
 *      Author: Jason
 */
#include "util.h"
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
