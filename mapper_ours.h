/*
 * util.h

 *
 *  Created on: 23/10/2014
 *      Author: Jason
 */
//2Bit hash table
extern HashTable *create_hash_table_2bit(uint32_t);
extern void put_2bit(HashTable *hash_table, uint64_t, uint32_t);
extern uint32_t get_2bit(HashTable*, uint64_t);
extern void iterate_and_add_2bit(LinkedList*, uint64_t,
		uint32_t);
extern uint64_t hash_seq(char*, short);

