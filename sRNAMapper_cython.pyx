#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "util.h"
#include "fasta.h" //Header file from fasta library
#include "mapper_ours.h"

#Timing Variables - These could be removed with the pre-processor somehow
cdef clock_t read_time = 0;
cdef clock_t write_time = 0;
cdef clock_t scan_time = 0;
cdef clock_t build_hashtable_time = 0;
cdef clock_t read_sRNAs_time = 0;
#More information about hash collisions and the total number of reads
cdef int num_collisions = 0;
cdef long num_sRNA_reads = 0;

