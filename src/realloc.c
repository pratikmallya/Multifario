#if HAVE_CONFIG_H
#include <multifarioConfig.h>
#endif

#undef realloc

#include <sys/types.h>

void *realloc (void*,size_t);

/* Allocate an N-byte block of memory from the heap. If N is zero, allocate a 1-byte block. */

void * rpl_realloc (void *ptr, size_t n) { if (n == 0) n = 1; return realloc (ptr,n); } 
