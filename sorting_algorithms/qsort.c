#include <stdio.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>

#include <x86intrin.h>

#include "sorting.h"


static int compare (const void *x, const void *y)
{
    /* TODO: comparison function to be used by qsort()*/

    /* cast x and y to uint64_t* before comparing */

  return (*(uint64_t*)x - *(uint64_t*)y);
    
}

void sequential_qsort_sort (uint64_t *T, const int size)
{

    /* TODO: sequential sorting based on libc qsort() function */

    qsort(T, size, sizeof(uint64_t), compare);

    return ;
}

/* 
   Merge two sorted chunks of array T!
   The two chunks are of size size
   First chunck starts at T[0], second chunck starts at T[size]
*/
void merge (uint64_t *T, const uint64_t size)
{
  uint64_t *X = (uint64_t *) malloc (2 * size * sizeof(uint64_t)) ;
  
  uint64_t i = 0 ;
  uint64_t j = size ;
  uint64_t k = 0 ;
  
  while ((i < size) && (j < 2*size))
    {
      if (T[i] < T [j])
	{
	  X [k] = T [i] ;
	  i = i + 1 ;
	}
      else
	{
	  X [k] = T [j] ;
	  j = j + 1 ;
	}
      k = k + 1 ;
    }

  if (i < size)
    {
      for (; i < size; i++, k++)
	{
	  X [k] = T [i] ;
	}
    }
  else
    {
      for (; j < 2*size; j++, k++)
	{
	  X [k] = T [j] ;
	}
    }
  
  memcpy (T, X, 2*size*sizeof(uint64_t)) ;
  free (X) ;
  
  return ;
}

void parallel_qsort_sort (uint64_t *T, const uint64_t size)
{

    /* TODO: parallel sorting based on libc qsort() function +
     * sequential merging */

    uint64_t size_of_chunks;
    int numt;

    #pragma omp parallel default(none) shared(T, size, size_of_chunks, numt)
    {
      int tid = omp_get_thread_num();
      numt = omp_get_num_threads();
      char flag = 0;

      // If number of threads are not propotional to the size of the array
      if(size%numt != 0)
      {
        size_of_chunks = size/numt + 1;
        flag = 1;
      }
      else
      {
        size_of_chunks = size/numt;
      }

      if(tid == (numt - 1) && flag)
      {
        qsort(T+(size_of_chunks*tid), size_of_chunks, sizeof(uint64_t), compare);
      }
      else
      {
        qsort(T+(size_of_chunks*tid), size_of_chunks, sizeof(uint64_t), compare);
      }        
    }

    // Sequential merge
    for(uint64_t i = 0; i <= size/2; i += 2*size_of_chunks)
    {
      merge(T+i, size_of_chunks);
    }
    merge(T, size/2);

    return;

}


void parallel_qsort_sort1 (uint64_t *T, const uint64_t size)
{

    /* TODO: parallel sorting based on libc qsort() function +
     * PARALLEL merging */

    uint64_t size_of_chunks;
    int numt;

    #pragma omp parallel default(none) shared(T, size, size_of_chunks, numt)
    {
      int tid = omp_get_thread_num();
      numt = omp_get_num_threads();
      char flag = 0;

      // If no of threads are not propotional to the size of the array
      if(size%numt != 0)
      {
        size_of_chunks = size/numt + 1;
        flag = 1;
      }
      else
      {
        size_of_chunks = size/numt;
      }

      if(tid == (numt - 1) && flag)
      {
        qsort(T+(size_of_chunks*tid), size_of_chunks, sizeof(uint64_t), compare);
      }
      else
      {
        qsort(T+(size_of_chunks*tid), size_of_chunks, sizeof(uint64_t), compare);
      }        
    } 

    // Parallel Merge
    #pragma omp parallel for default(none), shared(T, size, size_of_chunks)
    
    for(uint64_t i = 0; i <= size/2; i += 2*size_of_chunks)
    {
      merge(T+i, size_of_chunks);
    }
    merge(T, size/2);

    return;

}

int main (int argc, char **argv)
{
    uint64_t start, end;
    uint64_t av ;
    unsigned int exp ;

    /* the program takes one parameter N which is the size of the array to
       be sorted. The array will have size 2^N */
    if (argc != 2)
    {
        fprintf (stderr, "qsort.run N \n") ;
        exit (-1) ;
    }

    uint64_t N = 1 << (atoi(argv[1])) ;
    /* the array to be sorted */
    uint64_t *X = (uint64_t *) malloc (N * sizeof(uint64_t)) ;

    printf("--> Sorting an array of size %u\n",N);
#ifdef RINIT
    printf("--> The array is initialized randomly\n");
#endif
    

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++){
#ifdef RINIT
        init_array_random (X, N);
#else
        init_array_sequence (X, N);
#endif
        
      
        start = _rdtsc () ;
        
        sequential_qsort_sort (X, N) ;
     
        end = _rdtsc () ;
        experiments [exp] = end - start ;

        /* verifying that X is properly sorted */
#ifdef RINIT
        if (! is_sorted (X, N))
        {
            fprintf(stderr, "ERROR: the sequential sorting of the array failed\n") ;
            print_array (X, N) ;
            exit (-1) ;
	}
#else
        if (! is_sorted_sequence (X, N))
        {
            fprintf(stderr, "ERROR: the sequential sorting of the array failed\n") ;
            print_array (X, N) ;
            exit (-1) ;
	}
#endif
    }

    av = average_time() ;  

    double serial = (double)av/1000000;

    printf ("\n qsort serial \t\t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

  
    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
#ifdef RINIT
        init_array_random (X, N);
#else
        init_array_sequence (X, N);
#endif
        
        start = _rdtsc () ;
        
        parallel_qsort_sort (X, N) ;

     
        end = _rdtsc () ;
        experiments [exp] = end - start ;

        /* verifying that X is properly sorted */
#ifdef RINIT
        if (! is_sorted (X, N))
        {
            fprintf(stderr, "ERROR: the parallel sorting of the array failed\n") ;
            exit (-1) ;
	}
#else
        if (! is_sorted_sequence (X, N))
        {
            fprintf(stderr, "ERROR: the parallel sorting of the array failed\n") ;
            exit (-1) ;
	}
#endif
                
        
    }
    
    av = average_time() ;  
    printf ("\n qsort parallel (merge seq) \t\t %.2lf Mcycles\n\n", (double)av/1000000) ;
    
    


    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
#ifdef RINIT
        init_array_random (X, N);
#else
        init_array_sequence (X, N);
#endif
        
        start = _rdtsc () ;

        parallel_qsort_sort1 (X, N) ;
     
        end = _rdtsc () ;
        experiments [exp] = end - start ;

        /* verifying that X is properly sorted */
#ifdef RINIT
        if (! is_sorted (X, N))
        {
            fprintf(stderr, "ERROR: the parallel sorting of the array failed\n") ;
            exit (-1) ;
	}
#else
        if (! is_sorted_sequence (X, N))
        {
            fprintf(stderr, "ERROR: the parallel sorting of the array failed\n") ;
            exit (-1) ;
	}
#endif
                
        
    }
    
    av = average_time() ;  
    printf("\n qsort parallel \t\t %.2lf Mcycles\n\n", (double)av/1000000) ;
    
    

    /* before terminating, we run one extra test of the algorithm */
    uint64_t *Y = (uint64_t *) malloc (N * sizeof(uint64_t)) ;
    uint64_t *Z = (uint64_t *) malloc (N * sizeof(uint64_t)) ;

#ifdef RINIT
    init_array_random (Y, N);
#else
    init_array_sequence (Y, N);
#endif

    memcpy(Z, Y, N * sizeof(uint64_t));

    sequential_qsort_sort (Y, N) ;
    parallel_qsort_sort1 (Z, N) ;

    if (! are_vector_equals (Y, Z, N)) {
        fprintf(stderr, "ERROR: sorting with the sequential and the parallel algorithm does not give the same result\n") ;
        exit (-1) ;
    }


    free(X);
    free(Y);
    free(Z);
    
}
