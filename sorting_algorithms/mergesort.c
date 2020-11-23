#include <stdio.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>

#include <x86intrin.h>

#include "sorting.h"

/* 
   merge sort -- sequential, parallel -- 
*/

/* Function to merge the two haves arr[l..m] and arr[m+1..r] of array arr[] */
void merges(uint64_t arr[], uint64_t l, uint64_t m, uint64_t r)
{
    uint64_t i, j, k;
    uint64_t n1 = m - l + 1;
    uint64_t n2 =  r - m;
 
    /* create temp arrays */
    uint64_t L[n1], R[n2];
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void sequential_merge_sort (uint64_t *T, const uint64_t size)
{
    /* sequential implementation of merge sort */

    if(size == 2)
    {
        if(T[0] > T[1])
        {
            uint64_t tmp = T[0];
            T[0] = T[1];
            T[1] = tmp;
            return;
        }
        else
        {
            return;
        }        
    }

    for(int i = 0; i < 2; i++)
    {
        sequential_merge_sort(T + (i*size/2), size/2);
    }

    merges(T, 0, size/2-1, size-1);

    return;    
}

void parallel_merge_sort (uint64_t *T, const uint64_t size)
{
    /* parallel implementation of merge sort */

    if(size == 2)
    {
        if(T[0] > T[1])
        {
            uint64_t tmp = T[0];
            T[0] = T[1];
            T[1] = tmp;
            return;
        }
        else
        {
            return;
        }        
    }

    #pragma omp task shared(T, size)
    parallel_merge_sort(T, size/2);

    #pragma omp task shared(T, size)
    parallel_merge_sort(T + size/2, size/2);

    #pragma omp taskwait        

    merges(T, 0, size/2-1, size-1);

    return;
}

void parallel_merge_sort1 (uint64_t *T, const uint64_t size)
{
    /* optimized parallel implementation of merge sort */

    int curr_size;

   int left_start;
 
   for (curr_size=1; curr_size<size; curr_size = 2*curr_size)
   {       
       #pragma omp parallel for default(none), shared(T, size, curr_size, left_start)
       for (left_start=0; left_start<size-1; left_start += 2*curr_size)
       {
           int mid = left_start + curr_size - 1;
 
           int right_end = left_start + 2*curr_size - 1 < size-1? left_start + 2*curr_size - 1 : size-1;
 
           merges(T, left_start, mid, right_end);
       }
   }

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
        fprintf (stderr, "merge.run N \n") ;
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
        
        sequential_merge_sort (X, N) ;
     
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

    printf ("\n mergesort serial \t\t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

    double serial = (double)av/1000000;
  
    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
#ifdef RINIT
        init_array_random (X, N);
#else
        init_array_sequence (X, N);
#endif
        
        start = _rdtsc () ;

        parallel_merge_sort (X, N) ;
     
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
    printf ("\n mergesort parallel \t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
#ifdef RINIT
        init_array_random (X, N);
#else
        init_array_sequence (X, N);
#endif
        
        start = _rdtsc () ;

        parallel_merge_sort1 (X, N) ;
     
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
            fprintf(stderr, "ERROR: the parallel sorting(optimized) of the array failed\n") ;
            exit (-1) ;
	}
#endif
                
        
    }

     av = average_time() ;  
    printf ("\n mergesort parallel(optimized) \t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

    /* print_array (X, N) ; */

    /* before terminating, we run one extra test of the algorithm */
    uint64_t *Y = (uint64_t *) malloc (N * sizeof(uint64_t)) ;
    uint64_t *Z = (uint64_t *) malloc (N * sizeof(uint64_t)) ;

#ifdef RINIT
    init_array_random (Y, N);
#else
    init_array_sequence (Y, N);
#endif

    memcpy(Z, Y, N * sizeof(uint64_t));

    sequential_merge_sort (Y, N) ;
    parallel_merge_sort (Z, N) ;

    if (! are_vector_equals (Y, Z, N)) {
        fprintf(stderr, "ERROR: sorting with the sequential and the parallel algorithm does not give the same result\n") ;
        exit (-1) ;
    }


    free(X);
    free(Y);
    free(Z);
    
}
