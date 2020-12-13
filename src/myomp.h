// Header file to determine when OpenMP should be used
// Copied from the data.table package under the Mozilla public license 2.0
// https://github.com/Rdatatable/data.table/blob/master/configure
#ifdef _OPENMP
  // for machines with compiles with OpenMP support
  #include <omp.h>
  #if _OPENMP >= 201511
    #define monotonic_dynamic monotonic:dynamic // #4786
  #else
    #define monotonic_dynamic dynamic
  #endif
  #define MY_OPENMP              _OPENMP
#else
  // for machines with compilers void of OpenMP support
  #define omp_get_num_threads()  1
  #define omp_get_thread_num()   0
  #define omp_get_max_threads()  1
  #define omp_get_thread_limit() 1
  #define omp_get_num_procs()    1
  #define omp_set_nested(a)   // empty statement to remove the call
  #define omp_get_wtime()        0
  #define MY_OPENMP              0
#endif
