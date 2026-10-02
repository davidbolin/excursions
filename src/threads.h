#ifndef EXCURSIONS_THREADS_H
#define EXCURSIONS_THREADS_H

#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

// Number of threads to request for the parallel regions. With n_threads = 0
// this is the OpenMP default, which respects OMP_NUM_THREADS, and otherwise
// n_threads, at most the number of processors. In both cases it is at most
// OMP_THREAD_LIMIT. The runtime may still give fewer threads, so the parallel
// regions must use the size of the team they get. The number of threads is
// requested with num_threads() rather than omp_set_num_threads(), which would
// change the default for later calls and for other packages.
static inline int excursions_threads(int n_threads) {
#ifdef _OPENMP
  int nP = (n_threads <= 0) ? omp_get_max_threads()
                            : std::min(n_threads, omp_get_num_procs());
  nP = std::min(nP, omp_get_thread_limit());
  return std::max(nP, 1);
#else
  (void) n_threads;
  return 1;
#endif
}

#endif
