//
//  UtilitiesOMP.h
//  Gauss
//
//  Created by David Levin on 6/7/17.
//
//

#ifndef UtilitiesOMP_h
#define UtilitiesOMP_h

inline int omp_thread_count() {
    int n = 0;
#pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

#endif /* UtilitiesOMP_h */
