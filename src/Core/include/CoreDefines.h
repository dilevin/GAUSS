#ifndef _COREDEFINES
#define _COREDEFINES

/***********************************************************************
* Some useful preprocessor definitions
* Please add the approriate ifdef/endif blocks for various operating systems if you find any 
* incompatibilities
*/

#include <cstring>
#include <assert.h> 
#include <iostream>

#ifdef WIN32
    #define MEMCPY(dst, size, src, num) (memcpy_s((dst), (size), (src), (num)))
#else
    #define MEMCPY(dst, size, src, num) (memcpy((dst), (src), (num)))
#endif

#ifdef WIN32
    #define MEMSET
#else
    #define MEMSET(dst, val, num) (memset((dst), (val), (num)))
#endif

#ifdef DEBUG
    #define DEBUGPRINT(a) (std::cout<<a<<"\n")
#else
    #define DEBUGPRINT(a)
#endif

//A few useful, small methods
namespace Core
{
    
}

#endif