#ifndef _COREOSSPECIFIC_H
#define _COREOSSPECIFIC_H

/**
 * OS specific functions
 */
#include <iostream>
#include <string>
#include <vector>

#if defined(WIN32) || defined(_WIN32) || defined(_WIN64)
#define PATH_SEPARATOR '\\'
#else
#define PATH_SEPARATOR '/'
#endif

namespace Core
{
    std::string getDirectory(const std::string &pathName);
 
    std::string getFilename(const std::string &pathName);
    
    int getDirectoryListing(std::vector<std::string> &list, const std::string &directoryName);
}

#endif // COREOSSPECIFIC_H
