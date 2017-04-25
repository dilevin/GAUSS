#include "CoreOSSpecific.h"

#ifdef WIN32
//WINDOWS INCLUDES
#else
//*NIX includes
#include <dirent.h>
#endif

#ifdef WIN32

int Core::getDirectoryListing(std::vector<std::string> &list, const std::string &directoryName)
{
    std::cout<<"Windows getDirectoryListing not implemented \n";
    return 0;
}

#else

int Core::getDirectoryListing(std::vector<std::string> &list, const std::string &directoryName)
{
    struct dirent *dp;

    DIR *dfd = opendir(directoryName.c_str());

    list.clear();

    if(dfd != NULL) {
        while((dp = readdir(dfd)) != NULL)
        {
            list.push_back(dp->d_name);
        }

        closedir(dfd);
    } else {
        return 0;
    }

    return 1;
}

#endif

std::string Core::getDirectory(const std::string &pathName)
{
    //find first path seperator
    return pathName.substr(0,pathName.find_last_of(PATH_SEPARATOR));
    
}

std::string Core::getFilename(const std::string &pathName)
{
    return pathName.substr(pathName.find_last_of(PATH_SEPARATOR)+1);
}

