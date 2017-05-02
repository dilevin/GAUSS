//
//  Utilities.cpp
//  Gauss
//
//  Created by David Levin on 5/2/17.
//
//

#include <Utilities.h>

#include <chrono>
#include <ctime>

std::string Gauss::dataDir() {
    return std::string(DataDir(GAUSS_DATA_DIR));
}

std::string Gauss::timeStampString(std::string toStamp) {
    
    std::stringstream stampedString;
    
    //append date and time to this string
    std::chrono::time_point<std::chrono::high_resolution_clock> currentTime = std::chrono::high_resolution_clock::now();
    
    stampedString <<toStamp<<"_"<<currentTime.time_since_epoch().count();
    
    return stampedString.str();
    
}
