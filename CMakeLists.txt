cmake_minimum_required(VERSION 2.8)
SET(CMAKE_BUILD_TYPE "Debug")  
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")  
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")  
SET(EPICS_BASE "/home/ray/EPICS/epics-base-R3.14.12.7")
project( ReconstructBeamspot )
find_package( OpenCV REQUIRED )
include_directories(${EPICS_BASE}/include ${EPICS_BASE}/include/os/Linux)
link_directories(${EPICS_BASE}/lib/linux-x86_64)
set(SRC_LIST reconstructBeamspot.cpp precoimgPixel.cpp)
add_executable( ReconstructBeamspot reconstructBeamspot.cpp precoimgPixel.cpp fft.cpp)
target_link_libraries( ReconstructBeamspot ${OpenCV_LIBS} -lpthread -lm libca.so libCom.so)
