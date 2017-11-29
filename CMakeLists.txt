cmake_minimum_required(VERSION 2.8)
project(TCadToGeo)

set(CMAKE_CXX_STANDARD 11)

# Extend cmake module search path to provide FindTetgen.cmake.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# Find tetgen library, e.g. via LD_LIBRARY_PATH
find_package(Tetgen REQUIRED)
include_directories(${TETGEN_INCLUDE_DIRS})

# Find ROOT, Minimum Version 6
find_package(ROOT 6 REQUIRED)
message(STATUS "ROOT Version ${ROOT_VERSION} found in ${ROOT_root_CMD}")
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

# Build
root_generate_dictionary(G__TCadToGeo src/TCadToGeo.h src/TGeoTet.h MODULE TCadToGeo LINKDEF src/TCadToGeoLinkDef.h)
set(SOURCE_FILES src/TCadToGeo.cxx src/TGeoTet.cxx G__TCadToGeo.cxx)
add_library(TCadToGeo SHARED ${SOURCE_FILES})
## Geom is not included in ROOT_LIBRARIES ?!
target_link_libraries(TCadToGeo ${TETGEN_LIBRARIES} ${ROOT_LIBRARIES} Geom)

# Install
install(TARGETS TCadToGeo DESTINATION lib)
install(FILES ${CMAKE_BINARY_DIR}/libTCadToGeo_rdict.pcm ${CMAKE_BINARY_DIR}/libTCadToGeo.rootmap DESTINATION lib)
install(FILES src/TCadToGeo.h src/TGeoTet.h DESTINATION include)
