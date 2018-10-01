cmake_minimum_required(VERSION 2.6)
project(genomestore)

#set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_GENERATOR Visual Studio 15 2017 Win64)
#set(CMAKE_CXX_FLAGS "-O -std=c++1y")
# set(CMAKE_CXX_COMPILER "clang++")
# set(CMAKE_CXX_FLAGS "-O3 -std=c++1y")

set(SOURCE_DIR "./src")
set(LMDB_DIR "./lmdb/libraries/liblmdb/")
set(JSON_DIR "./json/include")
set(PYTHON_ROOT "D:/Python37")
if (WIN32)
    set(PYTHON_LIBRARY "${PYTHON_ROOT}/libs/python37.lib")
endif (WIN32)
set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include")
set(PYBIND_INCLUDE_DIR "D:/code/pybind11/include")
set(PYBIND_DIR "D:/code/pybind11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib")

#configure_file("${SOURCE_DIR}/base.hpp")

INCLUDE_DIRECTORIES(
    ${SOURCE_DIR}
    ${LMDB_DIR}
    ${JSON_DIR}
)

LINK_DIRECTORIES(
    )



message("Configuring core C++ binaries")
# file (GLOB MOODS_SRC 
#     "${SOURCE_DIR}/MOODS/*.h"
#     "${SOURCE_DIR}/MOODS/*.cpp"
#     )
# add_library(moods SHARED ${MOODS_SRC})
# INCLUDE_DIRECTORIES ("${SOURCE_DIR}/MOODS")
FIND_PACKAGE(PythonLibs 3.7)
add_library(llmdb ${LMDB_DIR}/midl.c ${LMDB_DIR}/lmdb.h ${LMDB_DIR}/mdb.c)
add_library(lmdb SHARED ${LMDB_DIR}/midl.c ${LMDB_DIR}/lmdb.h ${LMDB_DIR}/mdb.c)
add_library(lindexer ${SOURCE_DIR}/indexer.cpp ${SOURCE_DIR}/base.cpp)
add_library(indexer SHARED ${SOURCE_DIR}/indexer.cpp ${SOURCE_DIR}/base.cpp)
target_link_libraries(lindexer llmdb)
target_link_libraries(indexer llmdb)
# add_executable(mafslice.exe ${SOURCE_DIR}/main.cpp)
# target_link_libraries(mafslice.exe indexer boost_program_options)
IF(PYTHONLIBS_FOUND)
    message("PYTHON 3.7 FOUND, configuring python interface")
    # INCLUDE_DIRECTORIES ("${Boost_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES ( ${PYTHON_INCLUDE_DIR} )
    INCLUDE_DIRECTORIES ( ${PYBIND_INCLUDE_DIR} )
    add_subdirectory(${PYBIND_DIR} build)
    # SET(Boost_USE_STATIC_LIBS OFF)
    # SET(Boost_USE_MULTITHREADED ON)
    # SET(Boost_USE_STATIC_RUNTIME OFF)
    # FIND_PACKAGE(Boost 1.59.0 COMPONENTS python)
    # add_library(motifmapdb SHARED ${SOURCE_DIR}/motifmapdb.cpp ${SOURCE_DIR}/motifmapdb.hpp)
    # target_link_libraries(motifmapdb indexer moods ${Boost_LIBRARIES})
    # add_library(motifmapcore SHARED ${SOURCE_DIR}/motifmapcore.cpp)
    # target_link_libraries(motifmapcore motifmapdb ${Boost_LIBRARIES})
    # set_target_properties(motifmapcore PROPERTIES PREFIX "" )
    pybind11_add_module(genomestore ${SOURCE_DIR}/genomestore.cpp)
    target_link_libraries(genomestore PRIVATE llmdb lindexer)
ELSE()
    message("Error finding boost...")
ENDIF()
