cmake_minimum_required(VERSION 2.6)

project(genomestore)

set(SOURCE_DIR "./src")
set(LMDB_DIR "./lmdb/libraries/liblmdb/")
set(JSON_DIR "./json/include")
set(PYBIND_DIR "./pybind11")
set(PYBIND_INCLUDE_DIR "./pybind11/include")
IF(WIN32)
    set(CMAKE_GENERATOR "Visual Studio 15 2017 Win64")
    set(PYTHON_ROOT "D:/Python37")
    set(PYTHON_LIBRARY "${PYTHON_ROOT}/libs/python37.lib")
    set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include")
ELSEIF(APPLE)
    set(CMAKE_CXX_FLAGS "-O3 -std=c++17 -stdlib=libc++")
ELSEIF(UNIX AND NOT APPLE)
    set(PYTHON_LIBRARY   "/usr/lib/python3.6")
    set(PYTHON_INCLUDE_DIR   "/usr/include/python3.6m")
    set(CMAKE_CXX_FLAGS "-O3 -std=c++17 -fPIC -lm")
ENDIF()

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
IF(WIN32)
    target_link_libraries(lindexer llmdb)
    target_link_libraries(indexer llmdb)
ELSE()
    target_link_libraries(lindexer lmdb)
    target_link_libraries(indexer lmdb)
ENDIF()
# add_executable(mafslice.exe ${SOURCE_DIR}/main.cpp)
# target_link_libraries(mafslice.exe indexer boost_program_options)
IF(PYTHONLIBS_FOUND)
    message("PYTHON 3.7 FOUND, configuring python interface")
    INCLUDE_DIRECTORIES ( ${PYTHON_INCLUDE_DIR} )
    INCLUDE_DIRECTORIES ( ${PYBIND_INCLUDE_DIR} )
    add_subdirectory(${PYBIND_DIR} build)
    pybind11_add_module(genomestore ${SOURCE_DIR}/genomestore.cpp)
    IF(WIN32) # for windows only, compile against static
        target_link_libraries(genomestore PRIVATE llmdb lindexer)
    ELSE()
        target_link_libraries(genomestore PRIVATE lmdb indexer)
    ENDIF()
ELSE()
    message("Error finding pythonlib...")
ENDIF()
