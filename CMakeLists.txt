cmake_minimum_required(VERSION 3.7)
project(GreedyCut)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "the compiler is gcc")
    SET(CMAKE_CXX_FLAGS "-O2")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(STATUS "the compiler is MSVC")
    add_definitions(/D_USE_MATH_DEFINES)
endif()

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message(STATUS "Find OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option(USE_PARDISO "Use Pardiso" ON)
option(USE_MKL_PARDISO "Use MKL Pardiso" OFF)
option(USE_EIGEN "Use Eigen" OFF)

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "the OS is windows")
    include_directories("D:/LIBRARY/Eigen")
    include_directories("D:/LIBRARY/OpenMesh/src")

    if (USE_MKL_PARDISO)
        message(STATUS "Using MKL Pardiso")
        add_definitions(-DUSE_MKL_PARDISO)
        include_directories("D:/LIBRARY/mkl/compilers_and_libraries_2019.4.245/windows/mkl/include")
        link_directories("D:/LIBRARY/mkl/compilers_and_libraries_2019.4.245/windows/mkl/lib/intel64")
        SET(SOLVER_LIBRARY mkl_rt)
        SET(SOLVER MKLPardisoSolver.cpp)
    elseif(USE_PARDISO)
        message(STATUS "Using Pardiso")
        add_definitions(-DUSE_PARDISO)
        link_directories("D:/LIBRARY/Pardiso")
        SET(SOLVER_LIBRARY libpardiso600-WIN-X86-64.lib)
        SET(SOLVER PardisoSolver.cpp)
    elseif(USE_EIGEN)
        message(STATUS "Using Eigen")
        add_definitions(-DUSE_EIGEN)
        SET(SOLVER EigenLinSolver.cpp)
    endif()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        link_directories("D:/LIBRARY/OpenMesh/build_cmake/Build/lib")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        link_directories("D:/LIBRARY/OpenMesh/build_vs/Build/lib")
    endif()

elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
    message(STATUS "the OS is Linux")
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    message(STATUS "the OS is Darwin")
    
endif()

file (GLOB SRCFILES
    src/main.cpp
    src/Auxiliary.cpp
    src/PointSampling.cpp
    src/MeshCache.cpp
    src/MeshCut.cpp
    src/KPNewton.cpp
    src/Solver/${SOLVER}
    src/PointFinding.cpp
    src/GAP.cpp
    src/AddAuxiliaryPoint.cpp
)

add_executable(${PROJECT_NAME} ${SRCFILES})
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_link_libraries(${PROJECT_NAME} debug OpenMeshToolsd debug OpenMeshCored)
    target_link_libraries(${PROJECT_NAME} optimized OpenMeshTools optimized OpenMeshCore)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_link_libraries(${PROJECT_NAME}  OpenMeshTools OpenMeshCore)
endif()

target_link_libraries(${PROJECT_NAME} ${SOLVER_LIBRARY})