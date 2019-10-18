#设置cmake的最低版本
cmake_minimum_required(VERSION 3.7)
#设置项目名称
project(GreedyCut)
#识别使用的编译器
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "编译器是gcc")
    add_compile_options(-std=c++11)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(STATUS "编译器是msvc")
endif()

SET(CMAKE_CXX_FLAGS "-O2")
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("发现OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#包含mkl的头文件目录和库
include_directories("D:/LIBRARY/mkl/compilers_and_libraries_2019.4.245/windows/mkl/include")
link_directories("D:/LIBRARY/mkl/compilers_and_libraries_2019.4.245/windows/mkl/lib/intel64")

#包含Eigen的头文件和库
include_directories("D:/LIBRARY")
#包含OpenMesh的头文件目录和库
include_directories("D:/LIBRARY/OpenMesh/src")

#使用gcc编译的OpenMesh库
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    link_directories("D:/LIBRARY/OpenMesh/build_cmake/Build/lib")
#使用MSVC编译的OpenMesh库
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    link_directories("D:/LIBRARY/OpenMesh/build_vs/Build/lib")
endif()
link_libraries(OpenMeshTools)
link_libraries(OpenMeshCore)
file (GLOB SRCFILES
    src/main.cpp
    src/MeshDefinition/MeshDefinition.cpp
    src/MeshDefinition/MeshCache.cpp 
    src/Auxiliary/OptionReader.cpp
    src/Auxiliary/Algorithm.cpp
    src/PointSampling/PointSampling.cpp
    src/MeshCut/MeshCut.cpp
    src/KPNewton/KPNewton.cpp
    src/PardisoSolver/PardisoSolver.cpp
)

add_executable(GreedyCut ${SRCFILES})
target_link_libraries(GreedyCut OpenMeshTools)
target_link_libraries(GreedyCut OpenMeshCore)
target_link_libraries(GreedyCut mkl_rt)