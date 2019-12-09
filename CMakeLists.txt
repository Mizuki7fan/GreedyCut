#设置cmake的最低版本
cmake_minimum_required(VERSION 3.7)
#设置项目名称
project(GreedyCut)
#识别使用的编译器
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "编译器是gcc")
    SET(CMAKE_CXX_FLAGS "-O3")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(STATUS "编译器是msvc")
    #添加预处理器定义，/Db表示是MSVC的编译器
    add_definitions(/D_USE_MATH_DEFINES)
endif()

#需要并行加速，调用OpenMP
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message(STATUS "发现OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#指定线性方程组求解器，通过ON和OFF来控制
option(USE_PARDISO "Use Pardiso" ON)
option(USE_MKL_PARDISO "Use MKL Pardiso" OFF)

#检测操作系统并且进行配置
if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "操作系统是Windows")
    #指定Eigen的头文件目录
    include_directories("D:/LIBRARY/Eigen")
    #包含OpenMesh的头文件目录和库
    include_directories("D:/LIBRARY/OpenMesh/src")

    if (USE_MKL_PARDISO)
        message(STATUS "使用MKL Pardiso")
        add_definitions(-DUSE_MKL_PARDISO)
        include_directories("D:/LIBRARY/mkl/compilers_and_libraries_2019.4.245/windows/mkl/include")
        link_directories("D:/LIBRARY/mkl/compilers_and_libraries_2019.4.245/windows/mkl/lib/intel64")
        SET(PARDISO_LIBRARY mkl_rt)
        SET(PARDISO_SOLVER MKLPardisoSolver.cpp)
    elseif(USE_PARDISO)
        message(STATUS "使用Pardiso")
        add_definitions(-DUSE_PARDISO)
        link_directories("D:/LIBRARY/Pardiso/")
        SET(PARDISO_LIBRARY libpardiso600-WIN-X86-64.lib)
        SET(PARDISO_SOLVER PardisoSolver.cpp)
    endif()

    #使用gcc编译的OpenMesh库
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        link_directories("D:/LIBRARY/OpenMesh/build_cmake/Build/lib")
        #使用MSVC编译的OpenMesh库
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        link_directories("D:/LIBRARY/OpenMesh/build_vs/Build/lib")
    endif()

elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
    message(STATUS "操作系统是Linux")
endif()

file (GLOB SRCFILES
    src/main.cpp
    src/MeshDefinition/MeshDefinition.cpp
    src/MeshDefinition/MeshCache.cpp 
    src/Auxiliary/OptionReader.cpp
    src/Auxiliary/Algorithm.cpp
    src/PointSampling/PointSampling.cpp
    src/MeshCut/MeshCut.cpp
    src/KPNewton/KPNewton.cpp
    src/Solver/${PARDISO_SOLVER}
    src/PointFinding/PointFinding.cpp
    src/AddAuxiliaryPoint/AddAuxiliaryPoint.cpp
    src/GAP/GAP.cpp
)

add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} debug OpenMeshToolsd debug OpenMeshCored)
target_link_libraries(${PROJECT_NAME} optimized OpenMeshTools optimized OpenMeshCore)
target_link_libraries(${PROJECT_NAME} ${PARDISO_LIBRARY})