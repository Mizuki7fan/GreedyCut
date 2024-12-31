cmake_minimum_required(VERSION 3.7)
project(GreedyCut)

#设置不生成MinSizeRel;RelWithDebInfo
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

add_definitions(/D_USE_MATH_DEFINES)
SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置生成exe的存放位置
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_SOURCE_DIR}/bin/Debug")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/bin/Release")

# 读取设置系统环境变量
set (MyLib $ENV{MyLibPath})
if(NOT MyLib)
    message(FATAL_ERROR "没有在环境变量中设置本地库目录MyLibPath的值.")
else()
    message("本地库目录为: ${MyLib}")
endif()

#OpenMP设置
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    find_package(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message(STATUS "Find OpenMP")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OPENMP_PATH "C:\\Program Files\\LLVM\\lib\\clang\\19\\include")
    set(OPENMP_LIB_PATH "C:\\Program Files\\LLVM\\lib\\libomp.lib")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    include_directories(SYSTEM ${OPENMP_PATH})
endif()

SET(EIGEN_PATH "${MyLib}\\Eigen")
SET(OPENMESH_PATH "${MyLib}\\OpenMesh\\src")
SET(OPENMESH_LIB_PATH "${MyLib}\\OpenMesh\\build\\Build\\lib")
SET(MKL_PATH "C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2025.0\\include")
SET(MKL_LIB_PATH "C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2025.0\\lib")

option(USE_MKL "Use MKL Pardiso" OFF)
option(USE_EIGEN "Use Eigen" ON)

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "the OS is windows")
    include_directories(SYSTEM ${EIGEN_PATH})
    include_directories(SYSTEM ${OPENMESH_PATH})
    link_directories(${OPENMESH_LIB_PATH})

    if (USE_MKL)
        message(STATUS "Using MKL Pardiso")
        add_definitions(-DUSE_MKL_PARDISO)
        include_directories(SYSTEM ${MKL_PATH})
        link_directories(${MKL_LIB_PATH})
        SET(SOLVER_LIBRARY mkl_rt)
        SET(SOLVER_ MKLPardisoSolver)
    elseif(USE_EIGEN)
        message(STATUS "Using Eigen")
        add_definitions(-DUSE_EIGEN)
        SET(SOLVER_ EigenLinSolver)
    endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
    message(STATUS "the OS is Linux")
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    message(STATUS "the OS is Darwin")
endif()

file (GLOB MAIN_FILE
    src/main.cpp
    src/Auxiliary.cpp
    src/Auxiliary.h
    src/PointSampling.cpp
    src/PointSampling.h
    src/MeshCache.cpp
    src/MeshCache.h
    src/MeshCut.cpp
    src/MeshCut.h
    src/KPNewton.cpp
    src/KPNewton.h
    src/PointFinding.cpp
    src/PointFinding.h
    src/GAP.cpp
    src/GAP.h
    src/AddAuxiliaryPoint.cpp
    src/AddAuxiliaryPoint.h
)

file (GLOB SOLVER_FILE
    src/Solver/${SOLVER_}.cpp
    src/Solver/${SOLVER_}.h
    src/Solver/Solver.h
)

add_executable(${PROJECT_NAME} 
    ${MAIN_FILE}
    ${SOLVER_FILE}
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_link_libraries(${PROJECT_NAME} 
    debug OpenMeshToolsd
    debug OpenMeshCored
    optimized OpenMeshTools
    optimized OpenMeshCore)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_link_libraries(${PROJECT_NAME}  OpenMeshTools OpenMeshCore)
endif()

target_link_libraries(${PROJECT_NAME} ${SOLVER_LIBRARY})
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_link_libraries(${PROJECT_NAME} ${OPENMP_LIB_PATH})
endif()

# 设置编译器标志, 这句话要在target_link_libraries后面
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(${PROJECT_NAME} PRIVATE /W4 /WX)# /W4 较高的警告等级; /WX 所有警告都被视为错误 /wd4820 屏蔽C4820警告。
endif()
