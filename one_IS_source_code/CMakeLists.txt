cmake_minimum_required(VERSION 3.3.1)
project(firstproject)
set (CMAKE_CXX_STANDARD 11)


#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/localdisk/hornak/local/share/eigen3/cmake/")

set(Eigen3_DIR "/localdisk/hornak/local/share/eigen3/cmake/")
#set(Eigen3_DIR "/localdisk/hornak/local/EIGEN3/eigen")
#set(Eigen3_DIR "/localdisk/hornak/local/build_dir")


#SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "documentation string here")

set(DCMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE Release)

SET(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp")
set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}  -fPIC ${CMAKE_CXX_FLAGS}")
include_directories(src/)

file(GLOB SOURCES "src/*.cpp"   )

find_package(Eigen3 3.3 REQUIRED)
message("eigen3 include: ${EIGEN3_INCLUDE_DIR}")
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(/localdisk/hornak/local/include)

add_executable(firstproject main.cpp src/MTOC2.cpp src/Microtubule.cpp  src/ISCaptShrinkage.cpp src/KISS.cpp src/Cell.cpp src/Nucleus.cpp src/generalUtilities.cpp src/GeneralGeometry.cpp src/CellShape.cpp src/Cell_parametres.cpp src/ISCorticalSl.cpp src/Surface.cpp src/Dynein.cpp src/mtoc_simple.cpp
src/ideal_mtoc.cpp src/numerical_result.cpp src/Dynein_real.cpp src/mersenne.cpp src/Generator.cpp
)


install(TARGETS firstproject RUNTIME DESTINATION bin)

add_subdirectory(src/)
