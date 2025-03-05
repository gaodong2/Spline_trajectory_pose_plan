# cmake version
cmake_minimum_required(VERSION 3.12.4)

# project
project(optimal_planning)

# cmake config
set(CMAKE_BUILD_TYPE "Debug")
find_package(toppra REQUIRED)
include_directories(${toppra_INCLUDE_DIRS})

include_directories(PUBLIC ./include)


add_definitions(-std=c++11)

add_library(interface ./src/interface_toppra.cpp)
add_library(interface_shared SHARED ./src/interface_toppra.cpp)

add_library(kinematics ./src/kinematics_6axis_ur.cpp)
add_library(kinematics_shared SHARED ./src/kinematics_6axis_ur.cpp)

add_library(B_Spline ./src/B_Spline.cpp)
add_library(B_Spline_shared SHARED ./src/B_Spline.cpp)

add_library(trapezoidal_velocity ./src/trapezoidal_velocity.cpp)
add_library(trapezoidal_velocity_shared SHARED ./src/trapezoidal_velocity.cpp)

add_library(rotational_interpolation ./src/rotational_interpolation.cpp)
add_library(rotational_interpolation_shared SHARED ./src/rotational_interpolation.cpp)

add_library(SphericalQuadraticInterpolation ./src/SphericalQuadraticInterpolation.cpp)
add_library(SphericalQuadraticInterpolation_shared SHARED ./src/SphericalQuadraticInterpolation.cpp)

add_library(MultiPointsTrajectory ./src/MultiPointsTrajectory.cpp)
add_library(MultiPointsTrajectory_shared SHARED ./src/MultiPointsTrajectory.cpp)

add_executable(useHello ./src/main.cpp)
target_link_libraries(useHello interface_shared)
target_link_libraries(useHello kinematics_shared)
target_link_libraries(useHello toppra)
target_link_libraries(useHello B_Spline_shared)
target_link_libraries(useHello trapezoidal_velocity_shared)
target_link_libraries(useHello rotational_interpolation)
target_link_libraries(useHello SphericalQuadraticInterpolation)
target_link_libraries(useHello MultiPointsTrajectory)