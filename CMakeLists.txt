cmake_minimum_required(VERSION 3.20)

project(xraymcExample4 VERSION 0.0.1 DESCRIPTION "xraymc example 4" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Fetch xraymc library
# xraymc is rather lightweight and is designed to be downloaded at 
# configure time by cmake's FetchContent module

# Optional, make xraymc download material database from IAEA for you
set(XRAYMCLIB_EPICS_DOWNLOAD "ON")
# Optional, include option to save png files
set(XRAYMCLIB_USE_LOADPNG "ON")

# Include CMakes FetchContent module
include(FetchContent)
# Specify where to fint xraymc
FetchContent_Declare(
    libxraymc # name of cmake target (you can chose)
    GIT_REPOSITORY https://github.com/medicalphysics/XRayMClib.git # Git address
    GIT_TAG develop # Git branch name (master for stable version) 
)
FetchContent_MakeAvailable(libxraymc)


# Create an executable application
add_executable(xraymcexample4
	src/main.cpp   
    include/hollowsphere.hpp 
)

target_include_directories(xraymcexample4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(xraymcexample4 PRIVATE libxraymc)
# Adding data from xraymcCarm
xraymclib_add_physics_list(xraymcexample4)
