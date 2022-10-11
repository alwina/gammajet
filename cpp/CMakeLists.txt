cmake_minimum_required(VERSION 3.10)

# set the project name
project(gammajet)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add ROOT
find_package(ROOT)

# add OpenMP prallelization (MPI currently overkill, but may be necessary later)
find_package(OpenMP)

# RooUnfold is a manual install, and $ROOUNFOLD_DIR is a variable that doesn't naturally exist,
# so add to .bashrc or export it before compiling or whatever
# it only needs to exist for response_matrix

# the situation with yaml-cpp is similar to RooUnfold, but it needs to exist for all the scripts

# include directories
include_directories(${ROOT_INCLUDE_DIRS})
include_directories($ENV{HDF5_DIR}/include)
include_directories($ENV{ROOUNFOLD_DIR}/src)
include_directories($ENV{YAMLCPP}/include)
                    
# library directories
link_directories($ENV{ROOUNFOLD_DIR})
link_directories($ENV{HDF5_DIR}/lib)
link_directories($ENV{YAMLCPP}/build)

# add executables
add_executable(same_event same_event.cc)
add_executable(same_event_mc same_event_mc.cc)
add_executable(same_event_skimmedpp same_event_skimmedpp.cc)
add_executable(same_event_skimmedpp_mc same_event_skimmedpp_mc.cc)
add_executable(centrality_checks centrality_checks.cc)
add_executable(mixed_event mixed_event.cc)
add_executable(process_gjmc process_gjmc.cc)
add_executable(process_gjmc_pp process_gjmc_pp.cc)
add_executable(skim_pp skim_pp.cc)

# link libraries
target_link_libraries(same_event ${ROOT_LIBRARIES} yaml-cpp)
target_link_libraries(same_event_mc ${ROOT_LIBRARIES} yaml-cpp)
target_link_libraries(same_event_skimmedpp ${ROOT_LIBRARIES} yaml-cpp)
target_link_libraries(same_event_skimmedpp_mc ${ROOT_LIBRARIES} yaml-cpp)
target_link_libraries(centrality_checks ${ROOT_LIBRARIES} yaml-cpp)
target_link_libraries(mixed_event ${ROOT_LIBRARIES} hdf5_cpp hdf5 yaml-cpp OpenMP::OpenMP_CXX)
target_link_libraries(process_gjmc RooUnfold ${ROOT_LIBRARIES} yaml-cpp)
target_link_libraries(process_gjmc_pp RooUnfold ${ROOT_LIBRARIES} yaml-cpp)
target_link_libraries(skim_pp ${ROOT_LIBRARIES} yaml-cpp)
