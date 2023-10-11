cmake_minimum_required(VERSION 3.16...3.21 FATAL_ERROR)

project(Muography)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)	# src/
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)	# include/

add_executable(main main.cc ${sources} ${headers})
target_link_libraries(main ${Geant4_LIBRARIES})

add_custom_target(Muography DEPENDS main)