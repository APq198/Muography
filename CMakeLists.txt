cmake_minimum_required(VERSION 3.16...3.21 FATAL_ERROR)

project(Muography)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)	# src/
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)	# include/

add_executable(main main.cc ${sources} ${headers})
target_link_libraries(main ${Geant4_LIBRARIES})

file(GLOB MACRO_FILES
  "*.mac"
)
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

file(GLOB BASH_FILES
  "*.sh"
)
file(COPY ${BASH_FILES} DESTINATION ${PROJECT_BINARY_DIR})

add_custom_target(Muography DEPENDS main)