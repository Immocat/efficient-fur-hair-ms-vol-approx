add_library(homogeneous SHARED maxexp.h homogeneous.cpp)
add_library(heterogeneous SHARED heterogeneous.cpp)
add_library(hairheterogeneous SHARED hairheterogeneous.cpp)

target_link_libraries(homogeneous mitsuba-python)
target_link_libraries(heterogeneous mitsuba-python)
target_link_libraries(hairheterogeneous mitsuba-python)

set_target_properties(homogeneous PROPERTIES FOLDER "plugins")
set_target_properties(heterogeneous PROPERTIES FOLDER "plugins")
set_target_properties(hairheterogeneous PROPERTIES FOLDER "plugins")

add_library(hairmedium SHARED hairmedium.cpp hairmedium.h)
target_link_libraries(hairmedium mitsuba-python)
set_target_properties(hairmedium PROPERTIES FOLDER "plugins")
