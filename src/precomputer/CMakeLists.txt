add_library(hairvoxelizer SHARED hairvoxelizer.cpp hairvoxelizer.h)
target_link_libraries(hairvoxelizer mitsuba-python bsplinehair)
set_target_properties(hairvoxelizer PROPERTIES FOLDER "plugins")
target_include_directories(hairvoxelizer PRIVATE ${MITSUBA_ROOT_DIR})
target_compile_definitions(hairvoxelizer PRIVATE HairVoxel_EXPORTS)

add_library(hairvoxelizer2 SHARED hairvoxelizer2.cpp hairvoxelizer2.h)
target_link_libraries(hairvoxelizer2 mitsuba-python bsplinehair)
set_target_properties(hairvoxelizer2 PROPERTIES FOLDER "plugins")
target_include_directories(hairvoxelizer2 PRIVATE ${MITSUBA_ROOT_DIR})
target_compile_definitions(hairvoxelizer2 PRIVATE HairVoxelizer2_EXPORTS)

