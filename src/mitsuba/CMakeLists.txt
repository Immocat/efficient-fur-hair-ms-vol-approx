set(WMAIN_UTF8_SOURCE "$<$<PLATFORM_ID:Windows>:${CMAKE_CURRENT_SOURCE_DIR}/wmain_utf8.cpp>")

add_executable(mitsuba "${MITSUBA_ROOT_DIR}/include/mitsuba/mitsuba.h" "mitsuba.cpp" ${WMAIN_UTF8_SOURCE})
target_link_libraries(mitsuba getopt mitsuba-python ${XERCES_LIB_LIST} ${GL_LIB_LIST})

add_executable(mtssrv "mtssrv.cpp" ${WMAIN_UTF8_SOURCE})
target_link_libraries(mtssrv getopt mitsuba-python ${XERCES_LIB_LIST} ${GL_LIB_LIST})

add_executable(mtsutil "mtsutil.cpp" ${WMAIN_UTF8_SOURCE})
target_link_libraries(mtsutil getopt mitsuba-python ${XERCES_LIB_LIST} ${GL_LIB_LIST})

