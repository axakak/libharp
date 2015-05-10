set(harplib_SRCS c-rbf.cpp traceData.cpp)
add_library(harp ${harplib_SRCS})

#TODO: add compile features check, (i.e constexpr)
#target_compile_features(harp PRIVATE cxx_constexpr)

target_include_directories(harp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(harp ${YAMLCPP_LIBRARIES})
