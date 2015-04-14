set(harplib_SRCS c-rbf.cpp traceData.cpp)
add_library(libharp ${harplib_SRCS})

#TODO: add compile features, such as constexpr
#target_compile_features(libharp PRIVATE cxx_constexpr)

target_include_directories(libharp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(libharp ${YAMLCPP_LIBRARIES})
