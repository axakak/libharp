set(harptrain_SRCS crbfTrainer.cpp)
add_executable(harptrain ${harptrain_SRCS})
target_link_libraries(harptrain harp)


set(harpevaluate_SRCS crbfEvaluate.cpp)
add_executable(harpevaluate ${harpevaluate_SRCS})
target_link_libraries(harpevaluate harp)
