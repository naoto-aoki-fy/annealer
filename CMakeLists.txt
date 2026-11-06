cmake_minimum_required(VERSION 3.20)
project(annealer VERSION 0.1.0 LANGUAGES C CXX)

add_library(annealer INTERFACE)
add_library(annealer::annealer ALIAS annealer)
target_include_directories(annealer INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(annealer INTERFACE cxx_std_20)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(multivalued_example
    examples/main.cpp examples/user_model.cpp examples/c_factors.c)
target_link_libraries(multivalued_example PRIVATE annealer::annealer)
set_target_properties(multivalued_example PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED YES)

include(CTest)
if(BUILD_TESTING)
    add_executable(model_tests tests/model_tests.cpp)
    target_link_libraries(model_tests PRIVATE annealer::annealer)
    add_test(NAME model_tests COMMAND model_tests)
    add_test(NAME multivalued_example COMMAND multivalued_example)
    set_tests_properties(multivalued_example PROPERTIES
        PASS_REGULAR_EXPRESSION "Energy: -18.5")
endif()
