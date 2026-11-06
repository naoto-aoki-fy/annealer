cmake_minimum_required(VERSION 3.20)
project(annealer VERSION 0.1.0 LANGUAGES C CXX)

add_library(annealer INTERFACE)
add_library(annealer::annealer ALIAS annealer)
target_include_directories(annealer INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(annealer INTERFACE cxx_std_20)
set(CMAKE_CXX_EXTENSIONS OFF)
option(ANNEALER_ENABLE_MPI "Enable distributed MPI multi-start solver" OFF)
if(ANNEALER_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    add_library(annealer_mpi INTERFACE)
    add_library(annealer::mpi ALIAS annealer_mpi)
    target_link_libraries(annealer_mpi INTERFACE annealer::annealer MPI::MPI_CXX)
endif()

option(ANNEALER_ENABLE_IPO "Enable supported native link-time optimization" OFF)
if(ANNEALER_ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error LANGUAGES C CXX)
    if(NOT ipo_supported)
        message(FATAL_ERROR "IPO is unavailable: ${ipo_error}")
    endif()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_executable(multivalued_example
    examples/main.cpp examples/user_model.cpp examples/c_factors.c)
target_link_libraries(multivalued_example PRIVATE annealer::annealer)
set_target_properties(multivalued_example PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED YES)
foreach(demo binary_hubo nonlinear minimal gset_g1 dimacs_coloring maxsat set_covering qaplib tsplib jsplib_jobshop)
    add_executable(${demo}_example examples/${demo}.cpp)
    target_link_libraries(${demo}_example PRIVATE annealer::annealer)
endforeach()

option(ANNEALER_BUILD_BENCHMARKS "Build sequential reproducible benchmarks" OFF)
if(ANNEALER_BUILD_BENCHMARKS)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" bench_build_type)
    set(ANNEALER_BENCH_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${bench_build_type}}")
    configure_file(benchmarks/build_info.hpp.in generated/build_info.hpp @ONLY)
    add_executable(gset_benchmark benchmarks/gset_benchmark.cpp)
    target_include_directories(gset_benchmark PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
    target_link_libraries(gset_benchmark PRIVATE annealer::annealer)
    add_executable(serial_benchmark benchmarks/serial_benchmark.cpp)
    target_include_directories(serial_benchmark PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
    target_link_libraries(serial_benchmark PRIVATE annealer::annealer)
endif()

include(CTest)
if(BUILD_TESTING)
    if(ANNEALER_ENABLE_MPI)
        add_executable(mpi_tests tests/mpi_tests.cpp)
        target_link_libraries(mpi_tests PRIVATE annealer::mpi)
        foreach(ranks 1 2 4)
            add_test(NAME mpi_tests_${ranks}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi_tests> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(mpi_tests_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
        endforeach()
    endif()
    add_executable(jobshop_model_tests tests/jobshop_model_tests.cpp)
    target_link_libraries(jobshop_model_tests PRIVATE annealer::annealer)
    add_test(NAME jobshop_model_tests COMMAND jobshop_model_tests)
    add_test(NAME jobshop_cli_tests COMMAND ${CMAKE_COMMAND}
        "-DEXAMPLE=$<TARGET_FILE:jsplib_jobshop_example>"
        "-DWORK=${CMAKE_CURRENT_BINARY_DIR}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/jobshop_cli_tests.cmake")
    add_executable(jsplib_jobshop_tests tests/jsplib_jobshop_tests.cpp)
    target_compile_features(jsplib_jobshop_tests PRIVATE cxx_std_20)
    add_test(NAME jsplib_jobshop_tests COMMAND jsplib_jobshop_tests)
    add_executable(tsplib_tests tests/tsplib_tests.cpp)
    target_link_libraries(tsplib_tests PRIVATE annealer::annealer)
    add_test(NAME tsplib_tests COMMAND tsplib_tests)
    add_test(NAME tsplib_cli_tests COMMAND ${CMAKE_COMMAND}
        "-DEXAMPLE=$<TARGET_FILE:tsplib_example>"
        "-DWORK=${CMAKE_CURRENT_BINARY_DIR}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/tsplib_cli_tests.cmake")
    add_executable(qaplib_tests tests/qaplib_tests.cpp)
    target_link_libraries(qaplib_tests PRIVATE annealer::annealer)
    add_test(NAME qaplib_tests COMMAND qaplib_tests)
    add_test(NAME qaplib_cli_tests COMMAND ${CMAKE_COMMAND}
        "-DEXAMPLE=$<TARGET_FILE:qaplib_example>"
        "-DWORK=${CMAKE_CURRENT_BINARY_DIR}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/qaplib_cli_tests.cmake")
    add_executable(swap_tests tests/swap_tests.cpp)
    target_link_libraries(swap_tests PRIVATE annealer::annealer)
    add_test(NAME swap_tests COMMAND swap_tests)
    add_executable(set_covering_tests tests/set_covering_tests.cpp)
    target_link_libraries(set_covering_tests PRIVATE annealer::annealer)
    add_test(NAME set_covering_tests COMMAND set_covering_tests)
    add_test(NAME set_covering_cli_tests COMMAND ${CMAKE_COMMAND}
        "-DEXAMPLE=$<TARGET_FILE:set_covering_example>"
        "-DWORK=${CMAKE_CURRENT_BINARY_DIR}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/set_covering_cli_tests.cmake")
    add_executable(maxsat_tests tests/maxsat_tests.cpp)
    target_link_libraries(maxsat_tests PRIVATE annealer::annealer)
    add_test(NAME maxsat_tests COMMAND maxsat_tests)
    add_test(NAME maxsat_cli_tests COMMAND ${CMAKE_COMMAND}
        "-DEXAMPLE=$<TARGET_FILE:maxsat_example>"
        "-DWORK=${CMAKE_CURRENT_BINARY_DIR}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/maxsat_cli_tests.cmake")
    add_executable(coloring_tests tests/coloring_tests.cpp)
    target_link_libraries(coloring_tests PRIVATE annealer::annealer)
    add_test(NAME coloring_tests COMMAND coloring_tests)
    add_test(NAME coloring_cli_tests COMMAND ${CMAKE_COMMAND}
        "-DEXAMPLE=$<TARGET_FILE:dimacs_coloring_example>"
        "-DWORK=${CMAKE_CURRENT_BINARY_DIR}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/coloring_cli_tests.cmake")
    set(gset_cli_benchmark_arg "")
    if(ANNEALER_BUILD_BENCHMARKS)
        set(gset_cli_benchmark_arg "-DBENCHMARK=$<TARGET_FILE:gset_benchmark>")
    endif()
    add_test(NAME gset_cli_tests COMMAND ${CMAKE_COMMAND}
        "-DEXAMPLE=$<TARGET_FILE:gset_g1_example>"
        ${gset_cli_benchmark_arg} "-DWORK=${CMAKE_CURRENT_BINARY_DIR}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/gset_cli_tests.cmake")
    add_executable(instrumentation_tests tests/instrumentation_tests.cpp)
    target_link_libraries(instrumentation_tests PRIVATE annealer::annealer)
    add_test(NAME instrumentation_tests COMMAND instrumentation_tests)
    add_test(NAME gset_cli_error COMMAND gset_g1_example ignored --target-cut invalid)
    set_tests_properties(gset_cli_error PROPERTIES WILL_FAIL TRUE)
    if(ANNEALER_BUILD_BENCHMARKS)
        add_test(NAME gset_benchmark_cli_error COMMAND gset_benchmark ignored --budget 0)
        set_tests_properties(gset_benchmark_cli_error PROPERTIES WILL_FAIL TRUE)
    endif()
    add_executable(gset_tests tests/gset_tests.cpp)
    target_link_libraries(gset_tests PRIVATE annealer::annealer)
    add_test(NAME gset_tests COMMAND gset_tests)
    add_executable(model_tests tests/model_tests.cpp)
    target_link_libraries(model_tests PRIVATE annealer::annealer)
    add_test(NAME model_tests COMMAND model_tests)
    add_executable(execution_tests tests/execution_tests.cpp)
    target_link_libraries(execution_tests PRIVATE annealer::annealer)
    add_test(NAME execution_tests COMMAND execution_tests)
    add_executable(solver_tests tests/solver_tests.cpp)
    target_link_libraries(solver_tests PRIVATE annealer::annealer)
    add_test(NAME solver_tests COMMAND solver_tests)
    add_executable(prepared_tests tests/prepared_tests.cpp)
    target_link_libraries(prepared_tests PRIVATE annealer::annealer)
    add_test(NAME prepared_tests COMMAND prepared_tests)
    add_executable(integration_tests tests/integration_tests.cpp examples/user_model.cpp examples/c_factors.c)
    target_link_libraries(integration_tests PRIVATE annealer::annealer)
    add_test(NAME integration_tests COMMAND integration_tests)
    add_test(NAME multivalued_example COMMAND multivalued_example)
    set_tests_properties(multivalued_example PROPERTIES
        PASS_REGULAR_EXPRESSION "Energy: -18.5")
endif()
