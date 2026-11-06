file(WRITE "${WORK}/tiny.gset" "2 1\n1 2 3\n")
execute_process(COMMAND "${EXAMPLE}" "${WORK}/tiny.gset" --trials 1 --sweeps 0
    --target-cut 4 --trace-csv "${WORK}/trace.csv"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT output MATCHES "Target not reached")
    message(FATAL_ERROR "Example failed: ${output} ${error}")
endif()
file(READ "${WORK}/trace.csv" csv)
if(NOT csv MATCHES "kind,levels,sweeps" OR NOT csv MATCHES "history,200,0,30,")
    message(FATAL_ERROR "Unexpected history CSV: ${csv}")
endif()
execute_process(COMMAND "${EXAMPLE}" "${WORK}/tiny.gset" --trials 1 --levels 0 --target-cut -1
    RESULT_VARIABLE result OUTPUT_VARIABLE output)
if(NOT result EQUAL 0 OR NOT output MATCHES "Updates: 0")
    message(FATAL_ERROR "Initial target/zero levels failed: ${output}")
endif()
foreach(option IN ITEMS --target-cut --trace-csv --seed --unknown)
    execute_process(COMMAND "${EXAMPLE}" "${WORK}/tiny.gset" "${option}"
        RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
    if(result EQUAL 0)
        message(FATAL_ERROR "Missing value was accepted: ${option}")
    endif()
endforeach()
if(DEFINED BENCHMARK)
    execute_process(COMMAND "${BENCHMARK}" "${WORK}/tiny.gset" --levels 0 --budget 0.001
        --target-cut 3 --csv "${WORK}/budget.csv"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(NOT result EQUAL 0 OR NOT output MATCHES "Verified: yes")
        message(FATAL_ERROR "Benchmark failed: ${output} ${error}")
    endif()
    foreach(suffix IN ITEMS "" ".iterations.csv" ".summary.csv" ".metadata.csv")
        file(SIZE "${WORK}/budget.csv${suffix}" size)
        if(size LESS 10)
            message(FATAL_ERROR "Missing benchmark CSV content: ${suffix}")
        endif()
    endforeach()
    foreach(levels IN ITEMS "50," ",100" "x")
        execute_process(COMMAND "${BENCHMARK}" "${WORK}/tiny.gset" --levels "${levels}"
            RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
        if(result EQUAL 0)
            message(FATAL_ERROR "Invalid level list accepted: ${levels}")
        endif()
    endforeach()
endif()

execute_process(COMMAND "${EXAMPLE}" "${WORK}/tiny.gset" --mode serial --trials 1
    --levels 2 --sweeps 2 --profile --target-cut -1 --trace-csv "${WORK}/profile-trace.csv"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT output MATCHES "Profile seconds:" OR
   NOT output MATCHES "TTS seconds:" OR NOT output MATCHES "Verified: yes")
    message(FATAL_ERROR "Serial instrumentation failed: ${output} ${error}")
endif()
