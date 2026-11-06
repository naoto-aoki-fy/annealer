file(WRITE "${WORK}/jobshop-small.txt" "1 1 0 2\n")
execute_process(COMMAND "${EXAMPLE}" "${WORK}/jobshop-small.txt" --seed 42 --trials 2 --levels 2 --sweeps 1
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT output MATCHES "Seed: 43 Actual makespan: 2 Constraint violations: 0" OR
    NOT output MATCHES "Total energy: 2.*Feasible: yes")
    message(FATAL_ERROR "Incorrect schedule report: ${output} ${error}")
endif()
foreach(args "--memory-mib;0" "--memory-mib;18446744073709551615" "--trials;0"
    "--seed;18446744073709551615;--trials;2" "--initial-temperature;nan" "--final-temperature;31"
    "--levels;0" "--unknown;1" "--seed")
    execute_process(COMMAND "${EXAMPLE}" "${WORK}/jobshop-small.txt" ${args}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(result EQUAL 0)
        message(FATAL_ERROR "Accepted invalid options or memory budget: ${args}")
    endif()
endforeach()
