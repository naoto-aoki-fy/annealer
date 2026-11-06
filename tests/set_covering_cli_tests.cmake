file(WRITE "${WORK}/covering-small.txt" "2 3 5 2 8 2 1 2 2 2 3\n")
execute_process(COMMAND "${EXAMPLE}" "${WORK}/covering-small.txt"
    --seed 42 --trials 2 --levels 2 --sweeps 2
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT output MATCHES "Seed: 43" OR
   NOT output MATCHES "Selected-column cost: .*Uncovered rows: .*Total energy: .*Solve seconds: .*Feasible:" OR
   NOT output MATCHES "Best feasible cost:")
    message(FATAL_ERROR "Incorrect CLI report: ${output} ${error}")
endif()
foreach(args "--trials;0" "--seed;18446744073709551615;--trials;2"
    "--initial-temperature;nan" "--final-temperature;31" "--levels;0"
    "--sweeps;bad" "--unknown;1" "--seed")
    execute_process(COMMAND "${EXAMPLE}" "${WORK}/covering-small.txt" ${args}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(result EQUAL 0)
        message(FATAL_ERROR "Accepted invalid options: ${args}")
    endif()
endforeach()
file(WRITE "${WORK}/covering-bad.txt" "1 1 1 0\n")
execute_process(COMMAND "${EXAMPLE}" "${WORK}/covering-bad.txt" RESULT_VARIABLE result)
if(result EQUAL 0)
    message(FATAL_ERROR "Accepted uncoverable row")
endif()
