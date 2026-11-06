file(WRITE "${WORK}/tsp-small.dat" "TYPE: TSP\nDIMENSION: 2\nEDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION\n1 0 0\n2 3 4\n")
execute_process(COMMAND "${EXAMPLE}" "${WORK}/tsp-small.dat" --reference-cost 10
    --seed 42 --trials 2 --levels 2 --sweeps 2
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0 OR NOT output MATCHES "Seed: 43" OR
   NOT output MATCHES "Tour length: 10 Collisions: 0 Total energy: 10.*Valid tour: yes" OR
   NOT output MATCHES "Reference cost: 10 Objective gap: 0")
    message(FATAL_ERROR "Incorrect CLI output: ${output} ${error}")
endif()
foreach(args "--reference-cost;-1" "--reference-cost;nan" "--reference-cost;9007199254740993"
    "--seed;18446744073709551615;--trials;2" "--trials;0" "--levels;0"
    "--initial-temperature;nan" "--final-temperature;31" "--unknown;1" "--seed")
    execute_process(COMMAND "${EXAMPLE}" "${WORK}/tsp-small.dat" ${args}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(result EQUAL 0)
        message(FATAL_ERROR "Accepted invalid options: ${args}")
    endif()
endforeach()
