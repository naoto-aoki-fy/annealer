file(WRITE "${WORK}/coloring-triangle.col" "p edge 3 3\ne 1 2\ne 2 3\ne 3 1\n")
foreach(colors 1 3)
    execute_process(COMMAND "${EXAMPLE}" "${WORK}/coloring-triangle.col"
        --colors ${colors} --trials 1 --seed 42 --levels 10 --sweeps 10
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "CLI failed: ${error}")
    endif()
    if(colors EQUAL 1 AND NOT output MATCHES "Conflicting edges: 3.*Best energy: 3.*assignment contains conflicts")
        message(FATAL_ERROR "Incorrect conflict report: ${output}")
    endif()
    if(colors EQUAL 3 AND NOT output MATCHES "Conflicting edges: 0.*Best energy: 0.*valid K-coloring")
        message(FATAL_ERROR "Incorrect coloring report: ${output}")
    endif()
endforeach()
foreach(args "--colors;0" "--colors;-1" "--colors;2x" "--trials;0"
    "--seed;18446744073709551615;--trials;2" "--initial-temperature;nan"
    "--final-temperature;31" "--levels;-1" "--sweeps;bad" "--unknown;1" "--colors")
    execute_process(COMMAND "${EXAMPLE}" "${WORK}/coloring-triangle.col" ${args}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(result EQUAL 0)
        message(FATAL_ERROR "CLI accepted invalid options: ${args}")
    endif()
endforeach()
