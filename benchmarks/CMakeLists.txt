foreach(bench_name bench_collision bench_solver)
    add_executable(${bench_name} ${bench_name}.cpp)
    target_link_libraries(${bench_name} PRIVATE butter)
endforeach()
