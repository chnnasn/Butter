foreach(bench_name bench_collision bench_solver bench_2d bench_2d_solver_quality)
    add_executable(${bench_name} ${bench_name}.cpp)
    target_link_libraries(${bench_name} PRIVATE butter)
endforeach()
