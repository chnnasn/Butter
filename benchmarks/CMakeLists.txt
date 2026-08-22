foreach(bench_name bench_collision bench_solver bench_2d)
    add_executable(${bench_name} ${bench_name}.cpp)
    target_link_libraries(${bench_name} PRIVATE butter)
endforeach()
