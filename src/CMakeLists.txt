add_library(histo SHARED histo.c)
install(TARGETS histo DESTINATION lib)
install(FILES histo.h DESTINATION include)

add_executable(histo_test main.c)
target_link_libraries(histo_test histo)
