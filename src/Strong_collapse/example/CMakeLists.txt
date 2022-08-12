project(Strong_collapse_examples)

add_executable(strong_collapse_from_points strong_collapse_from_points.cpp)
if (TBB_FOUND)
  target_link_libraries(strong_collapse_from_points ${TBB_LIBRARIES})
endif()
add_test(NAME strong_collapse_from_points COMMAND $<TARGET_FILE:strong_collapse_from_points>)

add_executable(strong_collapse_rips_persistence_step_by_step strong_collapse_rips_persistence_step_by_step.cpp)
target_link_libraries(strong_collapse_rips_persistence_step_by_step ${Boost_PROGRAM_OPTIONS_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(strong_collapse_rips_persistence_step_by_step ${TBB_LIBRARIES})
endif()
add_test(NAME strong_collapse_rips_persistence_step_by_step COMMAND $<TARGET_FILE:strong_collapse_rips_persistence_step_by_step>)

add_executable(test_example test_example.cpp)
if (TBB_FOUND)
  target_link_libraries(test_example ${TBB_LIBRARIES})
endif()
