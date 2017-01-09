add_subdirectory(clustering)
add_subdirectory(geometry)
add_subdirectory(shapes)

SET(DGTAL_CONTRIB_TESTS
  testStatistics
  testModeller
  )

FOREACH(FILE ${DGTAL_CONTRIB_TESTS})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal DGtalIO)
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_PATH} )
ENDFOREACH(FILE)
