add_subdirectory(clustering)

SET(DGTAL_CONTRIB_TESTS
  testStatistics
  )

FOREACH(FILE ${DGTAL_CONTRIB_TESTS})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal DGtalIO)
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_PATH} )
ENDFOREACH(FILE)
