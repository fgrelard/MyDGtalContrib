add_subdirectory(clustering)
add_subdirectory(geometry)
add_subdirectory(shapes)
add_subdirectory(vcm)
add_subdirectory(hessian)

include_directories(${CMAKE_SOURCE_DIR}/)

SET(DGTAL_CONTRIB_TESTS
  testShapeDescriptor
  testModeller
  )

FOREACH(FILE ${DGTAL_CONTRIB_TESTS})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal DGtalIO)
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_PATH} )
ENDFOREACH(FILE)
