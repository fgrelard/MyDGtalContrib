add_subdirectory(path)

SET(DGTAL_CONTRIB_TESTS
  testCurve
  testCurveProcessor
  testCurveDecomposition
  testMedialAxis
  testDigitalPlaneProcessor
  testVCMAdjustableRadius
  testVCMOnDigitalSurfaceAdjustableRadius
  testMorphomaths
  testSphericalShellIntersection
  testSetProcessor
  testPointUtil
  testMultiPathThinner
  )

FOREACH(FILE ${DGTAL_CONTRIB_TESTS})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal DGtalIO)
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_PATH} )
ENDFOREACH(FILE)
