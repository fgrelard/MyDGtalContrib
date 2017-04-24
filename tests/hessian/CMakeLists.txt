SET(DGTAL_CONTRIB_TESTS
        testDerivativeOperator
        testGaussianDerivativeOperator
        testDiscreteHessianFunction
        testRecursiveGaussianDerivativeOperator
        testHessianRecursiveGaussian
        )

FOREACH (FILE ${DGTAL_CONTRIB_TESTS})
    add_executable(${FILE} ${FILE})
    target_link_libraries(${FILE} DGtal DGtalIO)
    set_target_properties(${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_PATH})
ENDFOREACH (FILE)
