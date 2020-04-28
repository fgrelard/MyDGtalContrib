SET(DGTAL_CONTRIB_EXAMPLES
        OrthogonalPlaneEstimation
        OrthogonalPlaneEstimationWithCurve
        PruningSkeletonOrthogonalPlanes
        SkeletonOrthogonalPlanes
        RecenterSkeletonPoints
        VisuSlicesFromSkeleton
        VCMTubularCriterion
        VCMBinaryTubularCriterion
        DeltaVCM
        DeltaDistanceVCMVesselness
        OrthogonalPlaneEvaluation
        Torsion
        )

FOREACH(FILE ${DGTAL_CONTRIB_EXAMPLES})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXAMPLE_OUTPUT_PATH} )
ENDFOREACH(FILE)
