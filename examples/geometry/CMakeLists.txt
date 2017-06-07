SET(DGTAL_CONTRIB_EXAMPLES
  EnsureConnectivity
  ShapeModelization
  FlipAxis
  HausdorffDistance
  HausdorffDistanceWithoutJunction
  )

if(WITH_ITK)
  SET(DGTAL_CONTRIB_EXAMPLES
    ${DGTAL_CONTRIB_EXAMPLES}
    DeltaDistance
    DeltaDistanceEdge
    VisuDeltaDistance
    HessianLocalAdaptative
    HessianLocalAdaptative3D
    MultiscaleVesselness
    MultiscaleVesselnessUniqueImage
    MyHessianLocalAdaptative3D
    BarycenterRegistration
    )
  
endif(WITH_ITK)

FOREACH(FILE ${DGTAL_CONTRIB_EXAMPLES})
  add_executable(${FILE} ${FILE})
  target_link_libraries(${FILE} DGtal DGtalIO ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY})
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXAMPLE_OUTPUT_PATH} )
ENDFOREACH(FILE)
