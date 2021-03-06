SET(DGTAL_CONTRIB_EXAMPLES
  EnsureConnectivity
  EnsureConnectivity2D
  ShapeModelization
  FlipAxis
  HausdorffDistance
  HausdorffDistanceWithoutJunction
  GeodesicLoops
  )

if(WITH_ITK)
  SET(DGTAL_CONTRIB_EXAMPLES
    ${DGTAL_CONTRIB_EXAMPLES}
    DeltaDistance
    DeltaDistanceEdge
    DeltaDistanceMax
    VisuDeltaDistance
    HessianLocalAdaptative
    HessianLocalAdaptative3D
    vcm-curvature
    # MultiscaleVesselness
    # MultiscaleVesselnessUniqueImage
    MyHessianLocalAdaptative3D
    BarycenterRegistration
    ExportToGraph
    SimpleThinning
    JaccardIndex
    ImageDifference
    ComparisonImages
    CannyEdgeDetector)
endif(WITH_ITK)

FOREACH(FILE ${DGTAL_CONTRIB_EXAMPLES})
  add_executable(${FILE} ${FILE})
  target_link_libraries(${FILE} DGtal ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY})
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXAMPLE_OUTPUT_PATH} )
ENDFOREACH(FILE)
