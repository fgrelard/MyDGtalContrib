SET(DGTAL_CONTRIB_EXAMPLES
  binvox2vol
  pgm2vol
  tre2vol
  vol2noff
  vol2skelab
  dicom2vol
  skelab2vol
  vol2dicom
  vol2pgm
  vol2surfacevol
  )

FOREACH(FILE ${DGTAL_CONTRIB_EXAMPLES})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal DGtalIO ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXAMPLE_OUTPUT_PATH} )
ENDFOREACH(FILE)