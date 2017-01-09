# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/florent/MyDGtalContrib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/florent/MyDGtalContrib/bin

# Include any dependencies generated for this target.
include tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/depend.make

# Include the progress variables for this target.
include tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/progress.make

# Include the compile flags for this target's objects.
include tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/flags.make

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o: tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/flags.make
tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o: ../tests/utils/geometry/testDigitalPlaneProcessor.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/florent/MyDGtalContrib/bin/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o"
	cd /home/florent/MyDGtalContrib/bin/tests/utils/geometry && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o -c /home/florent/MyDGtalContrib/tests/utils/geometry/testDigitalPlaneProcessor.cpp

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.i"
	cd /home/florent/MyDGtalContrib/bin/tests/utils/geometry && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/florent/MyDGtalContrib/tests/utils/geometry/testDigitalPlaneProcessor.cpp > CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.i

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.s"
	cd /home/florent/MyDGtalContrib/bin/tests/utils/geometry && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/florent/MyDGtalContrib/tests/utils/geometry/testDigitalPlaneProcessor.cpp -o CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.s

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.requires:
.PHONY : tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.requires

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.provides: tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.requires
	$(MAKE) -f tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/build.make tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.provides.build
.PHONY : tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.provides

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.provides.build: tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o

# Object files for target testDigitalPlaneProcessor
testDigitalPlaneProcessor_OBJECTS = \
"CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o"

# External object files for target testDigitalPlaneProcessor
testDigitalPlaneProcessor_EXTERNAL_OBJECTS =

tests/testDigitalPlaneProcessor: tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o
tests/testDigitalPlaneProcessor: tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/build.make
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libmpfr.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libgmp.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libCGAL_Core.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libCGAL.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libboost_thread.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libboost_system.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libpthread.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libDGtalIO.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libDGtal.so
tests/testDigitalPlaneProcessor: /usr/lib/libblas.so
tests/testDigitalPlaneProcessor: /usr/lib/liblapack.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libmpfr.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libgmpxx.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libgmp.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKDeprecated-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKDICOMParser-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOLSM-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOMesh-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKgiftiio-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOCSV-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOHDF5-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOMRC-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKOptimizersv4-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKReview-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOBMP-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOGDCM-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkgdcmMSFF-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkgdcmDICT-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkgdcmIOD-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkgdcmDSED-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkgdcmCommon-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOGIPL-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOJPEG-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOMeta-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIONIFTI-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKniftiio-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKznz-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIONRRD-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKNrrdIO-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOPNG-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkpng-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOTIFF-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitktiff-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkjpeg-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOVTK-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKLabelMap-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKQuadEdgeMesh-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKPolynomials-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKBiasCorrection-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKBioCell-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOBioRad-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOStimulate-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOSpatialObjects-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOXML-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKEXPAT-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKFEM-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKMetaIO-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKOptimizers-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOGE-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOSiemens-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOIPL-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOTransformHDF5-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkhdf5_cpp-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkhdf5-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkzlib-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOTransformInsightLegacy-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOTransformMatlab-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOTransformBase-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKKLMRegionGrowing-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKWatersheds-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKStatistics-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkNetlibSlatec-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKSpatialObjects-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKMesh-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKTransform-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKPath-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkopenjpeg-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKVideoIO-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKIOImageBase-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKVideoCore-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKVtkGlue-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKVTK-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKCommon-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkdouble-conversion-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitksys-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libITKVNLInstantiation-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkvnl_algo-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkvnl-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkv3p_netlib-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitknetlib-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libitkvcl-4.10.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkRenderingFreeTypeOpenGL-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkRenderingOpenGL-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkImagingHybrid-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkIOImage-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkDICOMParser-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkIOCore-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkmetaio-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkpng-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtktiff-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkjpeg-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/lib/libXNVCtrl.a
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libSM.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libICE.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libX11.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libXext.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libXt.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkRenderingFreeType-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkftgl-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkfreetype-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkzlib-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkInteractionStyle-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkRenderingCore-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkFiltersGeometry-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkFiltersExtraction-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkFiltersStatistics-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkImagingFourier-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkalglib-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkFiltersSources-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkFiltersGeneral-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkFiltersCore-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonComputationalGeometry-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkImagingSources-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkImagingCore-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonExecutionModel-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonDataModel-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonMisc-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonSystem-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtksys-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonTransforms-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonMath-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/local/lib/libvtkCommonCore-6.2.so.1
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libcairo.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libQtGui.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libQtXml.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libQtCore.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libQGLViewer.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libCGAL_Core.so
tests/testDigitalPlaneProcessor: /usr/local/lib/libCGAL.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libboost_thread.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libboost_system.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libpthread.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libGLU.so
tests/testDigitalPlaneProcessor: /usr/lib/x86_64-linux-gnu/libGL.so
tests/testDigitalPlaneProcessor: tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../testDigitalPlaneProcessor"
	cd /home/florent/MyDGtalContrib/bin/tests/utils/geometry && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testDigitalPlaneProcessor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/build: tests/testDigitalPlaneProcessor
.PHONY : tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/build

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/requires: tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/testDigitalPlaneProcessor.cpp.o.requires
.PHONY : tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/requires

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/clean:
	cd /home/florent/MyDGtalContrib/bin/tests/utils/geometry && $(CMAKE_COMMAND) -P CMakeFiles/testDigitalPlaneProcessor.dir/cmake_clean.cmake
.PHONY : tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/clean

tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/depend:
	cd /home/florent/MyDGtalContrib/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/florent/MyDGtalContrib /home/florent/MyDGtalContrib/tests/utils/geometry /home/florent/MyDGtalContrib/bin /home/florent/MyDGtalContrib/bin/tests/utils/geometry /home/florent/MyDGtalContrib/bin/tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/utils/geometry/CMakeFiles/testDigitalPlaneProcessor.dir/depend

