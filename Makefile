# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:
.PHONY : .NOTPARALLEL

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
CMAKE_SOURCE_DIR = /home/nicolas/dev/eclipse-workspace/Stage3A

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nicolas/dev/eclipse-workspace/Stage3A

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/nicolas/dev/eclipse-workspace/Stage3A/CMakeFiles /home/nicolas/dev/eclipse-workspace/Stage3A/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/nicolas/dev/eclipse-workspace/Stage3A/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named CancerClustering

# Build rule for target.
CancerClustering: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 CancerClustering
.PHONY : CancerClustering

# fast build rule for target.
CancerClustering/fast:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/build
.PHONY : CancerClustering/fast

src/correlationMatrix.o: src/correlationMatrix.cpp.o
.PHONY : src/correlationMatrix.o

# target to build an object file
src/correlationMatrix.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/correlationMatrix.cpp.o
.PHONY : src/correlationMatrix.cpp.o

src/correlationMatrix.i: src/correlationMatrix.cpp.i
.PHONY : src/correlationMatrix.i

# target to preprocess a source file
src/correlationMatrix.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/correlationMatrix.cpp.i
.PHONY : src/correlationMatrix.cpp.i

src/correlationMatrix.s: src/correlationMatrix.cpp.s
.PHONY : src/correlationMatrix.s

# target to generate assembly for a file
src/correlationMatrix.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/correlationMatrix.cpp.s
.PHONY : src/correlationMatrix.cpp.s

src/dataReader.o: src/dataReader.cpp.o
.PHONY : src/dataReader.o

# target to build an object file
src/dataReader.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/dataReader.cpp.o
.PHONY : src/dataReader.cpp.o

src/dataReader.i: src/dataReader.cpp.i
.PHONY : src/dataReader.i

# target to preprocess a source file
src/dataReader.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/dataReader.cpp.i
.PHONY : src/dataReader.cpp.i

src/dataReader.s: src/dataReader.cpp.s
.PHONY : src/dataReader.s

# target to generate assembly for a file
src/dataReader.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/dataReader.cpp.s
.PHONY : src/dataReader.cpp.s

src/heatMap.o: src/heatMap.cpp.o
.PHONY : src/heatMap.o

# target to build an object file
src/heatMap.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/heatMap.cpp.o
.PHONY : src/heatMap.cpp.o

src/heatMap.i: src/heatMap.cpp.i
.PHONY : src/heatMap.i

# target to preprocess a source file
src/heatMap.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/heatMap.cpp.i
.PHONY : src/heatMap.cpp.i

src/heatMap.s: src/heatMap.cpp.s
.PHONY : src/heatMap.s

# target to generate assembly for a file
src/heatMap.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/heatMap.cpp.s
.PHONY : src/heatMap.cpp.s

src/k_means.o: src/k_means.cpp.o
.PHONY : src/k_means.o

# target to build an object file
src/k_means.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/k_means.cpp.o
.PHONY : src/k_means.cpp.o

src/k_means.i: src/k_means.cpp.i
.PHONY : src/k_means.i

# target to preprocess a source file
src/k_means.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/k_means.cpp.i
.PHONY : src/k_means.cpp.i

src/k_means.s: src/k_means.cpp.s
.PHONY : src/k_means.s

# target to generate assembly for a file
src/k_means.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/k_means.cpp.s
.PHONY : src/k_means.cpp.s

src/lodePNG/lodepng.o: src/lodePNG/lodepng.cpp.o
.PHONY : src/lodePNG/lodepng.o

# target to build an object file
src/lodePNG/lodepng.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/lodePNG/lodepng.cpp.o
.PHONY : src/lodePNG/lodepng.cpp.o

src/lodePNG/lodepng.i: src/lodePNG/lodepng.cpp.i
.PHONY : src/lodePNG/lodepng.i

# target to preprocess a source file
src/lodePNG/lodepng.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/lodePNG/lodepng.cpp.i
.PHONY : src/lodePNG/lodepng.cpp.i

src/lodePNG/lodepng.s: src/lodePNG/lodepng.cpp.s
.PHONY : src/lodePNG/lodepng.s

# target to generate assembly for a file
src/lodePNG/lodepng.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/lodePNG/lodepng.cpp.s
.PHONY : src/lodePNG/lodepng.cpp.s

src/main.o: src/main.cpp.o
.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

src/stats.o: src/stats.cpp.o
.PHONY : src/stats.o

# target to build an object file
src/stats.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/stats.cpp.o
.PHONY : src/stats.cpp.o

src/stats.i: src/stats.cpp.i
.PHONY : src/stats.i

# target to preprocess a source file
src/stats.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/stats.cpp.i
.PHONY : src/stats.cpp.i

src/stats.s: src/stats.cpp.s
.PHONY : src/stats.s

# target to generate assembly for a file
src/stats.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/stats.cpp.s
.PHONY : src/stats.cpp.s

src/tests/correlationMatrix_test.o: src/tests/correlationMatrix_test.cpp.o
.PHONY : src/tests/correlationMatrix_test.o

# target to build an object file
src/tests/correlationMatrix_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/correlationMatrix_test.cpp.o
.PHONY : src/tests/correlationMatrix_test.cpp.o

src/tests/correlationMatrix_test.i: src/tests/correlationMatrix_test.cpp.i
.PHONY : src/tests/correlationMatrix_test.i

# target to preprocess a source file
src/tests/correlationMatrix_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/correlationMatrix_test.cpp.i
.PHONY : src/tests/correlationMatrix_test.cpp.i

src/tests/correlationMatrix_test.s: src/tests/correlationMatrix_test.cpp.s
.PHONY : src/tests/correlationMatrix_test.s

# target to generate assembly for a file
src/tests/correlationMatrix_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/correlationMatrix_test.cpp.s
.PHONY : src/tests/correlationMatrix_test.cpp.s

src/tests/dataReader_test.o: src/tests/dataReader_test.cpp.o
.PHONY : src/tests/dataReader_test.o

# target to build an object file
src/tests/dataReader_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/dataReader_test.cpp.o
.PHONY : src/tests/dataReader_test.cpp.o

src/tests/dataReader_test.i: src/tests/dataReader_test.cpp.i
.PHONY : src/tests/dataReader_test.i

# target to preprocess a source file
src/tests/dataReader_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/dataReader_test.cpp.i
.PHONY : src/tests/dataReader_test.cpp.i

src/tests/dataReader_test.s: src/tests/dataReader_test.cpp.s
.PHONY : src/tests/dataReader_test.s

# target to generate assembly for a file
src/tests/dataReader_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/dataReader_test.cpp.s
.PHONY : src/tests/dataReader_test.cpp.s

src/tests/general_test.o: src/tests/general_test.cpp.o
.PHONY : src/tests/general_test.o

# target to build an object file
src/tests/general_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/general_test.cpp.o
.PHONY : src/tests/general_test.cpp.o

src/tests/general_test.i: src/tests/general_test.cpp.i
.PHONY : src/tests/general_test.i

# target to preprocess a source file
src/tests/general_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/general_test.cpp.i
.PHONY : src/tests/general_test.cpp.i

src/tests/general_test.s: src/tests/general_test.cpp.s
.PHONY : src/tests/general_test.s

# target to generate assembly for a file
src/tests/general_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/general_test.cpp.s
.PHONY : src/tests/general_test.cpp.s

src/tests/k_means_test.o: src/tests/k_means_test.cpp.o
.PHONY : src/tests/k_means_test.o

# target to build an object file
src/tests/k_means_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/k_means_test.cpp.o
.PHONY : src/tests/k_means_test.cpp.o

src/tests/k_means_test.i: src/tests/k_means_test.cpp.i
.PHONY : src/tests/k_means_test.i

# target to preprocess a source file
src/tests/k_means_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/k_means_test.cpp.i
.PHONY : src/tests/k_means_test.cpp.i

src/tests/k_means_test.s: src/tests/k_means_test.cpp.s
.PHONY : src/tests/k_means_test.s

# target to generate assembly for a file
src/tests/k_means_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/k_means_test.cpp.s
.PHONY : src/tests/k_means_test.cpp.s

src/tests/lodePNG_test.o: src/tests/lodePNG_test.cpp.o
.PHONY : src/tests/lodePNG_test.o

# target to build an object file
src/tests/lodePNG_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/lodePNG_test.cpp.o
.PHONY : src/tests/lodePNG_test.cpp.o

src/tests/lodePNG_test.i: src/tests/lodePNG_test.cpp.i
.PHONY : src/tests/lodePNG_test.i

# target to preprocess a source file
src/tests/lodePNG_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/lodePNG_test.cpp.i
.PHONY : src/tests/lodePNG_test.cpp.i

src/tests/lodePNG_test.s: src/tests/lodePNG_test.cpp.s
.PHONY : src/tests/lodePNG_test.s

# target to generate assembly for a file
src/tests/lodePNG_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/lodePNG_test.cpp.s
.PHONY : src/tests/lodePNG_test.cpp.s

src/tests/stats_test.o: src/tests/stats_test.cpp.o
.PHONY : src/tests/stats_test.o

# target to build an object file
src/tests/stats_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/stats_test.cpp.o
.PHONY : src/tests/stats_test.cpp.o

src/tests/stats_test.i: src/tests/stats_test.cpp.i
.PHONY : src/tests/stats_test.i

# target to preprocess a source file
src/tests/stats_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/stats_test.cpp.i
.PHONY : src/tests/stats_test.cpp.i

src/tests/stats_test.s: src/tests/stats_test.cpp.s
.PHONY : src/tests/stats_test.s

# target to generate assembly for a file
src/tests/stats_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/stats_test.cpp.s
.PHONY : src/tests/stats_test.cpp.s

src/tests/unsupervisedNormalization_test.o: src/tests/unsupervisedNormalization_test.cpp.o
.PHONY : src/tests/unsupervisedNormalization_test.o

# target to build an object file
src/tests/unsupervisedNormalization_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/unsupervisedNormalization_test.cpp.o
.PHONY : src/tests/unsupervisedNormalization_test.cpp.o

src/tests/unsupervisedNormalization_test.i: src/tests/unsupervisedNormalization_test.cpp.i
.PHONY : src/tests/unsupervisedNormalization_test.i

# target to preprocess a source file
src/tests/unsupervisedNormalization_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/unsupervisedNormalization_test.cpp.i
.PHONY : src/tests/unsupervisedNormalization_test.cpp.i

src/tests/unsupervisedNormalization_test.s: src/tests/unsupervisedNormalization_test.cpp.s
.PHONY : src/tests/unsupervisedNormalization_test.s

# target to generate assembly for a file
src/tests/unsupervisedNormalization_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/unsupervisedNormalization_test.cpp.s
.PHONY : src/tests/unsupervisedNormalization_test.cpp.s

src/tests/utilities_test.o: src/tests/utilities_test.cpp.o
.PHONY : src/tests/utilities_test.o

# target to build an object file
src/tests/utilities_test.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/utilities_test.cpp.o
.PHONY : src/tests/utilities_test.cpp.o

src/tests/utilities_test.i: src/tests/utilities_test.cpp.i
.PHONY : src/tests/utilities_test.i

# target to preprocess a source file
src/tests/utilities_test.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/utilities_test.cpp.i
.PHONY : src/tests/utilities_test.cpp.i

src/tests/utilities_test.s: src/tests/utilities_test.cpp.s
.PHONY : src/tests/utilities_test.s

# target to generate assembly for a file
src/tests/utilities_test.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/tests/utilities_test.cpp.s
.PHONY : src/tests/utilities_test.cpp.s

src/unsupervisedNormalization.o: src/unsupervisedNormalization.cpp.o
.PHONY : src/unsupervisedNormalization.o

# target to build an object file
src/unsupervisedNormalization.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/unsupervisedNormalization.cpp.o
.PHONY : src/unsupervisedNormalization.cpp.o

src/unsupervisedNormalization.i: src/unsupervisedNormalization.cpp.i
.PHONY : src/unsupervisedNormalization.i

# target to preprocess a source file
src/unsupervisedNormalization.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/unsupervisedNormalization.cpp.i
.PHONY : src/unsupervisedNormalization.cpp.i

src/unsupervisedNormalization.s: src/unsupervisedNormalization.cpp.s
.PHONY : src/unsupervisedNormalization.s

# target to generate assembly for a file
src/unsupervisedNormalization.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/unsupervisedNormalization.cpp.s
.PHONY : src/unsupervisedNormalization.cpp.s

src/utilities.o: src/utilities.cpp.o
.PHONY : src/utilities.o

# target to build an object file
src/utilities.cpp.o:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/utilities.cpp.o
.PHONY : src/utilities.cpp.o

src/utilities.i: src/utilities.cpp.i
.PHONY : src/utilities.i

# target to preprocess a source file
src/utilities.cpp.i:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/utilities.cpp.i
.PHONY : src/utilities.cpp.i

src/utilities.s: src/utilities.cpp.s
.PHONY : src/utilities.s

# target to generate assembly for a file
src/utilities.cpp.s:
	$(MAKE) -f CMakeFiles/CancerClustering.dir/build.make CMakeFiles/CancerClustering.dir/src/utilities.cpp.s
.PHONY : src/utilities.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... CancerClustering"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... src/correlationMatrix.o"
	@echo "... src/correlationMatrix.i"
	@echo "... src/correlationMatrix.s"
	@echo "... src/dataReader.o"
	@echo "... src/dataReader.i"
	@echo "... src/dataReader.s"
	@echo "... src/heatMap.o"
	@echo "... src/heatMap.i"
	@echo "... src/heatMap.s"
	@echo "... src/k_means.o"
	@echo "... src/k_means.i"
	@echo "... src/k_means.s"
	@echo "... src/lodePNG/lodepng.o"
	@echo "... src/lodePNG/lodepng.i"
	@echo "... src/lodePNG/lodepng.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/stats.o"
	@echo "... src/stats.i"
	@echo "... src/stats.s"
	@echo "... src/tests/correlationMatrix_test.o"
	@echo "... src/tests/correlationMatrix_test.i"
	@echo "... src/tests/correlationMatrix_test.s"
	@echo "... src/tests/dataReader_test.o"
	@echo "... src/tests/dataReader_test.i"
	@echo "... src/tests/dataReader_test.s"
	@echo "... src/tests/general_test.o"
	@echo "... src/tests/general_test.i"
	@echo "... src/tests/general_test.s"
	@echo "... src/tests/k_means_test.o"
	@echo "... src/tests/k_means_test.i"
	@echo "... src/tests/k_means_test.s"
	@echo "... src/tests/lodePNG_test.o"
	@echo "... src/tests/lodePNG_test.i"
	@echo "... src/tests/lodePNG_test.s"
	@echo "... src/tests/stats_test.o"
	@echo "... src/tests/stats_test.i"
	@echo "... src/tests/stats_test.s"
	@echo "... src/tests/unsupervisedNormalization_test.o"
	@echo "... src/tests/unsupervisedNormalization_test.i"
	@echo "... src/tests/unsupervisedNormalization_test.s"
	@echo "... src/tests/utilities_test.o"
	@echo "... src/tests/utilities_test.i"
	@echo "... src/tests/utilities_test.s"
	@echo "... src/unsupervisedNormalization.o"
	@echo "... src/unsupervisedNormalization.i"
	@echo "... src/unsupervisedNormalization.s"
	@echo "... src/utilities.o"
	@echo "... src/utilities.i"
	@echo "... src/utilities.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

