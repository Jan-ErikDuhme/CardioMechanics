# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 4.0

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jd827/Software/CardioMechanics

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jd827/Software/CardioMechanics/_build

# Include any dependencies generated for this target.
include electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/compiler_depend.make

# Include the progress variables for this target.
include electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/progress.make

# Include the compile flags for this target's objects.
include electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/flags.make

electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/codegen:
.PHONY : electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/codegen

electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o: electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/flags.make
electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o: /Users/jd827/Software/CardioMechanics/electrophysiology/src/Lattice/kaVersion.cpp
electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o: electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jd827/Software/CardioMechanics/_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o"
	cd /Users/jd827/Software/CardioMechanics/_build/electrophysiology/src/CellModel && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o -MF CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o.d -o CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o -c /Users/jd827/Software/CardioMechanics/electrophysiology/src/Lattice/kaVersion.cpp

electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.i"
	cd /Users/jd827/Software/CardioMechanics/_build/electrophysiology/src/CellModel && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jd827/Software/CardioMechanics/electrophysiology/src/Lattice/kaVersion.cpp > CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.i

electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.s"
	cd /Users/jd827/Software/CardioMechanics/_build/electrophysiology/src/CellModel && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jd827/Software/CardioMechanics/electrophysiology/src/Lattice/kaVersion.cpp -o CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.s

kaVersion: electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/__/Lattice/kaVersion.cpp.o
kaVersion: electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/build.make
.PHONY : kaVersion

# Rule to build all files generated by this target.
electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/build: kaVersion
.PHONY : electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/build

electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/clean:
	cd /Users/jd827/Software/CardioMechanics/_build/electrophysiology/src/CellModel && $(CMAKE_COMMAND) -P CMakeFiles/kaVersion.dir/cmake_clean.cmake
.PHONY : electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/clean

electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/depend:
	cd /Users/jd827/Software/CardioMechanics/_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jd827/Software/CardioMechanics /Users/jd827/Software/CardioMechanics/electrophysiology/src/CellModel /Users/jd827/Software/CardioMechanics/_build /Users/jd827/Software/CardioMechanics/_build/electrophysiology/src/CellModel /Users/jd827/Software/CardioMechanics/_build/electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : electrophysiology/src/CellModel/CMakeFiles/kaVersion.dir/depend

