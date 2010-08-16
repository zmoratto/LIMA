# The Vision Workbench installation prefix (/usr/local by default)
# VWPREFIX = /usr/local # this is the default - below is what I think is correct dtj June 15th 2010
VWPREFIX = $(HOME)/projects/vision_workbench/build/x86_64_linux_gcc4.1
IRGPKG_DIR = /irg/packages/x86_64_linux_gcc4.1

# If you don’t have libtool installed, you can specify the full
# path to the libtool script in the base directory of your Vision
# Workbench source tree, e.g. $(HOME)/VisionWorkbench-1.0/libtool
LIBTOOL = libtool

# Compilation flags:
# -O3 turns on optimization, which you should almost always do
# -g enables debugging support
# -Wall turns on all compiler warnings
CXXFLAGS = -I$(VWPREFIX)/include -O3 -g -Wall

# Boost:
# The Vision Workbench header files require the boost headers. If
# boost is installed in a non-standord location, you may need
# to uncomment this line and insert the path to the boost headers.
# CXXFLAGS += -I<path to boost include dir>
CXXFLAGS += -I$(IRGPKG_DIR)/boost/include

# Linking flags:
# -lvw includes the Vision Workbench core libraries
LDFLAGS = -L$(VWPREFIX)/lib -lvw
LDFLAGS += -L$(IRGPKG_DIR)/boost/lib
LDFLAGS += -lvwCartography

# Object files:
# List the object files needed to build your program here.
OBJECTS = main.o

# Linking rule:
# Duplicate and modify this rule to build multiple programs.
main: $(OBJECTS)
$(LIBTOOL) --mode=link $(CXX) $(LDFLAGS) -o $@ $^ -g

# Clean-up rule:
clean:
	rm -f *.o *~ \#*
