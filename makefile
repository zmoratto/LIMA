VWPREFIX = /Users/anefian/projects/VisionWorkbench/build
ASPPREFIX = /Users/anefian/projects/StereoPipeline/build

LIBTOOL = libtool

#VW &ASP
CXXFLAGS = -I$(VWPREFIX)/include -O3 -g -Wall -Wextra -I$(ASPPREFIX)/include\
-Wno-unused-parameter
LDFLAGS = -L$(VWPREFIX)/lib -L$(ASPPREFIX)/lib -lvwCartography -lvwCamera -lvwStereo \
-lvwFileIO -lvwImage -lvwMath -lvwCore -lvwPhotometry -laspIsisIO

INCPATH = -I /opt/local/include
LIBPATH = -L/opt/local/lib/
LIBPATH += -L/usr/lib 
OPTIONS = -lpng -lboost_filesystem-mt -lboost_system-mt -lboost_thread-mt -lboost_program_options-mt


# ASP dependencies
CXXFLAGS +=  -I$(ISISROOT)/../noinstall/include -I$(ISISROOT)/../noinstall/include/QtCore -I$(ISISROOT)/3rdParty/include -I$(ISISROOT)/inc -DHAVE_CONFIG_H -DQT_SQL_LIB -DQT_OPENGL_LIB -DQT_GUI_LIB -DQT_NETWORK_LIB -DQT_CORE_LIB -DQT_SHARED -DDEBUG -fno-strict-aliasing

#needed for 32bit
CXXFLAGS += -arch i386
LDFLAGS += -arch i386

coregister: coregister.cc
	g++ -arch i386 -fopenmp $(CXXFLAGS) $(INCPATH) $(LIBPATH) $(LDFLAGS) $(OPTIONS) coregister.cc tracks.cc match.cc io.cc display.cc weights.cc featuresLOLA.cc -o coregister 

clean:
	rm coregister
#	rm *.o

