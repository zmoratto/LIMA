VWPREFIX = /Users/anefian/projects/VisionWorkbench/build
ASPPREFIX = /Users/anefian/projects/StereoPipeline/build

LIBTOOL = libtool

#VW &ASP
CXXFLAGS = -I$(VWPREFIX)/include -O3 -g -Wall -Wextra -I$(ASPPREFIX)/include\
-Wno-unused-parameter -DHAVE_CONFIG_H -DQT_SQL_LIB -DQT_OPENGL_LIB -DQT_GUI_LIB -DQT_NETWORK_LIB -DQT_CORE_LIB -DQT_SHARED -DDEBUG -fno-strict-aliasing
LDFLAGS = -L$(VWPREFIX)/lib -L$(ASPPREFIX)/lib -lvwCartography -lvwCamera -lvwStereo \
-lvwFileIO -lvwImage -lvwMath -lvwCore -lvwPhotometry -laspIsisIO

INCPATH = -I /opt/local/include
LIBPATH = -L/opt/local/lib/
LIBPATH += -L/usr/lib 
OPTIONS = -lpng -lboost_filesystem-mt -lboost_system-mt -lboost_thread-mt -lboost_program_options-mt

# ASP dependencies
CXXFLAGS_ISIS = -I$(ISISROOT)/../noinstall/include -I$(ISISROOT)/../noinstall/include/QtCore -I$(ISISROOT)/3rdParty/include -I$(ISISROOT)/inc

#needed for 32bit
CXXFLAGS += -arch i386
LDFLAGS += -arch i386


#lidar2image
lidar2img: lidar2img.cc
	g++ -arch i386 -fopenmp $(CXXFLAGS) $(CXXFLAGS_ISIS) $(INCPATH) $(LIBPATH) $(LDFLAGS) $(OPTIONS) coregister.cc tracks.cc match.cc io.cc display.cc weights.cc featuresLOLA.cc util.cc -o lidar2img 
#dem2dem
assembler: assembler.cc
	g++ -arch i386 $(CXXFLAGS) $(INCPATH) $(LIBPATH) $(LDFLAGS) $(OPTIONS)  io.cc icp.cc  assembler.cc -o assembler 

#lidar2dem
lidar2dem: lidar2dem.cc
	g++ -arch i386 $(CXXFLAGS) $(CXXFLAGS_ISIS) $(INCPATH) $(LIBPATH) $(LDFLAGS) $(OPTIONS)  icp.cc  io.cc tracks.cc match.cc display.cc weights.cc featuresLOLA.cc lidar2dem.cc -o lidar2dem 
clean:
	rm coregister
	rm assembler
	rm lidar2dem
        

