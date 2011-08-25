VWPREFIX = /Users/anefian/projects/VisionWorkbench/build
ASPPREFIX = /Users/anefian/projects/StereoPipeline/build
ISISSUPPORTDIR = ${ISISROOT}/../base/noinstall
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

CXXFLAGS_CAM2MAP = -arch i386 -I${ISISSUPPORTDIR}/include/naif  -I${ISISSUPPORTDIR}/include/QtCore -I${ISISSUPPORTDIR}/include/QtGui -I${ISISSUPPORTDIR}/include/QtNetwork -I${ISISSUPPORTDIR}/include/QtSql -I${ISISSUPPORTDIR}/include/QtSvg -I${ISISSUPPORTDIR}/include/QtXml -I${ISISSUPPORTDIR}/include/QtXmlPatterns -DQT_SHARED -I${ISISSUPPORTDIR}/include -I${ISISROOT}/inc
LDFLAGS_CAM2MAP  = -L$(ISISROOT)/3rdParty/lib -L$(ISISSUPPORTDIR)/lib -L$(ISISROOT)/lib -lisis3 -framework vecLib -lgsl -lqwt.5 -lgeos -lcspice -lxerces-c -lkdu_a63R -F$(ISISROOT)/3rdParty/lib -framework vecLib -framework QtCore -framework QtGui -framework QtNetwork -framework QtSql -framework QtSvg -framework QtXml -framework QtXmlPatterns

#needed for 32bit
CXXFLAGS += -arch i386
LDFLAGS += -arch i386

#lidar2image
lidar2img: lidar2img.cc
	g++ -arch i386 -fopenmp $(CXXFLAGS) $(CXXFLAGS_ISIS) $(CXXFLAGS_CAM2MAP) $(INCPATH) $(LIBPATH) $(LDFLAGS) $(LDFLAGS_CAM2MAP) $(OPTIONS) tracks.cc match.cc coregister.cc display.cc weights.cc featuresLOLA.cc util.cc lidar2img.cc map2cam.cc -o lidar2img 

#dem2dem
assembler: assembler.cc
	g++ -arch i386 $(CXXFLAGS) $(INCPATH) $(LIBPATH) $(LDFLAGS) $(OPTIONS)  coregister.cc icp.cc  assembler.cc -o assembler 

#dem2dem
dem2dem: dem2dem.cc
	g++ -arch i386 $(CXXFLAGS) $(INCPATH) $(LIBPATH) $(LDFLAGS) $(OPTIONS)  coregister.cc icp.cc  dem2dem.cc -o dem2dem 

#lidar2dem
lidar2dem: lidar2dem.cc
	g++ -arch i386 $(CXXFLAGS) $(CXXFLAGS_ISIS) $(CXXFLAGS_CAM2MAP) $(INCPATH) $(LIBPATH) $(LDFLAGS)  $(LDFLAGS_CAM2MAP) $(OPTIONS)  icp.cc  coregister.cc tracks.cc  util.cc lidar2dem.cc -o lidar2dem 


clean:
	rm lidar2img
	rm assembler
	rm lidar2dem
	rm dem2dem

