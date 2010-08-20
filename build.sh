#!/bin/sh

g++ coregister.cc tracks.cc match.cc io.cc display.cc -o coregister -I. -I/opt/local/include  -I /Users/anefian/opencv/cv/include -I /Users/anefian/opencv/cvaux/include -I /Users/anefian/opencv/cxcore/include -I /Users/anefian/opencv/otherlibs/highgui  -I/Users/anefian/projects/VisionWorkbench/build/include -L/opt/local/lib -L/Users/anefian/projects/VisionWorkbench/build/lib -lvwCore -lvwMath -lvwImage -lboost_thread-mt -lboost_program_options-mt -lboost_filesystem-mt -lboost_system-mt -lvwFileIO -lvwCartography -lvwPhotometry -lcv -lcvaux -lcxcore -lhighgui -framework VecLib

