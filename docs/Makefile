# GNU Makefile

OS := $(shell uname -s)

ALL: all
all: doc
docs: doc      
doc:
	pdflatex AlignmentTK.tex
	# Use 'open' on Mac systems to view the PDF
ifeq ($(OS),Darwin)
	open AlignmentTK.pdf
endif

clean:
	rm -rf AlignmentTK.aux AlignmentTK.log AlignmentTK.pdf
