CPLUSPLUS=g++
CFLAGS=-Wall -O3 --std=c++11
LIB=-lm
TEX=pdftex

all: bin doc
source:
	ctangle free.w - free.cpp
bin: source
	$(CPLUSPLUS) $(CFLAGS) free.cpp $(LIB) -o free
tex:
	cweave free.w
doc: tex
	$(TEX) free.tex
