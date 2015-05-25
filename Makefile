CPLUSPLUS=g++
CFLAGS=-Wall -O3 --std=c++11
LIB=-lm
TEX=pdftex
FILE=yukawa

all: bin doc
source:
	ctangle $(FILE).w - $(FILE).cpp
bin: source
	$(CPLUSPLUS) $(CFLAGS) $(FILE).cpp $(LIB) -o $(FILE)
tex:
	cweave $(FILE).w
doc: tex
	$(TEX) $(FILE).tex
