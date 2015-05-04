source:
	ctangle free.w - free.cpp
all: source
	g++ -g -Wall -O3 --std=c++11 free.cpp -o free
tex:
	cweave free.w
doc: tex
	pdftex free.tex
