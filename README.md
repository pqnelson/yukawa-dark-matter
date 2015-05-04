# Overview

This is a collection of programs to analyze Mark Wise and collaborator's
work on modeling dark matter as a Dirac fermion that's Yukawa-coupled to
a scalar field.

The original paper, linked in the CWEB, is at
[arXiv:1411.1772](http://arxiv.org/abs/1411.1772). **This program needs
CWEB, pdftex, and a C++ compiler!** You can set the variables in the
Makefile to your preferences, I use pdftex and the GNU C++ compiler
(because that's all I got).

This is written in CWEB, so it is its own documentation. To compile the
manual, simply type `make doc` in the console. To compile the program to
run, `make bin` will compile the code. If you type `make` or `make all`,
it will compile the executable *and* produce the documentation.

# License

Everything is under the MIT License, as noted in the LICENSE file and
the "License" section of the software.
