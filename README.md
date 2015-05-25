# Overview

This is a collection of programs to analyze Mark Wise and Yue Zhang's
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

# Algorithm

The basic algorithm can easily be generalized to ODEs with Neumann
boundary conditions. I actually haven't seen this approach taken
anywhere, and can be used for any elliptic problem with spherical
symmetry. Basically, if you have a problem which can be put in the form
`f'(x) = \int^{x}_{0}G(r,f(r))dr`, with `f'(0)=0`, you can use this
trick: `f'(x) = f'(x-h) + \int^{x}_{x-h}G(r,f(r))dr`. The numerical
routine takes advantage of this identity, which apparently Mathematica
cannot do.

# License

Everything is under the MIT License, as noted in the LICENSE file and
the "License" section of the software.
