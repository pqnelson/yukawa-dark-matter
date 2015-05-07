\def\arcsinh{\mathop{\rm ArcSinh}\nolimits}
\def\title{\uppercase{Dark Matter with Free Yukawa Scalar Field}}
\def\abs#1{|#1|}
\newcount\eqcount
\eqcount=0
\def\vec#1{{\bf #1}}
\def\slug{\hbox{\kern1.5pt\vrule width2.5pt height6pt depth1.5pt\kern1.5pt}}
\def\slugonright{\vrule width0pt\nobreak\hfill\slug}

\def\eqn{{\global\advance\eqcount by1}\eqno{(\the\eqcount)}}
\input macros

@* Include Files.
Before expounding anything, we need to first note some parameters and
headers. This is boring, skip to the next page.

@d _USE_MATH_DEFINES true
@d PI M_PI
@d PI_4 M_PI_4

@c
#include <ctime>
#include <cstdlib> /* for rand(), used in unit tests */
#include <stdexcept>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "types.hpp"
#include "newton.hpp"

@ We will be, throughout, using the maximum precision available on the
current machine. For most machines (using Intel or AMD processors),
that's 20 digits of precision. We also load the mathematical constants
from the C library.

@ @(types.hpp@>=
#ifndef __TYPES_HPP__
#define __TYPES_HPP__
typedef long double real;
typedef std::size_t index;
#endif 

@ {\bf License.}
Everything is under the MIT License (i.e., the ``please don't sue me if
this breaks your computer'' license), quoted here for pedantry:

\medbreak
{\narrower
The MIT License (MIT)

Copyright (c) 2015 by Alex Nelson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the ``Software''), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED ``AS IS'', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
\par}

@* Introduction, Review.
Following \arXiv{1411.1772}, we consider a fermionic field $\chi$
coupled to a scalar field $\phi$. Initially we consider $\phi$ massless
and non-self-interacting. For simplicity we assume spherical
symmetry. 

@ {\bf Lagrangian.} The system's Lagrangian density is
$$
{\cal L} = i\bar{\chi}\gamma^{\mu}\partial_{\mu}\chi
 - m_{\chi}\bar{\chi}\chi - g_{\chi}\bar{\chi}\chi\phi
 + {1\over 2}\partial\phi\cdot\partial\phi.\eqn{}
$$
We assume $g_{\chi}>0$.
This is hard, so instead of considering a fermionic field $\chi$, we
just consider a finite number $N$ of relativistic particles. 

The Lagrangian we work with instead is (setting $\partial_{t}\phi=0$ and
working in mostly-minuses convention)
$$
L = -\sum_{i}m(\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}-{1\over2}\int\nabla\phi\nabla\phi\,{\rm d}^{3}x.\eqn{}
$$
We have introduced the so-called ``scale-invariant mass''
$$
m(\vec{x}_{i}) = m_{\chi} + g_{\chi}\phi(\vec{x}_{i})\eqn{}
$$
where $\vec{x}_{i}(t)$ is the coordinate of the $i$-th $\chi$ particle.

@ {\bf Equations of Motion for Scalar Field.}
The equations of motion from such a Lagrangian for $\phi$ may be
written down:
$$
\nabla\cdot{\partial L\over {\partial\nabla\phi}}-{\partial L\over{\partial\phi}}=0.
$$
We see
$$
\nabla\cdot{\partial L\over {\partial\nabla\phi(\vec{x})}}=\nabla^{2}\phi(\vec{x})
$$
and
$$
{\partial L\over{\partial\phi(\vec{x})}}=\sum_{i}g_{\chi}\delta^{(3)}(\vec{x}-\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}.
$$
Thus we deduce the equations of motion for $\phi$:
$$
\nabla^{2}\phi(\vec{x})
=\sum_{i}g_{\chi}\delta^{(3)}(\vec{x}-\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}.\eqn{}
$$

@ {\bf Momentum.} The conjugate momentum for $\vec{x}_{i}(t)$ is quite simple:
$$
\vec{p}_{i} = m(\vec{x}_{i}){\dot{\vec{x}}_{i}\over\sqrt{1-\dot{\vec{x}}_{i}^{2}}}.\eqn{}
$$
Observe
$$
\vec{p}_{i}^{2} + m(\vec{x}_{i})^{2} = {m(\vec{x}_{i})^{2}\over{1-\dot{\vec{x}}_{i}^{2}}},\eqn{}
$$
and hence
$$
\sqrt{1-\dot{\vec{x}}_{i}^{2}} = {m(\vec{x}_{i})\over{\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}}}.\eqn{}
$$

@ Observe the field equations (4) change due to (7) as
$$
\nabla^{2}\phi(\vec{x})
=g_{\chi}\sum_{i}\delta^{3}(\vec{x}-\vec{x}_{i}) {m(\vec{x}_{i})\over
\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}}.
$$
This is a nonlinear differential equation we're trying to solve.

@ {\bf Hamiltonian.} We find the Hamiltonian from the Legendre transform as
$$
H = \sum_{i}\vec{p}_{i}\dot{\vec{x}}_{i} - L. \eqn{}
$$
We see
$$
\vec{p}_{i}\dot{\vec{x}}_{i} 
= m(\vec{x}_{i})
  {{\dot{\vec{x}}_{i}^{2}}\over{1-\dot{\vec{x}}_{i}^{2}}}
  \sqrt{1-\dot{\vec{x}}_{i}^{2}}.\eqn{}
$$
Hence
$$
\vec{p}_{i}\dot{\vec{x}}_{i} + m(\vec{x}_{i})\sqrt{1-\vec{\dot{x}}_{i}^{2}}
= m(\vec{x}_{i})
  {1\over{1-\dot{\vec{x}}_{i}^{2}}}\sqrt{1-\dot{\vec{x}}_{i}^{2}}.\eqn{}
$$
Therefore the kinetic term Legendre transforms into
$$
\vec{p}_{i}\cdot\dot{\vec{x}}_{i}
 + m(\vec{x}_{i})\sqrt{1-\vec{\dot{x}}_{i}^{2}}
=\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}.\eqn{}
$$
We can write
$$
H = \sum_{i}\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}} 
+ {1\over 2}\int\nabla\phi\cdot\nabla\phi\,{\rm d}^{3}x.\eqn{}
$$
Integration by parts
$$
H = \sum_{i}\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}
- {1\over 2}\int\phi\nabla^{2}\phi\,{\rm d}^{3}x.\eqn{}
$$
Then inserting the equations of motion
$$
H = \sum_{i}\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}} 
- {g_{\chi}\over 2}\sum_{i}\phi(\vec{x}_{i})
  {m(\vec{x}_{i})\over{\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}}}.\eqn{}
$$

@ {\it Remark. Nonzero Scalar Potential Term.}
If we had some nonzero potential $V(\phi)$ for the scalar field, how
would we change things? Well, we would have to modify Eq (2) as
$$
L = -\sum_{i}m(\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}+\int{-1\over2}\nabla\phi\nabla\phi-V(\phi)\,{\rm d}^{3}x
\eqno{(2')}
$$
This would produce an extra term in the $\partial L/\partial\phi$
equation, modifying Eq (4) as
$$
\nabla^{2}\phi(\vec{x})
=\sum_{i}g_{\chi}\delta^{(3)}(\vec{x}-\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}
-{\partial V(\phi)\over{\partial\phi(\vec{x})}}.
\eqno{(4')}
$$
This does not alter the definition of $\vec{p}_{i}$, but the Hamiltonian
changes in several ways: first, since we inserted the equations of
motion into (14), and second since the Lagrangian gains an extra term
for the potential $V(\phi)$. We adjust our Hamiltonian to become
$$
H = \sum_{i}\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}} 
- {g_{\chi}\over 2}\sum_{i}\phi(\vec{x}_{i})
  {m(\vec{x}_{i})\over{\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}}}
+\int\left(V(\phi)+{1\over 2}\phi(\vec{x}){\partial V(\phi)\over{\partial\phi(\vec{x})}}\right)
{\rm d}^{3}x.
\eqno{(14')}
$$

@ {\bf Continuum Limit.}
When $N\gg 1$, we may treat the collection of particles as a continuous
``field'', at least trading the sum over particles for integrals over
the phase space:
$$
\sum_{i}\to\int{\rm d}^{3}r\int{{{\rm d}^{3}p}\over{(2\pi)^{3}}}
f(\vec{r}, \vec{p}).\eqn{}
$$
As always, we assume spherical symmetry. We have the Fermi
momentum\footnote{${}^{1}$}{I honestly have not found a good reference
on what the ``Fermi momentum'' {\it is}. As best as I can tell, it's the
momenta for the ``most energized'' particle state in some configuration
of particles.} $p_{F}(r)$ which may depend on position, and the $\chi$
particles must be confined to a sphere of radius $R$. Hence we conclude
$$
f(\vec{r}, \vec{p}) = 2\theta(R-r)\theta(p_{F}(r)-p)\eqn{}
$$
where $\theta(-)$ is the Heaviside step function.

@ {\bf Number of $\chi$ Particles.}
We may consider the number of $\chi$ particles as
$$
N = \sum_{i}1 = {4\pi\over 3}\int^{R}_{0} p_{F}(r)^{3}r^{2}\,{\rm d}r.\eqn{}
$$
This parameter will be fixed, and given by the user.

@ {\bf Energy's Continuum Limit.}
So, we consider the Hamiltonian (14) and how it changes under the
continuum limit. We see
$$
\sum_{i}\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}
=\int^{R}_{0}r^{2}\int^{p_{F}(r)}_{0}\sqrt{p(r)^{2}+m(r)^{2}}p(r)^{2}\,{\rm
d}p\,{\rm d}r.\eqn{}
$$
A change of variables to $u=p/m$ gives us $m\,{\rm d}u = {\rm
d}p$. Hence we introduce the auxiliary functions
$$
h(z) = \int^{z}_{0}u^{2}\sqrt{1+u^{2}}{\rm d}u = {1\over 4}\left(i(z)+z^{3}\sqrt{1+z^{2}}\right)\eqn{}
$$
$$
i(z) = \int^{z}_{0}{u^{2}\over\sqrt{1+u^{2}}}\,{\rm d}u
     = {1\over 2}\left(z\sqrt{1+z^{2}} - \arcsinh(z)\right)\eqn{}
$$
to write 
$$
\sum_{i}\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}
=\int^{R}_{0}r^{2}m(r)^{4}h\bigl(p_{F}(r)/\abs{m(r)}\bigr)\,{\rm d}r.\eqn{}
$$
Likewise
$$
\sum_{i}\phi(\vec{x}_{i}){m(\vec{x}_{i})\over\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}}
=\int^{R}_{0}r^{2}m(r)^{3}\phi(r)i\bigl(p_{F}(r)/\abs{m(r)}\bigr)\,{\rm d}r.\eqn{}
$$
Hence we may write the total energy as
$$
E_{0}(R) = {4\pi\over 3}
\int^{R}_{0}r^{2}\left[
m(r)^{4}h\bigl(p_{F}(r)/\abs{m(r)}\bigr)
-{g_{\chi}\over 2}m(r)^{3}\phi(r)i\bigl(p_{F}(r)/\abs{m(r)}\bigr)
\right]{\rm d}r.\eqn{}
$$
If we added a nonzero $V(\phi)$ potential, then we'll need to add
$$
E_{{\rm int}}(R) = {4\pi\over 3}\int^{R}_{0}
r^{2}\left[
  V(\phi)+{1\over 2}\phi(r){\partial V(\phi)\over{\partial\phi(r)}}
\right]{\rm d}r\eqn{}
$$
to the total energy function.

@ @c real i(real z) {
  return 0.5*(z*hypot(1.0,z) - asinhl(z));
}

@ @c
/* "h" is used for our step size, */
/* so the h auxiliary function should be called something else */
real hUtil(real z) {
  real u = z*z;
  return 0.25*(i(z) + z*u*hypot(1.0,z));
}

@ {\bf Continuum Limit of Equations of Motion.}
Insider the nugget, when $r<R$, the equations of motion for $\phi$
changes to
$$
\nabla^{2}\phi(r) = {{{\rm d}^{2}}\over {{\rm d}r^{2}}}\phi(r)+{2\over
r}{{{\rm d}}\over {{\rm d}r}}\phi(r) =
{g_{\chi}\over\pi^{2}}m(r)^{3}i(p_{F}(r)/|m(r)|)-{\partial V(\phi)\over{\partial\phi(r)}}\eqn{}
$$
This is great, it's the equation we will be solving, but we have one
glaring problem: we never specified $p_{F}(r)$.

@ {\bf Momentum Ansatz.}
We expect the Fermi momentum to vanish at the surface of the nugget. So
$p_{F}(R)=0$ is a reasonable expectation. A one-parameter family of such
momenta (parametrized by $a\geq0$) would be:
$$
p_{F}(r) = \left({{3\pi\Gamma(4+3a)N}\over{8\Gamma(1+3a)}}\right)^{1/3}
\left(1-{r\over R}\right)^{a}{1\over R}.\eqn{}
$$
We will either set $a=0$ or $a=1/2$. Observe
$$
{\Gamma(4+3a)\over{\Gamma(1+3a)}}=3(1+a)(2+3a)(1+3a)
$$
hence
$$
p_{F}(r) = \left({{9\pi(1+a)(2+3a)(1+3a)N}\over{8}}\right)^{1/3}
\left(1-{r\over R}\right)^{a}{1\over R}.\eqn{}
$$
We note $9\pi/8 = (9/2)(\pi/4)$ for implementation purposes.

@ @c
real NewtonSolver::pF(index j, real a) {
  if (j >= length-1) return 0.0;
  real r = j*h;
  if (r >= R) return 0.0; /* avoid complex numbers! */
  return cbrt(4.5*N*PI_4*(1.0+a)*(2.0+3.0*a)*(1.0+3.0*a))*pow(1.0-(r/R),a)/R;
}


@* Summary of Algorithms, Big Picture.
The basic idea is to rewrite Eq $(25)$ as something of the form
$$G(\phi) = \nabla^{2}\phi-\sigma(\phi)$$
then try to find the $\phi$ which makes this function vanish. For
numerical purposes, we approximate $\nabla^{2}$ using the method of
finite differences. Then the problem naturally becomes a candidate for
Newton--Raphson iterations.

@ {\bf Problems, Solutions.}
It turns out the Jacobian matrix is {\it not} diagonal dominant, and due
to the boundary conditions {\it not} symmetric. This is quite
distressing from a numerical analyst's perspective, but there's still
hope! Instead of solving the system
$$
J(\phi^{(n)})\phi^{(n+1)} = \vec{b}(\phi^{(n)})\eqn{}
$$
we instead solve
$$
J(\phi^{(n)})^{T}J(\phi^{(n)})\phi^{(n+1)} = J(\phi^{(n)})^{T}\vec{b}(\phi^{(n)})\eqn{}
$$
where $J(\phi^{(n)})^{T}$ is the transpose of the Jacobian. We then get
a symmetric tridiagonal matrix on the left hand side.

@ {\bf Theorem.} {\it
The symmetric matrix $J^{T}J$ formed from the Jacobian matrix and its
transpose is positive-definite.}

\medbreak
\noindent{\it Proof.}
By the $LDL^{T}$ decomposition, we see $J=L^{T}$, and $D$ is
the identity matrix. Hence the desired matrix is, in fact, symmetric
positive-definite. \slugonright

\medbreak
(A simpler proof might have pointed out since the Jacobian is
invertible, i.e.\ ``regular'', then $JJ^{T}$ is positive-definite, and
hence $J^{T}J$ is positive-definite too.)

@ @(newton.hpp@>=
#include "types.hpp"
#ifndef __NEWTON_HPP__
#define __NEWTON_HPP__

class NewtonSolver {
private: @/
  @<NewtonSolver Fields@>@;
public: @/
  NewtonSolver(index length_, real R_, real alpha_, real massChi_, real N_, real massPhi_=0.0);
  ~NewtonSolver();
  real pF(index j, real a=0.50); /* Fermi Momentum */
  real energyDensity(real r, real phi);
  real energy();
  real source(index j);
  real partialSource(index j);
  real Laplacian(index i, index j);
  real BoundaryConditionAtOrigin(index i);
  real SurfaceBoundaryCondition(index j);
  real JacobianMatrix(index i, index j);
  void init();
  real a(index i);
  real b(index i);
  real c(index i);
  real d(index i);
  void CroutFactorization();
  void ThomasInvert(); /* Thomas' Algorithm */
  real residual();
  real norm();
  bool isDiagonalDominant();
  bool isSymmetric(); 
  void setR(real R_) { R = R_; }
  real currentR() const { return R; }
  unsigned iterateCount() const { return iterateCounter; }
  real* solution() const { return nextPhi; }
};@;
#endif

@ @<NewtonSolver Fields@>=
  index length;
  real R;
  real alpha; 
  real massChi;
  real N; /* number of fermions */
  real h; /* step size */
  real *phi;
  real *nextPhi;
  real *cPrime; /* to speed up inverting the Jacobian */
  real massPhi;
  real energy_;
  unsigned iterateCounter;

@ @c
NewtonSolver::NewtonSolver(index length_, real R_, real alpha_, real massChi_, real N_, real massPhi_) : 
  length(length_),@#
  R(R_),@#
  alpha(alpha_),@#
  massChi(massChi_),@#
  N(N_),@#
  h(R_/length_),@#
  phi(new real[length_]),@#
  nextPhi(new real[length_]),@#
  cPrime(new real[length_]),@#
  massPhi(massPhi_),@#
  energy_(0.0),@#
  iterateCounter(0) @# {
  index j;
  for (j=0; j<length; j++) {
    cPrime[j] = 0.0;
    phi[j] = nextPhi[j] = -1.0; 
  }
}

NewtonSolver::~NewtonSolver() {
  delete phi;
  delete nextPhi;
  delete cPrime;
}

@* Energy Density. 
Recall (\S\S13--17) the energy for the system will be
$$
E(N,R) = E_{0}(N,R) + E_{{\rm int}}(N,R)\eqn{}
$$
where we had $E_{0}(R)$ as
$$
E_{0}(R) = {4\pi\over 3}
\int^{R}_{0}r^{2}\left[
m(r)^{4}h\bigl(p_{F}(r)/\abs{m(r)}\bigr)
-{g_{\chi}\over 2}m(r)^{3}\phi(r)i\bigl(p_{F}(r)/\abs{m(r)}\bigr)
\right]{\rm d}r.\eqno{(23)}
$$
For some nonzero potential, we have
$$
E_{{\rm int}}(R) = {4\pi\over 3}\int^{R}_{0}
r^{2}\left[
  V(\phi)+{1\over 2}\phi(r){\partial V(\phi)\over{\partial\phi(r)}}
\right]{\rm d}r\eqno{(24)}
$$
We refer to the integrand of $E(N, R)$ as the ``Energy Density''.

@ Observe the energy density has extrema at $\partial_{R}E_{0}(N,R)=0$,
which would imply $p_{F}(R)=0$. This is how we get the boundary
condition on the Fermi momentum. But, if we add a mass to the scalar
field or a self-interacting contribution, then we change the boundary
condition $p_{F}(R)$ to something else.

{\bf To do:} Figure out $p_{F}(R)$ by finding
$\partial_{R}E_{0}(N,R)+V(\phi(R))=0$ where $V(\phi)$ is the potential
term for the scalar field from the Lagrangian density.

{\bf To Do:} Need to unit test the energy somehow...

@ @c
real NewtonSolver::energyDensity(real r, real phi) {
  real g = sqrt(4.0*PI*alpha);
  real scalingMass = massChi - g*phi;
  if (scalingMass == 0.0) throw std::logic_error("Divide by Zero");
  real p = pF(r)/fabs(scalingMass);
  real potentialTerm = pow(phi*massPhi, 2.0); 
  return (r*r*4.0*PI/3.0)*(potentialTerm + pow(scalingMass, 3.0)*(scalingMass*hUtil(p) - 0.5*g*phi*i(p)));
}

@ The energy will be computed numerically via Simpson's rule. This
means, if ${\cal E}_{j}$ is the energy density at node $j$, then
$$
E\approx{h\over 3}\sum_{j=0}^{N/2}{\cal E}_{2j-2}+4{\cal E}_{2j-1}+{\cal E}_{2j}.
$$
The error is $\bigO{h^{4}}$. We also cache the energy for speed.

{\bf To Do:} Consider using adaptive quadrature for performance. There's
no point in useless computation, might as well try to minimize it if at
all possible.

@ @c
real NewtonSolver::energy() {
  if (energy_==0.0) {
    real e = 0.0;
    index j;
    for(j = 1; j<length/2; j++) {
      e += 2*energyDensity(2*j*h, phi[2*j-2]);
      e += 4*energyDensity((2*j-1)*h, phi[2*j-1]);
    }
    energy_=(energyDensity(0.0, phi[0]) + energyDensity((length-1)*h, phi[length-1])+e)*h/3.0;
  }
  return energy_;
}

@* Laplacian. The Laplacian operator is a tad bit subtle, but in radial
coordinates it becomes
$$\nabla^{2}f(r) 
= {1\over r^{2}}{{\rm d}\over{{\rm d}r}}
  \left(r^{2}{{\rm d}\over{{\rm d}r}}\right)f(r)
={{\rm d}^{2}\over{{\rm d}r^{2}}}f(r)+{2\over r}f'(r)\eqn{}$$
The finite difference approximation gives us
$$f''(r)\mapsto{{f_{n+1}-2f_{n}+f_{n-1}}\over h^{2}}\eqn{}$$
and
$$f'(r)\mapsto{{f_{n+1}-f_{n-1}}\over{2h}}.\eqn{}$$
Adding these equations together gives
$$
f''(r) + {2\over r}f'(r)\mapsto
\left({1\over h^{2}}+{1\over{rh}}\right)f_{n+1}
-{2\over h^{2}}f_{n}+
\left({1\over {h^{2}}}-{{1}\over{hr}}\right)f_{n-1}.\eqn{}
$$
Hence, if $r_{n}=nh$, we have
$$\nabla^{2}f(r)\mapsto\left({n+1\over n}{1\over h^{2}}\right)f_{n+1}
-{2\over h^{2}}f_{n}
+\left({n-1\over n}{1\over h^{2}}\right)f_{n-1}.\eqn{}$$
The accuracy is, for the double derivative, of order $h^{4}f^{(4)}(c)/12
+ \bigO{h^{5}}$ for some $c\in[x-h,x+h]$; for the single derivative,
it is of order $h^{3}f^{(3)}(c)/3 + \bigO{h^{5}}$.

Numerically, we will delay dividing by $h^{2}$ until ``the last possible
minute''. This will allow us to take advantage of theorem 20 more easily.

@ @c
real NewtonSolver::Laplacian(index i, index j) {
  if (i == j)
    return -2.0;
  real r = i*h;
  if (i+1 == j)
    return 1.0 + (h/r);
  else if (i-1 == j)
    return 1.0 - (h/r);
  return 0.0;
}

@* Source Term.
Recall (\S17) the field equations resemble Poisson's equations 
$\nabla^{2}\phi = \sigma$. We therefore have
$$
\sigma(r) = -{{\partial V(\phi)}\over {\partial\phi}}
+{g_{\chi}\over\pi^{2}}m(r)^{3}i(p_{F}(r)/\abs{m(r)}).\eqn{}
$$
The discrete version of this is quite straightforward
$$
\sigma_{j} = -m_{\phi}^{2}\phi_{j} + {g_{\chi}\over\pi^{2}}m_{j}^{3}i(p_{F}(r_{j})/\abs{m_{j}})\eqn{}
$$
where $m_{j} = m_{\chi} + g_{\chi}\phi_{j}$ and $r_{j} = jh$.

@ @c
real NewtonSolver::source(index j) {
  real g = sqrt(4.0*alpha*PI);
  real scalingMass = massChi + g*phi[j];
  real p = pF(j);
  if (scalingMass == 0.0) throw std::logic_error("Divide by Zero");
  return -massPhi*massPhi*phi[j]+g*PI*pow(scalingMass/PI, 3.0)*i(p/fabs(scalingMass));
}

@ Observe the partial derivative of $\sigma$ with respect to the field
$\phi$ gives us
$$
{{\partial\sigma}\over{\partial\phi}} =
-{{\partial^{2}V}\over{\partial \phi^{2}}} + 
{g_{\chi}^{2}\over\pi^{2}}\left[3m(r)^{2}i(p_{F}(r)/\abs{m(r)})
-{p_{F}(r)^{2}\over\sqrt{m(r)^{2}+p_{F}(r)^{2}}}
\right]
$$
We can implement this precisely. Recall $\alpha=g_{\chi}^{2}/(4\pi)$.

@ @c real NewtonSolver::partialSource(index j) {
  real scalingMass = massChi + sqrt(4.0*PI*alpha)*phi[j];
  if (scalingMass == 0.0) throw std::logic_error("Divide by Zero");
  real fermiMomentum = pF(j);
  real p = fermiMomentum/fabs(scalingMass);
  real partialI = fermiMomentum*fermiMomentum*scalingMass/hypot(scalingMass,fermiMomentum);
  return -pow(massPhi, 2.0) + (alpha/PI_4)*(3.0*pow(scalingMass, 2.0)*i(p) - partialI);
}

@* Boundary Conditions.
The boundary conditions for the field are $\phi'(0)=0$ and, well, quite
complicated on the surface $r=R$. Lets focus first on the boundary
condition at the origin.

@ {\bf Origin Boundary Conditions.}
We see the boundary condition at the origin has the finite difference
approximation
$${1\over{2h}}(\phi_{1}-\phi_{-1})=0.\eqn{}$$
But what is $\phi_{-1}$? It's a ``ghost point'', which we can determine
via the field equation as
$$\phi_{1}-2\phi_{0}+\phi_{-1}=h^{2}\sigma_{0}\eqn{}$$
hence
$$\phi_{-1}=2\phi_{0}-\phi_{-1}+h^{2}\sigma_{0}.\eqn{}$$
This gives us the boundary condition
$${2\over h^{2}}\phi_{0}-{2\over h^{2}}\phi_{1}=\sigma_{0}.\eqn{}$$
We can implement this quite readily.

@ @c
real NewtonSolver::BoundaryConditionAtOrigin(index i) {
  if (i==0)      return 2.0;
  else if (i==1) return -2.0;
  else           return 0.0;
}

@ {\bf Surface Boundary Term.}
The boundary condition at the surface can likewise be solved. We
pretend outside the nugget, the scalar field acts like a Yukawa
potential. So that gives us
$$
\phi(r) = \phi(R){R\over r}e^{-m_{\phi}r}.\eqn{}
$$
Differentiation, then evaluating the result at $r=R$, gives us
$$\phi'(R) = -\phi(R)e^{-m_{\phi}R}\left({1\over R}+m_{\phi}\right).\eqn{}$$

@ {\bf Remark.}
This needs to be modified when we have a self-interacting term. Perhaps
one way out is to consider $\phi(\infty)=0$ as the boundary term? We
need to be sure it has the right asymptotic behaviour, but it's as good
as anything at the moment. (Luckily, we can put this problem on the back
burner, and focus on solving our free field situation.)

@ {\bf To Do.}
We should seriously consider instead the approximations
$$
\phi'(0)\approx {{-3\phi_{0} + 4\phi_{1}-\phi_{2}}\over{2h}}
+\bigO{h^{2}}
$$
and
$$
\phi'(R)\approx {{\phi_{N-2}-4\phi_{N-1}+3\phi_{N}}\over{2h}}+\bigO{h^{2}}
$$
These would produce suitably different results. We would have to, again,
modify the equations to make the resulting matrix tridiagonal\dots but
this would produce far more satisfactory results. The only cost is we
must assume $\phi(r)$ is sufficiently smooth. Observe the approximation
is
$$
-3\phi(x)+4\phi(x-h)-\phi(x-2h)=-2h\phi'(x)-{2h^{3}\over3}\phi'''(x)+\bigO{h^{4}}
$$
which is fairly decent, as good as our Laplacian approximation.

We would get for the boundary condition at the origin, after eliminating
the $\phi_{2}$ by adding $(h/2)$ times the next row of the matrix,
$$
\phi'(0)=0\mapsto 3{{-\phi_{0}+\phi_{1}}\over{2h}} = {-h\over2}\sigma_{1}
$$
and for the surface boundary condition
$$
\left[
{1\over{2h}}\left({N(N-1)\over{(N+1)(N-2)}}+3\right)
+\pmatrix{{\rm RHS\ of}\cr
{\rm Boundary}\cr
{\rm Condition}}
\right]\phi_{N} - {1\over h}\left[{N-1\over{N-2}}-2\right]\phi_{N-1}
= {-1\over 2}{{N-1}\over{N-2}}h\sigma_{N-1}.
$$

@
The finite difference approximation of this result is
$${{\phi_{N+1}-\phi_{N-1}}\over{2h}}=-\phi_{N}e^{-m_{\phi}R}\left({1\over R}+m_{\phi}\right).\eqn{}$$
Again, we must determine the ghost point $\phi_{N+1}$ in terms of prior
values through the field equation.

The field equations give us
$$\phi_{N+1} = {N\over {N+1}}h^{2}
\left[
  {2\over h^{2}}+{2\over R}e^{-m_{\phi}R}\left({1\over R}+m_{\phi}\right)
\right]
  \phi_{N}
-\left({{N-1}\over{N+1}}\right)\phi_{N-1}
+\left({N\over{N+1}}h^{2}\right)\sigma_{N}\eqn{}
$$
or simplifying a bit
$$\phi_{N+1} = {N\over {N+1}}h^{2}
\left[{2N\over {N+1}}
  + {{2N/R} \over {N+1}}e^{-m_{\phi}R}\left({1\over R}+m_{\phi}\right)
\right]
  \phi_{N}
-\left({{N-1}\over{N+1}}\right)\phi_{N-1}
+\left({N\over{N+1}}h^{2}\right)\sigma_{N}.\eqn{}
$$
Hence we get
$$
{1\over h}\left({N\over{N+1}}\right)\phi_{N-1}
-\left[
  {{N/h}\over{N+1}}
+\left(1 + {{N/R}\over{(N+1)h}}\right)\left({1\over R}+m_{\phi}\right)e^{-m_{\phi}R}
\right]\phi_{N}
=
{{hN/2}\over{N+1}}\sigma_{N}\eqn{}
$$
We want to have the right hand side be $\sigma_{N}$, so rearranging
terms gives us
$$
{2\over h^{2}}\phi_{N-1}
-\left[
  {2\over h^{2}}
+{{2Rh(N+1)+2N}\over{h^{2}R}}\left({1\over R}+m_{\phi}\right)e^{-m_{\phi}R}
\right]\phi_{N}
=
\sigma_{N}.\eqn{}
$$
Numerically, though, we delay dividing by $h^{2}$ until the ``last
possible moment''.

@ @c
real NewtonSolver::SurfaceBoundaryCondition(index j) {
  if (j<length-2) return 0.0;
  if (j==length-2) return 2.0;
  if (j==length-1) {
    return -2.0-((h*R*length-1+length)/(0.5*R))*(massPhi + 1.0/R)*exp(-massPhi*R);
  }
  else return 0.0;
}

@* Jacobian Matrix.
The Jacobian $J(\phi)_{ij} = M_{ij} - (\partial\sigma/\partial\phi)_{ij}(\phi)$ 
is diagonal dominant. We have three helper functions, which describe the
tridiagonal components of the ``symmetrized Jacobian'' $(J^{T}J)$:
$$ a_{i} = (J^{T}J)_{i,i-1}\eqn{}$$
$$ b_{i} = (J^{T}J)_{i,i} \eqn{}$$
$$ c_{i} = (J^{T}J)_{i,i+1} \eqn{}$$
We also have the modified right hand side of Newton's equations, since
we have to multiply by $J^{T}(\phi^{(n)})$.

@ @c
real NewtonSolver::JacobianMatrix(index i, index j) {
  if (i>=length || j>=length) throw std::out_of_range("Index out of bound");
  if (i==0) return BoundaryConditionAtOrigin(j);
  if (i==length-1) return SurfaceBoundaryCondition(j);
  
  if (j == i+1) {
    return Laplacian(i, j);
  }
  if (j == i-1) {
    return Laplacian(i,j);
  }
  if (j == i) return Laplacian(i, j)-0*h*h*partialSource(i);
  return 0.0; /* somehow someone asked for a zero entry */
}

@ We have some helper functions to access the off-diagonal and diagonal
components explicitly. While constructing this, we noted that
$J_{2,1}=0$ quite unexpectedly. This caused problems with
|CroutFactorization|, but enables us to write another algorithm that's
$\bigO{N}$ to invert the Jacobian.

@ Lets reflect about the system of equations. We have our system become
$$ J^{T}J\phi^{(n+1)}=J^{T}J\phi^{(n)}-J^{T}F(\phi^{(n)}).\eqn{}$$
The right hand side becomes
$$
{\rm RHS} = -M^{T}\left({\partial\sigma\over{\partial\phi}}\right)\phi^{(n)}
-M^{T}\sigma
+\left({\partial\sigma\over{\partial\phi}}\right)\sigma
+\left({\partial\sigma\over{\partial\phi}}\right)^{2}\phi^{(n)}\eqn{}
$$
and the left hand side
$$
{\rm LHS}_{j} = \left[
M^{T}M - \left({\partial\sigma\over{\partial\phi}}\right)^{T}M
- M^{T}\left({\partial\sigma\over{\partial\phi}}\right) +
\left({\partial\sigma\over{\partial\phi}}\right)^{T}\left({\partial\sigma\over{\partial\phi}}\right)
\right]\phi^{(n+1)}_{j}.\eqn{}
$$
We cheated a little since the source's Jacobian $\partial\sigma/\partial\phi$ 
is a diagonal matrix.

Observe further that, component-wise, $M^{T}M$ has:
$$
\pmatrix{
b_{0} & a_{1} &       & \dots\cr
c_{0} & b_{1} & a_{2} & \dots\cr
      & c_{1} & b_{2} & \dots\cr
      &       &       & \dots}
\pmatrix{
b_{0} & c_{0} &       & \dots\cr
a_{1} & b_{1} & c_{1} & \dots\cr
      & a_{2} & b_{2} & \dots\cr
      &       &       & \dots}
=
\pmatrix{
b_{0}^{2} + a_{1}^{2}   & b_{0}c_{0} + b_{1}a_{1}       &
 & \dots\cr
c_{0}b_{0}+b_{1}a_{1} & c_{0}^{2}+b_{1}^{2}+a_{2}^{2} & a_{2}b_{2}+b_{1}c_{1}
  & \dots\cr
                        &    c_{1}b_{1} + b_{2}a_{2}           &   &}
$$
Hence
$$
\tilde{a}_{j} = a_{j}b_{j}+c_{j-1}b_{j-1} = M_{j,j-1}M_{j,j} + M_{j-1,j}M_{j-1,j-1}\eqn{}
$$
$$
\tilde{b}_{j} = a_{j+1}^{2}+b_{j}^{2}+c_{j-1}^{2}=M_{j+1,j}^{2}+M_{j,j}^{2}+M_{j-1,j}^{2}\eqn{}
$$
and
$$
\tilde{c}_{j} = b_{j}c_{j}+a_{j+1}b_{j+1}=M_{j,j}M_{j,j+1}+M_{j+1,j}M_{j+1,j+1}\eqn{}
$$
describe the tridiagonal components of $(M^{T}M)$.

We also see
$$
\left({{\partial\sigma}\over{\partial\phi}}\right)^{T}M
=
\pmatrix{
\partial\sigma_{0} &                &                & \cr
               & \partial\sigma_{1} &                & \cr 
               &                & \partial\sigma_{2} & \cr
               &                &                &\ddots
}
\pmatrix{
b_{0} & c_{0} &       & \dots\cr
a_{1} & b_{1} & c_{1} & \dots\cr
      & a_{2} & b_{2} & \dots\cr
      &       &       & \ddots}
= \pmatrix{
b_{0}\partial\sigma_{0} & c_{0}\partial\sigma_{0} &       & \dots\cr
a_{1}\partial\sigma_{1} & b_{1}\partial\sigma_{1} & c_{1}\partial\sigma_{1} & \dots\cr
      & a_{2}\partial\sigma_{2} & b_{2}\partial\sigma_{2} & \dots\cr
      &       &       & \ddots}
$$
and
$$
M^{T}\left({{\partial\sigma}\over{\partial\phi}}\right)
=
\pmatrix{
b_{0} & a_{1} &       & \dots\cr
c_{0} & b_{1} & a_{2} & \dots\cr
      & c_{1} & b_{2} & \dots\cr
      &       &       & \ddots}
\pmatrix{
\partial\sigma_{0} &                &                & \cr
               & \partial\sigma_{1} &                & \cr 
               &                & \partial\sigma_{2} & \cr
               &                &                &\ddots
}
= \pmatrix{
b_{0}\partial\sigma_{0} & a_{1}\partial\sigma_{1} &       & \dots\cr
c_{0}\partial\sigma_{0} & b_{1}\partial\sigma_{1} & a_{2}\partial\sigma_{2} & \dots\cr
      & c_{1}\partial\sigma_{1} & b_{2}\partial\sigma_{2} & \dots\cr
      &       &       & \ddots}
$$

@ @c real NewtonSolver::a(index j) {
  if (j==0) return 0.0;
  real u = h*h;
  real MMtrans = (JacobianMatrix(j,j-1)*JacobianMatrix(j,j)+JacobianMatrix(j-1,j)*JacobianMatrix(j-1,j-1))/u;
  real partialSourceM = JacobianMatrix(j,j-1)*(partialSource(j)*u); /* $a_{j}\partial\sigma_{j}$ */
  partialSourceM += JacobianMatrix(j-1,j)*(partialSource(j-1)*u); /* $c_{j-1}\partial\sigma_{j-1}$ */
  return MMtrans-partialSourceM;
}
@ @c real NewtonSolver::b(index j) {
  if (j>length-1 || j<0) return 0.0;
  real MMtrans = pow(JacobianMatrix(j, j), 2.0);
  if (j<length-1) MMtrans += pow(JacobianMatrix(j+1,j), 2.0);
  if (j>0) MMtrans += pow(JacobianMatrix(j-1,j),2.0);
  real partialSourceM = (partialSource(j)*h*h*2)*JacobianMatrix(j,j);
  real sourceSquared = pow(partialSource(j)*h,2.0);
  return (MMtrans/(h*h))-partialSourceM+sourceSquared;
}
@ @c real NewtonSolver::c(index j) {
  if (j>=length-1) return 0.0;
  real u = h*h;
  real MMtrans = (JacobianMatrix(j,j+1)*JacobianMatrix(j,j)+JacobianMatrix(j+1,j+1)*JacobianMatrix(j+1,j))/u;
  real partialSourceM = JacobianMatrix(j,j+1)*(partialSource(j)*u); /* $c_{j}\partial\sigma_{j}$ */
  partialSourceM += JacobianMatrix(j+1,j)*(partialSource(j+1)*u); /* $a_{j+1}\partial\sigma_{j+1}$ */
  return MMtrans-partialSourceM;
}
@ @c real NewtonSolver::d(index j) {
  int ctr;
  real ret = 0.0;
  index start;
  if (j==0) start = 0;
  else {     start = j; start--; }
  for(ctr = start; ctr<length && ctr<j+1; ctr++) {
    ret += -JacobianMatrix(ctr,j)*(h*h*source(ctr)+h*h*partialSource(ctr)*phi[ctr]);
  }
  real p = partialSource(j);
  ret += (h*h*p)*(source(j)+p*phi[j]);
  return ret;
}

@ {\bf Unit Test.}
After all this, we should make certain that the matrix is symmetric with
some unit tests.

@ @c
bool LHSisSymmetricTests() {
  real alpha, massChi, R, N;
  NewtonSolver *f;
  alpha = 0.1; /* random numbers */
  massChi = 10.0;
  R = 0.4;
  N = 1e2;
  f = new NewtonSolver(10000, R, alpha, massChi, N);
  if (!f->isSymmetric()) return false;
  delete f;
  srand (time(NULL));
  for(int j=0; j<10; j++) {
    massChi = rand();
    alpha = rand();
    N = rand();
    R = rand();
    f = new NewtonSolver(10000, R, alpha, massChi, N);
    if (!f->isSymmetric()) return false;
    delete f;
  }
  return true;
}


@* Thomas' Algorithm.
When you have a tridiagonal matrix that's also diagonal dominant, you
can invert it in $\bigO{2N}$ operations instead of $\bigO{N^{3}}$
operations. Which is smashing. The basic idea is to do an LU
decomposition, then perform Gaussian elimination.

@ @c
void NewtonSolver::ThomasInvert() {
  iterateCounter++;
  index j;
  real m;
  std::swap(phi,nextPhi);

  @<Forward ...@>@;
  @<Backs...@>
}

@ @<Forward Sweep@>=
  cPrime[0] = c(0)/b(0);
  nextPhi[0] = phi[0]/b(0);
  for(j = 1; j < length-1; j++) {
    m = 1.0/(b(j) - a(j)*cPrime[j-1]);
    cPrime[j] = c(j)*m;
    nextPhi[j] = (d(j) - a(j)*nextPhi[j-1])*m;
  }

@ @<Backsubstitution@>=
  for(j = length - 2; j>=0; j--) {
    nextPhi[j] = nextPhi[j] - cPrime[j]*nextPhi[j+1];
    if (j==0) break;
  }

@* Crout Factorization.
Another LU-decomposition of a tridiagonal matrix which is not
necessarily diagonal dominant. Alas, the Jacobian matrix is rarely
diagonal dominant for our situation --- the three entries are roughly of
the same magnitudes.

If we write out the matrix multiplication $A=LU$ explicitly, we find:
$$
\left(\matrix{
a_{11} & a_{12} &      0 &    0   & \cdots & 0\cr
a_{21} & a_{22} & a_{23} &    0   & \ddots & \vdots\cr
   0   & a_{32} & a_{33} & a_{34} & \ddots & \vdots\cr
\vdots & \ddots & \ddots & \ddots & \ddots & \vdots\cr
0      & \cdots & \cdots & \cdots & a_{n,n-1} & a_{n,n}}
\right)
= \left(\matrix{
l_{11} & 0      & \dots & 0\cr
l_{21} & l_{22} & \ddots&\vdots\cr
0      & \ddots & \ddots&\vdots\cr
\vdots & \ddots & \ddots& 0\cr
0      & \dots  & l_{n,n-1} & l_{n,n}
}\right)
\left(\matrix{
1      & u_{12} & 0      & \dots & 0\cr
0      & 1      & u_{23} & \dots & 0\cr
\vdots & \ddots & \ddots & \ddots&\vdots\cr
\vdots & \ddots & \ddots & \ddots&\vdots\cr
0      & \dots  & \dots  & 0     &1}\right)
$$
Matrix multiplication will tell us:
$$a_{11} = l_{11}$$
$$a_{i,i-1} = l_{i,i-1}$$
for $i=2,3,\dots,n$;
$$a_{ii}=l_{i,i-1}u_{i-1,i} + l_{ii}$$
for $i=2,3,\dots,n$; and
$$
a_{i,i+1}=l_{i,i}u_{i,i+1}
$$
for $i=1,2,\dots,n-1$.


@ @c
void NewtonSolver::CroutFactorization() {
  iterateCounter++;
  index j;
  real l[2];
  real *u = cPrime;
  std::swap(phi,nextPhi);
  
  @<Solve Lower Triangular Part@>@;
  @<Solve Upper Triangular Part@>@;
}

@ @<Solve Lower...@>=
  std::cout<<"Solving the lower"<<std::endl;
  l[0] = 0.0;
  l[1]=b(0);
  if (l[1]==0.0)
    u[0]=0.0;
  else
    u[0] = c(0)/l[1];
  nextPhi[0] = d(0)/l[1];
  /* j-th row of L */
  for(j=1; j<length-1; j++) {
    l[0] = a(j);
    l[1] = b(j) - l[0]*u[j-1];
    if (l[1]==0.0) {
      std::cout<<"l["<<j<<"]["<<j<<"]=="<<l[0]<<std::endl;
    }
    u[j] = c(j)/l[1];
    nextPhi[j] = (d(j) - l[0]*nextPhi[j-1])/l[1];
  }
  /* last row */
  l[0] = a(length-1);
  l[1] = b(length-1) - l[0]*u[length-2];
  if (l[1]==0.0) {
    std::cout<<"l[N-1][N-1]=="<<l[0]<<std::endl;
  }
  nextPhi[length-1] = (d(length-1)-l[0]*nextPhi[length-2])/l[1];

@ @<Solve Upper...@>=
  std::cout<<"Solving the upper"<<std::endl;
  for(j=length-2; j>0; j--) {
    if (j==0) { std::cout<<"j==0"<<std::endl;
    std::cout<<"u[0] = "<<u[0]<<std::endl; }
    nextPhi[j]=nextPhi[j]-(u[j]*nextPhi[j+1]);
  }
  std::cout<<"We're done!!"<<std::endl;

@* Minimizing Energy.
We need to keep iteratively adjusting the radius $R$ until we minimize
the energy. How do we do this? Well, we should keep performing Newton's
method until the residual is ``good enough''. Then we perform {\it another}
Newton's method to adjust $R$ via the ``steepest descent''.

@ {\bf Residual.}
The residual for an approximate solution $\vec{x}^{*}$ to the system
$A\vec{x}=\vec{b}$ is precisely the vector
$\vec{r}=A\vec{x}^{*}-\vec{b}$. Its norm is a measure of error.

@ @c
real NewtonSolver::residual() {
  if (iterateCounter == 0) throw std::out_of_range("Must iterate before computing residual");
  real r = 0.0;
  real dr = 0.0;
  index j;
  for (j=1; j<length-1; j++) {
    dr=fabs((Laplacian(j,j-1)*phi[j-1]+Laplacian(j,j)*phi[j]+Laplacian(j,j+1)*phi[j+1])-h*h*source(j));
    r += dr*dr;
  }
  return sqrt(r);
}

@ @<Residual at Origin@>=
  r += pow((BoundaryConditionAtOrigin(0)*phi[0]+BoundaryConditionAtOrigin(1)*phi[1])-h*h*source(0),2.0);

@ @<Residual at Surface@>=
  r += pow((SurfaceBoundaryCondition(length-2)*phi[length-2]+SurfaceBoundaryCondition(length-1)*phi[length-1])-h*h*source(length-1),2.0);

@ @c
real NewtonSolver::norm() {
  if (iterateCounter == 0) throw std::out_of_range("Must iterate before computing residual");
  real r = 0.0;
  index j;
  for (j=0; j<length; j++) {
    r+=pow(nextPhi[j],2.0);
  }
  return sqrt(r);
}

@
We use a modified version of the method of steepest descent. Although
the fundamental theorem of calculus applies directly to the energy
integral, the integrand is always positive---hence the method will never
find a minima. We use a finite difference approximation of the
derivative of the energy with respect to $R$.

If the additional sliver of energy is negative, and if we want to
minimize the energy for the system, then we should increase $R$. (``A
little must be good, more must be better'' --- the motto of alcoholics
and gradient descent.)


@ @c
void run(NewtonSolver *f, real tol) {
  index i;
  real res_;
  f->ThomasInvert();
  for(i=0; 10>i; i++) {
    f->ThomasInvert();
    res_ = f->residual();
    std::cout<<"Residual: "<<res_<<std::endl;
    if (i>0 && res_<tol) break;
  }
}

@ @c
real nextR(index len, real R_, real alpha_, real massChi_, real N_, real massPhi_ = 0.0) {
  real E[2];
  real dE;
  NewtonSolver* f;
  real cbrtEps = cbrt(std::numeric_limits<real>::epsilon());
  real tol = cbrtEps*R_;
  f = new NewtonSolver(len+1, R_+cbrtEps, alpha_, massChi_, N_, massPhi_);
  run(f, tol);
  E[1] = f->energy();
  delete f;
  f = new NewtonSolver(len-1, R_-cbrtEps, alpha_, massChi_, N_, massPhi_);
  run(f, tol);
  E[0] = f->energy();
  dE = (E[1]-E[0])*0.5;
  std::cout<<"dE = "<<dE<<std::endl;
  while (fabs(dE)*5 >= R_) {
    dE *= 0.5;
  }
  delete f;
  std::cout<<">> dE* = "<<dE<<std::endl;
  return R_-dE;
}

@ @c
real norm(real *v, index length) {
  real r = 0.0;
  index j;
  for(j=0; j<length; j++) {
    r += v[j]*v[j];
  }
  return sqrt(r);
}

@ @c
bool NewtonSolver::isDiagonalDominant() {
  index j=0;

  for(j=3; j<length-1; j++) {
    if (fabs(b(j))<fabs(a(j))+fabs(c(j))) {
      std::cout<<"b("<<j<<") = "<<b(j)<<std::endl;
      std::cout<<"a("<<j<<")+c("<<j<<") = ";
      std::cout<<fabs(a(j))+fabs(c(j))<<std::endl;
      std::cout<<"Diff:"<<fabs(b(j))-fabs(a(j))-fabs(c(j))<<std::endl;
      std::cout<<"j = "<<j<<std::endl;
      return false;
    }
  }
  return true;
}
bool NewtonSolver::isSymmetric() {
  index j;
  for(j=0; j<length-1; j++) {
    if(fabs(c(j)-a(j+1))>1e-7) {
      std::cout<<"a("<<j+1<<") = "<<a(j+1)<<std::endl;
      std::cout<<"c("<<j<<") = "<<c(j)<<std::endl;
      return false;
    }
  }
  return true;
}

@ {\bf Initial Guess for Radius.}
Just a note to myself, the initial guess for the radius I picked
$\log_{10}(N)/m_{\chi}$ which appears to be as good as anything else.


@ @c void dumpData(NewtonSolver *f) {
  std::cout<<"\n>> R: "<<(f->currentR())<<std::endl;
  std::cout<<"Residual: "<<(f->residual())<<std::endl;
  std::cout<<"Norm:     "<<(f->norm())<<std::endl;
  std::cout<<"Momentum: "<<(f->pF(f->currentR(), 0.0))<<std::endl;
  std::cout<<"Energy: "<<(f->energy())<<std::endl;
}

@ @c
void minimizeR() { @#
  index j;
  index len;
  real N_ = 1e3;
  real alpha_ = 0.1;
  real massChi_ = 100.0;
  real R_ = 10*log10(N_)/massChi_;
  R_ = 0.62;
  real eps = std::numeric_limits<real>::epsilon();
  real cbrtEps = cbrt(eps);
  real tol;
  NewtonSolver *f;
  std::cout<<"Cube root of epsilon: "<<cbrtEps<<std::endl;
  std::cout<<"Epsilon: "<<eps<<std::endl;
  for(j=0; 80>j; j++) {
    tol = cbrtEps*R_;
    len = (index)(R_/cbrtEps); /* step size should be cuberoot of epsilon */
    if (len%2==1) len++; /* len should be even :) */
    f = new NewtonSolver(len, R_, alpha_, massChi_, N_);
    run(f, tol);
    dumpData(f);
    delete f;
    R_ = nextR(len, R_, alpha_, massChi_, N_);
    std::cout<<"Next R: "<<R_<<std::endl;
  }
}

@ @c void dumpEnergyPoints() {
  NewtonSolver *f;
  std::ofstream myfile;
  myfile.open("data.csv");
  real eps = std::numeric_limits<real>::epsilon();
  real cbrtEps = cbrt(std::numeric_limits<double>::epsilon());
  real sqrtEps = sqrt(eps);
  index len;
  real R_;
  real alpha_=0.1;
  real massChi_=100.0;
  real N_=1e3;
  real step = 0.02;
  index maxIter = (index)(2.0/step);
  for(index i=1; i<maxIter; i++) {
    /* sample 5 points around the location */
    for(int j=-5; j<5; j++) {
      R_ = i*step+j*sqrtEps;
      len = ((index)(R_/cbrtEps))+j*j;
      if (len%2==1) len++;
      std::cout<<"Number of data points: "<<len<<std::endl;
      f = new NewtonSolver(len, R_, alpha_, massChi_, N_);
      run(f, cbrtEps*R_);
      myfile<<std::setprecision(21)<<R_<<","<<(f->energy())<<std::endl;
      delete f;
    }
  }
}


@ @c void minimizeEnergyPoints() {
  NewtonSolver *f;
  real eps = std::numeric_limits<real>::epsilon();
  real cbrtEps = cbrt(std::numeric_limits<double>::epsilon());
  real sqrtEps = sqrt(eps);
  index len;
  real R_, r, dE;
  real rMin = -1.0;
  real rMax = 0.0;
  real eMin = 0.0;
  real eMax = 0.0;
  real alpha_=0.1;
  real massChi_=100.0;
  real N_=1e3;
  real step = 0.02;
  index maxIter = 20;
  R_ = 0.1;
  for(index i=0; i<maxIter; i++) {
    len = (index)(R_/cbrtEps);
    f = new NewtonSolver(len, R_, alpha_, massChi_, N_);
    run(f, cbrtEps*R_);
    std::cout<<"R: "<<R_<<std::endl;
    std::cout<<"energy: "<<(f->energy())<<std::endl;
    if (i==0 || (f->energy())<eMax) {
      eMax = f->energy();
      rMax = R_;
    } else {
      
      rMin = R_;
      eMin = f->energy();
    }
    if (rMin > 0.0) {
      R_ = 0.5*(rMin + rMax);
    } else {
      R_ += 1.0;
    }
    delete f;
  }
}

@ We use the 5-point stencil to approximate the change in energy. This
is a good $\bigO{h^5}$ approximation to the derivative.

@ @<Change in Energy@>=
  delete f;
  f = new NewtonSolver(len+2, R_+2*cbrtEps, alpha_, massChi_, N_);
  run(f, cbrtEps*R_);
  dE = -(f->energy());
  delete f;
  f = new NewtonSolver(len+2, R_+cbrtEps, alpha_, massChi_, N_);
  run(f, cbrtEps*R_);
  dE += 8*(f->energy());
  delete f;
  f = new NewtonSolver(len+2, R_-cbrtEps, alpha_, massChi_, N_);
  run(f, cbrtEps*R_);
  dE += -8*(f->energy());
  delete f;
  f = new NewtonSolver(len+2, R_-2*cbrtEps, alpha_, massChi_, N_);
  run(f, cbrtEps*R_);
  dE += (f->energy());
  delete f;
  dE = dE/(12*sqrtEps);
  std::cout<<">>> dE = "<<dE<<std::endl;

@ @c int main() {
  minimizeEnergyPoints();
  return 0;
}

@ {\bf Do not panic.}
The first residual printed indicates how badly the initial guess {\it was}.
It is usually large, since it compares the non-perturbative guess with
{\it the actual result}. Newton's method converges fairly quickly, each
iteration reducing the residual by a factor of about $10^{6}$. So
do not be shocked to find the sequence of residuals resemble: $10^5$,
$1$, $10^{-5}$, $10^{-11}$.

\vfill\eject