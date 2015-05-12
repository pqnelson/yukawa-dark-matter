\def\title{\uppercase{Yukawa Dark Matter}}

\input macros

@* Include Files.
Before expounding anything, we need to first note some parameters and
headers. This is boring, skip to the next page.

@d _USE_MATH_DEFINES true
@d PI M_PI
@d PI_4 M_PI_4
@d MAX(a,b) (a>b?a:b)
@d SQ(x) (x*x)
@d CUBE(x) (x*x*x)
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
#include <string>
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

@ @c real iDeriv(real z) {
  return z*z/hypot(1.0,z);
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

@ {\bf First-Order Equations of Motion.}
We claim that we may reduce the equations of motion to a first-order
differential equation:
$$
\phi'(r) = {1\over{r^{2}}}\int^{r}_{0}(r')^{2}\left[
{g_{\chi}\over \pi^{2}}m(r')^{3}i\left({{p_{F}(r')}\over\abs{m(r')}}\right)
-{{\partial V}\over{\partial\phi}}
\right]{\rm d}r'.
$$
We can see this from
$$
{{\rm d}\over{{\rm d}r}}[r^{2}\phi'(r)] = 
{g_{\chi}\over \pi^{2}}r^{2}m(r)^{3}i\left({{p_{F}(r)}\over\abs{m(r)}}\right)
-r^{2}{{\partial V}\over{\partial\phi}}
$$
which is precisely the equations of motion multipled by $r^{2}$. This
first-order ODE remains quite nonlinear.

Let $V(m)=V(m_{\chi}+g_{\chi}\phi)$. Then we may rewrite this as a
first-order differential equation in the scaling mass
$$
{{\rm d}\over{{\rm d}r}}m(r) = 
{g_{\chi}\over{r^{2}}}\int^{r}_{0}(r')^{2}\left[
{g_{\chi}\over \pi^{2}}m(r')^{3}i\left({{p_{F}(r')}\over\abs{m(r')}}\right)
-g_{\chi}{{\partial V}\over{\partial m(r)}}
\right]{\rm d}r'.
$$


@ {\bf Constraint: Positivity of Mass.}
We may implement a constraint condition, since
$$
m'(r) = g_{\chi}\phi'(r)
$$
we see
$$
m'(r)
=
{g_{\chi}\over{r^{2}}}\int^{r}_{0}(r')^{2}\left[
{g_{\chi}\over \pi^{2}}m(r')^{3}i\left({{p_{F}(r')}\over\abs{m(r')}}\right)
-{{\partial V}\over{\partial\phi}}
\right]{\rm d}r'.
$$
For the free situation $V=0$, we have $m'(r)>0$ and hence does not
change sign.


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
  real JacobianTransposeJacobian(index i, index j);
  void init();
  real a(index i);
  real b(index i);
  real c(index i);
  real d(index i);
  real aTilde(index i);
  real bTilde(index i);
  real cTilde(index i);
  real dTilde(index i);
  real F(index j);
  void CroutFactorization();
  void ThomasInvert(); /* Thomas' Algorithm */
  real residual(bool swapIter=false);
  real norm();
  bool isDiagonalDominant();
  bool isSymmetric(); 
  real relativeError();
  void setR(real R_) { R = R_; }
  real currentR() const { return R; }
  unsigned iterateCount() const { return iterateCounter; }
  real* solution() const { return nextPhi; }
  real get(index j) const { return phi[j]; }
  void swapPhi() { std::swap(phi, nextPhi); }
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
  real r;
  for (j=0; j<length; j++) {
    r = (j+1)*h;
    cPrime[j] = 0.0;
    nextPhi[j] = phi[j] = -alpha*1.5*(N/R)*(3-pow(r/R,2.0))*pow(massChi/pF(r),3.0)*i(pF(r)/massChi);
    if (j>0 && j<length-1) {
      if (partialSource(j)==0.0)
        std::cout<<"partialSource("<<j<<") is zero"<<std::endl;
    }
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
  if (i+1 == j) {
    if (i==0) return 3.0;
    return 1.0 + (h/r);
  }
  else if (i-1 == j) {
    if (i==0) return -3.0;
    return 1.0 - (h/r);
  }
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
  return -massPhi*massPhi*phi[j]+(g/(PI*PI))*pow(scalingMass, 3.0)*i(p/fabs(scalingMass));
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

@ {\bf Lemma.} {\it We see the derivative of $i'(p/m)$ term gives us}
$$
m(r)^{3}{{\partial i(p_{F}(r)/\abs{m(r)})}\over{\partial\phi(r)}}
= -g_{\chi}\left[{{p_{F}(r)}\over{\sqrt{m(r)^{2} + p_{F}(r)^{2}}}}
\right]p_{F}(r)\abs{m(r)}.
$$
\medbreak

\noindent{\it Proof.} We see
$$
i'(z) = {z^{2}\over{\sqrt{1+z^{2}}}}
$$
by definition of $i(z)$. Then we have
$$
{{\partial i(p_{F}(r)/\abs{m(r)})}\over{\partial\phi(r)}}
= i'(p_{F}(r)/\abs{m(r)}){{-p_{F}(r)}\over{m(r)^{2}}}{{\partial m(r)}\over{\partial\phi(r)}}.
$$
We claim
$$
{{\partial m(r)}\over{\partial\phi(r)}}=g_{\chi}.
$$
Then we get
$$
{{\partial i(p_{F}(r)/\abs{m(r)})}\over{\partial\phi(r)}}
= g_{\chi}\left[{{p_{F}(r)/\abs{m(r)}}\over{\sqrt{1 + \left({p_{F}(r)\over\abs{m(r)}}\right)^{2}}}}
\right]{{-p_{F}(r)}\over{m(r)^{2}}}.
$$
Rearranging terms gives
$$
{{\partial i(p_{F}(r)/\abs{m(r)})}\over{\partial\phi(r)}}
= g_{\chi}\left[{{p_{F}(r)}\over{\sqrt{m(r)^{2} + p_{F}(r)^{2}}}}
\right]{{-p_{F}(r)}\over{m(r)^{2}}}.
$$
Hence we find the desired result.\quad\slug

@ @c real NewtonSolver::partialSource(index j) {
  if (j==0) throw std::out_of_range("Attempting to get partialSource(0)");
  if (j>=length-1) throw std::out_of_range("Attempting to get partialSource(j>N-1)");
  real scalingMass = massChi + sqrt(4.0*PI*alpha)*phi[j];
  if (scalingMass == 0.0) throw std::logic_error("Divide by Zero");
  real fermiMomentum = pF(j);
  real p = fermiMomentum/fabs(scalingMass);
  
  real partialI = pow(fermiMomentum,2.0)*fabs(scalingMass)/hypot(scalingMass,fermiMomentum);
  if (isnan(phi[j])) {
    std::cout<<"phi["<<j<<"] = "<<phi[j]<<std::endl;
    throw std::invalid_argument("phi[j] is a NaN");
  }
  if (isnan(massChi)) throw std::invalid_argument("massChi is a NaN");
  if (isnan(scalingMass)) throw std::invalid_argument("scalingMass is a NaN");
  if (isnan(p)) throw std::invalid_argument("p is a NaN");
  if (isnan(partialI)) throw std::invalid_argument("partialI is a NaN");
  return -pow(massPhi, 2.0) + (alpha/PI_4)*(3.0*pow(scalingMass, 2.0)*i(p) - partialI);
}

@* Boundary Conditions.
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
must assume $\phi(r)$ is sufficiently smooth. 

@ {\bf Proposition.} {\it We have the following approximation hold}
$$
-3\phi(x)+4\phi(x-h)-\phi(x-2h)=-2h\phi'(x)-{2h^{3}\over3}\phi'''(x)+\bigO{h^{4}}
$$

\medbreak
\noindent{\it Proof.} We simply use the Taylor series expansion
$$
\phi(x-2h)=\phi(x)-2h\cdot\phi'(x)+{4h^{2}\over{2!}}\phi''(x)-{8h^{3}\over{3!}}\phi'''(x)+\bigO{h^{4}}
$$
$$
\phi(x-h)=\phi(x)-h\cdot\phi'(x)+{h^{2}\over{2!}}\phi''(x)-{h^{3}\over{3!}}\phi'''(x)+\bigO{h^{4}}
$$
hence the relation holds immediately.\quad\slug

@ {\bf Origin Boundary Conditions.}
We may find the boundary conditions at the origin as
$$
\phi'(0)=0\mapsto{{3\phi_{0}-4\phi_{1}+\phi_{2}}\over{2h}}=0.
$$
However, for numerical purposes, we strongly desire this equation to
involve only $\phi_{0}$ and $\phi_{1}$. (This way, we'd get a
tridiagonal matrix, and we could invoke a host of nice numerical methods
that are optimized for tridiagonal matrices.)

We would get for the boundary condition at the origin, after eliminating
the $\phi_{2}$ by adding $(h/2)$ times the next row of the matrix,
$$
\phi'(0)=0\mapsto 3{{-\phi_{0}+\phi_{1}}\over{2h}} = {-h\over2}\sigma_{1}.\eqn{}
$$
This presents the solution to the unstable boundary condition we
couldn't tackle with ghost points. For simplicity, we multiply through
by $-h/2$:
$$
{-6\over {h^{2}}}\phi_{0}+{6\over{h^{2}}}\phi_{1}=\sigma_{1}\eqn{}
$$
which nicely speeds things up. Again, delaying dividing by $h^{2}$ until
the last possible moment gives us:

@ @c
real NewtonSolver::BoundaryConditionAtOrigin(index i) {
  if (i==0) return -6.0;
  if (i==1) return 6.0;
  return 0.0;
}

@ {\bf Surface Boundary Condition.}
This is a bit more involved, so let us be careful. We will denote
$$
\phi'(R) = c_{0}\phi(R)
$$
as the boundary condition on the surface. Then we have
$$
{1\over{2h}}\phi_{N-2}-{4\over{2h}}\phi_{N-1}+{3\over{2h}}\phi_{N} = c_{0}\phi_{N}
$$
approximate the condition. To eliminate $\phi_{N-2}$, we examine the
equations of motion from the previous row of the system
$$
{N\over{N-1}}{1\over{h^{2}}}\phi_{N}
-{2\over{h^{2}}}\phi_{N-1}
+{{N-2}\over{N-1}}{1\over{h^{2}}}\phi_{N-2}
=\sigma_{N-1}.
$$
We manipulate this equation to have the coefficient of $\phi_{N-2}$ be unity
$$
{{N}\over{N-2}}{1\over{h^{2}}}\phi_{N}-{{N-1}\over{N-2}}{2\over{h^{2}}}\phi_{N-1}
+{1\over{h^{2}}}\phi_{N-2}
={{N-1}\over{N-2}}\sigma_{N-1}
$$
We multiply through by $h/2$:
$$
{{N}\over{N-2}}{1\over{2h}}\phi_{N}
-{{N-1}\over{N-2}}{1\over{h}}\phi_{N-1}
+{1\over{2h}}\phi_{N-2}
={{N-1}\over{N-2}}{h\over 2}\sigma_{N-1}.
$$
Now we take our surface boundary condition and subtract out this
manipulated equation of motion:
$$
\left(3-{{N}\over{N-2}}\right){1\over{2h}}\phi_{N}
+\left(-2+{{N-1}\over{N-2}}\right){1\over{h}}\phi_{N-1}
=c_{0}\phi_{N}-{{N-1}\over{N-2}}{h\over 2}\sigma_{N-1}.
$$
Great. Now to make this uniform with the rest of the system, we should
have the right hand side be just $\sigma_{N-1}$.

Our first step is to throw $c_{0}\phi_{N}$ to the other side:
$$
\left(3-{{N}\over{N-2}}-2h c_{0}\right){1\over{2h}}\phi_{N}
+\left(-2+{{N-1}\over{N-2}}\right){1\over{h}}\phi_{N-1}
=-{{N-1}\over{N-2}}{h\over 2}\sigma_{N-1}.
$$
Now we just multiply through by the inverse of the coefficient of
$\sigma_{N-1}$ to get
$$
-{{N-2}\over{N-1}}{2\over h}\left(3-{{N}\over{N-2}}-2h c_{0}\right){1\over{2h}}\phi_{N}
-{{N-2}\over{N-1}}{2\over h}\left(-2+{{N-1}\over{N-2}}\right){1\over{h}}\phi_{N-1}
=\sigma_{N-1}.
$$
Gathering terms first gives us
$$
-{{N-2}\over{N-1}}{2\over h}\left({{2N-6}\over{N-2}}-2h c_{0}\right){1\over{2h}}\phi_{N}
-{{N-2}\over{N-1}}{2\over h}\left({{4-N}\over{N-2}}\right){1\over{h}}\phi_{N-1}
=\sigma_{N-1}.
$$
Hence we find
$$
-\left({{2N-6}\over{N-1}}-2h{{N-2}\over{N-1}}c_{0}\right){1\over{h^{2}}}\phi_{N}
+\left({{N-4}\over{N-1}}\right){2\over{h^{2}}}\phi_{N-1}
=\sigma_{N-1}.
$$
We simply take
$$
c_{0} = -Re^{-m_{\phi}R}\left[m_{\phi}+{1\over{R^{2}}}\right]
$$
and we're done.

@ @c
real NewtonSolver::SurfaceBoundaryCondition(index j) {
  index N_ = length-1;
  if (j==N_-1) {
    return (((N_-4)*2.0)/((N_-1)*1.0));
  }
  else if (j==N_) {
    real c = -exp(-massPhi*R)*(R*massPhi + (1.0/R));
    real ret = -2.0*((1.0*(N_-3))/(1.0*(N_-1)));
    ret += (2.0*h*(N_-2)/(1.0*(N_-1)))*c;
    return ret;
  }
  return 0.0;
}

@ {\bf On Ghost Boundary Conditions.}
There is one situation of note I should mention. The naive finite
difference one might attempt to use (the central divided difference)
would involve points like $\phi_{-1}$ or $\phi_{N+1}$, which we don't
keep track of. One way out is to use the equations of motion to
eliminate these points. Sometimes this works. But not for the
spherically symmetric Laplacian---at least, not for the $\phi'(0)$
condition.

Why? Well, the spherically symmetric Laplacian has a
$r^{-1}\partial_{r}$ term. We would need
$$
\lim_{\varepsilon\to0}{{\phi(\varepsilon)-\phi(0)}\over{\varepsilon^{2}}}=0
$$
or more precisely
$\phi(\varepsilon)\approx\phi(0)+\bigO{\varepsilon^{3}}$. We cannot make
this assumption, it'd make the entire system ``gravitate'' towards the
$\phi(r)=0$ solution.

@* Jacobian Matrix.
The Jacobian $J(\phi)_{ij} = M_{ij} - (\partial\sigma/\partial\phi)_{ij}(\phi)$ 
is diagonal dominant. Mostly. We see for $2\leq i\leq N$, the matrix is
diagonal dominant, for $i=1$ it's weakly diagonal dominant, and for the
boundary condition at the origin\dots again weakly diagonal dominant.

It'd be nicer if the Jacobian were {\it strictly} diagonal dominant, or
symmetric and positive-definite. The latter condition may be established
if we consider instead $J^{T}J$\dots which is harder. For now, we
consider the vanilla Jacobian matrix.

@ @c
real NewtonSolver::JacobianMatrix(index i, index j) {
  if (i>=length || j>=length) throw std::out_of_range("Index out of bound");
  if (i==0) {
    real ret = BoundaryConditionAtOrigin(j)/h;
    if (j==1) ret+=h*partialSource(1);
    return ret;
  }
  if (i==length-1) {
    real ret = SurfaceBoundaryCondition(j)/h;
    if (j==length-2) ret+=h*partialSource(j);
    return ret;
  }

  if (j == i+1) {
    return Laplacian(i, j)/h;
  }
  if (j == i-1) {
    return Laplacian(i,j)/h;
  }
  if (j == i) {
    if (i==0) {
      std::cout<<"Somehow got to (i,j)=(0,0) illegally..."<<std::endl;
    }
    return Laplacian(i, j)/h-h*partialSource(i);
  }
  return 0.0; /* somehow someone asked for a zero entry */
}

@ We also include the matrix multiplied $J^{T}J$. This is quick, since
the Jacobian is tridiagonal. Observe further that, component-wise, $J^{T}J$ has:
$$
\pmatrix{
b_{0} & a_{1} &       & \dots\cr
c_{0} & b_{1} & a_{2} & \dots\cr
      & c_{1} & b_{2} & \dots\cr
      &       &       & \ddots}
\pmatrix{
b_{0} & c_{0} &       & \dots\cr
a_{1} & b_{1} & c_{1} & \dots\cr
      & a_{2} & b_{2} & \dots\cr
      &       &       & \ddots}
=
\pmatrix{
b_{0}^{2} + a_{1}^{2}   & b_{0}c_{0} + b_{1}a_{1}       &
 & \dots\cr
c_{0}b_{0}+b_{1}a_{1} & c_{0}^{2}+b_{1}^{2}+a_{2}^{2} & a_{2}b_{2}+b_{1}c_{1}
  & \dots\cr
                        &    c_{1}b_{1} + b_{2}a_{2}           &   &}
$$
Hence we conclude the components of the resulting matrix are
Hence
$$
\tilde{a}_{j} = a_{j}b_{j}+c_{j-1}b_{j-1} = J_{j,j-1}J_{j,j} + J_{j-1,j}J_{j-1,j-1}\eqn{}
$$
$$
\tilde{b}_{j} = a_{j+1}^{2}+b_{j}^{2}+c_{j-1}^{2}=J_{j+1,j}^{2}+J_{j,j}^{2}+J_{j-1,j}^{2}\eqn{}
$$
and
$$
\tilde{c}_{j} = b_{j}c_{j}+a_{j+1}b_{j+1}=J_{j,j}J_{j,j+1}+J_{j+1,j}J_{j+1,j+1}\eqn{}
$$

@ @c real NewtonSolver::JacobianTransposeJacobian(index i, index j) {
  real ret = 0.0;
  if (j==i) { /* $\tilde{b}_{i}$ */
    ret += pow(JacobianMatrix(j,j),2.0);
    if (j>0) ret += pow(JacobianMatrix(j-1,j-1), 2.0);
    if (j<length-1) ret += pow(JacobianMatrix(j+1,j+1),2.0);
  }
  else if (j==i-1) { /* $\tilde{a}_{i}$ */
    ret += JacobianMatrix(i,i)*JacobianMatrix(i,j); /* $b_{i}a_{i}$ */
    if (j<length-1) ret += JacobianMatrix(j,i)*JacobianMatrix(j,j);
  }
  else if (j==i+1) { /* $\tilde{c}_{i}$ */
    ret += JacobianMatrix(j,j)*JacobianMatrix(j,i); /* $b_{i+1}a_{i+1}$ */
    ret += JacobianMatrix(i,j)*JacobianMatrix(i,i);
  }
  return ret;
}

@ @c real NewtonSolver::aTilde(index j) {
  if (j==0) return 0.0;
  return JacobianTransposeJacobian(j,j-1);
}
@ @c real NewtonSolver::bTilde(index j) {
  if (j>length-1 || j<0) return 0.0;
  return JacobianTransposeJacobian(j,j);
}
@ @c real NewtonSolver::cTilde(index j) {
  if (j>=length-1) return 0.0;
  return JacobianTransposeJacobian(j,j+1);
}
@ @c real NewtonSolver::a(index j) {
  return aTilde(j);
  if (j>0) return JacobianMatrix(j, j-1);
  return 0.0;
}

@ @c real NewtonSolver::b(index j) {
  return bTilde(j);
  if (j<0 || j>length-1) return 0.0;
  return JacobianMatrix(j,j);
}

@ @c real NewtonSolver::c(index j) {
  return cTilde(j);
  if (j>=length-1) return 0.0;
  return JacobianMatrix(j,j+1);
}

@ The system we are trying to solve is
$$
\vec{F}(\phi) = (M\phi - \sigma(\phi))
$$
It will come in handy to implement this as a utility function.

@ @c
real NewtonSolver::F(index j) {
  real ret =0.0;
  if (j==0) {
    ret = (BoundaryConditionAtOrigin(0)*phi[0]);
    ret += (BoundaryConditionAtOrigin(1)*phi[1]);
    ret = (ret/h) - source(1)*h;
  } else if (j==length-1) {
    ret = (SurfaceBoundaryCondition(length-2)*phi[length-2]);
    ret += (SurfaceBoundaryCondition(length-1)*phi[length-1]);
    ret = (ret/h) - source(length-2)*h;
  } else if (j<length-1) {
    ret = Laplacian(j,j)*phi[j];
    if (j<length-1) ret += Laplacian(j,j+1)*phi[j+1];
    if (j>0) ret += Laplacian(j,j-1)*phi[j-1];
    ret = (ret/h) - source(j)*h;
  }
  return ret;
}

@ @c real NewtonSolver::d(index j) {
  return dTilde(j);
  real ret = JacobianMatrix(j,j)*phi[j]-F(j);
  if(j>0) ret+= JacobianMatrix(j,j-1)*phi[j-1];
  if(j<length-1) ret+=JacobianMatrix(j,j+1)*phi[j+1];
  return ret;
}

@ @c real NewtonSolver::dTilde(index j) {
  real ret = JacobianTransposeJacobian(j,j)*phi[j]-JacobianMatrix(j,j)*F(j);
  if (j>0)
    ret += JacobianTransposeJacobian(j,j-1)*phi[j-1]-JacobianMatrix(j-1,j)*F(j-1);
  if (j<length-1)
    ret += JacobianTransposeJacobian(j, j+1)*phi[j+1]-JacobianMatrix(j+1,j)*F(j+1);
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

@* Minimizing Energy.
We need to keep iteratively adjusting the radius $R$ until we minimize
the energy. How do we do this? Well, we should keep performing Newton's
method until the residual is ``good enough''. Then we perform {\it another}
Newton's method to adjust $R$ via the ``steepest descent''.

@ {\bf Residual.}
The residual for an approximate solution $\vec{x}^{*}$ to the system
$A\vec{x}=\vec{b}$ is precisely the vector
$\vec{r}=A\vec{x}^{*}-\vec{b}$. Its norm is a measure of error. We see
that this is precisely the norm of |NewtonSolver::F|.

@ @c
real NewtonSolver::residual(bool swapIter) {
  if (iterateCounter == 0) throw std::out_of_range("Must iterate before computing residual");
  if(swapIter) std::swap(phi,nextPhi);
  real r = 0.0;
  real dr = 0.0;
  index j;
  for (j=1; j<length-1; j++) {
    dr=fabs(F(j));
    r += dr*dr;
  }
  if(swapIter) std::swap(phi,nextPhi);
  return sqrt(r);
}

@ {\bf Relative Error.}
On the other hand, we can measure a notion of ``the vectors are
converging towards something''. This is $\norm{\phi^{(n+1)}-\phi^{(n)}}$,
or the ``relative error''.

@ @c real NewtonSolver::relativeError() {
  real ret = 0.0;
  for(index j=0; j<length; j++) {
    ret += pow(phi[j]-nextPhi[j],2.0);
  }
  return sqrt(ret);
}

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
  real res_, prevRes, relErr;
  f->ThomasInvert();
  res_ = f->residual(true);
  std::cout<<"Residual(0): "<<res_<<std::endl;
  if (res_ < 1.0) return;
  for(i=0; 100>i; i++) {
    prevRes = res_;
    f->ThomasInvert();
    relErr = f->relativeError();
    res_ = f->residual();
    if (i%10==0 || (res_<tol || relErr<tol)) {
      std::cout<<"Residual("<<(i+1)<<"): "<<res_<<std::endl;
      std::cout<<"Relative Error("<<i+1<<"): "<<relErr<<std::endl;
    }
    if (i>0 && (res_<tol || relErr<tol)) break;
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
  index maxIter = 100;
  R_ = 0.6;
  for(index i=0; i<40; i++) {
    len = (index)(R_/cbrtEps);
    f = new NewtonSolver(len, R_, alpha_, massChi_, N_);
    run(f, cbrtEps*R_);
    if (i==0 || (f->energy())<eMax) {
      eMax = f->energy();
      rMax = R_;
    }
    R_ += 0.002;
    std::cout<<"E("<<R_<<") = "<<(f->energy())<<std::endl;
    delete f;
  }
  std::cout<<"R minimized at: "<<rMax<<std::endl;
  std::cout<<"Energy: "<<eMax<<std::endl;
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

@ {\bf Accuracy, Convergence.}
Although the relative error starts small, and decreases quickly, the
residual remains unacceptably large. (Something like $\sim\bigO{10\cdot N}$ for
$N$ components in the vector $\phi_{j}$.) Consequently, this approach
does not appear to work: the system is just too nonlinear. Newton's
method converges, but to the wrong value.

\vfill\eject