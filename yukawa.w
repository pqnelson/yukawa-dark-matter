\def\title{\uppercase{Yukawa Dark Matter}}

\input macros

@s real double
@s index size_t
@s constexpr const
@s NoSolutionExistsException char
@s MassOutOfBoundsException char
@s string char
@s std int
@s bool int
@s complex int
@s ios int
@s iterator int
@s less int
@s map int
@s pair int
@s static_cast int
@s vector int
@f numeric_limits make_pair

\def\topofcontents{
 \null\vskip-35pt
 \vfill
 \vfill
 \centerline{\bf Numerical Simulations of Yukawa-Coupled Dark Matter}
 \bigskip
 \centerline{Alex Nelson\footnote{${}^{\dagger}$}{\eightrm The latest version of
 this program may be found at \homepage{}.}}
 \smallskip
 \vskip 4\bigskipamount
 \centerline{\vtop{\hsize=.8\hsize
    \noindent{\bf Abstract.}\quad
    We numerically investigate and reproduce recent results from Wise and Zhang
    (\arXiv{1411.1772}). The basic model couples fermionic matter to a
    massless scalar field. We extend the model to a massive
    scalar field. Furthermore, we also discuss the
    numerical analysis of the problem in considerable detail. A working
    program is included.}}
 \vskip.7in
 \vfill} % this material will start the table of contents page

% the next lines patch cwebmac so that the contents page comes out first
\let\oldcon=\con
\let\con=\bye

\newread\tocfile \openin\tocfile=\contentsfile\relax
\ifeof\tocfile
 \message{Please run this file again thru TeX to get the table of contents!}
\else \pageno=0 \titletrue {\let\end=\relax \oldcon} \fi

\def\lheader{\mainfont\the\pageno\eightrm\qquad\grouptitle\hfill
  ALEX NELSON\qquad
  \mainfont\topsecno} % top line on left-hand pages
\def\rheader{\mainfont\topsecno\eightrm\qquad DARK MATTER\hfill
  \eightrm\grouptitle\qquad\mainfont\the\pageno} % top line on right-hand pages

@* Introduction, Review.
Following Mark Wise and Yue Zhang's work in \arXiv{1411.1772}, we
@^Wise, Mark@>
@^Zhang, Yue@>
consider a fermionic field $\chi$ coupled to a scalar field
$\phi$. Initially we consider $\phi$ massless and
non-self-interacting. For simplicity we assume spherical symmetry. 


@c
@<Include ...@>@;

class YukawaDarkMatter {
private:
        @<Model Parameters@>@;
public:
        @<Model constructor@>@;
        @<Mutator Methods@>@;
        @<Accessor Methods@>@;
        @<Routines@>@;
};

@ {\bf Lagrangian.} The system's Lagrangian density is
$$
{\cal L} = i\bar{\chi}\gamma^{\mu}\partial_{\mu}\chi
 - m_{\chi}\bar{\chi}\chi - g_{\chi}\bar{\chi}\chi\phi
 + {1\over 2}\partial\phi\cdot\partial\phi.\eqn{}
$$
We assume $g_{\chi}>0$.
@^notation: $g_{\chi}$@>
This is hard, so instead of considering a fermionic field $\chi$, we
just consider a finite number $N$ of relativistic particles. 
@^Fermionic Field@>

The Lagrangian we work with instead is (setting $\partial_{t}\phi=0$ and
working in mostly-minuses convention)
@^Steady-State Approximation@>
$$
L = -\sum_{i}m(\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}-{1\over2}\int\nabla\phi\nabla\phi\,{\rm d}^{3}x.\eqn{}
$$\labelEqn\defnLagrangian
We have introduced the so-called ``scale-invariant mass''
@^notation: $m(r)$@>
@^Scale-Invariant Mass@>
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
\nabla\cdot{\partial L\over {\partial\nabla\phi(\vec{x})}}=-\nabla^{2}\phi(\vec{x})
$$
and
$$
{\partial L\over{\partial\phi(\vec{x})}}=-\sum_{i}g_{\chi}\delta^{(3)}(\vec{x}-\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}.
$$
Thus we deduce the equations of motion for $\phi$:
$$
\nabla^{2}\phi(\vec{x})
=\sum_{i}g_{\chi}\delta^{(3)}(\vec{x}-\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}.\eqn{}
$$\labelEqn\scalarFieldEqns

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
$$\labelEqn\hamiltonionDefn

@ {\bf Nonzero Scalar Potential Term.}
If we had some nonzero potential $V(\phi)$ for the scalar field, how
would we change things? Well, we would have to modify Eq (\defnLagrangian) as
$$
L = -\sum_{i}m(\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}+\int{-1\over2}\nabla\phi\nabla\phi-V(\phi)\,{\rm d}^{3}x
\eqno{(\defnLagrangian{}')}
$$
This would produce an extra term in the $\partial L/\partial\phi$
equation, modifying Eq (\scalarFieldEqns) as
$$
\nabla^{2}\phi(\vec{x})
=\sum_{i}g_{\chi}\delta^{(3)}(\vec{x}-\vec{x}_{i})\sqrt{1-\dot{\vec{x}}_{i}^{2}}
+{\partial V(\phi)\over{\partial\phi(\vec{x})}}.
\eqno{(\scalarFieldEqns{}')}
$$
This does not alter the definition of $\vec{p}_{i}$, but the Hamiltonian
changes in several ways: first, since we inserted the equations of
motion into (\hamiltonionDefn), and second since the Lagrangian gains an
extra term for the potential $V(\phi)$. We adjust our Hamiltonian to become
$$
H = \sum_{i}\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}} 
- {g_{\chi}\over 2}\sum_{i}\phi(\vec{x}_{i})
  {m(\vec{x}_{i})\over{\sqrt{\vec{p}_{i}^{2}+m(\vec{x}_{i})^{2}}}}
+\int\left(V(\phi)-{1\over 2}\phi(\vec{x}){\partial V(\phi)\over{\partial\phi(\vec{x})}}\right)
{\rm d}^{3}x.
\eqno{(\hamiltonionDefn{}')}
$$

@c 
real YukawaDarkMatter::scalarPotential(real phi) {
  if (m_scalarMass==0.0) return 0.0;
  return 0.5*SQ(phi*m_scalarMass);
}

real YukawaDarkMatter::partialScalarPotential(real phi) {
  if (m_scalarMass==0.0) return 0.0;
  return phi*SQ(m_scalarMass);
}

@ {\bf Continuum Limit.}
When $N\gg 1$, we may treat the collection of particles as a continuous
@^Continuum Limit@>
``field'', at least trading the sum over particles for integrals over
the phase space:
$$
\sum_{i}\to\int{\rm d}^{3}r\int{{{\rm d}^{3}p}\over{(2\pi)^{3}}}
f(\vec{r}, \vec{p}).\eqn{}
$$
As always, we assume spherical symmetry. We have the Fermi
@^Fermi Momentum@>
@^notation: $p_{F}(r)$@>
momentum\footnote{${}^{1}$}{I honestly have not found a good reference
on what the ``Fermi momentum'' {\it is}. As best as I can tell, it's the
momenta for the ``most energized'' particle state in some configuration
of particles.} $p_{F}(r)$ which may depend on position, and the $\chi$
particles must be confined to a sphere of radius $R$. Hence we conclude
$$
f(\vec{r}, \vec{p}) = 2\theta(R-r)\theta(p_{F}(r)-p)\eqn{}
$$
where $\theta(-)$ is the Heaviside step function.
@^Heaviside Step Function@>
@^notation: $\theta(-)$@>

@ {\bf Number of $\chi$ Particles.}
We may consider the number of $\chi$ particles as
$$
N = \sum_{i}1 = {4\pi\over 3}\int^{R}_{0} p_{F}(r)^{3}r^{2}\,{\rm d}r.\eqn{}
$$
This parameter will be fixed, and given by the user.
@^notation: $N$@>

@ {\bf Energy's Continuum Limit.}
So, we consider the Hamiltonian (\hamiltonionDefn) and how it changes
under the continuum limit. We see
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
  V(\phi)-{1\over 2}\phi(r){\partial V(\phi)\over{\partial\phi(r)}}
\right]{\rm d}r\eqn{}
$$
to the total energy function.\callthis\energyContinuumLimit

\proclaim Theorem. For large $z$, 
$$i(z)\approx {z^{2}\over 2}+ {1\over 4}(1-\ln(4)-2\ln(z))+\bigO{z^{-2}}$$
and 
$$h(z)\approx {{1-2\ln(4)-4\ln(z)}\over{32}}+{z^{2}\over 4}+{z^{4}\over 4}+\bigO{z^{-2}}.$$

\rmk
We will use this theorem when $0<m(r)<1/10$, and hence
$\infty>1/m(r)>10$. How good of an approximation is this? Well, we find
$$
\int^{1/10}_{0}\norm{h(w^{-1})-h_{{\rm approx}}(w^{-1})}^{2}\,{\rm
d}w\approx 1.94\cdot10^{-9}
$$
and likewise
$$
\int^{1/10}_{0}\norm{i(w^{-1})-i_{{\rm approx}}(w^{-1})}^{2}\,{\rm
d}w\approx 6.989653\cdot10^{-8}.
$$


@c namespace Util {
  real i(real z) {
       return 0.5*(z*sqrt(1.0+z*z) - asinh(z));
  }
  real iDeriv(real z) {
        return z*z/hypot(1.0,z);
  }
  real h(real z) {
    real u = z*z;
    return 0.25*(i(z) + z*u*hypot(1.0,z));
  }
}
const real LN_4 = 1.3862943611198906188344642429164L;
real YukawaDarkMatter::energyDensity(real mass, real r) {
  real phi = massToField(mass);
  real massTerm, freeTerm;
  massTerm = scalarPotential(phi)-0.5*phi*partialScalarPotential(phi);
  if (fabs(mass)<0.1) {
    real p = fermiMomentum(r);
    real logP = log(p);
    real logM = log(mass);
    real hTerm = 0.25*pow(p,4.0)+SQ(mass)*(0.25*SQ(p)+0.03125*SQ(mass)*(1.0-2.0*LN_4-4.0*logP+4.0*logM));
    real iTerm = mass*(0.5*SQ(p)+0.25*(1.0-LN_4-2.0*logP+2.0*logM));
    freeTerm = (hTerm - 0.5*coupling()*iTerm*phi);
    return 4.0*PI_3*SQ(r)*(freeTerm + massTerm);
  } else {
    real p = fermiMomentum(r)/fabs(mass);
    freeTerm = CUBE(mass)*(mass*Util::h(p) - 0.5*coupling()*phi*Util::i(p));
    return 4.0*PI_3*SQ(r)*(freeTerm + massTerm);
  }
}

@ {\bf Continuum Limit of Equations of Motion.}
Insider the nugget, when $r<R$, the equations of motion for $\phi$
changes to
$$
\nabla^{2}\phi(r) = {{{\rm d}^{2}}\over {{\rm d}r^{2}}}\phi(r)+{2\over
r}{{{\rm d}}\over {{\rm d}r}}\phi(r) =
{g_{\chi}\over\pi^{2}}m(r)^{3}i(p_{F}(r)/|m(r)|)+{\partial V(\phi)\over{\partial\phi(r)}}\eqn{}
$$
This is great, it's the equation we will be solving, but we have one
glaring problem: we never specified $p_{F}(r)$.
@^Continuum Limit@>

\rmk
Again, we use the approximation when $m(r)<1/10$ that
$i(p_{F}(r)/\abs{m(r)})$ is approximately cubic plus some logarithmic
terms. We did this once before with the energy density (\S\energyContinuumLimit).

@c
real YukawaDarkMatter::source(real mass, real r) {
  real iTerm;
  if (fabs(mass)<0.1) {
    real p = fermiMomentum(r);
    real logP = log(p);
    real logM = log(mass);
    iTerm = mass*(0.5*SQ(p)+0.25*(1.0-LN_4-2.0*logP+2.0*logM));
  } else {
    iTerm = CUBE(mass)*Util::i(fermiMomentum(r)/fabs(mass));
  }
  real fermionContribution = (coupling()/PI_SQ)*iTerm;
  return fermionContribution+partialScalarPotential(massToField(mass));
}

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
{\Gamma(4+3a)\over{\Gamma(1+3a)}}=3(1+a)(2+3a)(1+3a)\eqn{}
$$
hence
$$
p_{F}(r) = \left({{9\pi(1+a)(2+3a)(1+3a)N}\over{8}}\right)^{1/3}
\left(1-{r\over R}\right)^{a}{1\over R}.\eqn{}
$$
@^Fermi Momentum@>

@c
real YukawaDarkMatter::fermiMomentum(real r) {
  if (r > m_nuggetSize) return 0.0; /* avoid complex numbers! */
  if (m_a==0.5) return m_fermiMomentumCoef*sqrt(1.0-(r/m_nuggetSize));
  if (m_a==0.0) return m_fermiMomentumCoef;
  return m_fermiMomentumCoef*pow(1.0-(r/m_nuggetSize),m_a);
}

@ {\bf First-Order Equations of Motion.}
We claim that we may reduce the equations of motion to a first-order
differential equation:
$$
\phi'(r) = {1\over{r^{2}}}\int^{r}_{0}(r')^{2}\left[
{g_{\chi}\over \pi^{2}}m(r')^{3}i\left({{p_{F}(r')}\over\abs{m(r')}}\right)
+{{\partial V}\over{\partial\phi}}
\right]{\rm d}r'.\eqn{}
$$
We can see this from
$$
{{\rm d}\over{{\rm d}r}}[r^{2}\phi'(r)] = 
{g_{\chi}\over \pi^{2}}r^{2}m(r)^{3}i\left({{p_{F}(r)}\over\abs{m(r)}}\right)
+r^{2}{{\partial V}\over{\partial\phi}}\eqn{}
$$
which is precisely the equations of motion multipled by $r^{2}$. This
first-order ODE remains quite nonlinear.
\callthis\firstOrderFieldEqns

Let $V(m)=V(m_{\chi}+g_{\chi}\phi)$. Then we may rewrite this as a
first-order differential equation in the scaling mass
$$
{{\rm d}\over{{\rm d}r}}m(r) = 
{g_{\chi}\over{r^{2}}}\int^{r}_{0}(r')^{2}\left[
{g_{\chi}\over \pi^{2}}m(r')^{3}i\left({{p_{F}(r')}\over\abs{m(r')}}\right)
+{{\partial V}\over{\partial m(r)}}
\right]{\rm d}r'.\eqn{}
$$

@ {\bf Surface Boundary Conditions.}
For the massive scalar field, the boundary condition may be deduced from
$$
\phi(r) = -\phi(R)e^{-m_{\phi}r}{R\over r}\eqn{}
$$
for $r>R$. Thus we find
$$
\phi'(R) = -\phi(R)e^{-m_{\phi}R}\left({1\over R}+m_{\phi}\right)\eqn{}
$$
and for the massless case, we just take $m_{\phi}\to0$.
@^Boundary condition, surface@>

@c real YukawaDarkMatter::surfaceBoundaryCondition(real phi) {
   real mPhi = m_scalarMass;
   real R = m_nuggetSize;
   return -phi*exp(-mPhi*R)*(mPhi + (1.0/R));
}

@ {\bf Constraint: Positivity of Mass.}
We may implement a constraint condition, since
$$
m'(r) = g_{\chi}\phi'(r)\eqn{}
$$
we see
$$
m'(r)
=
{g_{\chi}\over{r^{2}}}\int^{r}_{0}(r')^{2}\left[
{g_{\chi}\over \pi^{2}}m(r')^{3}i\left({{p_{F}(r')}\over\abs{m(r')}}\right)
+{{\partial V}\over{\partial\phi}}
\right]{\rm d}r'.\eqn{}
$$
For the free situation $V=0$, we have $m'(r)>0$ and hence does not
change sign.
@^Constraint, Positive Mass@>
\callthis\constrainMassAsPositive

But for the massive situation $m_{\phi}\neq0$ or the self-interacting
scalar field, we cannot neglect the potential contribution. This
prevents us from generalizing the constraint to our particular
situation. However, the {\it physical} situation would have
$\abs{m(r)}\leq m_{\chi}$, so we just impose this constraint instead.

\rmk
If the sign of the effective mass changes, we have to carefully revisit
the ansatz for our Fermi momentum. Specifically, how it would behave at
that point.
@^Effective Mass, Sign Changes@>

@c
bool YukawaDarkMatter::isValidMass(real mass) {
     return (fabs(mass)<=m_fermionMass);
}

@ @<Model Param...@>=
        real m_fermionMass, m_nuggetSize, m_couplingConstant, m_fermionNumber;
        real m_scalarMass, m_a, m_fermiMomentumCoef;
@ @<Model constructor@>=
  YukawaDarkMatter(real massChi, real R, real coupling, real
        fermionNumber, real scalarMass=0.0, real a=0.0): 
  m_fermionMass(massChi), 
  m_nuggetSize(R),
  m_couplingConstant(coupling),
  m_fermionNumber(fermionNumber),
  m_scalarMass(scalarMass),
  m_a(a),
  m_fermiMomentumCoef(cbrt(9.0*PI_8*m_fermionNumber*(1.0+m_a)*(2.0+3.0*m_a)*(1.0+3.0*m_a))/m_nuggetSize)
   {}

@ @<Accessor Methods@>=
  real fermionMass() const { return m_fermionMass; }
  real nuggetSize() const { return m_nuggetSize; }
  real coupling() const { return m_couplingConstant; }
  real fermionNumber() const { return m_fermionNumber; }
  real scalarMass() const { return m_scalarMass; }

@ @<Mutator Methods@>=
  void setFermionNumber(real N) {
    if (!(N>0.0)) {
      throw std::logic_error("N must be positive");
    }
    m_fermionNumber = N;
  }
  void setNuggetSize(real R) {
    if (!(R>0.0)) {
      throw std::logic_error("R must be positive");
    }
    m_nuggetSize = R;
    m_fermiMomentumCoef=cbrt(9.0*PI_8*m_fermionNumber*(1.0+m_a)*(2.0+3.0*m_a)*(1.0+3.0*m_a))/R;
  }
  void setMomentumParameter(real a) {
    m_a = a;
  }

@ @<Routines@>=
  real scalarPotential(real phi);
  real partialScalarPotential(real phi);
  real fermiMomentum(real r);
  real energyDensity(real mass, real r);
  bool isValidMass(real mass);
  real surfaceBoundaryCondition(real phi);
  real source(real mass, real r);
  real scaleInvariantMass(real phi) const { return m_fermionMass + m_couplingConstant*phi; }
  real massToField(real mass) const { return (mass-m_fermionMass)/m_couplingConstant; }

@i algorithm.w

@ @<Include Header Files@>=
@<License@>@;
#include <limits>
#include <functional> /* for |std::function| template */
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#define VERBOSITY FATAL
#include "log.h"

#define VERSION "0.1 Beta"

@<Define Exception...@>@;
constexpr auto hbarC = 0.2; /* in GeV fm */
constexpr auto PI = 3.1415926535897932384626433832795L;
constexpr auto PI_4 = 0.78539816339744830961566084581988L;
constexpr auto PI_3 = 1.0471975511965977461542144610932L;
constexpr auto PI_8 = 0.39269908169872415480783042290994L;
constexpr auto FOUR_PI = 12.566370614359172953850573533118L;
constexpr auto SQRT_FOUR_PI = 3.5449077018110320545963349666823L;
constexpr auto PI_SQ = 9.8696044010893586188344909998762L;
typedef double real;
typedef std::size_t index;
real convertAlphaToCoupling(real alpha) { return SQRT_FOUR_PI*sqrt(alpha); }
real SQ(real x) { return x*x; }
real CUBE(real x) { return x*x*x; }

@ {\bf Exceptions.}
We have several different exceptions we need to define. One is if we
ever violate the bounds for the mass $0\leq m(r) \leq m_{\chi}$. This
should indicate we need start with a smaller initial guess $m_{0}$.

The other exception we must worry about is if we ever divide by
zero in |YukawaDarkMatter::source()| method, if
somehow $m(r)=0$. This {\it shouldn't} happen, and even if it did\dots
we'd start over since we violated the mass bound.

@<Define Exception Classes@>=
#include <stdexcept>
class MassOutOfBoundsException : public std::logic_error {
public:
  using std::logic_error::logic_error;
};

class NoSolutionExistsException : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

@* Test Cases.
We have three test cases prepared, which correspond to the $\alpha=0.1$
and constant Fermi momentum scenario. The results correspond to the
charts in the cited \.{arXiv} paper to within 7 digits. We may arbitrarily
extend the method to full machine precision, but unless we can narrow
the bounds of the initial guess\dots it will be slow if we demand
greater precision.

@c
YukawaDarkMatter* TestA() {
  real massChi = 100.0;
  real R = 0.160;
  real g = convertAlphaToCoupling(0.1);
  real N = sqrt(1e3);
  YukawaDarkMatter *model = new YukawaDarkMatter(massChi, R, g, N);
  return model;
}
YukawaDarkMatter* TestB() {
  real massChi = 100.0;
  real R = 0.62308;
  real g = convertAlphaToCoupling(0.1);
  real N = 1e3;
  YukawaDarkMatter *model = new YukawaDarkMatter(massChi, R, g, N);
  return model;
}
YukawaDarkMatter* TestC() {
  real massChi = 100.0;
  real R = 11.8776;
  real g = convertAlphaToCoupling(0.1);
  real N = 1e5;
  YukawaDarkMatter *model = new YukawaDarkMatter(massChi, R, g, N);
  return model;
}
YukawaDarkMatter* TestD() {
  real massChi = 100.0;
  real R = 0.5023;
  real g = convertAlphaToCoupling(0.01);
  real N = 1e3;
  YukawaDarkMatter *model = new YukawaDarkMatter(massChi, R, g, N);
  return model;
}

@i user.w

@ {\bf The main method.}
Now for the {\it pi\`ece de r\'esistance}: the |main()| method! For now,
I'll just examine certain test cases.

@c
void sampleMasses(Solver *solver) {
  real *m = solver->getMass();
  index length = (solver->getLength());
  int delta = length/10;
  for(int j = 0; j<10; j++) {
    std::cout<<"m["<<(j*(solver->dx())*delta)<<"] = "<<m[j*delta]<<std::endl;
  }
  std::cout<<"m["<<(length*(solver->dx()))<<"] = "<<m[length-1]<<std::endl;
}

void printVersionInfo() {
  std::cout<<"Welcome to the Yukawa Solver, v."<<VERSION<<std::endl;
}

void warnAboutSelfInteraction() {
  std::cout<<"This version handles massive scalar fields, but not self-interacting"<<std::endl;
  std::cout<<"scalar fields yet."<<std::endl;
}

void printStartupMessage() {
  std::cout<<std::endl;
  printVersionInfo();
  std::cout<<std::endl;
  warnAboutSelfInteraction();
  std::cout<<std::endl;
}

int main() {
  printStartupMessage();
  repl();
  return 0;
}

@* References.

\newdimen\refindent \refindent=27pt
\def\ref[#1] {\smallskip\noindent\hangindent\refindent
  \hbox to\refindent{\hfill[#1]\enspace}}

\ref[\Fr] Marco Frasca, ``Exact solutions of classical scalar field
equations''. {\sl J.Nonlin.Math.Phys.} {\bf 18} (2011)
291--297, \arXiv{0907.4053}. Discusses exact solution for $\phi^{4}$
classically. 
@^Frasca, Marco@>

\ref[\WZ] Mark Wise and Yue Zhang, ``Yukawa Bound States of a Large
Number of Fermions''. {\sl JHEP} {\bf 1502} no.2 (2015) 23, \arXiv{1411.1772}.
@^Wise, Mark@>
@^Zhang, Yue@>

@* License.
All sourcecode is developed under the MIT License (i.e., the ``Please
don't sue me if your computer blows up when running this code''
license):

\bigbreak
{\narrower
The {\sc MIT} License ({\sc MIT})

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
THE SOFTWARE.\par}

\let\oldSHC\SHC
\def\SHC#1{\noindent\ignorespaces$//\,${\eightrm #1}\smallbreak}
@<License@>=
// The MIT License (MIT)

// Copyright (c) 2015 by Alex Nelson

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

@* Index.