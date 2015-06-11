@** Algorithmic Structure.
The basic algorithm solves the first-order scalar field equations
(\S\firstOrderFieldEqns) in a pseudo-shooting manner. Since we know the
mass is always positive (\S\constrainMassAsPositive), we can use the
bisection method on the initial value for $m(r)$. We use the surface boundary
condition as the error of the approximation, and iterate accordingly.

How will this work? If we just use something simple like Euler's method,
we still have to compute a rather unpleasant integral. The trick is to
observe from
$$
m'(r) = {c_{0}\over{r^{2}}}\int^{r}_{0}(r')^{2}\sigma(r')\,{\rm d}r'\eqn{}
$$
where $\sigma(r)$ is the ``source function'' (i.e., right hand side of
the scalar field equations) we have
$$
m'(r) = {(r-h)^{2}\over{r^{2}}}m'(r-h) +
{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.\eqn{}
$$
We will translate this into a recursive algorithm involving $m_{n+1}$,
$m_{n}$, and $m_{n-1}$. First we must prove a few lemmas.

@c
class Solver {
private:
        @<Solver's Fields@>@;
        @<Solver's Internal Routines@>@;
public:
        @<Solver Constructor and Destructor@>@;
        void run();
        @<Solver's Accessor Methods@>@;
        @<Solver's Computational Routines@>@;
};

@ {\bf Lemma.}
{\it We have one approximation of the derivative}
$$
{{-f(x+2h)+4f(x+h)-3f(x)}\over{2h}}=f'(x) - {{4h^{2}}\over{3!}}f'''(x) + \bigO{h^{3}}.\eqn{}
$$
\callthis\rhsDerivative

\proof
We first find the Taylor expansion of $f(x+2h)$ as
$$
f(x+2h) = f(x) + 2hf'(x) + {4h^{2}\over{2!}}f''(x) + {8h^{3}\over{3!}}f'''(x)
+\bigO{h^{4}}\eqn{}
$$
and likewise for $f(x+h)$ we get
$$
f(x+h) = f(x) + hf'(x) + {h^{2}\over{2!}}f''(x) + {h^{3}\over{3!}}f'''(x)
+\bigO{h^{4}}.\eqn{}
$$
Then subtracting gives us
$$
f(x+2h)-4f(x+h) = -3f(x)-2hf'(x)+{4h^{3}\over{3!}}f'''(x)+\bigO{h^{4}}.\eqn{}
$$
We can then produce the result immediately.\quad\slug

\rmk
We can get an approximation of order $\bigO{h^{4}}$ if we subtract out
$f'''(x)$ instead of $f''(x)$, then use the equations of motion for
$f''(x)$. As a first pass, we will not use this approach here.

@ {\bf Lemma.}
{\it We have the approximation}
$$
{{f(x)-f(x-h)}\over h}=f'(x-h) + {h\over 2}f''(x-h)+\bigO{h^{2}}\eqn{}
$$

\proof
Immediate from the Taylor expansion for $f(x)$ about $x-h$.\quad\slug

@ {\bf Lemma.}
{\it We have the approximation}
$$
{{f(x+h)-f(x-h)}\over{2h}}=f'(x) + {h^{2}\over3}f'''(x)+\bigO{h^{3}}.\eqn{}
$$
{\it This will be used to approximate $m'(r)$.}
\callthis\lhsDerivative

\rmk
This is probably the biggest source of truncation error. But given that
we only want to use $m_{n-1}$, $m_{n}$, and $m_{n+1}$, this is the best
we can do under these constraints. If we need greater precision, we need
to swap this out for a different approximation.

@ {\bf Basic Algorithm.}
We can combine these lemmas to produce an approximation
$$
{{m(r+h)-m(r-h)}\over{2h}} = {(r-h)^{2}\over{r^{2}}}\left(
{{-m(r+h)+4m(r)-3m(r-h)}\over{2h}}
\right)
+ {c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.\eqn{}
$$
or equivalently
$$
{m(r+h)-m(r-h)} = {(r-h)^{2}\over{r^{2}}}\left(
{-m(r+h)+4m(r)-3m(r-h)}
\right)
+ 2h{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.\eqn{}
$$
Rearranging terms, we see
$$
(1+u)m(r+h)
= 4{(r-h)^{2}\over{r^{2}}}m(r)
+\left[1-3{(r-h)^{2}\over{r^{2}}}\right]m(r-h)
+ 2h{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'\eqn{}
$$
where $u=(r-h)^{2}/r^2$.
We are left with trying to approximate the integral expression before we
can rest.

@ {\bf Approximating the Integral.}
We approximate the integral using Simpson's rule. Recall...

\proclaim Simpson's Rule. If $f\in C^{4}(a,b)$, then we may approximate
the integral
$$
\int^{b}_{a}f(x)\,{\rm d}x=
{h\over6}\left(f(a) + f\left({{a+b}\over2}\right)+f(b)\right)
 - \left.{h^{4}\over{2880}}
    {{\rm d}^{4}\over{{\rm d}x^{4}}}f(x)\right\evalAt_{r=\xi}\eqn{}
$$
for some $\xi\in(a,b)$.
@^Simpson's Rule@>

We can apply this directly to our problem, to find:
$$
\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'=
{h\over6}\bigl((r-h)^{2}\sigma(r-h)
+ (r - 0.5h)^2\sigma(r-0.5h)
+ r^{2}\sigma(r)\bigr)
 - \left.{h^{4}\over{2880}}
    {{\rm d}^{4}\over{{\rm d}r^{4}}}(r^{2}\sigma(r))\right\evalAt_{r=\xi}\eqn{}
$$
for some $\xi\in(r-h,r)$.


@* Iterative Procedure.
We can combine the previous steps to conclude
$$
\eqalign{m(r+h)
=& {4u\over{1+u}}m(r)
+\left[{{1-3u}\over{1+u}}\right]m(r-h)\cr
&+ {c_{0}h^{2}\over{3\cdot r^{2}(1+u)}}\left[(r-h)^{2}\sigma(r-h) + (r + 0.5h)^2\sigma(r+0.5h) +
r^{2}\sigma(r)\right].\cr}\eqn{}
$$
Or, writing $r_{n}=nh$, and $m_{n}=m(r_{n})$, we get
$$
\eqalign{m_{n+1}
=& {4u\over{1+u}}m_{n}
+\left[{{1-3u}\over{1+u}}\right]m_{n-1}\cr
&+ {c_{0}h^{2}\over{3\cdot r_{n}^{2}(1+u)}}\left[r_{n-1}^{2}\sigma(r_{n-1})
+ \left({{r_{n} + r_{n-1}}\over 2}\right)^2\sigma\left({{r_{n} + r_{n-1}}\over 2}\right)
+ r_{n}^{2}\sigma(r_{n})\right].\cr}\eqn{}
$$

@c void Solver::iterate(index j) {
  if(j<2) return;
  real h = dx();
  real r = (j-1)*h;
  real rSq = SQ(r);
  real rPrime = (j-2)*h;
  real rPrimeSq = SQ(rPrime);
  real u = (rPrime/r);
  real c = (model->coupling())/rSq;
  real v = 1.0/(1.0+u);
  real massTerms = 4.0*u*v*m_mass[j-1] + v*(1.0-3.0*u)*m_mass[j-2];
  // Simpson's rule for quadrature
  real k[3];
  k[0] = rPrimeSq*(model->source(m_mass[j-2], rPrime));
  real m = 0.5*(m_mass[j-2] + m_mass[j-1]);
  u = 0.5*(r + rPrime);
  k[1] = SQ(u)*(model->source(m, u));
  k[2] = rSq*(model->source(m_mass[j-1], r));

  real integral = (h*ONE_SIXTH)*(k[0] + 4.0*k[1] + k[2]);

  m_mass[j] = massTerms + c*h*v*integral;
}

@ {\bf Solving the mass differential equation.}
We simply iterate through the $m_{n}$ for $n=2,\dots,N-1$ to solve the
differential equation at hand. If we get an invalid mass, we should
throw an exception indicating the bad state.

@c
void Solver::solveScalarField() {
  for(index j=0; j<length; j++) {
    iterate(j);
    if(!model->isValidMass(m_mass[j])) {
      throw MassOutOfBoundsException("");
    }
  }
}


@ {\bf Bisection Method.}
A far more cautious approach, because look: we know $0\leq m(0)\leq
m_{\chi}$, so why not just use the bisection method and iteratively
solve the system 30 times?
@^Bisection method@>

@c
void Solver::bisectionMethod() {
     real massLowerBound = 0.0;
     real massUpperBound = model->fermionMass();
     real m = 0.1*massUpperBound, res;
     for(index j=0; j<30+(index)(log2(model->fermionMass())+0.5); j++) {
       @<Set the initial condition@>@;
       try {
         solveScalarField();
         res = residual();
         if(res>1e-7) {
           massUpperBound = m;
         } else if (res<-1e-7) {
           massLowerBound = m;
         } else {
           break;
         }
       } catch (MassOutOfBoundsException e) {
         massUpperBound = m;
         j--;
       }
     }
}

@ {\bf Initial Guess.}
We can do some approximations to speed things up by improving the
initial guess. We can measure empirically the initial guess $m_{0}$
which solves various test cases, then extrapolate out some power
law. When the Fermi momentum is constant, we find the following
situation for $\alpha=0.1$:
$$
m_{0}(N=10^{3})=2.21073,\quad{\rm and}\quad m_{0}(N=10^{3/2})=72.2797.\eqn{}
$$
We guess a power law behavior $m_{0}(N) = c_{N}N^{b}$. We find
$$
b = {2\over 3}\log_{10}\left({{2.21073}\over{72.2797}}\right)=-1.00965.\eqn{}
$$
Hence $m_{0}\sim N^{-1}$, and $C_{N}\approx 2210$.

We can also figure out the same behaviour for
$m_{0}$'s dependence on the coupling $\alpha=4\pi g_{\chi}^{2}$. We find
empirically:
$$
m_{0}(\alpha=0.1,N=10^{3})=2.21073,\quad{\rm and}\quad
m_{0}(\alpha=0.01,N=10^{3})=72.7375.\eqn{}
$$
We again suppose $m_{0}(\alpha)=c_{\alpha}\alpha^{b'}$. We find
$$
b' = -\log_{10}\left({{72.7375}\over{2.21073}}\right)\approx-1.51722.\eqn{}
$$
So we find $m_{0}\sim\alpha^{-3/2}$. We find $c_{\alpha}\approx 0.07$.

@* Find Solution Satisfying the Boundary Conditions.
We now iteratively determine the solution. The method so far is
incredibly naive, readjusting the initial position based on the sign of
the residual of the surface boundary condition. One method to speed this
up is to come up with a linear approximation of the residual for the
first couple iterations, the solve this approximation for the mass which
would make the residual vanish. So instead of 30 or more iterations,
it'd boil down to---say---5 or so.
@^Shooting method@>

If $m_{0}$ produces a residual $\rho_{0}$, and $m_{1}$ produces residual
$\rho_{1}$, then we have the linear approximation
$$
\rho(m) = \rho_{0} + \left({{\rho_{1}-\rho_{0}}\over{m_{1}-m_{0}}}\right)(m-m_{0}).\eqn{}
$$
We want to find the $m$ such that $\rho(m)=0$. We find, by basic
algebra, that
$$
m=m_{0}-\rho_{0}\left({{m_{1}-m_{0}}\over{\rho_{1}-\rho_{0}}}\right).\eqn{}
$$
This gives us a way to find the solution faster (in theory).

@c
void Solver::shootingMethod() {
  const real TOL = 1e-12;
  real massLowerBound = 0.0;
  real massUpperBound = model->fermionMass();
  real m = 0.0, residual_ =0.0;
  real residualUpper = 10000.0;
  real residualLower = -10000.0;
  real masses[2];
  real res[2];
  int k=0;
  int MAX_BISECTION = 30;
  for(int j=0; j<2; j++) {
    @<Set the initial condition@>@;
    @<Solve the System Once@>@;
    masses[j] = m;
    residual_ = res[j] = residual();
    if(res[j]>0.0) {
      massUpperBound = m;
      residualUpper = res[j];
    }
    if(res[j]<0.0) {
      massLowerBound = m;
      residualLower = res[j];
    }
    if(j==0 && (massUpperBound==model->fermionMass())) {
      @<Adjust the bounds if needed@>@;
    }
  }
  for(k=0; k<12; k++) {
    @<Solve System with Extrapolated...@>@;
    @<Update the Guess from the Residual and Masses@>@;
    if (residual_>TOL && residual_<residualUpper) {
      residualUpper = residual_;
      massUpperBound = m;
    } else if (residual_<-TOL && residualLower < residual_) {
      massLowerBound = m;
      residualLower = residual_;
    } else {
      break;
    }
  }
}

@ @<Update the Guess from the Residual and Masses@>=
  LOG::trace<<"Updating the initial condition guess from the residual"<<std::endl;
  masses[0] = masses[1];
  res[0] = res[1];
  masses[1] = m;
  residual_ = res[1] = residual();
  @<Adjust the bounds if needed@>@;

@ Use the secant method to approximate the next guess for the initial
condition.

@<Solve System with Extrapolated Mass@>=
  LOG::trace<<"Solving the system with the new initial condition"<<std::endl;
  real deltaM = res[0]*((masses[1]-masses[0])/(res[1]-res[0]));
  m = masses[0] + (residual_<0.0 ? -1.0 : 1.0)*fabs(deltaM);
  @<Fallback to Bisection@>@;
  m_mass[0] = m;
  m_mass[1] = m;
  @<Solve the System Once@>@;

@ If something catastrophic happens, and our linear extrapolation takes
us to somewhere {\it worse} than our previous guesses, we should just
resort to the bisection method. We have an upper bound on the number of
times we can default to bisection, approximately 30 times will give us
precision to 7 digits or so.

@<Fallback to Bisection@>=
  if(!(massLowerBound<m && m<massUpperBound)) {
    if (MAX_BISECTION<0) break;
    MAX_BISECTION--;
    LOG::warn<<"Falling back to bisection method"<<std::endl;
    m = 0.5*(massLowerBound + massUpperBound);
    k--;
  }

@ @<Adjust the bounds if needed@>=
  LOG::trace<<"Checking to see if we need to adjust the bounds on the guess for the initial condition"<<std::endl;
  if(residual()>TOL && m<massUpperBound) {
    massUpperBound = m;
  } else if (residual()<-TOL && m>massLowerBound) {
    massLowerBound = m;
  } else if (fabs(residual())<TOL) {
    return;
  }

@ @<Solve the System Once@>=
  LOG::trace<<"Iteratively solving the system (slow)..."<<std::endl;
  while (true) {
    try {
      solveScalarField();
      break;
    } catch (MassOutOfBoundsException e) {
      LOG::info<<"Guessed too high (Out of bounds exception caught)"<<std::endl;
      if (fabs(m)<1e-21) {
        throw NoSolutionExistsException("Cannot find a nonzero solution");
      } else {
        massUpperBound = m;
      }
      @<Set the initial condition@>@;
    }
  }
  LOG::trace<<"Solution found!"<<std::endl;

@ @<Set the initial condition@>=
  m = 0.5*(massUpperBound + massLowerBound);
  if (m==0.0) {
    m += 1e-7;
  }
  if (false) {
    throw NoSolutionExistsException("Solver::shootingMethod() initial guess smaller than machine epsilon");
  }
  LOG::info<<"Guessing again with m = "<<m<<std::endl;
  m_mass[0] = m;
  m_mass[1] = m;

@ @c
void Solver::run() {
  LOG::trace<<"Solver::run() called..."<<std::endl;
  shootingMethod();
  LOG::trace<<"Solver::run() terminating..."<<std::endl;
}

@ We then consider the residual of the surface boundary conditions,
which indicates ``how far off'' we are.

Note, this is a good proxy for the overall residual since the uniqueness
theorem for differential equations guarantees a function satisfying the
boundary conditions to be uniquely determined. Since we have the
$\phi'(0)=0$ trivially solved, the surface boundary condition's residual
would indicate ``how far off'' we are from the honest solution.

@c
real Solver::residual() {
  real h = dx();
  real dr;
  real field[3];
  field[0] = model->massToField(m_mass[length-3]);
  field[1] = model->massToField(m_mass[length-2]);
  field[2] = model->massToField(m_mass[length-1]);      
  real derivativeOnSurface = (0.5*(field[0]-4.0*field[1]+3.0*field[2])/h);
  dr = derivativeOnSurface - model->surfaceBoundaryCondition(field[2]);
  return dr;
}  

@* Computing the Energy.
A critical component to our analysis is computing the energy for a given
radius. For the time being, we will just do the simplest thing to
program: Simpson's rule. We can skip computing the energy density at
$r=0$ since the integrand includes a factor of $r^{2}$, hence vanishes.
@^Energy, computing@>

@c
real Solver::computeEnergy() {
  real energy = 0.0;
  real h = dx();
  for(index j=1; j<(length/2)-1; j++) {
    energy += 4.0*(model->energyDensity(m_mass[2*j-1], (2*j-1)*h));
    energy += 2.0*(model->energyDensity(m_mass[2*j], 2*j*h));
  }
  energy += 4.0*(model->energyDensity(m_mass[length-2], (length-2)*h));
  energy += (model->energyDensity(m_mass[length-1], (length-1)*h));
  return energy*h*ONE_SIXTH;
}

@ {\bf Dump Energy.}
This is a simple method that computes the $(R,E)$ --- the nugget size
$R$, and corresponding energy $E$ --- ranging over $0\leq R\leq 5$
stepping $\Delta R$.

@c
void Solver::dumpEnergy() {
     real E = 0.0;
     real dR = 0.05;
     for(index j=1; j<30; j++) {
         real R = dR*j;
         model->setNuggetSize(R);
         run();
         E = computeEnergy();
         std::cout<<"E("<<R<<") = "<<E<<std::endl;
     }
}

@* Determining the Nugget Size.
We want to find the nugget size which minimizes the energy
(\S\energyContinuumLimit).

One approach is to just ``walk'' along the values of $R$, and determine
when the energy is increasing. Then go back to the last place where it
was decreasing, and walk smaller steps. Keep iterating until you're
satisfied. This is a ``bisection-like method'', and incredibly slow.
@^Nugget Size, determining@>

The clever approach is to sample energy points at the current nugget
size $R_{1}=R$, but also at $R_{0}=R-h$ and $R_{2}=R+h$ for some
$0<h\ll R$. Then construct a quadratic polynomial
$$
E(R) = E_{0} + E'_{0}(R-R_{0}) + {E''_{0}\over 2!}(R-R_{0})^{2}.
$$
Take its derivative, solve for $E'(R)=0$. This will give us our next
guess for $R$ which {\it extremizes} the energy. Since the energy
appears to be a ``parabolic'' function of nugget size, this would
produce the minima desired.

\proclaim Puzzle. When $m_{\phi}=0$, the energy as a function of nugget
size $E(R)$ is ``parabolic''. But is this still true for $m_{\phi}\neq0$?

@c
void Solver::findNuggetSize() {
  LOG::trace<<"Solver::findNuggetSize() called..."<<std::endl;
  real h = 0.1;
  @<Initial Guess for Nugget Size@>@;
  real nextR;
  real dE[2];
  real E[3];
  int silverBullets = 5;
  for(int j=0; j<20; j++) {
    try {
      h = initialH/(j+1.0);
      @<Sample Local Energy@>@;
      @<Compute Radial Derivatives of Energy@>@;
      @<Determine Next Guess for Nugget Size@>@;
      LOG::info<<std::setprecision(20)
               <<"R = "
               <<nextR<<std::endl;
      @<Terminate Nugget Iterative Loop if Good Enough@>@;
      R = nextR;
    } catch (NoSolutionExistsException e) {
      @<Handle no solution for given radius@>@;
    } catch (NegativeDistanceException e) {
      R = -10.0*R;
      LOG::info<<"Solver::findNuggetSize() Guessed a negative distance, trying again with R = "
               <<R<<std::endl;
      if (silverBullets>0) {
        silverBullets--;
        j=0;
      }
    }
  }
  LOG::info<<"Solver::findNuggetSize() Setting R = "<<std::setprecision(20)<<R<<std::endl;
  model->setNuggetSize(R);
  LOG::trace<<"Solver::findNuggetSize() terminating..."<<std::endl;
}

@ {\bf Todo: Pick next radius when we cannot find a solution.}
This is the biggest problem I need to think deeply about: when the
solver cannot find a solution for the given radius, how should this be
handled? We can pick a bigger radius (or a smaller one), or do something
else.
@^TODO: Minimize Energy Problem@>

@<Handle no solution for given radius@>=
      if (R<0.25) {
        R *= 2.0;
      } else {
        R += 0.25;
      }
      LOG::info<<"Solver::findNuggetSize() No such solution exists, guessing again with R = "
               <<R<<std::endl;
      j--;


@ We know that $R\sim 1/m_{\chi}$, but after a few numerical experiments
I found that $R\approx 16\sqrt{N}\hbar c/m_{\chi}$ is a decent
approximation when $m_{\phi}=0$.
\callthis\initialGuessForNuggetSize

@<Initial Guess for Nugget Size@>=
  real R = 16.0*sqrt(model->fermionNumber())*(hbarC/(model->fermionMass()));
  if (R<h) {
    h = 0.125*R;
  }
  real initialH = h;
  LOG::trace<<"Solver::findNuggetSize() Initial guess for R: "
            <<std::setprecision(20)
            <<R<<std::endl;

@ @<Sample Local Energy@>=
    for(int k=0; k<3; k++) {
      if (R<h) h = 0.25*R;
      model->setNuggetSize(R+((k-1)*h));
      run();
      E[k] = computeEnergy();
      LOG::info<<"Solver::findNuggetSize() E["<<k<<"] = "
               <<std::setprecision(20)
               <<E[k]<<std::endl;
    }

@ We approximate the radial derivatives |dE| by using the forward
difference. We take their average to get the centered difference and
store that as |derivative|.

@<Compute Radial Derivatives of Energy@>=
    dE[0] = (E[1]-E[0])/h;
    dE[1] = (E[2]-E[1])/h;
    real dSqE = (dE[1]-dE[0])/h;
    real derivative = 0.5*(dE[0]+dE[1]);
    
@ Based on whether the change in energy is increasing or decreasing, we
adjust our guesses accordingly. We use the local finite differences to
construct a quadratic polynomial
$$
E(R) = E_{0} + {{\Delta E_{0}}\over{\Delta R_{0}}}(R-R_{0})
+{1\over2}{{\Delta^{2}E_{0}}\over{\Delta^{2}R_{0}}}(R-R_{0})(R-R_{1})\eqn{}
$$
We take its derivative
$$
E'(R) = {{\Delta E_{0}}\over{\Delta R_{0}}}
+{1\over2}{{\Delta^{2}E_{0}}\over{\Delta^{2}R_{0}}}(2R-R_{0}-R_{1})\eqn{}
$$
then setting it to zero we find
$$
R={R_{0}+R_{1}\over2}-\left({{\Delta^{2}E_{0}}\over{\Delta^{2}R_{0}}}\right)^{-1}{{\Delta E_{0}}\over{\Delta R_{0}}}.\eqn{}
$$
This is how we determine our next guess.

@<Determine Next Guess for Nugget Size@>=
    if ((dE[0]>0.0 && dE[1]>0.0) || (dE[0]<0.0 && dE[1]<0.0)) {
      nextR = R - 0.5*h - (dE[0]>0.0 ? 1.0 : -1.0)*fabs(derivative/dSqE);
      if (nextR < 0.0) {
        nextR = 0.5*R;
      }
    } else if (dE[0]>0.0 && dE[1]<0.0) {
      nextR = R;
      LOG::error<<"Solver::findNuggetSize() Derivatives have incorrect signs..."<<std::endl;
    } else {
      h *= 0.5;
      nextR = R + h;
    }

@ @<Terminate Nugget Iterative Loop if Good Enough@>=
    if (fabs(nextR-R)<1e-14) {
      LOG::info<<" terminating approximation for R, next R too close"<<std::endl;
      R = nextR;
      break;
    }

@ {\bf Backmatter.}
Now that we are ``done'' with the algorithms, I will just dump the
various components of the |Solver| class here.

@ @<Solver's Fields@>=
  YukawaDarkMatter *model;
  real *m_mass;
  index length;

@ @<Solver's Internal Routines@>=
  void iterate(index j);
  void solveScalarField();
  void shootingMethod();
  void bisectionMethod();

@ @<Solver Constructor and Destructor@>=
  Solver(YukawaDarkMatter *m, index l):
    model(m),
    m_mass(new real[l]),
    length(l)
    {};
  ~Solver() { delete m_mass; }

@ @<Solver's Accessor Methods@>=
  real* getMass() const { return m_mass; }
  index getLength() const { return length; }
  real dx() const { return (model->nuggetSize())/length; }

@ @<Solver's Computational Routines@>=
  real residual();
  void dumpEnergy();
  void findNuggetSize();
  real computeEnergy();
  