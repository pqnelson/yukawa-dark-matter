@* Algorithmic Structure.
The basic algorithm solves the first-order scalar field equations
(\S\firstOrderFieldEqns) in a pseudo-shooting manner. Since we know the
mass is always positive (\S\constrainMassAsPositive), we can use the
bisection method on the initial value for $m(r)$. We use the surface boundary
condition as the error of the approximation, and iterate accordingly.

How will this work? If we just use something simple like Euler's method,
we still have to compute a rather unpleasant integral. The trick is to
observe from
$$
m'(r) = {c_{0}\over{r^{2}}}\int^{r}_{0}(r')^{2}\sigma(r')\,{\rm d}r'
$$
where $\sigma(r)$ is the ``source function'' (i.e., right hand side of
the scalar field equations) we have
$$
m'(r) = {(r-h)^{2}\over{r^{2}}}m'(r-h) +
{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.
$$
We will translate this into a recursive algorithm involving $m_{n+1}$,
$m_{n}$, and $m_{n-1}$. First we must prove a few lemmas.

@c
class Solver {
private:
        YukawaDarkMatter *model;
        real *m_mass;
        index length;
        void iterate(index j);
        void solveScalarField();
        void shootingMethod();
        void bisectionMethod();
public:
        Solver(YukawaDarkMatter *m, index l): model(m), m_mass(new
        real[l]), length(l) {};
        ~Solver() { delete model; delete m_mass; }
        void run();
        real residual();
        real* getMass() const { return m_mass; }
        index getLength() const { return length; }
        real dx() const { return (model->nuggetSize())/length; }
};
        
@ {\bf Lemma.}
{\it We have one approximation of the derivative}
$$
{{-f(x+2h)+4f(x+h)-3f(x)}\over{2h}}=f'(x) - {{4h^{2}}\over{3!}}f'''(x) + \bigO{h^{3}}.
$$
\callthis\rhsDerivative

\proof
We first find the Taylor expansion of $f(x+2h)$ as
$$
f(x+2h) = f(x) + 2hf'(x) + {4h^{2}\over{2!}}f''(x) + {8h^{3}\over{3!}}f'''(x)
+\bigO{h^{4}}
$$
and likewise for $f(x+h)$ we get
$$
f(x+h) = f(x) + hf'(x) + {h^{2}\over{2!}}f''(x) + {h^{3}\over{3!}}f'''(x)
+\bigO{h^{4}}.
$$
Then subtracting gives us
$$
f(x+2h)-4f(x+h) = -3f(x)-2hf'(x)+{4h^{3}\over{3!}}f'''(x)+\bigO{h^{4}}.
$$
We can then produce the result immediately.\quad\slug

\rmk
We can get an approximation of order $\bigO{h^{4}}$ if we subtract out
$f'''(x)$ instead of $f''(x)$, then use the equations of motion for
$f''(x)$. As a first pass, we will not use this approach here.

\rmk
In practice, this doesn't work too well. It seems there's a problem I
cannot quite identify, and I'm unwilling to waste time figuring it
out. I'll put this on the back burner, maybe come back to it.

@ {\bf Lemma.}
{\it We have the approximation}
$$
{{f(x)-f(x-h)}\over h}=f'(x-h) + {h\over 2}f''(x-h)+\bigO{h^{2}}
$$

\proof
Immediate from the Taylor expansion for $f(x)$ about $x-h$.\quad\slug

@ {\bf Lemma.}
{\it We have the approximation}
$$
{{f(x+h)-f(x-h)}\over{2h}}=f'(x) + {h^{2}\over3}f'''(x)+\bigO{h^{3}}.
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
{{m(r)-m(r-h)}\over{h}}
\right)
+ {c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.
$$
or equivalently
$$
{m(r+h)-m(r-h)} = 2{(r-h)^{2}\over{r^{2}}}\left(
{m(r)-m(r-h)}
\right)
+ 2h{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.
$$
Rearranging terms, we see
$$
m(r+h)
= 2{(r-h)^{2}\over{r^{2}}}m(r)
+\left[1-2{(r-h)^{2}\over{r^{2}}}\right]m(r-h)
+ 2h{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.
$$
We are left with trying to approximate the integral expression before we
can rest.

@ {\bf Approximating the Integral.}
We approximate the integral using Simpson's rule:
$$
\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'=
{h\over6}\bigl((r-h)^{2}\sigma(r-h) + (r + 0.5h)^2\sigma(r+0.5h) +
r^{2}\sigma(r)\bigr)
 - \left.{h^{4}\over{2880}}
    {{\rm d}^{4}\over{{\rm d}r^{4}}}(r^{2}\sigma(r))\right\evalAt_{r=\xi}
$$
for some $\xi\in(r-h,r)$.

@ {\bf Iterative Procedure.}
We can combine the previous steps to conclude
$$
m(r+h)
= 2{(r-h)^{2}\over{r^{2}}}m(r)
+\left[1-2{(r-h)^{2}\over{r^{2}}}\right]m(r-h)
+ h{c_{0}\over{r^{2}}}{h\over3}\bigl((r-h)^{2}\sigma(r-h) + (r + 0.5h)^2\sigma(r+0.5h) +
r^{2}\sigma(r)\bigr).
$$
Or, writing $r_{n}=nh$, and $m_{n}=m(r_{n})$, we get
$$
m_{n+1}
= {{r_{n-1}^{2}}\over{r_{n}^{2}}}2 m_{n}
+\left(1-2{{r_{n-1}^{2}}\over{r_{n}^{2}}}\right)m_{n-1}
+ h{c_{0}\over{r_{n}^{2}}}{h\over3}\bigl(r_{n-1}^{2}\sigma(r_{n-1})
 + \left({{r_{n} + r_{n-1}}\over 2}\right)^2\sigma\left({{r_{n} + r_{n-1}}\over 2}\right)
+ r_{n}^{2}\sigma(r_{n})\bigr).
$$

\rmk
When we modify the scalar field to include some nonzero $m_{\phi}$, we
just have to modify the |YukawaDarkMatter::source()| routine to include
the extra $-m_{\phi}^{2}\phi(r)$ term.

@c void Solver::iterate(index j) {
   if(j<2) return;
   real h = dx();
   real r = (j-1)*h;
   real rSq = SQ(r);
   real rPrime = (j-2)*h;
   real rPrimeSq = SQ(rPrime);
   real u = rPrime/r; 
   real c = (model->coupling())/rSq;

   real massTerms = 2.0*u*m_mass[j-1] + (1.0 - 2.0*u)*m_mass[j-2];
   // Simpson's rule for quadrature
   real k[3];
   k[0] = rPrimeSq*(model->source(m_mass[j-2], rPrime));
   real m = 0.5*(m_mass[j-2] + m_mass[j-1]);
   u = 0.5*(r + rPrime);
   k[1] = SQ(u)*(model->source(m, u));
   k[2] = rSq*(model->source(m_mass[j-1], r));

   real integral = (h/6.0)*(k[0] + 4.0*k[1] + k[2]);

   m_mass[j] = massTerms + c*h*integral;
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

@c
void Solver::bisectionMethod() {
     real massLowerBound = 0.0;
     real massUpperBound = model->fermionMass();
     real m, res;
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
m_{0}(N=10^{3})=2.21073,\quad{\rm and}\quad m_{0}(N=10^{3/2})=72.2797.
$$
We guess a power law behavior $m_{0}(N) = c_{N}N^{b}$. We find
$$
b = {2\over 3}\log_{10}\left({{2.21073}\over{72.2797}}\right)=-1.00965.
$$
Hence $m_{0}\sim N^{-1}$, and $C_{N}\approx 2210$.

We can also figure out the same behaviour for
$m_{0}$'s dependence on the coupling $\alpha=4\pi g_{\chi}^{2}$. We find
empirically:
$$
m_{0}(\alpha=0.1,N=10^{3})=2.21073,\quad{\rm and}\quad
m_{0}(\alpha=0.01,N=10^{3})=72.7375.
$$
We again suppose $m_{0}(\alpha)=c_{\alpha}\alpha^{b'}$. We find
$$
b' = -\log_{10}\left({{72.7375}\over{2.21073}}\right)\approx-1.51722.
$$
So we find $m_{0}\sim\alpha^{-3/2}$. We find $c_{\alpha}\approx 0.07$.

@ {\bf Shooting Method.}
We now iteratively determine the solution. The method so far is
incredibly naive, readjusting the initial position based on the sign of
the residual of the surface boundary condition. One method to speed this
up is to come up with a linear approximation of the residual for the
first couple iterations, the solve this approximation for the mass which
would make the residual vanish. So instead of 30 or more iterations,
it'd boil down to---say---5 or so.

If $m_{0}$ produces a residual $\rho_{0}$, and $m_{1}$ produces residual
$\rho_{1}$, then we have the linear approximation
$$
\rho(m) = \rho_{0} + \left({{\rho_{1}-\rho_{0}}\over{m_{1}-m_{0}}}\right)(m-m_{0}).
$$
We want to find the $m$ such that $\rho(m)=0$. We find, by basic
algebra, that
$$
m=m_{0}-\rho_{0}\left({{m_{1}-m_{0}}\over{\rho_{1}-\rho_{0}}}\right).
$$
This gives us a way to find the solution faster (in theory).

@c
void Solver::shootingMethod() {
  const real TOL = 1e-7;
  real massLowerBound = 0.0;
  real massUpperBound = model->fermionMass();
  real m = 0.0, residual_ =0.0;
  real masses[2];
  real res[2];
  int k=0;
  for(int j=0; j<2; j++) {
    @<Set the initial condition@>@;
    @<Solve the System Once@>@;
    masses[j] = m;
    residual_ = res[j] = residual();
    if(res[j]>0.0) massUpperBound = m;
    if(res[j]<0.0) massLowerBound = m;
    if(j==0 && (massUpperBound==model->fermionMass())) {
      @<Adjust the bounds if needed@>@;
    }
  }
  for(k=0; k<10; k++) {
    @<Solve System with Extrapolated...@>@;
    @<Update the Residual and Masses@>@;
  }
}

@ @<Update the Residual and Masses@>=
  masses[0] = masses[1];
  res[0] = res[1];
  masses[1] = m;
  residual_ = res[1] = residual();
  @<Adjust the bounds if needed@>@;

@ @<Solve System with Extrapolated Mass@>=
  m = masses[0] - res[0]*((masses[1]-masses[0])/(res[1]-res[0]));
  @<Fallback to Bisection@>@;
  m_mass[0] = m;
  m_mass[1] = m;
  @<Solve the System Once@>@;

@ If something catastrophic happens, and our linear extrapolation takes
us to somewhere {\it worse} than our previous guesses, we should just
resort to the bisection method.

@<Fallback to Bisection@>=
  if(!(massLowerBound<m && m<massUpperBound)) {
    std::cout<<"[WARN] Falling back to bisection method"<<std::endl; 
    m = 0.5*(massLowerBound + massUpperBound);
    k--;
  }

@ @<Adjust the bounds if needed@>=
  if(residual_>TOL && m<massUpperBound) {
    massUpperBound = m;
  } else if (residual_<-TOL && m>massLowerBound) {
    massLowerBound = m;
  } else if (fabs(residual_)<TOL) {
    break;
  }

@ @<Solve the System Once@>=
     do {
       try {
         solveScalarField();
         break;
       } catch (MassOutOfBoundsException e) {
         massUpperBound = m;
         @<Set the initial condition@>@;
       }
     } while(true);

@ @c
void Solver::run() {
     shootingMethod();
}

@ @<Set the initial condition@>=
  if(m>0.0) {
    m = 0.5*(massUpperBound + massLowerBound);
  } else {
    m = 1.0;
  }
  m_mass[0] = m;
  m_mass[1] = m;

@ We then consider the residual of the surface boundary conditions,
which indicates ``how far off'' we are.

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
  std::cout<<"Surface boundary condition gives: "<<dr<<std::endl;
  return dr;
}  
