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
        void trickyIterate(index j);
        void rungeKuttaIterate(index j);
        void solveScalarField();
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
{\it In particular, this will be useful to approximate $m'(r-h)$.}
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
{{-m(r+h)+4m(r)-3m(r-h)}\over{2h}}
\right)
+ {c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.
$$
or equivalently
$$
{m(r+h)-m(r-h)} = {(r-h)^{2}\over{r^{2}}}\left(
{-m(r+h)+4m(r)-3m(r-h)}
\right)
+ 2h{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.
$$
Rearranging terms, we see
$$
\left[1+{(r-h)^{2}\over{r^{2}}}\right]m(r+h)
= 4{(r-h)^{2}\over{r^{2}}}m(r)
+\left[1-3{(r-h)^{2}\over{r^{2}}}\right]m(r-h)
+ 2h{c_{0}\over{r^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'
$$
and thus
$$
m(r+h)
= {{(r-h)^{2}}\over{r^{2}+(r-h)^{2}}}4 m(r)
+{{r^{2}}\over{r^{2}+(r-h)^{2}}}\left[1-3{(r-h)^{2}\over{r^{2}}}\right]m(r-h)
+ 2h{c_{0}\over{r^{2}+(r-h)^{2}}}\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'.
$$
We are left with trying to approximate the integral expression before we
can rest.

@ {\bf Approximating the Integral.}
We see the truncation error, so far, is quadratic in $h$ (this is
because the left hand side's approximation of $m'(r)$ has
$h^{2}m'''(r)$ as the truncation error---c.f., \S\lhsDerivative). We
already have the coefficient of the integral term involve a factor of
$h/r^{2}$, modulo some arbitrary constant. But iteratively,
$r=r_{n}=nh$. So the coefficient is $h/r^{2}\sim\bigO{h^{-1}}$. We need
an approximation that is at least cubic in $h$ for the truncation
error.

The trapezoid rule is one such approximation:
$$
\int^{r}_{r-h}(r')^{2}\sigma(r')\,{\rm d}r'=
{h\over2}\bigl((r-h)^{2}\sigma(r-h) + r^{2}\sigma(r)\bigr) - {h^{3}\over{12}}\left(2\sigma(\xi)+4\xi\sigma'(\xi)+\xi^{2}\sigma''(\xi)\right)
$$
for some $\xi\in(r-h,r)$.

\rmk
For a higher order approximation, we may have to use Simpson's rule and
possibly some sort of Taylor expansion of $m(r - h/2)$ about $m(r-h)$ or
$m(r)$.

@ {\bf Iterative Procedure.}
We can combine the previous steps to conclude
$$
m(r+h)
= {{(r-h)^{2}}\over{r^{2}+(r-h)^{2}}}4 m(r)
+{{r^{2}}\over{r^{2}+(r-h)^{2}}}\left[1-3{(r-h)^{2}\over{r^{2}}}\right]m(r-h)
+ {h^{2}c_{0}\bigl((r-h)^{2}\sigma(r-h) + r^{2}\sigma(r)\bigr)\over{r^{2}+(r-h)^{2}}}.
$$
Or, writing $r_{n}=nh$, and $m_{n}=m(r_{n})$, we get
$$
m_{n+1}
= {{r_{n-1}^{2}}\over{r_{n}^{2}+r_{n-1}^{2}}}4 m_{n}
+{{r_{n}^{2}}\over{r_{n}^{2}+r_{n-1}^{2}}}\left[1-3{r_{n-1}^{2}\over{r_{n}^{2}}}\right]m_{n-1}
+ {h^{2}c_{0}\bigl(r_{n-1}^{2}\sigma(r_{n-1}) + r_{n}^{2}\sigma(r_{n})\bigr)\over{r_{n}^{2}+r_{n-1}^{2}}}.
$$

\rmk
When we modify the scalar field to include some nonzero $m_{\phi}$, we
will either (a) have to modify this routine, or (b) modify the
|source()| function to incorporate the changes.

@c void Solver::iterate(index j) {
   if(j<2) return;
   return rungeKuttaIterate(j);
}

void Solver::rungeKuttaIterate(index j) {
   real h = dx();
   real r = (j-1)*h;
   real rSq = SQ(r);
   real rPrime = (j-2)*h;
   real rPrimeSq = SQ(rPrime);
   real u = rPrime/r; // 1.0 - (2.0/j) + SQ(1.0/j);
   real c = 1.0*SQ(model->coupling())/(PI*PI*rSq);

   real massTerms = 2.0*u*m_mass[j-1] + (1.0 - 2.0*u)*m_mass[j-2]; // centered difference

   // Simpson's rule for quadrature
   real k[3];
   k[0] = h*rPrimeSq*(model->source(m_mass[j-2], rPrime));
   real m = 0.5*(m_mass[j-2] + m_mass[j-1]);
   u = 0.5*(r + rPrime);
   k[1] = 4.0*h*SQ(u)*(model->source(m, u));
   k[2] = h*rSq*(model->source(m_mass[j-1], r));

   m_mass[j] = massTerms + h*(k[0]+k[1]+k[2])*c/6.0;
}

void Solver::trickyIterate(index j) {
   real h = dx();
   real r = (j-1)*h;
   real rSq = SQ(r);
   real rPrime = (j-2)*h;
   real rPrimeSq = SQ(rPrime);
   real u =rSq + rPrimeSq;
   real g = model->coupling();
   
   real massTerms = (rPrimeSq/u)*4.0*m_mass[j-1];
   massTerms += (rSq/u)*(1.0-3.0*SQ(rPrime/r))*m_mass[j-2];

   real sourceTerms = h*rSq*(model->source(m_mass[j-1], r));
   
   real c = 2.0*SQ(g)/PI;   
   m_mass[j] = massTerms + (c*h*sourceTerms/u);
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
We now iteratively determine the solution.

@c
void Solver::run() {
     real massLowerBound = 0.0;
     real massUpperBound = model->fermionMass();
     real m, res;
     for(index j=0; j<30+(int)(log2(model->fermionMass())+0.5); j++) {
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

@ @<Set the initial condition@>=
  m = 0.5*(massUpperBound + massLowerBound);
  m_mass[0] = m;
  m_mass[1] = m;

@ We then consider the residual of the surface boundary conditions,
which indicates ``how far off'' we are.

@c
real Solver::residual() {
  real residual;
  real coef = (model->coupling()/SQ(PI));
  real h = dx();
  real u = SQ(h);
  real dr;
  real field[3];
  residual+=SQ(0.5*(m_mass[2]-4*m_mass[1]+3*m_mass[0])/u);
  for(index j=1; j<length-2; j++) {
    field[0] = model->massToField(m_mass[j-1]);
    field[1] = model->massToField(m_mass[j]);
    field[2] = model->massToField(m_mass[j+1]);
    real laplacian = (field[0]-2.0*field[1]+field[2])/u;
    dr = laplacian - coef*CUBE(m_mass[j])*Util::i((model->fermiMomentum(j*h))/fabs(m_mass[j]));
    if (isnan(dr)) {
      std::cout<<j<<std::endl;
    } else {
      residual += SQ(dr);
    }
  }
  field[0] = model->massToField(m_mass[length-3]);
  field[1] = model->massToField(m_mass[length-2]);
  field[2] = model->massToField(m_mass[length-1]);      
  real derivativeOnSurface = (0.5*(field[0]-4.0*field[1]+3.0*field[2])/h);
  dr = derivativeOnSurface + (field[2]/(model->nuggetSize()));
  std::cout<<"Surface boundary condition gives: "<<dr<<std::endl;
  return dr;
  residual +=SQ(dr);
  residual = sqrt(residual);
  std::cout<<"Residual: "<<residual<<std::endl;
  return residual;
}  
