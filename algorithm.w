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
We now iteratively determine the solution.

@c
void Solver::run() {
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

@ @<Set the initial condition@>=
  m = 0.5*(massUpperBound + massLowerBound);
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
