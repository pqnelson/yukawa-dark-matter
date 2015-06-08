@** User Interface.
I've decided to opt for a simple command line interface which will ask
the user for various options. This is nothing fancy, and really will
just be responsible for determining which routines to run.

@c
void plotBindingEnergy(YukawaDarkMatter* oldParameters);
void plotNuggetSize(YukawaDarkMatter* oldParameters);
void plotEffectiveMass(YukawaDarkMatter* oldParameters);
void determineNuggetSize(YukawaDarkMatter* oldParameters);
void plotNuggetSizeAgainstScalarMass(YukawaDarkMatter* oldParameters);
void updateVerbosity();
enum UserInput {
  USER_UNKNOWN = 0,@/
  USER_DETERMINE_NUGGET_SIZE = 1,@/
  USER_ENERGY_VS_FERMION_NUMBER = 2,@/
  USER_RADIUS_VS_FERMION_NUMBER = 3,@/
  USER_PLOT_EFFECTIVE_MASS = 4,@/
  USER_NUGGET_SIZE_VS_SCALAR_MASS = 5,@/
  USER_VERBOSITY = 6,@/
  USER_QUIT = 7
};

@ We begin by asking the user what they'd like to do. This is a simple
``Read-Eval-Print'' loop, hence we call the method |repl()|.
@^Read-Eval-Print@>

@c
void repl() {
  bool terminate = false;
  YukawaDarkMatter* model = nullptr;
  while(!terminate) {
    try {
      @<Print User Choices and Get Choice@>@;
      @<Evaluate User Input@>@;
    } catch (NoSolutionExistsException e) {
      std::cerr<<"Error: no such solution exists!"<<std::endl;
    } catch (...) {
      std::cerr<<"Unknown error caught!"<<std::endl;
    }
    delete model;
  }
  delete model;
}

@ @<Print User Choices and Get Choice@>=
  std::cout<<"(1) Determine the nugget size R for given inputs"<<std::endl;
  std::cout<<"(2) Plot Fermion Number against Binding Energy"<<std::endl;
  std::cout<<"(3) Plot Nugget size (R) against Binding Energy"<<std::endl;
  std::cout<<"(4) Plot the scale-invariant mass given the "
           <<"fermion mass and fermion number"<<std::endl;
  std::cout<<"(5) Plot the Nugget size (R) as a function of scalar mass"
           <<std::endl;
  std::cout<<"(6) Change warning level [Currently: "
           <<LOG::getVerbosity()<<"]"<<std::endl;
  std::cout<<"(7) Quit"<<std::endl;
  int userChoice;
  std::cin>>userChoice;

@ @<Evaluate User Input@>=
  switch(userChoice) {
    case USER_DETERMINE_NUGGET_SIZE:@#
      determineNuggetSize(model);
      break;
    case USER_ENERGY_VS_FERMION_NUMBER:@#
      plotBindingEnergy(model);
      break;
    case USER_RADIUS_VS_FERMION_NUMBER:
      break;
    case USER_PLOT_EFFECTIVE_MASS:@#
      plotEffectiveMass(model);
      break;
    case USER_NUGGET_SIZE_VS_SCALAR_MASS:@#
      plotNuggetSizeAgainstScalarMass(model);
      break;
    case USER_VERBOSITY:@#
      updateVerbosity();
      break;
    case USER_QUIT:@#
      std::cout<<"Goodbye"<<std::endl;
      terminate=true;
      break;
    default:@#
      std::cout<<"Invalid choice, "<<userChoice<<std::endl;
      break;
  }

@ Getting the verbosity is a fairly straightforward thing.

@c
void updateVerbosity() {
  std::cout<<"What level of debugging/tracing would you like?\n"
           <<"(1) Trace\n"
           <<"(2) Debug\n"
           <<"(3) Info\n"
           <<"(4) Warn\n"
           <<"(5) Error\n"
           <<"(6) Fatal\n"
           <<"(7) Quiet"<<std::endl;
  int verb;
  std::cin>>verb;
  LOG::setVerbosity(6-verb);
}

@* Getting Parameter Values.
Any of these requires the user to initialize the model with various
parameters. This will be used many times, so we make this its own
routine. Since in all cases we want to determine the nugget size, we
simply use our previous approximation (\S\initialGuessForNuggetSize)
that $R\approx 16\sqrt{N}\hbar c/m_{\chi}$.

@c
real getParam(const char* paramName, bool strictlyPositive=false) {
  real val;
  std::cout<<"Enter the "<<paramName<<": ";
  std::cin>>val;
  while((val<0.0)||(val==0.0&&strictlyPositive)) {
    std::cout<<"Invalid value for "<<paramName;
    std::cout<<", it must be ";
    if (strictlyPositive) {
      std::cout<<"strictly positive";
    } else {
      std::cout<<"non-negative";
    }
    std::cout<<". Please enter a new value: ";
    std::cin>>val;
  }
  return val;
}

YukawaDarkMatter* promptUserForParameters(YukawaDarkMatter *oldValues) {
  delete oldValues;
  real massChi, alpha, fermionNumber, scalarMass, a;
  massChi = getParam("fermion mass (m_{chi}) in GeV", true);
  std::cout<<"Enter the alpha value: "; 
  std::cin>>alpha;
  fermionNumber = getParam("fermion number (N)", true);
  scalarMass = getParam("scalar mass (m_{phi}) in GeV");
  real R = 16.0*sqrt(fermionNumber)*(hbarC/massChi);
  a = getParam("momentum parameter (a)");
  YukawaDarkMatter* model = new YukawaDarkMatter(massChi, R, convertAlphaToCoupling(alpha), fermionNumber, scalarMass, a);
  return model;
}

@ {\bf ``Hit Enter for Default Value''.}
This functionality is reused many places, so I thought I ought to
refactor it out as a function.

@c

template <typename T>
bool useDefaultValue(std::string msg, T defaultVal, T &returnValue) {
  std::cout<<msg
           <<" [Hit enter for "
           <<defaultVal
           <<"]"<<std::endl;
  std::cin.sync();
  char c;
  std::cin.get();
  std::cin.get(c);
  if (c=='\n') {
    std::cin.sync();
    returnValue = defaultVal;
    return true;
  }
  std::cin.putback(c);
  std::cin>>returnValue;
  return false;
}

@ Determining the nugget size is fairly straightforward, and will
just be printed to the screen. 

@c
index promptUserForLength() {
  index length, defaultLength = 1000000;
  useDefaultValue("How many subintervals would you like?",
                  defaultLength, length);
  return length;
}

void updateOldParameters(YukawaDarkMatter *oldParameters, YukawaDarkMatter *model) {
  std::swap(oldParameters, model);
  if (model) delete model;
}

void determineNuggetSize(YukawaDarkMatter* oldParameters) {
  YukawaDarkMatter* model = promptUserForParameters(oldParameters);
  index length = promptUserForLength();
  Solver *solver = new Solver(model, length);
  solver->findNuggetSize();
  solver->run();
  std::cout<<"Residual: "
           <<std::setprecision(20)
           <<(solver->residual())<<std::endl;
  std::cout<<"Nugget Size (R) = "
           <<std::setprecision(20)
           <<(model->nuggetSize())<<std::endl;
  delete solver;
  updateOldParameters(oldParameters, model);
}

@* Plot Effective Mass.
This is a bit of a lie, what we really do is produce data then pick out
a reasonable subcollection printed to a file which the user plots later
using Gnuplot (or something else).

@c
void plotEffectiveMass(YukawaDarkMatter *oldParameters) {
  YukawaDarkMatter* model = promptUserForParameters(oldParameters);
  index length = promptUserForLength();
  Solver *solver = new Solver(model, length);
  solver->findNuggetSize();
  solver->run();
  
  @<Determine Output File@>@;
  @<Print Solution to File@>@;
  @<Print Momentum to File@>@;
  data.close();
  delete solver;
  updateOldParameters(oldParameters, model);
}

@ @<Determine Output File@>=
  std::ofstream data;
  std::string filename = "";
  std::cout<<"What file would you like to print the data?"<<std::endl;
  std::cin>>filename;
  data.open(filename);
  data<<"# fermion mass = "<<(model->fermionMass())
      <<", alpha = "<<couplingToAlpha(model->coupling())
      <<", N = "<<(model->fermionNumber())<<"\n";

@ @<Print Solution to File@>=
  index MAX_ITER = 50;
  real dx = (model->nuggetSize())/(1.0*MAX_ITER);
  index step = (solver->getLength())/MAX_ITER;
  real *solution = solver->getMass();
  data<<"# First index block (index 0)\n";
  data<<"# r    m(r)\n";
  data<<"0.0"<<"    "<<(solution[0])<<"\n";
  for(index j=1; j<MAX_ITER; j++) {
    data<<(j*dx)<<"    "<<(solution[j*step])<<"\n";
  }
  data<<(model->nuggetSize())<<"    "<<(solution[(solver->getLength())-1])<<"\n";

@ @<Print Momentum to File@>=
  data<<"\n\n\n";
  data<<"# Second index block (index 1)\n";
  data<<"# r    p(r)\n";
  data<<"0.0"<<"    "<<(model->fermiMomentum(0.0))<<"\n";
  real r;
  for(index j=1; j<MAX_ITER; j++) {
    r = j*dx;
    data<<r<<"    "<<(model->fermiMomentum(r))<<"\n";
  }
  r = model->nuggetSize();
  data<<r<<"    "<<(model->fermiMomentum(r))<<"\n";

@* Plot binding energy against fermion number.
We will plot the binding energy as a function of the fermion number
$N$. The binding energy is defined as
$$
E_{{\rm bind}} := m_{\chi} - {E(N,R)\over N}.\eqn{}
$$
We will have the plot be logarithmic in $N$, so it will step as
$10^{k/2}$ for $k=2,\dots,k_{{\rm max}}$.

@c
real bindingEnergy(real energy, real N, real massChi) {
  return massChi - (energy/N);
}

@ @c
index promptUserForMaxFermionNumber() {
  index ret;
  std::cout<<"What is the log-base 10 of the maximum fermion number (N=10^k, k=?)?"<<std::endl;
  std::cout<<"k = ";
  std::cin>>ret;
  return ret;
}

@ @c
real getBindingEnergy(Solver* solver, real N, YukawaDarkMatter *model) {
  return bindingEnergy(solver->computeEnergy(), N, model->fermionMass());
}

@ We will end up producing a couple of log-plots, so it makes sense to
refactor the code to keep things DRY\footnote{${}^{*}$}{DRY = ``Don't
Repeat Yourself''. Repetition of code is bad form, a source of errors,
hinders robust products, causes endless suffering, etc.}

@c
typedef std::function<real(Solver*,real,YukawaDarkMatter*)> plotFn;

void logPlot(std::string label, YukawaDarkMatter* oldParameters, plotFn fn) {
  YukawaDarkMatter* model = promptUserForParameters(oldParameters);
  index length = promptUserForLength();
  index MAX_ITER = promptUserForMaxFermionNumber();
  @<Determine Output File@>@;
  Solver *solver = new Solver(model, length);
  data<<"# N    "<<label<<"\n";
  real N;
  for(index j=0; j<1+2*MAX_ITER; j++) {
    N = pow(10.0, 1.0+0.5*j);
    model->setFermionNumber(N);
    solver->findNuggetSize();
    solver->run();
    data<<N<<"    "<<(fn(solver, N, model))<<"\n";
  }
  data.close();
  delete solver;
  updateOldParameters(oldParameters, model);
}

@ @c
void plotBindingEnergy(YukawaDarkMatter* oldParameters) {
  logPlot("E", oldParameters, getBindingEnergy);
}

@* Plot Nugget Size against Fermion Number.
The last item which we will be graphing is the nugget size against the
fermion number $N$. This is remarkably similar to what we had done with
the binding energy against $N$.

@c
real getNuggetSize(Solver* solver, real N, YukawaDarkMatter *model) {
  return (model->nuggetSize());
}

void plotNuggetSize(YukawaDarkMatter* oldParameters) {
  logPlot("R", oldParameters, getNuggetSize);
}

@* Plotting Nugget Size as Function of Scalar Mass.
This actually comes in handy when trying to figure out how scalar mass
affects the nugget's stability.

@c
void plotNuggetSizeAgainstScalarMass(YukawaDarkMatter* oldParameters) {
  YukawaDarkMatter* model = promptUserForParameters(oldParameters);
  real maxScalarMass, dScalarMass;
  std::cout<<"What is the max scalar mass?"<<std::endl;
  std::cin>>maxScalarMass;
  std::cout<<"What is the step size for the scalar mass?"<<std::endl;
  std::cin>>dScalarMass;
  index MAX_ITER = (index)(maxScalarMass/dScalarMass);
  std::cout<<"Will have "<<MAX_ITER<<" iterations"<<std::endl;
  index length = promptUserForLength();
  @<Determine Output File@>@;
  Solver *solver = new Solver(model, length);
  data<<"# scalar mass    R    residual\n";
  real R, scalarMass;
  for(index j=0; j<1+MAX_ITER; j++) {
    try {
      scalarMass = j*dScalarMass;
      model->setScalarMass(scalarMass);
      solver->findNuggetSize();
      solver->run();
      R = model->nuggetSize();
      data<<scalarMass<<"    "<<R<<"    "<<(solver->residual())<<"\n";
      data.flush();
    } catch (NoSolutionExistsException e) {
      LOG::warn<<"Skipping scalar mass "
               <<scalarMass
               <<", no solution could be found"
               <<std::endl;
    }
  }
  data.close();
  delete solver;
  updateOldParameters(oldParameters, model);
}