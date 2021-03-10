/************************************/
/*PARTICLE SWARM OPTIMIZATION */
/************************************/
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;
// 
//  Change any of these parameters to match your needs 
//
# define POPSIZE 10
# define MAXGENS 1000
# define NVARS 5
# define PXOVER 0.8
# define PMUTATION 0.15
//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//
struct genotype
{
  double gene[NVARS];
  double fitness;
  double upper[NVARS];
  double lower[NVARS];
  double rfitness;
  double cfitness;
};

struct genotype population[POPSIZE+1];
struct genotype newpopulation[POPSIZE+1]; 

double target;
double zeroone(double x)
{
  if (x<=0.5){return 0.;}
  return 1;
}
int main ( int argc, char *argv[]);
void crossover ( int &seed );
void elitist ( );
void evaluate ( );
int i4_uniform_ab ( int a, int b, int &seed );
void initialize ( string filename, int &seed );
void keep_the_best ( );
void mutate ( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
void report ( int generation );
void selector ( int &seed );
void timestamp ( );
void Xover ( int one, int two, int &seed );

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <iostream>
/*CONSTANTS */
#define NOPS 10    /*NUMBER OF PARTICLES*/
#define LIMITL 256 /*MAXIMUM SIZE OF THE ARRAY FOR PARTICLE DATA*/
#define ILIMIT 150 /*MAXIMUM NUMBER OF ITERATION*/
#define SEED 32767 /*INITIALIZATION OF RANDOM NUMBER*/
#define W 0.7      /*INERTIA CONSTANT*/
#define C1 1.4     /*CONSTANT OF ATTRACTION FROM PERSONAL BEST*/
#define C2 1.4     /*CONSTANT OF ATTRACTION FROM GROUP BEST*/

#define DMAX 1 /*STRUCTURE FOR ARBITRARY COORDINATE*/
struct point{
  double x[DMAX];
};

typedef std::complex<double> Complex;

extern "C"
{
    void zheev_ ( const char& JOBZ, const char& UPLO,
		const int& N, Complex* A, const int& LDA,
		double* WW, Complex* WORK, const int& LWORK,
		double* RWORK,
		int& INFO, int JOBZlen, int UPLOlen );
  void zgetrf_(int *M, int *N, Complex *A, int *LDA, int *IPIV, int *INFO);
  void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
  void zhegv_(const int &ITYPE, const char &JOBZ, const char &UPLO,
              const int &N,
              Complex **A, const int &LDA,
              Complex **B, const int &LDB,
              double *WW, Complex *WORK, const int &LWORK,
              double *RWORK,
              int &INFO, int JOBZlen, int UPLOlen);
};

using namespace std;

/*STRUCTURE FOR PARTILE*/
struct particle
{
  struct point pos;     /*POSITION*/
  double value;         /*EVALUATED VALUE*/
  struct point v;       /*VELOCITY*/
  struct point bestpos; /*THE BEST OF EVER-REACHED POSITIONS*/
  double bestval;       /*THE BEST OF EVER-EVALUATED VALUES*/
};

/*PROTOTYPES OF FUNCTIONS */
void initps(struct particle ps[]);   /*INITIALIZATION OF PARTICLE SWARM*/
void printps(struct particle ps[]);  /*PRINT OUT OF PARTICLE SWARM DATA*/
double frand();                      /*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE*/
double calcval(double *x);           /*THE OBJECTIVE FUNCTION*/
double calcresidue(double *x);       /*THE RESIDUE OF THE OPTIMIZATION*/
void optimize(struct particle ps[]); /*OPTIMIZATION*/
void setgbest(struct particle ps[]); /*FIND THE GROUP BEST*/

/*GLOBAL VARIABLES*/
struct point gbestpos;       /*THE BEST POSITION IN THE GROUP*/
double gbestval;             /*THE BEST VALUE IN THE GROUP*/

/**********************************************/
/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE */
/**********************************************/
double doublerand(void)
{
  double result;
  while ((result = (double)rand() / RAND_MAX) >= 1)
    ;
  return result;
}

/*VARABLES TO CONSTRUCT A QUANTUM GRAPH*/

double lambda[10];
double L[10];
double A[10];
int B = 5;
double k;

double h00(double k)
{
  double cc = 0.;
  cc += lambda[0] / k;
  for (int i = 1; i <= B; i++)
  {
    cc += 1. / tan(k * L[i]);
  }
  return cc;
}

double hjj(double k, int j)
{
  double cc = 0.0;
  cc += lambda[j] / k;
  cc += 1. / tan(k * L[j]);
  return cc;
}

complex<double> h0j(double k, int j)
{
  complex<double>
   cc;
  cc = complex<double>(cos(A[j] * L[j]), -sin(A[j] * L[j]));
  cc /= sin(k * L[j]);
  cc *= -1;
  return cc;
}

double zetah(double k)
{
  complex<double> A[(B + 1) * (B + 1)];
  A[0] = h00(k);
  for (int i = 1; i <= B; i++)
  {
    A[0 + (B + 1) * i] = h0j(k, i);
    A[i + (B + 1) * 0] = conj(A[0 + (B + 1) * i]);
    A[i + (B + 1) * i] = hjj(k, i);
  }

  int m = B + 1;   // 行のサイズ
  int n = B + 1;   // 列のサイズ
  int lda = B + 1; // mと同じ値

  int info;                    // 計算が成功すれば0を返す
  int ipiv[B + 1];             // 要素数はm,nのうち小さい方とする
  int lwork = B + 1;           // nと同じ値
  complex<double> work[B + 1]; // 要素数はlworkと同じ値
  //printf("\n%+8.5lf %+8.5lf\n", real(A[0]), real(A[2]));
  //printf("%+8.5lf %+8.5lf\n", real(A[1]), real(A[3]));
  //printf("%+8.5lf %+8.5lf\n", imag(A[0]), imag(A[2]));
  //printf("%+8.5lf %+8.5lf\n", imag(A[1]), imag(A[3]));

  // LAPACKのdgetrfサブルーチンを呼んで、行列AをLU分解
  // 引数は全て参照渡し
  zgetrf_(&m, &n, A, &lda, ipiv, &info);

  // LU分解後の行列から逆行列を求める
  // 逆行列は元の配列Aに入る
  //dgetri_( &n, A, &lda, ipiv, work, &lwork, &info);

  //  printf("\n%+8.5lf %+8.5lf\n", real(A[0]), real(A[2]));
  //  printf("%+8.5lf %+8.5lf\n", real(A[1]), real(A[3]));
  //  printf("%+8.5lf %+8.5lf\n", imag(A[0]), imag(A[2]));
  //  printf("%+8.5lf %+8.5lf\n", imag(A[1]), imag(A[3]));
  double det = 1;
  for (int i = 0; i <= B; i++)
  {
    int j = i + (B + 1) * i;
    //printf("%+8.5lf %+8.5lf %d\n", real(A[j]), imag(A[j]),ipiv[i]);
    det *= real(A[j]);
    if (ipiv[i] != i + 1)
      det = -det;
  }
  return det;
  //out<<"det="<<det<<endl;
}

double zetak(double k, double *e)
{
  double cc = h00(k);

  for (int i = 1; i <= B; i++)
  {
//    cout<<L[i];
    double absh0j = abs(h0j(k, i));
    cc -= e[i-1]*absh0j * absh0j / hjj(k, i);
  }
//  cout<<endl;
  for (int i = 1; i <= B; i++)
  {
    cc *= hjj(k, i);
  }
  return cc;
}



double calceigen(double *XIN)
{
  complex<double> U[(B+1)*(B+1)];
  int N=B+1;
  double E[B+1];
  int i, j, info;
  complex<double> cwork[4*(B+1)];
  double rwork[4*(B+1)];
  U[0]=complex<double>(h00(XIN[0]),0.0);
  for (i=1;i<=B;i++)
  {
    
    U[0+(B+1)*i]=h0j(XIN[0],i);
    U[i+(B+1)*0]=conj(h0j(XIN[0],i));
    U[i+(B+1)*i]=Complex(hjj(XIN[0],i));
    
//    U[i+(B+1)*i]=i;
//    U[i+(B+1)*0]=conj(h0j(XIN[0],i));
//    U[i+(B+1)*i]=Complex(hjj(XIN[0],i));

  }
  for(int j=0;j<=B;j++){
    for(int i=0;i<=B;i++){
      cout<<U[i+j*(B+1)]<<" ";
    }
    cout<<endl;
  }
  zheev_( 'V', 'U', N, (Complex*)U, N, E, cwork, 4*N, rwork, info, 1, 1 );
  double det=1;
  for(i=0;i<=B;i++){
    det*=E[i];
    printf("%lf :",E[i]);
    for(j=0;j<=B;j++){
    cout<<U[j+i*(B+1)]<< " ";
  }
  cout<<endl;
  }

  printf( " det*det=%lf\n",det*det);
  return  0;
}



/*********************************************/
/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE*/
/*********************************************/
double frand(void)
{
  double result;
  while ((result = (double)rand() / RAND_MAX) >= 1)
    ;
  return result;
}

double findzerop(double *eg)
{
//  cout<<eg[0]<<eg[1]<<eg[2]<<eg[3]<<eg[4]<<endl;
  double v=zetak(0.04,eg);
  double ka=0.04;
  for (double k = 0.05; k <= 3; k += 0.01)
  {
    double vp=zetak(k,eg);
    if (v*vp<0.){
//      cout<<v<<" "<<vp<<" "<<k<<endl;
      return (k+ka)/2.;
    }
    ka=k;
    v=vp;
  }
  return 0;
}

//****************************************************************************80

int main (int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN supervises the genetic algorithm.
//
//  Discussion:
//
//    Each generation involves selecting the best 
//    members, performing crossover & mutation and then 
//    evaluating the resulting population, until the terminating 
//    condition is satisfied   
//
//    This is a simple genetic algorithm implementation where the 
//    evaluation function takes positive values only and the 
//    fitness of an individual is the same as the value of the 
//    objective function.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    MAXGENS is the maximum number of generations.
//
//    NVARS is the number of problem variables.
//
//    PMUTATION is the probability of mutation.
//
//    POPSIZE is the population size. 
//
//    PXOVER is the probability of crossover.                          
//
{
  if (argc!=0)
  {
    target=atof(argv[1]);
  }
  complex<double> a, b, c;
  a = complex<double>(3, 2);
  cout << a << endl;
  double eg[11];
  int seed0=123456798;
  for (int i = 0; i <= 10; i++)
  {
    L[i] = 1;
    L[i]= r8_uniform_ab(0.5,1.5,seed0);
    A[i] = 0;
    lambda[i] = r8_uniform_ab(0.,1.,seed0);
    eg[i]=1;
  }

  eg[1]=0;
  eg[2]=0;
  eg[3]=0;
  eg[4]=0;
  for (double k = 0.05; k <= 3; k += 0.01)
  {
  //  cout <<"K "<< k << " " << zetak(k,eg) << endl;
  }
  cout<<"zero="<<findzerop(eg)<<endl;
  //exit(1);

  string filename = "simple_ga_input.txt";
  int generation;
  int i;
  int seed;

  timestamp ( );
  cout << "\n";
  cout << "SIMPLE_GA:\n";
  cout << "  C++ version\n";
  cout << "  A simple example of a genetic algorithm.\n";

  if ( NVARS < 2 )
  {
    cout << "\n";
    cout << "  The crossover modification will not be available,\n";
    cout << "  since it requires 2 <= NVARS.\n";
  }

  seed = 123456789;

  initialize ( filename, seed );

  evaluate ( );

  keep_the_best ( );

  for ( generation = 0; generation < MAXGENS; generation++ )
  {
    selector ( seed );
    crossover ( seed );
    mutate ( seed );
    report ( generation );
    evaluate ( );
    elitist ( );
  }

  cout << "\n";
  cout << "  Best member after " << MAXGENS << " generations:\n";
  cout << "\n";

  for ( i = 0; i < NVARS; i++ )
  {
    cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
  }

  cout << "\n";
  cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "SIMPLE_GA:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";

  cout <<"TARGET:"<<target<<" BEST GENE: (";
  for ( i = 0; i < NVARS; i++ )
  {
    cout<<" " << population[POPSIZE].gene[i] << ",";
  }
  cout<<") BEST VALUE:"<<findzerop(&population[POPSIZE].gene[0]);
  cout<<endl;
  timestamp ( );

  return 0;
}
//****************************************************************************80

void crossover ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int FIRST, is a count of the number of members chosen.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int mem;
  int one;
  int first = 0;
  double x;

  for ( mem = 0; mem < POPSIZE; ++mem )
  {
    x = r8_uniform_ab ( a, b, seed );

    if ( x < PXOVER )
    {
      ++first;

      if ( first % 2 == 0 )
      {
        Xover ( one, mem, seed );
      }
      else
      {
        one = mem;
      }

    }
  }
  return;
}
//****************************************************************************80

void elitist ( )

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;

  best = population[0].fitness;
  worst = population[0].fitness;

  for ( i = 0; i < POPSIZE - 1; ++i )
  {
    if ( population[i+1].fitness < population[i].fitness )
    {

      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }

    }

  }
// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//
  if ( population[POPSIZE].fitness <= best )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[POPSIZE].gene[i] = population[best_mem].gene[i];
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
  }
  else
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[worst_mem].gene[i] = population[POPSIZE].gene[i];
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  } 

  return;
}
//****************************************************************************80

double calcval(double x,double y) 
{ 
double twopi=8.*atan(1.);
double A= exp(-0.2*pow(0.5*(x*x+y*y),0.5));
double B= exp(0.5*(cos(twopi*x)+cos(twopi*y)));

double x1=x-20;
double y1=y-20;
double A1= exp(-0.2*pow(0.5*(x1*x1+y1*y1),0.5));
double B1= exp(0.5*(cos(twopi*x1)+cos(twopi*y1)));
return 20*A + B;
return -20*A-B+20+exp(1.);
} 

double calcval0(double x,double y) 
{ 
double resultat=0;

FILE *f1;
/*
FILE *f1=fopen("objv.in.txt","w");
fprintf(f1,"%lf %lf\n",x,y);
fclose(f1);
*/

char command[4000];
sprintf(command,"sed \"s/\\(^theta_layer1=\\).*/\\1%12.8e/\" mkinpfile_tm1da.sh.INIT | sed \"s/\\(^theta_layer2=\\).*/\\1%12.8e/\"  > mkinpfile_tm1da.sh",x,y);
system(command);
//printf("%s\n",command);
//printf("start run.sh\n");
system("sh run.sh >/dev/null");	
f1=fopen("ECHOSIMIL.txt","r");
fscanf(f1,"%lf",&resultat);
fclose(f1);
printf("resultat=%12.8lf\n",resultat);
return -resultat ;

double twopi=8.*atan(1.);
double A= exp(-0.2*pow(0.5*(x*x+y*y),0.5));
double B= exp(0.5*(cos(twopi*x)+cos(twopi*y)));

double x1=x+10;
double y1=y+10;
double A1= exp(-0.2*pow(0.5*(x1*x1+y1*y1),0.5));
double B1= exp(0.5*(cos(twopi*x1)+cos(twopi*y1)));

return x1*x1+y1*y1;
return -20*A1-B1+20+exp(1.);
} 


double fitfunc(double x)
{
//  cout<<target<<endl;
   double y=x-target;
//   return x;
   return exp(-y*y);
}
void evaluate ( )

//****************************************************************************80
// 
//  Purpose:
//
//    EVALUATE implements the user-defined valuation function
//
//  Discussion:
//
//    Each time this is changed, the code has to be recompiled.
//    The current function is:  x[1]^2-x[1]*x[2]+x[3]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
{
  int member;
  int i;
  double x[NVARS+1];

  for ( member = 0; member < POPSIZE; member++ )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      x[i+1] = population[member].gene[i];
    } 
//    population[member].fitness = 1-( x[1] * x[1] ) - ( x[2] * x[2] ) ;
    population[member].fitness = fitfunc(findzerop(&x[1]));
//     population[member].fitness = +( x[1] ) + ( x[2]  ) + x[3];
//    cout<<x[1]<<" "<<x[2]<<" "<<calcval(x[1],x[2])<<endl;
//    population[member].fitness = calcval( x[1],x[2]) ;

  }
  return;
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void initialize ( string filename, int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    INITIALIZE initializes the genes within the variables bounds. 
//
//  Discussion:
//
//    It also initializes (to zero) all fitness values for each
//    member of the population. It reads upper and lower bounds 
//    of each variable from the input file `gadata.txt'. It 
//    randomly generates values between these bounds for each 
//    gene of each genotype in the population. The format of 
//    the input file `gadata.txt' is 
//
//      var1_lower_bound var1_upper bound
//      var2_lower_bound var2_upper bound ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int i;
  ifstream input;
  int j;
  double lbound;
  double ubound;

  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "INITIALIZE - Fatal error!\n";
    cerr << "  Cannot open the input file!\n";
    exit ( 1 );
  }
// 
//  Initialize variables within the bounds 
//
  for ( i = 0; i < NVARS; i++ )
  {
    input >> lbound >> ubound;

    for ( j = 0; j < POPSIZE; j++ )
    {
      population[j].fitness = 0;
      population[j].rfitness = 0;
      population[j].cfitness = 0;
      population[j].lower[i] = lbound;
      population[j].upper[i]= ubound;
      population[j].gene[i] = zeroone( r8_uniform_ab ( lbound, ubound, seed ) );
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void keep_the_best ( )

//****************************************************************************80
// 
//  Purpose:
//
//    KEEP_THE_BEST keeps track of the best member of the population. 
//
//  Discussion:
//
//    Note that the last entry in the array Population holds a 
//    copy of the best individual.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int CUR_BEST, the index of the best individual.
//
{
  int cur_best;
  int mem;
  int i;

  cur_best = 0;

  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    if ( population[POPSIZE].fitness < population[mem].fitness )
    {
      cur_best = mem;
      population[POPSIZE].fitness = population[mem].fitness;
    }
  }
// 
//  Once the best member in the population is found, copy the genes.
//
  for ( i = 0; i < NVARS; i++ )
  {
    population[POPSIZE].gene[i] = population[cur_best].gene[i];
  }

  return;
}
//****************************************************************************80

void mutate ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    MUTATE performs a random uniform mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random value 
//    between the lower and upper bounds of this variable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;

  for ( i = 0; i < POPSIZE; i++ )
  {
    for ( j = 0; j < NVARS; j++ )
    {
      x = r8_uniform_ab ( a, b, seed );
      if ( x < PMUTATION )
      {
        lbound = population[i].lower[j];
        ubound = population[i].upper[j];  
        population[i].gene[j] = zeroone(r8_uniform_ab ( lbound, ubound, seed ));
      }
    }
  }

  return;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

void report ( int generation )

//****************************************************************************80
// 
//  Purpose:
//
//    REPORT reports progress of the simulation. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double avg, the average population fitness.
//
//    Local, best_val, the best population fitness.
//
//    Local, double square_sum, square of sum for std calc.
//
//    Local, double stddev, standard deviation of population fitness.
//
//    Local, double sum, the total population fitness.
//
//    Local, double sum_square, sum of squares for std calc.
//
{
  double avg;
  double best_val;
  int i;
  double square_sum;
  double stddev;
  double sum;
  double sum_square;

  if ( generation == 0 )
  {
    cout << "\n";
    cout << "  Generation       Best            Average       Standard \n";
    cout << "  number           value           fitness       deviation \n";
    cout << "\n";
  }

  sum = 0.0;
  sum_square = 0.0;

  for ( i = 0; i < POPSIZE; i++ )
  {
    sum = sum + population[i].fitness;
    sum_square = sum_square + population[i].fitness * population[i].fitness;
  }

  avg = sum / ( double ) POPSIZE;
  square_sum = avg * avg * POPSIZE;
  stddev = sqrt ( ( sum_square - square_sum ) / ( POPSIZE - 1 ) );
  best_val = population[POPSIZE].fitness;
  cout << "  " << setw(8) << generation 
       << "  " << setw(14) << best_val 
       << "  " << setw(14) << avg 
       << "  " << setw(14) << stddev << "\n";

  return;
}
//****************************************************************************80

void selector ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    SELECTOR is the selection function.
//
//  Discussion:
//
//    Standard proportional selection for maximization problems incorporating 
//    the elitist model.  This makes sure that the best member always survives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;
//
//  Find the total fitness of the population.
//
  sum = 0.0;
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    population[mem].rfitness = population[mem].fitness / sum;
  }
// 
//  Calculate the cumulative fitness.
//
  population[0].cfitness = population[0].rfitness;
  for ( mem = 1; mem < POPSIZE; mem++ )
  {
    population[mem].cfitness = population[mem-1].cfitness +       
      population[mem].rfitness;
  }
// 
//  Select survivors using cumulative fitness. 
//
  for ( i = 0; i < POPSIZE; i++ )
  { 
    p = r8_uniform_ab ( a, b, seed );
    if ( p < population[0].cfitness )
    {
      newpopulation[i] = population[0];      
    }
    else
    {
      for ( j = 0; j < POPSIZE; j++ )
      { 
        if ( population[j].cfitness <= p && p < population[j+1].cfitness )
        {
          newpopulation[i] = population[j+1];
        }
      }
    }
  }
// 
//  Overwrite the old population with the new one.
//
  for ( i = 0; i < POPSIZE; i++ )
  {
    population[i] = newpopulation[i]; 
  }

  return;     
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void Xover ( int one, int two, int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    XOVER performs crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int i;
  int point;
  double t;
// 
//  Select the crossover point.
//
  point = i4_uniform_ab ( 0, NVARS - 1, seed );
//
//  Swap genes in positions 0 through POINT-1.
//
  for ( i = 0; i < point; i++ )
  {
    t                       = population[one].gene[i];
    population[one].gene[i] = population[two].gene[i];
    population[two].gene[i] = t;
  }

  return;
}


