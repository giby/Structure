/*Library of random number generating functions.

VERSION 9-25-99

Includes:
Randomize()              (seed random number generator)
RandomReal(low,high)     (uniform)
RandomInteger(low,high)  (uniform integer in [low,high])
rnd()                    (uniform in (0,1))
RGamma(n,lambda)         (gamma)
RDirichlet(a[],k,b[])    (dirichlet)
RPoisson(mean)           (poisson)
RNormal(mean,sd)        (Normal)
RExpon(mean)             (Exponential)
snorm()                  (Standard normal)
Binomial(n, p)           (Binomial rv)


extern void Randomize(void);  
extern double RandomReal(double low, double high);
extern int RandomInteger(int low, int high);
extern double rnd();
extern double RGamma(double n,double lambda);
extern void RDirichlet(const double * a, const int k, double * b);
extern long RPoisson(double mu);
extern double RExpon(double av);
extern double RNormal(double mu,double sd) ;
extern double fsign( double num, double sign );
extern double sexpo(void);
extern double snorm();
extern double genexp(double av);   
extern long ignpoi(double mean);  
extern long ignuin(int low, int high);   
extern double genunf(double low, double high);   
extern long Binomial(int n, double p)

MORE DETAILS BELOW:




 Random number functions from random.c by Eric Roberts 
 void Randomize(void);    
                  Seed the random number generator 
 double RandomReal(double low, double high);
                  Get a random number between low and high 
 int RandomInteger(int low, int high);
                  Get a random integer between low and high INCLUSIVE
 double rnd();
                  Uniform(0,1) random number generation


 Random number functions from Matthew Stephens 
   
 double RGamma(double n,double lambda);
                 gamma random generator from Ripley, 1987, P230 
 void RDirichlet(const double * a, const int k, double * b);
                 Dirichlet random generator
    a and b are arrays of length k, containing doubles.
    a is the array of parameters
    b is the output array, where b ~ Dirichlet(a)  


Functions from Brown and Lovato

 long RPoisson(double mu);
                  Poisson with parameter mu
 double RExpon(double av);
                  exponential with parameter av
 double RNormal(double mu,double sigsq) ;
                  Normal with mean mu, var sigsq; by JKP
  ---------Helper functions from Brown and Lovato
 double fsign( double num, double sign );
 double sexpo(void);
                  exponential with parameter 1
 double snorm(); 
                  standard normal N(0,1)  


 double genexp(double av);   return RExpon(av); 
 long ignpoi(double mean);   return RPoisson(mean); 
 long ignuin(int low, int high);    return RandomInteger(low,high);
 double genunf(double low, double high);   return RandomReal(low,high);
 */


#include "ran.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define OVERFLO 1e100
#define UNDERFLO 1e-100
/*#define OVERFLOW   1e100 */ /*take log of likelihood when product is more than this*/
/*#define UNDERFLOW  1e-100 */

/*==============================================*/
/*==============================================*/
/*==============================================*/

/* Random number functions (from random.c by Eric Roberts) */

/*Melissa modified in 1/08 so that it either calls srand with given seed or generates one*/
void Randomize(int RANDOMIZE, int *seed)
/* Seed the random number generator */
{   
  FILE *outfile;
  if (RANDOMIZE) 
    *seed = (int)time(NULL);
  srand(*seed);
  outfile = fopen("seed.txt", "a");
  fprintf(outfile, "%i\n", *seed);
  fclose(outfile);
  /*  srand((int) time(NULL) );  */
}
/*-------------------------------------*/
double RandomReal(double low, double high)
/* Get a random number between low and high */

{
  double d;
 
  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low) );
}
/*-------------------------------------*/
int RandomInteger(int low, int high)
/* Get a random integer between low and high INCLUSIVE*/
{
  int k;
  double d;
 
  d = (double) rand() / ((double) RAND_MAX + 1);
  k = (int) (d * (high - low + 1));
  return (low + k);
}

/*=======================================================*/
/*  Uniform(0,1) random number generation*/

double rnd()
{
  double value;

  do
    value = RandomReal(0.0,1.0);
  while ((value==0.0)||(value==1.0));

  return value;
}
/*-----------Gamma and dirichlet from Matt.----------*/
/* gamma random generator from Ripley, 1987, P230 */


double RGamma(double n,double lambda)
{
  double aa;
  double w;
  /*  int i; */

	double x=0.0;
	if(n<1)
	{
		const double E=2.71828182;
		const double b=(n+E)/E;
		double p=0.0;
		one: 
		p=b*rnd();
		if(p>1) goto two;
		x=exp(log(p)/n);
		if(x>-log(rnd())) goto one;
		goto three;
		two: 
		x=-log((b-p)/n);
		if (((n-1)*log(x))<log(rnd())) goto one;
		three:;	
	}
	else if(n==1.0)

	  /* exponential random variable, from Ripley, 1987, P230  */	
	{
		double a=0.0;
		double u,u0,ustar;
	ten:
		u=rnd();
		u0=u;
	twenty:
		ustar=rnd();
		if(u<ustar) goto thirty;
		u=rnd();
		if(u<ustar) goto twenty;
		a++;
		goto ten;
	thirty:
		return (a+u0)/lambda;
	}
	else
	{
		double static nprev=0.0;
		double static c1=0.0;
		double static c2=0.0;
		double static c3=0.0;
		double static c4=0.0;
		double static c5=0.0;
		double u1;
		double u2;
		if(n!=nprev)
		{
			c1=n-1.0;
			aa=1.0/c1;
			c2=aa*(n-1/(6*n));
			c3=2*aa;
			c4=c3+2;
			if(n>2.5) c5=1/sqrt(n);
		}
		four:
		u1=rnd();
		u2=rnd();
		if(n<=2.5) goto five;
		u1=u2+c5*(1-1.86*u1);
		if ((u1<=0) || (u1>=1)) goto four;
		five:
		w=c2*u2/u1;
		if(c3*u1+w+1.0/w < c4) goto six;
		if(c3*log(u1)-log(w)+w >=1) goto four;
		six:
		x=c1*w;		
		nprev=n;
	}	

	return x/lambda;
}


/*
double
LogRGamma (double n, double lambda)
{
  //double aa;
  //  double w;
  //  int i;
  double logx;
  //  return log(RGamma(n, lambda));
  if (n < 1)
  //this is the case we need to worry about underflow
  //copied code from down below but work with logx
  //instead of x
    {
      const double E = 2.71828182;
      const double b = (n + E) / E;
      double p = 0.0;
    one:
      p = b * rnd ();
      if (p > 1)
        goto two;
      logx =  log (p) / n;
      if (logx > log(-log (rnd ())))
        goto one;
      goto three;
    two:
      logx = log(-log (b - p)) -log(n);

      if (((n - 1) * logx) < log (rnd ()))
        goto one;
    three:
return logx-log(lambda);
}
else
//otherwise log the standard version 
return log(RGamma(n,lambda));
}*/




/* Melissa's version, adapted from an algorithm on wikipedia.  January 08 */
double LogRGamma(double n, double lambda)
{
  double v0, v[3], E=2.71828182, em, logem, lognm;
  int i;
  if (lambda!=1.0) {
    printf("lambda=%e!\n", lambda); exit(-1);
  }
  if (n >= 1.0) {
    return log(RGamma(n, lambda));
  }
  v0 = E/(E+n);
  while (1) {
    for (i=0; i<3; i++) {
      v[i] = rnd();
    }
    
    if (v[0] <= v0) {
      logem = 1.0/n*log(v[1]);
      em = exp(logem);
      lognm = log(v[2])+(n-1)*logem;
    } else {
      em = 1.0-log(v[1]);
      logem = log(em);
      lognm = log(v[2]) - em;
    }
    if (lognm <= (n-1)*logem - em) {
      return logem - log(lambda);
    }
  }
}



/*--------------------------------------*/

/* Dirichlet random generator
   a and b are arrays of length k, containing doubles.
   a is the array of parameters
   b is the output array, where b ~ Dirichlet(a)  
   */

void RDirichlet(const double * a, const int k, double * b)
{
int i;
	double sum=0.0;
	for(i=0;i<k;i++)
	{
		b[i]=RGamma(a[i],1);
		sum += b[i];
	}
	for(i=0;i<k;i++)
	{
		b[i] /= sum;
	}
}


/*This function returns both a logged and unlogged version
of the dirichlet function. Designed to cope with
underflows in the RGamma function.
made by Daniel
b is the output array and c is a logged version of b*/

void
LogRDirichlet (const double *a, const int k, double *b,double *c)
{
  int i;
  double sum = 0.0;
  double sum2;
  for (i = 0; i < k; i++) {
    c[i] = LogRGamma (a[i], 1);
    b[i]=exp(c[i]);
    sum += b[i];
  }
  
  /* patch added May 2007 to set gene frequencies equal if all draws from the Gamma distribution are very low. Ensures that P and logP remain defined in this rare event */
  if(sum<UNDERFLO) {
    for(i=0;i<k;i++) {
      b[i] = 1.0/(double)(k);
      c[i] = log(b[i]);
    }
  } else {
    sum2=log(sum);
    for (i = 0; i < k; i++) {
      c[i]-=sum2;
      b[i]/=sum;
    }
  }
}


/*---------------------------------------*/

long RPoisson(double mu)
/*
**********************************************************************
     long RPoissondouble mu)
                    GENerate POIsson random deviate
                              Function
     Generates a single random deviate from a Poisson
     distribution with mean AV.
                              Arguments
     av --> The mean of the Poisson distribution from which
            a random deviate is to be generated.
     RExpon <-- The random deviate.
                              Method
     Renames KPOIS from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.

     ----substituted rnd for ranf--JKP, 11/98------

     For details see:
               Ahrens, J.H. and Dieter, U.
               Computer Generation of Poisson Deviates
               From Modified Normal Distributions.
               ACM Trans. Math. Software, 8, 2
               (June 1982),163-179
**********************************************************************
**********************************************************************
                                                                      
                                                                      
     P O I S S O N  DISTRIBUTION                                      
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER GENERATION OF POISSON DEVIATES                
               FROM MODIFIED NORMAL DISTRIBUTIONS.                    
               ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. 
                                                                      
     (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  
                                                                      
**********************************************************************
      INTEGER FUNCTION RPOISSONIR,MU)
     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
             MU=MEAN MU OF THE POISSON DISTRIBUTION
     OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
     SEPARATION OF CASES A AND B
*/
{
extern double fsign( double num, double sign );
static double a0 = -0.5;
static double a1 = 0.3333333;
static double a2 = -0.2500068;
static double a3 = 0.2000118;
static double a4 = -0.1661269;
static double a5 = 0.1421878;
static double a6 = -0.1384794;
static double a7 = 0.125006;
static double muold = 0.0;
static double muprev = 0.0;
static double fact[10] = {
    1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
};
static long ignpoi,j,k,kflag,l,m;
static double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,
    t,u,v,x,xx,pp[35];

    if(mu == muprev) goto S10;
    if(mu < 10.0) goto S120;
/*
     C A S E  A. (RECALCULATION OF S,D,L IF MU HAS CHANGED)
*/
    muprev = mu;
    s = sqrt(mu);
    d = 6.0*mu*mu;
/*
             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
*/
    l = (long) (mu-1.1484);
S10:
/*
     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
*/
    g = mu+s*snorm();
    if(g < 0.0) goto S20;
    ignpoi = (long) (g);
/*
     STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
*/
    if(ignpoi >= l) return ignpoi;
/*
     STEP S. SQUEEZE ACCEPTANCE - Srnd(IR) FOR (0,1)-SAMPLE U
*/
    fk = (double)ignpoi;
    difmuk = mu-fk;
    u = rnd();  /*was ranf -- JKP*/
    if(d*u >= difmuk*difmuk*difmuk) return ignpoi;
S20:
/*
     STEP P. PREPARATIONS FOR STEPS Q AND H.
             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
*/
    if(mu == muold) goto S30;
    muold = mu;
    omega = 0.3989423/s;
    b1 = 4.166667E-2/mu;
    b2 = 0.3*b1*b1;
    c3 = 0.1428571*b1*b2;
    c2 = b2-15.0*c3;
    c1 = b1-6.0*b2+45.0*c3;
    c0 = 1.0-b1+3.0*b2-15.0*c3;
    c = 0.1069/mu;
S30:
    if(g < 0.0) goto S50;
/*
             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
*/
    kflag = 0;
    goto S70;
S40:
/*
     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
*/
    if(fy-u*fy <= py*exp(px-fx)) return ignpoi;
S50:
/*
     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
             (IF T <= -.6744 THEN PK < FK FOR ALL MU >=
