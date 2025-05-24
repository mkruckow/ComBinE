//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
#include <fstream>	//in-/output to/from file
#include "ComBinElib.h" // local header file

//using namespace std;

//global variable
t_ecc eccarray;	// eccarray: values for the integral in eq. 5.14 in Peters (1964)

double mergetime(double M1, double M2, double a0, double e0){
  static double betaconstants = 64.0/5.0*pow(G,3)/pow(c,5);
  double beta,c0,t;

  // after Peters (1964)
  beta = betaconstants*M1*M2*(M1+M2);
  if (e0==0.0){	//circular orbit
    t = pow(a0,4)/(4.0*beta);
    return t;
  }else{	//eccentric orbit
    c0 = a0*(1.0-e0*e0)*pow(e0,-12.0/19.0)*pow(1.0+121.0/304.0*e0*e0,-870.0/2299.0);
    t = 12.0/19.0*pow(c0,4)/beta*intmergetime(e0);
    return t;
  }
}

double intmergetime(double e0){
  int i;	//loop variable
  double eratio;	//eccentricity ratio

  if (eccarray.n==0){	//reread the mergetime table
    destroymergetime();
    readmergetime();
  }
  for (i=max(1,(int)(e0*(eccarray.n-1)-2));i<eccarray.n-1;i++){	//find eccentricity
    if (e0<eccarray.ecc[i]) break;
  }
  eratio = ratio(e0,eccarray.ecc[i],eccarray.ecc[i-1]);
  if ((eratio<0)||(eratio>1)){
    cerr << "#Error: e0=" << e0 << " eccarray.ecc[" << i << "]=" << eccarray.ecc[i] << " eccarray.ecc[" << i-1 << "]=" << eccarray.ecc[i-1] << " eratio=" << eratio << endl;
    screen = true;
    return 0.0;
  }else{
    return eccarray.integral[i]-eratio*(eccarray.integral[i]-eccarray.integral[i-1]);
  }
}

void readmergetime(){
  ifstream infile;	//variable for read table from file
  char text[10000];	//some text
  int err;		//errorcode

  //open mergertime table
  err = sprintf(text,"%smergertime.int",tablelocation);
  if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
  infile.open(text);
  if (infile.fail()){	//if table does not exist create it
    infile.close();
    newinttable();
    infile.open(text);
  }
  if (infile.fail()){
    cerr << "#Error: no mergertime.int";
  }else{
    //reserve memory
    eccarray.ecc = (double *)malloc(sizeof(double));
    eccarray.integral = (double *)malloc(sizeof(double));
    //check if memory allocation fails
    if (eccarray.ecc==NULL) cerr << "#Error: memory allocation failed: eccarray.ecc" << endl;
    if (eccarray.integral==NULL) cerr << "#Error: memory allocation failed: eccarray.integral" << endl;
    //initialize values
    eccarray.n = 0;
    do{
      //enlarge reserved memory
      eccarray.ecc = (double *)realloc(eccarray.ecc,(eccarray.n+1)*sizeof(double));
      eccarray.integral = (double *)realloc(eccarray.integral,(eccarray.n+1)*sizeof(double));
      //check if memory allocation fails
      if (eccarray.ecc==NULL) cerr << "#Error: memory reallocation failed: eccarray.ecc" << endl;
      if (eccarray.integral==NULL) cerr << "#Error: memory reallocation failed: eccarray.integral" << endl;
      //read in new values
      infile >> eccarray.ecc[eccarray.n];
      infile >> eccarray.integral[eccarray.n];
      //increase number of array entries
      eccarray.n++;
    }while (!infile.fail());
    //remove last entry
    eccarray.n--;
    eccarray.ecc = (double *)realloc(eccarray.ecc,(eccarray.n)*sizeof(double));
    eccarray.integral = (double *)realloc(eccarray.integral,(eccarray.n)*sizeof(double));
    //check if memory allocation fails
    if (eccarray.ecc==NULL) cerr << "#Error: memory reallocation failed: eccarray.ecc" << endl;
    if (eccarray.integral==NULL) cerr << "#Error: memory reallocation failed: eccarray.integral" << endl;
    //close lambda table
    infile.close();
  }
}

void destroymergetime(){
  free(eccarray.ecc);
  free(eccarray.integral);
  eccarray.n = 0;
}

void newinttable(){
  int n=10000;	//number of eccentricities
  int s=10000;	//number of supporting points per eccentricity interval
  int i;	//loop variable
  int imax=s*n;	//number of supporting points
  double sum=0.0;	//integralsum
  double delta=1.0/imax;	//interval size
  double ecc;	//current eccentricity
  ofstream outfile;	//enable writing to file the outputdata
  char text[10000];	//some text
  int err;		//errorcode

  cout << "create new mergetime table ...";
  err = sprintf(text,"%smergertime.int",tablelocation);
  if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
  outfile.open(text);	// open file with name "tables/mergertime.int" to write
  outfile << 0.0 << " " << 0.0 << endl;
  for (i=1;i<=imax;i++){
    ecc = delta*i-0.5*delta; // mid-sum
    sum += (pow(ecc,29.0/19.0)*pow((1.0+(121.0/304.0)*ecc*ecc),1181.0/2299.0)*pow(1.0-ecc*ecc,-1.5))*delta;
    if (i%s==0) outfile << delta*i << " " << sum << endl;
  }
  outfile.close();
  cout << " done" << endl;
}
