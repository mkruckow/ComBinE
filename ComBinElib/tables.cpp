//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
//#include <fstream>	//in-/output to/from file
#include "ComBinElib.h"	// local header file

//using namespace std;

//global variables containing the data from the tables
int nstar=0;		// nstar: number of stellar tracks in stararray
int nstar1=0;		// nstar1: number of stellar tracks in stararray1
int nstar2=0;		// nstar2: number of stellar tracks in stararray2
int nhe=0;		// nhe: number of naked helium tracks in hestararray
int nhe1=0;		// nhe1: number of naked helium tracks in hestararray1
t_HRD* stararray=NULL;	// stararray: tracks for stars of different masses
t_HRD* stararray1=NULL;	// stararray1: tracks for MS-stars of different masses
t_HRD* stararray2=NULL;	// stararray2: tracks for Core-Helium-burning-stars of different masses
t_He* hestararray=NULL;	// hestararray: tracks for He-stars of different masses
t_HRD* hestararray1=NULL;	// hestararray1: tracks for He-stars of different masses
t_lambda lambdaarray;	// lambdaarray: lambda values for different masses and radii

char tablelocation[1000] = "./tables/";	//location of the tables in the file system

void readstararray(int kind){
/*kind<-2: old tables; kind==-1 or -2: trk1 and trk2 files; 0<kind<10: H,nonrotating stars; 10<kind<20: H,rotating stars; 20<kind<30: He,nonrotating stars*/
  ifstream infile;	//variable for read table from file
  ifstream findfile;	//variable for found files
  string line;		//a line
  istringstream sline;	//stream of line
  char text[10000];	//some text
  int err;		//errorcode
  int i,j;		//loop variable
  t_HRD tempHRD;	//temporary variable for sort
//  int j,k;		//loop variables
//  double diff;		//difference to next value
/*  double x1,x2,y1,y2,a,b;	//auxiliary variables for extrapolation*/
  double lambdaG,lambdaB;	//lambda values without/with thermal and radiation energy
  int nstar0=0;		//nstar: number of stellar tracks in stararray
  int jump=0;		//jump: indicates the index where a jump in coremass occures because of convection
//  double minvalue=0.0;	//minimal reference value for correction of convective core
  t_HRD* stararray0=NULL;	// stararray: tracks for stars of different masses

  //reserve memory
  stararray0 = (t_HRD *)malloc(sizeof(t_HRD));
  //check if memory allocation fails
  if (stararray0==NULL) cerr << "#Error: memory allocation failed: stararray0" << endl;
  //initialize number of read star tracks
  nstar0 = 0;

  if (kind==1){
    // lists all trk-files in tables/MWstars in Find_MW.trk
    err = sprintf(text,"find %sMWstars/*.trk -type f > Find_MW.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in MWstars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_MW.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==-1){
    // lists all trk1-files in tables/MWstars in Find_MW.trk1
    err = sprintf(text,"find %sMWstars/*.trk1 -type f > Find_MW.trk1",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk1-tables in MWstars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_MW.trk1");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==-2){
    // lists all trk2-files in tables/MWstars in Find_MW.trk2
    err = sprintf(text,"find %sMWstars/*.trk2 -type f > Find_MW.trk2",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk2-tables in MWstars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_MW.trk2");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
//    if (stararray1==NULL) cerr << "#Error: can't create trk2-tables without trk1-tables" << endl;
  }else if (kind==2){
    // lists all trk-files in tables/LMCstars in Find_LMC.trk
    err = sprintf(text,"find %sLMCstars/*.trk -type f > Find_LMC.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in LMCstars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_LMC.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==3){
    // lists all trk-files in tables/SMCstars in Find_SMC.trk
    err = sprintf(text,"find %sSMCstars/*.trk -type f > Find_SMC.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in SMCstars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_SMC.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==4){
    // lists all trk-files in tables/IZw18stars in Find_IZw18.trk
    err = sprintf(text,"find %sIZw18stars/*.trk -type f > Find_IZw18.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in IZw18stars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_IZw18.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==11){
    // lists all trk-files in tables/MWrotstars in Find_MW_rot.trk
    err = sprintf(text,"find %sMWrotstars/*.trk -type f > Find_MW_rot.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in MWrotstars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_MW_rot.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==21){
    // lists all trk-files in tables/MWhestars in Find_MW_he.trk
    err = sprintf(text,"find %sMWhestars/*.trk -type f > Find_MW_he.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in MWhestars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_MW_he.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==22){
    // lists all trk-files in tables/LMChestars in Find_LMC_he.trk
    err = sprintf(text,"find %sLMChestars/*.trk -type f > Find_LMC_he.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in LMChestars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_LMC_he.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==23){
    // lists all trk-files in tables/SMChestars in Find_SMC_he.trk
    err = sprintf(text,"find %sSMChestars/*.trk -type f > Find_SMC_he.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in SMChestars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_SMC_he.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else if (kind==24){
    // lists all trk-files in tables/IZw18hestars in Find_IZw18_he.trk
    err = sprintf(text,"find %sIZw18hestars/*.trk -type f > Find_IZw18_he.trk",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find trk-tables in IZw18hestars, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find_IZw18_he.trk");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }else{
    if (kind>-3) cerr << "#Error: unconsidered kind: " << kind << endl;
    // lists all files in tables/star in Find.star
    err = sprintf(text,"find %sstar/ -type f > Find.star",tablelocation);
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    err = system(text);
    if (err!=0) cerr << "#Error: can't find tables in star, code: " << err << endl;
    //open file with names of all star tables
    err = sprintf(text,"Find.star");
    if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
    findfile.open(text);
  }

  do{
    //get name of next star table
    findfile.getline(text,10000);
    if (findfile.fail()){
      break;
    }else{
//      cout << text << endl;
      //open next star table
      infile.open(text);
      if (debug) cerr << "open " << text << endl;
      //enlarge reserved memory
      stararray0 = (t_HRD *)realloc(stararray0,(nstar0+1)*sizeof(t_HRD));
      //check if memory allocation fails
      if (stararray0==NULL) cerr << "#Error: memory reallocation failed: stararray0" << endl;
      //reserve memory
      stararray0[nstar0].m = (double *)malloc(sizeof(double));
      stararray0[nstar0].t = (double *)malloc(sizeof(double));
      stararray0[nstar0].r = (double *)malloc(sizeof(double));
      stararray0[nstar0].cm = (double *)malloc(sizeof(double));
      stararray0[nstar0].llum = (double *)malloc(sizeof(double));
      stararray0[nstar0].lteff = (double *)malloc(sizeof(double));
      stararray0[nstar0].lambda = (double *)malloc(sizeof(double));
      stararray0[nstar0].cc = (double *)malloc(sizeof(double));
      stararray0[nstar0].cf = (double *)malloc(sizeof(double));
      stararray0[nstar0].omega = (double *)malloc(sizeof(double));
      //check if memory allocation fails
      if (stararray0[nstar0].m==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].m" << endl;
      if (stararray0[nstar0].t==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].t" << endl;
      if (stararray0[nstar0].r==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].r" << endl;
      if (stararray0[nstar0].cm==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].cm" << endl;
      if (stararray0[nstar0].llum==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].llum" << endl;
      if (stararray0[nstar0].lteff==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].lteff" << endl;
      if (stararray0[nstar0].lambda==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].lambda" << endl;
      if (stararray0[nstar0].cc==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].cc" << endl;
      if (stararray0[nstar0].cf==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].cf" << endl;
      if (stararray0[nstar0].omega==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].omega" << endl;
      //initialize values
      stararray0[nstar0].TAMS = -1;
      stararray0[nstar0].n = 1;
      stararray0[nstar0].m[0] = -1000.0;
      stararray0[nstar0].t[0] = -1000.0;
      stararray0[nstar0].r[0] = -1000.0;
      stararray0[nstar0].cm[0] = -1000.0;
      stararray0[nstar0].llum[0] = -1000.0;
      stararray0[nstar0].lteff[0] = -1000.0;
      stararray0[nstar0].lambda[0] = -1000.0;
      stararray0[nstar0].cc[0] = -1000.0;
      stararray0[nstar0].cf[0] = -1000.0;
      stararray0[nstar0].omega[0] = -1000.0;
//      do{
      while(getline(infile,line)){
        sline.clear();
        sline.str(line);
        //enlarge reserved memory
        stararray0[nstar0].m = (double *)realloc(stararray0[nstar0].m,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "m";
        stararray0[nstar0].t = (double *)realloc(stararray0[nstar0].t,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "t";
        stararray0[nstar0].r = (double *)realloc(stararray0[nstar0].r,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "r";
        stararray0[nstar0].cm = (double *)realloc(stararray0[nstar0].cm,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "c";
        stararray0[nstar0].llum = (double *)realloc(stararray0[nstar0].llum,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "L";
        stararray0[nstar0].lteff = (double *)realloc(stararray0[nstar0].lteff,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "T";
        stararray0[nstar0].lambda = (double *)realloc(stararray0[nstar0].lambda,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "l";
        stararray0[nstar0].cc = (double *)realloc(stararray0[nstar0].cc,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "C";
        stararray0[nstar0].cf = (double *)realloc(stararray0[nstar0].cf,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "f";
        stararray0[nstar0].omega = (double *)realloc(stararray0[nstar0].omega,(stararray0[nstar0].n+1)*sizeof(double));
        if (debug) cerr << stararray0[nstar0].n << "o, ";
        //check if memory allocation fails
        if (stararray0[nstar0].m==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].m" << endl;
        if (stararray0[nstar0].t==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].t" << endl;
        if (stararray0[nstar0].r==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].r" << endl;
        if (stararray0[nstar0].cm==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cm" << endl;
        if (stararray0[nstar0].llum==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].llum" << endl;
        if (stararray0[nstar0].lteff==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].lteff" << endl;
        if (stararray0[nstar0].lambda==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].lambda" << endl;
        if (stararray0[nstar0].cc==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cc" << endl;
        if (stararray0[nstar0].cf==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cf" << endl;
        if (stararray0[nstar0].omega==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].omega" << endl;
        //read in new values
//        infile >> stararray0[nstar0].m[stararray0[nstar0].n];
        sline >> stararray0[nstar0].m[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].m[stararray0[nstar0].n] = fmax(10.0*accuracy,stararray0[nstar0].m[stararray0[nstar0].n-1]);
        if (stararray0[nstar0].m[stararray0[nstar0].n]<0){	//TAMS index detected
          if ((kind<0)||(kind>20)) cerr << "#Error: no TAMS expected, negative mass? value=" << stararray0[nstar0].m[stararray0[nstar0].n] << endl;
          stararray0[nstar0].TAMS = -(int)stararray0[nstar0].m[stararray0[nstar0].n];
          goto endfile;
        }
//        infile >> stararray0[nstar0].t[stararray0[nstar0].n];
        sline >> stararray0[nstar0].t[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].t[stararray0[nstar0].n] = fmax(0.0,stararray0[nstar0].t[stararray0[nstar0].n-1]);
//        infile >> stararray0[nstar0].r[stararray0[nstar0].n];
        sline >> stararray0[nstar0].r[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].r[stararray0[nstar0].n] = fmax(accuracy,stararray0[nstar0].r[stararray0[nstar0].n-1]);
//        infile >> stararray0[nstar0].cm[stararray0[nstar0].n];
        sline >> stararray0[nstar0].cm[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].cm[stararray0[nstar0].n] = fmax(accuracy,stararray0[nstar0].cm[stararray0[nstar0].n-1]);
        if ((stararray0[nstar0].n>1)&&(stararray0[nstar0].cm[stararray0[nstar0].n-1]*100.0<stararray0[nstar0].cm[stararray0[nstar0].n])){	//core mass increased by 2 orders of magnitude
          jump = stararray0[nstar0].n;
          if (debug) cout << "set jump=" << jump << " cm[" << stararray0[nstar0].n-1 << "]=" << stararray0[nstar0].cm[stararray0[nstar0].n-1] << "Msun cm[" << stararray0[nstar0].n << "]=" << stararray0[nstar0].cm[stararray0[nstar0].n] << "Msun" << endl;
        }
//        infile >> stararray0[nstar0].llum[stararray0[nstar0].n];
        sline >> stararray0[nstar0].llum[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].llum[stararray0[nstar0].n] = fmax(-100.0,stararray0[nstar0].llum[stararray0[nstar0].n-1]);
//        infile >> stararray0[nstar0].lteff[stararray0[nstar0].n];
        sline >> stararray0[nstar0].lteff[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].lteff[stararray0[nstar0].n] = fmax(-100.0,stararray0[nstar0].lteff[stararray0[nstar0].n-1]);
        if (kind>=-2){
//          infile >> lambdaG;
          sline >> lambdaG;
          if (sline.fail()) lambdaG = 0.0;
          if (lambdaG==0.0) lambdaG = 1.0e-99;
//          infile >> lambdaB;
          sline >> lambdaB;
          if (sline.fail()) lambdaB = 0.0;
          if (lambdaB==0.0) lambdaB = 1.0e-99;
        }else{
          lambdaG = abs(lambda_const);
          lambdaB = 2.0*lambdaG;
        }
        if ((alphaTH>=0.0)&&(alphaTH<=1.0)){	//interpolate, see eq. 6 Dewi&Tauris 2000
          stararray0[nstar0].lambda[stararray0[nstar0].n] = lambdaG*lambdaB/(alphaTH*(lambdaG-lambdaB)+lambdaB);
        }else{	//virial equilibrium for alphaTH=2
          stararray0[nstar0].lambda[stararray0[nstar0].n] = alphaTH*lambdaG;
        }
        if ((stararray0[nstar0].lambda[stararray0[nstar0].n]<0)&&(alphaTH!=-1)){
          if (screen){
            cerr << "#Warning: replace negative lambda with 1.0e+10: ";
            if ((kind>0)&&(kind<20)) cerr << "stararray";
            else if (kind==-1) cerr << "stararray1";
            else if (kind==-2) cerr << "stararray2";
            else if ((kind>20)&&(kind<30)) cerr << "hestararray1";
            else cerr << "old stararray";
            cerr << "[" << nstar0 << "].lambda[" << stararray0[nstar0].n << "]=" << stararray0[nstar0].lambda[stararray0[nstar0].n] << endl;
          }
          stararray0[nstar0].lambda[stararray0[nstar0].n] = 1.0e+10;
        }
//        infile >> stararray0[nstar0].cc[stararray0[nstar0].n];
        sline >> stararray0[nstar0].cc[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].cc[stararray0[nstar0].n] = fmax(0.01*accuracy,stararray0[nstar0].cc[stararray0[nstar0].n-1]);
//        infile >> stararray0[nstar0].cf[stararray0[nstar0].n];
        sline >> stararray0[nstar0].cf[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].cf[stararray0[nstar0].n] = fmax(0.0,stararray0[nstar0].cf[stararray0[nstar0].n-1]);
//        infile >> stararray0[nstar0].omega[stararray0[nstar0].n];
        sline >> stararray0[nstar0].omega[stararray0[nstar0].n];
        if (sline.fail()) stararray0[nstar0].omega[stararray0[nstar0].n] = fmax(0.0,stararray0[nstar0].omega[stararray0[nstar0].n-1]);
//        infile.ignore(10000,'\n');	//ignore rest of the line(all the characters until '\n')
//        cout << stararray0[nstar0].m[stararray0[nstar0].n] << stararray0[nstar0].t[stararray0[nstar0].n] << stararray0[nstar0].r[stararray0[nstar0].n] << stararray0[nstar0].cm[stararray0[nstar0].n] << stararray0[nstar0].llum[stararray0[nstar0].n] << stararray0[nstar0].lteff[stararray0[nstar0].n] << endl;
        if (infile.fail()){	//if an error while last reading occured, e.g. end of file is reached
          endfile:	//goto-marker
          if (debug) cerr << endl << infile.fail() << endl;
          //remove last entry
          stararray0[nstar0].m = (double *)realloc(stararray0[nstar0].m,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].t = (double *)realloc(stararray0[nstar0].t,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].r = (double *)realloc(stararray0[nstar0].r,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].cm = (double *)realloc(stararray0[nstar0].cm,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].llum = (double *)realloc(stararray0[nstar0].llum,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].lteff = (double *)realloc(stararray0[nstar0].lteff,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].lambda = (double *)realloc(stararray0[nstar0].lambda,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].cc = (double *)realloc(stararray0[nstar0].cc,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].cf = (double *)realloc(stararray0[nstar0].cf,(stararray0[nstar0].n)*sizeof(double));
          stararray0[nstar0].omega = (double *)realloc(stararray0[nstar0].omega,(stararray0[nstar0].n)*sizeof(double));
          //check if memory allocation fails
          if (stararray0[nstar0].m==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].m" << endl;
          if (stararray0[nstar0].t==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].t" << endl;
          if (stararray0[nstar0].r==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].r" << endl;
          if (stararray0[nstar0].cm==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cm" << endl;
          if (stararray0[nstar0].llum==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].llum" << endl;
          if (stararray0[nstar0].lteff==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].lteff" << endl;
          if (stararray0[nstar0].lambda==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].lambda" << endl;
          if (stararray0[nstar0].cc==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cc" << endl;
          if (stararray0[nstar0].cf==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cf" << endl;
          if (stararray0[nstar0].omega==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].omega" << endl;
          break;
        }
        //check for maximum values
/*        if (stararray0[nstar0].m[stararray0[nstar0].n]>stararray0[nstar0].m[0]){
          stararray0[nstar0].m[0] = stararray0[nstar0].m[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].t[stararray0[nstar0].n]>stararray0[nstar0].t[0]){
          stararray0[nstar0].t[0] = stararray0[nstar0].t[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].r[stararray0[nstar0].n]>stararray0[nstar0].r[0]){
          stararray0[nstar0].r[0] = stararray0[nstar0].r[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].cm[stararray0[nstar0].n]>stararray0[nstar0].cm[0]){
          stararray0[nstar0].cm[0] = stararray0[nstar0].cm[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].llum[stararray0[nstar0].n]>stararray0[nstar0].llum[0]){
          stararray0[nstar0].llum[0] = stararray0[nstar0].llum[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].lteff[stararray0[nstar0].n]>stararray0[nstar0].lteff[0]){
          stararray0[nstar0].lteff[0] = stararray0[nstar0].lteff[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].lambda[stararray0[nstar0].n]>stararray0[nstar0].lambda[0]){
          stararray0[nstar0].lambda[0] = stararray0[nstar0].lambda[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].cc[stararray0[nstar0].n]>stararray0[nstar0].cc[0]){
          stararray0[nstar0].cc[0] = stararray0[nstar0].cc[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].cf[stararray0[nstar0].n]>stararray0[nstar0].cf[0]){
          stararray0[nstar0].cf[0] = stararray0[nstar0].cf[stararray0[nstar0].n];
        }
        if (stararray0[nstar0].omega[stararray0[nstar0].n]>stararray0[nstar0].omega[0]){
          stararray0[nstar0].omega[0] = stararray0[nstar0].omega[stararray0[nstar0].n];
        }*/
        //increase number of array entries
        stararray0[nstar0].n++;
//      }while (!infile.fail());
      }
      if (((fabs(stararray0[nstar0].t[1])<accuracy)&&(stararray0[nstar0].cm[1]<=accuracy))||(kind==-2)){
        for (i=1;i<stararray0[nstar0].n;i++){
          stararray0[nstar0].m[i-1] = stararray0[nstar0].m[i];
          stararray0[nstar0].t[i-1] = stararray0[nstar0].t[i];
          stararray0[nstar0].r[i-1] = stararray0[nstar0].r[i];
          stararray0[nstar0].cm[i-1] = stararray0[nstar0].cm[i];
          stararray0[nstar0].llum[i-1] = stararray0[nstar0].llum[i];
          stararray0[nstar0].lteff[i-1] = stararray0[nstar0].lteff[i];
          stararray0[nstar0].lambda[i-1] = stararray0[nstar0].lambda[i];
          stararray0[nstar0].cc[i-1] = stararray0[nstar0].cc[i];
          stararray0[nstar0].cf[i-1] = stararray0[nstar0].cf[i];
          stararray0[nstar0].omega[i-1] = stararray0[nstar0].omega[i];
        }
        stararray0[nstar0].n--;
        stararray0[nstar0].TAMS--;
        jump--;
      }else{	//set values for t=0
        stararray0[nstar0].m[0] = stararray0[nstar0].m[1];
        stararray0[nstar0].t[0] = 0.0;
        stararray0[nstar0].r[0] = stararray0[nstar0].r[1];
        stararray0[nstar0].cm[0] = accuracy;
        stararray0[nstar0].llum[0] = stararray0[nstar0].llum[1];
        stararray0[nstar0].lteff[0] = stararray0[nstar0].lteff[1];
        stararray0[nstar0].lambda[0] = stararray0[nstar0].lambda[1];
        stararray0[nstar0].cc[0] = 0.01*accuracy;
        stararray0[nstar0].cf[0] = stararray0[nstar0].cf[1];
        stararray0[nstar0].omega[0] = stararray0[nstar0].omega[1];
        for (i=1;i<stararray0[nstar0].n;i++){	//check that new values at t=0 are the minimum values
          if (stararray0[nstar0].t[i]<stararray0[nstar0].t[0]) stararray0[nstar0].t[i] = stararray0[nstar0].t[0];
          if (stararray0[nstar0].cm[i]<stararray0[nstar0].cm[0]) stararray0[nstar0].cm[i] = stararray0[nstar0].cm[0];
          if (stararray0[nstar0].cc[i]<stararray0[nstar0].cc[0]) stararray0[nstar0].cc[i] = stararray0[nstar0].cc[0];
        }
      }
//      debug = (stararray0[nstar0].m[0]==5.0);
/*      if ((jump>0)&&((kind==-1)||((kind>0)&&(kind<20)))){	//convective core from comparision to detailed modelling of Pablo	//||(kind==21) not proven for He-stars
        minvalue = stararray0[nstar0].llum[0];	//set minvalue
        for (i=1;i<stararray0[nstar0].n;i++){	//find time index of time
          if (stararray0[nstar0].llum[i]<minvalue) minvalue = stararray0[nstar0].llum[i];	//set minvalue
          if (stararray0[nstar0].cm[i]/stararray0[nstar0].m[i]>0.01) break;
        }
        if (debug) cerr << "minvalue=" << minvalue << " i=" << i << endl;
        if (i<jump){
          i = jump;
          if (debug) cerr << "inimass=" << stararray0[nstar0].m[0] << "Msun TAMS=" << stararray0[nstar0].TAMS << " i=" << i << " jump=" << jump << endl;
        }
        if ((stararray0[nstar0].TAMS+3<stararray0[nstar0].n)&&(stararray0[nstar0].TAMS>0)){
          i = stararray0[nstar0].TAMS;	//if the TAMS is reached it is better to use its value as reference to MS core mass
          while (stararray0[nstar0].llum[i-1]>stararray0[nstar0].llum[i]){	//avoid decrease of Luminosity towards TAMS
            i--;
            if (i<1) break;
          }
        }
        if ((i<jump)&&(debug)) cerr << "inimass=" << stararray0[nstar0].m[0] << "Msun TAMS=" << stararray0[nstar0].TAMS << " i=" << i << " jump=" << jump << endl;
        //enlarge reserved memory
        stararray0[nstar0].m = (double *)realloc(stararray0[nstar0].m,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].t = (double *)realloc(stararray0[nstar0].t,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].r = (double *)realloc(stararray0[nstar0].r,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].cm = (double *)realloc(stararray0[nstar0].cm,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].llum = (double *)realloc(stararray0[nstar0].llum,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].lteff = (double *)realloc(stararray0[nstar0].lteff,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].lambda = (double *)realloc(stararray0[nstar0].lambda,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].cc = (double *)realloc(stararray0[nstar0].cc,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].cf = (double *)realloc(stararray0[nstar0].cf,(stararray0[nstar0].n+1)*sizeof(double));
        stararray0[nstar0].omega = (double *)realloc(stararray0[nstar0].omega,(stararray0[nstar0].n+1)*sizeof(double));
        //check if memory allocation fails
        if (stararray0[nstar0].m==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].m" << endl;
        if (stararray0[nstar0].t==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].t" << endl;
        if (stararray0[nstar0].r==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].r" << endl;
        if (stararray0[nstar0].cm==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cm" << endl;
        if (stararray0[nstar0].llum==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].llum" << endl;
        if (stararray0[nstar0].lteff==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].lteff" << endl;
        if (stararray0[nstar0].lambda==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].lambda" << endl;
        if (stararray0[nstar0].cc==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cc" << endl;
        if (stararray0[nstar0].cf==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].cf" << endl;
        if (stararray0[nstar0].omega==NULL) cerr << "#Error: memory reallocation failed: stararray0[nstar0].omega" << endl;
        //copy values
        for (j=stararray0[nstar0].n;j>1;j--){
          stararray0[nstar0].m[j] = stararray0[nstar0].m[j-1];
          stararray0[nstar0].t[j] = stararray0[nstar0].t[j-1];
          stararray0[nstar0].r[j] = stararray0[nstar0].r[j-1];
          stararray0[nstar0].cm[j] = stararray0[nstar0].cm[j-1];
          stararray0[nstar0].llum[j] = stararray0[nstar0].llum[j-1];
          stararray0[nstar0].lteff[j] = stararray0[nstar0].lteff[j-1];
          stararray0[nstar0].lambda[j] = stararray0[nstar0].lambda[j-1];
          stararray0[nstar0].cc[j] = stararray0[nstar0].cc[j-1];
          stararray0[nstar0].cf[j] = stararray0[nstar0].cf[j-1];
          stararray0[nstar0].omega[j] = stararray0[nstar0].omega[j-1];
        }
        stararray0[nstar0].m[1] = (1.0-1.0e-5)*stararray0[nstar0].m[0]+1.0e-5*stararray0[nstar0].m[2];
        stararray0[nstar0].t[1] = (1.0-1.0e-5)*stararray0[nstar0].t[0]+1.0e-5*stararray0[nstar0].t[2];
        stararray0[nstar0].r[1] = (1.0-1.0e-5)*stararray0[nstar0].r[0]+1.0e-5*stararray0[nstar0].r[2];
        stararray0[nstar0].cm[1] = (1.0-1.0e-5)*stararray0[nstar0].cm[0]+1.0e-5*stararray0[nstar0].cm[2];
        stararray0[nstar0].llum[1] = (1.0-1.0e-5)*stararray0[nstar0].llum[0]+1.0e-5*stararray0[nstar0].llum[2];
        stararray0[nstar0].lteff[1] = (1.0-1.0e-5)*stararray0[nstar0].lteff[0]+1.0e-5*stararray0[nstar0].lteff[2];
        stararray0[nstar0].lambda[1] = (1.0-1.0e-5)*stararray0[nstar0].lambda[0]+1.0e-5*stararray0[nstar0].lambda[2];
        stararray0[nstar0].cc[1] = (1.0-1.0e-5)*stararray0[nstar0].cc[0]+1.0e-5*stararray0[nstar0].cc[2];
        stararray0[nstar0].cf[1] = (1.0-1.0e-5)*stararray0[nstar0].cf[0]+1.0e-5*stararray0[nstar0].cf[2];
        stararray0[nstar0].omega[1] = (1.0-1.0e-5)*stararray0[nstar0].omega[0]+1.0e-5*stararray0[nstar0].omega[2];
        stararray0[nstar0].n++;
        stararray0[nstar0].TAMS++;
        i++;
        if (debug) cerr << "stararray0[nstar0].cm=(" << stararray0[nstar0].cm[0];
        for (j=1;j<i;j++){
//          stararray0[nstar0].cm[j] = stararray0[nstar0].t[j]/stararray0[nstar0].t[i]*stararray0[nstar0].cm[i];	//estimate He(-2<kind<20)/C(20<kind<30) ammount of convective core
          if (stararray0[nstar0].llum[j]>=minvalue){	//estimate He(-2<kind<20)/C(20<kind<30) ammount of convective core
            if (j==1){
              if (stararray0[nstar0].m[0]<15.0) stararray0[nstar0].cm[j] = fmax(stararray0[nstar0].cm[0],ratio(stararray0[nstar0].llum[j],minvalue,stararray0[nstar0].llum[i])*stararray0[nstar0].cm[i]);
              else if (stararray0[nstar0].m[0]<45.0) stararray0[nstar0].cm[j] = fmax(stararray0[nstar0].cm[0],pow(ratio(stararray0[nstar0].llum[j],minvalue,stararray0[nstar0].llum[i]),(45.0-stararray0[nstar0].m[0])/(45.0-15.0))*stararray0[nstar0].cm[i]);
              else stararray0[nstar0].cm[j] = fmax(stararray0[nstar0].cm[j],stararray0[nstar0].cm[i]);
            }else{
              if (stararray0[nstar0].m[0]<15.0) stararray0[nstar0].cm[j] = fmax(stararray0[nstar0].cm[0],ratio(stararray0[nstar0].llum[j],minvalue,stararray0[nstar0].llum[i])*stararray0[nstar0].cm[i]);
              else if (stararray0[nstar0].m[0]<45.0) stararray0[nstar0].cm[j] = fmax(stararray0[nstar0].cm[0],pow(ratio(stararray0[nstar0].llum[j],minvalue,stararray0[nstar0].llum[i]),(45.0-stararray0[nstar0].m[0])/(45.0-15.0))*stararray0[nstar0].cm[i]);
              else stararray0[nstar0].cm[j] = fmax(stararray0[nstar0].cm[0],(1.0-1.0e-5)*stararray0[nstar0].cm[i]);
            }
          }
          if (debug) cerr << ", " << stararray0[nstar0].cm[j];
          if (stararray0[nstar0].cm[j]<0.0) cerr << "#Error: negative core mass: inimass=" << stararray0[nstar0].m[0] << "Msun cm[" << j << "]=" << stararray0[nstar0].cm[j] << "Msun jump=" << jump << endl;
        }
        if (debug) cerr << ")" << endl;
        jump = 0;	//reset jump-index
      }*/
      //close star table
      infile.close();
/*      if ((stararray0[nstar0].inimass>1.7)&&(stararray0[nstar0].inimass<1.9)){ debug = true; cerr << "debug on: kind=" << kind << endl;}*/
      smoothtable(stararray0[nstar0].m, stararray0[nstar0].n,-1);
/*      if ((stararray0[nstar0].inimass>1.7)&&(stararray0[nstar0].inimass<1.9)) debug = false;*/
      smoothtable(stararray0[nstar0].t, stararray0[nstar0].n,1);
      smoothtable(stararray0[nstar0].r, stararray0[nstar0].n,0);
      smoothtable(stararray0[nstar0].cm, stararray0[nstar0].n,1);
      smoothtable(stararray0[nstar0].llum, stararray0[nstar0].n,0);
      smoothtable(stararray0[nstar0].lteff, stararray0[nstar0].n,0);
      smoothtable(stararray0[nstar0].lambda, stararray0[nstar0].n,0);
      smoothtable(stararray0[nstar0].cc, stararray0[nstar0].n,1);
      smoothtable(stararray0[nstar0].cf, stararray0[nstar0].n,0);
      smoothtable(stararray0[nstar0].omega, stararray0[nstar0].n,0);
      if (kind==-2) stararray0[nstar0].TAMS = 0;
      else if (stararray0[nstar0].TAMS<0) stararray0[nstar0].TAMS = stararray0[nstar0].n;
      else if (stararray0[nstar0].TAMS>stararray0[nstar0].n) cerr << "#Error: stararray0[" << nstar0 << "].TAMS=" << stararray0[nstar0].TAMS << ">stararray0[" << nstar0 << "].n=" << stararray0[nstar0].n << endl;
      //set initial mass for sorting
      if (kind==-2){
        for (i=0;i<nstar1;i++){
          if (fabs(stararray0[nstar0].m[0]-stararray1[i].m[stararray1[i].n-1])<accuracy) break;
        }
        if (i<nstar1){
          stararray0[nstar0].inimass = stararray1[i].inimass;
          for (j=1;j<stararray0[nstar0].n;j++){
            if (stararray0[nstar0].m[j]==stararray0[nstar0].m[0]) stararray0[nstar0].m[j] = stararray1[i].m[stararray1[i].n-1];
          }
          stararray0[nstar0].m[0] = stararray1[i].m[stararray1[i].n-1];
        }else{
          cerr << "#Error: no corresponding trk1 found" << endl;
        }
      }else{
        stararray0[nstar0].inimass = stararray0[nstar0].m[0];
      }
      for (i=nstar0;i>0;i--){	//sort stararray0 from low initial mass to high initial mass
        if (stararray0[i].inimass<stararray0[i-1].inimass){	//change order if new mass is lower than the one before
          tempHRD = stararray0[i];
          stararray0[i] = stararray0[i-1];
          stararray0[i-1] = tempHRD;
//          cout << "change: " << i << " with " << i-1 << endl;
        }else{
          break;
        }
      }
/*      if (nstar0>0){
        if (stararray0[nstar0].n!=stararray0[0].n) cerr << "#Error: star tables have diffent length" << endl;
      }*/
/*      if (nstar0==0){
        for (i=0;i<stararray0[nstar0].n;i++) cout << stararray0[nstar0].t[i] << "yr, ";
        cout << endl;
      }*/
/*      for (i=0;i<stararray0[nstar0].n-1;i++){
        j = 1;
        while (stararray0[nstar0].m[i]==stararray0[nstar0].m[i+j]){
          j++;
          if (i+j>=stararray0[nstar0].n) break;
        }
        if (i+j<stararray0[nstar0].n) diff = stararray0[nstar0].m[i+j]-stararray0[nstar0].m[i];
        else if(i>0) diff = stararray0[nstar0].m[i]-stararray0[nstar0].m[i-1];
        else diff = stararray0[nstar0].m[i]*1.0E-10;
        for (k=1;k<j;k++){
          stararray0[nstar0].m[i+k] += diff*((double)k)/((double)j);
        }
        i += j-1;
      }*/
/*      if (nstar0==0){
        for (i=0;i<stararray0[nstar0].n;i++) cout << stararray0[nstar0].t[i] << "yr, ";
        cout << endl;
      }*/
      //increase number of star tables
      nstar0++;
    }
  }while (!findfile.fail());
  //close find-file
  findfile.close();
/*  cout << "[0]: n=" << stararray0[0].n << " m[0]" << stararray0[0].m[0] << " t[0]" << stararray0[0].t[0] << " r[0]" << stararray0[0].r[0] << " core mass(cm)[0]" << stararray0[0].cm[0] << " log(lum/Lsun)(llum)[0]" << stararray0[0].llum[0] << " log(T_eff/K)(lteff)[0]" << stararray0[0].lteff[0] << " lambda[0]" << stararray0[0].lambda[0] << " carbon core mass(cc)[0]" << stararray0[0].cc[0] << " concentration factor(cf)[0]" << stararray0[0].cf[0] << " angular velocity(1/yr)(omega)[0]" << stararray0[0].omega[0] << endl;
  cout << "[0]: n=" << stararray0[1].n << " m[1]" << stararray0[0].m[1] << " t[1]" << stararray0[0].t[1] << " r[1]" << stararray0[0].r[1] << " core mass(cm)[1]" << stararray0[0].cm[1] << " log(lum/Lsun)(llum)[1]" << stararray0[0].llum[1] << " log(T_eff/K)(lteff)[1]" << stararray0[0].lteff[1] << " lambda[1]" << stararray0[0].lambda[1] << " carbon core mass(cc)[1]" << stararray0[0].cc[1] << " concentration factor(cf)[1]" << stararray0[0].cf[1] << " angular velocity(1/yr)(omega)[1]" << stararray0[0].omega[1] << endl;
  cout << "[0]: n=" << stararray0[0].n << " m[n-1]" << stararray0[0].m[stararray0[0].n-1] << " t[n-1]" << stararray0[0].t[stararray0[0].n-1] << " r[n-1]" << stararray0[0].r[stararray0[0].n-1] << " core mass(cm)[n-1]" << stararray0[0].cm[stararray0[0].n-1] << " log(lum/Lsun)(llum)[n-1]" << stararray0[0].llum[stararray0[0].n-1] << " log(T_eff/K)(lteff)[n-1]" << stararray0[0].lteff[stararray0[0].n-1] << " lambda[n-1]" << stararray0[0].lambda[stararray0[0].n-1] << " carbon core mass(cc)[n-1]" << stararray0[0].cc[stararray0[0].n-1] << " concentration factor(cf)[n-1]" << stararray0[0].cf[stararray0[0].n-1] << " angular velocity(1/yr)(omega)[n-1]" << stararray0[0].omega[stararray0[0].n-1] << endl;*/
/*  while (stararray0[nstar0-1].m[1]<1000.0){ //create extrapolation grid
    //enlarge reserved memory
    stararray0 = (t_HRD *)realloc(stararray0,(nstar0+1)*sizeof(t_HRD));
    //check if memory allocation fails
    if (stararray0==NULL) cerr << "#Error: memory reallocation failed: stararray0" << endl;
    //take lenght from previous table
    stararray0[nstar0].n = stararray0[nstar0-1].n;
    //reserve memory
    stararray0[nstar0].m = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].t = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].r = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].cm = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].llum = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].lteff = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].lambda = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].cc = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].cf = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    stararray0[nstar0].omega = (double *)malloc(stararray0[nstar0].n*sizeof(double));
    //check if memory allocation fails
    if (stararray0[nstar0].m==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].m" << endl;
    if (stararray0[nstar0].t==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].t" << endl;
    if (stararray0[nstar0].r==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].r" << endl;
    if (stararray0[nstar0].cm==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].cm" << endl;
    if (stararray0[nstar0].llum==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].llum" << endl;
    if (stararray0[nstar0].lteff==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].lteff" << endl;
    if (stararray0[nstar0].lambda==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].lambda" << endl;
    if (stararray0[nstar0].cc==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].cc" << endl;
    if (stararray0[nstar0].cf==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].cf" << endl;
    if (stararray0[nstar0].omega==NULL) cerr << "#Error: memory allocation failed: stararray0[nstar0].omega" << endl;
    //extrapolate values of the last point
    stararray0[nstar0].m[1] = 100.0*ceil(0.01*stararray0[nstar0-1].m[1])+100.0;
    x1 = stararray0[nstar0-2].m[1];
    x2 = stararray0[nstar0-1].m[1];
    stararray0[nstar0].m[stararray0[nstar0].n-1] = stararray0[nstar0-1].m[stararray0[nstar0-1].n-1]*stararray0[nstar0].m[1]/stararray0[nstar0-1].m[1];
    y1 = stararray0[nstar0-2].t[stararray0[nstar0-2].n-1]-2.5e+6;	//stay above 2.5Myr
    y2 = stararray0[nstar0-1].t[stararray0[nstar0-1].n-1]-2.5e+6;	//stay above 2.5Myr
    a = exp((log(y1)*log(x2)-log(y2)*log(x1))/(log(y1)-log(y2)));
    b = log(y1/y2)/log(x2/x1);
//    cout << "a=" << a << " b=" << b << endl;
    stararray0[nstar0].t[stararray0[nstar0].n-1] = 2.5e+6+pow(a/stararray0[nstar0].m[1],b);	//stay above 2.5Myr
    y1 = stararray0[nstar0-2].r[1]-Schwarzschildradius(x1);	//stay above Schwarzschildradius
    y2 = stararray0[nstar0-1].r[1]-Schwarzschildradius(x2);	//stay above Schwarzschildradius
    a = exp((log(y1)*log(x2)-log(y2)*log(x1))/(log(y1)-log(y2)));
    b = log(y1/y2)/log(x2/x1);
    stararray0[nstar0].r[stararray0[nstar0].n-1] = (Schwarzschildradius(stararray0[nstar0].m[1])+pow(a/stararray0[nstar0].m[1],b))*stararray0[nstar0-1].r[stararray0[nstar0-1].n-1]/stararray0[nstar0-1].r[1];	//stay above Schwarzschildradius
    stararray0[nstar0].cm[stararray0[nstar0].n-1] = stararray0[nstar0-1].cm[stararray0[nstar0-1].n-1]*stararray0[nstar0].m[1]/stararray0[nstar0-1].m[1];
    stararray0[nstar0].llum[stararray0[nstar0].n-1] = stararray0[nstar0-2].llum[stararray0[nstar0-2].n-1]+(stararray0[nstar0].m[1]-x1)/(x2-x1)*(stararray0[nstar0-1].llum[stararray0[nstar0-1].n-1]-stararray0[nstar0-2].llum[stararray0[nstar0-2].n-1]);
    stararray0[nstar0].lteff[stararray0[nstar0].n-1] = stararray0[nstar0-2].lteff[stararray0[nstar0-2].n-1]+(stararray0[nstar0].m[1]-x1)/(x2-x1)*(stararray0[nstar0-1].lteff[stararray0[nstar0-1].n-1]-stararray0[nstar0-2].lteff[stararray0[nstar0-2].n-1]);
    stararray0[nstar0].lambda[stararray0[nstar0].n-1] = stararray0[nstar0-2].lambda[stararray0[nstar0-2].n-1]+(stararray0[nstar0].m[1]-x1)/(x2-x1)*(stararray0[nstar0-1].lambda[stararray0[nstar0-1].n-1]-stararray0[nstar0-2].lambda[stararray0[nstar0-2].n-1]);
    stararray0[nstar0].cc[stararray0[nstar0].n-1] = stararray0[nstar0-1].cc[stararray0[nstar0-1].n-1]*stararray0[nstar0].m[1]/stararray0[nstar0-1].m[1];
    stararray0[nstar0].cf[stararray0[nstar0].n-1] = stararray0[nstar0-1].cf[stararray0[nstar0-1].n-1]*stararray0[nstar0].m[1]/stararray0[nstar0-1].m[1];
    stararray0[nstar0].omega[stararray0[nstar0].n-1] = stararray0[nstar0-1].omega[stararray0[nstar0-1].n-1]*stararray0[nstar0].m[1]/stararray0[nstar0-1].m[1];
    for (i=1;i<stararray0[nstar0].n;i++){ //calculate extrapolation by keeping the ratio to last point as for the previous table
      stararray0[nstar0].m[i] = stararray0[nstar0].m[stararray0[nstar0].n-1]*stararray0[nstar0-1].m[i]/stararray0[nstar0-1].m[stararray0[nstar0-1].n-1];
      stararray0[nstar0].t[i] = stararray0[nstar0].t[stararray0[nstar0].n-1]*stararray0[nstar0-1].t[i]/stararray0[nstar0-1].t[stararray0[nstar0-1].n-1];
      stararray0[nstar0].r[i] = stararray0[nstar0].r[stararray0[nstar0].n-1]*stararray0[nstar0-1].r[i]/stararray0[nstar0-1].r[stararray0[nstar0-1].n-1];
      stararray0[nstar0].cm[i] = stararray0[nstar0].cm[stararray0[nstar0].n-1]*stararray0[nstar0-1].cm[i]/stararray0[nstar0-1].cm[stararray0[nstar0-1].n-1];
      stararray0[nstar0].llum[i] = stararray0[nstar0].llum[stararray0[nstar0].n-1]*stararray0[nstar0-1].llum[i]/stararray0[nstar0-1].llum[stararray0[nstar0-1].n-1];
      stararray0[nstar0].lteff[i] = stararray0[nstar0].lteff[stararray0[nstar0].n-1]*stararray0[nstar0-1].lteff[i]/stararray0[nstar0-1].lteff[stararray0[nstar0-1].n-1];
      stararray0[nstar0].lambda[i] = stararray0[nstar0].lambda[stararray0[nstar0].n-1]*stararray0[nstar0-1].lambda[i]/stararray0[nstar0-1].lambda[stararray0[nstar0-1].n-1];
      stararray0[nstar0].cc[i] = stararray0[nstar0].cc[stararray0[nstar0].n-1]*stararray0[nstar0-1].cc[i]/stararray0[nstar0-1].cc[stararray0[nstar0-1].n-1];
      stararray0[nstar0].cf[i] = stararray0[nstar0].cf[stararray0[nstar0].n-1]*stararray0[nstar0-1].cf[i]/stararray0[nstar0-1].cf[stararray0[nstar0-1].n-1];
      stararray0[nstar0].omega[i] = stararray0[nstar0].omega[stararray0[nstar0].n-1]*stararray0[nstar0-1].omega[i]/stararray0[nstar0-1].omega[stararray0[nstar0-1].n-1];
    }
    //set values for t=0
    stararray0[nstar0].m[0] = stararray0[nstar0].m[1];
    stararray0[nstar0].t[0] = 0.0;
    stararray0[nstar0].r[0] = stararray0[nstar0].r[1];
    stararray0[nstar0].cm[0] = 0.0;
    stararray0[nstar0].llum[0] = stararray0[nstar0].llum[1];
    stararray0[nstar0].lteff[0] = stararray0[nstar0].lteff[1];
    stararray0[nstar0].lambda[0] = stararray0[nstar0].lambda[1];
    stararray0[nstar0].cc[0] = 0.0;
    stararray0[nstar0].cf[0] = stararray0[nstar0].cf[1];
    stararray0[nstar0].omega[0] = stararray0[nstar0].omega[1];
    //increase number of star tables
    nstar0++;
  }*/
  if (kind==-1){
    nstar1 = nstar0;
    stararray1 = stararray0;
    nstar = nstar0;
    stararray = stararray0;
/*    nstar0 = 0;
    stararray0 = NULL;*/
  }else if (kind==-2){
    nstar2 = nstar0;
    stararray2 = stararray0;
/*    nstar0 = 0;
    stararray0 = NULL;*/
  }else if ((20<kind)&&(kind<30)){
    nhe1 = nstar0;
    hestararray1 = stararray0;
/*    nstar0 = 0;
    stararray0 = NULL;*/
  }else{
    nstar = nstar0;
    stararray = stararray0;
    nstar1 = nstar0;
    stararray1 = stararray0;
    nstar2 = nstar0;
    stararray2 = stararray0;
  }
}

void readhestararray(){
  ifstream infile;	//variable for read table from file
  ifstream findfile;	//variable for found files
  char text[10000];	//some text
  int err;		//errorcode
  int i;		//loop variable
  t_He tempHe;	//temporary variable for sort
//  double x1,x2,y1,y2,a,b;	//auxiliary variables for extrapolation

  // lists all files in tables/star in Find.star
  err = sprintf(text,"find %she/ -type f > %sFind.he",tablelocation,tablelocation);
  if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
  err = system(text);
  if (err!=0) cerr << "#Error: can't find tables in he, code: " << err << endl;

  //reserve memory
  hestararray = (t_He *)malloc(sizeof(t_He));
  //check if memory allocation fails
  if (hestararray==NULL) cerr << "#Error: memory allocation failed: hestararray" << endl;
  //initialize number of read star tracks
  nhe = 0;

  //open file with names of all star tables
  err = sprintf(text,"%sFind.he",tablelocation);
  if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
  findfile.open(text);
  do{
    //get name of next star table
    findfile.getline(text,1000);
    if (findfile.fail()){
      break;
    }else{
//      cout << text << endl;
      //open he table
      infile.open(text);
      //enlarge reserved memory
      hestararray = (t_He *)realloc(hestararray,(nhe+1)*sizeof(t_He));
      //check if memory allocation fails
      if (hestararray==NULL) cerr << "#Error: memory reallocation failed: hestararray" << endl;
      //reserve memory
      hestararray[nhe].m = (double *)malloc(sizeof(double));
      hestararray[nhe].t = (double *)malloc(sizeof(double));
      hestararray[nhe].lteff = (double *)malloc(sizeof(double));
      hestararray[nhe].llum = (double *)malloc(sizeof(double));
      hestararray[nhe].lr = (double *)malloc(sizeof(double));
      hestararray[nhe].mw = (double *)malloc(sizeof(double));
      hestararray[nhe].cm = (double *)malloc(sizeof(double));
      hestararray[nhe].r = (double *)malloc(sizeof(double));
      //check if memory allocation fails
      if (hestararray[nhe].m==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].m" << endl;
      if (hestararray[nhe].t==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].t" << endl;
      if (hestararray[nhe].lteff==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].lteff" << endl;
      if (hestararray[nhe].llum==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].llum" << endl;
      if (hestararray[nhe].lr==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].lr" << endl;
      if (hestararray[nhe].mw==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].mw" << endl;
      if (hestararray[nhe].cm==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].cm" << endl;
      if (hestararray[nhe].r==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].r" << endl;
      //initialize values
      hestararray[nhe].n = 0;
/*      hestararray[nhe].m[0] = -1000.;
      hestararray[nhe].t[0] = -1000.;
      hestararray[nhe].lteff[0] = -1000.;
      hestararray[nhe].llum[0] = -1000.;
      hestararray[nhe].lr[0] = -1000.;
      hestararray[nhe].mw[0] = -1000.;
      hestararray[nhe].cm[0] = -1000.;
      hestararray[nhe].r[0] = -1000.;*/
      do{
        //enlarge reserved memory
        hestararray[nhe].m = (double *)realloc(hestararray[nhe].m,(hestararray[nhe].n+1)*sizeof(double));
        hestararray[nhe].t = (double *)realloc(hestararray[nhe].t,(hestararray[nhe].n+1)*sizeof(double));
        hestararray[nhe].lteff = (double *)realloc(hestararray[nhe].lteff,(hestararray[nhe].n+1)*sizeof(double));
        hestararray[nhe].llum = (double *)realloc(hestararray[nhe].llum,(hestararray[nhe].n+1)*sizeof(double));
        hestararray[nhe].lr = (double *)realloc(hestararray[nhe].lr,(hestararray[nhe].n+1)*sizeof(double));
        hestararray[nhe].mw = (double *)realloc(hestararray[nhe].mw,(hestararray[nhe].n+1)*sizeof(double));
        hestararray[nhe].cm = (double *)realloc(hestararray[nhe].cm,(hestararray[nhe].n+1)*sizeof(double));
        hestararray[nhe].r = (double *)realloc(hestararray[nhe].r,(hestararray[nhe].n+1)*sizeof(double));
        //check if memory allocation fails
        if (hestararray[nhe].m==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].m" << endl;
        if (hestararray[nhe].t==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].t" << endl;
        if (hestararray[nhe].lteff==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].lteff" << endl;
        if (hestararray[nhe].llum==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].llum" << endl;
        if (hestararray[nhe].lr==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].lr" << endl;
        if (hestararray[nhe].mw==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].mw" << endl;
        if (hestararray[nhe].cm==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].cm" << endl;
        if (hestararray[nhe].r==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].r" << endl;
        //read in new values
        infile >> hestararray[nhe].m[hestararray[nhe].n];
        infile >> hestararray[nhe].t[hestararray[nhe].n];
        infile >> hestararray[nhe].lteff[hestararray[nhe].n];
        infile >> hestararray[nhe].llum[hestararray[nhe].n];
        infile >> hestararray[nhe].lr[hestararray[nhe].n];
        infile >> hestararray[nhe].mw[hestararray[nhe].n];
        infile >> hestararray[nhe].cm[hestararray[nhe].n];
        infile.ignore(10000,'\n');	//ignore rest of the line(all the characters until '\n')
        if (infile.fail()){	//if an error while last reading occured, e.g. end of file is reached
          //remove last entry
          hestararray[nhe].m = (double *)realloc(hestararray[nhe].m,(hestararray[nhe].n)*sizeof(double));
          hestararray[nhe].t = (double *)realloc(hestararray[nhe].t,(hestararray[nhe].n)*sizeof(double));
          hestararray[nhe].lteff = (double *)realloc(hestararray[nhe].lteff,(hestararray[nhe].n)*sizeof(double));
          hestararray[nhe].llum = (double *)realloc(hestararray[nhe].llum,(hestararray[nhe].n)*sizeof(double));
          hestararray[nhe].lr = (double *)realloc(hestararray[nhe].lr,(hestararray[nhe].n)*sizeof(double));
          hestararray[nhe].mw = (double *)realloc(hestararray[nhe].mw,(hestararray[nhe].n)*sizeof(double));
          hestararray[nhe].cm = (double *)realloc(hestararray[nhe].cm,(hestararray[nhe].n)*sizeof(double));
          hestararray[nhe].r = (double *)realloc(hestararray[nhe].r,(hestararray[nhe].n)*sizeof(double));
          //check if memory allocation fails
          if (hestararray[nhe].m==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].m" << endl;
          if (hestararray[nhe].t==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].t" << endl;
          if (hestararray[nhe].lteff==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].lteff" << endl;
          if (hestararray[nhe].llum==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].llum" << endl;
          if (hestararray[nhe].lr==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].lr" << endl;
          if (hestararray[nhe].mw==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].mw" << endl;
          if (hestararray[nhe].cm==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].cm" << endl;
          if (hestararray[nhe].r==NULL) cerr << "#Error: memory reallocation failed: hestararray[nhe].r" << endl;
          break;
        }else{
          //calculate radius in Rsun
          hestararray[nhe].r[hestararray[nhe].n] = pow(10.0,hestararray[nhe].lr[hestararray[nhe].n]);
        }
/*        //check for maximum values
        if (hestararray[nhe].m[hestararray[nhe].n]>hestararray[nhe].m[0]){
          hestararray[nhe].m[0] = hestararray[nhe].m[hestararray[nhe].n];
        }
        if (hestararray[nhe].t[hestararray[nhe].n]>hestararray[nhe].t[0]){
          hestararray[nhe].t[0] = hestararray[nhe].t[hestararray[nhe].n];
        }
        if (hestararray[nhe].lteff[hestararray[nhe].n]>hestararray[nhe].lteff[0]){
          hestararray[nhe].lteff[0] = hestararray[nhe].lteff[hestararray[nhe].n];
        }
        if (hestararray[nhe].llum[hestararray[nhe].n]>hestararray[nhe].llum[0]){
          hestararray[nhe].llum[0] = hestararray[nhe].llum[hestararray[nhe].n];
        }
        if (hestararray[nhe].lr[hestararray[nhe].n]>hestararray[nhe].lr[0]){
          hestararray[nhe].lr[0] = hestararray[nhe].lr[hestararray[nhe].n];
        }
        if (hestararray[nhe].mw[hestararray[nhe].n]>hestararray[nhe].mw[0]){
          hestararray[nhe].mw[0] = hestararray[nhe].mw[hestararray[nhe].n];
        }
        if (hestararray[nhe].cm[hestararray[nhe].n]>hestararray[nhe].cm[0]){
          hestararray[nhe].cm[0] = hestararray[nhe].cm[hestararray[nhe].n];
        }
        if (hestararray[nhe].r[hestararray[nhe].n]>hestararray[nhe].r[0]){
          hestararray[nhe].r[0] = hestararray[nhe].r[hestararray[nhe].n];
        }*/
        //increase number of array entries
        hestararray[nhe].n++;
      }while (!infile.fail());
      //close he table
      infile.close();
      for (i=nhe;i>0;i--){	//sort hestararray from low mass to high mass
        if (hestararray[i].m[1]<hestararray[i-1].m[1]){	//change order if new mass if lower than the one before
          tempHe = hestararray[i];
          hestararray[i] = hestararray[i-1];
          hestararray[i-1] = tempHe;
//          cout << "change: " << i << " with " << i-1 << endl;
        }else{
          break;
        }
      }
      smoothtable(hestararray[nhe].t, hestararray[nhe].n,1);
      smoothtable(hestararray[nhe].llum, hestararray[nhe].n,0);
      smoothtable(hestararray[nhe].lteff, hestararray[nhe].n,0);
      smoothtable(hestararray[nhe].r, hestararray[nhe].n,0);
      smoothtable(hestararray[nhe].mw, hestararray[nhe].n,-1);
      smoothtable(hestararray[nhe].cm, hestararray[nhe].n,1);
/*      if (nhe>0){
        if (hestararray[nhe].n!=hestararray[0].n) cerr << "#Warning: naked helium tables have different length" << endl;
      }*/
      //increase number of naked helium tables
      nhe++;
    }
  }while (!findfile.fail());
  //close find-file
  findfile.close();
/*  while (hestararray[nhe-1].mw[1]<500.0){ //create extrapolation grid
    //enlarge reserved memory
    hestararray = (t_HRD *)realloc(hestararray,(nhe+1)*sizeof(t_HRD));
    //check if memory allocation fails
    if (hestararray==NULL) cerr << "#Error: memory reallocation failed: hestararray" << endl;
    //take lenght from previous table
    hestararray[nhe].n = hestararray[nhe-1].n;
    //reserve memory
    hestararray[nhe].mw = (double *)malloc(hestararray[nhe].n*sizeof(double));
    hestararray[nhe].t = (double *)malloc(hestararray[nhe].n*sizeof(double));
    hestararray[nhe].r = (double *)malloc(hestararray[nhe].n*sizeof(double));
    hestararray[nhe].cm = (double *)malloc(hestararray[nhe].n*sizeof(double));
    hestararray[nhe].llum = (double *)malloc(hestararray[nhe].n*sizeof(double));
    hestararray[nhe].lteff = (double *)malloc(hestararray[nhe].n*sizeof(double));
    //check if memory allocation fails
    if (hestararray[nhe].mw==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].mw" << endl;
    if (hestararray[nhe].t==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].t" << endl;
    if (hestararray[nhe].r==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].r" << endl;
    if (hestararray[nhe].cm==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].cm" << endl;
    if (hestararray[nhe].llum==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].llum" << endl;
    if (hestararray[nhe].lteff==NULL) cerr << "#Error: memory allocation failed: hestararray[nhe].lteff" << endl;
    //extrapolate values of the last point
    x1 = hestararray[nhe-2].mw[1];
    x2 = hestararray[nhe-1].mw[1];
    hestararray[nhe].mw[hestararray[nhe].n-1] = 2.0*hestararray[nhe-1].mw[hestararray[nhe-1].n-1];
    y1 = hestararray[nhe-2].t[hestararray[nhe-2].n-1]-2.5e+6;	//stay above 2.5Myr
    y2 = hestararray[nhe-1].t[hestararray[nhe-1].n-1]-2.5e+6;	//stay above 2.5Myr
    a = exp((x2*log(y1)-x1*log(y2))/(x2-x1));
    b = log(y1/y2)/(x2-x1);
    hestararray[nhe].t[hestararray[nhe].n-1] = 2.5e+6+a*exp(-b*2.0*x2);	//stay above 2.5Myr
    y1 = hestararray[nhe-2].r[1]-Schwarzschildradius(x1);	//stay above Schwarzschildradius
    y2 = hestararray[nhe-1].r[1]-Schwarzschildradius(x2);	//stay above Schwarzschildradius
    a = exp((x2*log(y1)-x1*log(y2))/(x2-x1));
    b = log(y1/y2)/(x2-x1);
    hestararray[nhe].r[hestararray[nhe].n-1] = (Schwarzschildradius(2.0*x2)+a*exp(-b*2.0*x2))*hestararray[nhe-1].r[hestararray[nhe-1].n-1]/hestararray[nhe-1].r[1];	//stay above Schwarzschildradius
    hestararray[nhe].cm[hestararray[nhe].n-1] = 2*hestararray[nhe-1].cm[hestararray[nhe-1].n-1];
    hestararray[nhe].llum[hestararray[nhe].n-1] = hestararray[nhe-2].llum[hestararray[nhe-2].n-1]+(2.0*x2-x1)/(x2-x1)*(hestararray[nhe-1].llum[hestararray[nhe-1].n-1]-hestararray[nhe-2].llum[hestararray[nhe-2].n-1]);
    hestararray[nhe].lteff[hestararray[nhe].n-1] = hestararray[nhe-2].lteff[hestararray[nhe-2].n-1]+(2.0*x2-x1)/(x2-x1)*(hestararray[nhe-1].lteff[hestararray[nhe-1].n-1]-hestararray[nhe-2].lteff[hestararray[nhe-2].n-1]);
    for (i=1;i<hestararray[nhe].n;i++){ //calculate extrapolation by keeping the ratio to last point as for the previous table
      hestararray[nhe].mw[i] = hestararray[nhe].mw[hestararray[nhe].n-1]*hestararray[nhe-1].mw[i]/hestararray[nhe-1].mw[hestararray[nhe-1].n-1];
      hestararray[nhe].t[i] = hestararray[nhe].t[hestararray[nhe].n-1]*hestararray[nhe-1].t[i]/hestararray[nhe-1].t[hestararray[nhe-1].n-1];
      hestararray[nhe].r[i] = hestararray[nhe].r[hestararray[nhe].n-1]*hestararray[nhe-1].r[i]/hestararray[nhe-1].r[hestararray[nhe-1].n-1];
      hestararray[nhe].cm[i] = hestararray[nhe].cm[hestararray[nhe].n-1]*hestararray[nhe-1].cm[i]/hestararray[nhe-1].cm[hestararray[nhe-1].n-1];
      hestararray[nhe].llum[i] = hestararray[nhe].llum[hestararray[nhe].n-1]*hestararray[nhe-1].llum[i]/hestararray[nhe-1].llum[hestararray[nhe-1].n-1];
      hestararray[nhe].lteff[i] = hestararray[nhe].lteff[hestararray[nhe].n-1]*hestararray[nhe-1].lteff[i]/hestararray[nhe-1].lteff[hestararray[nhe-1].n-1];
    }
    //set values for t=0
    hestararray[nhe].mw[0] = hestararray[nhe].mw[1];
    hestararray[nhe].t[0] = 0.0;
    hestararray[nhe].r[0] = hestararray[nhe].r[1];
    hestararray[nhe].cm[0] = 0.0;
    hestararray[nhe].llum[0] = hestararray[nhe].llum[1];
    hestararray[nhe].lteff[0] = hestararray[nhe].lteff[1];
    //increase number of star tables
    nhe++;
  }*/
}

void readlambdaarray(){
  ifstream infile;	//variable for read table from file
  char text[10000];	//some text
  int err;		//errorcode

  //open lambda table
  err = sprintf(text,"%slambda/ext.lambda",tablelocation);
  if (err<=0) cerr << "#Error: can't write to text, code: " << err << endl;
  infile.open(text);
  //reserve memory
  lambdaarray.m = (double *)malloc(sizeof(double));
  lambdaarray.r = (double *)malloc(sizeof(double));
  lambdaarray.lambda = (double *)malloc(sizeof(double));
  //check if memory allocation fails
  if (lambdaarray.m==NULL) cerr << "#Error: memory allocation failed: lambdaarray.m" << endl;
  if (lambdaarray.r==NULL) cerr << "#Error: memory allocation failed: lambdaarray.r" << endl;
  if (lambdaarray.lambda==NULL) cerr << "#Error: memory allocation failed: lambdaarray.lambda" << endl;
  //initialize values
  lambdaarray.n = 0;
/*  lambdaarray.m[0] = -1000.;
  lambdaarray.r[0] = -1000.;
  lambdaarray.lambda[0] = -1000.;*/
  do{
    //enlarge reserved memory
    lambdaarray.m = (double *)realloc(lambdaarray.m,(lambdaarray.n+1)*sizeof(double));
    lambdaarray.r = (double *)realloc(lambdaarray.r,(lambdaarray.n+1)*sizeof(double));
    lambdaarray.lambda = (double *)realloc(lambdaarray.lambda,(lambdaarray.n+1)*sizeof(double));
    //check if memory allocation fails
    if (lambdaarray.m==NULL) cerr << "#Error: memory reallocation failed: lambdaarray.m" << endl;
    if (lambdaarray.r==NULL) cerr << "#Error: memory reallocation failed: lambdaarray.r" << endl;
    if (lambdaarray.lambda==NULL) cerr << "#Error: memory reallocation failed: lambdaarray.lambda" << endl;
    //read in new values
    infile >> lambdaarray.m[lambdaarray.n];
    infile >> lambdaarray.r[lambdaarray.n];
    infile >> lambdaarray.lambda[lambdaarray.n];
    infile.ignore(10000,'\n');	//ignore rest of the line(all the characters until '\n')
/*    //check for maximum values
    if (lambdaarray.m[lambdaarray.n]>lambdaarray.m[0]){
      lambdaarray.m[0] = lambdaarray.m[lambdaarray.n];
    }
    if (lambdaarray.r[lambdaarray.n]>lambdaarray.r[0]){
      lambdaarray.r[0] = lambdaarray.r[lambdaarray.n];
    }
    if (lambdaarray.lambda[lambdaarray.n]>lambdaarray.lambda[0]){
      lambdaarray.lambda[0] = lambdaarray.lambda[lambdaarray.n];
    }*/
    //increase number of array entries
    lambdaarray.n++;
  }while (!infile.fail());
  //remove last entry
  lambdaarray.n--;
  lambdaarray.m = (double *)realloc(lambdaarray.m,(lambdaarray.n)*sizeof(double));
  lambdaarray.r = (double *)realloc(lambdaarray.r,(lambdaarray.n)*sizeof(double));
  lambdaarray.lambda = (double *)realloc(lambdaarray.lambda,(lambdaarray.n)*sizeof(double));
  //check if memory allocation fails
  if (lambdaarray.m==NULL) cerr << "#Error: memory reallocation failed: lambdaarray.m" << endl;
  if (lambdaarray.r==NULL) cerr << "#Error: memory reallocation failed: lambdaarray.r" << endl;
  if (lambdaarray.lambda==NULL) cerr << "#Error: memory reallocation failed: lambdaarray.lambda" << endl;
  //close lambda table
  infile.close();
}

void readtables(double metal, bool rotation){
  ifstream infile;	//variable for read table location from file

  infile.open("tablelocation.txt");	//open file with the name "tablelocation.txt"
  if (infile.is_open()) infile >> tablelocation;	//get location of the tables
  infile.close();		//close the input file

  destroytables();	//free old memory of all tables for initializing them new
/*  readstararray(-99);	//initializes and read in star tables*/
/*  readstararray(-1);	//initializes and read in Ms-star tables
  readstararray(-2);	//initializes and read in Core-Helium-burning-star tables*/
/*  stararray = stararray1;
  nstar = nstar1;
  stararray1 = NULL;
  nstar1 = 0;*/
//  readhestararray();	//initializes and read in hestar tables
  if (metal>0.007){	//Milky Way metallicity(~0.009)
    if (rotation) readstararray(11);	//initializes and read in MWrotstars tables
    else readstararray(1);	//initializes and read in MWstars tables
    readstararray(21);	//initializes and read in MWhestars tables
  }else if(metal>0.0035){	//LMC metallicity(~0.005)
    readstararray(2);	//initializes and read in LMCstars tables
    readstararray(22);	//initializes and read in LMChestars tables
  }else if(metal>0.001){	//SMC metallicity(~0.002)
    readstararray(3);	//initializes and read in SMCstars tables
    readstararray(23);	//initializes and read in SMChestars tables
  }else{	//IZw18 metallicity(~0.0002)
    readstararray(4);	//initializes and read in IZw18stars tables
    readstararray(24);	//initializes and read in IZw18hestars tables
  }
  readlambdaarray();	//initializes and read in lambda tables
}

void destroytables(){	//delete arrays and free their memory of all tables
  int i;

  //stararray
  if ((nstar>0)&&(stararray!=stararray1)&&(stararray!=stararray2)&&(stararray!=hestararray1)){
//    if (screen) cout << "free stararray[i].*" << endl;
    for (i=0;i<nstar;i++){
      free(stararray[i].m);
      free(stararray[i].t);
      free(stararray[i].r);
      free(stararray[i].cm);
      free(stararray[i].llum);
      free(stararray[i].lteff);
      free(stararray[i].lambda);
      free(stararray[i].cc);
      free(stararray[i].cf);
      free(stararray[i].omega);
    }
//    if (screen) cout << "free stararray" << endl;
    free(stararray);
  }
  nstar = 0;

  //stararray1
  if ((nstar1>0)&&(stararray1!=stararray2)){
//    if (screen) cout << "free stararray1[i].*" << endl;
    for (i=0;i<nstar1;i++){
      free(stararray1[i].m);
      free(stararray1[i].t);
      free(stararray1[i].r);
      free(stararray1[i].cm);
      free(stararray1[i].llum);
      free(stararray1[i].lteff);
      free(stararray1[i].lambda);
      free(stararray1[i].cc);
      free(stararray1[i].cf);
      free(stararray1[i].omega);
    }
//    if (screen) cout << "free stararray1" << endl;
    free(stararray1);
  }
  nstar1 = 0;

  //stararray2
  if (nstar2>0){
//    if (screen) cout << "free stararray2[i].*" << endl;
    for (i=0;i<nstar2;i++){
      free(stararray2[i].m);
      free(stararray2[i].t);
      free(stararray2[i].r);
      free(stararray2[i].cm);
      free(stararray2[i].llum);
      free(stararray2[i].lteff);
      free(stararray2[i].lambda);
      free(stararray2[i].cc);
      free(stararray2[i].cf);
      free(stararray2[i].omega);
    }
//    if (screen) cout << "free stararray2" << endl;
    free(stararray2);
  }
  nstar2 = 0;

  //hestararray1
  if (nhe1>0){
//    if (screen) cout << "free hestararray1[i].*" << endl;
    for (i=0;i<nhe1;i++){
      free(hestararray1[i].m);
      free(hestararray1[i].t);
      free(hestararray1[i].r);
      free(hestararray1[i].cm);
      free(hestararray1[i].llum);
      free(hestararray1[i].lteff);
      free(hestararray1[i].lambda);
      free(hestararray1[i].cc);
      free(hestararray1[i].cf);
      free(hestararray1[i].omega);
    }
//    if (screen) cout << "free hestararray1" << endl;
    free(hestararray1);
  }
  nhe1 = 0;

  //hestararray
  if (nhe>0){
//    if (screen) cout << "free hestararray[i].*" << endl;
    for (i=0;i<nhe;i++){
      free(hestararray[i].m);
      free(hestararray[i].t);
      free(hestararray[i].lteff);
      free(hestararray[i].llum);
      free(hestararray[i].lr);
      free(hestararray[i].mw);
      free(hestararray[i].cm);
      free(hestararray[i].r);
    }
  }
//  if (screen) cout << "free hestararray" << endl;
  free(hestararray);
  nhe = 0;

  //lambdaarray
//  if (screen) cout << "free lambdaarray.*" << endl;
  free(lambdaarray.m);
  free(lambdaarray.r);
  free(lambdaarray.lambda);
  lambdaarray.n = 0;
}

void smoothtable(double* table, int n, int monotonic){	//smooth the table to avoid steps because of the precision in input data
  int i,j,k;		//loop variables
  double diff=0.0;	//difference to next value
  double lastvalue=0.0;	//last value

  if (debug) cerr << "table=(";
  for (i=0;i<n-1;i++){
    if (debug) cerr << table[i] << ", ";
    j = 1;
    while (table[i]==table[i+j]){
      j++;
      if (i+j>=n) break;
    }
    if (j>1){
      if (monotonic==0){
        if (i+j<n) monotonic = table[i+j]-table[i];
        else if (i>0) monotonic = table[i]-table[i-1];
      }
      if (table[i]!=0) lastvalue = table[i];
      else if (i+j<n) lastvalue = table[i+j];
      else lastvalue = 0.0;
      if (monotonic<0) diff = -abs(lastvalue)*accuracy;
      else diff = abs(lastvalue)*accuracy;
      for (k=1;k<j;k++){
        table[i+k] += diff*((double)k)/((double)j);
        if (debug && (i+k<n-1)) cerr << table[i+k] << "[" << table[i+k]-table[i+k-1] << "=" << diff*((double)k)/((double)j) << "=" << diff << "*" << k << "/" << j << "], ";
      }
    }
    i += j-1;
  }
  if (debug) cerr << table[n-1] << ")" << endl << " last diff=" << diff << endl;
}
