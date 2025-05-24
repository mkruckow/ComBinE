//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
#include "ComBinElib.h" // local header file

//using namespace std;

/*saves for separate evolution*/
t_HRD savetrack;	//HRD track of a star
int savermax;		//rmax of a star track

void initsystem(t_system& system, double mp, double q, double a, double e, double t, double metal, double vp, double vs, double x, double y, double z, double vx, double vy, double vz, double rhostar){
  int i;	//index variable

  system.n = 1;
  //set up the two components
  if (mp<accuracy){
    initnonstar(system.prim, fmax(-mp,accuracy), t, metal, vp, 6);
  }else{
    initstar(system.prim, mp, t, metal, vp);
  }
  if (q<accuracy){
    initnonstar(system.sec, fmax(-q,accuracy), t, metal, vs, 6);
  }else{
    initstar(system.sec, mp*q, t, metal, vs);
  }
  if (separateevolution){
    savetrack = system.sec.track;
    savermax = system.sec.rmax;
    system.sec.track.n = 3;
    //remove not needed array part
    system.sec.track.m = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.t = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.r = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.cm = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.llum = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.lteff = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.lambda = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.cc = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.cf = (double *)malloc(system.sec.track.n*sizeof(double));
    system.sec.track.omega = (double *)malloc(system.sec.track.n*sizeof(double));
    //check if memory allocation fails
    if (system.sec.track.m==NULL) cerr << "#Error: memory allocation failed: system.sec.track.m" << endl;
    if (system.sec.track.t==NULL) cerr << "#Error: memory allocation failed: system.sec.track.t" << endl;
    if (system.sec.track.r==NULL) cerr << "#Error: memory allocation failed: system.sec.track.r" << endl;
    if (system.sec.track.cm==NULL) cerr << "#Error: memory allocation failed: system.sec.track.cm" << endl;
    if (system.sec.track.llum==NULL) cerr << "#Error: memory allocation failed: system.sec.track.llum" << endl;
    if (system.sec.track.lteff==NULL) cerr << "#Error: memory allocation failed: system.sec.track.lteff" << endl;
    if (system.sec.track.lambda==NULL) cerr << "#Error: memory allocation failed: system.sec.track.lambda" << endl;
    if (system.sec.track.cc==NULL) cerr << "#Error: memory allocation failed: system.sec.track.cc" << endl;
    if (system.sec.track.cf==NULL) cerr << "#Error: memory allocation failed: system.sec.track.cf" << endl;
    if (system.sec.track.omega==NULL) cerr << "#Error: memory allocation failed: system.sec.track.omega" << endl;
    //copy values
    system.sec.track.m[0] = savetrack.m[1];
    system.sec.track.t[0] = savetrack.t[1];
    system.sec.track.r[0] = savetrack.r[1];
    system.sec.track.cm[0] = savetrack.cm[1];
    system.sec.track.llum[0] = savetrack.llum[1];
    system.sec.track.lteff[0] = savetrack.lteff[1];
    system.sec.track.lambda[0] = savetrack.lambda[1];
    system.sec.track.cc[0] = savetrack.cc[1];
    system.sec.track.cf[0] = savetrack.cf[1];
    system.sec.track.omega[0] = savetrack.omega[1];
    //value check
    cout << "memory check: &savetrack=" << &savetrack << " &system.sec.track=" << &system.sec.track << endl;
    cout << "&savetrack.(n,m,t,r,cm,llum,lteff,cc,cf,omega)=\t\t" << &savetrack.n << ",\t" << savetrack.m << ",\t" << savetrack.t << ",\t" << savetrack.r << ",\t" << savetrack.cm << ",\t" << savetrack.llum << ",\t" << savetrack.lteff << ",\t" << savetrack.cc << ",\t" << savetrack.cf << ",\t" << savetrack.omega << "\n&system.sec.track.(n,m,t,r,cm,llum,lteff,lambda,cc,cf,omega)=\t" << &system.sec.track.n << ",\t" << system.sec.track.m << ",\t" << system.sec.track.t << ",\t" << system.sec.track.r << ",\t" << system.sec.track.cm << ",\t" << system.sec.track.llum << ",\t" << system.sec.track.lteff << ",\t" << system.sec.track.lambda << ",\t" << system.sec.track.cc << ",\t" << system.sec.track.cf << ",\t" << system.sec.track.omega << endl;
    cout << "savetrack.n=" << savetrack.n << " savetrack.(m,t,r,cm,llum,lteff,lambda,cc,cf,omega)=\n[(" << savetrack.m[0] << "," << savetrack.t[0] << "," << savetrack.r[0] << "," << savetrack.cm[0] << "," << savetrack.llum[0] << "," << savetrack.lteff[0] << "," << savetrack.lambda[0] << "," << savetrack.cc[0] << "," << savetrack.cf[0] << "," << savetrack.omega[0] << ")";
    for (i=1;i<savetrack.n;i++){
      cout << ", (" << savetrack.m[i] << "," << savetrack.t[i] << "," << savetrack.r[i] << "," << savetrack.cm[i] << "," << savetrack.llum[i] << "," << savetrack.lteff[i] << "," << savetrack.lambda[i] << "," << savetrack.cc[i] << "," << savetrack.cf[i] << "," << savetrack.omega[i] << ")";
      if (i<system.sec.track.n){	//copy values
        system.sec.track.m[i] = savetrack.m[i];
        system.sec.track.t[i] = savetrack.t[i];
        system.sec.track.r[i] = savetrack.r[i];
        system.sec.track.cm[i] = savetrack.cm[i];
        system.sec.track.llum[i] = savetrack.llum[i];
        system.sec.track.lteff[i] = savetrack.lteff[i];
        system.sec.track.lambda[i] = savetrack.lambda[i];
        system.sec.track.cc[i] = savetrack.cc[i];
        system.sec.track.cf[i] = savetrack.cf[i];
        system.sec.track.omega[i] = savetrack.omega[i];
      }
    }
    system.sec.track.t[2] = 9.0e+99;	//new value
    system.sec.rmax = 1;	//new value
    cout << "]" << endl << "system.sec.track.n=" << system.sec.track.n << " system.sec.track.(m,t,r,cm,llum,lteff,lambda)=\n[(" << system.sec.track.m[0] << "," << system.sec.track.t[0] << "," << system.sec.track.r[0] << "," << system.sec.track.cm[0] << "," << system.sec.track.llum[0] << "," << system.sec.track.lteff[0] << "," << system.sec.track.lambda[0] << "," << system.sec.track.cc[0] << "," << system.sec.track.cf[0] << "," << system.sec.track.omega[0] << ")";
    for (i=1;i<system.sec.track.n;i++){
      cout << ", (" << system.sec.track.m[i] << "," << system.sec.track.t[i] << "," << system.sec.track.r[i] << "," << system.sec.track.cm[i] << "," << system.sec.track.llum[i] << "," << system.sec.track.lteff[i] << "," << system.sec.track.lambda[i] << "," << system.sec.track.cc[i] << "," << system.sec.track.cf[i] << "," << system.sec.track.omega[i] << ")";
    }
    cout << "]" << endl;
  }
  //determine system mass in Msun
  system.M = (double *)malloc(sizeof(double));
  if (system.M==NULL) cerr << "#Error: memory allocation failed: system.M" << endl;	//check if memory allocation fails
  system.M[0] = system.prim.m[0]+system.sec.m[0];
  //set mass ratios
  system.qp = (double *)malloc(sizeof(double));
  if (system.qp==NULL) cerr << "#Error: memory allocation failed: system.qp" << endl;	//check if memory allocation fails
  system.qp[0] = 1.0/q;
  system.qs = (double *)malloc(sizeof(double));
  if (system.qs==NULL) cerr << "#Error: memory allocation failed: system.qs" << endl;	//check if memory allocation fails
  system.qs[0] = q;
  //set semi-major-axis in Rsun
  system.a = (double *)malloc(sizeof(double));
  if (system.a==NULL) cerr << "#Error: memory allocation failed: system.a" << endl;	//check if memory allocation fails
  system.a[0] = a;
  //calculate period in yr
  system.P = (double *)malloc(sizeof(double));
  if (system.P==NULL) cerr << "#Error: memory allocation failed: system.P" << endl;	//check if memory allocation fails
  system.P[0] = 2.0*M_PI*sqrt(cubic(system.a[0])/(G*system.M[0]));
  //determine the Roche-lobes in Rsun
  system.rp = (double *)malloc(sizeof(double));
  if (system.rp==NULL) cerr << "#Error: memory allocation failed: system.rp" << endl;	//check if memory allocation fails
  system.rp[0] = rocherad(system.a[0],system.qp[0]);
  system.rs = (double *)malloc(sizeof(double));
  if (system.rs==NULL) cerr << "#Error: memory allocation failed: system.rs" << endl;	//check if memory allocation fails
  system.rs[0] = rocherad(system.a[0],system.qs[0]);
  //set the eccentricity
  system.e = (double *)malloc(sizeof(double));
  if (system.e==NULL) cerr << "#Error: memory allocation failed: system.e" << endl;	//check if memory allocation fails
  system.e[0] = e;
  //set the peri-astron
  system.peri = (double *)malloc(sizeof(double));
  if (system.peri==NULL) cerr << "#Error: memory allocation failed: system.peri" << endl;	//check if memory allocation fails
  system.peri[0] = system.a[0]*(1.0-system.e[0]);
  //set the galactic x-position
  system.x = (double *)malloc(sizeof(double));
  if (system.x==NULL) cerr << "#Error: memory allocation failed: system.x" << endl;	//check if memory allocation fails
  system.x[0] = x;
  //set the galactic y-position
  system.y = (double *)malloc(sizeof(double));
  if (system.y==NULL) cerr << "#Error: memory allocation failed: system.y" << endl;	//check if memory allocation fails
  system.y[0] = y;
  //set the galactic z-position
  system.z = (double *)malloc(sizeof(double));
  if (system.z==NULL) cerr << "#Error: memory allocation failed: system.z" << endl;	//check if memory allocation fails
  system.z[0] = z;
  //set the galactic x-velocity
  system.vx = (double *)malloc(sizeof(double));
  if (system.vx==NULL) cerr << "#Error: memory allocation failed: system.vx" << endl;	//check if memory allocation fails
  system.vx[0] = vx;
  //set the galactic y-velocity
  system.vy = (double *)malloc(sizeof(double));
  if (system.vy==NULL) cerr << "#Error: memory allocation failed: system.vy" << endl;	//check if memory allocation fails
  system.vy[0] = vy;
  //set the galactic z-velocity
  system.vz = (double *)malloc(sizeof(double));
  if (system.vz==NULL) cerr << "#Error: memory allocation failed: system.vz" << endl;	//check if memory allocation fails
  system.vz[0] = vz;
  //set the density of surrounding stars
  system.rhostar = (double *)malloc(sizeof(double));
  if (system.rhostar==NULL) cerr << "#Error: memory allocation failed: system.rhostar" << endl;	//check if memory allocation fails
  system.rhostar[0] = rhostar;
  //set the phase index (for details see ComBinElib.h)
  system.phase = (int *)malloc(sizeof(int));
  if (system.phase==NULL) cerr << "#Error: memory allocation failed: system.phase" << endl;	//check if memory allocation fails
  if ((system.prim.stage[0]==6)&&(system.sec.stage[0]==6)) system.phase[0] = -100;	//no stars
  else if ((system.prim.stage[0]==6)||(system.sec.stage[0]==6)) system.phase[0] = -1;	//start with single star wind
  else system.phase[0] = 0;	//start with binary wind phase
  system.stagechange = 0;
  system.formation = 0;
  system.tgw = 0.0;
}

void initstar(t_star& star, double mass, double time, double metal, double vrot){
  int iup,ilow,j;	//index variables	//,jlow,jup
/*  int nlow=0,nup=0;	//index variables: number of points in ilow/iup=mass index*/
  double mratio;	//ratio between the two nearest evolutionary tracks
  double tratio;	//ratio between the two nearest times
/*  double rmax;	//maximal radius in Rsun
  double ratiolow,ratioup;	//mass ratio, ratios at lower/upper mass track
  double m,t,r,cm,llum,lteff,lambda;	//interpolated values of mass, age, radius, core mass, log10(luminosity), log10(effective temperature), lambda*/
  double EbindperG=0.0;	//envelope binding energy for interpolation

//  if (screen) cout << "Create star: mass(m)=" << mass << "Msun age(t)="<<time << "yr" << endl;
  //create evolutionary track of this star
  for (iup=1;iup<nstar1-1;iup++){	//find mass index
    if(stararray1[iup].m[0]>mass){	//stararray must be sorted from low mass to high mass (done when table is read in)
      break;
    }
  }
//  if (screen) cout << "iup=" << iup << " n_star=" << nstar << " stararray[iup].m[1]=" << stararray[iup].m[1] << "Msun" << endl;
  ilow = iup-1;
  //calculate the ratio between the two nearest evolutionary tracks
  mratio = ratio(mass,stararray1[iup].m[0],stararray1[ilow].m[0]);
//  if (screen) cout << "m_ratio=" << mratio << endl;
  star.last = 1;
  if (debug) cerr << "mass=" << mass << "Msun stararray1[" << ilow << "].n=" << stararray1[ilow].n << " stararray1[" << iup << "].n=" << stararray1[iup].n << " mratio=" << mratio << endl;

  getnewtrack(stararray1, ilow, iup, mratio, -1, star.track, star.rmax);
  //star.track.omega[1] = vrot/star.track.r[1];
  //find time index in the evolution
  for (j=star.last+1;j<star.track.n-1;j++){
    if (star.track.t[j]>time){
      break;
    }
  }

  //calculate time ratio
  tratio = ratio(time,star.track.t[j],star.track.t[j-1]);
  //determine stage
  star.stage = (int *)malloc(sizeof(int));
  if (star.stage==NULL) cerr << "#Error: memory allocation failed: star.stage" << endl;	//check if memory allocation fails
  if ((j==star.track.n-1)&&(tratio<=0)){	//star at end of hydrogen burning / start of helium burning
    star.stage[0] = 1;
    tratio = 0.0;
  }else if ((tratio>=0)&&(tratio<=1)){	//star at hydrogen burning
    star.stage[0] = 0;
  }else{
    cerr << "#Error: wrong time span: star.track.t[" << j << "]=" << star.track.t[j] << "yr time=" << time << "yr star.track.t[" << j-1 << "]=" << star.track.t[j-1] << "yr t_ratio=" << tratio << endl;
  }

  //interpolate linearly the mass(Msun)
  star.m = (double *)malloc(sizeof(double));
  if (star.m==NULL) cerr << "#Error: memory allocation failed: star.m" << endl;	//check if memory allocation fails
  star.m[0] = star.track.m[j]-tratio*(star.track.m[j]-star.track.m[j-1]);
  if (star.m[0]<0) cerr << "#Error: m=" << star.m[0] << "Msun t_ratio=" << tratio << " j=" << j << " star.track.m[j]=" << star.track.m[j] << "Msun star.track.m[j-1]=" << star.track.m[j-1] << "Msun" << endl;	//check if mass is negative (unphysical)
  //set time(yr)
  star.t = (double *)malloc(sizeof(double));
  if (star.t==NULL) cerr << "#Error: memory allocation failed: star.t" << endl;	//check if memory allocation fails
  star.t[0] = time;
  //interpolate linearly the radius(Rsun)
  star.r = (double *)malloc(sizeof(double));
  if (star.r==NULL) cerr << "#Error: memory allocation failed: star.r" << endl;	//check if memory allocation fails
  star.r[0] = star.track.r[j]-tratio*(star.track.r[j]-star.track.r[j-1]);
  if (star.r[0]<0) cerr << "#Error: r=" << star.r[0] << "Rsun t_ratio=" << tratio << " j=" << j << " star.track.r[j]=" << star.track.r[j] << "Rsun star.track.r[j-1]=" << star.track.r[j-1] << "Rsun" << endl;	//check if radius is negative (unphysical)
  //interpolate linearly the core mass(Msun)
  star.cm = (double *)malloc(sizeof(double));
  if (star.cm==NULL) cerr << "#Error: memory allocation failed: star.cm" << endl;	//check if memory allocation fails
  star.cm[0] = star.track.cm[j]-tratio*(star.track.cm[j]-star.track.cm[j-1]);
  if (star.cm[0]<0) cerr << "#Error: cm=" << star.cm[0] << "Msun t_ratio=" << tratio << " j=" << j << " star.track.cm[j]=" << star.track.cm[j] << "Msun star.track.cm[j-1]=" << star.track.cm[j-1] << "Msun" << endl;	//check if core mass is negative (unphysical)
  //interpolate linearly the luminosity(lg(lum/Lsun))
  star.llum = (double *)malloc(sizeof(double));
  if (star.llum==NULL) cerr << "#Error: memory allocation failed: star.llum" << endl;	//check if memory allocation fails
  star.llum[0] = star.track.llum[j]-tratio*(star.track.llum[j]-star.track.llum[j-1]);
  //interpolate linearly the effective temperature(lg(teff/K))
  star.lteff = (double *)malloc(sizeof(double));
  if (star.lteff==NULL) cerr << "#Error: memory allocation failed: star.lteff" << endl;	//check if memory allocation fails
  star.lteff[0] = star.track.lteff[j]-tratio*(star.track.lteff[j]-star.track.lteff[j-1]);
  //interpolate linearly the lambda
  star.lambda = (double *)malloc(sizeof(double));
  if (star.lambda==NULL) cerr << "#Error: memory allocation failed: star.lambda" << endl;	//check if memory allocation fails
//  star.lambda[0] = star.track.lambda[j]-tratio*(star.track.lambda[j]-star.track.lambda[j-1]);
  EbindperG = (star.track.m[j]*(star.track.m[j]-star.track.cm[j])/(star.track.lambda[j]*star.track.r[j])
              -tratio*(star.track.m[j]*(star.track.m[j]-star.track.cm[j])/(star.track.lambda[j]*star.track.r[j])
                      -star.track.m[j-1]*(star.track.m[j-1]-star.track.cm[j-1])/(star.track.lambda[j-1]*star.track.r[j-1])));	//interpolate envelope binding energy
  if ((star.m[0]-star.cm[0]==0.0)&&(EbindperG==0.0)) star.lambda[0] = ignorev;	//no interpolation
  else star.lambda[0] = star.m[0]*(star.m[0]-star.cm[0])/(EbindperG*star.r[0]);	//interpolate lambda	// E_bind=G*M*M_env / lambda*R
  //interpolate linearly the carbon core mass(Msun)
  star.cc = (double *)malloc(sizeof(double));
  if (star.cc==NULL) cerr << "#Error: memory allocation failed: star.cc" << endl;	//check if memory allocation fails
  star.cc[0] = star.track.cc[j]-tratio*(star.track.cc[j]-star.track.cc[j-1]);
  if (star.cc[0]<0) cerr << "#Error: cc=" << star.cc[0] << "Msun t_ratio=" << tratio << " j=" << j << " star.track.cc[j]=" << star.track.cc[j] << "Msun star.track.cc[j-1]=" << star.track.cc[j-1] << "Msun" << endl;	//check if core mass is negative (unphysical)
  //interpolate linearly the concentration factor
  star.cf = (double *)malloc(sizeof(double));
  if (star.cf==NULL) cerr << "#Error: memory allocation failed: star.cf" << endl;	//check if memory allocation fails
  star.cf[0] = star.track.cf[j]-tratio*(star.track.cf[j]-star.track.cf[j-1]);
  //set initial value
  star.omega = (double *)malloc(sizeof(double));
  if (star.omega==NULL) cerr << "#Error: memory allocation failed: star.omega" << endl;	//check if memory allocation fails
  star.omega[0] = vrot/star.r[0];	//converting rotational velocity into angular velocity
/*  star.track.omega[j-1] = vrot/star.r[0];	//converting rotational velocity into angular velocity
  star.track.omega[j] = star.track.omega[j-1];
  star.omega[0] = star.track.omega[j]-tratio*(star.track.omega[j]-star.track.omega[j-1]);*/
/*  star.RLO = false;
  star.HeRLO = false;*/
  star.inimass = star.track.inimass;
  star.metal = metal;
  star.SN.inimass = 0.0;
  star.SN.Hecoremass = 0.0;
  star.SN.COcoremass = 0.0;
  star.SN.remnantmass = 0.0;
  star.SN.w = 0.0;
  star.SN.theta = 0.0;
  star.SN.phi = 0.0;
  star.SN.vsys = 0.0;
  star.SN.startype = 0;
  star.SN.type = 0;
}

void initnonstar(t_star& star, double mass, double time, double metal, double vrot, int kind){
  double radius=accuracy, coremass=mass, llum=ignorev, lteff=ignorev, lambda=ignorev, carboncoremass=mass, concentrationfactor=ignorev;

//  if (screen) cout << "Create non star: mass(m)=" << mass << "Msun age(t)=" << time << "yr kind=" << kind << endl;
  //get radius of compact object in Rsun
  if (kind==5){
    radius = Schwarzschildradius(mass);
  }else if (kind==4){
    radius = NSradius(mass);
  }else if (kind==3){
    radius = WDradius(mass);
  }
  //create evolutionary track of this star
  star.last = 1;
  star.rmax = 1;
  //set initial mass of the new track
  star.track.inimass = mass;
  //set length for the arrays (maximal values and first and last track point)
  star.track.n = 3;
  //set TAMS index
  star.track.TAMS = 0;	//no TAMS
  //reserve memory
  star.track.m = (double *)malloc(star.track.n*sizeof(double));
  star.track.t = (double *)malloc(star.track.n*sizeof(double));
  star.track.r = (double *)malloc(star.track.n*sizeof(double));
  star.track.cm = (double *)malloc(star.track.n*sizeof(double));
  star.track.llum = (double *)malloc(star.track.n*sizeof(double));
  star.track.lteff = (double *)malloc(star.track.n*sizeof(double));
  star.track.lambda = (double *)malloc(star.track.n*sizeof(double));
  star.track.cc = (double *)malloc(star.track.n*sizeof(double));
  star.track.cf = (double *)malloc(star.track.n*sizeof(double));
  star.track.omega = (double *)malloc(star.track.n*sizeof(double));
  //check if memory allocation fails
  if (star.track.m==NULL) cerr << "#Error: memory allocation failed: star.track.m" << endl;
  if (star.track.t==NULL) cerr << "#Error: memory allocation failed: star.track.t" << endl;
  if (star.track.r==NULL) cerr << "#Error: memory allocation failed: star.track.r" << endl;
  if (star.track.cm==NULL) cerr << "#Error: memory allocation failed: star.track.cm" << endl;
  if (star.track.llum==NULL) cerr << "#Error: memory allocation failed: star.track.llum" << endl;
  if (star.track.lteff==NULL) cerr << "#Error: memory allocation failed: star.track.lteff" << endl;
  if (star.track.lambda==NULL) cerr << "#Error: memory allocation failed: star.track.lambda" << endl;
  if (star.track.cc==NULL) cerr << "#Error: memory allocation failed: star.track.cc" << endl;
  if (star.track.cf==NULL) cerr << "#Error: memory allocation failed: star.track.cf" << endl;
  if (star.track.omega==NULL) cerr << "#Error: memory allocation failed: star.track.omega" << endl;
  //set new evolutionary track
  //set first track point
  star.track.m[1] = mass;
  star.track.t[1] = time;
  star.track.r[1] = radius;
  star.track.cm[1] = coremass;
  star.track.llum[1] = llum;
  star.track.lteff[1] = lteff;
  star.track.lambda[1] = lambda;
  star.track.cc[1] = carboncoremass;
  star.track.cf[1] = concentrationfactor;
  star.track.omega[1] = vrot/star.track.r[1];
  //set last track point
  star.track.m[2] = mass;
  star.track.t[2] = 1.0e+99;
  star.track.r[2] = radius;
  star.track.cm[2] = coremass;
  star.track.llum[2] = llum;
  star.track.lteff[2] = lteff;
  star.track.lambda[2] = lambda;
  star.track.cc[2] = carboncoremass;
  star.track.cf[2] = concentrationfactor;
  star.track.omega[2] = vrot/star.track.r[2];
  //set maximal values
  star.track.m[0] = fmax(star.track.m[1],star.track.m[2]);
  star.track.t[0] = fmax(star.track.t[1],star.track.t[2]);
  star.track.r[0] = fmax(star.track.r[1],star.track.r[2]);
  star.track.cm[0] = fmax(star.track.cm[1],star.track.cm[2]);
  star.track.llum[0] = fmax(star.track.llum[1],star.track.llum[2]);
  star.track.lteff[0] = fmax(star.track.lteff[1],star.track.lteff[2]);
  star.track.lambda[0] = fmax(star.track.lambda[1],star.track.lambda[2]);
  star.track.cc[0] = fmax(star.track.cc[1],star.track.cc[2]);
  star.track.cf[0] = fmax(star.track.cf[1],star.track.cf[2]);
  star.track.omega[0] = fmax(star.track.omega[1],star.track.omega[2]);

  //set stage
  star.stage = (int *)malloc(sizeof(int));
  if (star.stage==NULL) cerr << "#Error: memory allocation failed: star.stage" << endl;	//check if memory allocation fails
  star.stage[0] = kind;
  //set the mass(Msun)
  star.m = (double *)malloc(sizeof(double));
  if (star.m==NULL) cerr << "#Error: memory allocation failed: star.m" << endl;	//check if memory allocation fails
  star.m[0] = mass;
  //set the time(yr)
  star.t = (double *)malloc(sizeof(double));
  if (star.t==NULL) cerr << "#Error: memory allocation failed: star.t" << endl;	//check if memory allocation fails
  star.t[0] = time;
  //set the radius(Rsun)
  star.r = (double *)malloc(sizeof(double));
  if (star.r==NULL) cerr << "#Error: memory allocation failed: star.r" << endl;	//check if memory allocation fails
  star.r[0] = radius;
  //set the core mass(Msun)
  star.cm = (double *)malloc(sizeof(double));
  if (star.cm==NULL) cerr << "#Error: memory allocation failed: star.cm" << endl;	//check if memory allocation fails
  star.cm[0] = coremass;
  //set the luminosity(lg(lum/Lsun))
  star.llum = (double *)malloc(sizeof(double));
  if (star.llum==NULL) cerr << "#Error: memory allocation failed: star.llum" << endl;	//check if memory allocation fails
  star.llum[0] = llum;
  //set the effective temperature(lg(teff/K))
  star.lteff = (double *)malloc(sizeof(double));
  if (star.lteff==NULL) cerr << "#Error: memory allocation failed: star.lteff" << endl;	//check if memory allocation fails
  star.lteff[0] = lteff;
  //set the lambda
  star.lambda = (double *)malloc(sizeof(double));
  if (star.lambda==NULL) cerr << "#Error: memory allocation failed: star.lambda" << endl;	//check if memory allocation fails
  star.lambda[0] = lambda;
  //set the carbon core mass(Msun)
  star.cc = (double *)malloc(sizeof(double));
  if (star.cc==NULL) cerr << "#Error: memory allocation failed: star.cc" << endl;	//check if memory allocation fails
  star.cc[0] = carboncoremass;
  //set the concentration factor
  star.cf = (double *)malloc(sizeof(double));
  if (star.cf==NULL) cerr << "#Error: memory allocation failed: star.cf" << endl;	//check if memory allocation fails
  star.cf[0] = concentrationfactor;
  //set the angular velocity
  star.omega = (double *)malloc(sizeof(double));
  if (star.omega==NULL) cerr << "#Error: memory allocation failed: star.omega" << endl;	//check if memory allocation fails
  star.omega[0] = vrot/star.r[0];	//converting rotational velocity into angular velocity
/*  star.RLO = false;
  star.HeRLO = false;*/
  star.inimass = star.track.inimass;
  star.metal = metal;
  star.SN.inimass = 0.0;
  star.SN.Hecoremass = 0.0;
  star.SN.COcoremass = 0.0;
  star.SN.remnantmass = 0.0;
  star.SN.w = 0.0;
  star.SN.theta = 0.0;
  star.SN.phi = 0.0;
  star.SN.vsys = 0.0;
  star.SN.startype = 0;
  star.SN.type = 0;
}

void getbinaryparameters(t_system& system, int nextphase){
  int n=system.n-1;	//last entry index
  //determine values for the system
  if (nextphase>=0){	//determine current system/star values if system did not crash and is a binary
    system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
    system.qp[n] = system.prim.m[n]/system.sec.m[n];	//get mass ratio
    system.qs[n] = system.sec.m[n]/system.prim.m[n];	//get mass ratio
    if (system.a[n]>=0.0) system.P[n] = 2.0*M_PI*sqrt(cubic(system.a[n])/(G*system.M[n]));	//get period out of semi-major-axis
    else system.P[n] = ignorev;
    system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
    system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
    system.peri[n] = system.a[n]*(1.0-system.e[n]);	//get peri-astro
  }else if (nextphase>-100){
    system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
    system.qp[n] = ignorev;	//set mass ratio
    system.qs[n] = ignorev;	//set mass ratio
    system.a[n] = ignorev;	//set semi-major axis
    system.e[n] = ignorev;	//set eccentricity
    system.P[n] = ignorev;	//set period out of semi-major-axis
    system.rp[n] = ignorev;	//set Roche-Lobe-radius of primary
    system.rs[n] = ignorev;	//set Roche-Lobe-radius of secondary
    system.peri[n] = ignorev;	//set peri-astro

  }
}

void newphase(t_system& system){
  int n=system.n-1;	//last entry index
  //enlarge reserved memory
  system.prim.m = (double *)realloc(system.prim.m,(system.n+1)*sizeof(double));
  system.prim.t = (double *)realloc(system.prim.t,(system.n+1)*sizeof(double));
  system.prim.r = (double *)realloc(system.prim.r,(system.n+1)*sizeof(double));
  system.prim.cm = (double *)realloc(system.prim.cm,(system.n+1)*sizeof(double));
  system.prim.llum = (double *)realloc(system.prim.llum,(system.n+1)*sizeof(double));
  system.prim.lteff = (double *)realloc(system.prim.lteff,(system.n+1)*sizeof(double));
  system.prim.lambda = (double *)realloc(system.prim.lambda,(system.n+1)*sizeof(double));
  system.prim.cc = (double *)realloc(system.prim.cc,(system.n+1)*sizeof(double));
  system.prim.cf = (double *)realloc(system.prim.cf,(system.n+1)*sizeof(double));
  system.prim.omega = (double *)realloc(system.prim.omega,(system.n+1)*sizeof(double));
  system.prim.stage = (int *)realloc(system.prim.stage,(system.n+1)*sizeof(int));
  system.sec.m = (double *)realloc(system.sec.m,(system.n+1)*sizeof(double));
  system.sec.t = (double *)realloc(system.sec.t,(system.n+1)*sizeof(double));
  system.sec.r = (double *)realloc(system.sec.r,(system.n+1)*sizeof(double));
  system.sec.cm = (double *)realloc(system.sec.cm,(system.n+1)*sizeof(double));
  system.sec.llum = (double *)realloc(system.sec.llum,(system.n+1)*sizeof(double));
  system.sec.lteff = (double *)realloc(system.sec.lteff,(system.n+1)*sizeof(double));
  system.sec.lambda = (double *)realloc(system.sec.lambda,(system.n+1)*sizeof(double));
  system.sec.cc = (double *)realloc(system.sec.cc,(system.n+1)*sizeof(double));
  system.sec.cf = (double *)realloc(system.sec.cf,(system.n+1)*sizeof(double));
  system.sec.omega = (double *)realloc(system.sec.omega,(system.n+1)*sizeof(double));
  system.sec.stage = (int *)realloc(system.sec.stage,(system.n+1)*sizeof(int));
  system.M = (double *)realloc(system.M,(system.n+1)*sizeof(double));
  system.qp = (double *)realloc(system.qp,(system.n+1)*sizeof(double));
  system.qs = (double *)realloc(system.qs,(system.n+1)*sizeof(double));
  system.a = (double *)realloc(system.a,(system.n+1)*sizeof(double));
  system.P = (double *)realloc(system.P,(system.n+1)*sizeof(double));
  system.rp = (double *)realloc(system.rp,(system.n+1)*sizeof(double));
  system.rs = (double *)realloc(system.rs,(system.n+1)*sizeof(double));
  system.e = (double *)realloc(system.e,(system.n+1)*sizeof(double));
  system.peri = (double *)realloc(system.peri,(system.n+1)*sizeof(double));
  system.x = (double *)realloc(system.x,(system.n+1)*sizeof(double));
  system.y = (double *)realloc(system.y,(system.n+1)*sizeof(double));
  system.z = (double *)realloc(system.z,(system.n+1)*sizeof(double));
  system.vx = (double *)realloc(system.vx,(system.n+1)*sizeof(double));
  system.vy = (double *)realloc(system.vy,(system.n+1)*sizeof(double));
  system.vz = (double *)realloc(system.vz,(system.n+1)*sizeof(double));
  system.rhostar = (double *)realloc(system.rhostar,(system.n+1)*sizeof(double));
  system.phase = (int *)realloc(system.phase,(system.n+1)*sizeof(int));
  //check if memory allocation fails
  if (system.prim.m==NULL) cerr << "#Error: memory reallocation failed: system.prim.m" << endl;
  if (system.prim.t==NULL) cerr << "#Error: memory reallocation failed: system.prim.t" << endl;
  if (system.prim.r==NULL) cerr << "#Error: memory reallocation failed: system.prim.r" << endl;
  if (system.prim.cm==NULL) cerr << "#Error: memory reallocation failed: system.prim.cm" << endl;
  if (system.prim.llum==NULL) cerr << "#Error: memory reallocation failed: system.prim.llum" << endl;
  if (system.prim.lteff==NULL) cerr << "#Error: memory reallocation failed: system.prim.lteff" << endl;
  if (system.prim.lambda==NULL) cerr << "#Error: memory reallocation failed: system.prim.lambda" << endl;
  if (system.prim.cc==NULL) cerr << "#Error: memory reallocation failed: system.prim.cc" << endl;
  if (system.prim.cf==NULL) cerr << "#Error: memory reallocation failed: system.prim.cf" << endl;
  if (system.prim.omega==NULL) cerr << "#Error: memory reallocation failed: system.prim.omega" << endl;
  if (system.prim.stage==NULL) cerr << "#Error: memory reallocation failed: system.prim.stage" << endl;
  if (system.sec.m==NULL) cerr << "#Error: memory reallocation failed: system.sec.m" << endl;
  if (system.sec.t==NULL) cerr << "#Error: memory reallocation failed: system.sec.t" << endl;
  if (system.sec.r==NULL) cerr << "#Error: memory reallocation failed: system.sec.r" << endl;
  if (system.sec.cm==NULL) cerr << "#Error: memory reallocation failed: system.sec.cm" << endl;
  if (system.sec.llum==NULL) cerr << "#Error: memory reallocation failed: system.sec.llum" << endl;
  if (system.sec.lteff==NULL) cerr << "#Error: memory reallocation failed: system.sec.lteff" << endl;
  if (system.sec.lambda==NULL) cerr << "#Error: memory reallocation failed: system.sec.lambda" << endl;
  if (system.sec.cc==NULL) cerr << "#Error: memory reallocation failed: system.sec.cc" << endl;
  if (system.sec.cf==NULL) cerr << "#Error: memory reallocation failed: system.sec.cf" << endl;
  if (system.sec.omega==NULL) cerr << "#Error: memory reallocation failed: system.sec.omega" << endl;
  if (system.sec.stage==NULL) cerr << "#Error: memory reallocation failed: system.sec.stage" << endl;
  if (system.M==NULL) cerr << "#Error: memory reallocation failed: system.M" << endl;
  if (system.qp==NULL) cerr << "#Error: memory reallocation failed: system.qp" << endl;
  if (system.qs==NULL) cerr << "#Error: memory reallocation failed: system.qs" << endl;
  if (system.a==NULL) cerr << "#Error: memory reallocation failed: system.a" << endl;
  if (system.P==NULL) cerr << "#Error: memory reallocation failed: system.P" << endl;
  if (system.rp==NULL) cerr << "#Error: memory reallocation failed: system.rp" << endl;
  if (system.rs==NULL) cerr << "#Error: memory reallocation failed: system.rs" << endl;
  if (system.e==NULL) cerr << "#Error: memory reallocation failed: system.e" << endl;
  if (system.peri==NULL) cerr << "#Error: memory reallocation failed: system.peri" << endl;
  if (system.x==NULL) cerr << "#Error: memory reallocation failed: system.x" << endl;
  if (system.y==NULL) cerr << "#Error: memory reallocation failed: system.y" << endl;
  if (system.z==NULL) cerr << "#Error: memory reallocation failed: system.z" << endl;
  if (system.vx==NULL) cerr << "#Error: memory reallocation failed: system.vx" << endl;
  if (system.vy==NULL) cerr << "#Error: memory reallocation failed: system.vy" << endl;
  if (system.vz==NULL) cerr << "#Error: memory reallocation failed: system.vz" << endl;
  if (system.rhostar==NULL) cerr << "#Error: memory reallocation failed: system.rhostar" << endl;
  if (system.phase==NULL) cerr << "#Error: memory reallocation failed: system.phase" << endl;
  //copy last values to new values
  system.prim.m[system.n] = system.prim.m[n];
  system.prim.t[system.n] = system.prim.t[n];
  system.prim.r[system.n] = system.prim.r[n];
  system.prim.cm[system.n] = system.prim.cm[n];
  system.prim.llum[system.n] = system.prim.llum[n];
  system.prim.lteff[system.n] = system.prim.lteff[n];
  system.prim.lambda[system.n] = system.prim.lambda[n];
  system.prim.cc[system.n] = system.prim.cc[n];
  system.prim.cf[system.n] = system.prim.cf[n];
  system.prim.omega[system.n] = system.prim.omega[n];
  system.prim.stage[system.n] = system.prim.stage[n];
  system.sec.m[system.n] = system.sec.m[n];
  system.sec.t[system.n] = system.sec.t[n];
  system.sec.r[system.n] = system.sec.r[n];
  system.sec.cm[system.n] = system.sec.cm[n];
  system.sec.llum[system.n] = system.sec.llum[n];
  system.sec.lteff[system.n] = system.sec.lteff[n];
  system.sec.lambda[system.n] = system.sec.lambda[n];
  system.sec.cc[system.n] = system.sec.cc[n];
  system.sec.cf[system.n] = system.sec.cf[n];
  system.sec.omega[system.n] = system.sec.omega[n];
  system.sec.stage[system.n] = system.sec.stage[n];
  system.M[system.n] = system.M[n];
  system.qp[system.n] = system.qp[n];
  system.qs[system.n] = system.qs[n];
  system.a[system.n] = system.a[n];
  system.P[system.n] = system.P[n];
  system.rp[system.n] = system.rp[n];
  system.rs[system.n] = system.rs[n];
  system.e[system.n] = system.e[n];
  system.peri[system.n] = system.peri[n];
  system.x[system.n] = 0.0;	//no displacement
  system.y[system.n] = 0.0;	//no displacement
  system.z[system.n] = 0.0;	//no displacement
  system.vx[system.n] = 0.0;	//no kick
  system.vy[system.n] = 0.0;	//no kick
  system.vz[system.n] = 0.0;	//no kick
  system.rhostar[system.n] = system.rhostar[n];
  system.phase[system.n] = system.phase[n];
  system.stagechange = 0;
  //check for sensible values
  if(isnan(system.prim.m[system.n])||((system.prim.m[system.n]<=0)&&(system.prim.m[system.n]!=ignorev))){ screen = true; cerr << endl << "#prim.m=" << system.prim.m[system.n] << "Msun";}
  if(isnan(system.prim.t[system.n])||((system.prim.t[system.n]<0)&&(system.prim.t[system.n]!=ignorev))){ screen = true; cerr << endl << "#prim.t=" << system.prim.t[system.n] << "yr";}
  if(isnan(system.prim.r[system.n])||((system.prim.r[system.n]<=0)&&(system.prim.r[system.n]!=ignorev))){ screen = true; cerr << endl << "#prim.r=" << system.prim.r[system.n] << "Rsun";}
  if(isnan(system.prim.cm[system.n])||((system.prim.cm[system.n]<0)&&(system.prim.cm[system.n]!=ignorev))||(system.prim.cm[system.n]>system.prim.m[system.n]+accuracy)){ screen = true; cerr << endl << "#prim.cm=" << system.prim.cm[system.n] << "Msun prim.m=" << system.prim.m[system.n] << "Msun";}
  if(isnan(system.prim.llum[system.n])){ screen = true; cerr << endl << "#prim.llum=" << system.prim.llum[system.n];}
  if(isnan(system.prim.lteff[system.n])){ screen = true; cerr << endl << "#prim.lteff=" << system.prim.lteff[system.n];}
  if(isnan(system.prim.lambda[system.n])){ screen = true; cerr << endl << "#prim.lambda=" << system.prim.lambda[system.n];}
  if(isnan(system.prim.cc[system.n])||((system.prim.cc[system.n]<0)&&(system.prim.cc[system.n]!=ignorev))||(system.prim.cc[system.n]>system.prim.m[system.n]+accuracy)){ screen = true; cerr << endl << "#prim.cc=" << system.prim.cc[system.n] << "Msun prim.cm=" << system.prim.cm[system.n] << "Msun prim.m=" << system.prim.m[system.n] << "Msun";}
  if(isnan(system.prim.cf[system.n])||((system.prim.cf[system.n]<0)&&(system.prim.cf[system.n]!=ignorev))){ screen = true; cerr << endl << "#prim.cf=" << system.prim.cf[system.n];}
  if(isnan(system.prim.omega[system.n])){ screen = true; cerr << endl << "#prim.omega=" << system.prim.omega[system.n] << "/yr omega_crit=" << sqrt(G*system.prim.m[system.n]/cubic(system.prim.r[system.n])) << "/yr prim.vrot=" << system.prim.omega[system.n]*system.prim.r[system.n]*rroteq(system.prim.omega[system.n]/sqrt(2.0/3.0*G*system.prim.m[system.n]/cubic(system.prim.r[system.n])))*1.0e-3*Rsun/yr << "km/s" ;}	//||(system.prim.omega[system.n]>sqrt(G*system.prim.m[system.n]/cubic(system.prim.r[system.n])))
//  if(isnan(system.prim.stage[system.n])){ screen = true; cerr << endl << "#prim.stage=" << system.prim.stage[system.n];}
  if(isnan(system.sec.m[system.n])||((system.sec.m[system.n]<=0)&&(system.sec.m[system.n]!=ignorev))){ screen = true; cerr << endl << "#sec.m=" << system.sec.m[system.n] << "Msun";}
  if(isnan(system.sec.t[system.n])||((system.sec.t[system.n]<0)&&(system.sec.t[system.n]!=ignorev))){ screen = true; cerr << endl << "#sec.t=" << system.sec.t[system.n] << "yr";}
  if(isnan(system.sec.r[system.n])||((system.sec.r[system.n]<=0)&&(system.sec.r[system.n]!=ignorev))){ screen = true; cerr << endl << "#sec.r=" << system.sec.r[system.n] << "Rsun";}
  if(isnan(system.sec.cm[system.n])||((system.sec.cm[system.n]<0)&&(system.sec.cm[system.n]!=ignorev))||(system.sec.cm[system.n]>system.sec.m[system.n]+accuracy)){ screen = true; cerr << endl << "#sec.cm=" << system.sec.cm[system.n] << "Msun sec.m=" << system.sec.m[system.n] << "Msun";}
  if(isnan(system.sec.llum[system.n])){ screen = true; cerr << endl << "#sec.llum=" << system.sec.llum[system.n];}
  if(isnan(system.sec.lteff[system.n])){ screen = true; cerr << endl << "#sec.lteff=" << system.sec.lteff[system.n];}
  if(isnan(system.sec.lambda[system.n])){ screen = true; cerr << endl << "#sec.lambda=" << system.sec.lambda[system.n];}
  if(isnan(system.sec.cc[system.n])||((system.sec.cc[system.n]<0)&&(system.sec.cc[system.n]!=ignorev))||(system.sec.cc[system.n]>system.sec.m[system.n]+accuracy)){ screen = true; cerr << endl << "#sec.cc=" << system.sec.cc[system.n] << "Msun sec.cm=" << system.sec.cm[system.n] << "Msun sec.m=" << system.sec.m[system.n] << "Msun";}
  if(isnan(system.sec.cf[system.n])||((system.sec.cf[system.n]<0)&&(system.sec.cf[system.n]!=ignorev))){ screen = true; cerr << endl << "#sec.cf=" << system.sec.cf[system.n];}
  if(isnan(system.sec.omega[system.n])){ screen = true; cerr << endl << "#sec.omega=" << system.sec.omega[system.n] << "/yr omega_crit=" << sqrt(G*system.sec.m[system.n]/cubic(system.sec.r[system.n])) << "/yr prim.vrot=" << system.sec.omega[system.n]*system.sec.r[system.n]*rroteq(system.sec.omega[system.n]/sqrt(2.0/3.0*G*system.sec.m[system.n]/cubic(system.sec.r[system.n])))*1.0e-3*Rsun/yr << "km/s";}	//||(system.sec.omega[system.n]>sqrt(G*system.sec.m[system.n]/cubic(system.sec.r[system.n])))
//  if(isnan(system.sec.stage[system.n])){ screen = true; cerr << endl << "#sec.stage=" << system.sec.stage[system.n];}
  if(isnan(system.M[system.n])||((system.M[system.n]<=0)&&(system.M[system.n]!=ignorev))){ screen = true; cerr << endl << "#sys.M=" << system.M[system.n] << "Msun";}
  if(isnan(system.qp[system.n])||((system.qp[system.n]<=0)&&(system.qp[system.n]!=ignorev))){ screen = true; cerr << endl << "#sys.qp=" << system.qp[system.n];}
  if(isnan(system.qs[system.n])||((system.qs[system.n]<=0)&&(system.qs[system.n]!=ignorev))){ screen = true; cerr << endl << "#sys.qs=" << system.qs[system.n];}
  if(isnan(system.a[system.n])||((system.a[system.n]<=0)&&(system.a[system.n]!=ignorev)&&(system.phase[system.n]>=0)&&(system.phase[system.n]<90)&&(system.phase[system.n]!=4))){ screen = true; cerr << endl << "#sys.a=" << system.a[system.n] << "Rsun";}
  if(isnan(system.P[system.n])||((system.P[system.n]<=0)&&(system.P[system.n]!=ignorev)&&(system.phase[system.n]>=0)&&(system.phase[system.n]<90)&&(system.phase[system.n]!=4))){ screen = true; cerr << endl << "#sys.P=" << system.P[system.n] << "yr";}
  if(isnan(system.rp[system.n])){ screen = true; cerr << endl << "#sys.rp=" << system.rp[system.n] << "Rsun";}
  if(isnan(system.rs[system.n])){ screen = true; cerr << endl << "#sys.rs=" << system.rs[system.n] << "Rsun";}
  if(isnan(system.e[system.n])||((system.e[system.n]<0)&&(system.e[system.n]!=ignorev))){ screen = true; cerr << endl << "#sys.e=" << system.e[system.n];}
//  if(isnan(system.peri[system.n])||((system.peri[system.n]<=0)&&(system.peri[system.n]!=ignorev))){ screen = true; cerr << endl << "#sys.peri=" << system.peri[system.n] << "Rsun";}
  if(isnan(system.x[system.n])){ screen = true; cerr << endl << "#sys.x=" << system.x[system.n] << "pc";}
  if(isnan(system.y[system.n])){ screen = true; cerr << endl << "#sys.y=" << system.y[system.n] << "pc";}
  if(isnan(system.z[system.n])){ screen = true; cerr << endl << "#sys.z=" << system.z[system.n] << "pc";}
  if(isnan(system.vx[system.n])){ screen = true; cerr << endl << "#sys.vx=" << system.vx[system.n] << "km/s";}
  if(isnan(system.vy[system.n])){ screen = true; cerr << endl << "#sys.vy=" << system.vy[system.n] << "km/s";}
  if(isnan(system.vz[system.n])){ screen = true; cerr << endl << "#sys.vz=" << system.vz[system.n] << "km/s";}
  if(isnan(system.rhostar[system.n])||((system.rhostar[system.n]<0)&&(system.rhostar[system.n]!=ignorev))){ screen = true; cerr << endl << "#sys.rhostar=" << system.rhostar[system.n] << "/pc^3";}
//  if(isnan(system.phase[system.n])){ screen = true; cerr << endl << "#sys.phase=" << system.phase[system.n];}
  //increase number of array entries
  system.n++;
}

void destroysystem(t_system& system){
  //primary star
  free(system.prim.m);
  free(system.prim.t);
  free(system.prim.r);
  free(system.prim.cm);
  free(system.prim.llum);
  free(system.prim.lteff);
  free(system.prim.lambda);
  free(system.prim.cc);
  free(system.prim.cf);
  free(system.prim.omega);
  free(system.prim.stage);
  free(system.prim.track.m);
  free(system.prim.track.t);
  free(system.prim.track.r);
  free(system.prim.track.cm);
  free(system.prim.track.llum);
  free(system.prim.track.lteff);
  free(system.prim.track.lambda);
  free(system.prim.track.cc);
  free(system.prim.track.cf);
  free(system.prim.track.omega);
  //secondary star
  free(system.sec.m);
  free(system.sec.t);
  free(system.sec.r);
  free(system.sec.cm);
  free(system.sec.llum);
  free(system.sec.lteff);
  free(system.sec.lambda);
  free(system.sec.cc);
  free(system.sec.cf);
  free(system.sec.omega);
  free(system.sec.stage);
  free(system.sec.track.m);
  free(system.sec.track.t);
  free(system.sec.track.r);
  free(system.sec.track.cm);
  free(system.sec.track.llum);
  free(system.sec.track.lteff);
  free(system.sec.track.lambda);
  free(system.sec.track.cc);
  free(system.sec.track.cf);
  free(system.sec.track.omega);
  //system
  free(system.M);
  free(system.qp);
  free(system.qs);
  free(system.a);
  free(system.P);
  free(system.rp);
  free(system.rs);
  free(system.e);
  free(system.peri);
  free(system.x);
  free(system.y);
  free(system.z);
  free(system.vx);
  free(system.vy);
  free(system.vz);
  free(system.rhostar);
  free(system.phase);
  system.n = 0;
}

void outphasename(int phase){
  if (phase==0) cout << " phase=" << phase << " wind mass loss(following evolutionary tracks)" << endl;
  else if (phase==12) cout << " phase=" << phase << " Roche-lobe-overflow from primary star to secondary" << endl;
  else if (phase==21) cout << " phase=" << phase << " Roche-lobe-overflow from secondary star to primary" << endl;
  else if (phase==13) cout << " phase=" << phase << " common envelope from primary star" << endl;
  else if (phase==23) cout << " phase=" << phase << " common envelope from secondary star" << endl;
  else if (phase==4) cout << " phase=" << phase << " merger" << endl;
  else if (phase==15) cout << " phase=" << phase << " supernova explosion of primary";
  else if (phase==25) cout << " phase=" << phase << " supernova explosion of secondary";
  else if (phase==16) cout << " phase=" << phase << " planetary nebula formation arround primary" << endl;
  else if (phase==26) cout << " phase=" << phase << " planetary nebula formation arround secondary" << endl;
  else if (phase==-1) cout << " phase=" << phase << " single star system" << endl;
  else if (phase==94) cout << " phase=" << phase << " (gravitational) merger of two compact objects" << endl;
  else if (phase==99) cout << " phase=" << phase << " both stars are remnants" << endl;
  else if (phase==-99) cout << " phase=" << phase << " single star is a remnant" << endl;
  else if (phase==-100) cout << " phase=" << phase << " no stars left" << endl;
  else cout << " phase=" << phase << " undefined" << endl;
}

void outstagename(int stage, bool newline){
  if (stage==0) cout << " stage=" << stage << " H-burning star";
  else if (stage==1) cout << " stage=" << stage << " He-burning star";
  else if (stage==2) cout << " stage=" << stage << " naked He-star";
  else if (stage==3) cout << " stage=" << stage << " white dwarf";
  else if (stage==4) cout << " stage=" << stage << " neutron star";
  else if (stage==5) cout << " stage=" << stage << " black hole";
  else if (stage==6) cout << " stage=" << stage << " no star";
  else if (stage==-2) cout << " stage=" << stage << " stripped star";
  else if (stage==-3) cout << " stage=" << stage << " exploding star";
  else cout << " stage=" << stage << " undefined";
  if (newline) cout << endl;
}
