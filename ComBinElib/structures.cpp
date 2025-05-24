//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
#include "ComBinElib.h" // local header file

//using namespace std;

/*saves for separate evolution*/
t_HRD savetrack;	//HRD track of a star
int savermax;		//rmax of a star track

void initsystem(t_system& system, double mp, double q, double a, double e, double t, double metal, double omegap, double omegas, double x, double y, double z, double vx, double vy, double vz, double rhostar){
  int i;	//index variable

  system.n = 1;
  //set up the two components
  initstar(system.prim, mp, t, metal, omegap);
  initstar(system.sec, mp*q, t, metal, omegas);
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
    //check if memory allocation fails
    if (system.sec.track.m==NULL) cerr << "#Error: memory allocation failed: system.sec.track.m" << endl;
    if (system.sec.track.t==NULL) cerr << "#Error: memory allocation failed: system.sec.track.t" << endl;
    if (system.sec.track.r==NULL) cerr << "#Error: memory allocation failed: system.sec.track.r" << endl;
    if (system.sec.track.cm==NULL) cerr << "#Error: memory allocation failed: system.sec.track.cm" << endl;
    if (system.sec.track.llum==NULL) cerr << "#Error: memory allocation failed: system.sec.track.llum" << endl;
    if (system.sec.track.lteff==NULL) cerr << "#Error: memory allocation failed: system.sec.track.lteff" << endl;
    if (system.sec.track.lambda==NULL) cerr << "#Error: memory allocation failed: system.sec.track.lambda" << endl;
    //copy values
    system.sec.track.m[0] = savetrack.m[1];
    system.sec.track.t[0] = savetrack.t[1];
    system.sec.track.r[0] = savetrack.r[1];
    system.sec.track.cm[0] = savetrack.cm[1];
    system.sec.track.llum[0] = savetrack.llum[1];
    system.sec.track.lteff[0] = savetrack.lteff[1];
    system.sec.track.lambda[0] = savetrack.lambda[1];
    //value check
    cout << "memory check: &savetrack=" << &savetrack << " &system.sec.track=" << &system.sec.track << endl;
    cout << "&savetrack.(n,m,t,r,cm,llum,lteff)=\t\t" << &savetrack.n << ",\t" << savetrack.m << ",\t" << savetrack.t << ",\t" << savetrack.r << ",\t" << savetrack.cm << ",\t" << savetrack.llum << ",\t" << savetrack.lteff << "\n&system.sec.track.(n,m,t,r,cm,llum,lteff,lambda)=\t" << &system.sec.track.n << ",\t" << system.sec.track.m << ",\t" << system.sec.track.t << ",\t" << system.sec.track.r << ",\t" << system.sec.track.cm << ",\t" << system.sec.track.llum << ",\t" << system.sec.track.lteff << ",\t" << system.sec.track.lambda << endl;
    cout << "savetrack.n=" << savetrack.n << " savetrack.(m,t,r,cm,llum,lteff,lambda)=\n[(" << savetrack.m[0] << "," << savetrack.t[0] << "," << savetrack.r[0] << "," << savetrack.cm[0] << "," << savetrack.llum[0] << "," << savetrack.lteff[0] << "," << savetrack.lambda[0] << ")";
    for (i=1;i<savetrack.n;i++){
      cout << ", (" << savetrack.m[i] << "," << savetrack.t[i] << "," << savetrack.r[i] << "," << savetrack.cm[i] << "," << savetrack.llum[i] << "," << savetrack.lteff[i] << "," << savetrack.lambda[i] << ")";
      if (i<system.sec.track.n){	//copy values
        system.sec.track.m[i] = savetrack.m[i];
        system.sec.track.t[i] = savetrack.t[i];
        system.sec.track.r[i] = savetrack.r[i];
        system.sec.track.cm[i] = savetrack.cm[i];
        system.sec.track.llum[i] = savetrack.llum[i];
        system.sec.track.lteff[i] = savetrack.lteff[i];
        system.sec.track.lambda[i] = savetrack.lambda[i];
      }
    }
    system.sec.track.t[2] = 9.0e+99;	//new value
    system.sec.rmax = 1;	//new value
    cout << "]" << endl << "system.sec.track.n=" << system.sec.track.n << " system.sec.track.(m,t,r,cm,llum,lteff,lambda)=\n[(" << system.sec.track.m[0] << "," << system.sec.track.t[0] << "," << system.sec.track.r[0] << "," << system.sec.track.cm[0] << "," << system.sec.track.llum[0] << "," << system.sec.track.lteff[0] << "," << system.sec.track.lambda[0] << ")";
    for (i=1;i<system.sec.track.n;i++){
      cout << ", (" << system.sec.track.m[i] << "," << system.sec.track.t[i] << "," << system.sec.track.r[i] << "," << system.sec.track.cm[i] << "," << system.sec.track.llum[i] << "," << system.sec.track.lteff[i] << "," << system.sec.track.lambda[i] << ")";
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
  system.P[0] = 2.0*M_PI*sqrt(pow(system.a[0],3)/(G*system.M[0]));
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
  system.phase[0] = 0;
  system.stagechange = 0;
  system.formation = 0;
  system.tgw = 0.0;
}

void initstar(t_star& star, double mass, double time, double metal, double omega){
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
  if ((star.m[0]-star.cm[0]==0.0)&&(EbindperG==0.0)) star.lambda[0] = -1.0e+99;	//no interpolation
  else star.lambda[0] = star.m[0]*(star.m[0]-star.cm[0])/(EbindperG*star.r[0]);	//interpolate lambda	// E_bind=G*M*M_env / lambda*R
  //interpolate linearly the carbon core mass(Msun)
  star.cc = (double *)malloc(sizeof(double));
  if (star.cc==NULL) cerr << "#Error: memory allocation failed: star.cc" << endl;	//check if memory allocation fails
  star.cc[0] = star.track.cc[j]-tratio*(star.track.cc[j]-star.track.cc[j-1]);
  if (star.cc[0]<0) cerr << "#Error: cc=" << star.cc[0] << "Msun t_ratio=" << tratio << " j=" << j << " star.track.cc[j]=" << star.track.cc[j] << "Msun star.track.cc[j-1]=" << star.track.cc[j-1] << "Msun" << endl;	//check if core mass is negative (unphysical)
  //set initial value
  star.omega = (double *)malloc(sizeof(double));
  if (star.omega==NULL) cerr << "#Error: memory allocation failed: star.omega" << endl;	//check if memory allocation fails
  star.omega[0] = omega;
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
  if(isnan(system.prim.m[system.n])||((system.prim.m[system.n]<=0)&&(system.prim.m[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#prim.m=" << system.prim.m[system.n] << "Msun";}
  if(isnan(system.prim.t[system.n])||((system.prim.t[system.n]<0)&&(system.prim.t[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#prim.t=" << system.prim.t[system.n] << "yr";}
  if(isnan(system.prim.r[system.n])||((system.prim.r[system.n]<=0)&&(system.prim.r[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#prim.r=" << system.prim.r[system.n] << "Rsun";}
  if(isnan(system.prim.cm[system.n])||((system.prim.cm[system.n]<0)&&(system.prim.cm[system.n]!=-1.0e+99))||(system.prim.cm[system.n]>system.prim.m[system.n]+accuracy)){ screen = true; cerr << endl << "#prim.cm=" << system.prim.cm[system.n] << "Msun prim.m=" << system.prim.m[system.n] << "Msun";}
  if(isnan(system.prim.llum[system.n])){ screen = true; cerr << endl << "#prim.llum=" << system.prim.llum[system.n];}
  if(isnan(system.prim.lteff[system.n])){ screen = true; cerr << endl << "#prim.lteff=" << system.prim.lteff[system.n];}
  if(isnan(system.prim.lambda[system.n])){ screen = true; cerr << endl << "#prim.lambda=" << system.prim.lambda[system.n];}
  if(isnan(system.prim.cc[system.n])||((system.prim.cc[system.n]<0)&&(system.prim.cc[system.n]!=-1.0e+99))||(system.prim.cc[system.n]>system.prim.m[system.n]+accuracy)){ screen = true; cerr << endl << "#prim.cc=" << system.prim.cc[system.n] << "Msun prim.cm=" << system.prim.cm[system.n] << "Msun prim.m=" << system.prim.m[system.n] << "Msun";}
  if(isnan(system.prim.omega[system.n])){ screen = true; cerr << endl << "#prim.omega=" << system.prim.omega[system.n] << "/yr";}
//  if(isnan(system.prim.stage[system.n])){ screen = true; cerr << endl << "#prim.stage=" << system.prim.stage[system.n];}
  if(isnan(system.sec.m[system.n])||((system.sec.m[system.n]<=0)&&(system.sec.m[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sec.m=" << system.sec.m[system.n] << "Msun";}
  if(isnan(system.sec.t[system.n])||((system.sec.t[system.n]<0)&&(system.sec.t[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sec.t=" << system.sec.t[system.n] << "yr";}
  if(isnan(system.sec.r[system.n])||((system.sec.r[system.n]<=0)&&(system.sec.r[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sec.r=" << system.sec.r[system.n] << "Rsun";}
  if(isnan(system.sec.cm[system.n])||((system.sec.cm[system.n]<0)&&(system.sec.cm[system.n]!=-1.0e+99))||(system.sec.cm[system.n]>system.sec.m[system.n]+accuracy)){ screen = true; cerr << endl << "#sec.cm=" << system.sec.cm[system.n] << "Msun sec.m=" << system.sec.m[system.n] << "Msun";}
  if(isnan(system.sec.llum[system.n])){ screen = true; cerr << endl << "#sec.llum=" << system.sec.llum[system.n];}
  if(isnan(system.sec.lteff[system.n])){ screen = true; cerr << endl << "#sec.lteff=" << system.sec.lteff[system.n];}
  if(isnan(system.sec.lambda[system.n])){ screen = true; cerr << endl << "#sec.lambda=" << system.sec.lambda[system.n];}
  if(isnan(system.sec.cc[system.n])||((system.sec.cc[system.n]<0)&&(system.sec.cc[system.n]!=-1.0e+99))||(system.sec.cc[system.n]>system.sec.m[system.n]+accuracy)){ screen = true; cerr << endl << "#sec.cc=" << system.sec.cc[system.n] << "Msun sec.cm=" << system.sec.cm[system.n] << "Msun sec.m=" << system.sec.m[system.n] << "Msun";}
  if(isnan(system.sec.omega[system.n])){ screen = true; cerr << endl << "#sec.omega=" << system.sec.omega[system.n] << "/yr";}
//  if(isnan(system.sec.stage[system.n])){ screen = true; cerr << endl << "#sec.stage=" << system.sec.stage[system.n];}
  if(isnan(system.M[system.n])||((system.M[system.n]<=0)&&(system.M[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sys.M=" << system.M[system.n] << "Msun";}
  if(isnan(system.qp[system.n])||((system.qp[system.n]<=0)&&(system.qp[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sys.qp=" << system.qp[system.n];}
  if(isnan(system.qs[system.n])||((system.qs[system.n]<=0)&&(system.qs[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sys.qs=" << system.qs[system.n];}
  if(isnan(system.a[system.n])||((system.a[system.n]<=0)&&(system.a[system.n]!=-1.0e+99)&&(system.phase[system.n]>=0)&&(system.phase[system.n]<99)&&(system.phase[system.n]!=4))){ screen = true; cerr << endl << "#sys.a=" << system.a[system.n] << "Rsun";}
  if(isnan(system.P[system.n])||((system.P[system.n]<=0)&&(system.P[system.n]!=-1.0e+99)&&(system.phase[system.n]>=0)&&(system.phase[system.n]<99)&&(system.phase[system.n]!=4))){ screen = true; cerr << endl << "#sys.P=" << system.P[system.n] << "yr";}
  if(isnan(system.rp[system.n])){ screen = true; cerr << endl << "#sys.rp=" << system.rp[system.n] << "Rsun";}
  if(isnan(system.rs[system.n])){ screen = true; cerr << endl << "#sys.rs=" << system.rs[system.n] << "Rsun";}
  if(isnan(system.e[system.n])||((system.e[system.n]<0)&&(system.e[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sys.e=" << system.e[system.n];}
//  if(isnan(system.peri[system.n])||((system.peri[system.n]<=0)&&(system.peri[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sys.peri=" << system.peri[system.n] << "Rsun";}
  if(isnan(system.x[system.n])){ screen = true; cerr << endl << "#sys.x=" << system.x[system.n] << "pc";}
  if(isnan(system.y[system.n])){ screen = true; cerr << endl << "#sys.y=" << system.y[system.n] << "pc";}
  if(isnan(system.z[system.n])){ screen = true; cerr << endl << "#sys.z=" << system.z[system.n] << "pc";}
  if(isnan(system.vx[system.n])){ screen = true; cerr << endl << "#sys.vx=" << system.vx[system.n] << "km/s";}
  if(isnan(system.vy[system.n])){ screen = true; cerr << endl << "#sys.vy=" << system.vy[system.n] << "km/s";}
  if(isnan(system.vz[system.n])){ screen = true; cerr << endl << "#sys.vz=" << system.vz[system.n] << "km/s";}
  if(isnan(system.rhostar[system.n])||((system.rhostar[system.n]<0)&&(system.rhostar[system.n]!=-1.0e+99))){ screen = true; cerr << endl << "#sys.rhostar=" << system.rhostar[system.n] << "/pc^3";}
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
  //secondary star
  free(system.sec.m);
  free(system.sec.t);
  free(system.sec.r);
  free(system.sec.cm);
  free(system.sec.llum);
  free(system.sec.lteff);
  free(system.sec.lambda);
  free(system.sec.cc);
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
  else if (phase==-1) cout << " phase=" << phase << " system destroyed" << endl;
  else if (phase==94) cout << " phase=" << phase << " gravitational merger of two compact objects" << endl;
  else if (phase==99) cout << " phase=" << phase << " both stars are remnants" << endl;
  else cout << " phase=" << phase << " undefined" << endl;
}

void outstagename(int stage, bool newline){
  if (stage==0) cout << " stage=" << stage << " H-burning star";
  else if (stage==1) cout << " stage=" << stage << " He-burning star";
  else if (stage==2) cout << " stage=" << stage << " naked He-star";
  else if (stage==3) cout << " stage=" << stage << " white dwarf";
  else if (stage==4) cout << " stage=" << stage << " neutron star";
  else if (stage==5) cout << " stage=" << stage << " black hole";
  else if (stage==-3) cout << " stage=" << stage << " exploding star";
  else cout << " stage=" << stage << " undefined";
  if (newline) cout << endl;
}
