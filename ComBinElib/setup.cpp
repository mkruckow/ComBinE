//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
//#include <fstream>	//in-/output to/from file
#include <random>	//provides the random generator functions
#include "ComBinElib.h" // local header file

//using namespace std;

mt19937_64 generator_mp,generator_q,generator_a,generator_kick,generator_phase,generator_eccentricity,generator_galactic;	//random number generator
mt19937_64 save_mp,save_q,save_a,save_kick,save_phase,save_eccentricity,save_galactic;	//saves of random number generator
//counts for usage of random number generators
long count_mp=0;	//counts primary mass generator
long count_q=0;	//counts mass ratio generator
long count_a=0;	//counts semi-major axis generator
long count_kick=0;	//counts kick parameter generator
long count_phase=0;	//counts phase parameter generator
long count_eccentricity=0;	//counts eccentricity generator
long count_galactic=0;	//counts galactic generator
//uniform_real_distribution<double> distribution(0.0,1.0);
int IMF=1;			// flag for initial mass function
int qdist=1;			// flag for initial mass ratio distribution
int adist=1;			// flag for initial semi-major axis distribution
int edist=0;			// flag for initial eccentricity distribution
int tdist=0;			// flag for initial age distribution
int Zdist=0;			// flag for initial metallicity distribution
int rotdist=0;			// flag for initial stellar spin distribution
int Rdist=0;			// flag for initial space distribution in a galaxy
int Vdist=0;			// flag for initial velocity distribution in a galaxy
int rhodist=0;			// flag for initial stellar density distribution in a galaxy

long newseed(int i, long s, long d){	//set a new seed for the random number generator
  if (i==0){	//primary mass
    if (s==0) s = generator_mp();	//get random value
    generator_mp.seed(s);	//set seed
    count_mp = 0;
    //create image of the generatorstatus
/*    save_mp.seed(s);
    save_mp.discard(d);*/
    save_mp = generator_mp;
    //make one step
    generator_mp.discard(d+1);	//discard d+1 random numbers
  }else if (i==1){	//mass ratio
    if (s==0) s = generator_q();	//get random value
    generator_q.seed(s);	//set seed
    count_q = 0;
    //create image of the generatorstatus
/*    save_q.seed(s);
    save_q.discard(d);*/
    save_q = generator_q;
    //make one step
    generator_q.discard(d+1);	//discard d+1 random numbers
  }else if (i==2){	//semi-major-axis
    if (s==0) s = generator_a();	//get random value
    generator_a.seed(s);	//set seed
    count_a = 0;
    //create image of the generatorstatus
/*    save_a.seed(s);
    save_a.discard(d);*/
    save_a = generator_a;
    //make one step
    generator_a.discard(d+1);	//discard d+1 random numbers
  }else if (i==3){	//kick parameter
    if (s==0) s = generator_kick();	//get random value
    generator_kick.seed(s);	//set seed
    count_kick = 0;
    //create image of the generatorstatus
/*    save_kick.seed(s);
    save_kick.discard(d);*/
    save_kick = generator_kick;
    //make one step
    generator_kick.discard(d);	//discard d+1 random numbers
  }else if (i==4){	//phase parameter
    if (s==0) s = generator_phase();	//get random value
    generator_phase.seed(s);	//set seed
    count_phase = 0;
    //create image of the generatorstatus
/*    save_phase.seed(s);
    save_phase.discard(d);*/
    save_phase = generator_phase;
    //make one step
    generator_phase.discard(d);	//discard d+1 random numbers
  }else if (i==5){	//eccentricity
    if (s==0) s = generator_eccentricity();	//get random value
    generator_eccentricity.seed(s);	//set seed
    count_eccentricity = 0;
    //create image of the generatorstatus
/*    save_eccentricity.seed(s);
    save_eccentricity.discard(d);*/
    save_eccentricity = generator_eccentricity;
    //make one step
    generator_eccentricity.discard(d);	//discard d+1 random numbers
  }else if (i==6){	//galactic
    if (s==0) s = generator_galactic();	//get random value
    generator_galactic.seed(s);	//set seed
    count_galactic = 0;
    //create image of the generatorstatus
/*    save_galactic.seed(s);
    save_galactic.discard(d);*/
    save_galactic = generator_galactic;
    //make one step
    generator_galactic.discard(d);	//discard d+1 random numbers
  }else{
    cerr << "#Error: no random number generator selected(i=" << i << ") to set new seed" << endl;
    screen = true;
    return -1;
  }
  if (debug) cerr << "$newseed(i=" << i << ",s=" << s << ",d=" << d << ")" << endl;
  return s;	//return seed
}

double ranv(int i){
  unsigned long r,m;
  if (i==0){	//primary mass
#pragma omp critical
{
    r = generator_mp();
    if (save_mp==generator_mp){
//      if (newseed(i,r+1,0)!=(long)r) r=generator_mp();
//      cout << generator_mp() << "=" << save_mp();
      newseed(i,r+1,0);
      r = generator_mp();
    }
    m = generator_mp.max();
    count_mp++;
}	//end critical
  }else if (i==1){	//mass ratio
#pragma omp critical
{
    r = generator_q();
    if (save_q==generator_q){
//      if (newseed(i,r+1,0)!=(long)r) r=generator_q();
      newseed(i,r+1,0);
      r = generator_q();
    }
    m = generator_q.max();
    count_q++;
}	//end critical
  }else if (i==2){	//semi-major-axis
#pragma omp critical
{
    r = generator_a();
    if (save_a==generator_a){
//      if (newseed(i,r+1,0)!=(long)r) r=generator_a();
      newseed(i,r+1,0);
      r = generator_a();
    }
    m = generator_a.max();
    count_a++;
}	//end critical
  }else if (i==3){	//kick parameter
#pragma omp critical
{
    r = generator_kick();
    if (save_kick==generator_kick){
//      if (newseed(i,r+1,0)!=(long)r) r=generator_kick();
      newseed(i,r+1,0);
      r = generator_kick();
    }
    m = generator_kick.max();
    count_kick++;
}	//end critical
  }else if (i==4){	//phase parameter
#pragma omp critical
{
    r = generator_phase();
    if (save_phase==generator_phase){
//      if (newseed(i,r+1,0)!=(long)r) r=generator_phase();
      newseed(i,r+1,0);
      r = generator_phase();
    }
    m = generator_phase.max();
    count_phase++;
}	//end critical
  }else if (i==5){	//eccentricity
#pragma omp critical
{
    r = generator_eccentricity();
    if (save_eccentricity==generator_eccentricity){
//      if (newseed(i,r+1,0)!=(long)r) r=generator_eccentricity();
      newseed(i,r+1,0);
      r = generator_eccentricity();
    }
    m = generator_eccentricity.max();
    count_eccentricity++;
}	//end critical
  }else if (i==6){	//galactic
#pragma omp critical
{
    r = generator_galactic();
    if (save_galactic==generator_galactic){
//      if (newseed(i,r+1,0)!=(long)r) r=generator_galactic();
      newseed(i,r+1,0);
      r = generator_galactic();
    }
    m = generator_galactic.max();
    count_galactic++;
}	//end critical
  }else{
    cerr << "#Error: no random number generator selected to get random value" << endl;
    screen = true;
    return -1;
  }
//  srand(rand()+addran);	//set next random seed
//  return (double)rand()/(RAND_MAX+1.0);	//generate uniformly random numbers in [0,1[
//  generator.seed(r+addran);	//set seed and add value to seed from the random number generator
  return (double)r/(double)m;	//generate uniformly random numbers in [0,1[
//  return distribution(generator);	//generate uniformly random numbers in [0,1[
}

void initial(t_system& system, double Mp_max, double Mp_min, double Ms_max, double Ms_min, double a_max, double a_min, double metal){ 
// Using a Monte Carlo drawing to find initial Mp, q and a
  const double alpha = -2.7+1.0;	//using a Salpeter-like IMF: N(m) ~ m^-2.7 (Scalo 1986; Kroupa, Tout & Gilmore 1993) + 1.0 from integration
  static double kIMF = (pow(Mp_max,alpha)-pow(Mp_min,alpha));	//normalize IMF to 1
  //using the canonical IMF: N(m) ~ m^-1.3 for 0.08Msun<=m<=  0.5Msun (Kroupa 2008)	//lower boundary used from user input
  //using the canonical IMF: N(m) ~ m^-2.3 for 0.5 Msun<=m<=150.0Msun (Kroupa 2008)	//upper boundary used from user input
  static double m0 = Mp_min;	//lower mass boundary in Msun (Kroupa 2008)
  static double m1 = fmin(fmax(0.5,Mp_min),Mp_max);	//intermediate mass boundary in Msun (Kroupa 2008)
  static double m2 = Mp_max;	//upper mass boundary in Msun (Kroupa 2008)
  const double alpha1 = -1.3+1.0;	//power of lower part (Kroupa 2008) + 1.0 from integration
  const double alpha2 = -2.3+1.0;	//power of upper part (Kroupa 2008) + 1.0 from integration
  const double k2prok1IMF = pow(m1,alpha1-alpha2);	//ratio of constants in front of the two parts from continuity condition: k1*0.5^alpha1=k2*0.5^alpha2 => k2/k1=0.5^(alpha1-alpha2)	// + 1.0 from integration cancels out
  static double int1IMF = (pow(m1,alpha1)-pow(m0,alpha1))/alpha1;	//integral of first part without common prefactor
  static double int2IMF = k2prok1IMF*(pow(m2,alpha2)-pow(m1,alpha2))/alpha2;	//integral of second part without common prefactor
  static double normIMF = int1IMF+int2IMF;	//normalize IMF to 1 = k1*(mmid^alpha1-mlow^alpha1) + k2*(mup^alpha2-mmid^alpha2)
  double value;
  double q_max, q_min;
  const double kappa = -0.1+1.0;	//using a Sana mass ratio distribution: f(q) ~ q^-0.1 (Sana et al. 2012) + 1.0 from integration
  static double kq;
  double a_min_prim, a_min_sec, ka;
  double P_min, P_max, P;
  const double pi = -0.55+1.0;	//using a Sana period distribution: f(P) ~ P^-0.55 (Sana et al. 2012) + 1.0 from integration
  static double kP;
  double e_max;
  const double eta = -0.45+1.0;	//using a Sana eccentricity distribution: f(e) ~ e^-0.45 (Sana et al. 2012) + 1.0 from integration
//  static double ke;
  const double deltaP=45.0;
//  const double etaP=2.5;
  double Mp, q, a, e, t, omegap, omegas, x, y, z, rho, phi, vx, vy, vz, ax, ay, az, rhostar;
  double Ms;

//Picking random masses for the primary star between Mp_min and Mp_max in Msun
  //NOTE: The IMF of primaries may differ from that of single stars, see
  //Mermillod or EDVDH p.480 in D. Vanbeveren (Ed., conf.proc.) 2001.
  if (((IMF<1)||(IMF>2))&&(IMF>-10)){
    cerr << "#Warning: No initial mass function spezified: use IMF=1" << endl;
    IMF = 1;
  }
  if (IMF==1){
    //use salpeter IMF
    Mp = pow(pow(Mp_min,alpha)+kIMF*ranv(0),1.0/alpha);	//generate random mass from IMF
  }else if (IMF==2){
    //use canonical IMF (Kroupa 2008)
    value = normIMF*ranv(0);
    //generate random mass from IMF
    if (value<int1IMF) Mp = pow(pow(m0,alpha1)+alpha1*value,1.0/alpha1);	//within lower part
    else Mp = pow(pow(m1,alpha2)+alpha2*(value-int1IMF)/k2prok1IMF,1.0/alpha2);	//within upper part
    if (((Mp>m1)&&(value<int1IMF))||((Mp<m1)&&(value>int1IMF))){
      cerr << "#Error: Mp=" << Mp << "Msun m1=" << m1 << "Msun value=" << value << " int1IMF=" << int1IMF << endl;
      screen = true;
    }
  }else if (IMF<-9){
    Mp = UserInitial_Mp(Mp_max, Mp_min);
  }else{
    cerr << "#Error: No initial mass function spezified! Mp set to Mp_max=" << Mp_max << endl;
    Mp = Mp_max;
  }
  if ((Mp<Mp_min)||(Mp>Mp_max)){
    cerr << "#Error: Mp/Msun=" << Mp << " not in [" << Mp_min << "," << Mp_max << "]" << endl;
    screen = true;
  }

//Picking the mass ratio between q_min(Ms_min,Mp) and q_max(Ms_max,Mp) where q=Ms/Mp
  q_max = fmin(Ms_max/Mp,1.0);
  q_min = fmax(Ms_min/Mp,0.0);
  if (((qdist<1)||(qdist>3))&&(qdist>-10)){
    cerr << "#Warning: No initial mass ratio distribution spezified: use qdist=1" << endl;
    qdist = 1;
  }
  if (qdist==1){
    //from Kuiper (1935), f(q)=1/(1+q)^2 
    //NOTE: This is the distribution you get if you randomly chop
    //a piece of something (anything) into two parts
    //- see EDVDH in D. Vanbeveren (2001).
    //normalizing: f(q)=1/(1+q)^2  =>  f(q)=2/(1+q)^2 for q between 0 and 1, see eq.(4) Kuiper (1935) in PASP 47,15
    kq = (1.0+q_min)*(1.0+q_max)/(q_max-q_min);	//normalize q-distribution to unity for q between q_min and q_max
    q = ((1.0+q_min)*kq/(kq-(1.0+q_min)*(1.0-ranv(1)))-1.0);	//generate random q from q-distribution
  }else if (qdist==2){
    //flat distribution f(q) = 1
    q = q_min+(q_max-q_min)*ranv(1);	//generate random q between q_min and q_max
  }else if (qdist==3){
    //mass ratio distribution from Sana et al.(2012) f(q) ~ q^-0.1
    kq = (pow(q_max,kappa)-pow(q_min,kappa));	//normalize q-distribution to unity for q between q_min and q_max
    q = pow(pow(q_min,kappa)+kq*ranv(1),1.0/kappa);	//generate random q between q_min and q_max
  }else if (qdist<-9){
    q = UserInitial_q(Ms_max, Ms_min, Mp);
  }else{
    cerr << "#Error: No initial mass ratio distribution spezified! q set to q_max=" << q_max << endl;
    q = q_max;
  }
  Ms = Mp*q;

//Picking random separations between a_min and a_max in Rsun
  if (((adist<1)||(adist>3))&&(adist>-10)){
    cerr << "#Warning: No initial semi-major axis distribution spezified: use adist=1" << endl;
    adist = 1;
  }
  a_min_prim = radius_ini(Mp)/rocherad(1.0,1.0/q);	//minimal separation that primary fills roche-lobe
  a_min_sec = radius_ini(q*Mp)/rocherad(1.0,q);	//minimal separation that secondary fills roche-lobe
  a_min = fmax(a_min,fmax(a_min_prim,a_min_sec));	//reset minimal seperation in Rsun
  if (adist==1){
    //Using flat logP-dist: N(a)=1/a (e.g. Abt)
    ka = log(a_max)-log(a_min);	//normalize a-distribution to unity for a between a_min and a_max
    a = exp(log(a_min)+ka*ranv(2));	//generate random a from a-distribution
  }else if (adist==2){
    //Using period distribution from Kroupa (2008)
    P_min = 2.0*M_PI*sqrt(pow(a_min,3)/(G*Mp*(1.0+q)))/day;	//get minimum Period in days
    P_max = 2.0*M_PI*sqrt(pow(a_max,3)/(G*Mp*(1.0+q)))/day;	//get maximum Period in days
    kP = log(pow(log10(P_max/P_min),2)/deltaP+1.0);	//get normalisation in original formular kp=2.0/eta
    P = pow(10,sqrt(deltaP*(exp(kP*ranv(2))-1.0)))*P_min;	//generate random a from P-distribution in days
    a = pow(pow(P*day/(2.0*M_PI),2)*(G*Mp*(1.0+q)),1.0/3.0);	//convert back to semi-major axis in Rsun
  }else if (adist==3){
    //period distribution from Sana et al.(2012) f(P) ~ P^-0.55
    P_min = 2.0*M_PI*sqrt(pow(a_min,3)/(G*Mp*(1.0+q)))/day;	//get minimum Period in days
    P_max = 2.0*M_PI*sqrt(pow(a_max,3)/(G*Mp*(1.0+q)))/day;	//get maximum Period in days
    kP = (pow(P_max,pi)-pow(P_min,pi));	//normalize P-distribution to unity for P between P_min and P_max
    P = pow(pow(P_min,pi)+kP*ranv(2),1.0/pi);	//generate random P between P_min and P_max
    a = pow(pow(P*day/(2.0*M_PI),2)*(G*Mp*(1.0+q)),1.0/3.0);	//convert back to semi-major axis in Rsun
  }else if (adist<-9){
    a = a_min-1.0;
    while ((a<a_min)||(a>a_max)){
      a = UserInitial_a(a_max, a_min, Mp, Ms);
      if ((a<a_min)||(a>a_max)) cerr << "#Warning: The user defined semi-major axis is out of range: a=" << a << "Rsun a_min=" << a_min << "Rsun a_max=" << a_max << "Rsun" << endl;
    }
  }else{
    cerr << "#Error: No initial semi-major axis distribution spezified! a set to a_max=" << a_max << endl;
    a = a_max;
  }
  //##Correlation between system mass and semi-major axis?

//Picking random eccentricities between 0 and 1
  e_max = 1.0-a_min/a;	//a_min is the minimal peri-astron separation -> maximum eccentricity
  if (((edist<0)||(edist>3))&&(edist>-10)){
    cerr << "#Warning: No initial eccentricity distribution spezified: use edist=0" << endl;
    edist = 0;
  }
  if (edist==0){
    //circular
    e = 0.0;
  }else if (edist==1){
    //thermal-equilibrium eccentricity distribution f(e)=2*e
    e = sqrt(ranv(5))*e_max;
  }else if (edist==2){
    //flat distribution f(e)=1
    e = ranv(5)*e_max;
  }else if (edist==3){
    //flat angular momentum-distribution g(L)=1 and e~sqrt(1-L*L)
    e = sqrt(1.0-pow(sqrt(1.0-e_max*e_max)+(1.0-sqrt(1.0-e_max*e_max))*ranv(5),2));
  }else if (edist==4){
    //eccentricity distribution from Sana et al.(2012) f(e) ~ e^-0.45
/*    ke = (pow(1.0,eta)-pow(0.0,eta));	//normalize e-distribution to unity for e between 0 and 1
    e = pow(pow(0.0,eta)+ke*ranv(1),1.0/eta);	//generate random e between 0 and 1*/
    e = pow(ranv(5),1.0/eta);	//generate random e between 0 and 1
  }else if (edist<-9){
    e = -1.0;
    while ((e<0.0)||(e>e_max)){
      e = UserInitial_e(e_max, Mp, Ms, a);
      if ((e<0.0)||(e>e_max)) cerr << "#Warning: The user defined eccentricity is out of range: e=" << e << " e_min=0.0 e_max=" << e_max << endl;
    }
  }else{
    cerr << "#Error: No initial eccentricity distribution spezified! e set to 0.0" << endl;
    e = 0.0;
  }
  //##Eccentricity distribution should be used, maybe with dependense on the semi-major-axis

  if (((tdist<0)||(tdist>0))&&(tdist>-10)){
    cerr << "#Warning: No initial age distribution spezified: use tdist=0" << endl;
    tdist = 0;
  }
  if (tdist==0){
    t = 0.0;	//at ZAMS
  }else if (tdist<-9){
    t = UserInitial_t(Mp, Ms, a, e);
  }else{
    cerr << "#Error: No initial age distribution spezified! t set to 0.0" << endl;
    t = 0.0;
  }

  if (((Zdist<-3)||(Zdist>0))&&(Zdist>-10)){
    cerr << "#Warning: No initial metallicity distribution spezified: use Zdist=0" << endl;
    Zdist = 0;
  }
  if (Zdist==0){
    metal = 0.009;	//set metallicity: Milky Way metallicity(~0.009)
  }else if (Zdist==-1){
    metal = 0.005;		//set metallicity: LMC(~0.005)
  }else if (Zdist==-2){
    metal = 0.002;		//set metallicity: SMC(~0.002)
  }else if (Zdist==-3){
    metal = 0.0002;		//set metallicity: IZw18(~0.0002)
  }else if (Zdist<-9){
    metal = UserInitial_Z(Mp, Ms, a, e, t);
  }else{
    cerr << "#Error: No initial metallicity distribution spezified! metallicity set to 0.02" << endl;
    metal = 0.02;
  }

  if (((rotdist<-1)||(rotdist>0))&&(rotdist>-10)){
    cerr << "#Warning: No initial stellar spin distribution spezified: use rotdist=0" << endl;
    rotdist = 0;
  }
  if (rotdist==-1){
    //synchronised
    omegap = 1.0/(2.0*M_PI*sqrt(pow(a,3)/(G*Mp*(1.0+q))));
    omegas = omegap;
  }else if (rotdist==0){
    omegap = 0.0;	//non rotating primary
    omegas = 0.0;	//non rotating secondary
  }else if (rotdist<-9){
    UserInitial_omega(omegap, omegas, Mp, Ms, a, e, t, metal);
  }else{
    cerr << "#Error: No initial stellar spin distribution spezified! omegap and omegas set to 0.0" << endl;
    omegap = 0.0;	//non rotating primary
    omegas = 0.0;	//non rotating secondary
  }

  if (((Rdist<-1)||(Rdist>1))&&(Rdist>-10)){
    cerr << "#Warning: No initial space distribution in a galaxy spezified: use Rdist=0" << endl;
    Rdist = 0;
  }
  if (Rdist==-1){
    //use solar values
    x = -8500.0;
    y =     0.0;
    z =     0.0;
  }else if (Rdist==0){
    x = 0.0;	//at galactic center
    y = 0.0;	//at galactic center
    z = 0.0;	//at galactic center
  }else if (Rdist==1){	//disk of MW_potential by Allen, C., Santillan, A., 1991, RevMexAA, 22, 255
    getMWposition(rho, phi, z);
    x = rho*cos(phi);	//convert
    y = rho*sin(phi);	//convert
    x *= 1000.0;	//convert from kpc to pc
    y *= 1000.0;	//convert from kpc to pc
    z *= 1000.0;	//convert from kpc to pc
  }else if (Rdist<-9){
    UserInitial_R(x, y, z, Mp, Ms, a, e, t, metal, omegap, omegas);
  }else{
    cerr << "#Error: No initial space distribution in a galaxy spezified! position set to center" << endl;
    x = 0.0;	//at galactic center
    y = 0.0;	//at galactic center
    z = 0.0;	//at galactic center
  }

  if (((Vdist<-1)||(Vdist>1))&&(Vdist>-10)){
    cerr << "#Warning: No initial velocity distribution in a galaxy spezified: use Vdist=0" << endl;
    Vdist = 0;
  }
  if (Vdist==-1){
    //use solar values
    vx =   0.0+10.0;
    vy = 220.0+15.0;
    vz =   0.0+ 7.0;
  }else if (Vdist==0){
    vx = 0.0;	//at rest
    vy = 0.0;	//at rest
    vz = 0.0;	//at rest
  }else if (Vdist==1){
    potential(ax, ay, az, x*0.001, y*0.001, z*0.001);	//get acceleration(pc/yr^2) at galactic position(kpc)
    vx = sqrt(fabs(ax*x+ay*y+az*z)/(1.0+ax*ax/(ay*ay)));	//velocity perpendicular to acceleration
    vy = -ax/ay*vx;	//velocity perpendicular to acceleration
    vz = 0.0;	//no vertial velocity
    vx *= pc/kms;	//convert velocity from pc/yr to km/s
    vy *= pc/kms;	//convert velocity from pc/yr to km/s
    vz *= pc/kms;	//convert velocity from pc/yr to km/s
//    cout << endl << "a={" << ax << ", " << ay << "," << az << "}\tr={" << x << ", " << y << "," << z << "}\tv={" << vx << ", " << vy << "," << vz << "}" << endl;
  }else if (Vdist<-9){
    UserInitial_V(vx, vy, vz, Mp, Ms, a, e, t, metal, omegap, omegas, x, y, z);
  }else{
    cerr << "#Error: No initial velocity distribution in a galaxy spezified! v set to 0.0" << endl;
    vx = 0.0;	//at rest
    vy = 0.0;	//at rest
    vz = 0.0;	//at rest
  }

  if (((rhodist<0)||(rhodist>0))&&(rhodist>-10)){
    cerr << "#Warning: No initial stellar density distribution in a galaxy spezified: use rhodist=0" << endl;
    rhodist = 0;
  }
  if (rhodist==0){
    rhostar = 0.0;	//field star
  }else if (rhodist<-9){
    rhostar = UserInitial_rho(Mp, Ms, a, e, t, metal, omegap, omegas, x, y, z, vx, vy, vz);
  }else{
    cerr << "#Error: No initial stellar density distribution in a galaxy spezified! star is set in the field" << endl;
    rhostar = 0.0;	//field star
  }

  if (screen) cout << "System parameters: M_prim(Mp)=" << Mp << "Msun mass ratio(q)=" << q << " semi-major-axis(a)=" << a << "Rsun eccentricity(e)=" << e << " age(t)=" << t << "yr metallicity=" << metal << " omega_prim(omegap)=" << omegap << "/yr omega_sec(omegas)=" << omegas << "/yr galactic position(x,y,z)=(" << x << ", " << y << ", " << z << ")pc galactic velocity(vx,vy,vz)=(" << vx << ", " << vy << ", " << vz << ")km/s stellar density(rhostar)=" << rhostar << "/pc^3" << endl;
  initsystem(system, Mp, q, a, e, t, metal, omegap, omegas, x, y, z, vx, vy, vz, rhostar);	//create new system with specified parameters
}

void getrandomdirection(double& theta, double& phi){	//get random direction on unit sphere
  theta = -1.0e9;
  phi = -1.0e9;
  if ((kickv<0)&& single){
    cout << "write kickangle(theta) in degree(<0:random value):" << endl;
    cin >> theta;
    theta *= M_PI/180.0;
    cout << "write kickangle(phi) in degree(<0:random value):" << endl;
    cin >> phi;
    phi *= M_PI/180.0;
  }
  if (theta<0){
    theta = acos(1.0-2.0*ranv(3));	//get theta in [0,pi[ uniform in cos(theta)
  }
  if (phi<0){
    phi = 2.0*M_PI*ranv(3);	//get phi uniform in [0,2pi[
  }
/*  theta = 162.251002/180.0*M_PI;
  phi = 22.748833/180.0*M_PI;*/
}

double getkickvelocity(double kickparameter){	//get kick velocity of SN
  double kick=kickv;
  if (kickv<0){
//    kickparameter *= 0.5;
    if (single){
      cout << "write kickvelocity(kickparameter=" << kickparameter << ") in km/s(<0:random value):" << endl;
      cin >> kick;
    }
    if (kick>=0){
      return kick;
    }else if (kickparameter<0){	//iron core collapse SN	//Hobbs et al. 2005 for kickparameter=265km/s
      return ranmaxwellian(-kickparameter);	//get value from maxwell distribution with rms of "kickparameter" km/s
    }else if (kickparameter>0){	//black hole formation(kickparameter=200km/s) / electron capture SN(kickparameter=50km/s)
      return kickparameter*ranv(3);	//random value of a flat distribution in km/s
    }else{
      return 0.0;
    }
  }else{	//fixed kickvelocity
    return kickv;
  }
}

double ranmaxwellian(double rms){	//returns a random value from a maxwellian distribution
  double a,max,factor,d,w;	//rms of the maxwell distribution, scale parameter of the maxwell distribution, maximum distribution value, prefactor, distribution value at w
  a = rms/sqrt(3.0);	//get scale parameter
  max = 4.0/(sqrt(2.0*M_PI)*a*exp(1.0));	//maximum value of the maxwell distribution (reached at w=sqrt(2)*a)
  factor = sqrt(2.0/M_PI)/(a*a*a);	//prefactor
  if (debug) cerr << " ==> a=" << a << " max=" << max << " factor=" << factor;
  do{	//generate velocity with maxwell distribution in [0,9*a[
    w = 9*a*ranv(3);	//guess velocity
    d = factor*w*w*exp(-w*w/(2.0*a*a));	//get value from maxwell distribution
    if (debug) cerr << " w=" << w << " d=" << d;
  }while (d<max*ranv(3));	//check probability of the guessed velocity
  if (debug) cerr << endl;
  return w;
}

double ranphase(double e){	//get a random phase angle for eccentric systems for a flat eccentric anomaly
  double phase = 2.0*M_PI*ranv(4);	//get phase(eccentric anomaly) uniform in [0,2pi[
  double y = (1.0+e)*ranv(4);
  while (y>1.0-e*cos(phase)){
    phase = 2.0*M_PI*ranv(4);	//get phase(eccentric anomaly) uniform in [0,2pi[
    y = (1.0+e)*ranv(4);
  }
/*  if (phase<M_PI){	//convert eccentric anomaly to true anomaly
    phase = acos((cos(phase)-e)/(1.0-e*cos(phase)));
  }else{
    phase = 2.0*M_PI-acos((cos(phase)-e)/(1.0-e*cos(phase)));
  }*/
  //writes distribution of the phase
  if (dist) phasedist << phase*180.0/M_PI << endl;
  return phase;
}

void getMWposition(double& rho, double& phi, double& z){	//get position disk of MW_potential by Allen, C., Santillan, A., 1991, RevMexAA, 22, 255
  const double a2 = 5.3178;	//in kpc
  const double b2 = 0.2500;	//in kpc
  const double factor = pow(a2+b2,3.0)*pow(b2,3.0)/(a2+3.0*b2);
  const double g = 200.0;
  const double g3 = pow(g,1.0/3.0);
  double d,zb;
  double ran, frho;
  do{
    ran = ranv(6);
    if (ran<2.0/3.0){
      rho = 1.5*g3*ran;
      frho = 1.0;
    }else{
      rho = g3*pow(3.0-3.0*ran,-0.5);
      frho = g*pow(rho,-3.0);
    }
//    rho = 1000.0*ranv(6);	//get rho uniform in [0,1000[
    z = 44.0*ranv(6)-22.0;		//get z uniform in [-22,22[
    zb = sqrt(z*z+b2*b2);
    d = (a2*rho*rho+(a2+3.0*zb)*(a2+zb)*(a2+zb))*pow(rho*rho+(a2+zb)*(a2+zb),-2.5)*pow(zb,-3.0)*factor;
  }while(d<frho*ranv(6));
  phi = 2.0*M_PI*ranv(6);	//get phi uniform in [0,2pi[
}
