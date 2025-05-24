//This is the header file for the library ComBinElib.a

#include <iostream>	//in-/output to/from screen
#include <fstream>	//in-/output to/from file
#include <iomanip>	//manipulating the in-/output
#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
#include <math.h>	//provides the mathematical functions

using namespace std;	//allows to skip std::

/* Header for ComBinE.cpp*/
#ifndef ComBinE_h
#define ComBinE_h
//constants
extern const double Msun;	//solar Mass in kg
extern const double Rsun;	//solar Radius in m
extern const double yr;		//year in s
extern const double Lsun;	//solar Luminosity in Msun*Rsun^2/(yr^3)
extern const double G;		//graviational constant in Rsun^3/(Msun*yr^2)
extern const double c;		//seed of light in Rsun/yr
extern const double MH;		//ionized hydrogen mass(=proton mass) in Msun
//extern const double MHe;	//ionized helium mass in Msun
extern const double sigmaT;	//Thomson cross section of eletrons in Rsun^2
extern const double cgsEnergy;	//conversion: 1erg = cgsEnergy*Msun*Rsun^2/(yr^2)
extern const double day;	//(sidereal) day in s converted to yr
extern const double kms;	//conversion: 1km/s = kms*Rsun/yr
extern const double au;		//from PDG 2014: astronomical unit in m converted to Rsun
extern const double pc;		//from PDG 2014: parsec in m converted to Rsun
extern const double rsmax;	//Schwarzschild radius of a 1000Msun star
extern const double tmax;	// maximal time to emit graviational waves in yr
extern const double WDmax;	//maximal mass of a white dwarf (Chandrasekhar limit)
extern const double NSmax;	//maximal mass of a neutron star
//global variables
extern long number;		// number of calculated systems
extern bool screen;		// screen: additional output
extern bool debug;		// debug: additional output for debuging
extern bool separateevolution;	//evolve the stars one after the other while the companion is unevolved (for comparison with the old programm)
extern bool galacticintegration;	// calculate the movement of the system in the galaxy
extern double accuracy;		// minimal accuracy of calculation
extern double alphaRLO;		// direct wind mass loss from donor star during RLO 
extern double beta_const;	// all material transfered to secondary (non-deg.) star is not accreted
extern double Gamma;		// Gamma^2 = a_r/a where a_r is the radius of a circum stellar ring
extern double delta;		// no mass loss to circum stellar ring
extern double alphaCE;		// efficiency of converting E_orb to E_kin env. during CE
extern double lambda_const;	// lambda value for CE
extern double qlimit;		// somewhat arbitrary (both MS and comp.obj. accretors)
extern double Mhece;		// if RLO and M_he < Mhece => CE! - see Dewi & Pols (2002) 
extern double kickv;		// kick velocity in km/s; if negative take random value; if equal 0 no kick, else asymmetric SN (3 comp. Maxwellian)
extern double alphaTH;		// factor to interpolate between lambda with/without thermal and radiation energy (0.0-1.0: interpolate, -1: take old lambdatable, else assume viral equilibrium)
//extern int addran;		// number which is added to new random seed
extern char output;		// output format
extern bool single;		// single: run only one system with maximal values of Mp, Ms, a and e=0
#pragma omp threadprivate(screen, debug)

#endif

/* Header for tables.cpp*/
#ifndef tables_h
#define tables_h
typedef struct{
  int n;	//number of entries in the arrays
  int TAMS;	//array index of the TAMS
  double* m;	//mass in Msun
  double* t;	//time in yr
  double* r;	//radius in Rsun
  double* cm;	//core mass in Msun
  double* llum;	//luminosity in lg(llum/Lsun)
  double* lteff;	//effective temperature in lg(lteff/K)
  double* lambda;	//lambda
  double inimass;	//initial mass of the track
}t_HRD;
typedef struct{
  int n;	//number of entries in the arrays
  double* m;	//mass in Msun
  double* t;	//time in yr
  double* lteff;	//effective temperature in lg(lteff/K)
  double* llum;	//luminosity in lg(llum/Lsun)
  double* lr;	//radius in lg(lr/Rsun)
  double* mw;	//mass_w in Msun
  double* cm;	//core mass in Msun
  double* r;	//radius in Rsun
}t_He;
typedef struct{
  int n;	//number of entries in the arrays
  double* m;	//mass in Msun
  double* r;	//radius in Rsun
  double* lambda;	//lambda
}t_lambda;
void readstararray(int kind);	//read in star tables
void readhestararray();	//read in hestar tables
void readlambdaarray();	//read in lambda tables
void readtables(double metal, bool rotation);	//read in all tables
void destroytables();	//delete arrays and free their memory of all tables
void smoothtable(double* table, int n, int monotonic);	//smooth the table to avoid steps because of the precision in input data
//global variables containing the data from the tables
extern int nstar;		// nstar: number of stellar tracks in stararray
extern int nstar1;		// nstar1: number of stellar tracks in stararray1
extern int nstar2;		// nstar2: number of stellar tracks in stararray2
extern int nhe;			// nhe: number of naked helium tracks in hestararray
extern int nhe1;		// nhe1: number of naked helium tracks in hestararray1
extern t_HRD* stararray;	// stararray: tracks for stars of different masses
extern t_HRD* stararray1;	// stararray1: tracks for MS-stars of different masses
extern t_HRD* stararray2;	// stararray2: tracks for Core-Helium-burning-stars of different masses
extern t_He* hestararray;	// hestararray: tracks for He-stars of different masses
extern t_HRD* hestararray1;	// hestararray1: tracks for He-stars of different masses
extern t_lambda lambdaarray;	// lambdaarray: lambda values for different masses and radii
extern char tablelocation[1000];	//location of the tables in the file system
#endif

/* Header for stuctures.cpp*/
#ifndef stuctures_h
#define stuctures_h
typedef struct{
  double inimass;	//initial mass of the star/he-star
  double Hecoremass;	//pre-SN He-core mass
  double COcoremass;	//pre-SN CO-core mass
  double remnantmass;	//post-SN remnant mass
  double w;	//kick velocity in km/s
  double theta;	//kick angle theta in degree
  double phi;	//kick angle phi in degree
  double vsys;	//systemic velocity change by the kick in km/s
  int startype;	//stellar stage (see stage in t_star) of star for inimass
  int type;	//type of supernova: 0 no supernova; 1 planetary nebula instead of supernova; 2 electron capture supernova; 3 iron core collapse supernova; 4 collapse to blackhole
}t_SN;
typedef struct{
  double* m;	//mass in Msun
  double* t;	//time in yr
  double* r;	//radius in Rsun
  double* cm;	//core mass in Msun
  double* llum;	//luminosity in log(llum/Lsun)
  double* lteff;	//effective temperature in log(lteff/K)
  double* lambda;	//lambda
  double* omega;	//angular velocity
  int* stage;	//evolutionary stage: 0 hydrogen burning, 1 helium burning, 2 naked helium burning, 3 WD, 4 NS, 5 BH, -2 become 2, -3 become 3/4/5
  t_HRD track;	//HRD track of this star
  int rmax;	//track index where the radius is maximal
  int last;	//track index where the last time in track is
/*  bool RLO;	//RLO happend
  bool HeRLO;	//HeRLO happend*/
  double inimass;	//initial mass of the track
  double metal;	//metallicity
  t_SN SN;	//supernova values
}t_star;
typedef struct{
  int n;	//number of entries in the arrays
  t_star prim;	//primary star
  t_star sec;	//secondary star
  double* M;	//system mass in Msun
  double* qp;	//mass ratio = mass of primary/mass of secondary
  double* qs;	//mass ratio = mass of secondary/mass of primary
  double* a;	//semi major axis in Rsun
  double* P;	//period in yr
  double* rp;	//Roche-lobe radius of primary in Rsun
  double* rs;	//Roche-lobe radius of secondary in Rsun
  double* e;	//eccentricity
  double* peri;	//peri-astron in Rsun
  double* x;	//galactic x-position in pc
  double* y;	//galactic y-position in pc
  double* z;	//galactic z-position in pc
  double* vx;	//galactic x-velocity in km/s
  double* vy;	//galactic y-velocity in km/s
  double* vz;	//galactic z-velocity in km/s
  double* rhostar;	//density of surrounding stars in 1/pc^3
  int* phase;	//indicates the phase of the system
  int stagechange;	//indicate if the stage of one component has changed
  int formation;	//indicate which formation channel is used
  double tgw;	//time to merge by gravitational wave emission
}t_system;
/*phases:
 0 wind mass loss for both stars
 12 Roche-lobe overflow from primary to secondary
 21 Roche-lobe overflow from secondary to primary
 13 common envelope primary fills Roche lobe
 23 common envelope secondary fills Roche lobe
 4 merge
 15 supernova/planetary nebula of primary
 25 supernova/planetary nebula of secondary
 16 planetary nebula of primary
 26 planetary nebula of secondary
 -1 destroyed
 94 gravitational merger/both stars are remnants
 99 end stage/both stars are remnants
*/
extern t_HRD savetrack;	//HRD track of a star
extern int savermax;		//rmax of a star track
void initsystem(t_system& system, double mp, double q, double a, double e, double t, double metal, double omegap, double omegas, double x, double y, double z, double vx, double vy, double vz, double rhostar);	//create new system with specified parameters
void initstar(t_star& star, double mass, double time, double metal, double omega);	//create new star with specified parameters
void newphase(t_system& system);	//add a new phase to the system and copy the values from last phase
void destroysystem(t_system& system);	//delete arrays and free their memory
void outphasename(int phase);	//writes the name of the phase
void outstagename(int stage, bool newline);	//writes the name of the stage
#endif

/* Header for setup.cpp*/
#ifndef setup_h
#define setup_h
long newseed(int i, long s, long d);	//set a new seed for the random number generator
double ranv(int i);	//generate uniformly random numbers in [0,1[
void initial(t_system& system, double Mp_max, double Mp_min, double Ms_max, double Ms_min, double a_max, double a_min, double metal);	//generate new system
void getrandomdirection(double& theta, double& phi);	//get random direction on unit sphere
double getkickvelocity(double kickparameter);	//get kick velocity of SN (kind: 1=EC-SN; 2=Fe CCSN; 3=BH-formation)
double ranmaxwellian(double rms);	//returns a random value from a maxwellian distribution
double ranphase(double e);	//get a random phase angle for eccentric systems
void getMWposition(double& rho, double& phi, double& z);	//get position disk of MW_potential by Allen, C., Santillan, A., 1991, RevMexAA, 22, 255
//counts for usage of random number generators
extern long count_mp;		//counts primary mass generator
extern long count_q;		//counts mass ratio generator
extern long count_a;		//counts semi-major axis generator
extern long count_kick;		//counts kick parameter generator
extern long count_phase;	//counts phase parameter generator
extern long count_eccentricity;	//counts eccentricity generator
extern long count_galactic;	//counts galactic generator
extern int IMF;			// flag for initial mass function
extern int qdist;		// flag for initial mass ratio distribution
extern int adist;		// flag for initial semi-major axis distribution
extern int edist;		// flag for initial eccentricity distribution
extern int tdist;		// flag for initial age distribution
extern int Zdist;		// flag for initial metallicity distribution
extern int rotdist;		// flag for initial stellar spin distribution
extern int Rdist;		// flag for initial space distribution in a galaxy
extern int Vdist;		// flag for initial velocity distribution in a galaxy
extern int rhodist;		// flag for initial stellar density distribution in a galaxy
#endif

/* Header for evolution.cpp*/
#ifndef evolution_h
#define evolution_h
void evolve(t_system& system);	//evolve system
int wind(t_system& system, double& t, double& r);	//calculate time for normal evolution with wind mass loss
int overflow(t_system& system);	//apply a Roche-Lobe Overflow
int commonenvelope(t_system& system);	//apply a Common Envelope
int merge(t_system& system);	//apply a merger
int supernova(t_system& system);	//apply a supernova explosion
void circularize(t_system& system); //apply circularization
double convectiveenvelopefactor(double metal, double inimass);	//getting the factor in radius compared to maximum radius to have a convective envelope
#endif

/* Header for tracks.cpp*/
#ifndef tracks_h
#define tracks_h
void updatetrack(t_star& star, int n);	//updates the evolutionary track after mass transfer
void updatestarstrack(t_star& star, int n, t_HRD* stararray0, int nstar0);
//void updateHetrack(t_star& star, int n);	//updates the evolutionary track of naked helium star (determined by mass and core mass)
int nextstage(t_star& star, int n);	//set star to next evolutionary stage (helium-burning, naked helium phase, ...)
double ratio(double value, double high, double low);	//calculates the ratio and consider zero interval width
double radius_ini(double mass);	//calculates the radius of the star with given mass at ZAMS
void getnewtrack(t_HRD* stararray, int ilow, int iup, double mratio, int jpos, t_HRD& newtrack, int& jrmax);	//creates new evolutionary track
void changetrack(t_star& star, int n, t_HRD newtrack);	//replace track after current point by the new track
long double getA(long double j11,long double j12,long double j21,long double j22,long double m11,long double m12,long double m21,long double m22,long double cm11,long double cm12,long double cm21,long double cm22);	//calculates auxiliary variable A
long double getB(long double j21,long double j22,long double m,long double m11,long double m12,long double cm,long double cm11,long double cm12);	//calculates auxiliary variable B
long double getC(long double j11,long double j12,long double m,long double m21,long double m22,long double cm,long double cm21,long double cm22);	//calculates auxiliary variable C
long double getr1(long double m,long double m11,long double m12,long double m21,long double m22,long double j11,long double j12,long double j21,long double j22,long double r);	//calculates ratio in upper track
long double getr2(long double m,long double m11,long double m12,long double m21,long double m22,long double j11,long double j12,long double j21,long double j22,long double r);	//calculates ratio in lower track
double getCOcore(double& McoreHe);	//calculate the CO core mass from the He core mass by taking the CO core of the last track point in the naked Helium star grid
#endif

/* Header for radii.cpp*/
#ifndef radii_h
#define radii_h
double rocherad(double a, double q);	//calculate Roche-lobe radius
double RLFT(t_system system, int i, double& r, double rochemin);	//calculate the Roche-lobe-filling time(RLFT)
double RLFT2(t_star star, t_star companion, int n, int imin, bool same, double acircM, double& r);	//calculate the Roche-lobe-filling time(RLFT) with bisection method
double Schwarzschildradius(double M);	//calculate radius of a Black Hole
double NSradius(double M);	//calculate radius of a Neutron Star
double WDradius(double M);	//calculate radius of a White Dwarf
#endif

/* Header for gravitationalwave.cpp*/
#ifndef gravitationalwave_h
#define gravitationalwave_h
typedef struct{
  int n;	//number of entries in the arrays
  double* ecc;	//eccentricity
  double* integral;	//integral value
}t_ecc;
double mergetime(double M1, double M2, double a0, double e0);	//calculated the time for a system to merge by gravitational wave radiation
double intmergetime(double e0);	//interpolates the integral value in eq. 5.14 in Peters (1964)
void readmergetime();	//read in mergetime from file
void destroymergetime();	//delete arrays and free their memory
void newinttable();	//solves integral in eq. 5.14 in Peters (1964) and write it to a file
extern t_ecc eccarray;	// eccarray: values for the integral in eq. 5.14 in Peters (1964)
#endif

/* Header for statistics.cpp*/
#ifndef statistics_h
#define statistics_h
typedef struct{
  int n;	//number of bins in the arrays
  int subs;	//number of subdivisions in the arrays
  long ctot;	//total number of counts
  long* c;	//count
  long** subc;	//sub count
  double* x;	//x values of bin borders
  double min;	//minimal value
  double max;	//maximal value
  bool logscale;	//linear or logscale
}t_hist;
void initoutfiles();	//initialize output files 
void writedataout(int parallel, long seed, double Mp_max, double Mp_min, double Ms_max, double Ms_min, double a_max, double a_min, int galintegrator, int galpotential, long done);	//write the data.out file
void closeoutfiles();	//close output files 
void inithist(t_hist& hist, int nbin, double min, double max, bool logscale, int subs);	//initialize a histogram
void addtohist(t_hist& hist, double value, int sub);	//adds the value to the histogram
void freehist(t_hist& hist);	//free the histogram
void histout(t_hist& hist, bool perbinsize, char* name);	//writes the histogram data to file
void adddata(t_system& system, long i);	//add the system to the statistic, e.g. counters, histograms
/*counters*/
//warnings
extern long cnotrackupdate;	//count of skipped track updates
//events
extern long cRLO;	//count of Roche-Lobe overflows
extern long cRLOA;	//count of Roche-Lobe overflows of Case A
extern long cRLOB;	//count of Roche-Lobe overflows of Case B/C
extern long cRLOBB;	//count of Roche-Lobe overflows of Case BB
extern long cCE;	//count of common envelopes
extern long cCEA;	//count of common envelopes of Case A
extern long cCEB;	//count of common envelopes of Case B/C
extern long cCEBB;	//count of common envelopes of Case BB
extern long cMerger;	//count of merger events
extern long cMergerRLO;	//count of merger events in Roche-Lobe overflows
extern long cMergerCE;	//count of merger events in common envelopes
extern long cMergerSN;	//count of merger events in supernovae
extern long csupernova;	//count of supernova events
extern long cNS;	//count of fromed neutron stars
extern long cBH;	//count of fromed black holes
extern long cplanetarynebula;	//count of planetary nebula events
extern long cWD;	//count of fromed white dwarfs
//endstage
extern long cdestroyed;	//count of destroyed systems
extern long cmerged;	//count of merged systems
extern long cexploded;	//count of systems which are derupted by a supernova
extern long cunknown;	//count of systems which are derupted by an unknown process
extern long cfinal;	//count of systems which reach endstage
extern long cWDWD;	//count of white dwarf-white dwarf systems
extern long cWDNS;	//count of white dwarf-neutron star systems
extern long cWDBH;	//count of white dwarf-black hole systems
extern long cNSWD;	//count of neutron star-white dwarf systems
extern long cNSNS;	//count of neutron star-neutron star systems
extern long cNSBH;	//count of neutron star-black hole systems
extern long cBHWD;	//count of black hole-white dwarf systems
extern long cBHNS;	//count of black hole-neutron star systems
extern long cBHBH;	//count of black hole-black hole systems
extern long cgw;	//count of systems which merge by gravitational waves within tmax
extern long cgwWDWD;	//count of white dwarf-white dwarf systems which merge by gravitational waves within tmax
extern long cgwWDNS;	//count of white dwarf-neutron star systems which merge by gravitational waves within tmax
extern long cgwWDBH;	//count of white dwarf-black hole systems which merge by gravitational waves within tmax
extern long cgwNSWD;	//count of neutron star-white dwarf systems which merge by gravitational waves within tmax
extern long cgwNSNS;	//count of neutron star-neutron star systems which merge by gravitational waves within tmax
extern long cgwNSBH;	//count of neutron star-black hole systems which merge by gravitational waves within tmax
extern long cgwBHWD;	//count of black hole-white dwarf systems which merge by gravitational waves within tmax
extern long cgwBHNS;	//count of black hole-neutron star systems which merge by gravitational waves within tmax
extern long cgwBHBH;	//count of black hole-black hole systems which merge by gravitational waves within tmax
extern long cformation[10];
extern long cgwformation[10];
/*output data*/
extern bool data;	//write outputdata
extern bool dist;	//writeout the data of Mp, q, a, e, sigma, w, theta, phi
extern bool grid;	//make a grid of Mp, Ms and a
extern bool plot;	//plot histogramm(s)
extern std::ofstream dataout;	//enable writing to file the outputdata
extern std::ofstream inidist;	//enable writing to file the initial distributions
extern std::ofstream supernovadist;	//enable writing to file the supernova distributions
extern std::ofstream phasedist;	//enable writing to file the phase distributions
extern std::ofstream findist;	//enable writing to file the final distributions
extern std::ofstream gridout;	//enable writing to file the grid output data
extern FILE* distributionout;	//enable writing to a file
//extern std::ofstream distributionout;	//enable writing to a file
/*histogramms*/
extern t_hist tmerge;	//histogramm of merger times
#endif

/* Header for galactic.cpp*/
#ifndef galactic_h
#define galactic_h
typedef void (*integratorfunction)(double* y, double time);	//definition of a integrator funktion
typedef void (*potentialfunction)(double& ax, double& ay, double& az, double x, double y, double z);	//definition of a potential funktion
extern integratorfunction integrator;	//contains the used integrator
extern potentialfunction potential;	//contains the used potential
void movesystem(t_system& system, int m, int n);	//calculates movement of the system
void rk4(double* y, double time);	//Runge-Kutta 4 integrator: vector y={x,y,z,vx,vy,vz}
void MW_potential(double& ax, double& ay, double& az, double x, double y, double z);	//calculates the acceleration in the MW potential
#endif

/* Header for user.cpp*/
#ifndef user_h
#define user_h
void First();	//this function is called directly after the variable declaration in the main
void AfterReadParameters(double& Mp_max, double& Mp_min, double& Ms_max, double& Ms_min, double& a_max, double& a_min, double& Mp, double& Ms, double& a, double& e, double& metal, double& vp, double& vs, double& x, double& y, double& z, double& vx, double& vy, double& vz, double& rhostar, int& parallel, int& galintegrator, int& galpotential);	//this function is called after the parameters are read in from the command line or data.in
void PreReadTables();	//this function is called directly before the tables are read in
void AfterReadTables();	//this function is called after the tables are read in
void SetSeed(long& seed, long* discard, int n);	//this function is called directly before the random seeds are set; n is the number of entries in discard
void PreMainLoop();	//this function is called before the main loop of calculating the systems starts
double UserInitial_Mp(double Mp_max, double Mp_min);	//here you can define your own initial primary mass function
double UserInitial_q(double Ms_max, double Ms_min, double Mp);	//here you can define your own initial mass ratio function
double UserInitial_a(double a_max, double a_min, double Mp, double Ms);	//here you can define your own initial semi-major axis function
double UserInitial_e(double e_max, double Mp, double Ms, double a);	//here you can define your own initial eccentricity
double UserInitial_t(double Mp, double Ms, double a, double e);	//here you can define your own initial age
double UserInitial_Z(double Mp, double Ms, double a, double e, double t);	//here you can define your own initial metallicity
void UserInitial_omega(double& omegap, double& omegas, double Mp, double Ms, double a, double e, double t, double metal);	//here you can define your own initial angular spin
void UserInitial_R(double& x, double& y, double& z, double Mp, double Ms, double a, double e, double t, double metal, double omegap, double omegas);	//here you can define your own initial position in a galactic potential
void UserInitial_V(double& vx, double& vy, double& vz, double Mp, double Ms, double a, double e, double t, double metal, double omegap, double omegas, double x, double y, double z);	//here you can define your own initial velocity in a galactic potential
double UserInitial_rho(double Mp, double Ms, double a, double e, double t, double metal, double omegap, double omegas, double x, double y, double z, double vx, double vy, double vz);	//here you can define your own initial stellar density
void AfterEvolution(t_system& system);	//this function is called after the system is evolved, system contains all the evolution information of the system
long ToDebug(long i);	//this function returns which system is printed with debug-output
void AfterMainLoop();	//this function is called after the main loop of calculating the systems finishes
extern FILE* distfile1;	//enable writing to a file
extern FILE* distfile2;	//enable writing to a file
extern FILE* distfile3;	//enable writing to a file
extern bool snapshots;	//snapshots
extern FILE* snapshot0;	//enable writing to a file
extern FILE* snapshot1;	//enable writing to a file
extern FILE* snapshot2;	//enable writing to a file
extern FILE* snapshot3;	//enable writing to a file
extern FILE* snapshot4;	//enable writing to a file
extern FILE* snapshot5;	//enable writing to a file
extern FILE* snapshot6;	//enable writing to a file
extern FILE* snapshot7;	//enable writing to a file
extern FILE* snapshot8;	//enable writing to a file
extern FILE* snapshot9;	//enable writing to a file
extern FILE* snapshot10;	//enable writing to a file
extern FILE* snapshot11;	//enable writing to a file
extern FILE* snapshot12;	//enable writing to a file
extern FILE* snapshot13;	//enable writing to a file
#endif


