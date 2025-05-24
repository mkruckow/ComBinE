// Kruckow & Tauris 2013-2018. version 1.1.0
// Based on Starburst by Tauris & Voss August 2001.
// Simulation of NS-NS Mergers (Gamma Ray Bursts)
//#include <iomanip>
//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
//#include <fstream>	//in-/output to/from file
#include <limits>	//provides numerical limits
#include <random>	//provides the random generator functions
#include <time.h>	//provides the function time
#include "ComBinElib/ComBinElib.h" // local header file

#include<omp.h>		//provides parallel programming routines

//using namespace std;

/* constants*/
const double Msun=1.9885E+30;	//from PDG 2014: solar Mass in kg
const double Rsun=6.9551E+8;	//from PDG 2014: solar Radius in m
const double yr=31558149.8;	//from PDG 2014: (sidereal) year in s
const double Lsun=3.828E+26*yr*yr*yr/(Msun*Rsun*Rsun);	//from PDG 2014: solar Luminosity in kg*m^2/(s^3) converted to Msun*Rsun^2/(yr^3)
const double G=6.67384E-11*Msun*yr*yr/(Rsun*Rsun*Rsun);	//from PDG 2014: graviational constant in m^3/(kg*s^2) converted to Rsun^3/(Msun*yr^2)
const double c=299792458.0*yr/Rsun;	//from PDG 2014: speed of light in m/s converted to Rsun/yr
const double MH=1.672621777E-27/Msun;	//from PDG 2014: ionized hydrogen mass(=proton mass) in kg converted to Msun
//const double MHe=4.002603/1.007825*MH;		//ionized helium mass in Msun
const double sigmaT=0.6652458734*1.0E-28/(Rsun*Rsun);	//from PDG 2014: Thomson cross section of electrons in barn(=1.0E-28m^2) converted to Rsun^2
const double cgsEnergy=1.0E-7*yr*yr/(Msun*Rsun*Rsun);	//convertion: 1erg = cgsEnergy*Msun*Rsun^2/(yr^2)
const double kms=1.0E+3*yr/Rsun;	//conversion: 1km/s = kms*Rsun/yr
const double day=((23*60+56)*60+4.09053)/yr;	//from PDG 2014: (sidereal) day in s converted to yr
const double au=1.495978707E+11/Rsun;	//from PDG 2014: astronomical unit in m converted to Rsun
const double pc=3.08567758149E+16/Rsun;	//from PDG 2014: parsec in m converted to Rsun

const double rsmax=2.0*G*1000.0/(c*c);	//Schwarzschild radius of a 1000Msun star
const double tmax=1.0E+10;	// maximal time to emit graviational waves in yr
const double WDmax=1.48;	//maximal mass of a non-rotating white dwarf (Chandrasekhar limit)	//1.37 -> 1.48
const double NSmax=2.5;		//maximal mass of a neutron star

/*global variables*/
long number=1000000;		// number of calculated systems
bool screen=false;		// additional output
bool debug=false;		// additional output for debuging
bool separateevolution=false;	// evolve the stars one after the other while the companion is unevolved (for comparison with the old programm)
bool galacticintegration=false;	// calculate the movement of the system in the galaxy
double accuracy=1.0e-10;	// minimal accuracy of calculation
double alphaRLO=0.20;		// 20% direct wind mass loss from donor star during RLO 
double beta_const=0.75;		// all material transfered to secondary (non-deg.) star is not accreted
double Gamma=2.0;		// gamma^2 = a_r/a where a_r is the radius of a circum stellar ring
double delta=0.0;		// no mass loss to circum stellar ring
double alphaCE=0.5;		// 50% efficiency of converting E_orb to E_kin env. during CE
double lambda_const=1.0;	// constant lambda for CE: try 1.0, 0.5 and 0.2
double qlimit=2.5;		// somewhat arbitrary (both MS and comp.obj. accretors)
double Mhece=-3.3;		// if RLO and M_he < Mhece => CE! - see Dewi & Pols (2002): CE for M_he in (2.8,3.2) and no CE in (3.6,6.4), no CE M_he>3.3 or M_he>3.8 in close orbits
double kickv=-1.0;		// kick velocity in km/s; if negative take random value: if equal 0 no kick, else asymmetric SN (3 comp. Maxwellian)
double alphaTH=0.5;		// factor to interpolate between lambda with/without thermal and radiation energy (0.0-1.0: interpolate, -1: take old lambdatable, else assume viral equilibrium)
char output='M';		// output format: T (tiny) or M (massive)
//int addran=0;			// number which is added to new random seed
bool single=false;		//single: run only one system with maximal values of Mp, Ms, a and e=0

int main(int argc, char *argv[])	//argc: count of arguments, argv: argument values
{
  ifstream input;		//enable reading from input file
  bool singledebug=false;	// additional output for debuging in single run
  bool wanted=false;		//search for a wanted system
// SET PHYSICS PARAMETERS:
  double Mp_max=100.0;		// 24 for NSNS
  double Mp_min=4.0;		// 10.0 to form NS/BH (else 0.8)
  double Ms_max=Mp_max;		// 15 for NSNS; Mp_max (or 2.0 to make BMSPs)
  double Ms_min=1.0;		// 3.0 to form NS or BH since RLO is needed (else CE)
  double a_max=10000.0;		// maximal semi-major-axis in Rsun
  double a_min=2.0;		// minimal semi-major-axis in Rsun
  double Mp=0.0, Ms=0.0, a=0.0;	// primary mass, secondary mass, semi-major-axis of one system
  double e=0.0, metal=0.02, vp=0.0, vs=0.0;	// eccentricity(circular), metallicity(solar value), primary rotational velocity(non-rotating), secondary rotational velocity(non-rotating)
  double x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, rhostar=0.0;	// galactic position(at center), galactic velocity(at rest), stellar density(field star)
  double savevx=0.0, savevy=0.0, savevz=0.0;	// save of galactic velocity change
  double stepsize=-20.0;		//step size in the grid

  t_system system;		//structure with all system parameters

/*  double age, t_thermal;	//old code
  double Mp[11]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};	//old code
  double Ms[11]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};	//old code
  double a[11]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; 	//old code
  double times[11]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};	//old code 
  double q,Rl,ecc,vCM1,vCM2,t_merge,maxdist,R2;	//old code
  double tau,ytry2,w1,w2,theta1,phi1,theta2,phi2;	//old code*/
  long seed=0;			//current seed of the random number generator
  long seeds[7]={-8809843958428921796,-4468887128648617210,-5430439732881999393,6759374438094590665,0,0,0};	//current seeds of the random number generators{0,0,0,0,0,0}
  long discard[7]={0,0,0,0,0,0,0};	//number of random numbers to discard{0,0,0,0,0,0}
  static long discard_mp=0;	//discards of primary mass generator
  static long discard_q=0;	//discards of mass ratio generator
  static long discard_a=0;	//discards of semi-major axis generator
  static long discard_kick=0;	//discards of kick parameter generator
  static long discard_phase=0;	//discards of phase parameter generator
  static long discard_eccentricity=0;	//discards of eccentricity generator
  static long discard_galactic=0;	//discards of galactic generator
//  long saveseeds[7]={0,0,0,0,0,0,0};	//old values of seeds to check if same system are calculated again
//  static long discardout[7]={discard[0],discard[1],discard[2],discard[3],discard[4],discard[5],discard[6]};	//for output discard
  bool nseed=false;		//negative input seed
/*  long* saveseed=(long *)malloc(sizeof(long));	//old values of seed to check if same system are calculated again
  long nsave=0;			//number of entries in the arrays above
  long incaddran=-1;		//addran should increased if incaddran is 0 or larger*/
  long i,j;			//loop variables
  long done=0;			//number of already calculated systems
  time_t t1, t2;		//times
  int n;			//last entry of the system
  int parallel=1;		//number of parallel threads (0: get maximal number automatically)
  int galintegrator=1, galpotential=0;	//flag for integrator and potential
  static char phasehistory[200]="0";	//number which tells the formation canal
/*  int Heq,SNq;	//old code

  int inform[11]={0,0,0,0,0,0,0,0,0,0,0};	//old code
  int count[11]={0,0,0,0,0,0,0,0,0,0,0};	//old code*/

  First();	//user interface function

  t1 = time(NULL);		//meassure time of the program
  if (output=='M') cerr << "Start: " << ctime(&t1) << endl;	//output start time of the program

/*  if (seeds==NULL){
    cerr << "#Error: memory allocation failed: seeds" << endl;
    return -1;
  }*/
//  for (j=0;j<numeric_limits<int>::max();j++) if (seeds[j]) cout << "seeds[" << j << "]=true" << endl;
// The random generator and parameters
  if (argc>1){			//check for random seed as input parameter
    seed = -atol(argv[1]);	//convert the first input parameter to the seed/input (see below)
    if (seed==11) singledebug = true;	//input = -11 run a single evolution in debug mode
    if (seed==4){		//input = -4 old version (evolve one star after the other one) with given parameters
      separateevolution = true;	//set evolution to happen first for primary and when this is finished evolve the secondary
      alphaTH = -1.0;		//take old lambda table
      seed = 1;			//evolve a single system
    }
    if (seed==3){		//input = -3 make a grid
      grid = true;
      if (argc>2){
        stepsize = atof(argv[2]);	//read in step size from input arguments
        if (stepsize<0) stepsize = 1.0/stepsize;	//negative stepsize is number of steps in log
        if (argc>3){
          seed = -atol(argv[3]);	//convert the thrid input parameter to the seed/input (see below)
        }
      }
    }
    if (seed==2){		//input = -2 read all values from file data.in
      input.open("data.in");	//open file with the name "data.in"
      if (!input.is_open()) cerr << "#Error: can't open data.in" << endl;
      input.ignore(10000,':');	//ignore the characters until ':', here "NUMBER OF PROCESSORS TO USE (<=0 AUTOMATIC):"
      input >> parallel;	//read in the number of parallel threads
      input.ignore(10000,':');	//ignore the characters until ':', here "SEED:"
      input >> seed;		//read in the seed
/*      input.ignore(10000,':');	//ignore the characters until ':', here "ADDRAN:"
      input >> addran;		//read in the addran*/
      input.ignore(10000,':');	//ignore the characters until ':', here "NUMBER:"
      input >> number;		//read in the number
      input.ignore(10000,':');	//ignore the characters until ':', here "MAXIMUM PRIMARY MASS (MSUN):"
      input >> Mp_max;		//read in the maximal mass of the primary star in Msun
      input.ignore(10000,':');	//ignore the characters until ':', here "MINIMUM PRIMARY MASS (MSUN):"
      input >> Mp_min;		//read in the minimal mass of the primary star in Msun
      input.ignore(10000,':');	//ignore the characters until ':', here "MAXIMUM SECONDARY MASS (MSUN):"
      input >> Ms_max;		//read in the maximal mass of the secondary star in Msun
      input.ignore(10000,':');	//ignore the characters until ':', here "MINIMUM SECONDARY MASS (MSUN):"
      input >> Ms_min;		//read in the minimal mass of the secondary star in Msun
      input.ignore(10000,':');	//ignore the characters until ':', here "MAXIMUM SEMI-MAJOR AXIS (RSUN):"
      input >> a_max;		//read in the maximal semi-major-axis in Rsun
      input.ignore(10000,':');	//ignore the characters until ':', here "MINIMUM SEMI-MAJOR AXIS (RSUN):"
      input >> a_min;		//read in the minimal semi-major-axis in Rsun
      input.ignore(10000,':');	//ignore the characters until ':', here "WIND MASS LOSS DURING RLO (%):"
      input >> alphaRLO;	//read in the direct wind mass loss from donor star during RLO
      input.ignore(10000,':');	//ignore the characters until ':', here "RLO MASS EJECTION PARAMETER (%):"
      input >> beta_const;	//read in the not accreted material during RLO
      input.ignore(10000,':');	//ignore the characters until ':', here "RLO RADIUS PARAMETER OF CIRCUMBINARY TORUS:"
      input >> Gamma;		//read in the not accreted material during RLO
      input.ignore(10000,':');	//ignore the characters until ':', here "RLO MASS-FRACTION TO CIRCUMBINARY TORUS (%):"
      input >> delta;		//read in the material going to circum stellar ring during RLO
      input.ignore(10000,':');	//ignore the characters until ':', here "CE EFFICIENCY PARAMETER (%):"
      input >> alphaCE;		//read in the efficiency of converting E_orb to E_kin env. during CE
      input.ignore(10000,':');	//ignore the characters until ':', here "LAMBDA:"
      input >> lambda_const;	//read in the LAMBDA
      input.ignore(10000,':');	//ignore the characters until ':', here "QLIMIT:"
      input >> qlimit;		//read in the QLIMIT
      input.ignore(10000,':');	//ignore the characters until ':', here "MHECE:"
      input >> Mhece;		//read in the MHECE
      input.ignore(10000,':');	//ignore the characters until ':', here "KICK VELOCITY (KM/S; <0 TAKE RANDOM VALUE):"
      input >> kickv;		//read in the kick velocity
      input.ignore(10000,':');	//ignore the characters until ':', here "LAMBDA INTERPOLATION FACTOR (0.0-1.0 INTERPOLATE; <0.0 OLD TABLES; >1.0 CONSTANT FAKTOR OF LAMBDA_G[2.0 VIRIAL EQULIBRIUM]):"
      input >> alphaTH;		//read in the factor to interpolate between lambda with/without thermal and radiation energy
      input.ignore(10000,':');	//ignore the characters until ':', here "MOTION INTEGRATOR (=0 NO GALACTIC MOTION; =1 RUNGEKUTTA4):"
      input >> galintegrator;	//read in the flag of galactic integrator
      input.ignore(10000,':');	//ignore the characters until ':', here "GALACTIC POTENTIAL (=0 NO GALACTIC MOTION; =1 MW-POTENTIAL BY ALLEN&SANTILLAN 1991):"
      input >> galpotential;	//read in the flag of galactic potential
      input.ignore(10000,':');	//ignore the characters until ':', here "OUTPUT FORMAT (=M MASSIVE OUTPUT; =T TINY OUTPUT):"
      input >> output;		//read in the output format
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL MASS FUNCTION (=1 SALPETER 1955; =2 KROUPA 2008):"
      input >> IMF;		//read in the flag of initial mass function
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL MASS-RATIO DISTRIBUTION (=1 KUIPER 1935; =2 FLAT):"
      input >> qdist;		//read in the flag of initial mass ratio distribution
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL SEMI-MAJOR AXIS DISTRIBUTION (=1 ABT 1983):"
      input >> adist;		//read in the flag of initial semi-major axis distribution
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL ECCENTRICITY DISTRIBUTION (=0 CIRCULAR; =1 THERMAL(HEGGIE 1975); =2 FLAT; =3 FLAT IN ANGULAR MOMENTUM):"
      input >> edist;		//read in the flag of initial eccentricity distribution
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL AGE DISTRIBUTION (=0 ALL AT ZAMS):"
      input >> tdist;		//read in the flag of initial age distribution
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL METALLICITY DISTRIBUTION (=0 SOLAR METALLICITY):"
      input >> Zdist;		//read in the flag of initial metallicity distribution
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL ROTATION DISTRIBUTION (=-1 SYNCHRONISED; =0 NON-ROTATING):"
      input >> rotdist;		//read in the flag of initial stellar spin distribution
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL POSITION IN THE GALAXY (=-1 SUN; =0 AT GALACTIC CENTER):"
      input >> Rdist;		//read in the flag of initial space distribution in a galaxy
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL VELOCITY IN THE GALAXY (=-1 SUN; =0 AT REST):"
      input >> Vdist;		//read in the flag of initial velocity distribution in a galaxy
      input.ignore(10000,':');	//ignore the characters until ':', here "INITIAL STELLAR DENSITY (=0 FIELD STAR):"
      input >> rhodist;		//read in the flag of initial stellar density distribution in a galaxy
      input.close();		//close the input file
      if (seed<0) nseed = true;	//seed was negative before
      else seed = -seed;	//switch seed to negative value;
      if (Ms_max>Mp_max){
        Ms_max = Mp_max;	//secondary mass must be lower than the primary
        cout << "The maximal mass of the secondary is set to the maximal mass of the primary: " << Ms_max << "Msun" << endl;
      }
      if (Mp_max<=Mp_min){
        Mp_min = Mp_max-accuracy;	//minimum mass must be lower than the maximum mass
        cout << "The minimum mass of the primary is set to: " << Mp_min << "Msun" << endl;
      }
      if (Ms_max<=Ms_min){
        Ms_min = Ms_max-accuracy;	//minimum mass must be lower than the maximum mass
        cout << "The minimum mass of the secondary is set to: " << Ms_min << "Msun" << endl;
      }
      if (Ms_min>Mp_min){
        Ms_min = Mp_min;	//secondary mass must be lower than the primary
        cout << "The minimal mass of the secondary is set to the minimal mass of the primary: " << Ms_min << "Msun" << endl;
      }
      if (a_max<=a_min){
        a_min = a_max-accuracy;	//minimum semi-major-axis must be lower than the maximum semi-major-axis
        cout << "The minimum semi-major-axis is set to: " << a_min << "Rsun" << endl;
      }
      alphaRLO /= 100.0;		//convert the per cent input
      beta_const /= 100.0;	//convert the per cent input
      delta /= 100.0;		//convert the per cent input
      alphaCE /= 100.0;		//convert the per cent input
      if ((output!='M')&&(output!='T')){
        output = 'M';	//output format set to stardard
        cout << "The output format is set to: " << output << endl;
      }
    }
    if ((seed==1)||singledebug){	//input = -1 run a single evolution with following parameters
      single = true;		//evolve a single system
/*      Mp = Mp_max;		//set primary mass in Msun
      Ms = Ms_max;		//set secondary mass in Msun
      a = a_max;		//set semi-major-axis in Rsun	//pow(pow((12.0/365.25)/(2.0*M_PI),2)*(G*(Mp_max+Ms_max)),1.0/3.0);
//      kickv = 50.0;		//set kickvelocity in km/s
      e = 0.0;			//set eccentricity: circular
      metal = 0.02;			//set metallicity: solar metallicity
      vp = 0.0;			//set rotation velocity of primary: non-rotating
      vs = 0.0;			//set rotation velocity of secondary: non-rotating
      x = 0.0;			//set galactic x-position to center
      y = 0.0;			//set galactic y-position to center
      z = 0.0;			//set galactic z-position to center
      vx = 0.0;			//set galactic x-velocity at rest
      vy = 0.0;			//set galactic y-velocity at rest
      vz = 0.0;			//set galactic z-velocity at rest
      rhostar = 0.0;		//set stellar density: field star*/
      seed = newseed(3,0,0);	//initialize the random number generator for the kick parameters
      seed = newseed(4,0,0);	//initialize the random number generator for the phase
//      Zdist = -3;
//      seed = 13579;		//set specific seed
      if (argc>2){
        Mp = atof(argv[2]);	//read in primary mass in Msun from input arguments
        if (argc>3){
          Ms = atof(argv[3]);	//read in secondary mass in Msun from input arguments
          if (argc>4){
            a = atof(argv[4]);	//read in semi-major-axis in Rsun from input arguments
            if (argc>5){
              seed = atol(argv[5]);	//read in random seed from input arguments
              if (seed<0) nseed = true;	//seed was negative before
              else seed = -seed;	//switch seed to negative value;
              if (argc>6){
                kickv = atof(argv[6]);	//read in kickvelocity in km/s(<0 is random) from input arguments
              }
            }
          }
        }
      }
    }
  }
  if (seed>=0) seed = -time(NULL);	//if seed is not specified take seed from current time (seconds since 1970)

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
    metal = UserInitial_Z(0.0, 0.0, 0.0, 0.0, 0.0);
  }else{
    cerr << "#Error: No initial metallicity distribution spezified! metallicity set to 0.02" << endl;
    metal = 0.02;
  }

  AfterReadParameters(Mp_max, Mp_min, Ms_max, Ms_min, a_max, a_min, Mp, Ms, a, e, metal, vp, vs, x, y, z, vx, vy, vz, rhostar, parallel, galintegrator, galpotential);	//user interface function

  if (galintegrator==1) integrator = rk4;	//set galactic integrator
  else if (galintegrator!=0){	//check if galactic integrator has a valid value
    integrator = rk4;
    cerr << "#Warning: Invalid integrator spezified: use RungeKutta 4 integrator" << endl;
  }
  if (galpotential==1) potential = MW_potential;	//set galactic potential
  else if (galpotential!=0){	//check if galactic potential has a valid value
    galacticintegration = false;
    cerr << "#Warning: Invalid potential spezified: use no potential" << endl;
  }
  if ((galintegrator==0)||(galpotential==0)) galacticintegration = false;	//check if galactic integration is possible
  else galacticintegration = true;

  PreReadTables();		//user interface function

  readtables(metal, vp>0);	//initializes and read in all tables
  readmergetime();		//initializes and read mergetime table

  AfterReadTables();		//user interface function

  if (debug){			// check tables
    cerr << "nstar " << nstar << endl;
    cerr << "stararray" << endl;
    for(i=0;i<nstar;i++){		//go through all masses in stararray
      cerr << "[" << i << "]: inimass=" << stararray[i].inimass << " TAMS=" << stararray[i].TAMS << endl;	//write out initial mass
      cerr << "[" << i << "]: n=" << stararray[i].n << " m[0]" << stararray[i].m[0] << " t[0]" << stararray[i].t[0] << " r[0]" << stararray[i].r[0] << " cm[0]" << stararray[i].cm[0] << " llum[0]" << stararray[i].llum[0] << " lteff[0]" << stararray[i].lteff[0] << " lambda[0]" << stararray[i].lambda[0] << " cc[0]" << stararray[i].cc[0] << endl;	//write out maximum values
      for(j=1;j<min(stararray[i].n-1,999);j++){
        cerr << "[" << i << "]: j=" << j << " m[j]" << stararray[i].m[j] << " t[j]" << stararray[i].t[j] << " r[j]" << stararray[i].r[j] << " cm[j]" << stararray[i].cm[j] << " llum[j]" << stararray[i].llum[j] << " lteff[j]" << stararray[i].lteff[j] << " lambda[j]" << stararray[i].lambda[j] << " cc[j]" << stararray[i].cc[j] << endl;	//write out track point
//        cerr << " cc[j]-cc[j-1]" << stararray[i].cc[j]-stararray[i].cc[j-1] << endl;
      }
//      if ((stararray[i].TAMS>2)&&(stararray[i].n-1>stararray[i].TAMS)) cerr << "[" << i << "]: n=" << stararray[i].n << " m[TAMS]" << stararray[i].m[stararray[i].TAMS] << " t[TAMS]" << stararray[i].t[stararray[i].TAMS] << " r[TAMS]" << stararray[i].r[stararray[i].TAMS] << " cm[TAMS]" << stararray[i].cm[stararray[i].TAMS] << " llum[TAMS]" << stararray[i].llum[stararray[i].TAMS] << " lteff[TAMS]" << stararray[i].lteff[stararray[i].TAMS] << " lambda[TAMS]" << stararray[i].lambda[stararray[i].TAMS] << " cc[TAMS]" << stararray[i].cc[stararray[i].TAMS] << endl;	//write out track point at TAMS
      cerr << "[" << i << "]: n=" << stararray[i].n << " m[n-1]" << stararray[i].m[stararray[i].n-1] << " t[n-1]" << stararray[i].t[stararray[i].n-1] << " r[n-1]" << stararray[i].r[stararray[i].n-1] << " cm[n-1]" << stararray[i].cm[stararray[i].n-1] << " llum[n-1]" << stararray[i].llum[stararray[i].n-1] << " lteff[n-1]" << stararray[i].lteff[stararray[i].n-1] << " lambda[n-1]" << stararray[i].lambda[stararray[i].n-1] << " cc[n-1]" << stararray[i].cc[stararray[i].n-1] << endl;		//write out last track point
    }
    if (stararray2!=stararray){
      cerr << "nstar2 " << nstar2 << endl;
      cerr << "stararray2" << endl;
      for(i=0;i<nstar2;i++){		//go through all masses in stararray2
        cerr << "[" << i << "]: inimass=" << stararray2[i].inimass << " TAMS=" << stararray2[i].TAMS << endl;	//write out initial mass
        cerr << "[" << i << "]: n=" << stararray2[i].n << " m[0]" << stararray2[i].m[0] << " t[0]" << stararray2[i].t[0] << " r[0]" << stararray2[i].r[0] << " cm[0]" << stararray2[i].cm[0] << " llum[0]" << stararray2[i].llum[0] << " lteff[0]" << stararray2[i].lteff[0] << " lambda[0]" << stararray2[i].lambda[0] << " cc[0]" << stararray2[i].cc[0] << endl;	//write out maximum values
        for(j=1;j<min(stararray2[i].n-1,999);j++){
          cerr << "[" << i << "]: j=" << j << " m[j]" << stararray2[i].m[j] << " t[j]" << stararray2[i].t[j] << " r[j]" << stararray2[i].r[j] << " cm[j]" << stararray2[i].cm[j] << " llum[j]" << stararray2[i].llum[j] << " lteff[j]" << stararray2[i].lteff[j] << " lambda[j]" << stararray2[i].lambda[j] << " cc[j]" << stararray2[i].cc[j] << endl;	//write out track point
        }
//        if ((stararray2[i].TAMS>2)&&(stararray2[i].n-1>stararray2[i].TAMS)) cerr << "[" << i << "]: n=" << stararray2[i].n << " m[TAMS]" << stararray2[i].m[stararray2[i].TAMS] << " t[TAMS]" << stararray2[i].t[stararray2[i].TAMS] << " r[TAMS]" << stararray2[i].r[stararray2[i].TAMS] << " cm[TAMS]" << stararray2[i].cm[stararray2[i].TAMS] << " llum[TAMS]" << stararray2[i].llum[stararray2[i].TAMS] << " lteff[TAMS]" << stararray2[i].lteff[stararray2[i].TAMS] << " lambda[TAMS]" << stararray2[i].lambda[stararray2[i].TAMS] << " cc[TAMS]" << stararray2[i].cc[stararray2[i].TAMS] << endl;	//write out track point at TAMS
        cerr << "[" << i << "]: n=" << stararray2[i].n << " m[n-1]" << stararray2[i].m[stararray2[i].n-1] << " t[n-1]" << stararray2[i].t[stararray2[i].n-1] << " r[n-1]" << stararray2[i].r[stararray2[i].n-1] << " cm[n-1]" << stararray2[i].cm[stararray2[i].n-1] << " llum[n-1]" << stararray2[i].llum[stararray2[i].n-1] << " lteff[n-1]" << stararray2[i].lteff[stararray2[i].n-1] << " lambda[n-1]" << stararray2[i].lambda[stararray2[i].n-1] << " cc[n-1]" << stararray2[i].cc[stararray2[i].n-1] << endl;		//write out last track point
      }
    }
    cerr << "nhe1 " << nhe1 << endl;
    cerr << "hestararray1" << endl;
    for(i=0;i<nhe1;i++){		//go through all masses in hestararray1
      cerr << "[" << i << "]: inimass=" << hestararray1[i].inimass << " TAMS=" << hestararray1[i].TAMS << endl;	//write out initial mass
      cerr << "[" << i << "]: n=" << hestararray1[i].n << " m[0]" << hestararray1[i].m[0] << " t[0]" << hestararray1[i].t[0] << " r[0]" << hestararray1[i].r[0] << " cm[0]" << hestararray1[i].cm[0] << " llum[0]" << hestararray1[i].llum[0] << " lteff[0]" << hestararray1[i].lteff[0] << " lambda[0]" << hestararray1[i].lambda[0] << " cc[0]" << hestararray1[i].cc[0] << endl;	//write out maximum values
      for(j=1;j<min(hestararray1[i].n-1,999);j++){
        cerr << "[" << i << "]: j=" << j << " m[j]" << hestararray1[i].m[j] << " t[j]" << hestararray1[i].t[j] << " r[j]" << hestararray1[i].r[j] << " cm[j]" << hestararray1[i].cm[j] << " llum[j]" << hestararray1[i].llum[j] << " lteff[j]" << hestararray1[i].lteff[j] << " lambda[j]" << hestararray1[i].lambda[j] << " cc[j]" << hestararray1[i].cc[j] << endl;	//write out track point
      }
//      if ((hestararray1[i].TAMS>2)&&(hestararray1[i].n-1>hestararray1[i].TAMS)) cerr << "[" << i << "]: n=" << hestararray1[i].n << " m[TAMS]" << hestararray1[i].m[hestararray1[i].TAMS] << " t[TAMS]" << hestararray1[i].t[hestararray1[i].TAMS] << " r[TAMS]" << hestararray1[i].r[hestararray1[i].TAMS] << " cm[TAMS]" << hestararray1[i].cm[hestararray1[i].TAMS] << " llum[TAMS]" << hestararray1[i].llum[hestararray1[i].TAMS] << " lteff[TAMS]" << hestararray1[i].lteff[hestararray1[i].TAMS] << " lambda[TAMS]" << hestararray1[i].lambda[hestararray1[i].TAMS] << " cc[TAMS]" << hestararray1[i].cc[hestararray1[i].TAMS] << endl;	//write out track point at TAMS
      cerr << "[" << i << "]: n=" << hestararray1[i].n << " m[n-1]" << hestararray1[i].m[hestararray1[i].n-1] << " t[n-1]" << hestararray1[i].t[hestararray1[i].n-1] << " r[n-1]" << hestararray1[i].r[hestararray1[i].n-1] << " cm[n-1]" << hestararray1[i].cm[hestararray1[i].n-1] << " llum[n-1]" << hestararray1[i].llum[hestararray1[i].n-1] << " lteff[n-1]" << hestararray1[i].lteff[hestararray1[i].n-1] << " lambda[n-1]" << hestararray1[i].lambda[hestararray1[i].n-1] << " cc[n-1]" << hestararray1[i].cc[hestararray1[i].n-1] << endl;		//write out last track point
    }
    cerr << "nhe " << nhe << endl;
    cerr << "hestararray" << endl;
    for(i=0;i<nhe;i++){		//go through all masses in hestararray
      cerr << "[" << i << "]: n=" << hestararray[i].n << " m[0]" << hestararray[i].m[0] << " t[0]" << hestararray[i].t[0] << " lteff[0]" << hestararray[i].lteff[0] << " llum[0]" << hestararray[i].llum[0] << " lr[0]" << hestararray[i].lr[0] << " mw[0]" << hestararray[i].mw[0] << " cm[0]" << hestararray[i].cm[0] << endl;	//write out first track point
      cerr << "[" << i << "]: n=" << hestararray[i].n << " m[1]" << hestararray[i].m[1] << " t[1]" << hestararray[i].t[1] << " lteff[1]" << hestararray[i].lteff[1] << " llum[1]" << hestararray[i].llum[1] << " lr[1]" << hestararray[i].lr[1] << " mw[1]" << hestararray[i].mw[1] << " cm[1]" << hestararray[i].cm[1] << endl;	//write out second track point
      cerr << "[" << i << "]: n=" << hestararray[i].n << " m[n-1]" << hestararray[i].m[hestararray[i].n-1] << " t[n-1]" << hestararray[i].t[hestararray[i].n-1] << " lteff[n-1]" << hestararray[i].lteff[hestararray[i].n-1] << " llum[n-1]" << hestararray[i].llum[hestararray[i].n-1] << " lr[n-1]" << hestararray[i].lr[hestararray[i].n-1] << " mw[n-1]" << hestararray[i].mw[hestararray[i].n-1] << " cm[n-1]" << hestararray[i].cm[hestararray[i].n-1] << endl;	//write out last track point
    }
    cerr << "lambdaarray" << endl;	//print lambda array
    cerr << "n=" << lambdaarray.n << " m[0]" << lambdaarray.m[0] << " r[0]" << lambdaarray.r[0] << " lambda[0]" << lambdaarray.lambda[0] << endl;	//write out first track point
    cerr << "n=" << lambdaarray.n << " m[1]" << lambdaarray.m[1] << " r[1]" << lambdaarray.r[1] << " lambda[1]" << lambdaarray.lambda[1] << endl;	//write out second track point
    cerr << "n=" << lambdaarray.n << " m[n-1]" << lambdaarray.m[lambdaarray.n-1] << " r[n-1]" << lambdaarray.r[lambdaarray.n-1] << " lambda[n-1]" << lambdaarray.lambda[lambdaarray.n-1] << endl;	//write out last track point
  }

  if (single){			//evolve a single system
    cout << "Single" << endl;
    number = 1;
  }
  if (grid){			//run a grid starting at the maximal values
    cout << "Grid" << endl;
    screen = true;		//enable screen output
    Mp = Mp_max;		//set primary mass in Msun
    Ms = Ms_max;		//set secondary mass in Msun	//Ms = fmin(Ms_max,0.999*Mp);
    a = a_max;			//set semi-major-axis in Rsun
//    a_min = 10.0;		//set lower limit of semi-major-axis in Rsun
//    number = 1000000;		//maximal number of grid points
  }
  if (number==1) screen = true;	//enable screen output for a single run
  if (screen && single) cout << "seed=" << seed << endl;	//write out the random seed
  if (!nseed) seed = abs(seed);	//convert to positve/initial seed
  SetSeed(seed,discard,7);	//user interface function to manipulate discard
  discard_mp = discard[0];	//discards of primary mass generator
  discard_q = discard[1];	//discards of mass ratio generator
  discard_a = discard[2];	//discards of semi-major axis generator
  discard_kick = discard[3];	//discards of kick parameter generator
  discard_phase = discard[4];	//discards of phase parameter generator
  discard_eccentricity = discard[5];	//discards of eccentricity generator
  discard_galactic = discard[6];	//discards of galactic generator
  if (debug) cerr << "discard={" << discard_mp << "," << discard_q << "," << discard_a << "," << discard_kick << "," << discard_phase << "," << discard_eccentricity << "," << discard_galactic << "}" << endl;
  for (i=0; i<7; i++){
    seeds[i] = seed+i;		//set seeds
    if (newseed(i,seeds[i],discard[i])!=seeds[i]){	//initializes the random number generator with seed
      seeds[i] = newseed(i,0,0);	//if seed setting fails get a new seed
      cerr << "#Warning: generate new seed(i): " << seeds[i] << endl;
    }
//    saveseeds[i]=seeds[i];	//save seeds
  }
//  srand(seed);		//initializes the random number generator with seed

  initoutfiles();		//prepares the output files to write
/*  if (data){			//write parameters to data.out
    dataout << "#random seed\tnumber of systems\tprimary mass(max,min)\tsecondary mass(max,min)\tsemi-major-axis(max,min)\twind loss in RLO\tefficency of CE\tget lambda with alpha_th\tmass ratio limit(RLO-CE)\tmass limit(RLO-CE)\tkick velocity" << endl;
    dataout << seed << "\t" << number << "\t(" << Mp_max << "," << Mp_min << ")\t(" << Ms_max << "," << Ms_min << ")\t(" << a_max << "," << a_min << ")\t" << alphaRLO << "\t" << alphaCE << "\t" << alphaTH << "\t" << qlimit << "\t" << Mhece << "\t" << kickv << endl;
  }*/

  if (plot) inithist(tmerge,85,1.0e-5,1.0e+80,true,10);	//initialise histogramm

  if (parallel<=0) parallel = omp_get_num_procs();	//gets number of available processors
  omp_set_num_threads(parallel);	//set number of threads
  if (parallel>1) cout << "run " << parallel << "threads" << endl;

  PreMainLoop();
  t1 = time(NULL);		//meassure time of the program
/*main loop to run (number) systems*/
#pragma omp parallel for private(system, i, j, t2, n) firstprivate(wanted, Mp, Ms, a, e, metal, vp, vs, x, y, z, rhostar, seed, seeds, discard_mp, discard_q, discard_a, discard_kick, discard_phase, discard_eccentricity, discard_galactic, phasehistory, data, dist, grid, plot) shared(count_mp, count_q, count_a, count_kick, count_phase, count_eccentricity, count_galactic) reduction(+ : cnotrackupdate, cRLO, cRLOA, cRLOB, cRLOBB, cCE, cCEA, cCEB, cCEBB, cMerger, cMergerRLO, cMergerCE, cMergerSN, csupernova, cNS, cBH, cplanetarynebula, cWD, cdestroyed, cmerged, cexploded, cunknown, cfinal, cWDWD, cWDNS, cWDBH, cNSWD, cNSNS, cNSBH, cBHWD, cBHNS, cBHBH, cgw, cgwWDWD, cgwWDNS, cgwWDBH, cgwNSWD, cgwNSNS, cgwNSBH, cgwBHWD, cgwBHNS, cgwBHBH)
//screen, debug are threadprivate
//*cformation, *cgwformation, tmerge inside critical ranges
//addran, discardout at the moment not used
  for (i=1; i<=number; i++){
    if ((i==ToDebug(i))||singledebug){	//debug specified system
      debug = true;		//set debug output
      screen = true;		//set screen output
    }else{
      debug = false;		//unset debug output
    }
    if (screen) cout << endl;
#pragma omp critical
{
    if (((number>99)&&(done%(number/100)==0))||((number>100000000)&&(done%1000000==0))||((grid)&&(done%100000==0))){	//progress: output at each per cent or every million systems
      t2 = time(NULL);		//meassure time of the program
      writedataout(parallel, seed, Mp_max, Mp_min, Ms_max, Ms_min, a_max, a_min, galintegrator, galpotential, done);	//writes the current data to data.out
      if (output=='M') cerr << " Number: " << done << " speed:" << done/difftime(t2,t1) << "systems/second [" << done/(difftime(t2,t1)+1.0) << "," << done/(difftime(t2,t1)-1.0) << "]" << endl;	//output of number and speed
    }else if (number<100){
      if (output=='M') cerr << " Number: " << done << endl;	//output of number
    }
    done++;
}	//end critical
/*create new system/change parameters*/
    if (single){		//run single system
      initsystem(system, Mp, Ms/Mp, a, e, 0.0, metal, vp, vs, x, y, z, vx, vy, vz, rhostar);	//create new system with specified parameters (at ZAMS: t=0)
    }else if (grid){		//run grid
      initsystem(system, 0.001*round(1000*Mp), round(1000*Ms)/round(1000*Mp), 0.001*round(1000*a), e, 0.0, metal, vp, vs, x, y, z, vx, vy, vz, rhostar);	//create new system with specified parameters (at ZAMS: t=0)
      seed = newseed(3,i,0);	//set new seed for kick parameters
      seed = newseed(4,i,0);	//set new seed for phase
//      srand(i); 		//set new seed
      if (stepsize<0){		//logarithmic grid
        if (Ms/pow(10,-stepsize)>=Ms_min) Ms /= pow(10,-stepsize);	//use log-steps per dex in secondary mass
        else{
          if (Mp/pow(10,-stepsize)>=Mp_min) Mp /= pow(10,-stepsize);	//use log-steps per dex in primary mass
          else{
            if (a/pow(10,-stepsize)>=a_min) a /= pow(10,-stepsize);	//use log-steps per dex in semi-major-axis
            else if (i!=number){	//grid finished
//              cout << "grid end: " << i << endl;
              i = number;		//reduce number of systems to calculate
            }
            Mp = Mp_max;		//get back to maximum in primary mass
          }
          Ms = Mp;		//get back to maximum in secondary mass	//Ms = 0.999*Mp;
        }
      }else{		//linear grid
        if (Ms-stepsize>=Ms_min) Ms -= stepsize;	//use linear-steps per dex in secondary mass
        else{
          if (Mp-stepsize>=Mp_min) Mp -= stepsize;	//use linear-steps per dex in primary mass
          else{
            if (a-stepsize>=a_min) a -= stepsize;	//use linear-steps per dex in semi-major-axis
            else if (i!=number){	//grid finished
//              cout << "grid end: " << i << endl;
              i = number;		//reduce number of systems to calculate
            }
            Mp = Mp_max;		//get back to maximum in primary mass
          }
          Ms = Mp;		//get back to maximum in secondary mass	//Ms = 0.999*Mp;
        }
      }
      if (i==ToDebug(i)){		//debug specified system	//1133546755
        debug = true;		//set debug output
        screen = true;		//set screen output
      }else{
        debug = false;		//unset debug output
      }
    }else{
      initial(system, Mp_max, Mp_min, Ms_max, Ms_min, a_max, a_min, metal);	//create new random system (circular, at ZAMS, with non-rotating stars, at solar metallicity)
    }
    if ((i<10)||(i==0)) screen = true; else screen = debug;	//screen output for a specified amount of systems in a monte carlo run	//number-5	//5244
/*    if (number>1000000){
      if ((saveseeds[0]==seeds[0])&&(saveseeds[1]==seeds[1])&&(saveseeds[2]==seeds[2])&&(i>1)){	//check for repetion of random seed sequence
        for (j=0;j<4;j++){
          seeds[j]=newseed(j,seeds[j]+1,0);	//set new random seeds
          saveseeds[j]=seeds[j];		//save seeds
        }
        cerr << "#Warning new initial-seeds={" << seeds[0] << "," << seeds[1] << "," << seeds[2] << "}" << endl;
      }
      if ((saveseeds[3]==seeds[3])&&(i>1)&&(discardout[3]!=discard[3]+count_kick)){	//check for repetion of random seed sequence
        seeds[3]=newseed(3,seeds[3]+1,0);	//set new random seed
        saveseeds[3]=seeds[3];			//save seed
        cerr << "#Warning new kick-seed=" << seeds[3] << endl;
      }	UPDATE*/
/*      for (j=0;j<nsave;j++){	//check if system occured before
        if (saveseed[j]==seed){
          if (incaddran>=0) cerr << "#Error: Seed twice found: seed=" << seed << " saveseed[" << j << "]=" << saveseed[j] << " saveseed[" << incaddran << "]=" << saveseed[incaddran] << " nsave=" << nsave << endl;
          incaddran = j;
        }
      }
      if (incaddran>=0){
        addran++;	//increment addran
        cout << "addran=" << addran << endl;
        nsave = 0;	//start array new
        incaddran = -1;
      }
      if (i%1000000==0){
        //change reserved memory
        saveseed = (long *)realloc(saveseed,(nsave+1)*sizeof(long));
        //check if memory allocation fails
        if (saveseed==NULL) cerr << "#Error: memory reallocation failed: saveseed" << endl;
        //save current value
        saveseed[nsave] = seed;
        //increase number of entries in the array
        nsave++;
      }
    }*/
    if (screen && (!single)) cout << "i=" << i << endl;	//write out the number of the system
    if (screen && (output=='T')) cout << "Initial:" << endl << "  Mp Ms a: " << system.prim.m[system.n-1] << " " << system.sec.m[system.n-1] << " " << system.a[system.n-1] << endl;	//write out the initial conditions of the system
    newphase(system);		//create new phase to start evolution
/*evolve the system*/
    while ((system.phase[system.n-1]>=0)&&(system.phase[system.n-1]<99)&&(system.n<100)){	//evolve the system till the end or the system evolved through more than 100 phases
      evolve(system);
//      system.phase[system.n-1]=99;
//      cout << system.formation << endl;
      if ((system.prim.stage[system.n-1]>2)&&(system.sec.track.t[2]==9.0e+99)&& separateevolution){	//start evolution of secondary in seperate evolution
        //free memory
        free(system.sec.track.m);
        free(system.sec.track.t);
        free(system.sec.track.r);
        free(system.sec.track.cm);
        free(system.sec.track.llum);
        free(system.sec.track.lteff);
        free(system.sec.track.lambda);
        free(system.sec.track.cc);
//        cout << "memory check: &savetrack=" << &savetrack << " &system.sec.track=" << &system.sec.track << endl;
//        cout << "&savetrack.(n,m,t,r,cm,llum,lteff,lambda,cc)=\t\t" << &savetrack.n << ",\t" << savetrack.m << ",\t" << savetrack.t << ",\t" << savetrack.r << ",\t" << savetrack.cm << ",\t" << savetrack.llum << ",\t" << savetrack.lteff << ",\t" << savetrack.lambda << ",\t" << savetrack.cc << "\n&system.sec.track.(n,m,t,r,cm,llum,lteff,lambda,cc)=\t" << &system.sec.track.n << ",\t" << system.sec.track.m << ",\t" << system.sec.track.t << ",\t" << system.sec.track.r << ",\t" << system.sec.track.cm << ",\t" << system.sec.track.llum << ",\t" << system.sec.track.lteff << ",\t" << system.sec.track.lambda << ",\t" << system.sec.track.cc << endl;
        system.sec.track = savetrack;	//"copy" saved track
        system.sec.rmax = savermax;	//copy saved rmax
//        cout << "memory check: &savetrack=" << &savetrack << " &system.sec.track=" << &system.sec.track << endl;
//        cout << "&savetrack.(n,m,t,r,cm,llum,lteff,lambda,cc)=\t\t" << &savetrack.n << ",\t" << savetrack.m << ",\t" << savetrack.t << ",\t" << savetrack.r << ",\t" << savetrack.cm << ",\t" << savetrack.llum << ",\t" << savetrack.lteff << ",\t" << savetrack.lambda << ",\t" << savetrack.cc << "\n&system.sec.track.(n,m,t,r,cm,llum,lteff,lambda,cc)=\t" << &system.sec.track.n << ",\t" << system.sec.track.m << ",\t" << system.sec.track.t << ",\t" << system.sec.track.r << ",\t" << system.sec.track.cm << ",\t" << system.sec.track.llum << ",\t" << system.sec.track.lteff << ",\t" << system.sec.track.lambda << ",\t" << system.sec.track.cc << endl;
//        cout << endl << "system.sec.track.n=" << system.sec.track.n << " system.sec.track.(m,t,r,cm,llum,lteff,lambda)=\n[(" << system.sec.track.m[0] << "," << system.sec.track.t[0] << "," << system.sec.track.r[0] << "," << system.sec.track.cm[0] << "," << system.sec.track.llum[0] << "," << system.sec.track.lteff[0] << "," << system.sec.track.lambda[0] << "," << system.sec.track.cc[0] << ")";
        for (j=0;j<system.sec.track.n;j++){	//add time of primary evolution to secondary track
          system.sec.track.t[j] += system.prim.t[system.n-1];
//          cout << ", (" << system.sec.track.m[j] << "," << system.sec.track.t[j] << "," << system.sec.track.r[j] << "," << system.sec.track.cm[j] << "," << system.sec.track.llum[j] << "," << system.sec.track.lteff[j] << "," << system.sec.track.lambda[j] << "," << system.sec.track.cc[j] << ")";
        }
//        cout << "]" << endl;
        // circularization
        system.a[system.n-2] *= 1.0-pow(system.e[system.n-2],2);	//calculate semi-major axis with angular momentum conservation
        system.e[system.n-2] = 0.0;	//set eccentricity to 0
        cout << "circularization: a=" << system.a[system.n-1] << " e=" << system.e[system.n-1] << " => a=" << system.a[system.n-2] << " e=" << system.e[system.n-2] << endl;
      }
    }
/*write out system*/
    n = system.n-1;		//get last index
    if (system.formation==0){	//system gets destroyed before one star reached finial stage
      system.formation = 10*system.prim.stage[n]+system.sec.stage[n];
    }else if (system.formation<10){	//system gets destroyed before second star reached finial stage
      if ((system.formation==system.prim.stage[n])||(system.prim.stage[n]<0)) system.formation = 10*system.formation+system.sec.stage[n];	//primary reached finial stage
      else if ((system.formation==system.sec.stage[n])||(system.sec.stage[n]<0)) system.formation = 10*system.formation+system.prim.stage[n];	//secondary reached finial stage
      else{			//no star reached finial stage
        cerr << "#Warning: wrong formation value" << endl;
        screen = true;		//enable screen output
      }
    }else if (system.formation>99){	//system had more than two supernovea
      if (system.formation%10==system.prim.stage[n]) system.formation = 10*system.sec.stage[n]+system.prim.stage[n];
      else if (system.formation%10==system.sec.stage[n]) system.formation = 10*system.prim.stage[n]+system.sec.stage[n];
      else cerr << "#Error: system.formation=" << system.formation << " system.prim.stage[" << n << "]=" << system.prim.stage[n] << " system.sec.stage[" << n << "]=" << system.sec.stage[n] << endl;
    }
    if (screen){		//output
      if (output=='M') cout << " end";
      cout << endl;
    }
    if (system.n==100) screen = true;	//the system does not get to an endstage
    if ((system.prim.stage[n]==3)&&(system.prim.m[n]>WDmax)){	//the primary is a too massive white dwarf
      cerr << "#Error: too massive white dwarf(primary):" << system.prim.m[n] << endl;
      screen = true;		//enable screen output
    }
    if ((system.sec.stage[n]==3)&&(system.sec.m[n]>WDmax)){	//the secondary is a too massive white dwarf
      cerr << "#Error: too massive white dwarf(secondary):" << system.sec.m[n] << endl;
      screen = true;		//enable screen output
    }
    if ((system.prim.stage[n]==4)&&(system.prim.m[n]>NSmax)){	//the primary is a too massive neutron star
      cerr << "#Error: too massive neutron star(primary):" << system.prim.m[n] << endl;
      screen = true;		//enable screen output
    }
    if ((system.sec.stage[n]==4)&&(system.sec.m[n]>NSmax)){	//the secondary is a too massive neutron star
      cerr << "#Error: too massive neutron star(secondary):" << system.sec.m[n] << endl;
      screen = true;		//enable screen output
    }
//    if ((system.prim.stage[n]==4)&&(system.sec.stage[n]==4)&&(system.phase[n]==99)&&(system.phase[n-2]%10!=5)) wanted = true;
//    if ((system.prim.stage[n]==4)&&(system.sec.stage[n]==3)&&(system.phase[n]==99)&&(system.formation==34)) wanted = true;	//WDNS system
//    if ((system.phase[n-1]==4)&&((system.phase[n-2]==12)||(system.phase[n-2]==21))) wanted = true;	//merger after RLO
//    if ((system.prim.stage[n]==4)&&(system.sec.stage[n]==4)&&(system.phase[n]==99)&&(system.prim.stage[n-3]==3)) wanted = true;	//NSNS from WD implosion
//    if ((system.prim.stage[n]==4)&&(system.sec.stage[n]==4)&&(system.phase[n]==99)&&(system.prim.m[0]<7)) wanted = true;	//too low mass to get NSNS
//    if ((system.prim.stage[n]==4)&&(system.sec.stage[n]==4)&&(system.phase[n]==99)&&(system.sec.m[0]<3.4)) wanted = true;	//too low mass to get NSNS
//    if ((system.prim.stage[n]==4)&&(system.sec.stage[n]==4)&&(system.phase[n]==99)&&(system.phase[max(n-1,0)]==25)&&(system.phase[max(n-3,0)]==23)) wanted = true;
//    if ((system.phase[n]==99)&&((system.prim.m[n]>50.0)||(system.sec.m[n]>50.0))) wanted = true;	//very massive BHs
//    if ((system.phase[n]==99)&&(system.prim.stage[n]==5)&&(system.sec.stage[n]==5)&&(system.prim.m[n]<10.0)&&(system.sec.m[n]<10.0)) wanted = true;	//low mass BHs
//    if ((system.formation==45)&&(system.phase[n]==99)) wanted = true;	//NSBH
//    if ((system.phase[n]==99)&&(system.e[n]<=0.10)&&(system.e[n]>=0.08)&&(system.P[n]/day<=11)&&(system.P[n]/day>=9)&&(abs(system.prim.stage[n]-system.sec.stage[n])<=1)&&(system.prim.stage[n]+system.sec.stage[n]<=8)&&(system.prim.stage[n]+system.sec.stage[n]>=7)) wanted = true;	//WDNS system
    if (wanted) screen = true;
//    if ((system.prim.stage[n]==5)&&(system.sec.stage[n]==5)&&(system.phase[n]>90)&&(system.prim.m[n]>35)&&(system.prim.m[n]<37)&&(system.sec.m[n]>28)&&(system.sec.m[n]<30)) screen = true;	//LIGO BHBHmerger
    adddata(system,i);		//add system to statistics
/*    if ((wanted)&&(system.tgw>tmax)){
      screen = false;
      wanted = false;
    }*/
    if (snapshots &&(galintegrator!=0)&&(galpotential!=0)){
      if ((system.phase[n]==99)&&((system.formation==34)||(system.formation==43))&&((system.prim.m[n]<0.3)||(system.sec.m[n]<0.3))&&(system.prim.t[n]<1.3e+10)){	// (NS+WD[<0.3Msun])
        galacticintegration = true;
        if (system.prim.m[n]<0.3){
          fprintf(snapshot0,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", system.prim.t[0], system.x[0], system.y[0], system.z[0], system.vx[0], system.vy[0], system.vz[0], system.prim.m[0], system.prim.stage[0], system.sec.m[0], system.sec.stage[0]);
        }else{
          fprintf(snapshot0,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", system.prim.t[0], system.x[0], system.y[0], system.z[0], system.vx[0], system.vy[0], system.vz[0], system.sec.m[0], system.sec.stage[0], system.prim.m[0], system.prim.stage[0]);
        }
        fflush(snapshot0);	//update file
        if(system.prim.t[n]<1.0e+11){
          newphase(system);		//create new phase to start evolution
          n=system.n-1;
          system.phase[n] = 100;
          system.prim.t[n] = 1.0e+11;
          system.sec.t[n] = 1.0e+11;
        }
      }else{
        galacticintegration = false;
      }
    }
    if ((system.phase[n]==99)&&(system.tgw<tmax)){
      newphase(system);		//create new phase to start evolution
      n=system.n-1;
      system.phase[n] = 94;
      system.a[n] = -1.0e+99;
      system.P[n] = -1.0e+99;
      system.e[n] = -1.0e+99;
      system.rp[n] = -1.0e+99;
      system.rs[n] = -1.0e+99;
      system.prim.t[n] = system.prim.t[n-1]+system.tgw;
      system.sec.t[n] = system.prim.t[n-1]+system.tgw;
    }
    sprintf(phasehistory,"0");	//start phase history with phase 0
    system.stagechange = 0;	//use stagechange as index for movement calculation
    for (j=1; j<system.n; j++){
      if (j<n) sprintf(phasehistory,"%s%d",phasehistory,system.phase[j]);	//add next phase to phase history
    //write to distributionout # open/close file in statistics
      if ((distributionout!=NULL)&&(system.phase[j]==99)&&(system.prim.stage[j]>3)&&(system.sec.stage[j]>3)){	//write specified system (NS/BS binaries) to distribution output
        fprintf(distributionout,"%12.6e\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%14.7e\t%4d\t%12.5f\t%10.6f  \t%10.6f\t%9.3f\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%12.5f\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.tgw, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi, system.prim.SN.vsys, system.sec.SN.type, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi, system.sec.SN.vsys);
/*      if ((distributionout!=NULL)&&(system.phase[j]==99)&& snapshots && galacticintegration){	//write snapshot systems to distribution output
        fprintf(distributionout,"%12.6e\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%14.7e\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%4d\t%9.3f\t%10.6f  \t%10.6f\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.tgw, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi, system.sec.SN.type, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi);*/
/*      if ((distributionout!=NULL)&&(system.phase[j]==99)&&((system.formation==34)||(system.formation==44))){	//write specified system to distribution output
        fprintf(distributionout,"%12.6e\t%9.3f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%14.7e\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%13.6f\n", system.prim.t[j], system.P[j]/day, system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.tgw, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi, system.sec.SN.type, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi, sqrt(system.vx[j]*system.vx[j]+system.vy[j]*system.vy[j]+system.vz[j]*system.vz[j]));*/
/*      if ((distributionout!=NULL)&&(system.P[j]/day>=0.1)&&(system.phase[j]!=-1)&&(system.phase[j-1]==15)&&(system.prim.stage[j]==4)&&(system.sec.m[j]<1.0)){	//write specified system to distribution output (low mass MS star as companion of J1755-2550	&&(system.P[j]/day<=1000.0)
        fprintf(distributionout,"%12.6e\t%9.3f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%14.7e\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%13.6f\n", system.prim.t[j], system.P[j]/day, system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.tgw, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi, system.sec.SN.type, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi, sqrt(system.vx[j]*system.vx[j]+system.vy[j]*system.vy[j]+system.vz[j]*system.vz[j]));*/
/*      if ((system.phase[j]==15)&&(system.prim.stage[j-1]==2)&&(system.prim.m[j]>2.99)&&(system.prim.m[j]<3.01)&&(system.prim.stage[j+1]==4)&&(system.prim.r[j]<system.rp[j])&&(system.sec.stage[j-1]==0)&&(system.sec.r[j-1]<system.rs[j])){	//write specified system to distribution output	//&&(system.phase[j+1]==0)
        fprintf(distributionout,"%12.6e\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%10.4f\t%9.4f\t%10.6f\t%10.6f\t%13.4f\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j-1], system.sec.stage[j-1],phasehistory, system.a[j+1], system.e[j+1], system.prim.m[j+1], system.sec.m[j+1], system.a[j+1]*(1.0-pow(system.e[j+1],2)));*/
/*      if ((distributionout!=NULL)&&(system.phase[j]==99)&&(system.prim.stage[j]==4)&&(system.sec.stage[j]==4)){	//write specified system to distribution output
        fprintf(distributionout,"%12.6e\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%14.7e\t%13.6f\t%13.6f\t%13.6f\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%13.6f\t%13.6f\t%13.6f\t%4d\t%9.3f\t%10.6f  \t%10.6f\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.tgw, system.prim.SN.Hecoremass, system.prim.SN.COcoremass, system.prim.SN.remnantmass, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi, system.sec.SN.Hecoremass, system.sec.SN.COcoremass, system.sec.SN.remnantmass, system.sec.SN.type, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi);*/
        fflush(distributionout);	//update file
//        if (system.a[j]<1.0) screen=true;	//enable screen output
      }
      if (galacticintegration){
        if ((system.vx[j]!=0.0)||(system.vy[j]!=0.0)||(system.vz[j]!=0.0)){	//save kickvelocity
          savevx = system.vx[j];
          savevy = system.vy[j];
          savevz = system.vz[j];
          movesystem(system,system.stagechange,j);	//moves the system in space
          system.stagechange = j;	//save where last movement was calculated
          system.vx[j] += savevx;
          system.vy[j] += savevy;
          system.vz[j] += savevz;
        }else if (j<=n){	//no kick occured
          movesystem(system,system.stagechange,j);	//moves the system in space
          system.stagechange = j;	//save where last movement was calculated
        }else{	//no movement calculated -> no output
          system.x[j] = -1.0e+99;
          system.y[j] = -1.0e+99;
          system.z[j] = -1.0e+99;
          system.vx[j] = -1.0e+99;
          system.vy[j] = -1.0e+99;
          system.vz[j] = -1.0e+99;
        }
      }
      if ((system.rhostar[0]==0.0)&&(system.rhostar[j]==0.0)) system.rhostar[j] = -1.0e+99;	//field star -> no output
      if ((system.prim.m[j]==0)||(system.sec.m[j]==0)) cerr << "#Error: zero mass" << endl;	//check for zero mass
      if (system.M[j-1]+accuracy<system.M[j]){	//... whether the primary gained mass
        cerr << "#Error: mass increase of the system: dm=" << system.M[j]-system.M[j-1] << "Msun" << endl;
        screen = true;	//enable screen output
        debug = true;	//enable debugging
//        i=number;
      }
      if (system.phase[j-1]==0){	//ckeck in wind mass loss phase ...
        if (system.prim.m[j-1]+accuracy<system.prim.m[j]){	//... whether the primary gained mass
          cerr << "#Error: wind mass increase of primary: dm=" << system.prim.m[j]-system.prim.m[j-1] << "Msun" << endl;
          debug = true;	//enable debugging
        }
        if (system.sec.m[j-1]+accuracy<system.sec.m[j]){	//... whether the secondary gained mass
          cerr << "#Error: wind mass increase of secondary: dm=" << system.sec.m[j]-system.sec.m[j-1] << "Msun" << endl;
          debug = true;	//enable debugging
        }
      }
    }
#pragma omp critical
{
    if ((output=='M')&&(((screen)&&((debug)||(wanted)||(system.phase[n]==99)||(i>100)||(i==1))))){	//write out the evolution of the system	//||((system.prim.stage[n]!=3)&&(system.sec.stage[n]!=3)&&(system.phase[n]==99)&&(i<100000))
      if (i==1) cout << endl << "seeds={" << seeds[0] << "," << seeds[1] << "," << seeds[2] << "," << seeds[3] << "," << seeds[4] << "," << seeds[5] << "," << seeds[6] << "}";
      for (j=0; j<system.n; j++){	//go through the evolution
        if (j==0){		//initial output
          if (!single) cout << endl << "discard={" << discard_mp << "," << discard_q << "," << discard_a << "," << discard_kick << "," << discard_phase << "," << discard_eccentricity << "," << discard_galactic << "}";
          cout << endl << "Initial ";
        }
        if ((j==n)||(system.phase[j]==4)) cout << "Final ";		//final output
        cout << "system [" << j << "]:";	//write out count of phase
        outphasename(system.phase[j]);	//write out number and name of phase
        if (system.phase[j]==15) cout << "(type=" << system.prim.SN.type << "): kickvelocity=" << system.prim.SN.w << "km/s kickangles{theta,phi}={" << system.prim.SN.theta << "," << system.prim.SN.phi << "}degree" << endl;	//write primary SN-kick
        if (system.phase[j]==25) cout << "(type=" << system.sec.SN.type << "): kickvelocity=" << system.sec.SN.w << "km/s kickangles{theta,phi}={" << system.sec.SN.theta << "," << system.sec.SN.phi << "}degree" << endl;	//write secondary SN-kick
        cout << " sys :";	//write out system parameters
        if (system.M[j]!=-1.0e+99) cout << " M=" << system.M[j] << "Msun";	//write out total mass
        if ((system.qp[j]!=-1.0e+99)||(system.qs[j]!=-1.0e+99)) cout << " q=1/" << system.qp[j] << "=" << system.qs[j];	//write out mass ratio
        if ((system.a[j]!=-1.0e+99)&&(system.phase[j]!=-1)) cout << " a=" << system.a[j] << "Rsun=" << system.a[j]/au << "AU";	//write out semi-major-axis
        if ((system.P[j]!=-1.0e+99)&&(system.phase[j]!=-1)) cout << " P=" << system.P[j] << "yr=" << system.P[j]/day << "days";	//write out orbital period
        if (((system.rp[j]!=-1.0e+99)||(system.rs[j]!=-1.0e+99))&&(system.phase[j]!=-1)) cout << " roche{prim/sec}={" << system.rp[j] << "/" << system.rs[j] << "}Rsun";	//write out Roche-lobe radii
        if (system.e[j]!=-1.0e+99) cout << " e=" << system.e[j];	//write out eccentricity
        if (galacticintegration){
          cout << endl << "       ";
          if ((system.x[j]!=-1.0e+99)||(system.y[j]!=-1.0e+99)||(system.z[j]!=-1.0e+99)) cout << " position{x,y,z}={" << system.x[j] << "," << system.y[j] << "," << system.z[j] << "}pc";	//write out galactic position
          if ((system.vx[j]!=-1.0e+99)||(system.vy[j]!=-1.0e+99)||(system.vz[j]!=-1.0e+99)) cout << " velocity{x,y,z}={" << system.vx[j] << "," << system.vy[j] << "," << system.vz[j] << "}km/s";	//write out galactic velocity
          if ((system.rhostar[j]!=-1.0e+99)&&(rhodist!=0)) cout << " rhostar=" << system.rhostar[j] << "/pc^3";	//write out stellar density
        }
        cout << endl << " prim:";	//write out parameters of the primary
        if (system.prim.m[j]!=-1.0e+99) cout << " m=" << system.prim.m[j] << "Msun";	//write out mass
        if (system.prim.t[j]!=-1.0e+99) cout << " t=" << system.prim.t[j] << "yr";	//write out age
        if (system.prim.r[j]!=-1.0e+99) cout << " r=" << system.prim.r[j] << "Rsun";	//write out radius
        if (system.prim.cm[j]!=system.prim.m[j]) cout << " core mass(cm)=" << system.prim.cm[j] << "Msun";	//write out core mass
        if (system.prim.llum[j]!=-1.0e+99) cout << " lg(lum/Lsun)=" << system.prim.llum[j];	//write out luminosity
        if (system.prim.lteff[j]!=-1.0e+99) cout << " lg(T_eff/K)=" << system.prim.lteff[j];	//write out effective temperature
        if ((system.prim.lambda[j]!=-1.0e+99)&&(alphaTH!=-1.0)) cout << " lambda=" << system.prim.lambda[j];	//write out lambda
        if (system.prim.cc[j]!=system.prim.cm[j]) cout << " carbon core mass(cc)=" << system.prim.cc[j] << "Msun";	//write out carbon core mass
        if ((system.prim.omega[j]!=-1.0e+99)&&((system.prim.omega[j]!=0.0)||(j==0))) cout << " omega=" << system.prim.omega[j] << "/yr";	//write out angular velocity
        outstagename(system.prim.stage[j],true);	//write out number and name of stage
        cout << " sec :";	//write out parameters of the secondary
        if (system.sec.m[j]!=-1.0e+99) cout << " m=" << system.sec.m[j] << "Msun";	//write out mass
        if (system.sec.t[j]!=-1.0e+99) cout << " t=" << system.sec.t[j] << "yr";	//write out age
        if (system.sec.r[j]!=-1.0e+99) cout << " r=" << system.sec.r[j] << "Rsun";	//write out radius
        if (system.sec.cm[j]!=system.sec.m[j]) cout << " core mass(cm)=" << system.sec.cm[j] << "Msun";	//write out core mass
        if (system.sec.llum[j]!=-1.0e+99) cout << " lg(lum/Lsun)=" << system.sec.llum[j];	//write out luminosity
        if (system.sec.lteff[j]!=-1.0e+99) cout << " lg(T_eff/K)=" << system.sec.lteff[j];	//write out effective temperature
        if ((system.sec.lambda[j]!=-1.0e+99)&&(alphaTH!=-1.0)) cout << " lambda=" << system.sec.lambda[j];	//write out lambda
        if (system.sec.cc[j]!=system.sec.cm[j]) cout << " carbon core mass(cc)=" << system.sec.cc[j] << "Msun";	//write out carbon core mass
        if ((system.sec.omega[j]!=-1.0e+99)&&((system.sec.omega[j]!=0.0)||(j==0))) cout << " omega=" << system.sec.omega[j] << "/yr";	//write out angular velocity
        outstagename(system.sec.stage[j],true);	//write out number and name of stage

        if (system.phase[j]==4) j = n;	//skip post-merger phase
        cout << endl;
      }
      cout << "Further final parameters of system number " << i << ":" << endl;	//write out further parameters
      cout << " sys : " << system.formation;	//write out formation order
      outstagename(system.formation/10,false);	//write out number and name of final stage of first remnant
      outstagename(system.formation%10,true);	//write out number and name of final stage of second remnant
/*      if (system.prim.RLO||system.sec.RLO||system.prim.HeRLO||system.sec.HeRLO){	//check for Roche-lobe-overflow
        if (system.prim.RLO||system.sec.RLO) cout << " Roche-lobe-overflow happend";
        if ((system.prim.RLO||system.sec.RLO)&&(system.prim.HeRLO||system.sec.HeRLO)) cout << ",";
        if (system.prim.HeRLO||system.sec.HeRLO) cout << " Roche-lobe-overflow happend form Helium star";
        cout << endl;
      }*/
      cout << " prim:" << " n_track=" << system.prim.track.n << "\tsec :" << " n_track=" << system.sec.track.n << endl;	//write out count of track points
    }
    if ((wanted)&&(screen)) i = number;
//    debug = true;
    if (debug){	//write out track data for debuging
      cerr << "primary:" << endl << "j\tmass(Msun)\tage(yr)\tradius(Rsun)\tcore mass(Msun)\tlg(lum/Lsun)\tlg(T_eff/K)\tlambda\tcarbon core(Msun)" << endl;	//write head of primary track-table
      for(j=0;j<system.prim.track.n;j++){	//write track-table of the primary
        cerr << j << "\t" << system.prim.track.m[j] << "\t" << system.prim.track.t[j] << "\t" << system.prim.track.r[j] << "\t" << system.prim.track.cm[j] << "\t" << system.prim.track.llum[j] << "\t" << system.prim.track.lteff[j] << "\t" << system.prim.track.lambda[j] << "\t" << system.prim.track.cc[j] << "\t" << endl;
      }
      cerr << "secondary:" << endl << "j\tmass(Msun)\tage(yr)\tradius(Rsun)\tcore mass(Msun)\tlg(lum/Lsun)\tlg(T_eff/K)\tlambda\tcarbon core(Msun)" << endl;	//write head of secondary track-table
      for(j=0;j<system.sec.track.n;j++){	//write track-table of the secondary
        cerr << j << "\t" << system.sec.track.m[j] << "\t" << system.sec.track.t[j] << "\t" << system.sec.track.r[j] << "\t" << system.sec.track.cm[j] << "\t" << system.sec.track.llum[j] << "\t" << system.sec.track.lteff[j] << "\t" << system.sec.track.lambda[j] << "\t" << system.sec.track.cc[j] << "\t" << endl;
      }
    }

}	//end critical

    AfterEvolution(system);	//user interface function

    destroysystem(system);	//delete data of current system to continue with the next one
//    seed = rand();	//get next random seed
//    srand(seed);	//set next random seed
//    for (j=0;j<5;j++) seeds[j] = newseed(j,0,0);	//get/set next random seed
#pragma omp critical
{
    discard_mp = discard[0]+count_mp;	//write count into discard_mp
    discard_q = discard[1]+count_q;	//write count into discard_q
    discard_a = discard[2]+count_a;	//write count into discard_a
    discard_kick = discard[3]+count_kick;	//write count into discard_kick
    discard_phase = discard[4]+count_phase;	//write count into discard_phase
    discard_eccentricity = discard[5]+count_eccentricity;	//write count into discard_eccentricity
    discard_galactic = discard[6]+count_galactic;	//write count into discard_galactic
    seed = seeds[0];	//take seed
}	//end critical
//    if (i%5000000==0) addran++;	//increase number of loops in random seed sequence
  }
  t2 = time(NULL);		//meassure time of the program
  AfterMainLoop();
  writedataout(parallel, seed, Mp_max, Mp_min, Ms_max, Ms_min, a_max, a_min, galintegrator, galpotential, done);	//writes the current data to data.out
  if (output=='M') cerr << " Number: " << done << " speed:" << done/difftime(t2,t1) << "systems/second [" << done/(difftime(t2,t1)+1.0) << "," << done/(difftime(t2,t1)-1.0) << "]" << endl;	//output of number and speed
  if (plot){
    histout(tmerge, false, (char *)"Tmerge");	//write the histogramm data to file
    freehist(tmerge);	//free the histogramm memory
  }
  closeoutfiles();	//close the output files
  destroytables();	//free the memory for the stellar evolution tables
  destroymergetime();	//free the memory for the integration table

  if (cnotrackupdate>0) cerr << "#Warning: " << cnotrackupdate << "track(s) not updated, because star out of grid" << endl;	//give summarized warning of skipped track updates

  t2 = time(NULL);	//get end time of the calculation
  if (output=='M') cerr << "End: " << ctime(&t2) << endl;	//output end time

  return 0;		//program finishes without error report
}

