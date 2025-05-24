// user interface
// this file contains functions which are called in the ComBinE code
// only this file should be changed by users of ComBinE to extract data or make modifications

#include "ComBinElib.h" // local header file

/*define new histograms*/
//t_hist initialmass;
//t_hist initialperiod;
FILE* distfile1=NULL;	//enable writing to a file
FILE* distfile2=NULL;	//enable writing to a file
FILE* distfile3=NULL;	//enable writing to a file
FILE* distfile4=NULL;	//enable writing to a file
FILE* distfile5=NULL;	//enable writing to a file
FILE* distfile6=NULL;	//enable writing to a file
FILE* distfile7=NULL;	//enable writing to a file
bool snapshots=false;	//snapshots
FILE* snapshot0=NULL;	//enable writing to a file
FILE* snapshot1=NULL;	//enable writing to a file
FILE* snapshot2=NULL;	//enable writing to a file
FILE* snapshot3=NULL;	//enable writing to a file
FILE* snapshot4=NULL;	//enable writing to a file
FILE* snapshot5=NULL;	//enable writing to a file
FILE* snapshot6=NULL;	//enable writing to a file
FILE* snapshot7=NULL;	//enable writing to a file
FILE* snapshot8=NULL;	//enable writing to a file
FILE* snapshot9=NULL;	//enable writing to a file
FILE* snapshot10=NULL;	//enable writing to a file
FILE* snapshot11=NULL;	//enable writing to a file
FILE* snapshot12=NULL;	//enable writing to a file
FILE* snapshot13=NULL;	//enable writing to a file
FILE* Bep=NULL;		//enable writing to a file
FILE* Bes=NULL;		//enable writing to a file
FILE* Be=NULL;		//enable writing to a file
FILE* nathan=NULL;	//enable writing to a file
FILE* deM=NULL;		//enable writing to a file

void First(){	//this function is called directly after the variable declaration in the main
//  dist = true;
  return;
}

void AfterReadParameters(double& Mp_max, double& Mp_min, double& Ms_max, double& Ms_min, double& a_max, double& a_min, double& Mp, double& Ms, double& a, double& e, double& metal, double& vp, double& vs, double& x, double& y, double& z, double& vx, double& vy, double& vz, double& rhostar, int& parallel, int& galintegrator, int& galpotential){	//this function is called after the parameters are read in from the command line or data.in
/*calulation parameters*/
  number = number;		//change number of systems to calculate
  screen = false;		//enable/disable screen-output
  debug = false;		//enable/disable debug-output
  separateevolution = false;	//enable/disable to evolve the stars one after the other while the companion is unevolved (for comparison with the old programm)
  accuracy = 1.0e-10;		//change minimal accuracy of calculation; tested only for default: 1.0e-10
  alphaRLO = alphaRLO;		//change direct wind mass loss from donor star during RLO 
  beta_const = beta_const;	//change not accreted material
  Gamma = Gamma;		//change gamma^2 = a_r/a where a_r is the radius of a circum stellar ring
  delta = delta;		//change mass loss transfered to circum stellar ring
  alphaCE = alphaCE;		//change efficiency of converting E_orb to E_kin during CE
  lambda_const = lambda_const;	//change constant lambda
  qlimit = qlimit;		//change mass ratio limit to distinguish between RLO and CE
  Mhece = Mhece;		//change He star mass limit to distinguish between RLO and CE
  kickv = kickv;		//change kick velocity in km/s; if negative take random value: if equal 0 no kick, else asymmetric SN (3 comp. Maxwellian)
  alphaTH = alphaTH;		//change factor to interpolate between lambda with/without thermal and radiation energy (0.0-1.0: interpolate, -1: take old lambdatable, else assume viral equilibrium)
  output = output;		//change output format: T (tiny) or M (massive)
  parallel = parallel;		//change number of parallel threads
  galintegrator = galintegrator;	//change galactic integrator
  galpotential = galpotential;	//change galactic potential
/*parameters for a population synthesis run*/
  Mp_max = Mp_max;	//change maximal primary mass in Msun
  Mp_min = Mp_min;	//change minimal primary mass in Msun
  Ms_max = Ms_max;	//change maximal secondary mass in Msun
  Ms_min = Ms_min;	//change minimal secondary mass in Msun
  a_max = a_max;	//change maximal semi-major-axis in Rsun
  a_min = a_min;	//change minimal semi-major-axis in Rsun
/*parameters for a single run*/
  if (Mp==-1.0) Mp = Mp_max;		//set primary mass in Msun
  if (Ms==-1.0) Ms = Ms_max;		//set secondary mass in Msun
  if (a==-1.0) a = a_max;		//set semi-major-axis in Rsun
  e = 0.0;		//set eccentricity: circular
  metal = metal;	//set metallicity: Milky Way metallicity(~0.009), LMC(~0.005), SMC(~0.002), IZw18(~0.0002)
  vp = vp;		//set rotation velocity of primary: non-rotating
  vs = vs;		//set rotation velocity of secondary: non-rotating
  x = 0.0;		//set galactic x-position to center
  y = 0.0;		//set galactic y-position to center
  z = 0.0;		//set galactic z-position to center
  vx = 0.0;		//set galactic x-velocity at rest
  vy = 0.0;		//set galactic y-velocity at rest
  vz = 0.0;		//set galactic z-velocity at rest
  rhostar = 0.0;	//set stellar density: field star
  return;
}

void PreReadTables(){	//this function is called directly before the tables are read in
  debug = false;	//enable/disable debug-output when tables are read
  return;
}

void AfterReadTables(){	//this function is called after the tables are read in
  debug = false;	//enable/disable debug-output to print an overview of the read tables
  return;
}

void SetSeed(long& seed, long* discard, int n){	//this function is called directly before the random seeds are set; n is the number of entries in discard
  debug = false;
  seed = seed;	//change seed value
/*change the discard values to reproduce a specific system*/	//50303333,50303333,50303333,35778845,137556,0
  if (n>0) discard[0] = 0;	//discards of primary mass generator	//50253834
  if (n>1) discard[1] = 0;	//discards of mass ratio generator	//50253834
  if (n>2) discard[2] = 0;	//discards of semi-major axis generator	//50253834
  if (n>3) discard[3] = 0;	//discards of kick parameter generator	//35743237
  if (n>4) discard[4] = 0;	//discards of phase parameter generator	//137428
  if (n>5) discard[5] = 0;	//discards of eccentricity generator
  if (n>6) discard[6] = 0;	//discards of galactic generator
  if (n>7) discard[7] = 0;	//discards of rotation generator
  return;
}

long ToDebug(long i){	//this function returns which system is printed with debug-output
/*expamples for return: 0 no system; 1(,...,number) first system, i all systems*/
  return 0;	//number of system with debug-output
}

void PreMainLoop(){	//this function is called before the main loop of calculating the systems starts
/*use this function to create new histograms, fill them in the function AfterEvolution and write it out and delete it in the function AfterMainLoop*/
//  inithist(initialmass,90,0.1,100.0,true,1);	//initialise histogramm
//  inithist(initialperiod,90,1.0E-3,1.0E+6,true,1);	//initialise histogramm
//  distfile1=fopen("HMXB.csv","w");	// open file with name "HMXB.csv" to write
  if (distfile1!=NULL){
    fprintf(distfile1,"#current\t         \t       \t          \t          \t          \t          \tinitial   \t       \t          \t          \tevolution     \tfirst SN\n");
    fprintf(distfile1,"#t in yr\ta in Rsun\t      e\tMp in Msun\tMs in Msun\tRp in Rsun\tRs in Rsun\t a in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tphases\ttype\tw in km/s\ttheta in deg\tphi in deg\n");
//    fprintf(distfile1,"#initial   \t          \t         \tkick     \t            \t          \tPost-SN   \t          \t         \t       \t            \tevolution    \n");
//    fprintf(distfile1,"#Mp in Msun\tMs in Msun\ta in Rsun\tw in km/s\ttheta in deg\tphi in deg\tMp in Msun\tMs in Msun\ta in Rsun\t      e\tvsys in km/s\tmass transfer\n");
  }
//  distfile2=fopen("pre-SN.csv","w");	// open file with name "pre-SN.csv" to write
  if (distfile2!=NULL){
    fprintf(distfile2,"#SN-progenitor\t             \t             \t     \t             \t             \t    \t         \t            \t          \tcompanion\t     \tevolution\n");
    fprintf(distfile2,"#M(H) in Msun \tM(He) in Msun\tM(CO) in Msun\tstage\tM_ini in Msun\tM(NS) in Msun\ttype\tw in km/s\ttheta in deg\tphi in deg\tM in Msun\tstage\tphases   \n");
  }
//  distfile3=fopen("Galactic.csv","w");	// open file with name "Galactic.csv" to write
  if (distfile3!=NULL){
    fprintf(distfile3,"#current\t          \t          \t          \t          \tinitial   \t       \t          \t          \tevolution     \tpost evolution\tprimary SN\t   \t            \t          \tsecondary SN\t \t            \t          \tgalactic distance in pc\tgalactic velocity in km/s\n");
    fprintf(distfile3,"#t in yr\tMp in Msun\tMs in Msun\tRp in Rsun\tRs in Rsun\t a in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tphases\tt_merge in yr \ttype\tw in km/s\ttheta in deg\tphi in deg\ttype\tw in km/s\ttheta in deg\tphi in deg\t  initial      final   \t   initial       final   \n");
  }
//  distfile4=fopen("LB-1.csv","w");	// open file with name "LB-1.csv" to write
  if (distfile4!=NULL){
    fprintf(distfile4,"#current\t         \t       \t          \t          \t          \t          \tinitial   \t       \t          \t          \tevolution     \tfirst SN\n");
    fprintf(distfile4,"#t in yr\ta in Rsun\t      e\tMp in Msun\tMs in Msun\tRp in Rsun\tRs in Rsun\t a in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tphases\ttype\tw in km/s\ttheta in deg\tphi in deg\n");
  }
//  distfile5=fopen("BHBH_SN.csv","w");	// open file with name "BHBH_SN.csv" to write
  if (distfile5!=NULL){
    fprintf(distfile5,"#pre 2nd SN\t         \t       \t          \t          \t      \t2nd SN     \t           \t           \t         \t            \t          \t            \tpost      \tevolution    \t          \tinitial   \t       \t          \t          \n");
    fprintf(distfile5,"# r in Rsun\ta in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tCHe in Msun\tCCO in Msun\tMBH in Msun\tw in km/s\ttheta in deg\tphi in deg\tvsys in km/s\tt_GW in yr\tt_delay in yr\tphases    \t a in Rsun\t      e\tMp in Msun\tMs in Msun\n");
  }
//  distfile6=fopen("WD+He_merge.csv","w");	// open file with name "WD+He_merge.csv" to write
  if (distfile6!=NULL){
    fprintf(distfile6,"#He+WD       \t          \t       \t          \t          \t      \t          \tinitial   \t       \t          \t          \tevolution\n");
    fprintf(distfile6,"# age in yr  \t a in Rsun\t   e   \tMp in Msun\tMs in Msun\tstages\tt_GW in yr\t a in Rsun\t   e   \tMp in Msun\tMs in Msun\tphases\n");
  }
  distfile7=fopen("ELMWD+MS.csv","w");	// open file with name "ELMWD+MS.csv" to write
  if (distfile7!=NULL){
    fprintf(distfile7,"#ELMWD+MS    \t          \t          \t          \t          \t      \tinitial   \t          \t          \tevolution\n");
    fprintf(distfile7,"# age in yr  \t a in Rsun\tMc in Msun\tM2 in Msun\tR2 in Rsun\tstages\t a in Rsun\tMp in Msun\tMs in Msun\tmerge\tphases\n");
  }
  if(snapshots){
    snapshot0 =fopen("snapshot0.csv","w");	// open file with name "snapshot0.csv" to write
    snapshot1 =fopen("snapshot1.csv","w");	// open file with name "snapshot1.csv" to write
    snapshot2 =fopen("snapshot2.csv","w");	// open file with name "snapshot2.csv" to write
    snapshot3 =fopen("snapshot3.csv","w");	// open file with name "snapshot3.csv" to write
    snapshot4 =fopen("snapshot4.csv","w");	// open file with name "snapshot4.csv" to write
    snapshot5 =fopen("snapshot5.csv","w");	// open file with name "snapshot5.csv" to write
    snapshot6 =fopen("snapshot6.csv","w");	// open file with name "snapshot6.csv" to write
    snapshot7 =fopen("snapshot7.csv","w");	// open file with name "snapshot7.csv" to write
    snapshot8 =fopen("snapshot8.csv","w");	// open file with name "snapshot8.csv" to write
    snapshot9 =fopen("snapshot9.csv","w");	// open file with name "snapshot9.csv" to write
    snapshot10=fopen("snapshot10.csv","w");	// open file with name "snapshot10.csv" to write
    snapshot11=fopen("snapshot11.csv","w");	// open file with name "snapshot11.csv" to write
    snapshot12=fopen("snapshot12.csv","w");	// open file with name "snapshot12.csv" to write
    snapshot13=fopen("snapshot13.csv","w");	// open file with name "snapshot13.csv" to write
  }
  if (snapshot0!=NULL){
    fprintf(snapshot0,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot0,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot1!=NULL){
    fprintf(snapshot1,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot1,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot2!=NULL){
    fprintf(snapshot2,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot2,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot3!=NULL){
    fprintf(snapshot3,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot3,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot4!=NULL){
    fprintf(snapshot4,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot4,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot5!=NULL){
    fprintf(snapshot5,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot5,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot6!=NULL){
    fprintf(snapshot6,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot6,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot7!=NULL){
    fprintf(snapshot7,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot7,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot8!=NULL){
    fprintf(snapshot8,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot8,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot9!=NULL){
    fprintf(snapshot9,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot9,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot10!=NULL){
    fprintf(snapshot10,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot10,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot11!=NULL){
    fprintf(snapshot11,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot11,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot12!=NULL){
    fprintf(snapshot12,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot12,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
  if (snapshot13!=NULL){
    fprintf(snapshot13,"#age in yr\t     galactic position in pc    \t    galactic velocity in km/s   \t      WD      \t   companion  \n");
    fprintf(snapshot13,"#    t    \t     x          y          z    \t    vx         vy         vz    \tm in Msun type\tm in Msun type\n");
  }
//  Bep = fopen("Bep.txt","w");
//  Bes = fopen("Bes.txt","w");
//  Be = fopen("Be.txt","w");
//  nathan = fopen("nathan.txt","w");
//  deM = fopen("deMx.txt","w");
  return;
}

double UserInitial_Mp(double Mp_max, double Mp_min){	//here you can define your own initial primary mass function
  const double alpha = 0.1*IMF+1.0;	//using IMF as single powerlaw with power alpha: N(m) ~ m^alpha + 1.0 from integration
  static double kIMF = (pow(Mp_max,alpha)-pow(Mp_min,alpha));	//normalize IMF to 1

  //use IMF as single powerlaw
  return pow(pow(Mp_min,alpha)+kIMF*ranv(0),1.0/alpha);	//generate random mass from IMF
}

double UserInitial_q(double Ms_max, double Ms_min, double Mp){	//here you can define your own initial mass ratio function
  return Ms_max/Mp;
}

double UserInitial_a(double a_max, double a_min, double Mp, double Ms){	//here you can define your own initial semi-major axis function
  return a_max;
}

double UserInitial_e(double e_max, double Mp, double Ms, double a){	//here you can define your own initial eccentricity
  return 0.0;
}

double UserInitial_t(double Mp, double Ms, double a, double e){	//here you can define your own initial age
  return 0.0;
}

double UserInitial_Z(double Mp, double Ms, double a, double e, double t){	//here you can define your own initial metallicity
  return 0.02;
}

void UserInitial_vrot(double& vp, double& vs, double Mp, double Ms, double a, double e, double t, double metal){	//here you can define your own initial rotational spin
  vp = 0.0;
  vs = 0.0;
  return;
}

void UserInitial_R(double& x, double& y, double& z, double Mp, double Ms, double a, double e, double t, double metal, double vp, double vs){	//here you can define your own initial position in a galactic potential
  x = 0.0;
  y = 0.0;
  z = 0.0;
  return;
}

void UserInitial_V(double& vx, double& vy, double& vz, double Mp, double Ms, double a, double e, double t, double metal, double vp, double vs, double x, double y, double z){	//here you can define your own initial velocity in a galactic potential
  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  return;
}

double UserInitial_rho(double Mp, double Ms, double a, double e, double t, double metal, double vp, double vs, double x, double y, double z, double vx, double vy, double vz){	//here you can define your own initial stellar density
  return 0.0;
}

void AfterEvolution(t_system& system){	//this function is called after the system is evolved, system contains all the evolution information of the system, see ComBinElib.h for the structure information of the type t_system
  long j,js,jp,jj;
  int mem,st1,st2,cas;
  jp=1;js=1;
  double t;
  double rahier;
  double alter,lu1,ro1,ra1,om1,ma1,lu2,ro2,ra2,om2,ma2,M0,Om0,Om;
  static char phasehistory[200]="0";	//number which tells the formation canal
//  addtohist(initialmass,system.prim.m[0],0);	//add data to histogramm
//  addtohist(initialperiod,system.P[0],0);	//add data to histogramm
  sprintf(phasehistory,"%d",system.phase[0]);	//start phase history with first phase
  for (j=1; j<system.n; j++){
    if (j<system.n) sprintf(phasehistory,"%.189s%d",phasehistory,system.phase[j]);	//add next phase to phase history
    if ((distfile1!=NULL)&&(system.phase[j]>=0)&&(system.M[j]>10.0)&&(system.phase[j-1]==15)&&(system.sec.stage[j]<3)){	//write specified system to distfile1 output
      fprintf(distfile1,"%12.6e\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%4d\t%9.3f\t%10.6f  \t%10.6f\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi);
      fflush(distfile1);	//update file
    }else if ((distfile1!=NULL)&&(system.phase[j]>=0)&&(system.M[j]>10.0)&&(system.phase[j-1]==25)&&(system.prim.stage[j]<3)){	//write specified system to distfile2 output
      fprintf(distfile1,"%12.6e\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%4d\t%9.3f\t%10.6f  \t%10.6f\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.sec.SN.type, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi);
      fflush(distfile1);	//update file
    }
    if ((distfile3!=NULL)&&(system.phase[j]==94)&&(system.prim.stage[j]>3)&&(system.sec.stage[j]>3)){	//write specified system to distfile3 output
      fprintf(distfile3,"%12.6e\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%14.7e\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%4d\t%9.3f\t%10.6f  \t%10.6f\t%11.3f %11.3f\t%12.3f %12.3f\n", system.prim.t[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.tgw, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi, system.sec.SN.type, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi, sqrt(system.x[0]*system.x[0]+system.y[0]*system.y[0]+system.z[0]*system.z[0]), sqrt(system.x[j]*system.x[j]+system.y[j]*system.y[j]+system.z[j]*system.z[j]), sqrt(system.vx[0]*system.vx[0]+system.vy[0]*system.vy[0]+system.vz[0]*system.vz[0]), sqrt(system.vx[j]*system.vx[j]+system.vy[j]*system.vy[j]+system.vz[j]*system.vz[j]));
      fflush(distfile3);	//update file
    }
    if ((distfile4!=NULL)&&(system.phase[j]>=0)&&(system.a[j]>200.0)&&(system.a[j]<400.0)&&(system.prim.stage[j]==5)&&(system.sec.stage[j]<3)&&(system.sec.m[j]>4.0)&&(system.sec.m[j]<8.0)){	//write specified system to distfile4 output
      fprintf(distfile4,"%12.6e\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%s\t%4d\t%9.3f\t%10.6f  \t%10.6f\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.r[j], system.sec.r[j], system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], system.prim.stage[j], system.sec.stage[j],phasehistory, system.prim.SN.type, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi);
      fflush(distfile4);	//update file
    }
    if ((distfile5!=NULL)&&(system.phase[system.n-1]==94)&&(system.prim.t[system.n-1]<14.0e+9)){	//write specified system to distfile5 output
      if ((system.prim.stage[j]==5)&&(system.sec.stage[system.n-1]==5)&&(system.phase[j]==25)){
        fprintf(distfile5,"%11.6f\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%11.6f\t%11.6f\t%11.6f\t%9.3f\t%12.6f\t%10.6f\t%12.3f\t%10.3e\t%13.6e\t%s\t%10.4f\t%.5f\t%10.6f\t%10.6f\n", system.sec.SN.r_ini, system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.stage[j], system.sec.stage[j-1], system.sec.SN.Hecoremass, system.sec.SN.COcoremass, system.sec.SN.remnantmass, system.sec.SN.w, system.sec.SN.theta, system.sec.SN.phi, system.sec.SN.vsys, system.tgw, system.prim.t[system.n-1], phasehistory, system.a[0], system.e[0], system.prim.m[0], system.sec.m[0]);
        fflush(distfile5);	//update file
      }else if ((system.sec.stage[j]==5)&&(system.prim.stage[system.n-1]==5)&&(system.phase[j]==15)){
        fprintf(distfile5,"%11.6f\t%9.4f\t%.5f\t%10.6f\t%10.6f\t%2d /%2d\t%11.6f\t%11.6f\t%11.6f\t%9.3f\t%12.6f\t%10.6f\t%12.3f\t%10.3e\t%13.6e\t%s\t%10.4f\t%.5f\t%10.6f\t%10.6f\n", system.prim.SN.r_ini, system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.stage[j-1], system.sec.stage[j], system.prim.SN.Hecoremass, system.prim.SN.COcoremass, system.prim.SN.remnantmass, system.prim.SN.w, system.prim.SN.theta, system.prim.SN.phi, system.prim.SN.vsys, system.tgw, system.prim.t[system.n-1], phasehistory, system.a[0], system.e[0], system.prim.m[0], system.sec.m[0]);
        fflush(distfile5);	//update file
      }
    }
    if ((distfile6!=NULL)&&(system.phase[j]==0)&&(system.prim.stage[j]+system.sec.stage[j]==5)&&(system.prim.stage[j]<=3)&&(system.sec.stage[j]<=3)&&(system.M[j]>1.35)){	//write specified system to distfile6 output
      t = mergetime(system.prim.m[j],system.sec.m[j],system.a[j],system.e[j]);
      if (t<5.0e+8){
        if ((system.prim.stage[j]==2)&&(system.prim.m[j]>=0.4)&&(system.prim.m[j]<=1.1)){
          fprintf(distfile6,"%13.6e\t%10.4f\t%7.5f\t%10.6f\t%10.6f\t%2d /%2d\t%10.3e\t%10.4f\t%7.5f\t%10.6f\t%10.6f\t%s\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.stage[j], system.sec.stage[j], t, system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], phasehistory);
          fflush(distfile6);	//update file
        }else if ((system.sec.stage[j]==2)&&(system.sec.m[j]>=0.4)&&(system.sec.m[j]<=1.1)){
          fprintf(distfile6,"%13.6e\t%10.4f\t%7.5f\t%10.6f\t%10.6f\t%2d /%2d\t%10.3e\t%10.4f\t%7.5f\t%10.6f\t%10.6f\t%s\n", system.prim.t[j], system.a[j], system.e[j], system.prim.m[j], system.sec.m[j], system.prim.stage[j], system.sec.stage[j], t, system.a[0], system.e[0], system.prim.m[0], system.sec.m[0], phasehistory);
          fflush(distfile6);	//update file
        }
      }
    }
    if ((distfile7!=NULL)&&(system.phase[j-1]%10==3)&&(system.prim.stage[j]*system.sec.stage[j]==0)&&(system.prim.stage[j]<3)&&(system.sec.stage[j]<3)&&(system.prim.t[j]<1.3e+10)){	//write specified system to distfile7 output
      if ((system.phase[j-1]==13)&&(system.prim.cm[j-1]<0.5)){
        fprintf(distfile7,"%13.6e\t%10.4f\t%10.6f\t%10.6f\t%10.4f\t%2d /%2d\t%10.4f\t%10.6f\t%10.6f\t%d\t%s\n", system.prim.t[j], system.a[j], system.prim.cm[j-1], system.sec.m[j], system.sec.r[j], system.prim.stage[j], system.sec.stage[j], system.a[0], system.prim.m[0], system.sec.m[0], system.phase[j], phasehistory);
        fflush(distfile7);	//update file
      }else if ((system.phase[j-1]==23)&&(system.sec.cm[j-1]<0.5)){
        fprintf(distfile7,"%13.6e\t%10.4f\t%10.6f\t%10.6f\t%10.4f\t%2d /%2d\t%10.4f\t%10.6f\t%10.6f\t%d\t%s\n", system.prim.t[j], system.a[j], system.sec.cm[j-1], system.prim.m[j], system.prim.r[j], system.prim.stage[j], system.sec.stage[j], system.a[0], system.prim.m[0], system.sec.m[0], system.phase[j], phasehistory);
        fflush(distfile7);	//update file
      }
    }
  }
  if(Bep!=NULL) fprintf(Bep,"%10.5g\t%10.5g\t%10.5g",system.prim.omega[0]*system.prim.r[0]*Rsun/yr/1.0e3,system.prim.llum[0],system.prim.m[0]);
  for (t=1.0e5; t<=1.1e9; t=t*cbrt(10.0)){
    while ((system.prim.track.t[jp]<t)&&(jp+1<system.prim.track.n)){
      jp++;
    }
    rahier = (system.prim.track.t[jp]-t)/(system.prim.track.t[jp]-system.prim.track.t[jp-1]);
    if(Bep!=NULL) fprintf(Bep,"\t%10.5g\t%10.5g\t%10.5g",(system.prim.track.omega[jp]-rahier*(system.prim.track.omega[jp]-system.prim.track.omega[jp-1]))*(system.prim.track.r[jp]-rahier*(system.prim.track.r[jp]-system.prim.track.r[jp-1]))*Rsun/yr/1.0e3,(system.prim.track.llum[jp]-rahier*(system.prim.track.llum[jp]-system.prim.track.llum[jp-1])),(system.prim.track.m[jp]-rahier*(system.prim.track.m[jp]-system.prim.track.m[jp-1])));
  }
  if(Bep!=NULL) fprintf(Bep,"\n");
  if(Bes!=NULL) fprintf(Bes,"%10.5g\t%10.5g\t%10.5g",system.sec.omega[0]*system.sec.r[0]*Rsun/yr/1.0e3,system.sec.llum[0],system.sec.m[0]);
  for (t=1.0e5; t<=1.1e9; t=t*cbrt(10.0)){
    while ((system.sec.track.t[js]<t)&&(js+1<system.sec.track.n)){
      js++;
    }
    rahier = (system.sec.track.t[js]-t)/(system.sec.track.t[js]-system.sec.track.t[js-1]);
    if(Bes!=NULL) fprintf(Bes,"\t%10.5g\t%10.5g\t%10.5g",(system.sec.track.omega[js]-rahier*(system.sec.track.omega[js]-system.sec.track.omega[js-1]))*(system.sec.track.r[js]-rahier*(system.sec.track.r[js]-system.sec.track.r[js-1]))*Rsun/yr/1.0e3,(system.sec.track.llum[js]-rahier*(system.sec.track.llum[js]-system.sec.track.llum[js-1])),(system.sec.track.m[js]-rahier*(system.sec.track.m[js]-system.sec.track.m[js-1])));
  }
  if(Bes!=NULL) fprintf(Bes,"\n");
  
  for (j=1; j<system.n; j++){
    if(Be!=NULL){
      if (((system.phase[j-1]==12)&&(system.phase[j]==0))||((system.phase[j-1]==13)&&(system.phase[j]==0))){
        fprintf(Be,"%10.5g\t%10.5g\t%10.5g\n",2.0*M_PI*system.a[j]*Rsun*1.0e-3/system.P[j]/yr,system.prim.omega[j]*system.prim.r[j]*Rsun/yr/1.0e3,system.sec.omega[j]*system.sec.r[j]*Rsun/yr/1.0e3);
      }
    }
  }
  
  //~ alter = 1.8e7;//*ranv(6); 
  alter = 1.0e8*ranv(6);
  for (jp=1; jp<system.prim.track.n; jp++){
    if (system.prim.track.t[jp] > alter) break;
  }
  for (js=1; js<system.sec.track.n; js++){
    if (system.sec.track.t[js] > alter) break;
  }
  for (jj=0; jj<system.n; jj++){
    if (system.prim.t[jj] > alter) break;
  }
  //js,jp,jj are in the step AFTER alter -> jj-1 and jj etc are limits
  
  mem = 0;
  
  if(system.sec.stage[0]==6) mem = 6;
  
  if(system.phase[1]==12 && system.phase[2]==0 && jj-1==1) mem = 120; //Case A RLO
  if(system.phase[1]==12 && system.phase[2]==0 && jj-1<1) mem = 2;
  if(system.phase[1]==12 && system.phase[2]==0 && jj-1>1) mem = 12;
  if(system.phase[1]==12 && system.phase[2]==0) cas = 0;
  
  if(system.phase[2]==13 && system.phase[3]==0 && jj-1<1) mem = 2; //Case A CE
  if(system.phase[2]==13 && system.phase[3]==0 && jj-1>2) mem = 12;
  if(system.phase[2]==13 && system.phase[3]==0) cas = 1;
  
  if(system.phase[2]==12 && system.phase[3]==0 && jj-1<2) mem = 2; //Case B RLO
  if(system.phase[2]==12 && system.phase[3]==0 && jj-1>2) mem = 12;
  if(system.phase[2]==12 && system.phase[3]==0) cas = 2;
  
  if(system.phase[3]==13 && system.phase[4]==0 && jj-1<2) mem = 2; //Case B CE
  if(system.phase[3]==13 && system.phase[4]==0 && jj-1>3) mem = 12;
  if(system.phase[3]==13 && system.phase[4]==0) cas = 3;
  
  if(system.phase[3]==4 && jj-1 >= 3 && system.prim.stage[jj-1]==0) mem = 4; //Case A Merger
  
  st1 = system.prim.stage[jj-1];
  st2 = system.sec.stage[jj-1];
  
  if(mem!=0){
    M0 = system.M[jj-1];
    Om0 = 2.*M_PI/system.P[jj-1];
    ma1 = system.prim.track.m[jp] - (system.prim.track.t[jp]-alter)/(system.prim.track.t[jp]-system.prim.track.t[jp-1])*(system.prim.track.m[jp]-system.prim.track.m[jp-1]);
    ma2 = system.sec.track.m[js] - (system.sec.track.t[js]-alter)/(system.sec.track.t[js]-system.sec.track.t[js-1])*(system.sec.track.m[js]-system.sec.track.m[js-1]);
    ra1 = system.prim.track.r[jp] - (system.prim.track.t[jp]-alter)/(system.prim.track.t[jp]-system.prim.track.t[jp-1])*(system.prim.track.r[jp]-system.prim.track.r[jp-1]);
    ra2 = system.sec.track.r[js] - (system.sec.track.t[js]-alter)/(system.sec.track.t[js]-system.sec.track.t[js-1])*(system.sec.track.r[js]-system.sec.track.r[js-1]);
    lu1 = pow(10.0,system.prim.track.llum[jp] - (system.prim.track.t[jp]-alter)/(system.prim.track.t[jp]-system.prim.track.t[jp-1])*(system.prim.track.llum[jp]-system.prim.track.llum[jp-1]));
    lu2 = pow(10.0,system.sec.track.llum[js] - (system.sec.track.t[js]-alter)/(system.sec.track.t[js]-system.sec.track.t[js-1])*(system.sec.track.llum[js]-system.sec.track.llum[js-1]));
    
    Om = Om0*(ma1+ma2)/M0;
    om1 = Om - (Om-system.prim.track.omega[jp-1]) * pow((Om-system.prim.track.omega[jp])/(Om-system.prim.track.omega[jp-1]), (alter-system.prim.track.t[jp-1])/(system.prim.track.t[jp]-system.prim.track.t[jp-1]));
    om2 = Om - (Om-system.sec.track.omega[js-1]) * pow((Om-system.sec.track.omega[js])/(Om-system.sec.track.omega[js-1]), (alter-system.sec.track.t[js-1])/(system.sec.track.t[js]-system.sec.track.t[js-1]));
    ro1 = om1*ra1*rroteq(om1/sqrt(2.0*G*ma1/3.0/cubic(ra1))) * 1.0e-3*Rsun/yr;
    ro2 = om2*ra2*rroteq(om2/sqrt(2.0*G*ma2/3.0/cubic(ra2))) * 1.0e-3*Rsun/yr;
    //~ if(lu1>lu2){
      if((lu1>1000.0)&&(ro1>0.0)&&(st1==0)&&(deM!=NULL)) fprintf(deM,"%3i\t%2i\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\n",mem,cas,lu1,ro1,alter,system.prim.m[0],system.sec.m[0],system.a[0]);
    //~ }
    //~ if(lu2>lu1){
      if((lu2>1000.0)&&(ro2>0.0)&&(st2==0)&&(deM!=NULL)) fprintf(deM,"%3i\t%2i\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\n",mem,cas,lu2,ro2,alter,system.prim.m[0],system.sec.m[0],system.a[0]);
    //~ }
  }
    
      
      //~ mem1 = 0;
      //~ if(system.sec.stage[0]==6) mem1 = 6;
      //~ for(jj=0; system.prim.t[jj]<alter && jj<system.n; jj++){
        //~ if(system.phase[jj] == 12) mem1 = 12;
        //~ if(system.phase[jj] == 21) mem1 = 21;
        //~ if(system.phase[jj] == 4) mem1 = 4;
      //~ }
      //~ if(jj==system.n && mem1==0) mem1 = 6;
      //~ if(system.phase[jj-1] == 12) mem1 = 120;
      //~ if(system.phase[jj-1] == 21) mem1 = 210;
      //~ st1 = system.prim.stage[jj-1];
      //~ while(jj<system.n && mem1 == 0){
        //~ if(system.phase[jj] == 12 || system.phase[jj] == 21) mem1 = 2;
        //~ jj++;
      //~ }
      
      

  
  return;
}

void AfterMainLoop(){	//this function is called after the main loop of calculating the systems finishes
//  histout(initialmass, false, (char *)"initialmass");	//write the histogramm data to file
//  freehist(initialmass);	//free the histogramm memory
//  histout(initialperiod, false, (char *)"initialperiod");	//write the histogramm data to file
//  freehist(initialperiod);	//free the histogramm memory
  if (distfile1!=NULL) fclose(distfile1);
  if (distfile2!=NULL) fclose(distfile2);
  if (distfile3!=NULL) fclose(distfile3);
  if (distfile4!=NULL) fclose(distfile4);
  if (distfile5!=NULL) fclose(distfile5);
  if (distfile6!=NULL) fclose(distfile6);
  if (distfile7!=NULL) fclose(distfile7);
  if (snapshot0!=NULL) fclose(snapshot0);
  if (snapshot1!=NULL) fclose(snapshot1);
  if (snapshot2!=NULL) fclose(snapshot2);
  if (snapshot3!=NULL) fclose(snapshot3);
  if (snapshot4!=NULL) fclose(snapshot4);
  if (snapshot5!=NULL) fclose(snapshot5);
  if (snapshot6!=NULL) fclose(snapshot6);
  if (snapshot7!=NULL) fclose(snapshot7);
  if (snapshot8!=NULL) fclose(snapshot8);
  if (snapshot9!=NULL) fclose(snapshot9);
  if (snapshot10!=NULL) fclose(snapshot10);
  if (snapshot11!=NULL) fclose(snapshot11);
  if (snapshot12!=NULL) fclose(snapshot12);
  if (snapshot13!=NULL) fclose(snapshot13);
  if (Bep!=NULL) fclose(Bep);
  if (Bes!=NULL) fclose(Bes);
  if (Be!=NULL) fclose(Be);
  if (nathan!=NULL) fclose(nathan);
  if (deM!=NULL) fclose(deM);
  return;
}


