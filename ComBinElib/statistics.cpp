//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
//#include <fstream>	//in-/output to/from file
//#include <iomanip>	//manipulating the in-/output
#include "ComBinElib.h"	// local header file

//using namespace std;

/*counters*/
//warnings
long cnotrackupdate=0;	//count of skipped track updates
//events
long cRLO=0;	//count of Roche-Lobe overflows
long cRLOA=0;	//count of Roche-Lobe overflows of Case A
long cRLOB=0;	//count of Roche-Lobe overflows of Case B/C
long cRLOBB=0;	//count of Roche-Lobe overflows of Case BB
long cCE=0;	//count of common envelopes
long cCEA=0;	//count of common envelopes of Case A
long cCEB=0;	//count of common envelopes of Case B/C
long cCEBB=0;	//count of common envelopes of Case BB
long cMerger=0;	//count of merger events
long cMergerRLO=0;	//count of merger events in Roche-Lobe overflows
long cMergerCE=0;	//count of merger events in common envelopes
long cMergerSN=0;	//count of merger events in supernovae
long csupernova=0;	//count of supernova events
long cNS=0;	//count of fromed neutron stars
long cBH=0;	//count of fromed black holes
long cplanetarynebula=0;	//count of planetary nebula events
long cWD=0;	//count of fromed white dwarfs
//endstage
long cdestroyed=0;	//count of destroyed systems
long cmerged=0;	//count of merged systems
long cexploded=0;	//count of systems which are derupted by a supernova
long cunknown=0;	//count of systems which are derupted by an unknown process
long cfinal=0;	//count of systems which reach endstage
long cWDWD=0;	//count of white dwarf-white dwarf systems
long cWDNS=0;	//count of white dwarf-neutron star systems
long cWDBH=0;	//count of white dwarf-black hole systems
long cNSWD=0;	//count of neutron star-white dwarf systems
long cNSNS=0;	//count of neutron star-neutron star systems
long cNSBH=0;	//count of neutron star-black hole systems
long cBHWD=0;	//count of black hole-white dwarf systems
long cBHNS=0;	//count of black hole-neutron star systems
long cBHBH=0;	//count of black hole-black hole systems
long cgw=0;	//count of systems which merge by gravitational waves within tmax
long cgwWDWD=0;	//count of white dwarf-white dwarf systems which merge by gravitational waves within tmax
long cgwWDNS=0;	//count of white dwarf-neutron star systems which merge by gravitational waves within tmax
long cgwWDBH=0;	//count of white dwarf-black hole systems which merge by gravitational waves within tmax
long cgwNSWD=0;	//count of neutron star-white dwarf systems which merge by gravitational waves within tmax
long cgwNSNS=0;	//count of neutron star-neutron star systems which merge by gravitational waves within tmax
long cgwNSBH=0;	//count of neutron star-black hole systems which merge by gravitational waves within tmax
long cgwBHWD=0;	//count of black hole-white dwarf systems which merge by gravitational waves within tmax
long cgwBHNS=0;	//count of black hole-neutron star systems which merge by gravitational waves within tmax
long cgwBHBH=0;	//count of black hole-black hole systems which merge by gravitational waves within tmax
long cformation[10]={0,0,0,0,0,0,0,0,0,0};
long cgwformation[10]={0,0,0,0,0,0,0,0,0,0};

/*output data*/
bool data=true;	//write outputdata
bool dist=false;	//writeout the data of Mp, q, a, e, sigma, w, theta, phi, phase angle
bool grid=false;	//make a grid of Mp, Ms and a
bool plot=true;	//plot histogram(s)
ofstream dataout;	//enable writing to file the outputdata
ofstream inidist;	//enable writing to file the initial distributions
ofstream supernovadist;	//enable writing to file the supernova distributions
ofstream phasedist;	//enable writing to file the phase distributions
ofstream findist;	//enable writing to file the final distributions
ofstream gridout;	//enable writing to file the grid output data
FILE* distributionout=NULL;	//enable writing to a file
bool formationchange=false;

/*histograms*/ 
t_hist tmerge;

int getformation(int f){
  if (f==33) return 1;
  else if (f==34) return 2;
  else if (f==35) return 3;
  else if (f==43) return 4;
  else if (f==44) return 5;
  else if (f==45) return 6;
  else if (f==53) return 7;
  else if (f==54) return 8;
  else if (f==55) return 9;
  else return 0;
}

const char * getheader(int i){
  if (i==1) return "WD-WD";
  else if (i==2) return "WD-NS";
  else if (i==3) return "WD-BH";
  else if (i==4) return "NS-WD";
  else if (i==5) return "NS-NS";
  else if (i==6) return "NS-BH";
  else if (i==7) return "BH-WD";
  else if (i==8) return "BH-NS";
  else if (i==9) return "BH-BH";
  else return "unknown";
}

void initoutfiles(){	//initialize output files 
  if (dist){
    inidist.open("initial.dist");	// open file with name "initial.dist" to write
    inidist << "#Mp in Msun\tq\ta in Rsun\te" << endl;
    supernovadist.open("supernova.dist");	// open file with name "supernova.dist" to write
    supernovadist << "#sigma in km/s\tw in km/s\ttheta in degree\tphi in degree" << endl;
    phasedist.open("phase.dist");	// open file with name "phase.dist" to write
    phasedist << "#phase in degree" << endl;
    findist.open("final.dist");	// open file with name "final.dist" to write
    findist << "#Mp in Msun\tq\ta in Rsun\te\tphase\tprimstage\tsecstage\ttmerge in yr" << endl;
  }
  distributionout=fopen("distribution.csv","w");	// open file with name "distribution.csv" to write
  if (distributionout!=NULL){
    fprintf(distributionout,"#current\t         \t       \t          \t          \t          \t          \tinitial   \t       \t          \t          \tevolution     \tpost evolution\tprimary SN\t   \t            \t          \tsecondary SN\t \t            \t          \n");
    fprintf(distributionout,"#t in yr\ta in Rsun\t      e\tMp in Msun\tMs in Msun\tRp in Rsun\tRs in Rsun\t a in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tphases\tt_merge in yr \ttype\tw in km/s\ttheta in deg\tphi in deg\tvsys in km/s\ttype\tw in km/s\ttheta in deg\tphi in deg\tvsys in km/s\n");
//    fprintf(distributionout,"#t in yr\tP in days\t      e\tMp in Msun\tMs in Msun\tRp in Rsun\tRs in Rsun\t a in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tphases\tt_merge in yr \ttype\tw in km/s\ttheta in deg\tphi in deg\ttype\tw in km/s\ttheta in deg\tphi in deg\tv_sys in km/s\n");	//J1755-2550
/*    fprintf(distributionout,"#current\t         \t       \t          \t          \t          \t          \tinitial   \t       \t          \t          \tevolution     \tpost SN\n");
    fprintf(distributionout,"#t in yr\ta in Rsun\t      e\tMp in Msun\tMs in Msun\tRp in Rsun\tRs in Rsun\t a in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tphases\t a in Rsun\t        e\tMp in Msun\tMs in Msun\tacirc in Rsun\n");*/
/*    fprintf(distributionout,"#current\t         \t       \t          \t          \t          \t          \tinitial   \t       \t          \t          \tevolution     \tpost evolution\tprimary SN   \t             \t             \t    \t         \t            \t          \tsecondary SN \t             \t             \t    \t         \t            \t          \n");
    fprintf(distributionout,"#t in yr\ta in Rsun\t      e\tMp in Msun\tMs in Msun\tRp in Rsun\tRs in Rsun\t a in Rsun\t      e\tMp in Msun\tMs in Msun\tstages\tphases\tt_merge in yr \tM(He) in Msun\tM(CO) in Msun\tMremn in Msun\ttype\tw in km/s\ttheta in deg\tphi in deg\tM(He) in Msun\tM(CO) in Msun\tMremn in Msun\ttype\tw in km/s\ttheta in deg\tphi in deg\n");*/
    fflush(distributionout);	//update file
  }
  if (grid){
    gridout.open("grid.csv");	// open file with name "grid.csv" to write
    //write parameters to grid.csv
    gridout << "\tinitial\t \t \tfinal" << endl << "index\tMp(Msun)\tMs(msun)\ta(Rsun)\tphase\tprimstage\tsecstage\tMp(Msun)\tMs(Msun)\ta(Rsun)\ttmerge(yr)" << endl;
  }
}

void writedataout(int parallel, long seed, double Mp_max, double Mp_min, double Ms_max, double Ms_min, double a_max, double a_min, int galintegrator, int galpotential, long done){ 
  if (data){
    dataout.open("data.out");	// open file with name "data.out" to write
    //write input parameters to data.out
    dataout << "#parallel\tseed\tnumber\t(Mp_max,Mp_min)\t(Ms_max,Ms_min)\t(a_max,a_min)\talphaRLO\tbeta_const\tGamma\tdelta\talphaCE\tlambda_const\tqlimit\tMhece\tkickv\talphaTH\tgalintegrator\tgalpotential\toutput\tIMF\tqdist\tadist\tedist\ttdist\tZdist\trotdist\tRdist\tVdist\trhodist" << endl;
    dataout << parallel << "\t" << seed << "\t" << number << "\t(" << Mp_max << "," << Mp_min << ")\t(" << Ms_max << "," << Ms_min << ")\t(" << a_max << "," << a_min << ")\t" << alphaRLO << "\t" << beta_const << "\t" << Gamma << "\t" << delta << "\t" << alphaCE << "\t" << lambda_const << "\t" << qlimit << "\t" << Mhece << "\t" << kickv << "\t" << alphaTH << "\t" << galintegrator << "\t" << galpotential << "\t" << output << "\t" << IMF << "\t" << qdist << "\t" << adist << "\t" << edist << "\t" << tdist << "\t" << Zdist << "\t" << rotdist << "\t" << Rdist << "\t" << Vdist << "\t" << rhodist << endl;
    //write counters to data.out
    dataout << "#counts RLO: total\tcase A\tcase B/C\tcase BB" << endl;
    dataout << cRLO << "\t" << cRLOA << "\t" << cRLOB << "\t" << cRLOBB << endl;
    dataout << "#counts CE: total\tcase A\tcase B/C\tcase BB" << endl;
    dataout << cCE << "\t" << cCEA << "\t" << cCEB << "\t" << cCEBB << endl;
    dataout << "#counts merger: total\tin RLO\tin CE\tin SN" << endl;
    dataout << cMerger << "\t" << cMergerRLO << "\t" << cMergerCE << "\t" << cMergerSN << endl;
    dataout << "#counts SN/PN: SN\tNS\tBH\tPN\tWD" << endl;
    dataout << csupernova << "\t" << cNS << "\t" << cBH << "\t" << cplanetarynebula << "\t" << cWD << endl;
    dataout << "#number of systems\tsurvive(SS)\tdestroyed(DS)" << endl;
    dataout << done << "\t" << cfinal << "\t" << cdestroyed << endl;
    dataout << "#number of SS\tWD-WD\tWD-NS\tWD-BH\tNS-WD\tNS-NS\tNS-BH\tBH-WD\tBH-NS\tBH-BH\tunknown" << endl;
//    dataout << cfinal << "\t" << cWDWD << "\t" << cWDNS << "\t" << cWDBH << "\t" << cNSWD << "\t" << cNSNS << "\t" << cNSBH << "\t" << cBHWD << "\t" << cBHNS << "\t" << cBHBH << endl;
    dataout << cfinal << "\t" << cformation[getformation(33)] << "\t" << cformation[getformation(34)] << "\t" << cformation[getformation(35)] << "\t" << cformation[getformation(43)] << "\t" << cformation[getformation(44)] << "\t" << cformation[getformation(45)] << "\t" << cformation[getformation(53)] << "\t" << cformation[getformation(54)] << "\t" << cformation[getformation(55)] << "\t" << cformation[0] << endl;
    dataout << "#SS t_gw<" << tmax << "yr\tWD-WD\tWD-NS\tWD-BH\tNS-WD\tNS-NS\tNS-BH\tBH-WD\tBH-NS\tBH-BH\tunknown" << endl;
//    dataout << cgw << "\t" << cgwWDWD << "\t" << cgwWDNS << "\t" << cgwWDBH << "\t" << cgwNSWD << "\t" << cgwNSNS << "\t" << cgwNSBH << "\t" << cgwBHWD << "\t" << cgwBHNS << "\t" << cgwBHBH << endl;
    dataout << cgw << "\t" << cgwformation[getformation(33)] << "\t" << cgwformation[getformation(34)] << "\t" << cgwformation[getformation(35)] << "\t" << cgwformation[getformation(43)] << "\t" << cgwformation[getformation(44)] << "\t" << cgwformation[getformation(45)] << "\t" << cgwformation[getformation(53)] << "\t" << cgwformation[getformation(54)] << "\t" << cgwformation[getformation(55)] << "\t" << cgwformation[0] << endl;
    dataout << "#number of DS\tmerged sys.\tdisrupted sys.\tunknown" << endl;
    dataout << cdestroyed << "\t" << cmerged << "\t" << cexploded << "\t" << cunknown << endl;
    dataout.close();	//close the data.out-file
  }
}

void closeoutfiles(){	//close output files 
  if (dist){	//close distribution outputs
    inidist.close();
    supernovadist.close();
    phasedist.close();
    findist.close();
  }
  if (distributionout!=NULL) fclose(distributionout);
  if (grid){	//close grid output
    gridout.close();
  }
}

void inithist(t_hist& hist, int nbin, double min, double max, bool logscale, int subs){	//initialize a histogram
  int i,j;	//loop variables

  //check for a positive natural number in nbin and subs
  if (nbin<1) nbin = 1;
  if (subs<1) subs = 1;

  hist.n = nbin;
  hist.subs = subs;
  //reserve memory
  hist.c = (long *)malloc(nbin*sizeof(long));
  hist.x = (double *)malloc((nbin+1)*sizeof(double));
  hist.subc = (long **)malloc(subs*sizeof(long*));
  //check if memory allocation fails
  if (hist.c==NULL) cerr << "#Error: memory allocation failed: hist.c" << endl;
  if (hist.x==NULL) cerr << "#Error: memory allocation failed: hist.x" << endl;
  if (hist.subc==NULL) cerr << "#Error: memory allocation failed: hist.subc" << endl;
  for (j=0;j<subs;j++){
    //reserve memory
    hist.subc[j] = (long *)malloc(nbin*sizeof(long));
    //check if memory allocation fails
    if (hist.subc[j]==NULL) cerr << "#Error: memory allocation failed: hist.subc[" << j << "]" << endl;
  }
  //set initial values
  hist.min = -1.0e-99;
  hist.max = hist.min;
  hist.ctot = 0;
  hist.logscale = logscale;
  for (i=0;i<nbin;i++){
    hist.c[i] = 0;
    if (logscale) hist.x[i] = pow(10.0,log10(min)+(double)i*(log10(max)-log10(min))/(double)(nbin));
    else hist.x[i] = min+(double)i*(max-min)/(double)(nbin);
    for (j=0;j<subs;j++) hist.subc[j][i] = 0;
  }
  hist.x[nbin] = max;
}

void addtohist(t_hist& hist, double value, int sub){	//adds the value to the histogram
  int i;	//loop variable

#pragma omp critical
{
  if (hist.ctot==0){	//first 
    hist.min = value;
    hist.max = hist.min;
  }else{	//update minimum and maximum values
    if (hist.min>value) hist.min = value;
    else if (hist.max<value) hist.max = value;
  }
  //find bin (anything above the range contains to last bin and all below the range contains to first bin)
  for (i=hist.n-1;i>0;i--){
    if (hist.x[i]<value) break;
  }
  hist.c[i]++;	//increase bin entry
  hist.ctot++;
  if ((0<=sub)&&(sub<hist.subs)) hist.subc[sub][i]++;
}	//end critical
}

void freehist(t_hist& hist){	//free the histogram
  int j;	//loop variable

  hist.n = 0;
  for (j=0;j<hist.subs;j++) free(hist.subc[j]);
  free(hist.c);
  free(hist.x);
  hist.subs = 0;
  free(hist.subc);
}

void histout(t_hist& hist, bool perbinsize, char* name){	//writes the histogram data to file
  ofstream histfile;	//enable writing to file the output file of the histogram
  int i,j;	//loop variables
  double xmin, xmax;	//minimum and maximum value of the x-axis
  double inversebinsize;	//1.0/bin size
  char filename[1000];

  sprintf(filename,"%s.hist",name);
  histfile.open(filename);	// open file to write the histogram
  histfile << "#Bins:" << hist.n << "\tSubs:" << hist.subs << endl;	// write number of bins
  if (&hist==&tmerge) histfile << "#tgw[yr] ";	//header of table
  else histfile << "#x\t";	//header of table
  xmin = fmin(hist.min,hist.x[0]);	//determine minimal x-value
  xmax = fmax(hist.max,hist.x[hist.n]);	//determine maximal x-value
  if (perbinsize){
    hist.x[0] = xmin;	//adjust histogram border
    hist.x[hist.n] = xmax;	//adjust histogram border
    histfile << "count/bin size";	//header of table
    for (j=0;j<hist.subs;j++){
      if (&hist==&tmerge) histfile << "\t" << getheader(j) << "/bin size";	//header of table
      else histfile << "\tsubcount" << j << "/bin size";	//header of table
    }
  }else{
    histfile << "counts";	//header of table
    for (j=0;j<hist.subs;j++){
      if (&hist==&tmerge) histfile << "\t" << getheader(j);	//header of table
      else histfile << "\tsubcount" << j;	//header of table
    }
  }
  /*write histogram data*/
  for (i=0;i<hist.n;i++){
    histfile << endl << hist.x[i] << "\t";	//write x-value
    if (perbinsize){
      inversebinsize = 1.0/(hist.x[i+1]-hist.x[i]);	//determine the inverse of the bin size
      histfile << hist.c[i]*inversebinsize;	//write y-value
      for (j=0;j<hist.subs;j++) histfile << "\t" << hist.subc[j][i]*inversebinsize;	//write suby-values
    }else{
      histfile << hist.c[i];	//write y-value
      for (j=0;j<hist.subs;j++) histfile << "\t" << hist.subc[j][i];	//write suby-values
    }
  }
  histfile << endl << hist.x[hist.n];
  histfile.close();	//close the output file of the histogram
}

void adddata(t_system& system, long i){
  int n=system.n-1;
  int sub=0;

  //counters
  if (system.phase[n]==-1){	//system is destroyed
    cdestroyed++;
    if (system.phase[n-1]==4) cmerged++;	//system has merged ((system.phase[n-1]==4)||((system.phase[n-1]==-1)&&(system.phase[n-2]==4)))
    else if ((system.phase[n-1]%10==5)||((system.phase[n-1]==-1)&&(system.phase[n-2]%10==5))) cexploded++;	//system is desrupted by a supernova
    else {
      cunknown++;
      screen = true;
    }
  }else if (system.phase[n]==99){	//system get to endstage
    cfinal++;
#pragma omp critical
{
    cformation[getformation(system.formation)]++;
}	//end critical
    if (getformation(system.formation)==0) screen = true;
    if (system.prim.stage[n]==3){	//primary is a white dwarf
      if (system.sec.stage[n]==3){	//secondary is a white dwarf
        cWDWD++;
        sub = 1;
        if ((system.formation!=33)&& formationchange){
          cout << "33->" << system.formation << endl;
          screen = true;
        }
      }else if (system.sec.stage[n]==4){	//secondary is a neutron star
        cWDNS++;
        sub = 2;
        if ((system.formation!=34)&& formationchange){
          cout << "34->" << system.formation << endl;
          screen = true;
        }
      }else if (system.sec.stage[n]==5){	//secondary is a black hole
        cWDBH++;
        sub = 3;
        if ((system.formation!=35)&& formationchange){
          cout << "35->" << system.formation << endl;
          screen = true;
        }
      }
    }else if (system.prim.stage[n]==4){	//primary is a neutron star
      if (system.sec.stage[n]==3){	//secondary is a white dwarf
        cNSWD++;
        sub = 4;
        if ((system.formation!=43)&& formationchange){
          cout << "43->" << system.formation << endl;
          screen = true;
        }
      }else if (system.sec.stage[n]==4){	//secondary is a neutron star
        cNSNS++;
        sub = 5;
        if ((system.formation!=44)&& formationchange){
          cout << "44->" << system.formation << endl;
          screen = true;
        }
      }else if (system.sec.stage[n]==5){	//secondary is a black hole
        cNSBH++;
        sub = 6;
        if ((system.formation!=45)&& formationchange){
          cout << "45->" << system.formation << endl;
          screen = true;
        }
      }
    }else if (system.prim.stage[n]==5){	//primary is a black hole
      if (system.sec.stage[n]==3){	//secondary is a white dwarf
        cBHWD++;
        sub = 7;
        if ((system.formation!=53)&& formationchange){
          cout << "53->" << system.formation << endl;
          screen = true;
        }
      }else if (system.sec.stage[n]==4){	//secondary is a neutron star
        cBHNS++;
        sub = 8;
        if ((system.formation!=54)&& formationchange){
          cout << "54->" << system.formation << endl;
          screen = true;
        }
      }else if (system.sec.stage[n]==5){	//secondary is a black hole
        cBHBH++;
        sub = 9;
        if ((system.formation!=55)&& formationchange){
          cout << "55->" << system.formation << endl;
          screen = true;
        }
      }
    }
/*    if (system.prim.stage[n]==3){	//primary is a white dwarf
      if (system.sec.stage[n]==3) cWDWD++;	//secondary is a white dwarf
      else if (system.sec.stage[n]==4) cWDNS++;	//secondary is a neutron star
      else if (system.sec.stage[n]==5) cWDBH++;	//secondary is a black hole
    }else if (system.prim.stage[n]==4){	//primary is a neutron star
      if (system.sec.stage[n]==3) cNSWD++;	//secondary is a white dwarf
      else if (system.sec.stage[n]==4) cNSNS++;	//secondary is a neutron star
      else if (system.sec.stage[n]==5) cNSBH++;	//secondary is a black hole
    }else if (system.prim.stage[n]==5){	//primary is a black hole
      if (system.sec.stage[n]==3) cBHWD++;	//secondary is a white dwarf
      else if (system.sec.stage[n]==4) cBHNS++;	//secondary is a neutron star
      else if (system.sec.stage[n]==5) cBHBH++;	//secondary is a black hole
    }*/
    system.tgw = mergetime(system.prim.m[n],system.sec.m[n],system.a[n],system.e[n]);
//    if (tgw<1.0) screen = true;
    if (screen) cout << "mergertime = " << system.tgw << "yr" << endl;
    if (system.tgw<tmax){	//system merges due to gravitational wave radiation before tmax
      cgw++;
#pragma omp critical
{
      cgwformation[getformation(system.formation)]++;
}	//end critical
      if (system.prim.stage[n]==3){	//primary is a white dwarf
        if (system.sec.stage[n]==3){	//secondary is a white dwarf
          cgwWDWD++;
//          sub = 1;
        }else if (system.sec.stage[n]==4){	//secondary is a neutron star
          cgwWDNS++;
//          sub = 2;
        }else if (system.sec.stage[n]==5){	//secondary is a black hole
          cgwWDBH++;
//          sub = 3;
        }
      }else if (system.prim.stage[n]==4){	//primary is a neutron star
        if (system.sec.stage[n]==3){	//secondary is a white dwarf
          cgwNSWD++;
//          sub = 4;
        }else if (system.sec.stage[n]==4){	//secondary is a neutron star
          cgwNSNS++;
//          sub = 5;
        }else if (system.sec.stage[n]==5){	//secondary is a black hole
          cgwNSBH++;
//          sub = 6;
        }
      }else if (system.prim.stage[n]==5){	//primary is a black hole
        if (system.sec.stage[n]==3){	//secondary is a white dwarf
          cgwBHWD++;
//          sub = 7;
        }else if (system.sec.stage[n]==4){	//secondary is a neutron star
          cgwBHNS++;
//          sub = 8;
        }else if (system.sec.stage[n]==5){	//secondary is a black hole
          cgwBHBH++;
//          sub = 9;
        }
      }
    }
//    if (plot) addtohist(tmerge,system.tgw,sub);
    if (plot) addtohist(tmerge,system.tgw,getformation(system.formation));
    if (debug) cerr << endl << "sub=" << sub << endl;
  }

  //distribution output
  if (dist){
    //writes initial/final primary mass, mass ratio, semi-major axis, eccentricity
    inidist << system.prim.m[0] << "\t" << system.qs[0] << "\t" << system.a[0] << "\t" << system.e[0] << endl;
    if (system.prim.m[n]>system.sec.m[n])
      findist << system.prim.m[n] << "\t" << system.qs[n];
    else
      findist << system.sec.m[n] << "\t" << system.qp[n];
    findist << "\t" << system.a[n] << "\t" << system.e[n] << "\t" << system.phase[n] << "\t" << system.prim.stage[n] << "\t" << system.prim.stage[n] << "\t" << system.tgw << endl;
  }
  //grid output
  if (grid){
    if (system.phase[n]==99) gridout << i << "\t" << system.prim.m[0] << "\t" << system.sec.m[0] << "\t" << system.a[0] << "\t" << system.phase[n] << "\t" << system.prim.stage[n] << "\t" << system.sec.stage[n] << "\t" << system.prim.m[n] << "\t" << system.sec.m[n] << "\t" << system.a[n] << "\t" << system.tgw << endl;
    else if (system.a[n]>0) gridout << i << "\t" << system.prim.m[0] << "\t" << system.sec.m[0] << "\t" << system.a[0] << "\t" << system.phase[n] << "\t" << system.prim.stage[n] << "\t" << system.sec.stage[n] << "\t" << system.prim.m[n] << "\t" << system.sec.m[n] << "\t" << system.a[n] << endl;
    else gridout << i << "\t" << system.prim.m[0] << "\t" << system.sec.m[0] << "\t" << system.a[0] << "\t" << system.phase[n] << "\t" << system.prim.stage[n] << "\t" << system.sec.stage[n] << "\t" << system.prim.m[n] << "\t" << system.sec.m[n] << endl;
  }
/*
//  Testing GWR timescale:
//  Mp[10]=1.441;
//  Ms[10]=1.387;
//  a[10]=2.801;
//  ecc=0.61713;
    t_merge=mergetime(Mp[10],Ms[10],a[10],ecc);
    if (screen) cout << endl << "TIME MERGING (yr): " << t_merge << endl;
    maxdist=((times[10]-times[5])*vCM1 + t_merge*vCM2)*3.1557E+7/3.08E+13;
    if (screen) cout << " max. dist. SN1+SN2 (pc): " << maxdist << endl << endl;
//  break; // Finish with one system
//  if (data && inform[5]==1 && inform[10]==1) dataout << endl << Mp[0] << " " << Ms[0] << " " << a[0] << " "
//    << Mp[10] << " " << Ms[10] << " " << a[10] << endl 
//    << inform[1] << " " << inform[2] << " " << inform[3] << " " 
//    << inform[4] << " " << inform[5] << " " << inform[6] << " " 
//    << inform[7] << " " << inform[8] << " " << inform[9] << " " 
//    << inform[10] << " " << t_merge << " " << vCM2 << " " << ecc << " " << maxdist;
//if (data && (inform[5]==1 || inform[5]==2) && (inform[10]==1 || inform[10]==2) && (inform[6]==1))
  if (data && (inform[5]==1 || inform[5]==2) && (inform[10]==1 || inform[10]==2)) 	//NS/BH+NS/BH
//if (data && (inform[5]==2) && (inform[10]==0))   // BH+WD
//if (data && (inform[5]==2) && (inform[10]==1))   // BH+NS
//if (data && (inform[5]==1) && (inform[10]==1))   // NS+NS
//if (data && (inform[5]==1) && (inform[10]==2))   // NS+BH
//if (data && (inform[5]==2) && (inform[10]==2))   // BH+BH
//if (data && (inform[5]==1 && inform[10]==1)) 
//if (data && (inform[5]==1 && inform[10]==1)) dataout << endl << a[10] << " " << ecc; 
//{  
// dataout << endl << Mp[0] << " " << Ms[0] << " " << a[0] << " " << a[10] << " " << ecc << " " << inform[6];
// cout << endl << "Ms Ms_AO MHe Mhece: " << Ms[0] << " " << Ms[5] << " " << Ms[7] << " " << Mhece;
// if (Ms[7]<4.0) cout << endl << "Ms Ms_AO MHe Mhece: " << Ms[0] << " " << Ms[5] << " " << Ms[7] << " " << Mhece;
   dataout << endl << Mp[0] << " " << Ms[0] << " " << a[0] << " " << a[10] << " " << ecc << " " 
           << t_merge << " " << maxdist << " " << vCM1 << " " << vCM2 << " "
           << inform[5] << " " << inform[10] << " " << Mp[10] << " " << Ms[10] << " " << inform[1];
//
//}
//  if (data && inform[5]==1 && inform[10]==1 && a[10] < 1.0 && ecc > 0.55) 
//    {
//       cout << endl << "Mp Ms a0 a ecc: " << Mp[0] << " " << Ms[0] << " " << a[0] << " " << a[10] << " " << ecc;
//    }
*/
  return;
}
