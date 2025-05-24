//#include <iostream>	//in-/output to/from screen
//#include <math.h>	//provides the mathematical functions
#include "ComBinElib.h" // local header file

//using namespace std;

void evolve(t_system& system){	//evolve system
  int n=system.n-1;	//last entry index
  int jp=1,js=1;	//entry index (p,s denote primary/secondary)
  int nextphase=0;	//next phase
  double t=0.0;	//time
  double r=0.0;	//radius
  double ratiop,ratios;	//ratio between the two nearest times/radii (p,s denote primary/secondary)
  double EbindperG=0.0;	//envelope binding energy for interpolation

  if (((system.prim.t[n]>system.prim.track.t[system.prim.track.n-1])||(system.sec.t[n]>system.sec.track.t[system.sec.track.n-1]))&&(system.phase[n]>=0)&&(system.phase[n]<99)){	//check if both stars are not evolved above their maximal age
    cerr << "#Warning: prim.t[" << n << "]=" << system.prim.t[n] << "yr prim.track.t[" << system.prim.track.n-1 << "]=" << system.prim.track.t[system.prim.track.n-1] << "yr sec.t[" << n << "]=" << system.sec.t[n] << "yr sec.track.t[" << system.sec.track.n-1 << "]=" << system.sec.track.t[system.sec.track.n-1] << "yr" << endl;	//write warning
    screen = true;	//enable screen output
  }
  if (debug) cerr << "prim.inimass=" << system.prim.inimass << "Msun sec.inimass=" << system.sec.inimass << "Msun ";	//additional output for debugging

  if (system.phase[n]==0){	//evolution in wind mass loss phase
    if (screen && (output=='M')) cout << "wind: ";
    t = system.prim.t[n];	//get current simulation time
    nextphase = wind(system,t,r);	//check till when the stars evolve without interaction, kind of interaction will be returned
    system.phase[n] = nextphase;	//save next phase
    if ((nextphase==12)||(nextphase==21)||(nextphase%10==3)||(nextphase==4)){	//if stars interact/are very close together
      circularize(system);	//apply circularization
    }
    if (nextphase!=-1){	//determine current system/star values if system did not crash
      //determine/interpolate values for the primary
      if (((nextphase==12)||(nextphase==4))&&(t>0.0)&&(r>0.0)){	//use radius for interpolation to avoid equal times
        for (jp=system.prim.last+1;jp<system.prim.track.n-1;jp++){	//find time index of r
          if ((fmax(system.prim.track.r[jp-1],system.prim.track.r[jp])>r)&&(fmin(system.prim.track.r[jp-1],system.prim.track.r[jp])<=r)&&(system.prim.track.t[jp]>=t)) break;	//&&(jp>system.prim.last)
//          if ((system.prim.track.r[jp]>r)&&(system.prim.track.r[jp-1]<=r)) break;
          if ((fmax(system.prim.track.r[jp-1],system.prim.track.r[jp])>r)&&(fmin(system.prim.track.r[jp-1],system.prim.track.r[jp])<=r)&&(debug)) cerr << " t=" << t << "yr system.prim.track.t[" << jp << "]=" << system.prim.track.t[jp] << "yr system.prim.track.t[" << jp-1 << "]=" << system.prim.track.t[jp-1] << "yr diff=" << system.prim.track.t[jp]-system.prim.track.t[jp-1] << "yr" << endl;	//additional output for debugging
        }
        ratiop = ratio(r,system.prim.track.r[jp],system.prim.track.r[jp-1]);	//calculate ratio in radius for the interpolation
        if (debug) cerr << " r=" << r << "Rsun system.prim.track.r[" << jp << "]=" << system.prim.track.r[jp] << "Rsun system.prim.track.r[" << jp-1 << "]=" << system.prim.track.r[jp-1] << "Rsun diff=" << system.prim.track.r[jp]-system.prim.track.r[jp-1] << "Rsun system.prim.track.n-1=" << system.prim.track.n-1 << " ratio_p=" << ratiop << endl;	//additional output for debugging
        if ((ratiop<0.0)||(ratiop>1.0)){	//out of interpolation range
          cerr << "#Error: radius out of track: r=" << r << "Rsun system.prim.track.r[" << jp << "]=" << system.prim.track.r[jp] << "Rsun system.prim.track.r[" << jp-1 << "]=" << system.prim.track.r[jp-1] << "Rsun diff=" << system.prim.track.r[jp]-system.prim.track.r[jp-1] << "Rsun system.prim.track.n-1=" << system.prim.track.n-1 << " ratio_p=" << ratiop << " system.prim.last=" << system.prim.last << endl;	//write error
          screen = true;	//enable screen output
          if (ratiop<0.0) ratiop = 0.0; else ratiop = 1.0;	//go back to nearest border of the interpolation range
        }
      }else{	//use time for interpolation
        for (jp=system.prim.last;jp<system.prim.track.n-1;jp++){	//find time index of t
          if (system.prim.track.t[jp]>=t) break;
        }
        ratiop = ratio(t,system.prim.track.t[jp],system.prim.track.t[jp-1]);	//calculate ratio in time for the interpolation
        if (debug) cerr << " t=" << t << "yr system.prim.track.t[" << jp << "]=" << system.prim.track.t[jp] << "yr system.prim.track.t[" << jp-1 << "]=" << system.prim.track.t[jp-1] << "yr diff=" << system.prim.track.t[jp]-system.prim.track.t[jp-1] << "yr system.prim.track.n-1=" << system.prim.track.n-1 << " ratio_p=" << ratiop << endl;	//additional output for debugging
        if ((ratiop<0.0)||(ratiop>1.0)){	//out of interpolation range
          cerr << endl << "#Error: time out of track: t=" << t << "yr system.prim.track.t[" << jp << "]=" << system.prim.track.t[jp] << "yr system.prim.track.t[" << jp-1 << "]=" << system.prim.track.t[jp-1] << "yr diff=" << system.prim.track.t[jp]-system.prim.track.t[jp-1] << "yr system.prim.track.n-1=" << system.prim.track.n-1 << " ratio_p=" << ratiop << " system.prim.last=" << system.prim.last << endl;	//write error
          screen = true;	//enable screen output
          if (ratiop<0.0) ratiop = 0.0; else ratiop = 1.0;	//go back to nearest border of the interpolation range
        }
      }
      system.prim.t[n] = t;	//set new time
      system.prim.m[n] = system.prim.track.m[jp]-ratiop*(system.prim.track.m[jp]-system.prim.track.m[jp-1]);	//interpolate mass
      system.prim.r[n] = system.prim.track.r[jp]-ratiop*(system.prim.track.r[jp]-system.prim.track.r[jp-1]);	//interpolate radius
      system.prim.cm[n] = system.prim.track.cm[jp]-ratiop*(system.prim.track.cm[jp]-system.prim.track.cm[jp-1]);	//interpolate core mass
      system.prim.llum[n] = system.prim.track.llum[jp]-ratiop*(system.prim.track.llum[jp]-system.prim.track.llum[jp-1]);	//interpolate luminosity
      system.prim.lteff[n] = system.prim.track.lteff[jp]-ratiop*(system.prim.track.lteff[jp]-system.prim.track.lteff[jp-1]);	//interpolate effective temperature
//      system.prim.lambda[n] = system.prim.track.lambda[jp]-ratiop*(system.prim.track.lambda[jp]-system.prim.track.lambda[jp-1]);	//interpolate lambda
      EbindperG = system.prim.track.m[jp]*(system.prim.track.m[jp]-system.prim.track.cm[jp])/(system.prim.track.lambda[jp]*system.prim.track.r[jp])
                  -ratiop*(system.prim.track.m[jp]*(system.prim.track.m[jp]-system.prim.track.cm[jp])/(system.prim.track.lambda[jp]*system.prim.track.r[jp])
                          -system.prim.track.m[jp-1]*(system.prim.track.m[jp-1]-system.prim.track.cm[jp-1])/(system.prim.track.lambda[jp-1]*system.prim.track.r[jp-1]));	//interpolate current binding energy of the envelope
      if ((system.prim.m[n]-system.prim.cm[n]==0.0)&&(EbindperG==0.0)) system.prim.lambda[n] = -1.0e+99;	//no-interpolation
      else system.prim.lambda[n] = system.prim.m[n]*(system.prim.m[n]-system.prim.cm[n])/(system.prim.r[n]*EbindperG);	//get lambda from interpolated E_bind	// E_bind=G*M*M_env / lambda*R
      if(isnan(system.prim.lambda[n])){ screen = true; cerr << endl << "#system.prim.lambda[n]=" << system.prim.lambda[n] << " system.prim.m[n]=" << system.prim.m[n] << "Msun system.prim.cm[n]=" << system.prim.cm[n] << "Msun jp=" << jp << endl << "EbindperG=" << EbindperG << "Msun^2/Rsun=" << system.prim.track.m[jp] << "*(" << system.prim.track.m[jp] << "-" << system.prim.track.cm[jp] << ")/(" << system.prim.track.lambda[jp] << "*" << system.prim.track.r[jp] << ")-" << ratiop << "*(" << system.prim.track.m[jp] << "*(" << system.prim.track.m[jp] << "-" << system.prim.track.cm[jp] << ")/(" << system.prim.track.lambda[jp] << "*" << system.prim.track.r[jp] << ")-" << system.prim.track.m[jp-1] << "*(" << system.prim.track.m[jp-1] << "-" << system.prim.track.cm[jp-1] << ")/(" << system.prim.track.lambda[jp-1] << "*" << system.prim.track.r[jp-1] << "))Msun^2/Rsun";}
      system.prim.last = jp-1;	//remember position of the interpolation in the track array
      if (system.prim.r[n]<0) cerr << "#Error: system.prim.r[" << n << "]=" << system.prim.r[n] << "Rsun system.prim.track.r[" << jp << "]=" << system.prim.track.r[jp] << "Rsun system.prim.track.r[" << jp-1 << "]=" << system.prim.track.r[jp-1] << "Rsun ratio_p=" << ratiop << endl;	//check for negative radius	//write error
      //determine/interpolate values for the secondary
      if ((nextphase==21)&&(t>0.0)&&(r>0.0)){	//use radius for interpolation to avoid equal times
        for (js=system.sec.last+1;js<system.sec.track.n-1;js++){	//find time index of r
          if ((fmax(system.sec.track.r[js-1],system.sec.track.r[js])>r)&&(fmin(system.sec.track.r[js-1],system.sec.track.r[js])<=r)&&(system.sec.track.t[js]>=t)) break;	//&&(js>system.sec.last)
//          if (system.sec.track.r[js]>r) break;
        }
        ratios = ratio(r,system.sec.track.r[js],system.sec.track.r[js-1]);	//calculate ratio in radius for the interpolation
        if (debug) cerr << " r=" << r << "Rsun system.sec.track.r[" << js << "]=" << system.sec.track.r[js] << "Rsun system.sec.track.r[" << js-1 << "]=" << system.sec.track.r[js-1] << "Rsun diff=" << system.sec.track.r[js]-system.sec.track.r[js-1] << "Rsun system.sec.track.n-1=" << system.sec.track.n-1 << " ratio_s=" << ratios << endl;	//additional output for debugging
        if ((ratios<0.0)||(ratios>1.0)){	//out of interpolation range
          cerr << "#Error: radius out of track: r=" << r << "Rsun system.sec.track.r[" << js << "]=" << system.sec.track.r[js] << "Rsun system.sec.track.r[" << js-1 << "]=" << system.sec.track.r[js-1] << "Rsun diff=" << system.sec.track.r[js]-system.sec.track.r[js-1] << "Rsun system.sec.track.n-1=" << system.sec.track.n-1 << " ratio_s=" << ratios << " system.sec.last=" << system.sec.last << endl;	//write error
          screen = true;	//enable screen output
          if (ratios<0.0) ratios = 0.0; else ratios = 1.0;	//go back to nearest border of the interpolation range
        }
      }else{	//use time for interpolation
        for (js=system.sec.last;js<system.sec.track.n-1;js++){	//find time index of t
          if (system.sec.track.t[js]>=t) break;
        }
        ratios = ratio(t,system.sec.track.t[js],system.sec.track.t[js-1]);	//calculate ratio in time for the interpolation
        if (debug) cerr << " t=" << t << "yr system.sec.track.t[" << js << "]=" << system.sec.track.t[js] << "yr system.sec.track.t[" << js-1 << "]=" << system.sec.track.t[js-1] << "yr diff=" << system.sec.track.t[js]-system.sec.track.t[js-1] << "yr system.sec.track.n-1=" << system.sec.track.n-1 << " ratio_s=" << ratios << endl;	//additional output for debugging
        if ((ratios<0.0)||(ratios>1.0)){	//out of interpolation range
          cerr << endl << "#Error: time out of track: t=" << t << "yr system.sec.track.t[" << js << "]=" << system.sec.track.t[js] << "yr system.sec.track.t[" << js-1 << "]=" << system.sec.track.t[js-1] << "yr diff=" << system.sec.track.t[js]-system.sec.track.t[js-1] << "yr system.sec.track.n-1=" << system.sec.track.n-1 << " ratio_s=" << ratios << " system.sec.last=" << system.sec.last << endl;	//write error
          screen = true;	//enable screen output
          if (ratios<0.0) ratios = 0.0; else ratios = 1.0;	//go back to nearest border of the interpolation range
        }
      }
      system.sec.t[n] = t;	//set new time
      system.sec.m[n] = system.sec.track.m[js]-ratios*(system.sec.track.m[js]-system.sec.track.m[js-1]);	//interpolate mass
      system.sec.r[n] = system.sec.track.r[js]-ratios*(system.sec.track.r[js]-system.sec.track.r[js-1]);	//interpolate radius
      system.sec.cm[n] = system.sec.track.cm[js]-ratios*(system.sec.track.cm[js]-system.sec.track.cm[js-1]);	//interpolate core mass
      system.sec.llum[n] = system.sec.track.llum[js]-ratios*(system.sec.track.llum[js]-system.sec.track.llum[js-1]);	//interpolate luminosity
      system.sec.lteff[n] = system.sec.track.lteff[js]-ratios*(system.sec.track.lteff[js]-system.sec.track.lteff[js-1]);	//interpolate effective temperature
//      system.sec.lambda[n] = system.sec.track.lambda[js]-ratios*(system.sec.track.lambda[js]-system.sec.track.lambda[js-1]);	//interpolate lambda
      EbindperG = system.sec.track.m[js]*(system.sec.track.m[js]-system.sec.track.cm[js])/(system.sec.track.lambda[js]*system.sec.track.r[js])
                  -ratiop*(system.sec.track.m[js]*(system.sec.track.m[js]-system.sec.track.cm[js])/(system.sec.track.lambda[js]*system.sec.track.r[js])
                          -system.sec.track.m[js-1]*(system.sec.track.m[js-1]-system.sec.track.cm[js-1])/(system.sec.track.lambda[js-1]*system.sec.track.r[js-1]));	//interpolate current binding energy of the envelope
      if ((system.sec.m[n]-system.sec.cm[n]==0.0)&&(EbindperG==0.0)) system.sec.lambda[n] = -1.0e+99;	//no-interpolation
      else system.sec.lambda[n] = system.sec.m[n]*(system.sec.m[n]-system.sec.cm[n])/(system.sec.r[n]*EbindperG);	//get lambda from interpolated E_bind	// E_bind=G*M*M_env / lambda*R
      if(isnan(system.sec.lambda[n])){ screen = true; cerr << endl << "#system.sec.lambda[n]=" << system.sec.lambda[n] << " system.sec.m[n]=" << system.sec.m[n] << "Msun system.sec.cm[n]=" << system.sec.cm[n] << "Msun js=" << js << endl << "EbindperG=" << EbindperG << "Msun^2/Rsun=" << system.sec.track.m[js] << "*(" << system.sec.track.m[js] << "-" << system.sec.track.cm[js] << ")/(" << system.sec.track.lambda[js] << "*" << system.sec.track.r[js] << ")-" << ratiop << "*(" << system.sec.track.m[js] << "*(" << system.sec.track.m[js] << "-" << system.sec.track.cm[js] << ")/(" << system.sec.track.lambda[js] << "*" << system.sec.track.r[js] << ")-" << system.sec.track.m[js-1] << "*(" << system.sec.track.m[js-1] << "-" << system.sec.track.cm[js-1] << ")/(" << system.sec.track.lambda[js-1] << "*" << system.sec.track.r[js-1] << "))Msun^2/Rsun";}
      system.sec.last = js-1;	//remember position of the interpolation in the track array
      if (system.sec.r[n]<0) cerr << "#Error: system.sec.r[" << n << "]=" << system.sec.r[n] << "Rsun system.sec.track.r[" << js << "]=" << system.sec.track.r[js] << "Rsun system.sec.track.r[" << js-1 << "]=" << system.sec.track.r[js-1] << "Rsun ratio_s=" << ratios << endl;	//check for negative radius	//write error
      //determine values for the system
      system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
      system.qp[n] = system.prim.m[n]/system.sec.m[n];	//get mass ratio
      system.qs[n] = system.sec.m[n]/system.prim.m[n];	//get mass ratio
      system.a[n] = system.a[n]*(system.M[n-1])/(system.M[n]);	//get semi-major-axis, increased by wind mass loss; use simplyfied eq.(B6) in Soberman et al. (1997)
//      system.e[n] = system.e[n-1];	//get eccentricity
      system.P[n] = 2.0*M_PI*sqrt(pow(system.a[n],3)/(G*system.M[n]));	//get period out of semi-major-axis
      system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
      system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
//      system.peri[n] = system.a[n]*(1.0-system.e[n]);	//get peri-astro
      if ((((nextphase==12)&&(fabs(system.prim.r[n]-system.rp[n])>0.01))||((nextphase==21)&&(fabs(system.sec.r[n]-system.rs[n])>0.01)))&&(t>0.0)){	//check if the Roche-lobe filling star does so and if the time is in the future
        cerr << "#Error: RL-radii=(" << system.rp[n] << ", " << system.rs[n] << ")Rsun stellar radii=(" << system.prim.r[n] << ", " << system.sec.r[n] << ")Rsun ratiop=" << ratiop << " ratios=" << ratios << endl << "system.prim.track.t[" << jp-1 << "]=" << system.prim.track.t[jp-1] << "yr t=" << t << "yr system.prim.track.t[" << jp << "]=" << system.prim.track.t[jp] << "yr system.sec.track.t[" << js-1 << "]=" << system.sec.track.t[js-1] << "yr t=" << t << "yr system.sec.track.t[" << js << "]=" << system.sec.track.t[js] << "yr system.prim.m[" << n << "]=" << system.prim.m[n] << "Msun system.sec.m[" << n << "]=" << system.sec.m[n] << "Msun" << endl;	//write error
        screen = true;	//enable screen output
      }
    }
  }else if ((system.phase[n]==12)||(system.phase[n]==21)){	//Roche-lobe overflow
    if (screen && (output=='M')){
      if (system.phase[n]==12) cout << "RLO prim->sec: ";
      else cout << "RLO sec->prim: ";
    }
    nextphase = overflow(system);	//apply Roche-lobe overflow and get next kind of interaction
    system.phase[n] = nextphase;	//save next phase
    if (nextphase!=-1){	//determine current system/star values if system did not crash
      //determine values for the system
      system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
      system.qp[n] = system.prim.m[n]/system.sec.m[n];	//get mass ratio
      system.qs[n] = system.sec.m[n]/system.prim.m[n];	//get mass ratio
      system.e[n] = system.e[n-1];	//get eccentricity
//      circularize(system);	// circularization
      system.P[n] = 2.0*M_PI*sqrt(pow(system.a[n],3)/(G*system.M[n]));	//get period out of semi-major-axis
      system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
      system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
      system.peri[n] = system.a[n]*(1.0-system.e[n]);	//get peri-astro
      if ((nextphase%10!=3)&&((system.rp[n]<system.prim.r[n])||(system.rs[n]<system.sec.r[n]))) system.phase[n] = 4;	//if both stars fill their Roche-lobe --> merge
//      if (system.rp[n]<system.prim.r[n]) system.phase[n] = 12;	//if primary star fill its Roche-lobe --> Roche-lobe from primary to secondary
//      if (system.rs[n]<system.sec.r[n]) system.phase[n] = 21;	//if secondary star fill its Roche-lobe --> Roche-lobe from secondary to primary
    }
    if (nextphase%10!=3){	//if system does not go into common envelope
#pragma omp critical
{
      cRLO++;	//increase counter of Roche-lobe overflows
      if (system.phase[n-1]==12){
        if (system.prim.stage[n-1]==0) cRLOA++;	//increase counter of case A Roche-lobe overflows
        else if (system.prim.stage[n-1]==1) cRLOB++;	//increase counter of case B/C Roche-lobe overflows
        else if (system.prim.stage[n-1]==2) cRLOBB++;	//increase counter of case BB Roche-lobe overflows
      }else if (system.phase[n-1]==21){
        if (system.sec.stage[n-1]==0) cRLOA++;	//increase counter of case A Roche-lobe overflows
        else if (system.sec.stage[n-1]==1) cRLOB++;	//increase counter of case B/C Roche-lobe overflows
        else if (system.sec.stage[n-1]==2) cRLOBB++;	//increase counter of case BB Roche-lobe overflows
      }
      if (system.phase[n]==4) cMergerRLO++;	//increase counter of mergers by Roche-lobe overflows
}	//end critical
      if (screen && (output=='T')){
        cout << endl << "post RLO -- time: " << system.prim.t[n] << endl;
        cout << "  Mp Ms a: " << system.prim.m[n] << " " << system.sec.m[n] << " " << system.a[n] << endl;
      }
    }
  }else if (system.phase[n]%10==3){	//Common Envelope
    if (screen && (output=='M')) cout << "CE: ";
    nextphase = commonenvelope(system);	//apply Common Envelope and get next kind of interaction
    system.phase[n] = nextphase;	//save next phase
    if (nextphase!=-1){	//determine current system/star values if system did not crash
      //determine values for the system
      system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
      system.qp[n] = system.prim.m[n]/system.sec.m[n];	//get mass ratio
      system.qs[n] = system.sec.m[n]/system.prim.m[n];	//get mass ratio
      system.e[n] = system.e[n-1];	//get eccentricity
//      circularize(system);	// circularization
      system.P[n] = 2.0*M_PI*sqrt(pow(system.a[n],3)/(G*system.M[n]));	//get period out of semi-major-axis
      system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
      system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
      system.peri[n] = system.a[n]*(1.0-system.e[n]);	//get peri-astro
    }
#pragma omp critical
{
    cCE++;	//increase counter of common envelopes
    if (system.phase[n-1]==13){
      if (system.prim.stage[n-1]==0) cCEA++;	//increase counter of case A common envelopes
      else if (system.prim.stage[n-1]==1) cCEB++;	//increase counter of case B/C common envelopes
      else if (system.prim.stage[n-1]==2) cCEBB++;	//increase counter of case BB common envelopes
    }else if (system.phase[n-1]==23){
      if (system.sec.stage[n-1]==0) cCEA++;	//increase counter of case A common envelopes
      else if (system.sec.stage[n-1]==1) cCEB++;	//increase counter of case B/C common envelopes
      else if (system.sec.stage[n-1]==2) cCEBB++;	//increase counter of case BB common envelopes
    }
    if (system.phase[n]==4) cMergerCE++;	//increase counter of mergers by common envelopes
}	//end critical
    if (screen && (output=='T')){
      cout << "  Mp Ms a: " << system.prim.m[n-1] << " " << system.sec.m[n-1] << " " << system.a[n-1] << endl;
      cout << "post CE -- time: " << system.prim.t[n] << endl;
      cout << "  Mp Ms a: " << system.prim.m[n] << " " << system.sec.m[n] << " " << system.a[n] << endl;
    }
  }else if (system.phase[n]==4){	//merger
    if (screen){
      if (output=='M') cout << "merge: ";
      if (output=='T') cout << "merger" << endl;
    }
    nextphase = merge(system);	//apply merger and get next kind of interaction
    system.phase[n] = nextphase;	//save next phase
    if (nextphase!=-1){	//determine current system/star values if system did not crash
      //determine values for the system
      system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
      system.qp[n] = system.prim.m[n]/system.sec.m[n];	//get mass ratio
      system.qs[n] = system.sec.m[n]/system.prim.m[n];	//get mass ratio
      system.e[n] = system.e[n-1];	//get eccentricity
//      circularize(system);	// circularization
      system.P[n] = 2.0*M_PI*sqrt(pow(system.a[n],3)/(G*system.M[n]));	//get period out of semi-major-axis
      system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
      system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
      system.peri[n] = system.a[n]*(1.0-system.e[n]);	//get peri-astro
    }
#pragma omp critical
{
    cMerger++;	//increase counter of mergers
}	//end critical
  }else if (system.phase[n]%10==5){	//supernova/planetary nebula explosion
    if (screen && (output=='M')) cout << "supernova/planetary nebula: ";
    nextphase = supernova(system);	//apply supernova/planetary nebula and get next kind of interaction	//counter cplanetarynebula and csupernova increased within supernova()
    system.phase[n] = nextphase;	//save next phase
    if (nextphase!=-1){	//determine current system/star values if system did not crash
      //determine values for the system
      system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
      system.qp[n] = system.prim.m[n]/system.sec.m[n];	//get mass ratio
      system.qs[n] = system.sec.m[n]/system.prim.m[n];	//get mass ratio
//      circularize(system);	// circularization
      if (system.a[n]>=0.0) system.P[n] = 2.0*M_PI*sqrt(pow(system.a[n],3)/(G*system.M[n]));	//get period out of semi-major-axis
      else system.P[n] = -1.0e+99;
      system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
      system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
    }
    if (((system.prim.stage[n]==3)&&(system.prim.stage[n-1]==-3))||((system.sec.stage[n]==3)&&(system.sec.stage[n-1]==-3))){
      if (screen && (output=='T')){
        cout << "pre PN -- time: " << system.prim.t[n-1] << endl;
        cout << "  Mp Ms a: " << system.prim.m[n-1] << " " << system.sec.m[n-1] << " " << system.a[n-1] << endl;
        cout << "post PN -- time e: " << system.prim.t[n] << " " << system.e[n] << endl;
        cout << "  Mp Ms a: " << system.prim.m[n] << " " << system.sec.m[n] << " " << system.a[n] << endl;
      }
    }else{
      if (screen && (output=='T')){
        cout << "post SN -- time e: " << system.prim.t[n] << " " << system.e[n] << endl;
        cout << "  Mp Ms a: " << system.prim.m[n] << " " << system.sec.m[n] << " " << system.a[n] << endl;
      }
    }
#pragma omp critical
{
    if (system.phase[n]==4) cMergerSN++;	//increase counter of mergers by supernovae
}	//end critical
  }else{	//code crashed
    if (screen && (output=='M')) cout << "crash: ";
    system.phase[n] = -1;
  }

  if ((system.phase[n]>=0)&&(system.phase[n]<99)){	//if system evolves further
    if (screen && (output=='M')) cout << " => ";
    if ((system.phase[n]!=4)&&(system.phase[n]%10!=5)){	//not a merger and explosion not jet indicated
      if (system.prim.stage[n]==-3) system.phase[n] = 15;	//if primary star will explode to get a WD, NS or BH
      else if (system.sec.stage[n]==-3) system.phase[n] = 25;	//if secondary star will explode to get a WD, NS or BH
    }
  }
  if (system.prim.r[n]<0){	//check for non-negative radius
    cerr << endl << "#Error: system.prim.r[" << n << "]=" << system.prim.r[n] << "Rsun" << endl;	//write error
    screen = true;	//enable screen output
  }
  if (system.sec.r[n]<0){	//check for non-negative radius
    cerr << endl << "#Error: system.sec.r[" << n << "]=" << system.sec.r[n] << "Rsun" << endl;	//write error
    screen = true;	//enable screen output
  }
  if ((system.a[n]<0)&&(system.phase[n]!=4)&&(system.phase[n]!=-1)){	//check for non-negative semi-major-axis
    cerr << endl << "#Error: system.a[" << n << "]=" << system.a[n] << "Rsun system.phase[" << n << "]=" << system.phase[n] << endl;	//write error
    screen = true;	//enable screen output
  }
  if (((system.phase[n]!=system.phase[n-1])||(system.stagechange==1)||(system.prim.stage[n]==-3)||(system.sec.stage[n]==-3))&&(system.phase[n]>=0)&&(system.phase[n]<99)){
    newphase(system);	//prepare new phase
  }
}

int wind(t_system& system, double& t, double& r){
  int n=system.n-1;	//last entry index
  double trochep=1.0e+99,troches=1.0e+99;	//times for RLO form primary/secondary
  double rrochep=1.0e+99,rroches=1.0e+99;	//radius at RLO time form primary/secondary
  double rp,rs;	//RL-radius of primary/secondary at end of the evolutionary track if companion has fixed mass
  double tp=1.0e+99,ts=1.0e+99;	//times before stage change of primary/secondary

  rp = fmin(system.rp[n],rocherad(system.a[n]*(system.M[n])/(system.prim.track.m[system.prim.rmax]+system.sec.m[n])*(1.0-pow(system.e[n],2)),system.prim.track.m[system.prim.track.n-1]/system.sec.m[n]));
  rs = fmin(system.rs[n],rocherad(system.a[n]*(system.M[n])/(system.sec.track.m[system.sec.rmax]+system.prim.m[n])*(1.0-pow(system.e[n],2)),system.sec.track.m[system.sec.track.n-1]/system.prim.m[n]));
//  if (screen) cout << "Roche-lobe-check: system.prim.track.r[" << system.prim.rmax << "]=" << system.prim.track.r[system.prim.rmax] << "Rsun system.rp[" << n << "]=" << system.rp[n] << "Rsun system.sec.track.r[" << system.sec.rmax << "]=" << system.sec.track.r[system.sec.rmax] << "Rsun system.rs[" << n << "]=" << system.rs[n] << "Rsun" << endl;
  if (debug) cerr << "Roche-lobe-check: system.prim.track.r[" << system.prim.rmax << "]=" << system.prim.track.r[system.prim.rmax] << "Rsun rp=" << rp << "Rsun system.sec.track.r[" << system.sec.rmax << "]=" << system.sec.track.r[system.sec.rmax] << "Rsun rs=" << rs << "Rsun system.a[" << n << "]=" << system.a[n] << "Rsun system.e[" << n << "]=" << system.e[n] << " system.prim.m[" << n << "]=" << system.prim.m[n] << "Msun system.sec.m[" << n << "]=" << system.sec.m[n] << "Msun" << endl;
  //get maximal times
  if (system.prim.stage[n]==0) tp = system.prim.track.t[system.prim.track.TAMS];	//maximal time is at TAMS
  else tp = system.prim.track.t[system.prim.track.n-1];	//maximal time is at the end of the track
  if (system.sec.stage[n]==0) ts = system.sec.track.t[system.sec.track.TAMS];	//maximal time is at TAMS
  else ts = system.sec.track.t[system.sec.track.n-1];	//maximal time is at the end of the track

//  if ((system.prim.track.r[system.prim.rmax]>system.rp[n])||(system.sec.track.r[system.sec.rmax]>system.rs[n])){
  if ((system.prim.track.r[system.prim.rmax]>rp)||(system.sec.track.r[system.sec.rmax]>rs)){
//    if (system.prim.track.r[system.prim.rmax]>system.rp[n]){
    if (system.prim.track.r[system.prim.rmax]>rp){
      if (debug) cerr << "primary Roch-Lobe filling time:" << endl;
      trochep = RLFT(system,1,rrochep,rp);
    }
//    if (system.sec.track.r[system.sec.rmax]>system.rs[n]){
    if (system.sec.track.r[system.sec.rmax]>rs){
      if (debug) cerr << "secondary Roch-Lobe filling time:" << endl;
      troches = RLFT(system,2,rroches,rs);
    }
    if (((trochep==system.prim.t[n])||(troches==system.sec.t[n]))&&(system.sec.t[n]==0)) cerr << endl << "#Warning: Roche-lobe-overflow setup: prim.r[0]=" << system.prim.r[0] << "Rsun rp[0]=" << system.rp[0] << "Rsun sec.r[0]=" << system.sec.r[0] << "Rsun rs[0]=" << system.rs[0] << "Rsun" << endl;	//write warning
    if (((trochep==system.prim.t[n])||(troches==system.sec.t[n]))&& screen && (output=='M')) cout << "Instantaneous Roche-lobe overflow: ";
    if (screen&&((trochep<1.0e+99)||(troches<1.0e+99))){
      if (output=='M') cout << "Roche-lobe-filling-time:";
      if (trochep<1.0e+99){
        if (output=='M') cout << " t_roche_p=" << trochep << "yr r_roche_p=" << rrochep << "Rsun";
      }
      if (troches<1.0e+99){
        if (output=='M') cout << " t_roche_s=" << troches << "yr r_roche_s=" << rroches << "Rsun";
      }
    }
    if (trochep<troches){	//Roche-lobe overflow from primary
      t = trochep;
      r = rrochep;
      if ((trochep<1.0e+99)&&((trochep>tp)||(trochep>ts))){	//earlier checks failed, time is out of track
        cerr << endl << "#Error: t_roche_p=" << trochep << "yr tp=" << tp << "yr ts=" << ts << "yr" << endl;	//write error
        screen = true;	//enable screen output
        return -1;	//system crashed
      }
      return 12;
    }else if (trochep>troches){	//Roche-lobe overflow from secondary
      t = troches;
      r = rroches;
      if ((troches<1.0e+99)&&((troches>tp)||(troches>ts))){	//earlier checks failed, time is out of track
        cerr << endl << "#Error: t_roche_s=" << troches << "yr tp=" << tp << "yr ts=" << ts << "yr" << endl;	//write error
        screen = true;	//enable screen output
        return -1;	//system crashed
      }
      return 21;
    }else if (trochep==1.0e+99){	//no Roche-lobe overflow
      if ((tp<=ts)&&(system.stagechange==0)){	//primary ends its evolutionary track first
        t = tp;
        system.stagechange=nextstage(system.prim,n);	//end of evolutionary track = next stage
//        cout << " prim.stage[n]=" << system.prim.stage[n];
      }
      if ((tp>=ts)&&(system.stagechange==0)){	//secondary ends its evolutionary track first
        t = ts;
        system.stagechange=nextstage(system.sec,n);	//end of evolutionary track = next stage
//        cout << " sec.stage[n]=" << system.sec.stage[n];
      }
      if (debug) cerr << " t=" << t << "yr tp=" << tp << "yr ts=" << ts << "yr system.prim.stage[n]=" << system.prim.stage[n] << " system.sec.stage[n]=" << system.sec.stage[n] << endl;
      if (system.prim.stage[n]==-3) return 15;
      else if (system.sec.stage[n]==-3) return 25;
      else return 0;
    }else{	//trochep=troches: merger
      t = trochep;
      r = rrochep;
      return 4;
    }
  }else{	//no Roche-lobe overflow
    if (screen && (output=='M')) cout << "No Roche-lobe overflow possible";
    if ((tp<=ts)&&(system.stagechange==0)){	//primary ends its evolutionary track first
      t = tp;
      system.stagechange=nextstage(system.prim,n);	//end of evolutionary track = next stage
//      cout << " prim.stage[n]=" << system.prim.stage[n];
    }
    if ((tp>=ts)&&(system.stagechange==0)){	//secondary ends its evolutionary track first
      t = ts;
      system.stagechange=nextstage(system.sec,n);	//end of evolutionary track = next stage
//      cout << " sec.stage[n]=" << system.sec.stage[n];
    }
    if (debug) cerr << " t=" << t << "yr tp=" << tp << "yr ts=" << ts << "yr system.prim.stage[n]=" << system.prim.stage[n] << " system.sec.stage[n]=" << system.sec.stage[n] << endl;
    if (system.prim.stage[n]==-3) return 15;
    else if (system.sec.stage[n]==-3) return 25;
    else return 0;
  }
}

int overflow(t_system& system){
  int n=system.n-1;	//last entry index
  int jt;		//time index
  t_star* donor;	//donor
  t_star* accretor;	//accretor
  double Atidal=1.0;	// neglecting tidal effects - see Soberman et al (1997) eq.(14)
  double alpha=alphaRLO;	// direct wind mass loss from donor star during RLO 
  double beta=beta_const;	// all not accreted material transfered to secondary (non-deg.) star is
  double epsilon=1.0-alphaRLO-beta_const-delta;	//RLO efficency, eq.(B1) in Soberman et al. (1997)
  double A5,B5,C5;	//see Soberman et al. (1997)
//  double Aquadratic,Bquadratic,Cquadratic,x1,x2;	//0=Aquadratic*x^2+Bquadratic*x+Cquadratic -> x1,2=0.5*(-Bquadratic+-sqrt(Bquadratic^2-4*Aquadratic*Cquadratic))/Aquadratic
  double q=1.0,nq=1.0;	//old/new mass ratio=M_{donor}/M_{accretor}
  double tthermal;	//thermal/nuclear timescale of the donor	//,tnuclear
  double tthermal2;	//thermal timescale of the accretor
  double tRLO;	//time estimate for the RLO
  double Mdot=0.0,Medd=0.0;	//mass transfer rate and its Eddington limit
  double E_acc=0.0,E_nuc=0.0;	//specific energy of accretion/nuclear buring of transfered material
  double ratiot;	//ratio in time
  double newdonormass;	//the new mass of the donor after the mass transfer
  double Mcore,Menvelope,MenvelopeCO=0.03,Porbini;	//iron core mass after case BB RLO; helium/CO envelope mass after case BB RLO; orbital period after last event before case BB RLO
  double EbindperG=0.0;	//envelope binding energy for interpolation
  double McoreFeSi;	//FeSi core mass

  if (system.phase[n]==12){	//RLO from primaray to secondary
    donor = &system.prim;
    accretor = &system.sec;
    q = system.qp[n];	//q=M_donor/M_accretor
  }else if (system.phase[n]==21){	//RLO from secondary to primaray
    donor = &system.sec;
    accretor = &system.prim;
    q = system.qs[n];	//q=M_donor/M_accretor
  }else{	//undefined phase
    cerr << "#Error: overflow: system.phase[" << n << "]=" << system.phase[n] << endl;	//write error
    return -1;
  }
  if (debug) cout << endl << "convective envelope at " << convectiveenvelopefactor(donor[0].metal, donor[0].inimass)*donor[0].track.r[donor[0].rmax] << "Rsun" << endl;
//  if ((donor[0].stage[n]<2)&&(q<qlimit)&&(system.P[n]>0.008)&&(donor[0].r[n]<convectiveenvelopefactor(donor[0].metal, donor[0].inimass)*donor[0].track.r[donor[0].rmax])){	//RLO from hydrogen-rich star	//orbital period should be >~ 3day cf. Pablo Marchant	//evelope should be non-convective (30% for Z_MW and 70% for Z_IZw18)
  if ((donor[0].stage[n]<2)&&(q<qlimit)&&(system.P[n]>0.008)&&((donor[0].r[n]<convectiveenvelopefactor(donor[0].metal, donor[0].inimass)*donor[0].track.r[donor[0].rmax])||(q<1.5))){	//RLO from hydrogen-rich star	//orbital period should be >~ 3day cf. Pablo Marchant	//evelope should be non-convective (30% for Z_MW and 70% for Z_IZw18) or q<1.5
//    donor[0].RLO = true;	//transfer hydrogen
//    Medd = accretor[0].r[n]*Rsun*1.5E-12;	//Eddington-limit: Medd = R_acc/(10km) * 1.5E-8 M_sun/yr
    Medd = accretor[0].r[n]*4.0*M_PI*MH*c/(sigmaT);	//Eddington-limit: Medd = 4*pi*(mass of accreted object)*c/(count*Thomson cross section of radiative influenced particle per accreted object) here: accreted material=hydrogen, radiative influenced particle=electron
    if (accretor[0].stage[n]>2){	//accretor is compact object
      E_acc = G*accretor[0].m[n]/accretor[0].r[n];	//epsilon_acc = G*M/R
      if (accretor[0].stage[n]==5) E_acc *= 0.846;	//the accreation efficency of a black hole should be corrected from eta=0.5 to eta=0.423 -> factor 0.423/0.5=0.846
      if (accretor[0].stage[n]==3){	//accreted onto white drawf --> burn to helium
        E_nuc = (1.0-4.002603/(1.007825*4.0))*c*c;	//epsilon_nuc = deltaM/M*c^2 = (1-(m_{reaction product}*A_{particle})/(m_{particle}*A_{reaction product}))*c^2, here A_{particle}=1
      }else{	//accreted onto NS or black hole --> burn to iron
        E_nuc = (1.0-55.934940/(1.007825*56.0))*c*c;	//epsilon_nuc = deltaM/M*c^2 = (1-(m_{reaction product}*A_{particle})/(m_{particle}*A_{reaction product}))*c^2, here A_{particle}=1
      }
      if (debug) cerr << endl << "E_acc=" << E_acc/(1.0E+3*Msun*cgsEnergy) << "erg/g E_nuc=" << E_nuc/(1.0E+3*Msun*cgsEnergy) << "erg/g";
      Medd = 4.0*M_PI*G*accretor[0].m[n]*MH*c/(sigmaT*(E_acc+E_nuc));	//M_edd=L_edd/(epsilon_acc+epsilon_nuc) and L=4*pi*G*M_accretor*m_paricle*c/(count_{electrons per particle}*sigma_T), here count_{electrons per particle}=1
    }
    tthermal = G*pow(donor[0].m[n],2.0)/(donor[0].r[n]*pow(10.0,donor[0].llum[n])*Lsun);	//calculate thermal timescale: t_thermal=G(in Rsun^3/(Msun*yr^2))*M(in Msun)^2/(R(in Rsun)*L(in Msun*Rsun^2/yr^3)) in yr
    if (debug) cerr << endl << "G=" << G << "Rsun^3/(Msun*yr^2) m_donor=" << donor[0].m[n] << "Msun r_donor=" << donor[0].r[n] << "Rsun L_donor=" << pow(10.0,donor[0].llum[n]) << "Lsun tthermal=" << tthermal << "yr" << endl;
    tRLO = 3.0*tthermal;	// our best estimate in general
//    newdonormass = donor[0].cm[n];	//new mass will be former core mass
    for (jt=donor[0].last;jt<donor[0].track.n-1;jt++){	//find time index of time
      if (donor[0].track.t[jt]>donor[0].t[n]+tRLO) break;
    }
    ratiot = ratio(donor[0].t[n]+tRLO,donor[0].track.t[jt],donor[0].track.t[jt-1]);	//calculate ratio in time for the interpolation
    newdonormass = donor[0].track.cm[jt]-ratiot*(donor[0].track.cm[jt]-donor[0].track.cm[jt-1]);	//interpolate core mass
  }else if ((donor[0].stage[n]==2)&&(((accretor[0].stage[n]<3)&&(donor[0].m[max(n-1,0)]>Mhece)&&(Mhece>=0.0))||((accretor[0].stage[n]<3)&&(q<qlimit)&&(Mhece<0.0))||((accretor[0].stage[n]>=3)&&(system.P[n]/day>0.07)))){	//RLO from naked helium star, see Dewi & Pols (2003), Tauris et al.(2015)
//    donor[0].HeRLO = true;	//transfer helium
    if (donor[0].cm[n]<0.01) return 4;	//no significant core -> merge
//    Medd = accretor[0].r[n]*Rsun*3.0E-12;	//Eddington-limit: Medd = R_acc/10km * 2 * 1.5E-8 M_sun/yr
    Medd = accretor[0].r[n]*4.0*M_PI*4.002603/1.007825*MH*c/(2.0*sigmaT);	//Eddington-limit: Medd = 4*pi*(mass of accreted object)*c/(count*Thomson cross section of radiative influenced particle per accreted object) here: accreted object=helium, radiative influenced particle=2 electrons
    if (accretor[0].stage[n]>2){	//accretor is compact object
      E_acc = G*accretor[0].m[n]/accretor[0].r[n];	//epsilon_acc = G*M/R
      if (accretor[0].stage[n]==5) E_acc *= 0.846;	//the accreation efficency of a black hole should be corrected from eta=0.5 to eta=0.423 -> factor 0.423/0.5=0.846
      if (accretor[0].stage[n]==3){	//accreted onto white drawf --> helium accretion, no burning
        E_nuc = 0.0;	//epsilon_nuc = deltaM/M*c^2 = (1-(m_{reaction product}*A_{particle})/(m_{particle}*A_{reaction product}))*c^2, here A_{particle}=1
      }else{	//accreted onto NS or black hole --> burn to iron
        E_nuc = (1.0-55.934940*4.0/(4.002603*56.0))*c*c;	//epsilon_nuc = deltaM/M*c^2 = (1-(m_{reaction product}*A_{particle})/(m_{particle}*A_{reaction product}))*c^2, here A_{particle}=1
      }
      if (debug||screen) cerr << endl << "E_acc=" << E_acc/(1.0E+3*Msun*cgsEnergy) << "erg/g E_nuc=" << E_nuc/(1.0E+3*Msun*cgsEnergy) << "erg/g";
      Medd = 4.0*M_PI*G*accretor[0].m[n]*4.002603/1.007825*MH*c/(2.0*sigmaT*(E_acc+E_nuc));	//M_edd=L_edd/(epsilon_acc+epsilon_nuc) and L=4*pi*G*M_accretor*m_paricle*c/(count_{electrons per particle}*sigma_T), here count_{electrons per particle}=1
    }
/*    tnuclear = donor[0].t[n]-donor[0].t[max(n-2,0)];	//calculate nuclear timescale
    tRLO = 0.02*tnuclear;	// about 2% of t_nuclear He-star*/
    tthermal = G*pow(donor[0].m[n],2.0)/(donor[0].r[n]*pow(10.0,donor[0].llum[n])*Lsun);	//calculate thermal timescale: t_thermal=G(in Rsun^3/(Msun*yr^2))*M(in Msun)^2/(R(in Rsun)*L(in Msun*Rsun^2/yr^3)) in yr
    tRLO = tthermal;	// our best estimate in general	(Tauris in prep.)
//    alpha = 0.0;	// no wind massloss
    if (donor[0].inimass<=3.5){
      Porbini = 2.0*M_PI*sqrt(pow(donor[0].r[n]/rocherad(1.0,donor[0].m[n]/1.35)*(donor[0].m[n]+1.35)/(donor[0].inimass+1.35),3)/(G*(donor[0].inimass+1.35)))/day;	//initial orbital period in days: get seperation if companion is a 1.35Msun NS, undo wind mass loss to initial mass, get initial period
      Mcore = (1/(400.0*Porbini)+0.49)*donor[0].inimass-(0.016/Porbini-0.106);	//Tauris,Langer,Podsiadlowski 2015 eq.(3)
      Mcore = fmax(donor[0].cm[n],Mcore);	//final core is at least core at onset of RLO
      if (Porbini<2.0) Menvelope = 0.18*pow(Porbini,0.45)*(log(pow(Mcore,4))-1.05);	//Tauris,Langer,Podsiadlowski 2015 eq.(4)
      else Menvelope = Mcore*(log(pow(Porbini,-0.2))+1.0)+log(sqrt(Porbini))-1.6;	//Tauris,Langer,Podsiadlowski 2015 eq.(4), corrected by -0.1Msun because of Fig.16
      if (Menvelope<accuracy){	//avoid negative total envelope masses
        MenvelopeCO = fmax(accuracy,0.03-accuracy+Menvelope);
        Menvelope = accuracy;
      }
      if (screen && (output=='M')) cerr << " Porbini=" << Porbini << "day McoreCO=" << Mcore << "Msun MenvelopeHe=" << Menvelope << "Msun MenvelopeCO=" << MenvelopeCO << "Msun ";
      newdonormass = fmin(Mcore+Menvelope+MenvelopeCO,fmax(donor[0].m[n-1]-0.1,Mcore));	//Tauris,Langer,Podsiadlowski 2015: M_total=McoreCO+MenvelopeHe+MenvelopeCO with MenvelopeCO=0.03 Msun, minimum lost mass in case BB RLO is set to 0.1Msun while core is smaller
      donor[0].cm[n] = Mcore;	//set new core mass
    }else{
      //see supernova
      McoreFeSi = 1.05667*pow(donor[0].cm[n],0.453707);	//fit to Timmes et al. 1995
      newdonormass = (-1.0+sqrt(1.0+0.336*McoreFeSi))/0.168;	//solution of M_NS=M_core-M_NS^bind from Tauris,Langer,Podsiadlowski 2015 where M_NS^bind=0.084*M_NS^2/Msun from Lattimer & Yahil 1989 eq.(10) by converting energy to mass in Msun

      newdonormass = fmax(donor[0].cm[n]+0.1,newdonormass);	//new mass will be former core mass+0.1 to correct for burning during the mass transfer phase
    }
  }else{	//common envelope
/*    if (donor[0].stage[n]<2) donor[0].RLO = true;	//lose hydrogen
    else if (donor[0].stage[n]==2) donor[0].HeRLO = true;	//lose helium*/
    if (system.phase[n]==12) return 13;
    else if (system.phase[n]==21) return 23;
    else return -1;	//undefined phase
  }
  if (screen && (output=='T')){
    if (system.phase[n-1]==12) cout << "pre RLO p->s -- time(yr) Rl(Rsun): " << system.prim.t[n-1] << " " << system.rp[n-1] << endl;
    else cout << "pre RLO s->p -- time Rl: " << system.sec.t[n-1] << " " << system.rs[n-1] << endl;
    cout << "  Mp Ms a: " << system.prim.m[n-1] << " " << system.sec.m[n-1] << " " << system.a[n-1] << endl << "  ";
  }
  //check that life time is not exceeded
  if (donor[0].t[n]+tRLO>donor[0].track.t[donor[0].track.n-1]) tRLO = donor[0].track.t[donor[0].track.n-1]-donor[0].t[n]-0.001;	//maximal time for donor
  if (accretor[0].t[n]+tRLO>accretor[0].track.t[accretor[0].track.n-1]) tRLO = accretor[0].track.t[accretor[0].track.n-1]-accretor[0].t[n]-0.001;	//maximal time for accretor
  epsilon=1.0-alpha-beta-delta;	//RLO efficency, eq.(B1) in Soberman et al. (1997)
  Mdot = (donor[0].m[n]-newdonormass)/tRLO;
  if (debug||screen) cerr << endl << "Medd=" << Medd << "Msun/yr Mdot=" << Mdot << "Msun/yr=(" << donor[0].m[n] << "-" << newdonormass << ")/" << tRLO << "Msun/yr epsilon=" << epsilon << " n=" << n;
  if ((Mdot*epsilon>Medd)){	// check if all transfered material can be accreted onto remnant	&&(accretor[0].stage[n]>2)
    epsilon = Medd/Mdot;
    beta = 1.0-alpha-delta-epsilon;
    if (screen && (output=='M')) cout << "transfer at Eddington-limit: Medd=" << Medd << "Msun/yr Mdot=" << Mdot << "Msun/yr ";
    if (accretor[0].stage[n]<=2){	//common envelope
      if (system.phase[n]==12) return 13;
      else if (system.phase[n]==21) return 23;
      else return -1;	//undefined phase
    }
  }
  if (epsilon>1.0-1.0e-5){	// for smooth calculations
    beta += epsilon-(1.0-1.0e-5);	// conserve epsilon=1.0-alpha-beta-delta
    epsilon = 1.0-1.0e-5;
  }
  donor[0].t[n] += tRLO;	//add time for the overflow
  accretor[0].t[n] += tRLO;	//add time for the overflow
  nq = newdonormass/accretor[0].m[n];	//new mass ratio
  if (screen) cout << " Atidal=" << Atidal << " alpha=" << alpha << " Gamma=" << Gamma << " delta=" << delta << " beta=" << beta << " epsilon=" << epsilon << " nq=" << nq << " q=" << q << " tRLO=" << tRLO << "yr";
  A5 = Atidal*alpha+Gamma*delta;	// eq.(B8) in Soberman et al. (1997)
  B5 = (Atidal*alpha+beta)/(1.0-epsilon);	// eq.(B9) in Soberman et al. (1997)
  C5 = Gamma*delta*(1.0-epsilon)/epsilon + Atidal*alpha*epsilon/(1.0-epsilon) + beta/(epsilon*(1.0-epsilon));	// eq.(B10) in Soberman et al. (1997)
/*  //get nearest orbit during RLO from minimum in eq.(B6) in Soberman et al. (1997) with A5,B5,C5>=0: x3=donor[0].m[n] (for A5>1.5), x4=-accretor[0].m[n]/epsilon (for B5>1.5), x1 and x2 from quadratic equation
  Aquadratic = -epsilon*(epsilon-1.0)*(2.0*C5+3.0);
  Bquadratic = accretor[0].m[n]*(2.0*(epsilon-1.0)*A5+2.0*B5-2.0*epsilon*C5-5.0*epsilon+1.0)+epsilon*donor[0].m[n]*(2.0*(epsilon-1.0)*A5+2.0*B5-4.0*C5+2.0*epsilon*C5+epsilon-5.0);
  Cquadratic = 2.0*accretor[0].m[n]*accretor[0].m[n]*(A5-1.0)+accretor[0].m[n]*donor[0].m[n]*(2.0*(1.0+epsilon)*A5-2.0*B5+2.0*epsilon*C5+epsilon-1.0)+2.0*epsilon*donor[0].m[n]*donor[0].m[n]*(A5-B5+C5+1.0);
  x1 = 0.5*(-Bquadratic+sqrt(Bquadratic*Bquadratic-4*Aquadratic*Cquadratic))/Aquadratic;	//solution 1
  x2 = 0.5*(-Bquadratic-sqrt(Bquadratic*Bquadratic-4*Aquadratic*Cquadratic))/Aquadratic;	//solution 2
  if (debug||screen) cout << " x1=" << x1 << "Msun x2=" << x2 << "Msun";
  if ((x1>0)&&(newdonormass+x1<donor[0].m[n])){
    nq = (donor[0].m[n]-x1)/(accretor[0].m[n]+x1*epsilon);	//new mass ratio
    if (debug||screen) cout << " rochemin=" << rocherad(system.a[n]*(pow((nq/q),(2.0*A5-2.0))*pow((1.0+nq)/(1.0+q),(1.0-2.0*B5))*pow((1.0+epsilon*nq)/(1.0+epsilon*q),(2.0*C5+3.0))),1.0/nq) << "Rsun";
    if (accretor[0].r[n]>rocherad(system.a[n]*(pow((nq/q),(2.0*A5-2.0))*pow((1.0+nq)/(1.0+q),(1.0-2.0*B5))*pow((1.0+epsilon*nq)/(1.0+epsilon*q),(2.0*C5+3.0))),1.0/nq)) return 4;
  }
  if ((x2>0)&&(newdonormass+x2<donor[0].m[n])){
    nq = (donor[0].m[n]-x2)/(accretor[0].m[n]+x2*epsilon);	//new mass ratio
    if (debug||screen) cout << " rochemin=" << rocherad(system.a[n]*(pow((nq/q),(2.0*A5-2.0))*pow((1.0+nq)/(1.0+q),(1.0-2.0*B5))*pow((1.0+epsilon*nq)/(1.0+epsilon*q),(2.0*C5+3.0))),1.0/nq) << "Rsun";
    if (debug) cout << " rochemin=" << accretor[0].r[n]/(pow((nq/q),(2.0*A5-2.0))*pow((1.0+nq)/(1.0+q),(1.0-2.0*B5))*pow((1.0+epsilon*nq)/(1.0+epsilon*q),(2.0*C5+3.0))) << "Rsun";
  }*/
  accretor[0].m[n] += (donor[0].m[n]-newdonormass)*epsilon;	//assume donor lost all its envelope => add envelope mass * efficency to accretormass
// add this also to core mass
  if (accretor[0].stage[n]>2) accretor[0].cm[n] += (donor[0].m[n]-newdonormass)*epsilon;	//if the transfer is onto a remnant
  else{
    tthermal2 = G*pow(accretor[0].m[n-1],2.0)/(accretor[0].r[n-1]*pow(10.0,accretor[0].llum[n-1])*Lsun);	//calculate thermal timescale: t_thermal=G(in Rsun^3/(Msun*yr^2))*M(in Msun)^2/(R(in Rsun)*L(in Msun*Rsun^2/yr^3)) in yr
    if ((donor[0].stage[n]==2)&&(accretor[0].stage[n]<1)) accretor[0].cm[n] += (donor[0].m[n]-newdonormass)*tRLO/tthermal2*epsilon;	//if helium is transfered to a star with a helium core; the mass which adds to the core is limited by the thermal time scale of the accretor
    if (donor[0].stage[n]>2) accretor[0].cm[n] += (donor[0].m[n]-newdonormass)*tRLO/tthermal2*epsilon;	//if heavier elements are transfered; the mass which adds to the core is limited by the thermal time scale of the accretor
  }
  donor[0].m[n] = newdonormass;	//assume donor lost all its envelope => new mass = previous core mass
  system.a[n] *= pow((nq/q),(2.0*A5-2.0))*pow((1.0+nq)/(1.0+q),(1.0-2.0*B5))*pow((1.0+epsilon*nq)/(1.0+epsilon*q),(2.0*C5+3.0));	// eq.(B6) in Soberman et al. (1997)
  if ((accretor[0].stage[n]>=0)&&(accretor[0].stage[n]<=2)){	//calculate new track of the accretor
    if ((donor[0].stage[n]<2)&&(accretor[0].stage[n]==2)){	//hydrogen is accreted onto a naked helium star
      accretor[0].stage[n] = 0;	//get back to a hydrogen burning star
      accretor[0].cm[n] = accretor[0].m[n-1]; //old mass gets new He-core mass
      updatetrack(accretor[0],n);
/*      if (accretor[0].stage[n]<0){	//no hydrogen burning track -> restore old values
        accretor[0].m[n] = accretor[0].m[n-1];	//take old mass
        accretor[0].cm[n] = accretor[0].cm[n-1];	//take old core mass
        accretor[0].stage[n] = 2;	//take old stage
      }*/
    }else{
      updatetrack(accretor[0],n);
    }
//    if (screen) cout << "accretor[0].m[" << n << "]=" << accretor[0].m[n] << "Msun accretor[0].track.m[" << accretor[0].last << "]=" << accretor[0].track.m[accretor[0].last] << "Msun" << endl << "accretor[0].cm[" << n << "]=" << accretor[0].cm[n] << "Msun accretor[0].track.cm[" << accretor[0].last << "]=" << accretor[0].track.cm[accretor[0].last] << "Msun" << endl << "accretor[0].t[" << n << "]=" << accretor[0].t[n] << "yr accretor[0].track.t[" << accretor[0].last << "]=" << accretor[0].track.t[accretor[0].last] << "yr" << endl;
    if ((accretor[0].m[n]==accretor[0].m[n-1])&&(accretor[0].r[n]==accretor[0].r[n-1])&&(accretor[0].cm[n]==accretor[0].cm[n-1])){	//no trackupdate
      for (jt=accretor[0].last;jt<accretor[0].track.n-1;jt++){	//find time index of time
        if (accretor[0].track.t[jt]>accretor[0].t[n]) break;
      }
      ratiot = ratio(accretor[0].t[n],accretor[0].track.t[jt],accretor[0].track.t[jt-1]);	//calculate ratio in time for the interpolation
      accretor[0].m[n] = accretor[0].track.m[jt]-ratiot*(accretor[0].track.m[jt]-accretor[0].track.m[jt-1]);	//interpolate mass
      accretor[0].r[n] = accretor[0].track.r[jt]-ratiot*(accretor[0].track.r[jt]-accretor[0].track.r[jt-1]);	//interpolate radius
      accretor[0].cm[n] = accretor[0].track.cm[jt]-ratiot*(accretor[0].track.cm[jt]-accretor[0].track.cm[jt-1]);	//interpolate core mass
      accretor[0].llum[n] = accretor[0].track.llum[jt]-ratiot*(accretor[0].track.llum[jt]-accretor[0].track.llum[jt-1]);	//interpolate luminosity
      accretor[0].lteff[n] = accretor[0].track.lteff[jt]-ratiot*(accretor[0].track.lteff[jt]-accretor[0].track.lteff[jt-1]);	//interpolate effective temperature
//      accretor[0].lambda[n] = accretor[0].track.lambda[jt]-ratiot*(accretor[0].track.lambda[jt]-accretor[0].track.lambda[jt-1]);	//interpolate lambda
      EbindperG = accretor[0].track.m[jt]*(accretor[0].track.m[jt]-accretor[0].track.cm[jt])/(accretor[0].track.lambda[jt]*accretor[0].track.r[jt])
                  -ratiot*(accretor[0].track.m[jt]*(accretor[0].track.m[jt]-accretor[0].track.cm[jt])/(accretor[0].track.lambda[jt]*accretor[0].track.r[jt])
                          -accretor[0].track.m[jt-1]*(accretor[0].track.m[jt-1]-accretor[0].track.cm[jt-1])/(accretor[0].track.lambda[jt-1]*accretor[0].track.r[jt-1]));	//interpolate current binding energy of the envelope
      if ((accretor[0].m[n]-accretor[0].cm[n]==0.0)&&(EbindperG==0.0)) accretor[0].lambda[n] = -1.0e+99;	//no interpolation
      else accretor[0].lambda[n] = accretor[0].m[n]*(accretor[0].m[n]-accretor[0].cm[n])/(EbindperG*accretor[0].r[n]);	//get lambda from interpolated Ebind	// E_bind=G*M*M_env / lambda*R
      if(isnan(accretor[0].lambda[n])){ screen = true; cerr << endl << "#accretor[0].lambda[n]=" << accretor[0].lambda[n];}
      accretor[0].last = jt-1;	//remember position of the interpolation in the track array
    }else{
      accretor[0].r[n] = accretor[0].track.r[accretor[0].last];
      accretor[0].llum[n] = accretor[0].track.llum[accretor[0].last];
      accretor[0].lteff[n] = accretor[0].track.lteff[accretor[0].last];
      accretor[0].lambda[n] = accretor[0].track.lambda[accretor[0].last];
    }
  }else if ((accretor[0].stage[n]==3)&&(accretor[0].cm[n]<WDmax)){	//WD below Chandrasekhar limit(stable WD)
    accretor[0].r[n] = WDradius(accretor[0].cm[n]);	//remnant radius
    system.stagechange = max(system.stagechange,nextstage(accretor[0],n));
  }else if ((accretor[0].stage[n]==4)&&(accretor[0].cm[n]<=NSmax)){	//NS below SN-limit(stable NS)
    accretor[0].r[n] = NSradius(accretor[0].cm[n]);	//remnant radius
    system.stagechange = max(system.stagechange,nextstage(accretor[0],n));
  }else if (accretor[0].stage[n]==5){	//stable BH
    accretor[0].r[n] = Schwarzschildradius(accretor[0].cm[n]);	//remnant radius
    system.stagechange = max(system.stagechange,nextstage(accretor[0],n));
  }else{	//accretor explodes
    accretor[0].stage[n] = -3;
  }
  if (donor[0].stage[n]<2){	//star become naked helium star
    donor[0].stage[n] = -2;
    system.stagechange = max(system.stagechange,nextstage(donor[0],n));
    donor[0].r[n] = donor[0].track.r[donor[0].last];
    donor[0].cm[n] = donor[0].track.cm[donor[0].last];
//    if (screen) cout << "donor.track.mass=" << donor[0].track.m[donor[0].last] << "Msun donor.mass=" << donor[0].m[n] << "Msun";
    donor[0].llum[n] = donor[0].track.llum[donor[0].last];
    donor[0].lteff[n] = donor[0].track.lteff[donor[0].last];
    donor[0].lambda[n] = donor[0].track.lambda[donor[0].last];
    while ((donor[0].t[n]>donor[0].track.t[donor[0].track.n-1])&&(donor[0].stage[n]!=-3)) system.stagechange = max(system.stagechange,nextstage(donor[0],n));	//check if donor gets older than its current evolutionary track
    return 0;
  }else{	//donor explodes
    donor[0].stage[n] = -3;
    donor[0].r[n] = WDradius(donor[0].m[n]);
  }
  if ((donor[0].stage[n]==-3)&&(system.phase[n]==12)) return 15;		//primary explodes
  else if ((donor[0].stage[n]==-3)&&(system.phase[n]==21)) return 25;		//secondary explodes
  else if ((accretor[0].stage[n]==-3)&&(system.phase[n]==21)) return 15;	//primary explodes
  else if ((accretor[0].stage[n]==-3)&&(system.phase[n]==12)) return 25;	//secondary explodes
  else return -1;
}

int commonenvelope(t_system& system){
  int n=system.n-1;	//last entry index
  int i,j;	//index variables
  int ilow=0,iup=0,jlow=0,jup=0;	//index variables: ilow/iup=mass index, jlow/jup=radius index
  t_star* star;	//star
  t_star* companion;	//companion
  double mratio,ratiolow,ratioup,tratio;	//mass ratio, ratios at lower/upper mass track, time ratio
  double lambda=0.0,lambdalow,lambdaup;	//lambda value, lambda values at lower/upper mass track
  double tCE=1000.0;	//time estimate for the common envelope: assume 1000yr
  double E_acc=0.0,E_nuc=0.0,E_bind=1.0;	//E_acc: releast energy by accreting material onto a remnant, E_nuc: nuclear burning energy of accreted material, E_bind: binding energy
  double M_edd=0.0;	//Eddington limited transfered mass
  double A_preburning=1.0, m_preburning=1.007825, A_postburning=56.0, m_postburning=55.934940, c_preburning=1.0;	//atomic number and mass before and after burning (initial: hyrogen -> iron); count of electrons per particle (hydrogen)
  double EbindperG=0.0;	//envelope binding energy for interpolation

  if (system.rp[n]<=system.prim.r[n]*(1.0+1.0e-4)){	//primaray fills Roche-lobe
    star = &system.prim;
    companion = &system.sec;
  }else if (system.rs[n]<=system.sec.r[n]*(1.0+1.0e-4)){	//secondary fills Roche-lobe
    star = &system.sec;
    companion = &system.prim;
  }else if (n>1){	//find star which filled its Roche-lobe last
    if (system.phase[n-1]==13){	//primaray filled Roche-lobe before
      star = &system.prim;
      companion = &system.sec;
    }else if (system.phase[n-1]==23){	//secondary filled Roche-lobe before
      star = &system.sec;
      companion = &system.prim;
    }else{	//no star filled Roche-lobe before
      cerr << "#Error: common envelope: system.rp[" << n << "]=" << system.rp[n] << "Rsun system.prim.r[" << n << "]=" << system.prim.r[n] << "Rsun system.rs[" << n << "]=" << system.rs[n] << "Rsun system.sec.r[" << n << "]=" << system.sec.r[n] << "Rsun" << endl;	//write error
      screen = true;	//enable screen output
      return -1;
    }
  }else{	//no star fills Roche-lobe
    cerr << "#Error: common envelope: system.rp[" << n << "]=" << system.rp[n] << "Rsun system.prim.r[" << n << "]=" << system.prim.r[n] << "Rsun system.rs[" << n << "]=" << system.rs[n] << "Rsun system.sec.r[" << n << "]=" << system.sec.r[n] << "Rsun" << endl;	//write error
    screen = true;	//enable screen output
    return -1;
  }
  if (star[0].m[n]==star[0].cm[n]) return 4;	//no envelope to eject -> merge
//  if ((star[0].stage[n]==2)||(companion[0].stage[n]==2)){	// He-star CE
//  if (star[0].stage[n]==2){	// He-star CE
//    lambda = lambda_const;
//  }else{
    if (alphaTH<0.0){
      iup = lambdaarray.n-1;
      while (lambdaarray.m[iup-1]==lambdaarray.m[lambdaarray.n-1]) iup--;
      jup = lambdaarray.n-iup;
      for (i=0;i<lambdaarray.n;i++){
        if (lambdaarray.m[i]<star[0].inimass){	//find largest mass which is lower than initial mass of the star
          if (lambdaarray.m[i]==lambdaarray.m[ilow]){
            jlow++;	//increse number of the lower mass
          }else{
            ilow = i;	//indicatest the first position of the lower mass
            jlow = 1;	//set number of the lower mass
          }
        }else{	//find mass which is larger than initial mass of the star
          if (jup==lambdaarray.n-iup){
            iup = i;	//indicatest the first position of the upper mass
            jup = 1;	//set number of the upper mass
          }else if (lambdaarray.m[i]==lambdaarray.m[iup]){
            jup++;	//increse number of the upper mass
          }else{
            break;
          }
        }
      }
      mratio = ratio(star[0].inimass,lambdaarray.m[iup],lambdaarray.m[ilow]);	//calculate ratio between masses
      for (j=1;j<jlow;j++){
        if (lambdaarray.r[ilow+j]>star[0].r[n]){	//find radius of the star at lower mass
          break;
        }
      }
      if (j==jlow) j--;	//extrapolate to larger radii
      jlow = j;
      ratiolow = ratio(star[0].r[n],lambdaarray.r[ilow+jlow],lambdaarray.r[ilow+jlow-1]);	//calculate ratio in radius at lower mass
      lambdalow = lambdaarray.lambda[ilow+jlow]-ratiolow*(lambdaarray.lambda[ilow+jlow]-lambdaarray.lambda[ilow+jlow-1]);	//interpolate/extrapolate lambda at lower mass
      for (j=1;j<jup;j++){
        if (lambdaarray.r[iup+j]>star[0].r[n]){	//find radius of the star at upper mass
          break;
        }
      }
      if (j==jup) j--;	//extrapolate to larger radii
      jup = j;
      ratioup = ratio(star[0].r[n],lambdaarray.r[iup+jup],lambdaarray.r[iup+jup-1]);	//calculate ratio in radius at upper mass
      lambdaup = lambdaarray.lambda[iup+jup]-ratioup*(lambdaarray.lambda[iup+jup]-lambdaarray.lambda[iup+jup-1]);	//interpolate/extrapolate lambda at upper mass
      lambda = lambdaup-mratio*(lambdaup-lambdalow);	//interpolate/extrapolate lambda between masses
      lambda = 2.5*lambda; // APPROX! to account for M_core being a little larger than defined by X < 0.10
    }else{
      lambda = star[0].lambda[n];
    }
//  }
  if (screen){
    if (output=='M') cout << "alpha_th=" << alphaTH << " lambda=" << lambda;
    if (output=='T'){
      if (star[0].stage[n]==2) cout << "pre He-star";
      else cout << "pre";
      cout << " CE -- time lambda Rl: " << star[0].t[n-1] << " " << lambda << " " << star[0].r[n-1] << endl;
    }
  }
  if (lambda<=0.0){	//extrapolated lambda which lead to zero or negative semi-major-axis
    system.a[n] = 1.0e-99;
  }else{	//determine semi-major-axis change
    system.a[n] *= (fmax(star[0].cm[n],1.0e-10)*companion[0].m[n]/star[0].m[n])/(companion[0].m[n]+2.0*(star[0].m[n]-star[0].cm[n])/(alphaCE*lambda*star[0].r[n]/system.a[n]));	//consider a minimal mass of 10^{-10}Msun as left remnant
  }
  if (star[0].t[n]+tCE>=star[0].track.t[star[0].track.n-1]) tCE = star[0].track.t[star[0].track.n-1]-star[0].t[n];	//maximal time for donor	//(1.0+accuracy)*
  if (companion[0].t[n]+tCE>=companion[0].track.t[companion[0].track.n-1]) tCE = companion[0].track.t[companion[0].track.n-1]-companion[0].t[n];	//maximal time for accretor	//(1.0+accuracy)*
  if (companion[0].stage[n]>2){	//companion is compact object
//    E_bind = 6.67E-8*star[0].m[n]*(star[0].m[n]-star[0].cm[n])*pow(1.989E+33,2)/(lambda*star[0].r[n]*6.96E+10); // E_bind=G*M*M_env / lambda*R
    E_bind = G*star[0].m[n]*(star[0].m[n]-star[0].cm[n])/(lambda*star[0].r[n])/cgsEnergy; // E_bind=G*M*M_env / lambda*R
// APPROX to include the released grav. bind. energy of material accreted onto NS/BH in CE: E_acc=tau*L_edd.
//  where L_edd for hydrogen is given by e.g. Heuvel, Saas-Fee p.313 or King "Acc. Power")
    if (star[0].stage[n]==2){	// accretion of helium
//      E_acc = companion[0].m[n]*7.6E+48;	//tau=1000 yr
      A_preburning = 4.0;
      m_preburning = 4.002603;
      c_preburning = 2.0;
    }else{	// accretion of hydrogen
//      E_acc = companion[0].m[n]*3.8E+48;	//tau=1000 yr
      A_preburning = 1.0;
      m_preburning = 1.007825;
      c_preburning = 1.0;
    }
    if (companion[0].stage[n]==3){	//accreted onto white drawf --> burn to helium
      A_postburning = 4.0;
      m_postburning = 4.002603;
    }else{	//accreted onto NS or black hole --> burn to iron
      A_postburning = 56.0;
      m_postburning = 55.934940;
    }
    E_acc = G*companion[0].m[n]/companion[0].r[n];	//epsilon_acc = G*M/R
    if (companion[0].stage[n]==5) E_acc *= 0.846;	//the accreation efficency of a black hole should be corrected from eta=0.5 to eta=0.423 -> factor 0.423/0.5=0.846
    E_nuc = (1.0-m_postburning*A_preburning/(m_preburning*A_postburning))*c*c;	//epsilon_nuc = deltaM/M*c^2 = (1-(m_{reaction product}*A_{particle})/(m_{particle}*A_{reaction product}))*c^2
    if (debug) cerr << endl << "E_acc=" << E_acc/(1.0E+3*Msun*cgsEnergy) << "erg/g E_nuc=" << E_nuc/(1.0E+3*Msun*cgsEnergy) << "erg/g";
    M_edd = 4.0*M_PI*G*companion[0].m[n]*m_preburning/1.007825*MH*c/(c_preburning*sigmaT*(E_acc+E_nuc))*tCE;	//M_edd=Mdot*t_CE with Mdot=L_edd/(epsilon_acc+epsilon_nuc) and L=4*pi*G*M_accretor*m_paricle*c/(count_{electrons per particle}*sigma_T)
    E_acc *= M_edd/cgsEnergy;	//get energy procuced by Eddington accretion in erg
    E_nuc *= M_edd/cgsEnergy;	//get energy procuced by Eddington accretion in erg
/*    M_edd = 4.0*M_PI*companion[0].r[n]*m_preburning/1.007825*MH*c/(sigmaT)*tCE;	//M_edd=Mdot*t_CE with Mdot=L_edd*R_accretor/(G*M_accretor)=4*pi*R_accretor*m_paricle*c/(count_{electrons per particle}*sigma_T)
    E_acc = 4.0*M_PI*G*companion[0].m[n]*m_preburning/1.007825*MH*c/(sigmaT)*tCE/cgsEnergy;	//E_acc=L_edd*t_CE with L_edd=4*pi*G*M_accretor*m_paricle*c/(count_{electrons per particle}*sigma_T)
    E_nuc = (1.0-m_postburning*A_preburning/(m_preburning*A_postburning))*M_edd*c*c/cgsEnergy;	//E_nuc=deltaM*c*c with deltaM=M_edd^particle-M_edd^{reaction product} where M_edd^{reaction product}=n_{reaction product}*m_{reaction product} and n_{reaction product}=n_particle*A_particle/A_{reaction product} ==> deltaM=(1-(m_{reaction product}*A_particle)/(m_particle*A_{reaction product}))*M_edd^particle
    if (companion[0].stage[n]==5) E_acc *= 0.846;	//the accreation efficency of a black hole should be corrected from eta=0.5 to eta=0.423 -> factor 0.423/0.5=0.846*/
    if (screen && (output=='M')) cout << " M_acc_edd=" << M_edd << "Msun E_acc=" << E_acc << "erg E_nuc=" << E_nuc << "erg E_bind=" << E_bind << "erg";
    if (screen && (output=='T')) cout << "  M_acc_edd E_acc E_nuc E_bind: " << M_edd << " " << E_acc << " " << E_nuc << " " << E_bind << endl;
    if (M_edd>star[0].m[n]-star[0].cm[n]){	//if Eddington limit is higher than mass of the common envelope scale down to mass of the common envelope
      E_acc *= (star[0].m[n]-star[0].cm[n])/M_edd;
      E_nuc *= (star[0].m[n]-star[0].cm[n])/M_edd;
      tCE *= (star[0].m[n]-star[0].cm[n])/M_edd;
      M_edd = star[0].m[n]-star[0].cm[n];
      if ((screen)&&(output=='M')) cout << " super Eddington accretion limited";
    }
//  system.a[n] = system.a[n]*(1.0+E_acc/E_bind); // simple interpolation formula, get maximal the initial separation
    //if E_acc>E_bind: semi-major-axis cannot exceed the increase by wind mass loss; use simplyfied eq.(B6) in Soberman et al. (1997) as upper limit
    system.a[n] = fmin(system.a[n]*(1.0+(E_acc+E_nuc)/E_bind),system.a[n-1]*(companion[0].m[n-1]+star[0].m[n-1])/(companion[0].m[n]+star[0].cm[n])); // simple interpolation formula, get maximal the initial separation with wind mass loss
    //add accreted mass to companion
    companion[0].m[n] += M_edd;
    companion[0].cm[n] += M_edd;
    if ((companion[0].stage[n]==3)&&(companion[0].cm[n]<WDmax)){	//WD below Chandrasekhar limit(stable WD)
      companion[0].r[n] = WDradius(companion[0].cm[n]);	//remnant radius
      system.stagechange = max(system.stagechange,nextstage(companion[0],n));
    }else if ((companion[0].stage[n]==4)&&(companion[0].cm[n]<=NSmax)){	//NS below SN-limit(stable NS)
      companion[0].r[n] = NSradius(companion[0].cm[n]);	//remnant radius
      system.stagechange = max(system.stagechange,nextstage(companion[0],n));
    }else if (companion[0].stage[n]==5){	//stable BH
      companion[0].r[n] = Schwarzschildradius(companion[0].cm[n]);	//remnant radius
      system.stagechange = max(system.stagechange,nextstage(companion[0],n));
    }else{	//companion explodes
      companion[0].stage[n] = -3;
    }
  }
//  cout << endl << "a_ini-a_fin=" << system.a[n-1]-system.a[n] << " a_fin/a_ini=" << system.a[n]/system.a[n-1] << endl;
  //update times
  star[0].t[n] += tCE;
  companion[0].t[n] += tCE;
  if (star[0].t[n]>star[0].track.t[star[0].track.n-1]){
    cerr << "#Warning: star[0].t[" << n << "]=" << star[0].t[n] << "yr set to " << star[0].track.t[star[0].track.n-1] << "yr" << endl;
    star[0].t[n] = star[0].track.t[star[0].track.n-1];	//maximal time for donor
  }
  if (companion[0].t[n]>companion[0].track.t[companion[0].track.n-1]){
    cerr << "#Warning: companion[0].t[" << n << "]=" << companion[0].t[n] << "yr set to " << companion[0].track.t[companion[0].track.n-1] << "yr" << endl;
    companion[0].t[n] = companion[0].track.t[companion[0].track.n-1];	//maximal time for accretor
  }
  while (((companion[0].t[n]>companion[0].track.t[companion[0].track.n-1])||((companion[0].stage[n]==0)&&(companion[0].t[n]>companion[0].track.t[companion[0].track.TAMS])))&&(companion[0].stage[n]!=-3)){
    system.stagechange=max(system.stagechange,nextstage(companion[0],n));	//check if companion gets older than its current evolutionary track (consider particularly the TAMS)
  }
  if (companion[0].stage[n]!=-3){
    //determine/interpolate values for the companion
    for (j=companion[0].last;j<companion[0].track.n-1;j++){	//find time index of companion[0].t[n]
      if (companion[0].track.t[j]>companion[0].t[n]) break;
    }
    tratio = ratio(companion[0].t[n],companion[0].track.t[j],companion[0].track.t[j-1]);	//calculate ratio in time
//    if (screen) cout << "companion[0].t[" << n << "]=" << companion[0].t[n] << "yr companion[0].track.t[" << j << "]=" << companion[0].track.t[j] << "yr companion[0].track.t[" << j-1 << "]=" << companion[0].track.t[j-1] << "yr diff=" << companion[0].track.t[j]-companion[0].track.t[j-1] << "yr companion[0].track.n-1=" << companion[0].track.n-1 << " t_ratio=" << tratio << endl;
    if ((tratio<0)||(tratio>1)){
      cerr << endl << "#Error: time out of track: companion[0].t[" << n << "]=" << companion[0].t[n] << "yr companion[0].track.t[" << j << "]=" << companion[0].track.t[j] << "yr companion[0].track.t[" << j-1 << "]=" << companion[0].track.t[j-1] << "yr diff=" << companion[0].track.t[j]-companion[0].track.t[j-1] << "yr companion[0].track.n-1=" << companion[0].track.n-1 << " t_ratio=" << tratio << endl;	//write error
      screen = true;
      if (tratio<0) tratio = 0.0; else tratio = 1.0;
    }
    companion[0].m[n] = companion[0].track.m[j]-tratio*(companion[0].track.m[j]-companion[0].track.m[j-1]);	//interpolate mass
    companion[0].r[n] = companion[0].track.r[j]-tratio*(companion[0].track.r[j]-companion[0].track.r[j-1]);	//interpolate radius
    companion[0].cm[n] = companion[0].track.cm[j]-tratio*(companion[0].track.cm[j]-companion[0].track.cm[j-1]);	//interpolate core mass
    companion[0].llum[n] = companion[0].track.llum[j]-tratio*(companion[0].track.llum[j]-companion[0].track.llum[j-1]);	//interpolate luminosity
    companion[0].lteff[n] = companion[0].track.lteff[j]-tratio*(companion[0].track.lteff[j]-companion[0].track.lteff[j-1]);	//interpolate effective temperature
//    companion[0].lambda[n] = companion[0].track.lambda[j]-tratio*(companion[0].track.lambda[j]-companion[0].track.lambda[j-1]);	//interpolate lambda
    EbindperG = companion[0].track.m[j]*(companion[0].track.m[j]-companion[0].track.cm[j])/(companion[0].track.lambda[j]*companion[0].track.r[j])
                -tratio*(companion[0].track.m[j]*(companion[0].track.m[j]-companion[0].track.cm[j])/(companion[0].track.lambda[j]*companion[0].track.r[j])
                        -companion[0].track.m[j-1]*(companion[0].track.m[j-1]-companion[0].track.cm[j-1])/(companion[0].track.lambda[j-1]*companion[0].track.r[j-1]));	//interpolate current binding energy of the envelope
    if ((companion[0].m[n]-companion[0].cm[n]==0.0)&&(EbindperG==0.0)) companion[0].lambda[n] = -1.0e+99;	//no interpolation
    else companion[0].lambda[n] = companion[0].m[n]*(companion[0].m[n]-companion[0].cm[n])/(EbindperG*companion[0].r[n]);	//get lambda from interpolated Ebind	// E_bind=G*M*M_env / lambda*R
    if(isnan(companion[0].lambda[n])){ screen = true; cerr << endl << "#companion[0].lambda[n]=" << companion[0].lambda[n];}
    companion[0].last = j-1;
    if (companion[0].last>=companion[0].rmax){	//update rmax
      companion[0].rmax = min(companion[0].last+1,companion[0].track.n-1);
      if (debug) cerr << endl << "companion[0].rmax=" << companion[0].rmax;
      for (j=companion[0].last+1;j<companion[0].track.n-1;j++){
        if (companion[0].track.r[j]>companion[0].track.r[companion[0].rmax]){
          if ((j<=companion[0].track.TAMS)||(companion[0].last>=companion[0].track.TAMS)) companion[0].rmax = j;
        }
      }
      if (debug) cerr << " --> companion[0].rmax=" << companion[0].rmax << endl;
    }
    if (companion[0].r[n]<0) cerr << "#Error: companion[0].r[" << n << "]=" << companion[0].r[n] << "Rsun companion[0].track.r[" << j << "]=" << companion[0].track.r[j] << "Rsun companion[0].track.r[" << j-1 << "]=" << companion[0].track.r[j-1] << "Rsun t_ratio=" << tratio << endl;	//write error
  }
  star[0].m[n] = fmax(star[0].cm[n],1.0e-10);	//star mass reduced to core mass, but minimal to 10^(-10)Msun
//  star[0].r[n] = 0.212*pow(star[0].m[n],0.654);	// Rasmus fit to Onno's ZAMS He-stars. 
  //see also T.Wettig, G.E. Brown/new astronomy 1 (1996) 17-34, (A.17)
  if (star[0].stage[n]<2){	//star become naked helium star
    star[0].stage[n] = -2;
    system.stagechange = max(system.stagechange,nextstage(star[0],n));
    star[0].r[n] = star[0].track.r[star[0].last];
    star[0].cm[n] = star[0].track.cm[star[0].last];
    star[0].llum[n] = star[0].track.llum[star[0].last];
    star[0].lteff[n] = star[0].track.lteff[star[0].last];
    star[0].lambda[n] = star[0].track.lambda[star[0].last];
    if ((rocherad(system.a[n],star[0].m[n]/companion[0].m[n])<star[0].r[n])||(star[0].m[n]<0.1)){	// merging (via instant case A RLO)
//    if ((rocherad(system.a[n],star[0].m[n]/companion[0].m[n])<pow(10.0,-0.820+0.911*log10(star[0].m[n])-0.156*log10(star[0].m[n])*log10(star[0].m[n])))||(star[0].m[n]<0.1)){	// merging (via instant case A RLO) use R_{He}^{ZAMS} from Crowther(2007): lg(R/Rsun)=-0.820+0.911*lg(M/Msun)-0.156*lg(M/Msun)**2
      return 4;
    }else if (rocherad(system.a[n],companion[0].m[n]/star[0].m[n])<companion[0].r[n]){	// merging (via instant case A RLO from secondary! star)
      return 4;
    }else{
      return 0;
    }
  }else if (star[0].stage[n]==2){	//star will explode
    star[0].stage[n] = -3;
    star[0].r[n] = WDradius(star[0].m[n]);
    if ((rocherad(system.a[n],star[0].m[n]/companion[0].m[n])<star[0].r[n])||(star[0].m[n]<0.1)){	// merging (via instant case A RLO)
      return 4;
    }else if (rocherad(system.a[n],companion[0].m[n]/star[0].m[n])<companion[0].r[n]){	// merging (via instant case A RLO from secondary! star)
      return 4;
    }else{
      if (system.phase[n-1]==13) return 15;	//primary explodes
      else if (system.phase[n-1]==23) return 25;	//secondary explodes
      else return -1;	//undefined
    }
  }else{
    return 99;
  }
}

int merge(t_system& system){
  int n=system.n-1;	//last entry index

  system.a[n] = -1.0e+99;
  return -1;
}

int supernova(t_system& system){

/*
  This subroutine calculates the effects on the binary system of a
  supernova explosion of one of the stars. It is assumed that all CO-cores 
  more massive than 1.7 M_sun (roughly corresponding to a Helium star mass
  of 2.8 M_sun, depending on P_orb) will explode and leave a neutron star.
  See Dewi et al. (2002) and Pols (2002) for further details.
  CO-cores more massive than 4.0 M_sun are assumed to leave black holes.
  WD masses are given by:
  For binary ZAMS stars, this corresponds to a mass of 10 M_sun (8-9 caseC 
  10-12 caseB) and 8 M_sun for very wide binaries (and single stars)
  -see Bhattacharya & van den Heuvel p.39, and for further details see
  van den Heuvel, Accretion-driven stellar X-ray sources.
  It is assumed that the collapse of the CO-core takes
  place instantaneously compared with the binary's orbital period, and
  that the mass of a neutron star is equal to 1.3 M_sun. We also assume
  that the mass of the companion star is not affected by the impact of the
  supernova shell - Tauris & Takens (1998) and Wheeler, Lecar & McKee (1975).
  The pre-supernova orbit is assumed to be circular. 
  If there are no kicks then systems which eject at least half their total
  mass (mu<0.5) become unbound, while those that eject less that half their
  total mass (mu>0.5) remain bound (virial theorem!), Dewey & Cordes p.783.
  Including kicks we use the formulaes of J.G.Hills (1983) for bound systems.
*/

  int n=system.n-1;	//last entry index
  t_star* star;	//star
  t_star* companion;	//companion
  double v,vcore;	//velocity in km/s	//,factorCO=0.0
  double BH_limit=30.0,NS_limit=NSmax,CO_limit=0.0;	//mass limits to form a black hole(BH), a neutron star(NS) or a carbon-oxygen white dwarf(CO-WD)
  double NS_CC_limit=NSmax;	//mass limits to get a core collape SN(CC)	//NS_EC_limit,	//an electron capture SN(EC) or 
  double PI_limit_low=140.0,PI_limit_up=260.0;	//mass limits for a pair-instability "SN"
  double theta=0.0,phi=0.0,w=0.0;	//velocity kick parameter: theta & phi: angles, w: maxwell velocity
  double phase=0.0;	//phase angle of an eccentric system
  double wx,wy,wz;	//velocity kick parameter: x,y,z component of maxwell velocity
  double vCM_x,vCM_y,vCM_z,vCM;
//  double vCM_90,angle;
//  double delta_M,bracket21=1.0,theta_crit,vfactor;
//  double Mreduced,L_orb,E_orb;
//  double R2;
  double M_shell=0.0;	//mass of ejected shell
  double r_ini,r_comp;	//initial separation of the system, initial radius of the companion star
  double E_SN,v_eject,v_esc,psi,x_crit,FF,eta_x,vim;
  double mm,xi,PP,QQ,RR,SS,v15,v15x,v15y,v15z,v25,v25x,v25y,v25z,tan_gamma;
  double McoreFeSi, McoreCO, MlayerHe=0.0;	//FeSi core mass, CO core mass, He layer mass
  double fHe=0.80;	//fraction of the contributing He layer	//put to 80% because of wind loss during late burning
  double Porbini;	//iron core mass after case BB RLO; helium envelope mass after case BB RLO; orbital period after last event before case BB RLO	//Mcore,Menvelope,
  int oldstage=0;
  int oldphase=0;
  int whichstar=0;
  bool getWD=true;
  long i;

//  cout  << "system.prim.stage[" << n-1 << "]=" << system.prim.stage[n-1] << " system.sec.stage[" << n-1 << "]=" << system.sec.stage[n-1] << " system.prim.stage[" << n << "]=" << system.prim.stage[n] << " system.sec.stage[" << n << "]=" << system.sec.stage[n] << endl;
  if (debug) cerr << "n=" << n << endl;
  if ((system.prim.stage[n]==-3)&&(system.phase[n]!=25)){	//primary explodes
    star = &system.prim;
    companion = &system.sec;
    whichstar = 1;
  }else if (system.sec.stage[n]==-3){	//secondary explodes
    star = &system.sec;
    companion = &system.prim;
    whichstar = 2;
  }else{
    cerr << "#Error: no star for supernova/planetary nebula selected: system.prim.stage[" << n << "]=" << system.prim.stage[n] << " system.sec.stage[" << n << "]=" << system.sec.stage[n] << endl;	//write error
    screen = true;	//enable screen output
    return -1;
  }
  r_comp = companion[0].r[n];	//set companion radius
  if (n>=2){	//get stage and phase before the SN
    oldphase = system.phase[n-2];
    oldstage = star[0].stage[n-1];
    if ((oldstage==-3)&&(oldphase==whichstar+((whichstar%2)+1)*10)){	//if accretor explodes
      if (companion[0].stage[n-2]<2) oldstage = 1;	//hydrogen accreted -> hydrogen star explodes
	//otherwise star type did not change during mass transfer
    }
    if ((oldstage>2)||(oldstage<0)) oldstage = star[0].stage[n-2];	//get previous stage
    Porbini = 2.0*M_PI*sqrt(pow(star[0].r[n-2]/rocherad(1.0,star[0].m[n-2]/1.35)*(star[0].m[n-2]+1.35)/(star[0].inimass+1.35),3)/(G*(star[0].inimass+1.35)))/day;	//initial orbital period in days: get seperation if companion is a 1.35Msun NS, undo wind mass loss to initial mass, get initial period
    if ((oldstage==-3)&&(n>=3)){	//get stage and phase before previous SN
      oldphase = system.phase[n-3];
      oldstage = star[0].stage[n-3];
    }
    if (debug){ cerr << "oldstage=" << oldstage;
cerr << " oldphase=" << oldphase << endl;}
//    if (((system.phase[n-2]==21)&&(system.sec.stage[n]==-3))||((system.phase[n-2]==12)&&(system.prim.stage[n]==-3))) r_comp = companion[0].r[n-2];	//needed if much mass is added to core during transfer
  }else{
    Porbini = system.P[0]/day;	//initial orbital period in days
  }
  if (debug) cerr << "Porbini=" << Porbini << "day" << endl;
  // set limits to decide between WD, NS and BH
  if ((!separateevolution)&&(oldstage>2)){	//a remnat acceted mass
    McoreCO = star[0].m[n];	//save mass of the CO core
    if (oldstage==3){	//was a white dwarf(WD)
      if (star[0].m[n]>WDmax){	//star exceeds Chandrasekhar limit(maximal mass of a WD)
        return -1;
        getWD = false;
        star[0].stage[n] = 5;	//set stage to black hole
        star[0].m[n] *= 0.5;
        star[0].cm[n] = star[0].m[n];
        w = getkickvelocity(200);	//get kick velocity (uniform [0:200[)
        getrandomdirection(theta, phi);	//get kick direction
        star[0].r[n] = Schwarzschildradius(star[0].m[n]);	//get remnant radius in Rsun
        M_shell = star[0].m[n-1]-star[0].m[n];	//get mass of the ejected shell in Msun
      }else star[0].stage[n] = oldstage;	//WD stays as a WD
    }
    if (oldstage==4){	//was a neutron star(NS)
//      if (star[0].m[n]>NSmax) BH_limit = star[0].inimass;	//star exceeds maximal mass of a NS
      if (star[0].m[n]>NSmax) McoreCO = -star[0].m[n];	//star exceeds maximal mass of a NS
      else star[0].stage[n] = oldstage;	//NS stays as a NS
    }
    if (oldstage==5) star[0].stage[n] = oldstage;	//black hole(BH) stays as a BH
    if (star[0].stage[n]!=-3){	//no supernova
      if (debug) cerr << "no supernova: star[0].stage[" << n << "]=" << star[0].stage[n] << endl;
      if (companion[0].stage[n]>2) return 99;	//both stars are at their final stage
      return 0;
    }
  }else if(oldstage==2){	//initial He-star mass
    PI_limit_up=133.0;	//upper mass limit for pair-instability (Heger & Woosley 2002, ApJ)
    PI_limit_low=64.0;	//lower mass limit for pair-instability (Heger & Woosley 2002, ApJ)
    BH_limit = 14.0;	//mass limit to get a black hole in Msun
    NS_limit = 2.55;	//mass limit to get a neutron star in Msun
    NS_CC_limit = 2.75;	//mass limit to get a core collapse SN in Msun
    CO_limit = 0.0;	//mass limit to get a carbon-oxygen white dwarf in Msun
    McoreCO = star[0].cm[n];	//save mass of the CO core
    MlayerHe = star[0].m[n]-McoreCO;	//get the He mass in Msun
    if (oldphase!=0){	//SN after a mass transfer: Tauris,Langer,Podsiadlowski 2015 eq.(8)
      if (screen && (output=='M')) cerr << "Porbini=" << Porbini << "day ";
      NS_limit = 1.0/(61.0*Porbini)+2.58;	//threshold according to eq.(10)
      if ((star[0].inimass>=NS_limit)&&(McoreCO<1.37)) McoreCO = 1.37;
      NS_CC_limit = 2.74*pow(Porbini,-0.00209/Porbini);	//threshold according to eq.(9)
      if ((star[0].inimass>=NS_CC_limit)&&(McoreCO<1.435)) McoreCO = 1.435;
    }
//    NS_EC_limit = NS_limit;	//mass limit to get an electron capture SN in Msun
  }else{	//initial H-star mass
    PI_limit_up=260.0;	//upper mass limit for pair-instability (Heger & Woosley 2002, ApJ)
    PI_limit_low=140.0;	//lower mass limit for pair-instability (Heger & Woosley 2002, ApJ)
    BH_limit = 30.0;	//mass limit to get a black hole in Msun
    NS_limit = 8.3;	//mass limit to get a neutron star in Msun
    NS_CC_limit = 9.0;	//mass limit to get a core collapse SN in Msun
    CO_limit = 0.0;	//mass limit to get a carbon-oxygen white dwarf in Msun
    if (star[0].m[n]==star[0].cm[n]){
      MlayerHe = star[0].track.cm[0];
      McoreCO = fmin(star[0].cm[n],getCOcore(MlayerHe));	//get mass of the CO core depending the maximum mass of the He core in Msun
      MlayerHe = MlayerHe-McoreCO;	//get mass of the He layer in Msun
    }else{
      MlayerHe = star[0].cm[n];
      McoreCO = getCOcore(MlayerHe);	//get mass of the CO core depending on the mass of the He core in Msun
      MlayerHe = MlayerHe-McoreCO;	//get mass of the He layer in Msun
    }
  }
  //save supernova data
  star[0].SN.startype = oldstage;
  star[0].SN.inimass = star[0].inimass;
  star[0].SN.Hecoremass = MlayerHe+McoreCO;
  star[0].SN.COcoremass = McoreCO;
  if (debug||screen) cerr << "BH_limit=" << BH_limit << "Msun NS_limit=" << NS_limit << "Msun NS_CC_limit=" << NS_CC_limit << "Msun CO_limit=" << CO_limit << "Msun star[0].inimass=" << star[0].inimass << "Msun McoreCO=" << McoreCO << "Msun MlayerHe=" << MlayerHe << "Msun" << endl;
//  if (star[0].inimass>=BH_limit){	//Black hole(BH) formation
  if ((McoreCO>=6.5)||(McoreCO<0)){	//Black hole(BH) formation
#pragma omp critical
{
    cBH++;	//increase counter of black holes
}	//end critical
    if (McoreCO<0){
      //screen = true;
      McoreCO = -McoreCO;
    }
    star[0].stage[n] = 5;	//set stage to black hole
    if (oldstage!=star[0].stage[n]) system.formation = 10*system.formation+star[0].stage[n];	//write formation channel
    star[0].m[n] = 0.8*(McoreCO+fHe*MlayerHe);	//get remnant mass in Msun
    if ((star[0].inimass>PI_limit_low)&&(star[0].inimass<PI_limit_up)) star[0].m[n] = 1e-99;	//pair-instability SN getting no remnant (Heger & Woosley 2002, ApJ)
    w = getkickvelocity(200.0);	//get kick velocity (uniform [0:200[)
    getrandomdirection(theta, phi);	//get kick direction
    star[0].r[n] = Schwarzschildradius(star[0].m[n]);	//get remnant radius in Rsun
    M_shell = star[0].m[n-1]-star[0].m[n];	//get mass of the ejected shell in Msun
    star[0].SN.remnantmass = star[0].m[n];	//save post-SN remnant mass
    star[0].SN.type = 4;	//set supernova type to collape to blackhole
//  }else if (star[0].inimass>=NS_limit){	//Neutron star(NS) formation
  }else if (McoreCO>=1.37){	//Neutron star(NS) formation
#pragma omp critical
{
    cNS++;	//increase counter of neutron stars
}	//end critical
    star[0].stage[n] = 4;	//set stage to neutron star
    if (oldstage!=star[0].stage[n]) system.formation = 10*system.formation+star[0].stage[n];	//write formation channel
//    if (star[0].inimass>=NS_CC_limit){	//iron core collapse SN
    if (McoreCO>=1.435){	//iron core collapse SN
      if ((oldstage==2)&&(oldphase==whichstar*10+(whichstar%2)+1)){	//stripped He-star with stable mass transfer(Tauris,Langer,Podsiadlowski 2015)	//&&(companion[0].stage[n]>=4)&&(companion[0].stage[n]<=5)&&(star[0].inimass<=3.5)
        McoreFeSi = 1.05667*pow(star[0].cm[n],0.453707);	//fit to Timmes et al. 1996
//        McoreFeSi = 1.3*pow(star[0].cm[n],0.35);	//test1 above Timmes et al. 1996
//        McoreFeSi = pow(star[0].cm[n]-1.4,0.3)+0.9;	//test2 above Timmes et al. 1996
        star[0].m[n] = (-1.0+sqrt(1.0+0.336*McoreFeSi))/0.168;	//solution of M_NS=M_core-M_NS^bind from Tauris,Langer,Podsiadlowski 2015 where M_NS^bind=0.084*M_NS^2/Msun from Lattimer & Yahil 1989 eq.(10) by converting energy to mass in Msun
        if (star[0].m[n]>0.99*NSmax){
//          screen = true;
          star[0].m[n] = 0.99*NSmax;
        }
//        Mcore = (1/(400.0*Porbini)+0.49)*star[0].inimass-(0.016/Porbini-0.106);	//Tauris,Langer,Podsiadlowski 2015 eq.(3)
//        star[0].m[n] = fmin((-1.0+sqrt(1.0+0.336*Mcore))/0.168,0.99*NSmax);	//solution of M_NS=M_core-M_NS^bind from Tauris,Langer,Podsiadlowski 2015 where M_NS^bind=0.084*M_NS^2/Msun from Lattimer & Yahil 1989 eq.(10) by converting energy to mass in Msun
//        Menvelope = star[0].m[n-1]-Mcore-0.03;	//correct envelope mass if needed
//        system.M[n-1] = star[0].m[n-1]+companion[0].m[n-1];	//correct pre SN mass in Msun
//        if ((Porbini<3.0)&&(companion[0].stage[n]>3)) screen = true;
//        if (star[0].m[n-1]<star[0].m[n]){
//          if (debug) cerr << "correct star[0].m[n-1]=" << star[0].m[n-1] << "Msun to " << star[0].m[n]+1.0e-10 << "Msun" << endl;
//          star[0].m[n-1] = star[0].m[n]+1.0e-10;
//          system.M[n-1] = system.prim.m[n-1]+system.sec.m[n-1];
//        }
      }else{
        star[0].m[n] = fmin(0.23*McoreCO+0.83,star[0].m[n-1]);	//get remnant mass in Msun
      }
      if (oldstage==2){
        if (ranv(5)<0.8){
          if ((oldphase==whichstar*10+(whichstar%2)+1)||(oldphase==whichstar*10+3)){
            if ((companion[0].stage[n]>=4)&&(companion[0].stage[n]<=5)){
              w = getkickvelocity(-30.0*sqrt(3.0));	//get kick velocity after ultra-stripping (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 30km/s [try 60km/s])
            }else w = getkickvelocity(-60.0*sqrt(3.0));	//get kick velocity after He-stripping (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 60km/s)
          }else w = getkickvelocity(-120.0*sqrt(3.0));	//get kick velocity for He-star (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 120km/s)
        }else w = getkickvelocity(-200.0*sqrt(3.0));	//20% of high kick (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 200km/s)
      }else w = getkickvelocity(-265.0*sqrt(3.0));	//get kick velocity of isolated star (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 265km/s, cf. Hobbs et al. 2005)
//      }else w = getkickvelocity(-265);	//get kick velocity of isolated star (maxwell: rms of 265km/s, cf. Hobbs et al. 2005)
      star[0].SN.remnantmass = star[0].m[n];	//save post-SN remnant mass
      star[0].SN.type = 3;	//set supernova type to iron core collape
    }else{	//electron capture SN
      star[0].m[n] = fmin(1.24+0.8*(McoreCO-1.37),star[0].m[n-1]);	//get remnant mass in Msun
      w = getkickvelocity(50.0);	//get kick velocity (uniform [0:50[)
/*      if (oldstage==2){			//try large kicks like CC
        if (ranv(5)<0.8){
          if ((oldphase==whichstar*10+(whichstar%2)+1)||(oldphase==whichstar*10+3)){
            if ((companion[0].stage[n]>=4)&&(companion[0].stage[n]<=5)){
              w = getkickvelocity(-30.0*sqrt(3.0));	//get kick velocity after ultra-stripping (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 30km/s [try 60km/s])
            }else w = getkickvelocity(-60.0*sqrt(3.0));	//get kick velocity after He-stripping (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 60km/s)
          }else w = getkickvelocity(-120.0*sqrt(3.0));	//get kick velocity for He-star (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 120km/s)
        }else w = getkickvelocity(-200.0*sqrt(3.0));	//20% of high kick (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 200km/s)
      }else w = getkickvelocity(-265.0*sqrt(3.0));	//get kick velocity of isolated star (maxwell: projected 1D-rms[3D-rms/sqrt(3)] of 265km/s, cf. Hobbs et al. 2005)*/
      star[0].SN.remnantmass = star[0].m[n];	//save post-SN remnant mass
      star[0].SN.type = 2;	//set supernova type to electron capture
    }
    getrandomdirection(theta, phi);	//get kick direction
    star[0].r[n] = NSradius(star[0].m[n]);	//get remnant radius in Rsun
    M_shell = star[0].m[n-1]-star[0].m[n];	//get mass of the ejected shell in Msun
    if (distfile2!=NULL){	//write specified system to distfile2 output
      fprintf(distfile2,"%14.5f\t%13.5f\t%13.5f\t%3d  \t%13.5f\t%13.5f\t%2d  \t%9.3f\t%10.6f  \t%10.6f\t%9.5f\t%3d  \t", star[0].m[n-1], star[0].SN.Hecoremass, star[0].SN.COcoremass, star[0].SN.startype, star[0].SN.inimass, star[0].SN.remnantmass, star[0].SN.type, w, theta*180.0/M_PI, phi*180.0/M_PI, companion[0].m[n], companion[0].stage[n]);
      for (i=0;i<n;i++) fprintf(distfile2,"%d",system.phase[i]);
      fprintf(distfile2,"\n");
      fflush(distfile2);	//update file
    }
  }else if(getWD){	//White dwarf(WD) formation
#pragma omp critical
{
    cWD++;	//increase counter of white dwarfs
}	//end critical
    system.phase[n-1]++;	//set current phase to planetary nebula
    star[0].stage[n] = 3;	//set stage to white dwarf
    if (oldstage!=star[0].stage[n]) system.formation = 10*system.formation+star[0].stage[n];	//write formation channel
/*    if ((oldstage==2)||(star[0].inimass>=CO_limit)) star[0].m[n] = McoreCO;	//get remnant mass in Msun of CO white dwarf
    else star[0].m[n] = star[0].cm[n];	//get remnant mass in Msun of He white dwarf*/
    star[0].m[n] = fmin(fmax(McoreCO,0.6),star[0].m[n]);	//get remnant mass in Msun of white dwarf being a CO-WD with m>=0.6
    if (star[0].m[n]>WDmax){
      if (!debug) cout << "BH_limit=" << BH_limit << "Msun NS_limit=" << NS_limit << "Msun NS_CC_limit=" << NS_CC_limit << "Msun CO_limit=" << CO_limit << "Msun star[0].inimass=" << star[0].inimass << "Msun McoreCO=" << McoreCO << "Msun" << endl;
      cerr << "#Error: white dwarf mass=" << star[0].m[n] << "Msun McoreCO=" << McoreCO << "Msun star[0].cm[n]=" << star[0].cm[n] << "Msun" << endl;
      screen = true;
    }
    star[0].r[n] = WDradius(star[0].m[n]);	//get remnant radius in Rsun
    system.M[n] = star[0].m[n]+companion[0].m[n];
    system.a[n] = system.a[n-1]*(system.M[n-1])/(system.M[n]);
    system.peri[n] = system.a[n]*(1.0-system.e[n]);
    //do not use values any longer
    star[0].cm[n] = star[0].m[n];
    star[0].llum[n] = -1.0e+99;
    star[0].lteff[n] = -1.0e+99;
    star[0].lambda[n] = -1.0e+99;
    system.stagechange = max(system.stagechange,nextstage(star[0],n));
#pragma omp critical
{
    cplanetarynebula++;	//increase counter of planetary nebulae
}	//end critical
    star[0].SN.remnantmass = star[0].m[n];	//save post-SN remnant mass
    star[0].SN.type = 1;	//set supernova type to planetary nebula
    if (rocherad(system.a[n],companion[0].m[n]/star[0].m[n])<r_comp) return 4;	//system merges
    if (companion[0].stage[n]>2) return 99;	//both stars are at their final stage
    return 0;
  }
//  w *= 0.5;
/*  if ((oldstage==2)&&(oldphase!=0)) star[0].cm[n] = star[0].m[n];
  // set limits to decide between WD, NS and BH
  if ((!separateevolution)&&(oldstage>2)){	//a remnat acceted mass
    BH_limit = NSmax;
    NS_limit = WDmax;	//Chandrasekhar limit
    CO_limit = 0.0;
    if ((oldstage==3)&&(star[0].cm[n]<NS_limit)) star[0].stage[n] = oldstage;	//WD stays as a WD
    if ((oldstage==4)&&(star[0].cm[n]<=BH_limit)) star[0].stage[n] = oldstage;	//NS stays as a NS
    if (oldstage==5) star[0].stage[n] = oldstage;	//BH stays as a BH
    if (star[0].stage[n]!=-3) return 0;	//no supernova
  }else if (star[0].HeRLO){	//Roche-lobe overflow from Helium star happend
    BH_limit = 4.0;
    NS_limit = 1.7;
    CO_limit = 1.1;
  }else if(star[0].RLO){	//Roche-lobe overflow from normal star happend
    BH_limit = 5.5;
    NS_limit = 2.8;
    CO_limit = 2.0;
  }else{	//no Roche-lobe overflow
    BH_limit = 20.0;  // ZAMS of roughly 25 M_sun
    NS_limit = 10.0;  // ZAMS of roughly 11 M_sun
    CO_limit = 0.0;   
    if (star[0].inimass>=24.0) BH_limit = 0.99*star[0].cm[n];
    if (star[0].inimass>=10.0) NS_limit = 0.99*star[0].cm[n]; 
  }
  if ((oldstage==2)&&(oldphase!=0)){	//Tauris,Langer,Podsiadlowski 2015
    NS_limit = 1.0/61.0*(day/system.P[n])+2.58;	//threshold according to eq.(8) but uses mass at onset of the supernova
  }
  // apply supernova depending on the kind of remnant
  if (star[0].cm[n]>BH_limit){	//Black hole(BH)
    star[0].stage[n] = 5;
    if (oldstage!=star[0].stage[n]) system.formation = 10*system.formation+star[0].stage[n];
//    delta_M = 0.3*star[0].cm[n];	//lost mass: APPROX. for small mass BH
    M_shell = 0.16*star[0].cm[n];	//get mass of the shell with assuming a BH E_bind of 20%: M_shell=M1_core-(M1_remnant*1.20)=M1_core-(0.7*M1_core*1.20)=0.16*M1_core
    star[0].m[n] = 0.7*star[0].cm[n];	//remnant mass
    star[0].r[n] = Schwarzschildradius(star[0].m[n]);	//remnant radius
    //get random direction on unit sphere
    if (kickv!=0.0){	//if supernova should be asymmetric
      w = getkickvelocity(2);	//get kick velocity
      getrandomdirection(theta, phi);	//get kick direction
      if (kickv<0.0) w *= 1.3/star[0].m[n];
    }
//    w = 0.0;	//no kick
  }else if (star[0].cm[n]>=NS_limit){	//Neutron star(NS)
    star[0].stage[n] = 4;
    if (oldstage!=star[0].stage[n]) system.formation = 10*system.formation+star[0].stage[n];
    //remnant mass
    if ((oldstage==2)&&(oldphase!=0)&&(star[0].track.t[star[0].track.n-1]-star[0].t[n]<1.0)){	//Tauris,Langer,Podsiadlowski 2015	//&&(companion[0].stage[n]==4)
      if (star[0].m[n]>2.74*pow(system.P[n]/day,-0.00209*(day/system.P[n]))){	//threshold according to eq.(7) but uses mass at onset of the supernova	//FeCCSN
        Mcore = (day/(400.0*system.P[n])+0.49)*star[0].inimass-(0.016*day/system.P[n]-0.106);	//Tauris,Langer,Podsiadlowski 2015 eq.(3)
        star[0].m[n] = fmin((-1.0+sqrt(1.0+0.336*Mcore))/0.168,0.9*NSmax);	//solution of M_NS=M_core-M_NS^bind from Tauris,Langer,Podsiadlowski 2015 where M_NS^bind=0.084*M_NS^2/Msun from Lattimer & Yahil 1989 eq.(10) by converting energy to mass in Msun
      }else star[0].m[n] = 1.25;	//EC SN
    }else star[0].m[n] = 1.3;
    star[0].r[n] = NSradius(star[0].m[n]);	//remnant radius
// NOTE!!!! @Matthias: This is the assumed gravitational mass which is 
//          about 10% smaller that its rest mass which is used in the eqs. below.
//          A correction shold be incorporated sometime.
//    delta_M = star[0].cm[n]-star[0].m[n];	//lost mass
    M_shell = fmax(1.0e-99,star[0].cm[n]-(star[0].m[n]+0.084*pow(star[0].m[n],2))); // including E_bind from Lattimer & Yahil (1989)
    if (kickv!=0.0){	//if supernova should be asymmetric
      w = getkickvelocity(2);	//get kick velocity
      getrandomdirection(theta, phi);	//get kick direction
    }
  }else{	// White dwarf(WD)
    system.phase[n-1]++;
    if (oldstage==2){	//Tauris,Langer,Podsiadlowski 2015
      if (system.P[n]/day<=0.07) star[0].m[n] = 0.9;
      else star[0].m[n] = 1.2;
    }else if (star[0].cm[n]>CO_limit){	//CO-WD
      if (star[0].HeRLO){
        star[0].m[n] = 1.3-0.33*(1.7-star[0].cm[n]);
      }else if (star[0].RLO){
        star[0].m[n] = 1.3-0.29*(2.7-star[0].cm[n]);
      }else{
        star[0].m[n] = fmin(1.3-0.06*(NS_limit-star[0].cm[n]),star[0].cm[n]);
      }
    }else{
      if (star[0].HeRLO) star[0].m[n] = star[0].cm[n];
      else if (star[0].RLO){	//Test to avoid M_WD>1.1 and get M_WD~M_core for small M_core
        factorCO = fmin(CO_limit-1.1,2.2-CO_limit);
        star[0].m[n] = star[0].cm[n]+(factorCO-CO_limit+1.1)*pow(star[0].cm[n]/CO_limit,2)-factorCO*pow(star[0].cm[n]/CO_limit,3);
      }
//      else if (star[0].RLO) star[0].m[n] = 0.55*star[0].cm[n]; 
//      star[0].m[n] = star[0].cm[n];  // WD TESTING !!!!!!!  to avoid Mcore=0.35 -> 0.18 M WD
    }
    star[0].r[n] = WDradius(star[0].m[n]);	//remnant radius
    star[0].stage[n] = 3;
    if (oldstage!=star[0].stage[n]) system.formation = 10*system.formation+star[0].stage[n];
//    delta_M = fmax(star[0].cm[n]-star[0].m[n],0.0);	//determine lost core mass, minimal is no mass lost
//    M_shell=delta_M;
    system.a[n] = system.a[n-1]*(star[0].cm[n]+companion[0].m[n])/(star[0].m[n]+companion[0].m[n]);
    system.peri[n] = system.a[n]*(1.0-system.e[n]);
    //do not use values any longer
    star[0].cm[n] = star[0].m[n];
    star[0].llum[n] = -1.0e+99;
    star[0].lteff[n] = -1.0e+99;
    star[0].lambda[n] = -1.0e+99;
    system.stagechange = max(system.stagechange,nextstage(star[0],n));
    if (rocherad(system.a[n],companion[0].m[n]/star[0].m[n])<r_comp) return 4;	//system merges
    if (companion[0].stage[n]>2) return 99;	//both stars are at their final stage
    return 0;
  }*/
//  if (screen && (output=='M')) cout << "w=" << w/kms << "km/s theta=" << theta*180.0/M_PI << "degree phi=" << phi*180.0/M_PI << "degree";
  if (screen && (output=='M')) cout << "w=" << w << "km/s theta=" << theta*180.0/M_PI << "degree phi=" << phi*180.0/M_PI << "degree";
  if (screen && (output=='T')){
    cout << "pre SN -- time: " << system.prim.t[n-1] << endl;
    cout << "  Mp Ms a: " << system.prim.m[n-1] << " " << system.sec.m[n-1] << " " << system.a[n-1] << endl;
    cout << "  (w,theta,phi)=(" << w << ", " << theta*180.0/M_PI << ", " << phi*180.0/M_PI << ")" << endl;
  }
  //save supernova data
  star[0].SN.w = w;
  star[0].SN.theta = theta*180.0/M_PI;
  star[0].SN.phi = phi*180.0/M_PI;
#pragma omp critical
{
  csupernova++;	//increase counter of supernovae
}	//end critical
  //spherical to cartesian coordinates
  wx = w*cos(theta);
  wy = w*sin(theta)*cos(phi);
  wz = w*sin(theta)*sin(phi);
//  if (debug) cerr << endl << "wx=" << wx/kms << "km/s wy=" << wy/kms << "km/s wz=" << wz/kms << "km/s";
  if (debug) cerr << endl << "wx=" << wx << "km/s wy=" << wy << "km/s wz=" << wz << "km/s";
  // pre-SN rel. velocity: ref. frame fixed on comp. star
  if(system.e[n]<1.0E-5){	//circular orbit
//    v = sqrt(G*(star[0].cm[n]+companion[0].m[n])/system.a[n]);
//    v = sqrt(G*(star[0].cm[n]+companion[0].m[n])/system.a[n])/kms;
    v = sqrt(G*(star[0].m[n-1]+companion[0].m[n])/system.a[n])/kms;
    r_ini = system.a[n];
  }else{	//eccentric orbit at phase angle phase
    phase = ranphase(system.e[n]);
//    if (count_mp<1000) screen = true;	//enable screen output
//    v = sqrt((1.+pow(system.e[n],2)+2.*system.e[n]*cos(phase))*G*(star[0].cm[n]+companion[0].m[n])/((1.-pow(system.e[n],2))*system.a[n]));
//    v = sqrt((1.+pow(system.e[n],2)+2.*system.e[n]*cos(phase))*G*(star[0].cm[n]+companion[0].m[n])/((1.-pow(system.e[n],2))*system.a[n]))/kms;
    v = sqrt((1.+pow(system.e[n],2)+2.*system.e[n]*cos(phase))*G*(star[0].m[n-1]+companion[0].m[n])/((1.-pow(system.e[n],2))*system.a[n]))/kms;
    r_ini = (1.-pow(system.e[n],2))/(1.+system.e[n]*cos(phase))*system.a[n];
//    if (screen && (output=='M')) cout << " phase angle=" << phase*180.0/M_PI << "degree v=" << v/kms << "km/s r_ini=" << r_ini << "Rsun";
    if (screen && (output=='M')) cout << " phase angle=" << phase*180.0/M_PI << "degree v=" << v << "km/s r_ini=" << r_ini << "Rsun";
  }
  // pre-SN velocity of exploding core in c.m. ref. frame
//  vcore = v*(companion[0].m[n])/(star[0].cm[n]+companion[0].m[n]);
  vcore = v*(companion[0].m[n])/(star[0].m[n-1]+companion[0].m[n]);
//  if (debug) cerr << endl << "v=" << v/kms << "km/s r_ini=" << r_ini << "Rsun vcore=" << vcore/kms << "km/s";
  if (debug) cerr << endl << "v=" << v << "km/s r_ini=" << r_ini << "Rsun vcore=" << vcore << "km/s";
  // Impact on companion star of SN shell
  // The formula below is Tauris & Taken's eq.(5) + Sect.3.1.2
  // See also Wheeler, Lecar & McKee (1975) and Tauris+(2014) on HVS stars
  // The SN explosion energy is typically 1-10 for for SNe Ib/c. May depend on M1_core...?
//  E_SN = 1.23E+51*cgsEnergy;	// SN explosion energy (can be up to 10 foe), here 1.23 foe
  E_SN = 1.23E+51;	// SN explosion energy (can be up to 10 foe), here 1.23 foe
  M_shell = fmax(M_shell,1.0e-10);	//mass has to be positive
//  v_eject = sqrt(2.*E_SN/M_shell);	// velocity of the ejected shell: typically ~10^4 km/s. It can vary...; assuming all the SN explosion energy goes into the velocity of the shell
  v_eject = sqrt(2.*E_SN*cgsEnergy/M_shell)/kms;	// velocity of the ejected shell: typically ~10^4 km/s. It can vary...; assuming all the SN explosion energy goes into the velocity of the shell
  // The factor 5/3 below is because v_esc is not at the surface, see sec.3.1.2 Tauris & Takens(1998), eq.10 from Wheeler, Lecar & McKee(1975)
//  v_esc = sqrt(2.*G*(5./3.)*companion[0].m[n]/r_comp);	// escape velocity of striped material: typically ~ 800 km/s
  v_esc = sqrt(2.*G*(5./3.)*companion[0].m[n]/r_comp)/kms;	// escape velocity of striped material: typically ~ 800 km/s
  psi = max(pow(r_comp/(2.*r_ini),2)*(M_shell/companion[0].m[n])*(v_eject/v_esc-1.),0.0);	// eq.58 from Tauris & Takens(1998)
//  psi = pow(r_comp/(2.*r_ini),2)*(M_shell/companion[0].m[n])*(v_eject/v_esc-1.);	// eq.58 from Tauris & Takens(1998)
  // Note: psi=6.*(1./r_ini)**(2.) only approx for 1 M_sun star companion; eq.59 from Tauris & Takens(1998)
  if ((psi<0)&&(companion[0].stage[n]<3)){	//lowest psi is 0.0 check it!
    cerr << endl << "#Error: negative psi=" << psi;	//write error
    screen = true;	//enable screen output
//    debug = true;	//enable debug output
  }
//  if (debug) cerr << endl << "M_shell=" << M_shell << "Msun r_comp=" << r_comp << "Rsun v_eject=" << v_eject/kms << "km/s v_esc=" << v_esc/kms << "km/s psi=" << psi;
  if (debug) cerr << endl << "M_shell=" << M_shell << "Msun star[0].m[n-1]=" << star[0].m[n-1] << "Msun star[0].m[n]=" << star[0].m[n] << "Msun r_comp=" << r_comp << "Rsun v_eject=" << v_eject << "km/s v_esc=" << v_esc << "km/s psi=" << psi;
  // Fit to Table 1(Tauris & Takens(1998)/Wheeler, Lecar & McKee(1975)):
  x_crit = 1.09-0.654*pow(psi,0.159);
  FF = pow(fmin(0.432*log(psi+1.01),0.999),1.3);	//F_strip+F_ablate
  eta_x = 0.2;
  vim = eta_x*v_eject*pow((r_comp/(2.*r_ini)),2)*(M_shell/companion[0].m[n])*pow(x_crit,2)*(1.+log(2.*v_eject/v_esc))/(1.-FF);	//eq.5 from Tauris & Takens(1998)
  if (companion[0].stage[n]>2){	//companion is a compact object -> no shell impact
    x_crit = 1.;
    FF = 0.;
    vim = 0.;
//  }else{
//    cout << endl << "psi=" << psi << " x_crit=" << x_crit << " FF=" << FF << " vim=" << vim << "km/s";
  }
//  if (debug) cerr << endl << "x_crit=" << x_crit << " FF=" << FF << " vim=" << vim/kms << "km/s";
  if (debug) cerr << endl << "x_crit=" << x_crit << " FF=" << FF << " vim=" << vim << "km/s";
  companion[0].m[n] *= (1.-FF);	//new companion mass; FF*M2 is stripped and ablated
  if (companion[0].m[n]<companion[0].cm[n]){	//consider loss of core material
    if (debug) cerr << endl << "companion[0].m[n]=" << companion[0].m[n] << "Msun < companion[0].cm[n]=" << companion[0].cm[n] << "Msun";
    if ((companion[0].stage[n]>=0)&&(companion[0].stage[n]<2)){	//star become naked helium star
      companion[0].stage[n] = -2;
      system.stagechange = max(system.stagechange,nextstage(companion[0],n));
      companion[0].r[n] = companion[0].track.r[companion[0].last];
      companion[0].cm[n] = companion[0].track.cm[companion[0].last];
      companion[0].llum[n] = companion[0].track.llum[companion[0].last];
      companion[0].lteff[n] = companion[0].track.lteff[companion[0].last];
      companion[0].lambda[n] = companion[0].track.lambda[companion[0].last];
//      screen = true;
    }else{	//companion will explode too
      companion[0].cm[n] = fmin(companion[0].cm[n],companion[0].m[n]);	//reduce core mass
      companion[0].stage[n] = -3;
    }
  }else if (companion[0].stage[n]<3){
    if (debug) cerr << endl << "companion[0].m[n-1]=" << companion[0].m[n-1] << "Msun -> companion[0].m[n]=" << companion[0].m[n] << "Msun" << endl;
    updatetrack(companion[0],n);
    if (companion[0].m[n]!=companion[0].m[n-1]){
      companion[0].r[n] = companion[0].track.r[companion[0].last];
      companion[0].llum[n] = companion[0].track.llum[companion[0].last];
      companion[0].lteff[n] = companion[0].track.lteff[companion[0].last];
      companion[0].lambda[n] = companion[0].track.lambda[companion[0].last];
    }
  }
  mm = (star[0].m[n]+companion[0].m[n])/(star[0].m[n-1]+companion[0].m[n-1]);	//ratio of total finial and initial mass; eq.22 from Tauris & Takens(1998)
  xi = (pow(v,2)+pow(w,2)+pow(vim,2)+2.*w*(v*cos(theta)-vim*sin(theta)*cos(phi)))/(mm*pow(v,2));	//eq.23 with eq.16 from Tauris & Takens(1998)
  if (debug) cerr << endl << "mm=" << mm << " xi=" << xi;
  PP = 1.-2.*mm+pow(w/v,2)+pow(vim/v,2)+2.*(w/pow(v,2))*(v*cos(theta)-vim*sin(theta)*cos(phi));	//eq.44 from Tauris & Takens(1998)
//  PP = mm*(xi-2);
  QQ = 1.+(PP/mm)-pow(w*sin(theta)*cos(phi)-vim,2)/(mm*pow(v,2));	//eq.45 from Tauris & Takens(1998)
  RR = ((sqrt(PP)/(mm*v))*(w*sin(theta)*cos(phi)-vim)-(PP/mm)-1.)*(star[0].m[n]+companion[0].m[n])/companion[0].m[n];	//eq.46 from Tauris & Takens(1998)
  SS = (1.+(PP/mm)*(QQ+1.))*(star[0].m[n]+companion[0].m[n])/companion[0].m[n];	//eq.47 from Tauris & Takens(1998)
  if (debug) cerr << endl << "PP=" << PP << " QQ=" << QQ << " RR=" << RR << " SS=" << SS;

  // The formula below is taken from J.G. Hills (21) and are ONLY valid when NO shel impact is included!
//  if (w!=0) bracket21 = (1.0-2.0*(delta_M/(star[0].cm[n]+companion[0].m[n]))-pow((w/v),2))*(0.5*(v/w));	//*(v/v));
//  if (bracket21 < -1.) bracket21 = -1.0;
//  if (bracket21 > 1.) bracket21 = 1.0;
//  theta_crit = acos(bracket21);
//  if (debug) cerr << endl << "delta_M=" << delta_M << "Msun star[0].cm[n]=" << star[0].cm[n] << "Msun companion[0].m[n]=" << companion[0].m[n] << "Msun v=" << v/kms << "km/s bracket21=" << bracket21 << endl;
//  if (debug) cerr << endl << "delta_M=" << delta_M << "Msun star[0].cm[n]=" << star[0].cm[n] << "Msun companion[0].m[n]=" << companion[0].m[n] << "Msun v=" << v << "km/s bracket21=" << bracket21 << endl;
//  if (theta < theta_crit) return -1;	// The binary is disrupted in SN

  // Calculation of the periastron distance:
  system.a[n] = -(r_ini*mm/PP);	//eq.26 or part of eq.57 from Tauris & Takens(1998)
  system.e[n] = sqrt(fabs(SS*companion[0].m[n]/(star[0].m[n]+companion[0].m[n])));	//part of eq.57 from Tauris & Takens(1998)
  system.peri[n] = system.a[n]*(1.0-system.e[n]);	//eq.57 from Tauris & Takens(1998)
  if(xi>2.){	//system is disrupted
//  if(PP>0.){	//system is disrupted
    v15x = wx*((1./RR)+1.)+((1./RR)+companion[0].m[n-1]/(star[0].m[n-1]+companion[0].m[n-1]))*v;	//eq.51 from Tauris & Takens(1998)
    v15y = wy*(1.-1./SS)+vim/SS+QQ*sqrt(PP)*v/SS;	//eq.52 from Tauris & Takens(1998)
    v15z = wz*((1./RR)+1.);	//eq.53 from Tauris & Takens(1998)
//    cout << endl << "v15: (" << v15x/kms << "," << v15y/kms << "," << v15z/kms << ")km/s" << endl;
//    cout << endl << "v15: (" << v15x << "," << v15y << "," << v15z << ")km/s" << endl;
    v15 = sqrt(v15x*v15x+v15y*v15y+v15z*v15z);
    v25x = -wx/((companion[0].m[n]/star[0].m[n])*RR)-(1./((companion[0].m[n]/star[0].m[n])*RR)+star[0].m[n-1]/(star[0].m[n-1]+companion[0].m[n-1]))*v;	//eq.54 from Tauris & Takens(1998)
    v25y = wy/((companion[0].m[n]/star[0].m[n])*SS)+(1.-1./((companion[0].m[n]/star[0].m[n])*SS))*vim -QQ*sqrt(PP)*v/((companion[0].m[n]/star[0].m[n])*SS);	//eq.55 from Tauris & Takens(1998)
    v25z = -wz/((companion[0].m[n]/star[0].m[n])*RR);	//eq.56 from Tauris & Takens(1998)
//    cout << endl << "v25: (" << v25x/kms << "," << v25y/kms << "," << v25z/kms << ")km/s" << endl;
//    cout << endl << "v25: (" << v25x << "," << v25y << "," << v25z << ")km/s" << endl;
    v25 = sqrt(v25x*v25x+v25y*v25y+v25z*v25z);
    // check if (otherwise) disrupted systems merges on inbound leg:
    // eq.(17): if tan(gamma)>0 <=> gamma>0 inbound leg (Fig.2)
    tan_gamma = (w*sin(theta)*cos(phi)-vim)/sqrt(pow(v+w*cos(theta),2)+pow(w*sin(theta)*sin(phi),2));
    if (debug) cerr << endl << "tan(gamma)=" << tan_gamma;
    //do not use values any longer
//    system.a[n] = fabs(system.a[n]);
    star[0].cm[n] = star[0].m[n];
    star[0].llum[n] = -1.0e+99;
    star[0].lteff[n] = -1.0e+99;
    star[0].lambda[n] = -1.0e+99;
    if (tan_gamma>0.){
//      if (rocherad(system.a[n],companion[0].m[n]/star[0].m[n])<companion[0].r[n]) return 4;	//system merges
//      if (rocherad(system.a[n]*(1.0-pow(system.e[n],2)),companion[0].m[n]/star[0].m[n])<companion[0].r[n]) return 4;	//system merges	//use semi-major-axis after circularization
      if (system.peri[n]<companion[0].r[n]) return 4;	//system merges
    }
//    if (screen && (output=='M')) cout << endl << "Ejection speed of neutron star/black hole: " << v15/kms << "km/s\tEjection speed of companion star: " << v25/kms << "km/s" << endl;
    if (screen && (output=='M')) cout << endl << "Ejection speed of neutron star/black hole: " << v15 << "km/s\tEjection speed of companion star: " << v25 << "km/s" << endl;
    if (screen && (output=='T')) cout << "  disrupted -- w,theta,phi: " << w << " " << theta*180./M_PI << " " << phi*180./M_PI << endl;
    if (companion[0].stage[n]==-3){	//companion will explode too
      //determine values for the system
      system.M[n] = system.prim.m[n]+system.sec.m[n];	//get system mass
      system.qp[n] = system.prim.m[n]/system.sec.m[n];	//get mass ratio
      system.qs[n] = system.sec.m[n]/system.prim.m[n];	//get mass ratio
//      circularize(system);	// circularization
//      system.P[n] = 2.0*M_PI*sqrt(pow(system.a[n],3)/(G*system.M[n]));	//get period out of semi-major-axis
      system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
      system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
      system.phase[n] = -1;
      newphase(system);	//prepare new phase
//      cout  << "system.prim.stage[" << n << "]=" << system.prim.stage[n] << " system.sec.stage[" << n << "]=" << system.sec.stage[n] << " system.prim.stage[" << n+1 << "]=" << system.prim.stage[n+1] << " system.sec.stage[" << n+1 << "]=" << system.sec.stage[n+1] << endl;
      system.phase[n] = 0*supernova(system)-1;
      system.phase[n+1] = -1;
//      cout  << "system.phase[" << n << "]=" << system.phase[n] << " system.phase[" << n+1 << "]=" << system.phase[n+1] << endl;
    }
    return -1;
  }else{	//system remains bound
    // See eqs.(13-14) in J.G.Hills
    //vfactor = pow((w/v),2)+2.0*(w/v)*cos(theta);	//*(v/v)
    //system.a[n] = system.a[n-1]*(1.0-delta_M/(star[0].cm[n]+companion[0].m[n-1]))/(1.0-2.0*delta_M/(star[0].cm[n]+companion[0].m[n])-vfactor);	// J.G. Hills (13)
    //if (debug) cerr << endl << "vfactor=" << vfactor << " delta_M=" << delta_M << "Msun";
    // The pre-SN orbit is assumed circular. The momentum of star 1 is
    // equal (and opposite) to the momentum of star 2. Because star 1 loses
    // mass in the explosion, the system will receive a recoil velocity.
    // Also the orbit becomes eccentric.
    //Mreduced=star[0].m[n]*companion[0].m[n]/(star[0].m[n]+companion[0].m[n]);
    //L=(r)x(p) => L=a_i*Mreduced*vrel_perpendicular.
    // vrel.^2=v^2 + w^2 +2*v*w*cos(theta) J.G.Hills(14)
    //L_orb=(system.a[n-1]*Mreduced)*sqrt(pow(v+w*cos(theta),2)+pow(w*sin(theta)*sin(phi),2));
    //L_orb=(system.a[n-1]*Mreduced)*sqrt(pow(v+w*cos(theta),2)+pow(w*sin(theta)*sin(phi),2))*kms;
    //E_orb=-G*star[0].m[n]*companion[0].m[n]/system.a[n-1]+0.5*Mreduced*(pow(v,2)+pow(w,2)+2.0*v*w*cos(theta));
    //E_orb=-G*star[0].m[n]*companion[0].m[n]/system.a[n-1]+0.5*Mreduced*(pow(v,2)+pow(w,2)+2.0*v*w*cos(theta))*kms*kms;
    // The equation below is equivalent to the one above.
    //E_orb=-G*M1_remnant*M2_f/(2.*a_f) !!!
    //system.e[n]=sqrt(1.0+(2.0*E_orb*pow(L_orb,2))/(Mreduced*pow(G*star[0].m[n]*companion[0].m[n],2)));
    //system.peri[n] = system.a[n]*(1.0-system.e[n]);

    //calculate kick of the system(center of mass)
    vCM_x = (star[0].m[n]*(vcore+w*cos(theta))-star[0].m[n-1]*vcore)/(star[0].m[n]+companion[0].m[n]);
    vCM_y = star[0].m[n]*w*sin(theta)*cos(phi)/(star[0].m[n]+companion[0].m[n]);
    vCM_z = star[0].m[n]*w*sin(theta)*sin(phi)/(star[0].m[n]+companion[0].m[n]);
    vCM = sqrt(vCM_x*vCM_x+vCM_y*vCM_y+vCM_z*vCM_z);
    //save supernova data
    star[0].SN.vsys = vCM;
    /*// Calculate vCM perpendicular to the plane of the old orbit:
    vCM_90 = (star[0].m[n]*w*sin(theta)*sin(phi))/(star[0].m[n]+companion[0].m[n]);
    // The angle between old and new orbital plane: // J.G. Hills (24,26)
    angle = atan((w*sin(theta)*sin(phi))/sqrt(pow((v+w*cos(theta)),2)+pow(w*sin(theta)*cos(phi),2)));
    //angle=(angle/(2.*pi))*360.;
    if (debug) cerr << endl << "vCM=" << vCM/kms << "km/s vCM_90=" << vCM_90/kms << "km/s angle=" angle*180.0/M_PI << "degree";
    if (debug) cerr << endl << "vCM=" << vCM << "km/s vCM_90=" << vCM_90 << "km/s angle=" angle*180.0/M_PI << "degree";*/
    getrandomdirection(theta, phi);	//get kick direction of the system in the galaxy
    //change system velocity
    system.vx[n] += vCM*cos(theta);
    system.vy[n] += vCM*sin(theta)*cos(phi);
    system.vz[n] += vCM*sin(theta)*sin(phi);
/*    if (star[0].stage[n]==4){
      screen = true;
      cout << "vCM=" << vCM << "km/s theta=" << theta*180.0/M_PI << "degree phi=" << phi*180.0/M_PI << "degree system velocity{x,y,z}[" << n << "]={" << system.vx[n] << "," << system.vy[n] << "," << system.vz[n] << "}km/s" << endl;
    }*/

    //do not use values any longer
    star[0].cm[n] = star[0].m[n];
    star[0].llum[n] = -1.0e+99;
    star[0].lteff[n] = -1.0e+99;
    star[0].lambda[n] = -1.0e+99;
    system.stagechange = max(system.stagechange,nextstage(star[0],n));
    // Check coalesence:
    //if (system.peri[n]<companion[0].r[n]) return 4;
//    if (rocherad(system.a[n]*(1.0-pow(system.e[n],2)),companion[0].m[n]/star[0].m[n])<companion[0].r[n]) return 4;	//system merges	//use semi-major-axis after circularization
    if (rocherad(system.a[n]*(1.0-system.e[n]),companion[0].m[n]/star[0].m[n])<companion[0].r[n]) return 4;	//system merges	//use distance at periastron
    if (companion[0].stage[n]==-3){	//companion will explode too
      if (whichstar==1) return 25;	//secondary will explode too
      else if (whichstar==2) return 15;	//primary will explode too
      else return -1;	//undefined
    }
    if (companion[0].stage[n]>2) return 99;	//both stars are at their final stage
    return 0;
  }
}

void circularize(t_system& system){ // circularization
  int n=system.n-1;	//last entry index
  double ae;
  if (!separateevolution){
    if ((system.e[n]>0.0)&&(system.e[n]<=1.0)){
      while (system.phase[n-1]==0){
        ae = system.a[n]*(1.0-pow(system.e[n],2));	//calculate semi-major axis with angular momentum conservation and assume fixed masses
        if (screen && (system.phase[n-2]!=0)){
          if (output=='M') cout << " => circularization at system[" << n-1 << "]: a=" << system.a[n-1] << "Rsun e=" << system.e[n-1] << " a_new=" << system.a[n-1]*(1.0-pow(system.e[n-1],2)) << "Rsun";
          if (output=='T') cout << "circularization -- a: " << system.a[n-1] << " -> " << system.a[n-1]*(1.0-pow(system.e[n-1],2)) << endl;
        }
        system.e[n] = 0.0;	//set eccentricity to 0
        system.a[n] = ae;	//calculate semi-major axis with angular momentum conservation and assume fixed masses	\\/(1.0-pow(system.e[n],2))
        system.P[n] = 2.0*M_PI*sqrt(pow(system.a[n],3)/(G*system.M[n]));	//get period out of semi-major-axis
        system.rp[n] = rocherad(system.a[n],system.qp[n]);	//get Roche-Lobe-radius of primary
        system.rs[n] = rocherad(system.a[n],system.qs[n]);	//get Roche-Lobe-radius of secondary
        system.peri[n] = system.a[n];	//get peri-astro (after circularization e=0 => peri-astro=semi-major-axis)
        n--;
        if (n==0) break;	//no circularization of the initial system
      }
    }
  }
}

double convectiveenvelopefactor(double metal, double inimass){	//getting the factor in radius compared to maximum radius to have a convective envelope
  if (metal>0.007){	//Milky Way metallicity(~0.009)
//    return 0.4;
    if (inimass<1.816) return 0.01*(-1.208+1.796*inimass);
    else if (inimass<19) return 0.01*(2.054+20.22*pow(inimass-1.816,0.23));
    else return 0.01*84.26;
  }else if(metal>0.0035){	//LMC metallicity(~0.005)
    return 0.01*50;	//guess
  }else if(metal>0.001){	//SMC metallicity(~0.002)
    return 0.01*60;	//guess
  }else{	//IZw18 metallicity(~0.0002)
    return 0.01*70;
  }
}
