// rotation.cpp
// This subroutine calculates the influence of stellar rotation

//#include <iostream>	//in-/output to/from screen
//#include <math.h>	//provides the mathematical functions
#include "ComBinElib.h" // local header file

/*global variables*/
int tides=0;	//flag for tides

void RotationTides(t_system& system, double t){	//calculates the change of the stellar rotation indiuced by tides
  int n=system.n-1;	//last entry index
  int i,ip=0,is=0;	//loop
//  int jp=min(system.prim.last+1,n),js=min(system.sec.last+1,n);	//entry index (p,s denote primary/secondary)
  double ratioop,ratioos;	//ratio for primary/secondary track
  double xp,xs; 	//Mass loss factor, Langer 89
  double mpfcor,msfcor;	//Final Masses of Primary/Seconday CORrigated by enhanced mass los
  double superdmp=0.0,superdms=0.0;	//extra mass loss in case of supercritical rotation
  double saveomp,saveoms;	//angular velocity of primary and secondary
  //old values
  double Mold,saveOm;	//system: mass, mean orbital angular velocity
  double mpold,cpold,Rpold,Lpold,cfpold,ompold;	//primary: mass, core mass, radius, luminosity, concentration factor, angular velocity
  double msold,csold,Rsold,Lsold,cfsold,omsold;	//secondary: mass, core mass, radius, luminosity, concentration factor, angular velocity
  //new values
  double Mnew;	//system: mass
  double mpnew,cpnew,Rpnew,Lpnew,cfpnew,ompnew;	//primary: mass, core mass, radius, luminosity, concentration factor, angular velocity
  double msnew,csnew,Rsnew,Lsnew,cfsnew,omsnew;	//secondary: mass, core mass, radius, luminosity, concentration factor, angular velocity
  //time scales
  double tsyncp=0.0, tsyncs=0.0;	//synchronisation time scale of primary and secondary
  double tsyncpold=0.0,tsyncsold=0.0;	//old synchronisation time scale of primary and secondary
  double tsyncpnew=0.0,tsyncsnew=0.0;	//new synchronisation time scale of primary and secondary
  double tconv,fconv,tkc;
  double lastt,deltat;	//last time and time between two track points
  double anaus = 1.0;

  if ((t<0)||(tides==0)) return;	//without tides
  if (debug) cerr << endl << "t=" << t << "yr system.prim.t[" << n-1 << "]=" << system.prim.t[n-1] << "yr system.sec.t[" << n-1 << "]=" << system.sec.t[n-1] << "yr" << endl;

  //begin spin evolution algorithm during wind mass loss
  for (i=1; i<system.prim.track.n; i++){
    if (system.prim.track.t[i]>system.prim.t[n-1]){
      ip = i-1;
      break;
    }
  }
  if (ip==0) ip = system.prim.track.n-1;	//no track point found, use last track point
  for (i=1; i<system.sec.track.n; i++){
    if (system.sec.track.t[i]>system.sec.t[n-1]){
      is = i-1;
      break;
    }
  }
  if (is==0) is = system.sec.track.n-1;	//no track point found, use last track point
  if (debug) cerr << "ip=" << ip << " prim.last=" << system.prim.last << " is=" << is << " sec.last=" << system.sec.last << endl;

  //set first values
  Mold = system.M[n-1];
  mpold = system.prim.m[n-1];
  msold = system.sec.m[n-1];
  cpold = system.prim.cm[n-1];
  csold = system.sec.cm[n-1];
  Rpold = system.prim.r[n-1];
  Rsold = system.sec.r[n-1];
  Lpold = pow(10.0,system.prim.llum[n-1]);
  Lsold = pow(10.0,system.sec.llum[n-1]);
  cfpold = system.prim.cf[n-1];
  cfsold = system.sec.cf[n-1];
  ompold = fmax(system.prim.track.omega[ip],0.0);
  omsold = fmax(system.sec.track.omega[is],0.0);

  saveomp = fmax(system.prim.omega[n-1],0.0);
  saveoms = fmax(system.sec.omega[n-1],0.0);
  saveOm = fmax(2.0*M_PI/system.P[n-1],0.0);
  mpfcor = mpold;
  msfcor = msold;
  if (debug) cerr << "saveomp=" << saveomp << "/yr saveoms=" << saveoms << "/yr saveOm=" << saveOm << "/yr" << endl;
  if (system.prim.track.t[ip]==system.prim.t[n-1]) system.prim.track.omega[ip] = system.prim.omega[n-1];
  if (system.sec.track.t[is]==system.sec.t[n-1]) system.sec.track.omega[is] = system.sec.omega[n-1];

  //rotationally enhanced mass loss, eq.(3) in Langer98
  xp = pow(1.0-saveomp/omegacrit(mpold,Rpold,0.0),-0.43);
  xs = pow(1.0-saveoms/omegacrit(msold,Rsold,0.0),-0.43);
  lastt = system.prim.t[n-1];
  if (debug) cerr << "xp=" << xp << " xs=" << xs << endl;

  if (debug) cerr << "t=" << t << "yr system.prim.track.t[" << ip+1 << "]=" << system.prim.track.t[ip+1] << "yr system.sec.track.t[" << is+1 << "]=" << system.sec.track.t[is+1] << "yr" << endl;
  while (((system.prim.track.t[ip+1]<=t)||(system.sec.track.t[is+1]<=t))&&(ip+1<system.prim.track.n)&&(is+1<system.sec.track.n)){
    if (system.prim.track.t[ip+1] < system.sec.track.t[is+1]){	//next time on primary's track
      ip++;
      if (debug) cerr << "ip=" << ip;
      //get new values (interpolate for secondary)
      deltat = system.prim.track.t[ip] - lastt;
      ratioos = ratio(system.prim.track.t[ip],system.sec.track.t[is+1],system.sec.track.t[is]);
      if (debug) cerr << " deltat=" << deltat << "yr ratioos=" << ratioos << endl;
      mpnew = system.prim.track.m[ip];
      msnew = system.sec.track.m[is+1]+ratioos*(system.sec.track.m[is]-system.sec.track.m[is+1]);
      cpnew = system.prim.track.cm[ip];
      csnew = system.sec.track.cm[is+1]+ratioos*(system.sec.track.cm[is]-system.sec.track.cm[is+1]);
      Mnew = mpnew+msnew;
      Rpnew = system.prim.track.r[ip];
      Rsnew = system.sec.track.r[is+1]+ratioos*(system.sec.track.r[is]-system.sec.track.r[is+1]);
      Lpnew = pow(10.0,system.prim.track.llum[ip]);
      Lsnew = pow(10.0,system.sec.track.llum[is+1]+ratioos*(system.sec.track.llum[is]-system.sec.track.llum[is+1]));
      cfpnew = system.prim.track.cf[ip];
      cfsnew = system.sec.track.cf[is+1]+ratioos*(system.sec.track.cf[is]-system.sec.track.cf[is+1]);
      ompnew = fmax(system.prim.track.omega[ip],0.0);
      saveOm = square(Mnew/Mold)*saveOm;
      if ((ompold>0.0)&&(ompnew>0.0)){	//take relative change of angular momentum for single star evolution into account
        saveomp = saveomp * ompnew/ompold;
        if (debug) cerr << "saveomp=" << saveomp << "/yr ompold=" << ompold << "/yr ompnew=" << ompnew << "/yr" << endl;
      }
      //~ saveomp = saveomp*square(Rpold/Rpnew)*(cfpold*mpold-2.0/3.0*(mpold-mpnew))/(cfpnew*mpnew);
      //~ saveoms = saveoms*square(Rsold/Rsnew)*(cfsold*msold-2.0/3.0*(msold-msnew))/(cfsnew*msnew);
      if (system.prim.stage[n-1]!=6){
        saveomp = saveomp * cfpold/cfpnew*square(Rpold/Rpnew) * pow(mpnew/mpold,xp/0.75/(cfpold+cfpnew)-1.0);
        if (debug) cerr << "saveomp=" << saveomp << "/yr cfpold/cfpnew=" << cfpold/cfpnew << " Rpold/Rpnew=" << Rpold/Rpnew << " mpnew/mpold=" << mpnew/mpold << endl;
      }
      if (system.sec.stage[n-1]!=6){
        saveoms = saveoms * cfsold/cfsnew*square(Rsold/Rsnew) * pow(msnew/msold,xs/0.75/(cfsold+cfsnew)-1.0);
        if (debug) cerr << "saveoms=" << saveoms << "/yr cfsold/cfsnew=" << cfsold/cfsnew << " Rsold/Rsnew=" << Rsold/Rsnew << " msnew/msold=" << msnew/msold << endl;
      }
      if (system.phase[n-1]>=0){	//system is a binary
        //primary
        if (system.prim.stage[n-1]==0){	//H burning star
          //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
          tsyncpold = anaus/omegacrit(mpold,Rpold,0.0) * cfpold * square(mpold/msold) * pow(1.0+msold/mpold,-5.0/6.0) * 8.67064e+5 * pow(mpold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rpold,8.5);
          tsyncpnew = anaus/omegacrit(mpnew,Rpnew,0.0) * cfpnew * square(mpnew/msnew) * pow(1.0+msnew/mpnew,-5.0/6.0) * 8.67064e+5 * pow(mpnew,-2.84) * pow(cbrt(G*Mnew/square(saveOm))/Rpnew,8.5);
          if (debug) cerr << "anaus=" << anaus << " mpnew=" << mpnew << "Msun Rpnew=" << Rpnew << "Rsun cfpnew=" << cfpnew << " msnew=" << msnew << "Msun saveOm=" << saveOm << "/yr" << endl;
        }else if (system.prim.stage[n-1]==1){	//He burning star
          tconv = 0.4311 *cbrt((mpold-cpold)*square(Rpold)*0.125/Lpold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * mpold/(mpold-cpold);	//eq.(30) in Hurley,Tout,Pols02
          tsyncpold = tkc/3.0 * square(mpold/msold) * cfpold * square(G*Mold/square(saveOm)/cubic(Rpold));	//eq.(27) in Hurley,Tout,Pols02
          tconv = 0.4311 *cbrt((mpnew-cpnew)*square(Rpnew)*0.125/Lpnew);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * mpnew/(mpnew-cpnew);	//eq.(30) in Hurley,Tout,Pols02
          tsyncpnew = tkc/3.0 * square(mpnew/msnew) * cfpnew * square(G*Mnew/square(saveOm)/cubic(Rpnew));	//eq.(27) in Hurley,Tout,Pols02
        }else if (system.prim.stage[n-1]==2){	//naked He star
          //no tides
          tsyncpold = 1.0e99;
          tsyncpnew = 1.0e99;
        }else{
          //no tides
          tsyncpold = 1.0e99;
          tsyncpnew = 1.0e99;
        }
        //secondary
        if (system.sec.stage[n-1]==0){	//H burning star
          //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
          tsyncsold = anaus/omegacrit(msold,Rsold,0.0) * cfsold * square(msold/mpold) * pow(1.0+mpold/msold,-5.0/6.0) * 8.67064e+5 * pow(msold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rsold,8.5);
          tsyncsnew = anaus/omegacrit(msnew,Rsnew,0.0) * cfsnew * square(msnew/mpnew) * pow(1.0+mpnew/msnew,-5.0/6.0) * 8.67064e+5 * pow(msnew,-2.84) * pow(cbrt(G*Mnew/square(saveOm))/Rsnew,8.5);
          if (debug) cerr << "anaus=" << anaus << " msnew=" << msnew << "Msun Rsnew=" << Rsnew << "Rsun cfsnew=" << cfsnew << " mpnew=" << mpnew << "Msun saveOm=" << saveOm << "/yr" << endl;
        }else if (system.sec.stage[n-1] == 1){	//He burning star
          tconv = 0.4311 *cbrt((msold-csold)*square(Rsold)*0.125/Lsold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * msold/(msold-csold);	//eq.(30) in Hurley,Tout,Pols02
          tsyncsold = tkc/3.0 * square(msold/mpold) * cfsold * square(G*Mold/square(saveOm)/cubic(Rsold));	//eq.(27) in Hurley,Tout,Pols02
          tconv = 0.4311 *cbrt((msnew-csnew)*square(Rsnew)*0.125/Lsnew);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * msnew/(msnew-csnew);	//eq.(30) in Hurley,Tout,Pols02
          tsyncsnew = tkc/3.0 * square(msnew/mpnew) * cfsnew * square(G*Mnew/square(saveOm)/cubic(Rsnew));	//eq.(27) in Hurley,Tout,Pols02
        }else if (system.sec.stage[n-1] == 2){	//naked He star
          //no tides
          tsyncsold = 1.0e99;
          tsyncsnew = 1.0e99;
        }else{
          //no tides
          tsyncsold = 1.0e99;
          tsyncsnew = 1.0e99;
        }
        if (debug) cerr << "tsyncpold=" << tsyncpold << "yr tsyncpnew=" << tsyncpnew << "yr  tsyncsold=" << tsyncsold << "yr tsyncsnew=" << tsyncsnew << "yr" << endl;
        tsyncs = tsyncsnew;
        tsyncp = tsyncpnew;
        //solving eq.(27) in Hurley,Tout,Pols02
        saveomp = (saveomp-saveOm) * exp(-deltat/tsyncp) + saveOm;
        saveoms = (saveoms-saveOm) * exp(-deltat/tsyncs) + saveOm;
        if (debug) cerr << "saveomp=" << saveomp << "/yr saveoms=" << saveoms << "/yr saveOm=" << saveOm << "/yr deltat=" << deltat << "yr tsyncp=" << tsyncp << "yr tsyncs=" << tsyncs << "yr" << endl;
      }
      //enhanced mass loss by super-critical rotation
      superdmp=0.0;
      superdms=0.0;
      if (saveomp > omegacrit(mpnew,1.5*Rpnew,0.0)){
        superdmp = (0.99-saveomp/omegacrit(mpnew,Rpnew*rroteq(0.99),GammaEdd(Lpnew,mpnew,system.prim.stage[n-1])))*(cfpold+cfpnew)*0.75;
        saveomp = 0.99*omegacrit(mpnew,1.5*Rpnew,0.0);
        if (debug) cerr << "saveomp=" << saveomp << "/yr mpnew=" << mpnew << "Msun Rpnew=" << Rpnew << "Rsun" << endl;
      }
      if (saveoms > omegacrit(msnew,1.5*Rsnew,0.0)){
        superdms = (0.99-saveoms/omegacrit(msnew,Rsnew*rroteq(0.99),GammaEdd(Lsnew,msnew,system.sec.stage[n-1])))*(cfsold+cfsnew)*0.75;
        saveoms = 0.99*omegacrit(msnew,1.5*Rsnew,0.0);
        if (debug) cerr << "saveoms=" << saveoms << "/yr msnew=" << msnew << "Msun Rsnew=" << Rsnew << "Rsun" << endl;
      }
      system.prim.track.omega[ip] = saveomp;
      mpfcor += xp*(mpnew-mpold)+superdmp;
      msfcor += xs*(msnew-msold)+superdms;
      //~ cout << "p " << xp*(mpnew-mpold)+superdmp << " " << mpnew << " " << mpold << " " << endl;
      //~ cout << "s " << xs*(msnew-msold)+superdms << " " << msnew << " " << msold << " " << system.sec.track.m[is] << " " << system.sec.track.m[is+1] << " " << ratioos << endl;
      //cout << "s* " << msnew << " " << msold << " " << system.sec.track.m[is+1] << " " << system.sec.track.m[is] << " " << ratioos << endl;
      //~ cout << "p: " << ip << " " << saveomp << endl;

      //save old values
      Mold = Mnew;
      mpold = mpnew;
      msold = msnew;
      cpold = cpnew;
      csold = csnew;
      Rpold = Rpnew;
      Rsold = Rsnew;
      Lpold = Lpnew;
      Lsold = Lsnew;
      cfpold = cfpnew;
      cfsold = cfsnew;
      ompold = ompnew;
      //rotationally enhanced mass loss, eq.(3) in Langer98
      xp = pow(1.0-saveomp/omegacrit(mpold,Rpold,0.0),-0.43);
      xs = pow(1.0-saveoms/omegacrit(msold,Rsold,0.0),-0.43);
      lastt = system.prim.track.t[ip];
      if (debug) cerr << "xp=" << xp << " xs=" << xs << endl;
    }else if (system.prim.track.t[ip+1] > system.sec.track.t[is+1]){	//next time on secondary's track
      is++;
      if (debug) cerr << "is=" << is;
      //get new values (interpolate for primary)
      deltat = system.sec.track.t[is] - lastt;
      ratioop = ratio(system.sec.track.t[is],system.prim.track.t[ip+1],system.prim.track.t[ip]);
      if (debug) cerr << " deltat=" << deltat << "yr ratioop=" << ratioop << endl;
      mpnew = system.prim.track.m[ip+1]+ratioop*(system.prim.track.m[ip]-system.prim.track.m[ip+1]);
      msnew = system.sec.track.m[is];
      cpnew = system.prim.track.cm[ip+1]+ratioop*(system.prim.track.cm[ip]-system.prim.track.cm[ip+1]);
      csnew = system.sec.track.cm[is];
      Mnew = mpnew+msnew;
      Rpnew = system.prim.track.r[ip+1]+ratioop*(system.prim.track.r[ip]-system.prim.track.r[ip+1]);
      Rsnew = system.sec.track.r[is];
      Lpnew = pow(10.0,system.prim.track.llum[ip+1]+ratioop*(system.prim.track.llum[ip]-system.prim.track.llum[ip+1]));
      Lsnew = pow(10.0,system.sec.track.r[is]);
      cfpnew = system.prim.track.cf[ip+1]+ratioop*(system.prim.track.cf[ip]-system.prim.track.cf[ip+1]);
      cfsnew = system.sec.track.cf[is];
      omsnew = fmax(system.sec.track.omega[is],0.0);
      saveOm = square(Mnew/Mold)*saveOm;
      if ((omsold>0.0)&&(omsnew>0.0)){	//take relative change of angular momentum for single star evolution into account
        saveoms = saveoms * omsnew/omsold;
        if (debug) cerr << "saveoms=" << saveoms << "/yr omsold=" << omsold << "/yr omsnew=" << omsnew << "/yr" << endl;
      }
      //~ saveomp = saveomp*square(Rpold/Rpnew)*(cfpold*mpold-2.0/3.0*(mpold-mpnew))/(cfpnew*mpnew);
      //~ saveoms = saveoms*square(Rsold/Rsnew)*(cfsold*msold-2.0/3.0*(msold-msnew))/(cfsnew*msnew);
      if (system.prim.stage[n-1]!=6){
        saveomp = saveomp * cfpold/cfpnew*square(Rpold/Rpnew) * pow(mpnew/mpold,xp/0.75/(cfpold+cfpnew)-1.0);
        if (debug) cerr << "saveomp=" << saveomp << "/yr cfpold/cfpnew=" << cfpold/cfpnew << " Rpold/Rpnew=" << Rpold/Rpnew << " mpnew/mpold=" << mpnew/mpold << endl;
      }
      if (system.sec.stage[n-1]!=6){
        saveoms = saveoms * cfsold/cfsnew*square(Rsold/Rsnew) * pow(msnew/msold,xs/0.75/(cfsold+cfsnew)-1.0);
        if (debug) cerr << "saveoms=" << saveoms << "/yr cfsold/cfsnew=" << cfsold/cfsnew << " Rsold/Rsnew=" << Rsold/Rsnew << " msnew/msold=" << msnew/msold << endl;
      }
      if (system.phase[n-1]>=0){	//system is a binary
        //primary
        if (system.prim.stage[n-1]==0){	//H burning star
          //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
          tsyncpold = anaus/omegacrit(mpold,Rpold,0.0) * cfpold * square(mpold/msold) * pow(1.0+msold/mpold,-5.0/6.0) * 8.67064e+5 * pow(mpold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rpold,8.5);
          tsyncpnew = anaus/omegacrit(mpnew,Rpnew,0.0) * cfpnew * square(mpnew/msnew) * pow(1.0+msnew/mpnew,-5.0/6.0) * 8.67064e+5 * pow(mpnew,-2.84) * pow(cbrt(G*Mnew/square(saveOm))/Rpnew,8.5);
          if (debug) cerr << "anaus=" << anaus << " mpnew=" << mpnew << "Msun Rpnew=" << Rpnew << "Rsun cfpnew=" << cfpnew << " msnew=" << msnew << "Msun saveOm=" << saveOm << "/yr" << endl;
        }else if (system.prim.stage[n-1]==1){	//He burning star
          tconv = 0.4311 *cbrt((mpold-cpold)*square(Rpold)*0.125/Lpold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * mpold/(mpold-cpold);	//eq.(30) in Hurley,Tout,Pols02
          tsyncpold = tkc/3.0 * square(mpold/msold) * cfpold * square(G*Mold/square(saveOm)/cubic(Rpold));	//eq.(27) in Hurley,Tout,Pols02
          tconv = 0.4311 *cbrt((mpnew-cpnew)*square(Rpnew)*0.125/Lpnew);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * mpnew/(mpnew-cpnew);	//eq.(30) in Hurley,Tout,Pols02
          tsyncpnew = tkc/3.0 * square(mpnew/msnew) * cfpnew * square(G*Mnew/square(saveOm)/cubic(Rpnew));	//eq.(27) in Hurley,Tout,Pols02
        }else if (system.prim.stage[n-1]==2){	//naked He star
          //no tides
          tsyncpold = 1.0e99;
          tsyncpnew = 1.0e99;
        }else{
          //no tides
          tsyncpold = 1.0e99;
          tsyncpnew = 1.0e99;
        }
        //secondary
        if (system.sec.stage[n-1]==0){	//H burning star
          //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
          tsyncsold = anaus/omegacrit(msold,Rsold,0.0) * cfsold * square(msold/mpold) * pow(1.0+mpold/msold,-5.0/6.0) * 8.67064e+5 * pow(msold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rsold,8.5);
          tsyncsnew = anaus/omegacrit(msnew,Rsnew,0.0) * cfsnew * square(msnew/mpnew) * pow(1.0+mpnew/msnew,-5.0/6.0) * 8.67064e+5 * pow(msnew,-2.84) * pow(cbrt(G*Mnew/square(saveOm))/Rsnew,8.5);
          if (debug) cerr << "anaus=" << anaus << " msnew=" << msnew << "Msun Rsnew=" << Rsnew << "Rsun cfsnew=" << cfsnew << " mpnew=" << mpnew << "Msun saveOm=" << saveOm << "/yr" << endl;
        }else if (system.sec.stage[n-1] == 1){	//He burning star
          tconv = 0.4311 *cbrt((msold-csold)*square(Rsold)*0.125/Lsold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * msold/(msold-csold);	//eq.(30) in Hurley,Tout,Pols02
          tsyncsold = tkc/3.0 * square(msold/mpold) * cfsold * square(G*Mold/square(saveOm)/cubic(Rsold));	//eq.(27) in Hurley,Tout,Pols02
          tconv = 0.4311 *cbrt((msnew-csnew)*square(Rsnew)*0.125/Lsnew);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * msnew/(msnew-csnew);	//eq.(30) in Hurley,Tout,Pols02
          tsyncsnew = tkc/3.0 * square(msnew/mpnew) * cfsnew * square(G*Mnew/square(saveOm)/cubic(Rsnew));	//eq.(27) in Hurley,Tout,Pols02
        }else if (system.sec.stage[n-1] == 2){	//naked He star
          //no tides
          tsyncsold = 1.0e99;
          tsyncsnew = 1.0e99;
        }else{
          //no tides
          tsyncsold = 1.0e99;
          tsyncsnew = 1.0e99;
        }
        if (debug) cerr << "tsyncpold=" << tsyncpold << "yr tsyncpnew=" << tsyncpnew << "yr  tsyncsold=" << tsyncsold << "yr tsyncsnew=" << tsyncsnew << "yr" << endl;
        tsyncs = tsyncsnew;
        tsyncp = tsyncpnew;
        //solving eq.(27) in Hurley,Tout,Pols02
        saveomp = (saveomp-saveOm) * exp(-deltat/tsyncp) + saveOm;
        saveoms = (saveoms-saveOm) * exp(-deltat/tsyncs) + saveOm;
        if (debug) cerr << "saveomp=" << saveomp << "/yr saveoms=" << saveoms << "/yr saveOm=" << saveOm << "/yr deltat=" << deltat << "yr tsyncp=" << tsyncp << "yr tsyncs=" << tsyncs << "yr" << endl;
      }
      //enhanced mass loss by super-critical rotation
      superdmp=0.0;
      superdms=0.0;
      if (saveomp > omegacrit(mpnew,1.5*Rpnew,0.0)){
        superdmp = (0.99-saveomp/omegacrit(mpnew,Rpnew*rroteq(0.99),GammaEdd(Lpnew,mpnew,system.prim.stage[n-1])))*(cfpold+cfpnew)*0.75;
        saveomp = 0.99*omegacrit(mpnew,1.5*Rpnew,0.0);
      }
      if (saveoms > omegacrit(msnew,1.5*Rsnew,0.0)){
        superdms = (0.99-saveoms/omegacrit(msnew,Rsnew*rroteq(0.99),GammaEdd(Lsnew,msnew,system.sec.stage[n-1])))*(cfsold+cfsnew)*0.75;
        saveoms = 0.99*omegacrit(msnew,1.5*Rsnew,0.0);
      }
      if (debug) cerr << "saveoms=" << saveoms << "/yr msnew=" << msnew << "Msun Rsnew=" << Rsnew << "Rsun" << endl;
      system.sec.track.omega[is] = saveoms;
      mpfcor += xp*(mpnew-mpold)+superdmp;
      msfcor += xs*(msnew-msold)+superdms;
      //~ cout << "p " << xp*(mpnew-mpold)+superdmp << " " << mpnew << " " << mpold << " " << system.prim.track.m[ip] << " " << system.prim.track.m[ip+1] << " " << ratioop << endl;
      //~ cout << "s " << xs*(msnew-msold)+superdms << " " << msnew << " " << msold << " " << endl;

      //~ cout << "s: " << is << " " << saveoms << endl;

      //save old values
      Mold = Mnew;
      mpold = mpnew;
      msold = msnew;
      cpold = cpnew;
      csold = csnew;
      Rpold = Rpnew;
      Rsold = Rsnew;
      Lpold = Lpnew;
      Lsold = Lsnew;
      cfpold = cfpnew;
      cfsold = cfsnew;
      omsold = omsnew;
      //rotationally enhanced mass loss, eq.(3) in Langer98
      xp = pow(1.0-saveomp/omegacrit(mpold,Rpold,0.0),-0.43);
      xs = pow(1.0-saveoms/omegacrit(msold,Rsold,0.0),-0.43);
      lastt = system.sec.track.t[is];
      if (debug) cerr << "xp=" << xp << " xs=" << xs << endl;
    }else if (system.prim.track.t[ip+1] == system.sec.track.t[is+1]){	//next time on primary's and secondary's track
      ip++;
      is++;
      if (debug) cerr << "ip=" << ip << " is=" << is;
      //get new values (interpolate for primary)
      deltat = system.prim.track.t[ip] - lastt;
      if (debug) cerr << "deltat=" << deltat << "yr" << endl;
      mpnew = system.prim.track.m[ip];
      msnew = system.sec.track.m[is];
      cpnew = system.prim.track.cm[ip];
      csnew = system.sec.track.cm[is];
      Mnew = mpnew+msnew;
      Rpnew = system.prim.track.r[ip];
      Rsnew = system.sec.track.r[is];
      Lpnew = pow(10.0,system.prim.track.llum[ip]);
      Lsnew = pow(10.0,system.sec.track.llum[is]);
      cfpnew = system.prim.track.cf[ip];
      cfsnew = system.sec.track.cf[is];
      ompnew = fmax(system.prim.track.omega[ip],0.0);
      omsnew = fmax(system.sec.track.omega[is],0.0);
      saveOm = square(Mnew/Mold)*saveOm;
      if ((ompold>0.0)&&(ompnew>0.0)){	//take relative change of angular momentum for single star evolution into account
        saveomp = saveomp * ompnew/ompold;
        if (debug) cerr << "saveomp=" << saveomp << "/yr ompold=" << ompold << "/yr ompnew=" << ompnew << "/yr" << endl;
      }
      if ((omsold>0.0)&&(omsnew>0.0)){	//take relative change of angular momentum for single star evolution into account
        saveoms = saveoms * omsnew/omsold;
        if (debug) cerr << "saveoms=" << saveoms << "/yr omsold=" << omsold << "/yr omsnew=" << omsnew << "/yr" << endl;
      }
      //~ saveomp = saveomp*square(Rpold/Rpnew)*(cfpold*mpold-2.0/3.0*(mpold-mpnew))/(cfpnew*mpnew);
      //~ saveoms = saveoms*square(Rsold/Rsnew)*(cfsold*msold-2.0/3.0*(msold-msnew))/(cfsnew*msnew);
      if (system.prim.stage[n-1]!=6){
        saveomp = saveomp * cfpold/cfpnew*square(Rpold/Rpnew) * pow(mpnew/mpold,xp/0.75/(cfpold+cfpnew)-1.0);
        if (debug) cerr << "saveomp=" << saveomp << "/yr cfpold/cfpnew=" << cfpold/cfpnew << " Rpold/Rpnew=" << Rpold/Rpnew << " mpnew/mpold=" << mpnew/mpold << endl;
      }
      if (system.sec.stage[n-1]!=6){
        saveoms = saveoms * cfsold/cfsnew*square(Rsold/Rsnew) * pow(msnew/msold,xs/0.75/(cfsold+cfsnew)-1.0);
        if (debug) cerr << "saveoms=" << saveoms << "/yr cfsold/cfsnew=" << cfsold/cfsnew << " Rsold/Rsnew=" << Rsold/Rsnew << " msnew/msold=" << msnew/msold << endl;
      }
      if (system.phase[n-1]>=0){	//system is a binary
        //primary
        if (system.prim.stage[n-1]==0){	//H burning star
          //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
          tsyncpold = anaus/omegacrit(mpold,Rpold,0.0) * cfpold * square(mpold/msold) * pow(1.0+msold/mpold,-5.0/6.0) * 8.67064e+5 * pow(mpold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rpold,8.5);
          tsyncpnew = anaus/omegacrit(mpnew,Rpnew,0.0) * cfpnew * square(mpnew/msnew) * pow(1.0+msnew/mpnew,-5.0/6.0) * 8.67064e+5 * pow(mpnew,-2.84) * pow(cbrt(G*Mnew/square(saveOm))/Rpnew,8.5);
          if (debug) cerr << "anaus=" << anaus << " mpnew=" << mpnew << "Msun Rpnew=" << Rpnew << "Rsun cfpnew=" << cfpnew << " msnew=" << msnew << "Msun saveOm=" << saveOm << "/yr" << endl;
        }else if (system.prim.stage[n-1]==1){	//He burning star
          tconv = 0.4311 *cbrt((mpold-cpold)*square(Rpold)*0.125/Lpold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * mpold/(mpold-cpold);	//eq.(30) in Hurley,Tout,Pols02
          tsyncpold = tkc/3.0 * square(mpold/msold) * cfpold * square(G*Mold/square(saveOm)/cubic(Rpold));	//eq.(27) in Hurley,Tout,Pols02
          tconv = 0.4311 *cbrt((mpnew-cpnew)*square(Rpnew)*0.125/Lpnew);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * mpnew/(mpnew-cpnew);	//eq.(30) in Hurley,Tout,Pols02
          tsyncpnew = tkc/3.0 * square(mpnew/msnew) * cfpnew * square(G*Mnew/square(saveOm)/cubic(Rpnew));	//eq.(27) in Hurley,Tout,Pols02
        }else if (system.prim.stage[n-1]==2){	//naked He star
          //no tides
          tsyncpold = 1.0e99;
          tsyncpnew = 1.0e99;
        }else{
          //no tides
          tsyncpold = 1.0e99;
          tsyncpnew = 1.0e99;
        }
        //secondary
        if (system.sec.stage[n-1]==0){	//H burning star
          //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
          tsyncsold = anaus/omegacrit(msold,Rsold,0.0) * cfsold * square(msold/mpold) * pow(1.0+mpold/msold,-5.0/6.0) * 8.67064e+5 * pow(msold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rsold,8.5);
          tsyncsnew = anaus/omegacrit(msnew,Rsnew,0.0) * cfsnew * square(msnew/mpnew) * pow(1.0+mpnew/msnew,-5.0/6.0) * 8.67064e+5 * pow(msnew,-2.84) * pow(cbrt(G*Mnew/square(saveOm))/Rsnew,8.5);
          if (debug) cerr << "anaus=" << anaus << " msnew=" << msnew << "Msun Rsnew=" << Rsnew << "Rsun cfsnew=" << cfsnew << " mpnew=" << mpnew << "Msun saveOm=" << saveOm << "/yr" << endl;
        }else if (system.sec.stage[n-1] == 1){	//He burning star
          tconv = 0.4311 *cbrt((msold-csold)*square(Rsold)*0.125/Lsold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * msold/(msold-csold);	//eq.(30) in Hurley,Tout,Pols02
          tsyncsold = tkc/3.0 * square(msold/mpold) * cfsold * square(G*Mold/square(saveOm)/cubic(Rsold));	//eq.(27) in Hurley,Tout,Pols02
          tconv = 0.4311 *cbrt((msnew-csnew)*square(Rsnew)*0.125/Lsnew);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
          fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
          tkc = 10.5 * tconv/fconv * msnew/(msnew-csnew);	//eq.(30) in Hurley,Tout,Pols02
          tsyncsnew = tkc/3.0 * square(msnew/mpnew) * cfsnew * square(G*Mnew/square(saveOm)/cubic(Rsnew));	//eq.(27) in Hurley,Tout,Pols02
        }else if (system.sec.stage[n-1] == 2){	//naked He star
          //no tides
          tsyncsold = 1.0e99;
          tsyncsnew = 1.0e99;
        }else{
          //no tides
          tsyncsold = 1.0e99;
          tsyncsnew = 1.0e99;
        }
        if (debug) cerr << "tsyncpold=" << tsyncpold << "yr tsyncpnew=" << tsyncpnew << "yr  tsyncsold=" << tsyncsold << "yr tsyncsnew=" << tsyncsnew << "yr" << endl;
        tsyncs = tsyncsnew;
        tsyncp = tsyncpnew;
        //solving eq.(27) in Hurley,Tout,Pols02
        saveomp = (saveomp-saveOm) * exp(-deltat/tsyncp) + saveOm;
        saveoms = (saveoms-saveOm) * exp(-deltat/tsyncs) + saveOm;
        if (debug) cerr << "saveomp=" << saveomp << "/yr saveoms=" << saveoms << "/yr saveOm=" << saveOm << "/yr deltat=" << deltat << "yr tsyncp=" << tsyncp << "yr tsyncs=" << tsyncs << "yr" << endl;
      }
      //enhanced mass loss by super-critical rotation
      superdmp=0.0;
      superdms=0.0;
      if (saveomp > omegacrit(mpnew,1.5*Rpnew,0.0)){
        superdmp = (0.99-saveomp/omegacrit(mpnew,Rpnew*rroteq(0.99),GammaEdd(Lpnew,mpnew,system.prim.stage[n-1])))*(cfpold+cfpnew)*0.75;
        saveomp = 0.99*omegacrit(mpnew,1.5*Rpnew,0.0);
      }
      if (saveoms > omegacrit(msnew,1.5*Rsnew,0.0)){
        superdms = (0.99-saveoms/omegacrit(msnew,Rsnew*rroteq(0.99),GammaEdd(Lsnew,msnew,system.sec.stage[n-1])))*(cfsold+cfsnew)*0.75;
        saveoms = 0.99*omegacrit(msnew,1.5*Rsnew,0.0);
      }
      if (debug) cerr << "saveomp=" << saveomp << "/yr mpnew=" << mpnew << "Msun Rpnew=" << Rpnew << "Rsun" << endl;
      if (debug) cerr << "saveoms=" << saveoms << "/yr msnew=" << msnew << "Msun Rsnew=" << Rsnew << "Rsun" << endl;
      system.prim.track.omega[ip] = saveomp;
      system.sec.track.omega[is] = saveoms;
      mpfcor += xp*(mpnew-mpold)+superdmp;
      msfcor += xs*(msnew-msold)+superdms;
      //~ cout << "p " << xp*(mpnew-mpold)+superdmp << endl;
      //~ cout << "s " << xs*(msnew-msold)+superdms << endl;

      //save old values
      Mold = Mnew;
      mpold = mpnew;
      msold = msnew;
      cpold = cpnew;
      csold = csnew;
      Rpold = Rpnew;
      Rsold = Rsnew;
      Lpold = Lpnew;
      Lsold = Lsnew;
      cfpold = cfpnew;
      cfsold = cfsnew;
      ompold = ompnew;
      omsold = omsnew;
      //rotationally enhanced mass loss, eq.(3) in Langer98
      xp = pow(1.0-saveomp/omegacrit(mpold,Rpold,0.0),-0.43);
      xs = pow(1.0-saveoms/omegacrit(msold,Rsold,0.0),-0.43);
      lastt = system.prim.track.t[ip];
      if (debug) cerr << "xp=" << xp << " xs=" << xs << endl;
    }
    //~ cout << tsyncp << " " << mpold << endl;
    //~ cout << check << " " << t << " " << system.prim.track.t[ip] << " " << system.sec.track.t[is] << endl;
  }
  ompnew = fmax(system.prim.omega[n],0.0);
  omsnew = fmax(system.sec.omega[n],0.0);
  saveOm = square(system.M[n]/Mold)*saveOm;
  if ((ompold>0.0)&&(ompnew>0.0)){	//take relative change of angular momentum for single star evolution into account
    saveomp = saveomp * ompnew/ompold;
    if (debug) cerr << "saveomp=" << saveomp << "/yr ompold=" << ompold << "/yr ompnew=" << ompnew << "/yr" << endl;
  }
  if ((omsold>0.0)&&(omsnew>0.0)){	//take relative change of angular momentum for single star evolution into account
    saveoms = saveoms * omsnew/omsold;
    if (debug) cerr << "saveoms=" << saveoms << "/yr omsold=" << omsold << "/yr omsnew=" << omsnew << "/yr" << endl;
  }
  deltat = t - lastt;
  //~ saveomp = saveomp*square(Rpold/system.prim.r[n])*(cfpold*mpold-2.0/3.0*(mpold-system.prim.m[n]))/(system.prim.cf[n]*system.prim.m[n]);
  if (system.prim.stage[n-1]!=6){
    saveomp = saveomp * cfpold/system.prim.cf[n]*square(Rpold/system.prim.r[n]) * pow(system.prim.m[n]/mpold,xp/0.75/(cfpold+system.prim.cf[n])-1.0);
    if (debug) cerr << "saveomp=" << saveomp << "/yr cfpold/system.prim.cf[n]=" << cfpold/system.prim.cf[n] << " Rpold/system.prim.r[n]=" << Rpold/system.prim.r[n] << " system.prim.m[n]/mpold=" << system.prim.m[n]/mpold << endl;
  }
  if (system.sec.stage[n-1]!=6){
    saveoms = saveoms * cfsold/system.sec.cf[n]*square(Rsold/system.sec.r[n]) * pow(system.sec.m[n]/msold,xs/0.75/(cfsold+system.sec.cf[n])-1.0);
    if (debug) cerr << "saveoms=" << saveoms << "/yr cfsold/system.sec.cf[n]=" << cfsold/system.sec.cf[n] << " Rsold/system.sec.r[n]=" << Rsold/system.sec.r[n] << " system.sec.m[n]/msold=" << system.sec.m[n]/msold << endl;
  }
  if (system.phase[n-1]>=0){	//system is a binary
    //primary
    if (system.prim.stage[n-1]==0){	//H burning star
      //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
      tsyncpold = anaus/omegacrit(mpold,Rpold,0.0) * cfpold * square(mpold/msold) * pow(1.0+msold/mpold,-5.0/6.0) * 8.67064e+5 * pow(mpold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rpold,8.5);
      tsyncpnew = anaus/omegacrit(system.prim.m[n],system.prim.r[n],0.0) * system.prim.cf[n] * square(system.prim.m[n]/system.sec.m[n]) * pow(1.0+system.sec.m[n]/system.prim.m[n],-5.0/6.0) * 8.67064e+5 * pow(system.prim.m[n],-2.84) * pow(cbrt(G*system.M[n]/square(saveOm))/system.prim.r[n],8.5);
    }else if (system.prim.stage[n-1]==1){	//He burning star
      tconv = 0.4311 *cbrt((mpold-cpold)*square(Rpold)*0.125/Lpold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
      fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
      tkc = 10.5 * tconv/fconv * mpold/(mpold-cpold);	//eq.(30) in Hurley,Tout,Pols02
      tsyncpold = tkc/3.0 * square(mpold/msold) * cfpold * square(G*Mold/square(saveOm)/cubic(Rpold));	//eq.(27) in Hurley,Tout,Pols02
      tconv = 0.4311 *cbrt((system.prim.m[n]-system.prim.cm[n])*square(system.prim.r[n])*0.125/pow(10.0,system.prim.llum[n]));	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
      fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
      tkc = 10.5 * tconv/fconv * system.prim.m[n]/(system.prim.m[n]-system.prim.cm[n]);	//eq.(30) in Hurley,Tout,Pols02
      tsyncpnew = tkc/3.0 * square(system.prim.m[n]/system.sec.m[n]) * system.prim.cf[n] * square(G*system.M[n]/square(saveOm)/cubic(system.prim.r[n]));	//eq.(27) in Hurley,Tout,Pols02
    }else if (system.prim.stage[n-1]==2){	//naked He star
      //no tides
      tsyncpold = 1.0e99;
      tsyncpnew = 1.0e99;
    }else{
      //no tides
      tsyncpold = 1.0e99;
      tsyncpnew = 1.0e99;
    }
    //secondary
    if (system.sec.stage[n-1]==0){	//H burning star
      //syncronisation time scale, eq.(44) in Hurley,Tout,Pols02
      tsyncsold = anaus/omegacrit(msold,Rsold,0.0) * cfsold * square(msold/mpold) * pow(1.0+mpold/msold,-5.0/6.0) * 8.67064e+5 * pow(msold,-2.84) * pow(cbrt(G*Mold/square(saveOm))/Rsold,8.5);
      tsyncsnew = anaus/omegacrit(system.sec.m[n],system.sec.r[n],0.0) * system.sec.cf[n] * square(system.sec.m[n]/system.prim.m[n]) * pow(1.0+system.prim.m[n]/system.sec.m[n],-5.0/6.0) * 8.67064e+5 * pow(system.sec.m[n],-2.84) * pow(cbrt(G*system.M[n]/square(saveOm))/system.sec.r[n],8.5);
    }else if (system.sec.stage[n-1] == 1){	//He burning star
      tconv = 0.4311 *cbrt((msold-csold)*square(Rsold)*0.125/Lsold);	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
      fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
      tkc = 10.5 * tconv/fconv * msold/(msold-csold);	//eq.(30) in Hurley,Tout,Pols02
      tsyncsold = tkc/3.0 * square(msold/mpold) * cfsold * square(G*Mold/square(saveOm)/cubic(Rsold));	//eq.(27) in Hurley,Tout,Pols02
      tconv = 0.4311 *cbrt((system.sec.m[n]-system.sec.cm[n])*square(system.sec.r[n])*0.125/pow(10.0,system.sec.llum[n]));	//eq.(31) in Hurley,Tout,Pols02 with Renv = 0.5*R
      fconv = fmin(1.0,square(M_PI/fabs(saveOm-saveomp)/tconv));	//eq.(32) in Hurley,Tout,Pols02 with use of eq.(33)
      tkc = 10.5 * tconv/fconv * system.sec.m[n]/(system.sec.m[n]-system.sec.cm[n]);	//eq.(30) in Hurley,Tout,Pols02
      tsyncsnew = tkc/3.0 * square(system.sec.m[n]/system.prim.m[n]) * system.sec.cf[n] * square(G*system.M[n]/square(saveOm)/cubic(system.sec.r[n]));	//eq.(27) in Hurley,Tout,Pols02
    }else if (system.sec.stage[n-1] == 2){	//naked He star
      //no tides
      tsyncsold = 1.0e99;
      tsyncsnew = 1.0e99;
    }else{
      //no tides
      tsyncsold = 1.0e99;
      tsyncsnew = 1.0e99;
    }
    if (debug) cerr << "tsyncpold=" << tsyncpold << "yr tsyncpnew=" << tsyncpnew << "yr  tsyncsold=" << tsyncsold << "yr tsyncsnew=" << tsyncsnew << "yr" << endl;
    tsyncs = tsyncsnew;
    tsyncp = tsyncpnew;
    //solving eq.(27) in Hurley,Tout,Pols02
    saveomp = fmin((saveomp-saveOm) * exp(-deltat/tsyncp) + saveOm,0.99*omegacrit(system.prim.m[n],system.prim.r[n],0.0));
    saveoms = fmin((saveoms-saveOm) * exp(-deltat/tsyncs) + saveOm,0.99*omegacrit(system.sec.m[n],system.sec.r[n],0.0));
    if (debug) cerr << "saveomp=" << saveomp << "/yr saveoms=" << saveoms << "/yr saveOm=" << saveOm << "/yr deltat=" << deltat << "yr tsyncp=" << tsyncp << "yr tsyncs=" << tsyncs << "yr" << endl;
  }
  if (system.prim.stage[n-1]!=6) system.prim.omega[n] = saveomp;
  if (system.sec.stage[n-1]!=6) system.sec.omega[n] = saveoms;

  if ((mpfcor < 0.9*system.prim.m[n])&& screen){
    cout << " mpfcor/system.prim.m[" << n << "]=" << mpfcor/system.prim.m[n] << endl;
  }
  if ((msfcor < 0.9*system.sec.m[n])&& screen){
    cout << " mpfcor/system.sec.m[" << n << "]=" << msfcor/system.sec.m[n] << endl;
  }
  //end spin

  return;
}

double impactJ(double ra, double md, double ma, double a){	//get the specific angular momentum transfered on the accretor during RLO; see Lubow&Shu 1975, Ulrich&Burger 1976
  double q = ma/md;	//ma=mass accretor, md=mass donor, ra=radius accretor
  double rm;		//minimal radius for an accretion disk
  if (q<0.0667) rm = 0.0425*a*pow(0.0667+0.0667*0.0667,0.25);	//edge-values
  else if (q>15.0) rm = 0.0425*a*pow(15.0+15.0*15.0,0.25);
  else rm = 0.0425*a*pow(q+q*q,0.25);
  if (ra>rm) return sqrt(1.7*G*ma*rm);	//ballisitc impact of the stream on accretor
  else return sqrt(G*ma*ra);		//accretion stream collides with itself, Keplerian accretion disk
}

double rroteq(double om){
  return 1.0;	//pow(10., om*(1.7539e-2*tan(1.422*om) + 4.511e-2*sin(om)));
}

double omegacrit(double M, double R, double Gamma){
  return sqrt(G*M*(1.0-Gamma)/R)/R;	//see e.g. Langer 1998, eq. (2): v_crit^2=G*M*(1-Gamma)/R
}
/*double omegacr(double M, double R, double L, int stage){
  double kappa;	//opacity
  L = pow(10.0,L);
  if(stage == 0 || stage == 1) kappa = 0.2*(1.0+0.74)*(1.0e3*Msun)/(1.0e4*Rsun*Rsun);
  else if(stage == 2) kappa = 0.2*(1.0e3*Msun)/(1.0e4*Rsun*Rsun);
  else kappa = 0.;
  return sqrt(M*G/pow(3.0*R/2.0,3)*(1.0-kappa*L*Lsun/(4.0*M_PI*c*G*M)));
}*/

double GammaEdd(double L, double M, int stage){
  double kappa;
  L = pow(10.0,L);
  if((stage == 0)||(stage == 1)) kappa = 0.2*(1.0+0.74)*(1.0e3*Msun)/(1.0e4*Rsun*Rsun);	//kappa=0.2*(1+X) cm^2/g
  else if(stage == 2) kappa = 0.2*(1.0e3*Msun)/(1.0e4*Rsun*Rsun);
  else kappa = 0.0;
  return kappa*L*Lsun/(4.0*M_PI*c*G*M);
}

double getOmegaMassloss(t_star& star, int n){
  /*double xp;
  if (star.stage[n-1]!=6){
    xp = pow(1.0-star.omega[n-1]*star.r[n-1]/sqrt(star.m[n-1]*G/star.r[n-1]),-0.43);
    if (debug) cerr << "star.omega[n-1]=" << star.omega[n-1] << "/yr star.cf[n-1]/star.track.cf[star.last]=" << star.cf[n-1]/star.track.cf[star.last] << " star.r[n-1]/star.track.r[star.last]=" << star.r[n-1]/star.track.r[star.last] << " star.track.m[star.last]/star.m[n-1]=" << star.track.m[star.last]/star.m[n-1] << " xp=" << xp << endl;
    return star.omega[n-1] * star.cf[n-1]/star.track.cf[star.last]*square(star.r[n-1]/star.track.r[star.last]) * pow(star.track.m[star.last]/star.m[n-1],xp/0.75/(star.cf[n-1]+star.track.cf[star.last])-1.0);
  }else*/{
    return star.omega[n-1];
  }
} 
