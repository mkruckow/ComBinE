//      radii.cpp
//      This subroutine calculates the critical Roche-lobe radius, when
//      a and q is given, Eggleton (1983). q=M_donor/M_accretor

//#include <iostream>	//in-/output to/from screen
//#include <math.h>	//provides the mathematical functions
#include "ComBinElib.h" // local header file

//using namespace std; // must be used with iostream to enable write output. TT2013

double rocherad(double a, double q){
//       cout << endl << "The Rochelobe radius is: " 
//      << a*0.49*pow(q,2./3.)/(0.6*pow(q,2./3.)+log(1.+pow(q,1./3.)));
//       cout << endl << "for a0=" << a << " and q=" << q;
//     Note here, q=M_donor/M_accretor:
  double cbrtq=cbrt(q);
  return a*0.49*square(cbrtq)/(0.6*square(cbrtq)+log(1.0+cbrtq));
}

double RLFT(t_system system, int index, double& r, double rochemin){
//RLFT: Roche-lobe-filling time
  t_star star;	//star
  t_star companion;	//companion star
  int n=system.n-1;	//last entry index
  int cloop=0,cjump=0;	//counter of the loop, counter of jumps
  int j,jr=0;	//time index
  int jmax;	//maximal index of companion
  double roche;	//Roche lobe in Rsun
  double rocheold,rocheolder=0.0;	//last and before last Roche lobe in Rsun
  double rratio=1.0;	//ratio between the two nearest radii used for the star
  double t,told;	//time and last approximated time in yr
  double tratio;	//ratio between the two nearest times used for the companion
  double nstarmass,ncompanionmass,nq;	//new star/companion mass in Msun, new mass ratio
  double na;	//new semi-major-axis in Rsun
  double acirc=system.a[n]*(1.0-square(system.e[n]));	//circuarised semi-major-axis
  double ts=1.0e+99,tc=1.0e+99;	//times before stage change of star/companion
//  double na,a;	//new and last semi-major-axis in Rsun
//  double rad;	//radius in Rsun
//  double nperi;	//new periastron in Rsun
  bool constantradius=false;	//if the stellar radius is constant -> Roche lobe is used for interpolation
  bool extrapolation=false;	//if radius is extrapolated

  if (index==1){	//use primary star
    star = system.prim;
    companion = system.sec;
    roche = system.rp[n];
  }else{	//use secondary star
    star = system.sec;
    companion = system.prim;
    roche = system.rs[n];
  }
  //get maximal times
  if (star.stage[n]==0) ts = star.track.t[star.track.TAMS];	//maximal time is at TAMS
  else ts = star.track.t[star.track.n-1];	//maximal time is at the end of the track
  if (companion.stage[n]==0) jmax = companion.track.TAMS;	//maximal index is at TAMS
  else jmax = companion.track.n-1;	//maximal index is at the end of the track
  tc = companion.track.t[jmax];	//maximal time of the companion

  if (star.r[n]>roche){	//instantaneous Roche-lobe overflow
    r = star.r[n];
    return star.t[n];
  }
  if ((roche>star.track.r[star.rmax])&&(rochemin<star.track.r[star.rmax])) roche = fmax(rochemin,1.001*star.r[n]);
  if (star.last>=star.rmax){	//check if time of maximal radius is in the past
    if (debug) cerr << "star.last=" << star.last << " star.rmax=" << star.rmax << endl;
    return 1.0e+99;	//no Roche-lobe overflow
  }
  na = acirc;
//  nperi = system.peri[n];
  t = star.t[n];
  if (debug){
    cerr << "roche=" << roche << "Rsun rochemin=" << rochemin << "Rsun na=" << na << "Rsun t=" << t << "yr star.last=" << star.last << " star.rmax=" << star.rmax << " ts=" << ts << "yr";
    if (star.stage[n]==0) cerr << " TAMS:" << star.track.TAMS;
    else cerr << " end:" << star.track.n-1;
    cerr << " tc=" << tc << "yr";
    if (companion.stage[n]==0) cerr << " TAMS:" << companion.track.TAMS;
    else cerr << " end:" << companion.track.n-1;
    cerr << endl;
  }
  if (t!=companion.t[n]) cerr << "#Error: System is not synchronous." << endl;
  told = t;
  rocheold = roche;
  do{
    if (star.track.r[star.rmax]<fmin(roche,rochemin)){	//check if maximal stellar radius can fill the minimal Roche lobe
      if (debug) cerr << "star.track.r[star.rmax=" << star.rmax << "]=" << star.track.r[star.rmax] << "Rsun roche=" << roche << "Rsun rochemin=" << rochemin << "Rsun" << endl;
      return 1.0e+99;	//no Roche-lobe overflow
    }
    if ((fabs(rocheolder-roche)<accuracy*roche)&&(cloop>10)){	//avoid alternating
      if (debug) cerr << "roche=" << roche << "Rsun rocheold=" << rocheold << "Rsun rocheolder=" << rocheolder << "Rsun |rocheolder-roche|=" << fabs(rocheolder-roche) << endl;
      if ((star.track.r[star.rmax]<roche)&&((rratio>1)||(rratio<0))){	//check if maximal stellar radius can fill Roche lobe
        return 1.0e+99;	//no Roche-lobe overflow
      }
      for (j=star.last+1;j<star.rmax;j++){	//find time index of Roche-lobe overflow
        if ((star.track.r[j]>roche)&&(star.track.r[j]>star.track.r[j-1])&&(star.track.t[j]>=star.t[n])){
          break;
        }
      }
      return RLFT2(star,companion,n,min(j,jr)-1,(j==jr),acirc*system.M[n],r);
    }
    if ((cloop>1000)&&(cloop%10000>100)) roche = rocheold+1000.0/cloop*(roche-rocheold);	//take value between roche and rocheold	//take mean of roche and rocheold to get closer to fix point
    if (cloop>100000){	//if loop is getting endless
      if ((star.track.r[star.rmax]<roche)&&((rratio>1)||(rratio<0))){	//check if maximal stellar radius can fill Roche lobe
        return 1.0e+99;	//no Roche-lobe overflow
      }
      break;	//avoid alternating
    }
    cloop++;
    told = t;
    rocheolder = rocheold;
    rocheold = roche;
    for (j=star.last+1;j<star.rmax;j++){	//find time index of Roche-lobe overflow
      if ((star.track.r[j]>roche)&&(star.track.r[j]>star.track.r[j-1])&&(star.track.t[j]>=star.t[n])){
        if (cloop<2000) break;
        if (j>jr-2) break;
        else if (cjump<5){
          cjump++;
          break;
        }	//else if (cjump==5) cjump++;
      }
    }
    jr = j;
    if (star.track.r[j]==star.track.r[j-1]){
      rratio *= ratio(star.track.r[j],rocherad(acirc*(system.M[n])/(companion.m[n]+star.track.m[j]),star.track.m[j]/companion.m[n]),roche);	//calculate ratio in radius
      if (debug) cerr << "star.track.r[" << j << "]=" << star.track.r[j] << "Rsun rocherad(" << n << ")=" << rocherad(acirc*(system.M[n])/(companion.m[n]+star.track.m[j]),star.track.m[j]/companion.m[n]) << "Rsun roche=" << roche << "Rsun" << endl;
      constantradius = true;
    }else{
      rratio = ratio(roche,star.track.r[j],star.track.r[j-1]);	//calculate ratio in radius
    }
    if (isnan(rratio)){
      cerr << "#Error: r_ratio=" << rratio << " roche=" << roche << "Rsun star.track.r[j]=" << star.track.r[j] << "Rsun star.track.r[j-1]=" << star.track.r[j-1] << "Rsun j=" << j << endl;
      screen = true;
    }
    if ((rratio<0)&&(j==star.rmax)){	//star does not reach Roche lobe any longer
//&&(((roche!=system.rp[n])&&(roche!=system.rs[n]))||(roche==rochemin))
//      if (rocherad(acirc*(system.M[n])/(companion.m[n]+star.track.m[j]),star.track.m[j]/companion.m[n])>star.track.r[j])
      if (debug) cerr << "rratio=" << rratio << " j=star.rmax=" << star.rmax << endl;
      if (extrapolation) return 1.0e+99;	//no Roche-lobe overflow
//      else if (companion.stage[n]==0) rratio = fmax(0.0,ratio(companion.track.t[companion.track.TAMS],star.track.t[j],star.track.t[j-1]));
      else rratio = fmax(0.0,ratio(tc,star.track.t[j],star.track.t[j-1]));
      extrapolation = true;
      if (debug){
        cerr << "rratio_new=" << rratio;
        if (companion.stage[n]==0) cerr << " rratio(companion.track.t[companion.track.TAMS])=" << ratio(companion.track.t[companion.track.TAMS],star.track.t[j],star.track.t[j-1]) << endl;
        else cerr << " rratio(companion.track.t[companion.track.n-1])=" << ratio(companion.track.t[companion.track.n-1],star.track.t[j],star.track.t[j-1]) << endl;
      }
    }
    nstarmass = star.track.m[j]-rratio*(star.track.m[j]-star.track.m[j-1]);	//calculate mass when radius = roche
    if ((roche>star.track.r[j])&&(rochemin>star.track.r[j])){	//star does not reach Roche lobe any longer
      cerr << "#Error: roche=" << roche << "Rsun star.track.r[" << j << "]=" << star.track.r[j] << "Rsun star.track.r[" << star.rmax << "]=" << star.track.r[star.rmax] << "Rsun" << endl;
      return 1.0e+99;	//no Roche-lobe overflow
    }
    if (nstarmass<1.0e-10){
      cerr << "#Error: new star mass=" << nstarmass << "Msun star.m[" << n << "]=" << star.m[n] << "Msun star.track.m[j]=" << star.track.m[j] << "Msun star.track.m[j-1]=" << star.track.m[j-1] << "Msun j=" << j << " last=" << star.last << " star.rmax=" << star.rmax << " star.track.n-1=" << star.track.n-1 << " r_ratio=" << rratio << " roche=" << roche << "Rsun star.track.r[j]=" << star.track.r[j] << "Rsun star.track.r[j-1]=" << star.track.r[j-1] << "Rsun" << endl;
      screen = true;
    }
//    rad = star.track.r[j]-rratio*(star.track.r[j]-star.track.r[j-1]);
    t = star.track.t[j]-rratio*(star.track.t[j]-star.track.t[j-1]);	//calculate time when radius = roche
//    t = star.track.t[j-1]+(1.0-rratio)*(star.track.t[j]-star.track.t[j-1]);	//calculate time when radius = roche
    if (t<star.t[n]){
      if (screen) cerr << "#Warning: t=" << t << "yr star.t[" << n << "]=" << star.t[n] << "yr star.track.t[j]=" << star.track.t[j] << "yr star.track.t[j-1]=" << star.track.t[j-1] << "yr j=" << j << " last=" << star.last << " star.rmax=" << star.rmax << " star.track.n-1=" << star.track.n-1 << " r_ratio=" << rratio << " roche=" << roche << "Rsun star.track.r[j]=" << star.track.r[j] << "Rsun star.track.r[j-1]=" << star.track.r[j-1] << "Rsun cloop=" << cloop << " cjump=" << cjump << endl;
//      return RLFT2(star,companion,n,star.last,false,acirc*system.M[n],r);
    }
    if (debug) cerr << "j(star)=" << j << " rratio=" << rratio << " star.track.r[j-1]=" << star.track.r[j-1] << "Rsun star.track.r[j]=" << star.track.r[j] << "Rsun t=" << t << "yr nstarmass=" << nstarmass << "Msun";
    if (t>ts){	//star does not reach Roche lobe any longer
      cerr << endl << "#Error: t=" << t << "yr star.track.t[" << j << "]=" << star.track.t[j] << "yr star.track.t[" << j-1 << "]=" << star.track.t[j-1] << "yr ts=" << ts << "yr r_ratio=" << rratio << endl;
      return 1.0e+99;	//no Roche-lobe overflow
    }
    for (j=companion.last+1;j<companion.track.n-1;j++){	//find time index of t
      if (companion.track.t[j]>t){
        break;
      }
    }
    tratio = ratio(t,companion.track.t[j],companion.track.t[j-1]);	//calculate ratio in time
    if (t>tc){	//companion does not reach time t
      if (debug) cerr << " ts=" << ts << "yr t=" << t << "yr tc=" << tc << "yr t-tc=" << t-tc << "yr t-ts=" << t-ts << "yr tratio=" << tratio << endl;
      ncompanionmass = companion.track.m[jmax];
      for (j=star.last+1;j<star.rmax;j++){	//find time index of maximal time of the companion
        if (star.track.t[j]>tc) break;
      }
      tratio = ratio(tc,star.track.t[j],star.track.t[j-1]);	//calculate ratio in time
      nstarmass = star.track.m[j]-tratio*(star.track.m[j]-star.track.m[j-1]);	//calculate mass of the star at maximal time of the companion
      na = acirc*(system.M[n])/(nstarmass+ncompanionmass);	//calculate new semi-major-axis
      nq = nstarmass/ncompanionmass;	//calculate new mass ratio
      roche = rocherad(na,nq);	//calculate new Roche lobe
      if (debug) cerr << "star.rmax=" << star.rmax << " star.track.t[" << j << "]=" << star.track.t[j] << "yr tc=" << tc << "yr companion.track.m[" << jmax << "]=" << companion.track.m[jmax] << "Msun tratio=" << tratio << " nstarmass=" << star.track.m[j]-tratio*(star.track.m[j]-star.track.m[j-1]) << "Msun roche=" << roche << "Rsun r_star=" << star.track.r[j]-tratio*(star.track.r[j]-star.track.r[j-1]) << "Rsun" << endl;
      if (roche<star.track.r[j]-tratio*(star.track.r[j]-star.track.r[j-1])) return RLFT2(star,companion,n,star.last,false,acirc*system.M[n],r);
      else return 1.0e+99;	//no Roche-lobe overflow
    }
    if (isnan(tratio)){
      cerr << endl << "#Error: t_ratio=" << tratio << " t=" << t << "yr companion.track.t[j]=" << companion.track.t[j] << "yr companion.track.t[j-1]=" << companion.track.t[j-1] << "yr j=" << j << endl;
      screen = true;
    }
    ncompanionmass = companion.track.m[j]-tratio*(companion.track.m[j]-companion.track.m[j-1]);	//calculate mass at time
    if (ncompanionmass<=0){
      cerr << endl << "#Error: new companion mass=" << ncompanionmass << "Msun companion.m[" << n << "]=" << companion.m[n] << "Msun companion.track.m[j]=" << companion.track.m[j] << "Msun companion.track.m[j-1]=" << companion.track.m[j-1] << "Msun j=" << j << " last=" << companion.last << " companion.track.n-1=" << companion.track.n-1 << " t_ratio=" << tratio << " t=" << t << "yr companion.track.t[j]=" << companion.track.t[j] << "yr companion.track.t[j-1]=" << companion.track.t[j-1] << "yr" << endl;
      screen = true;
    }
    na = acirc*(system.M[n])/(nstarmass+ncompanionmass);	//calculate new semi-major-axis
//    nperi = system.peri[n]*(system.M[n])/(nstarmass+ncompanionmass);	//calculate new periastron
    nq = nstarmass/ncompanionmass;	//calculate new mass ratio
    roche = rocherad(na,nq);	//calculate new Roche lobe
    if (isnan(roche)){
      cerr << endl << "#Error: Roche lobe=" << roche << "Rsun na=" << na << "Rsun new star mass=" << nstarmass << "Msun new companion mass=" << ncompanionmass << "Msun nq=" << nq << " r_ratio=" << rratio << " t_ratio=" << tratio << endl;
// nperi=" << nperi << "Rsun
      screen = true;
    }
    if (debug) cerr << " j(comp)=" << j << " tratio=" << tratio << " companion.track.t[j-1]=" << companion.track.t[j-1] << "yr companion.track.t[j]=" << companion.track.t[j] << "yr ncompanionmass=" << ncompanionmass << "Msun na=" << na << "Rsun nq=" << nq << " roche=" << roche << "Rsun rocheold=" << rocheold << "Rsun" << endl;
  }while ((fabs(t-told)>1.0)||(fabs(rocheold-roche)>0.0001)||(rratio<0.0)||(rratio>1.0)||((!constantradius)&&((roche<fmin(star.track.r[jr-1],star.track.r[jr]))||(roche>fmax(star.track.r[jr-1],star.track.r[jr]))))||((constantradius)&&(fabs(roche-star.track.r[jr])>0.0001)));	//iterate until time is depermined to 1yr precision and Roche Lobe is depermined to 0.0001Rsun
//  if ((fabs(roche-rad)>0.01)||(rratio<0.0)||(rratio>1.0)||(tratio<0.0)||(tratio>1.0)) cout << "roche=" << roche << "Rsun radius=" << rad << "Rsun t=" << t << "yr rratio=" << rratio << " tratio=" << tratio << endl;
  if (debug) cerr << "t-told=" << t-told << "yr rocheold-roche=" << rocheold-roche << "Rsun roche=" << roche << "Rsun" << endl;
  if (t<star.t[n]){	//check that RLO is not in the past
    r = star.r[n];
//    return star.t[n];
    return 1.0e+99;	//no Roche-lobe overflow
  }else{
    if (constantradius) r = -roche;
    else r = roche;
    return t;
  }
}

double RLFT2(t_star star, t_star companion, int n, int imin, bool same, double acircM, double& r){
  int i,j;	//time index of the star and the companion
  int iold,jold;	//old time index of the star and the companion
  int jmax;	//maximal index of companion
  double tratio;	//ratio between the two nearest times
  double t;	//time in yr
  double tlow, tup;	//time range in yr
  double ms,mc;	//star and companion mass in Msun
  double rs,roche;	//stars radius and Roche lobe in Rsun
  double diff;

  if (companion.stage[n]==0) jmax = companion.track.TAMS;	//maximal index is at TAMS
  else jmax = companion.track.n-1;	//maximal index is at the end of the track
  i = imin;
  t = fmax(star.track.t[i],star.t[n]);	//set time to time of the star
  tratio = ratio(t,star.track.t[i+1],star.track.t[i]);	//calculate ratio in time
  ms = star.track.m[i+1]-tratio*(star.track.m[i+1]-star.track.m[i]);	//set mass of the star
  rs = star.track.r[i+1]-tratio*(star.track.r[i+1]-star.track.r[i]);	//set radius of the star
  for (j=companion.last+1;j<jmax;j++){	//find time index of companion
    if (companion.track.t[j]>t){
      break;
    }
  }
  tratio = ratio(t,companion.track.t[j],companion.track.t[j-1]);	//calculate ratio in time
  mc = companion.track.m[j]-tratio*(companion.track.m[j]-companion.track.m[j-1]);	//calculate mass of the companion at time
  roche = rocherad(acircM/(ms+mc),ms/mc);	//calculate Roche lobe
  diff = rs-roche;	//calculate the difference between real and expacted radius
  j--;
  if (debug) cerr << "i=" << i << " j=" << j << " t=" << t << "yr star.track.t[i]=" << star.track.t[i] << "yr star.t[" << n << "]=" << star.t[n] << "yr tratio=" << tratio << " m_star=" << ms << "Msun m_comp=" << mc << "Msun r_star=" << rs << "Rsun roche=" << roche << "Rsun diff=" << diff << "Rsun" << endl;
  do{	//search for point, where radius=roche lobe
    iold = i;
    jold = j;
    tlow = t;
    diff = rs-roche;	//calculate the difference between real and expacted radius
    if (star.track.t[min(i+1,star.rmax)]<companion.track.t[min(j+1,jmax)]){
      i++;	//increase star index
      if (i>star.rmax){
        cerr << "#Error: i=" << i << ">star.rmax=" << star.rmax << endl;
        screen = true;
        return 1.0e+99;	//no Roche-lobe overflow
      }
      t = fmax(star.track.t[i],star.t[n]);	//set time to time of the star
      ms = star.track.m[i];	//set mass of the star
      rs = star.track.r[i];	//set radius of the star
      tratio = ratio(t,companion.track.t[j+1],companion.track.t[j]);	//calculate ratio in time
      mc = companion.track.m[j+1]-tratio*(companion.track.m[j+1]-companion.track.m[j]);	//calculate mass of the companion at time
      if (same){
        roche = rocherad(acircM/(ms+mc),ms/mc);	//calculate Roche lobe
        if (diff*(rs-roche)>0) return 1.0e+99;	//no Roche-lobe overflow
        break;
      }
    }else{
      j++;	//increase companion index
      if (j>jmax){
        cerr << "#Error: j=" << j << ">jmax=" << jmax << endl;
        screen = true;
        return 1.0e+99;	//no Roche-lobe overflow
      }
      t = companion.track.t[j];	//set time to time of the companion
      mc = companion.track.m[j];	//set mass of the companion
      tratio = ratio(t,star.track.t[i+1],star.track.t[i]);	//calculate ratio in time
      ms = star.track.m[i+1]-tratio*(star.track.m[i+1]-star.track.m[i]);	//calculate mass of the star at time
      rs = star.track.r[i+1]-tratio*(star.track.r[i+1]-star.track.r[i]);	//calculate radius of the star at time
    }
    roche = rocherad(acircM/(ms+mc),ms/mc);	//calculate Roche lobe
    if (debug) cerr << "i=" << i << " j=" << j << " t=" << t << "yr m_star=" << ms << "Msun m_comp=" << mc << "Msun r_star=" << rs << "Rsun roche=" << roche << "Rsun diff=" << diff << "Rsun" << endl;
  }while(diff*(rs-roche)>0);	//check for zero in the difference
  if (debug) cerr << "i=" << i << " j=" << j << " t=" << t << "yr star.track.t[i]=" << star.track.t[i] << "yr star.t[" << n << "]=" << star.t[n] << "yr tratio=" << tratio << " m_star=" << ms << "Msun m_comp=" << mc << "Msun r_star=" << rs << "Rsun roche=" << roche << "Rsun diff=" << diff << "Rsun" << endl;
  tup = t;
  t = 0.5*(tlow+tup);
  do{
    tratio = ratio(t,star.track.t[iold+1],star.track.t[iold]);	//calculate ratio in time
    ms = star.track.m[iold+1]-tratio*(star.track.m[iold+1]-star.track.m[iold]);	//calculate mass of the star at time
    rs = star.track.r[iold+1]-tratio*(star.track.r[iold+1]-star.track.r[iold]);	//calculate radius of the star at time
    tratio = ratio(t,companion.track.t[jold+1],companion.track.t[jold]);	//calculate ratio in time
    mc = companion.track.m[jold+1]-tratio*(companion.track.m[jold+1]-companion.track.m[jold]);	//calculate mass of the companion at time
    roche = rocherad(acircM/(ms+mc),ms/mc);	//calculate Roche lobe
    if (debug) cerr << "tlow=" << tlow << "yr tup=" << tup << "yr t=" << t << "yr m_star=" << ms << "Msun m_comp=" << mc << "Msun r_star=" << rs << "Rsun roche=" << roche << "Rsun diff=" << diff << "Rsun" << endl;
    if (rs==roche) break;
    else if (diff*(rs-roche)>0){
      diff = rs-roche;
      tlow = t;
    }else tup = t;
    t = 0.5*(tlow+tup);
    if (tup==tlow){
      screen = true;
      if (screen) cerr << "#Warning: no solution in RLFT2" << endl;
      return 1.0e+99;	//no Roche-lobe overflow
    }else if ((t==tlow)||(tup==t)){
      screen = true;
      if (screen) cerr << "#Warning: accuracy not reached in RLFT2: rs=" << rs << "Rsun roche=" << roche << "Rsun diff=" << diff << "Rsun 1.0-roche/rs=" << 1.0-roche/rs << "=" << diff/rs << endl;
      if (debug) cerr << "          r(tlow=" << tlow << "yr)=" << star.track.r[iold+1]-ratio(tlow,star.track.t[iold+1],star.track.t[iold])*(star.track.r[iold+1]-star.track.r[iold])
                      << "Rsun r(t=" << t << "yr)=" << star.track.r[iold+1]-ratio(t,star.track.t[iold+1],star.track.t[iold])*(star.track.r[iold+1]-star.track.r[iold])
                      << "Rsun r(tup=" << tup << "yr)=" << star.track.r[iold+1]-ratio(tup,star.track.t[iold+1],star.track.t[iold])*(star.track.r[iold+1]-star.track.r[iold]) << "Rsun" << endl;
      break;
    }
  }while((tup-tlow>0.1)||(fabs(1.0-roche/rs)>accuracy));
  r = roche;
  return t;
}

double RLOAtime(double Md, double Ma, double P, double Dt){  //getting the duration of the RLO Case A from donor mass, accretor mass and period, based on fitting by Joana Kramer
  //~ double q = Ma/Md;
  //~ P = log10(P*365.25);
  //~ cout << Dt + pow(Dt,1.3667)*pow(10.,-2.7887) << endl;
  //~ return Dt + pow(Dt,1.3667)*pow(10.,-2.7887);
  
  Md = log10(Md);
  P = log10(P*365.25);
  if(Md > 1.6) Md = 1.6;
  if(Md < 1.0) Md = 1.0;
  if(P > 0.619 + 1.868*square(Md-0.957)) P = 0.619 + 1.868*square(Md-0.957);
  return pow(10.0, 13.457 - Md*10.128 - P*6.245 + Md*P*19.341 + Md*Md*2.856 - P*P*19.227 - Md*Md*P*11.647 + Md*P*P*15.273 + Md*Md*Md*0.757 - P*P*P*2.175);
}

double RLOAmass(double Md, double Ma, double P){  //getting the donor mass after the RLO Case A from donor mass, accretor mass and period, based on fitting by Joana Kramer
  double q = Ma/Md;
  Md = log10(Md);
  P = log10(P*365.25);
  //cout << pow(10.0, -1.7 - 0.57*Md + 3.95*P + 1.6*q + 2.3*Md*Md - 1.25*P*P - 0.98*q*q - 2.66*Md*P - 0.2*Md*q - 1.5*P*q - 0.65*Md*Md*Md + 0.246*P*P*P + 0.20*q*q*q + 0.60*Md*Md*P - 0.2*Md*Md*q + 0.20*P*P*Md + 0.24*P*P*q + 0.3*q*q*Md + 0.08*q*q*P + 0.53*Md*P*q) << endl;
  return pow(10.0, -1.7 - 0.57*Md + 3.95*P + 1.6*q + 2.3*Md*Md - 1.25*P*P - 0.98*q*q - 2.66*Md*P - 0.2*Md*q - 1.5*P*q - 0.65*Md*Md*Md + 0.246*P*P*P + 0.20*q*q*q + 0.60*Md*Md*P - 0.2*Md*Md*q + 0.20*P*P*Md + 0.24*P*P*q + 0.3*q*q*Md + 0.08*q*q*P + 0.53*Md*P*q);
}

double Schwarzschildradius(double M){	//calculate radius of a Black Hole
  return 2.0*G*M/(c*c);	//remnant radius = 2.0*G*M/c^2
}

double NSradius(double M){	//calculate radius of a Neutron Star
  return 1.0E+4/Rsun;	//remnant radius = 10km
}

double WDradius(double M){	//calculate radius of a White Dwarf
  return 0.013/cbrt(M);	//non-relativistic WD
}
