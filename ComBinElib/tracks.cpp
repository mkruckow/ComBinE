//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
#include <stdio.h>
#include "ComBinElib.h" // local header file

//using namespace std;

void updatetrack(t_star& star, int n){
  if (star.stage[n]==0) updatestarstrack(star, n, stararray1, nstar1, false);	//hydrogen burning star
  else if (star.stage[n]==1) updatestarstrack(star, n, stararray2, nstar2, false);	//core helium burning star
  else if (star.stage[n]==2) updatestarstrack(star, n, hestararray1, nhe1, false);	//naked helium star
//  else if (star.stage[n]==2) updateHetrack(star,n);	//naked helium star
  else if (star.stage[n]==-3) return;	//get remnant
  else if ((star.stage[n]>6)||(star.stage[n]<3)) cerr << "#Error: no track update: star.stage[n]=" << star.stage[n] << endl;	//is not remnant and a star
}

void updatestarstrack(t_star& star, int n, t_HRD* stararray0, int nstar0, bool SN){
  int i;		//index variable
  int ilow, iup;	//index of lower/upper track
  int jlow=1, jup=1;	//index in lower/upper track
  int jmax=1, jmaxlow;	//maximal value for jup/jlow
  int nlow, nup;	//number of track points in lower/upper track
  int cnosol=0;		//count no solutions
  double mratio;	//ratio in track of the wanted mass
  double cmup=0.0;	//core mass at lower/upper track at the wanted mass	//cmlow=0.0, 
  long double j11=0,j12,j21=0,j22;	//j-position at jlow-1, jlow, jup-1, jup	//j-position=normalized position of two neaburing track
  double maxr=0.0;	//maximal radius
  long double A=0.0,B=0.0,C=0.0;	//auxiliary variables
  long double r=0.0,r1=0.0,r2=0.0;	//auxiliary variables: ratios
  t_HRD newtrack;	//new track part
  bool trackadaption=false;	//if the track is adapted
  long double totalup=0.0, totallow=0.0;	//total length in upper/lower track
  long double clup=0.0, cllow=0.0, clupold=0.0, cllowold=0.0;	//current length in upper/lower track

  if (debug) cerr << endl << "search for: star.m[" << n << "]=" << star.m[n] << "Msun and star.cm[" << n << "]=" << star.cm[n] << "Msun" <<endl;
  //find grid position of first track guess
  for (iup=1;iup<nstar0;iup++){
    if (stararray0[iup].m[0]>star.m[n]){	//find track with the first possible solution
      break;
    }
  }
  if (iup>=nstar0){	//too large mass
//    goto notrackupdate;
    iup = nstar0-1;
    cerr << "#Warning: reduce mass from star.m[" << n << "]=" << star.m[n] << "Msun to " << stararray0[iup].m[0] << "Msun" << endl;
    star.m[n] = stararray0[iup].m[0];	//reduce mass
    trackadaption = true;
/*    cnotrackupdate++;
    if (screen) cerr << "#Warning: no track update, because star out of grid: star.m[" << n << "]=" << star.m[n] << "Msun and star.cm[" << n << "]=" << star.cm[n] << "Msun" << endl;
      undophase(star, n, true);
    return;*/
  }
  if (debug) cerr << "star.m[" << n << "]=" << star.m[n] << "Msun in stararray0[" << iup-1 << "].m[0]=" << stararray0[iup-1].m[0] << "Msun stararray0[" << iup << "].m[0]=" << stararray0[iup].m[0] << "Msun" << endl;
  mratio = ratio(star.m[n],stararray0[iup].m[0],stararray0[iup-1].m[0]);	//calculate ratio between first track points
  cmup = stararray0[iup].cm[0]-mratio*(stararray0[iup].cm[0]-stararray0[iup-1].cm[0]);	//calculate lowest possible core mass between first track points
  if (debug) cerr << "at first track point: cmup=" << cmup << "Msun star.cm[" << n << "]-cmup=" << star.cm[n]-cmup << "Msun mratio=" << mratio << endl;
  if (star.cm[n]<=cmup){	//no significant core
    ilow = iup-1;	//set ilow to privious track of iup
    mratio = ratio(star.m[n],stararray0[iup].m[0],stararray0[ilow].m[0]);
    if (debug) cerr << "star.cm[n]=" << star.cm[n] << "Msun<=cmup=" << cmup << "Msun mratio=" << mratio << endl;
    nlow = stararray0[ilow].n-1;	//get number of track points at ilow
    nup = stararray0[iup].n-1;	//get number of track points at iup
//    goto goalongtrack;
  }else{
    iup--;	//go one step back
    //find nearest tracks
    while (star.cm[n]>cmup){
      if (iup<nstar0-1){	//while search is in the grid
        iup++;	//go to next track
      }else{	//leaves grid: too large mass
/*        cerr << "#Error: core mass not reached: star.cm[" << n << "]=" << star.cm[n] << "Msun cmlow=" << cmlow << "Msun cmup=" << cmup << "Msun" << endl << "#Error: core mass not reached: cmlow=" << cmlow << "Msun stararray0[" << iup-1 << "].cm[" << jlow-1 << "]=" << stararray0[iup-1].cm[jlow-1] << "Msun stararray0[" << iup-1 << "].cm[" << jlow << "]=" << stararray0[iup-1].cm[jlow] << "Msun" << endl << "#Error: core mass not reached: cmup=" << cmup << "Msun stararray0[" << iup << "].cm[" << jup-1 << "]=" << stararray0[iup].cm[jup-1] << "Msun stararray0[" << iup << "].cm[" << jup << "]=" << stararray0[iup].cm[jup] << "Msun mratio=" << mratio << endl;
        screen = true;*/
        if (debug) cerr << "iup=" << iup << ">=nstar0-1=" << nstar0-1 << endl;
//        goto notrackupdate;
        iup = nstar0-1;
        //find position of core mass in track iup
        for (jup=1;jup<stararray0[iup].n;jup++){
          if (stararray0[iup].cm[jup]>star.cm[n]){
            break;
          }
        }
        if (jup>=stararray0[iup].n){	//too large core mass
          //find position of maximum core mass in track iup
          for (jup=1;jup<stararray0[iup].n-1;jup++){
            if (stararray0[iup].cm[jup+1]<stararray0[iup].cm[jup]){
              break;
            }
          }
          cerr << "#Warning: reduce core mass from star.cm[" << n << "]=" << star.cm[n] << "Msun to " << stararray0[iup].cm[jup] << "Msun" << endl;
          star.cm[n] = stararray0[iup].cm[jup];	//reduce core mass
        }
        mratio = ratio(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]);	//get coremassratio
        cmup = (stararray0[iup].m[jup-1]-mratio*(stararray0[iup].m[jup-1]-stararray0[iup].m[jup]))*(1.0-accuracy);	//calculate mass at core mass position at track iup
        cerr << "#Warning: reduce mass from star.m[" << n << "]=" << star.m[n] << "Msun to " << cmup << "Msun" << endl;
        star.m[n] = cmup;	//reduce mass
        trackadaption = true;
//        screen = true;
/*        cnotrackupdate++;
        if (screen) cerr << "#Warning: no track update, because star out of grid: star.m[" << n << "]=" << star.m[n] << "Msun and star.cm[" << n << "]=" << star.cm[n] << "Msun" << endl;
          undophase(star, n, true);
        return;*/
      }
      //save values of privious track
      jlow = jup;
//      cmlow = cmup;
      //find position in track iup
      for (jup=1;jup<stararray0[iup].n;jup++){
        if ((stararray0[iup].m[jup]<star.m[n])||(stararray0[iup].cm[jup]>star.cm[n])){
          break;
        }
      }
      if ((stararray0[iup].cm[jup]>star.cm[n])&&(stararray0[iup].m[jup]>star.m[n])){	//core mass reached at track iup
        mratio = 1.0;	//set massratio
        cmup = stararray0[iup].cm[jup-1]-mratio*(stararray0[iup].cm[jup-1]-stararray0[iup].cm[jup]);	//calculate core mass at mass position at track iup
        ilow = iup-1;
        for (jlow=1;jlow<stararray0[ilow].n;jlow++){
          if (stararray0[ilow].m[jlow]<star.m[n]){
            break;
          }
        }
      }else if (jup<stararray0[iup].n){	//mass reached at track iup
        mratio = ratio(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup]);	//get massratio
        if ((mratio<0)||(mratio>1)){
          cerr << "#Error: mratio out of range: mratio=" << mratio << " star.m[" << n << "]=" << star.m[n] << "Msun stararray0[" << iup << "].m[" << jup-1 << "]=" << stararray0[iup].m[jup-1] << "Msun stararray0[" << iup << "].m[" << jup << "]=" << stararray0[iup].m[jup] << "Msun" << endl;	//check if ratio is in [0,1]
          screen = true;
        }
        cmup = stararray0[iup].cm[jup-1]-mratio*(stararray0[iup].cm[jup-1]-stararray0[iup].cm[jup]);	//calculate core mass at mass position at track iup
      }else{	//mass not reached at track iup
        jup = stararray0[iup].n-1;
        break;
      }
      if (debug) cerr << "star.m[" << n << "]=" << star.m[n] << "Msun in stararray0[" << iup << "].m[" << jup-1 << "]=" << stararray0[iup].m[jup-1] << "Msun stararray0[" << iup << "].m[" << jup << "]=" << stararray0[iup].m[jup] << "Msun cmup=" << cmup << "Msun mratio=" << mratio << endl;
    }
//    goalongtrack:	//goto-marker
    ilow = iup-1;	//set ilow to privious track of iup
    nexttrack:	//goto-marker
    if (debug) cerr << "star.m[" << n << "]=" << star.m[n] << "Msun in stararray0[" << ilow << "].m[" << jlow-1 << "]=" << stararray0[ilow].m[jlow-1] << "Msun stararray0[" << ilow << "].m[" << jlow << "]=" << stararray0[ilow].m[jlow] << "Msun stararray0[" << iup << "].m[" << jup-1 << "]=" << stararray0[iup].m[jup-1] << "Msun stararray0[" << iup << "].m[" << jup << "]=" << stararray0[iup].m[jup] << "Msun" << endl << "star.cm[" << n << "]=" << star.cm[n] << "Msun in stararray0[" << ilow << "].cm[" << jlow-1 << "]=" << stararray0[ilow].cm[jlow-1] << "Msun stararray0[" << ilow << "].cm[" << jlow << "]=" << stararray0[ilow].cm[jlow] << "Msun stararray0[" << iup << "].cm[" << jup-1 << "]=" << stararray0[iup].cm[jup-1] << "Msun stararray0[" << iup << "].cm[" << jup << "]=" << stararray0[iup].cm[jup] << "Msun" << endl;

    nlow = stararray0[ilow].n-1;	//get number of track points at ilow
    nup = stararray0[iup].n-1;	//get number of track points at iup
    jmaxlow = nlow;	//save maximum of jlow
    jmax = nup;	//save maximum of jup
    if (jlow>1) jlow--;	//go one step back
    if (jup>1) jup--;	//go one step back
    totallow = stararray0[ilow].t[nlow];	//get total length of ilow
    totalup = stararray0[iup].t[nup];	//get total length of iup
    cllow = stararray0[ilow].t[jlow];	//get current length at jlow in ilow
    cllowold = stararray0[ilow].t[jlow-1];	//get previous length at jlow in ilow
    clup = stararray0[iup].t[jup];	//get current length at jup in iup
    clupold = stararray0[iup].t[jup-1];	//get previous length at jup in iup
    j11 = clupold*totallow;	//get lower j at iup so that j at ilow and iup are the same at the beginning and the end of the tracks
    j12 = clup*totallow;	//get upper j at iup so that j at ilow and iup are the same at the beginning and the end of the tracks
    j21 = cllowold*totalup;	//get lower j at ilow so that j at ilow and iup are the same at the beginning and the end of the tracks
    j22 = cllow*totalup;	//get upper j at ilow so that j at ilow and iup are the same at the beginning and the end of the tracks
    if (debug) cerr << "findj: j11 = " << j11 << " j12 = " << j12 << " j21 = " << j21 << " j22 = " << j22 << " totallow = " << totallow << " totalup = " << totalup << " cllowold = " << cllowold << " cllow = " << cllow << " clupold = " << clupold << " clup = " << clup << " jlow = " << jlow << " jup = " << jup << " jmaxup = " << jmax << " jmaxlow = " << jmaxlow << endl;
    if(j11<j21){	//j11<j21
      for (jlow=1;jlow<nlow;jlow++){	//find relative length of clupold in ilow to get the starting jlow
        if (stararray0[ilow].t[jlow]*totalup>=j11) break;
      }
    }else{	//j11>=j21
      for (jup=1;jup<nup;jup++){	//find relative length of cllowold in iup to get the starting jup
        if (stararray0[iup].t[jup]*totallow>=j21) break;
      }
    }
/*    if(jlow>1) jlow--;
    cllowold = stararray0[ilow].t[jlow-1];	//get previous length at jlow in ilow
    cerr << "cllow=" << cllowold*totalup << "yr^2" << endl;
    for (jup=1;jup<nup;jup++){	//find relative length of cllowold in iup to get the starting jup
      cerr << "t(i=" << iup << ",j=" << jup+1 << ")=" << stararray0[iup].t[jup+1] << "yr -> clup=" << stararray0[iup].t[jup+1]*totallow << "yr^2" << endl;
      if (stararray0[iup].t[jup+1]*totallow>cllowold*totalup) break;
    }
    clup = stararray0[iup].t[min(jup+1,nup)];	//get next length at jup in iup
    clupold = stararray0[iup].t[min(jup,nup-1)];	//get current length at jup in iup
    if (jup>=jmax){
      jup = jmax;
      for (jmax=1;jmax<nup;jmax++){	//find relative length of cllow in iup to get the jmax
        if (stararray0[iup].t[jmax]*totallow>cllow*totalup) break;
      }
      clup = stararray0[iup].t[min(jup+1,nup)];	//get next length at jup in iup
      clupold = stararray0[iup].t[min(jup,nup-1)];	//get current length at jup in iup
      for (jlow=1;jlow<nlow;jlow++){	//find relative length of clupold in ilow to get the starting jlow
        if (stararray0[ilow].t[jlow+1]*totalup>clupold*totallow) break;
      }
    }*/
    //first digit of j??, c?? and m??: 1=up, 2=low
    //second digit of j??, c?? and m??: 1=jup/jlow-1, 2=jup/jlow
    cllowold = stararray0[ilow].t[jlow-1];	//get previous length at jlow in ilow
    j21 = cllowold*totalup;	//get lower j at ilow so that j at ilow and iup are the same at the beginning and the end of the tracks
    cllow = stararray0[ilow].t[jlow];	//get current length at jlow in ilow
    j22 = cllow*totalup;	//get upper j at ilow so that j at ilow and iup are the same at the beginning and the end of the tracks
    clupold = stararray0[iup].t[jup-1];	//get previous length at jup in iup
    j11 = clupold*totallow;	//get lower j at iup so that j at ilow and iup are the same at the beginning and the end of the tracks
    clup = stararray0[iup].t[jup];	//get current length at jup in iup
    j12 = clup*totallow;	//get upper j at iup so that j at ilow and iup are the same at the beginning and the end of the tracks
    if (debug) cerr << "guess: j11 = " << j11 << " j12 = " << j12 << " j21 = " << j21 << " j22 = " << j22 << " totallow = " << totallow << " totalup = " << totalup << " cllowold = " << cllowold << " cllow = " << cllow << " clupold = " << clupold << " clup = " << clup << " jlow = " << jlow << " jup = " << jup << " jmaxup = " << jmax << " jmaxlow = " << jmaxlow << endl;
    //find first position which includes the wanted core mass
    while (fmax(fmax(stararray0[iup ].cm[jup -1],stararray0[iup ].cm[jup   ]),fmax(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow  ]))<star.cm[n]){
      if (j22<j12){	//go to next point at ilow
        jlow++;
        cllowold = cllow;
        cllow = stararray0[ilow].t[jlow];
        j21 = cllowold*totalup;
        j22 = cllow*totalup;
      }else if (j22>j12){	//go to next point at iup
        jup++;
        clupold = clup;
        clup = stararray0[iup].t[jup];
        j11 = clupold*totallow;
        j12 = clup*totallow;
      }else{	//go to next point at ilow and iup
        jlow++;
        jup++;
        cllowold = cllow;
        cllow = stararray0[ilow].t[jlow];
        clupold = clup;
        clup = stararray0[iup].t[jup];
        j21 = cllowold*totalup;
        j22 = cllow*totalup;
        j11 = clupold*totallow;
        j12 = clup*totallow;
      }
      if (jlow>nlow){
        jlow = nlow;
        if (jup>nup) jup = nup;
        break;
      }
      if (jup>nup){
        jup = nup;
        break;
      }
    }
    if (debug) cerr << "start: j11 = " << j11 << " j12 = " << j12 << " j21 = " << j21 << " j22 = " << j22 << " totallow = " << totallow << " totalup = " << totalup << " cllowold = " << cllowold << " cllow = " << cllow << " clupold = " << clupold << " clup = " << clup << " jlow = " << jlow << " jup = " << jup << " jmaxup = " << jmax << " jmaxlow = " << jmaxlow << endl;
    mratio = -1;	//set undefined value for the ratio
    //calculate solutions for the ratio between the two nearest tracks
    while ((jup<=jmax)&&(jlow<=jmaxlow)){
      if (debug) cerr << "j11 = " << j11 << " j12 = " << j12 << " j21 = " << j21 << " j22 = " << j22 << " totallow = " << totallow << " totalup = " << totalup << " cllowold = " << cllowold << " cllow = " << cllow << " clupold = " << clupold << " clup = " << clup << " jlow = " << jlow << " jup = " << jup << endl
                      << "c11 = " << stararray0[iup ].cm[jup -1] << " c12 = " << stararray0[iup ].cm[jup   ] << " c21 = " << stararray0[ilow].cm[jlow-1] << " c22 = " << stararray0[ilow].cm[jlow  ] << endl
                      << "m11 = " << stararray0[iup ].m[jup -1] << " m12 = " << stararray0[iup ].m[jup   ] << " m21 = " << stararray0[ilow].m[jlow-1] << " m22 = " << stararray0[ilow].m[jlow  ] << endl
                      << "m11-m12=" << stararray0[iup ].m[jup -1]-stararray0[iup ].m[jup   ] << " m21-m22=" << stararray0[ilow].m[jlow-1]-stararray0[ilow].m[jlow  ] << endl
                      << "c11-c12=" << stararray0[iup ].cm[jup -1]-stararray0[iup ].cm[jup   ] << " c21-c22=" << stararray0[ilow].cm[jlow-1]-stararray0[ilow].cm[jlow  ] << endl;
      if (fmin(fmin(stararray0[iup ].cm[jup -1],stararray0[iup ].cm[jup   ]),fmin(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow  ]))>star.cm[n]){	//check if solution is possible
        if (fmax(fmax(stararray0[iup ].m[jup -1],stararray0[iup ].m[jup   ]),fmax(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow  ]))<star.m[n]){	//check if solution is possible in this track
          if (iup<nstar0-1){	//while search is in the grid
            iup++;
            ilow++;
          }else{
            star.stage[n] = -3;	//let star outside the gird explode
            screen = true;
            if (debug||true) cerr << "end of grid and too small core mass: iup=" << iup << endl;
            goto notrackupdate;
          }
          for (jlow=1;jlow<stararray0[ilow].n;jlow++){
            if (stararray0[ilow].m[jlow]<star.m[n]){
              break;
            }
          }
          for (jup=1;jup<stararray0[iup].n;jup++){
            if (stararray0[iup].m[jup]<star.m[n]){
              break;
            }
          }
          goto nexttrack;
        }
        cnosol++;
        if (debug) cerr << "go to next point: too small core mass" << endl;
        goto nextpoint;
/*        cerr << "#Error: no solution found: cmmin" << endl;
        screen = true;
        break;*/
      }
      if (fmax(fmax(stararray0[iup ].cm[jup -1],stararray0[iup ].cm[jup   ]),fmax(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow  ]))<star.cm[n]){	//check if solution is possible
        if (fmax(fmax(stararray0[iup ].m[jup -1],stararray0[iup ].m[jup   ]),fmax(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow  ]))<star.m[n]){	//check if solution is possible in this track
          if (iup<nstar0-1){	//while search is in the grid
            iup++;
            ilow++;
          }else{
            star.stage[n] = -3;	//let star outside the gird explode
            //screen = true;
            if (debug) cerr << "end of grid and too large core mass: iup=" << iup << endl;
            goto notrackupdate;
          }
          for (jlow=1;jlow<stararray0[ilow].n;jlow++){
            if (stararray0[ilow].m[jlow]<star.m[n]){
              break;
            }
          }
          for (jup=1;jup<stararray0[iup].n;jup++){
            if (stararray0[iup].m[jup]<star.m[n]){
              break;
            }
          }
          goto nexttrack;
        }
        cnosol++;
        if (debug) cerr << "go to next point: too large core mass" << endl;
        goto nextpoint;
/*        if ((stararray0[iup ].cm[jup -1]<stararray0[iup ].cm[jup   ])&&(stararray0[ilow].cm[jlow-1]<stararray0[ilow].cm[jlow  ])){	//core mass is not decreasing
          cerr << "#Error: no solution found: cmmax" << endl;
          screen = true;
          break;
        }*/
      }
      if (fmin(fmin(stararray0[iup ].m[jup -1],stararray0[iup ].m[jup   ]),fmin(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow  ]))>star.m[n]){	//check if solution is possible
        cnosol++;
        if (debug) cerr << "go to next point: too small mass" << endl;
        goto nextpoint;
/*        if ((stararray0[iup ].m[jup -1]>stararray0[iup ].m[jup   ])&&(stararray0[ilow].m[jlow-1]>stararray0[ilow].m[jlow  ])){	//core mass is not increasing
          cerr << "#Error: no solution found: mmin" << endl;
          screen = true;
          break;
        }*/
      }
      if (fmax(fmax(stararray0[iup ].m[jup -1],stararray0[iup ].m[jup   ]),fmax(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow  ]))<star.m[n]){	//check if solution is possible
        cnosol++;
        if (debug) cerr << "go to next point: too large mass" << endl;
        goto nextpoint;
/*        cerr << "#Error: no solution found: mmax" << endl;
        screen = true;
        break;*/
      }
      //calculate auxiliary variables
      A = getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]);
      B = getB(j21,j22,star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]);
      C = getC(j11,j12,star.m[n],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],star.cm[n],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]);
/*      A = ((long double)j22-(long double)j21)*((long double)stararray0[iup ].cm[jup -1]*(long double)stararray0[iup ].m[jup   ] - (long double)stararray0[iup ].cm[jup   ]*(long double)stararray0[iup ].m[jup -1])	//=(j22-j21)*(c11*m12-c12*m11)
        + ((long double)j12-(long double)j22)*((long double)stararray0[iup ].cm[jup -1]*(long double)stararray0[ilow].m[jlow-1] - (long double)stararray0[ilow].cm[jlow-1]*(long double)stararray0[iup ].m[jup -1])	//+(j12-j22)*(c11*m21-c21*m11)
        + ((long double)j22-(long double)j11)*((long double)stararray0[iup ].cm[jup   ]*(long double)stararray0[ilow].m[jlow-1] - (long double)stararray0[ilow].cm[jlow-1]*(long double)stararray0[iup ].m[jup   ])	//+(j22-j11)*(c12*m21-c21*m12)
        + ((long double)j12-(long double)j21)*((long double)stararray0[ilow].cm[jlow  ]*(long double)stararray0[iup ].m[jup -1] - (long double)stararray0[iup ].cm[jup -1]*(long double)stararray0[ilow].m[jlow  ])	//+(j12-j21)*(c22*m11-c11*m22)
        + ((long double)j21-(long double)j11)*((long double)stararray0[ilow].cm[jlow  ]*(long double)stararray0[iup ].m[jup   ] - (long double)stararray0[iup ].cm[jup   ]*(long double)stararray0[ilow].m[jlow  ])	//+(j21-j11)*(c22*m12-c12*m22)
        + ((long double)j12-(long double)j11)*((long double)stararray0[ilow].cm[jlow-1]*(long double)stararray0[ilow].m[jlow  ] - (long double)stararray0[ilow].cm[jlow  ]*(long double)stararray0[ilow].m[jlow-1]);	//+(j12-j11)*(c21*m22-c22*m21)
      B = ((long double)j22-(long double)j21)*((long double)star.cm[n]*((long double)stararray0[iup].m[jup-1]-(long double)stararray0[iup].m[jup]) + (long double)stararray0[iup].cm[jup-1]*((long double)stararray0[iup].m[jup]-(long double)star.m[n]) + (long double)stararray0[iup].cm[jup]*((long double)star.m[n]-(long double)stararray0[iup].m[jup-1]));	//=(j22-j21)*[c*(m11-m12)+c11*(m12-m)+c12*(m-m11)]
      C = ((long double)j12-(long double)j11)*((long double)star.cm[n]*((long double)stararray0[ilow].m[jlow-1]-(long double)stararray0[ilow].m[jlow]) + (long double)stararray0[ilow].cm[jlow-1]*((long double)stararray0[ilow].m[jlow]-(long double)star.m[n]) + (long double)stararray0[ilow].cm[jlow]*((long double)star.m[n]-(long double)stararray0[ilow].m[jlow-1]));	//=(j12-j11)*[c*(m21-m22)+c21*(m22-m)+c22*(m-m21)]*/
/*      A = (j22-j21)*(stararray0[iup ].cm[jup -1]*stararray0[iup ].m[jup   ] - stararray0[iup ].cm[jup   ]*stararray0[iup ].m[jup -1])	//=(j22-j21)*(c11*m12-c12*m11)
        + (j12-j22)*(stararray0[iup ].cm[jup -1]*stararray0[ilow].m[jlow-1] - stararray0[ilow].cm[jlow-1]*stararray0[iup ].m[jup -1])	//+(j12-j22)*(c11*m21-c21*m11)
        + (j22-j11)*(stararray0[iup ].cm[jup   ]*stararray0[ilow].m[jlow-1] - stararray0[ilow].cm[jlow-1]*stararray0[iup ].m[jup   ])	//+(j22-j11)*(c12*m21-c21*m12)
        + (j12-j21)*(stararray0[ilow].cm[jlow  ]*stararray0[iup ].m[jup -1] - stararray0[iup ].cm[jup -1]*stararray0[ilow].m[jlow  ])	//+(j12-j21)*(c22*m11-c11*m22)
        + (j21-j11)*(stararray0[ilow].cm[jlow  ]*stararray0[iup ].m[jup   ] - stararray0[iup ].cm[jup   ]*stararray0[ilow].m[jlow  ])	//+(j21-j11)*(c22*m12-c12*m22)
        + (j12-j11)*(stararray0[ilow].cm[jlow-1]*stararray0[ilow].m[jlow  ] - stararray0[ilow].cm[jlow  ]*stararray0[ilow].m[jlow-1]);	//+(j12-j11)*(c21*m22-c22*m21)
      B = (j22-j21)*(star.cm[n]*(stararray0[iup].m[jup-1]-stararray0[iup].m[jup]) + stararray0[iup].cm[jup-1]*(stararray0[iup].m[jup]-star.m[n]) + stararray0[iup].cm[jup]*(star.m[n]-stararray0[iup].m[jup-1]));	//=(j22-j21)*[c*(m11-m12)+c11*(m12-m)+c12*(m-m11)]
      C = (j12-j11)*(star.cm[n]*(stararray0[ilow].m[jlow-1]-stararray0[ilow].m[jlow]) + stararray0[ilow].cm[jlow-1]*(stararray0[ilow].m[jlow]-star.m[n]) + stararray0[ilow].cm[jlow]*(star.m[n]-stararray0[ilow].m[jlow-1]));	//=(j12-j11)*[c*(m21-m22)+c21*(m22-m)+c22*(m-m21)]*/
      if (debug) cerr << "A=" << A  << " B=" << B << "=" << getB(j21,j22,star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]) << " C=" << C << "=" << getC(j11,j12,star.m[n],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],star.cm[n],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << endl;
/*      if (debug) cerr << "Variation: (1+acc)*m, (1+acc)*m11, (1+acc)*m12, (1+acc)*m21, (1+acc)*m22, (1+acc)*cm, (1+acc)*c11, (1+acc)*c12, (1+acc)*c21, (1+acc)*c22" << endl
                      << "A = " << A << ", " << getA(j11,j12,j21,j22,(1.0+accuracy)*stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", " <<
                                                getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],(1.0+accuracy)*stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", " <<
                                                getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],stararray0[iup].m[jup],(1.0+accuracy)*stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", " <<
                                                getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],(1.0+accuracy)*stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", " <<
                                   A << ", " << getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],(1.0+accuracy)*stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", " <<
                                                getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],(1.0+accuracy)*stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", " <<
                                                getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],(1.0+accuracy)*stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", " <<
                                                getA(j11,j12,j21,j22,stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],(1.0+accuracy)*stararray0[ilow].cm[jlow]) << endl
                      << "B = " << getB(j21,j22,(1.0+accuracy)*star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]) << ", "
                      << getB(j21,j22,star.m[n],(1.0+accuracy)*stararray0[iup].m[jup-1],stararray0[iup].m[jup],star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]) << ", "
                      << getB(j21,j22,star.m[n],stararray0[iup].m[jup-1],(1.0+accuracy)*stararray0[iup].m[jup],star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]) << ", " << B << ", " << B << ", "
                      << getB(j21,j22,star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],(1.0+accuracy)*star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]) << ", "
                      << getB(j21,j22,star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],star.cm[n],(1.0+accuracy)*stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]) << ", "
                      << getB(j21,j22,star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],star.cm[n],stararray0[iup].cm[jup-1],(1.0+accuracy)*stararray0[iup].cm[jup]) << ", " << B << ", " << B << endl
                      << "C = " << getC(j11,j12,(1.0+accuracy)*star.m[n],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],star.cm[n],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", "
      << C << ", " << C << ", " << getC(j11,j12,star.m[n],(1.0+accuracy)*stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],star.cm[n],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", "
                                << getC(j11,j12,star.m[n],stararray0[ilow].m[jlow-1],(1.0+accuracy)*stararray0[ilow].m[jlow],star.cm[n],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", "
                                << getC(j11,j12,star.m[n],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],(1.0+accuracy)*star.cm[n],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", "
      << C << ", " << C << ", " << getC(j11,j12,star.m[n],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],star.cm[n],(1.0+accuracy)*stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]) << ", "
                                << getC(j11,j12,star.m[n],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],star.cm[n],stararray0[ilow].cm[jlow-1],(1.0+accuracy)*stararray0[ilow].cm[jlow]) << endl;*/
      if ((A==0.0)||(fabs(A)<fabs(B)*accuracy)){	//linear solution of 0=(C-B-A)*r+B
        r = (long double)B/((long double)A+(long double)B-(long double)C);
//        if (fmin(fabs(stararray0[iup].m[jup-1]-stararray0[iup].m[jup]),fabs(stararray0[ilow].m[jlow-1]-stararray0[ilow].m[jlow]))>fmin(fabs(stararray0[iup].cm[jup-1]-stararray0[iup].cm[jup]),fabs(stararray0[ilow].cm[jlow-1]-stararray0[ilow].cm[jlow]))){	//use masses
        if (fmin(fmax(stararray0[iup].m[jup-1],stararray0[iup].m[jup])/fmin(stararray0[iup].m[jup-1],stararray0[iup].m[jup]),fmax(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow])/fmin(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow]))>fmin(fmax(stararray0[iup].cm[jup-1],stararray0[iup].cm[jup])/fmin(stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]),fmax(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow])/fmin(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]))){	//use masses
          r1 = getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
          r2 = getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
//          r1 = (((long double)star.m[n]-(long double)stararray0[iup].m[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].m[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].m[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in upper track: {(m-m11*(1-r))*(j22-j21)+r*[m22*(j21-j11)-m21*(j22-j11)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
//          r2 = (((long double)star.m[n]-(long double)stararray0[ilow].m[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].m[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].m[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in lower track: {(m-m21*r)*(j12-j11)+(1-r)*[m12*(j11-j21)-m11*(j12-j21)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
          if (((r1<=0.0)&&(r1+accuracy>0.0))||((r2<=0.0)&&(r2+accuracy>0.0))){	//value not trustable => use core masses
            r1 = getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
            r2 = getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
          }
        }else{	//use core masses
          r1 = getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
          r2 = getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
//          r1 = (((long double)star.cm[n]-(long double)stararray0[iup].cm[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].cm[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].cm[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in upper track: {(c-c11*(1-r))*(j22-j21)+r*[c22*(j21-j11)-c21*(j22-j11)]}/{(c22-c21)*r*(j12-j11)+(c12-c11)*(1-r)*(j22-j21)}
//          r2 = (((long double)star.cm[n]-(long double)stararray0[ilow].cm[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].cm[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].cm[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in lower track: {(c-c21*r)*(j12-j11)+(1-r)*[c12*(j11-j21)-c11*(j12-j21)]}/{(c22-c21)*r*(j12-j11)+(c12-c11)*(1-r)*(j22-j21)}
          if (((r1<=0.0)&&(r1+accuracy>0.0))||((r2<=0.0)&&(r2+accuracy>0.0))){	//value not trustable => use masses
            r1 = getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
            r2 = getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
          }
        }
        if ((0.0<=r)&&(r<=1.0)&&(0.0<=r1)&&(r1<=1.0)&&(0.0<=r2)&&(r2<=1.0)){	//solution found
          if (debug) cout << "linear solution: r=" << r << " r1=" << r1 << " r2=" << r2 << endl;
          mratio = r;	//save ratio
          break;
        }
        if (debug) cerr << "linear: r=" << r << " r1=" << r1 << "=" << getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r) << "=" << getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r) << " r2=" << r2 << "=" << getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r) << "=" << getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r) << endl;
//        if (debug) cerr << "linear: r=" << r << " r1=" << r1 << "=" << (((long double)star.m[n]-(long double)stararray0[iup].m[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].m[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].m[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << "=" << (((long double)star.cm[n]-(long double)stararray0[iup].cm[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].cm[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].cm[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << " r2=" << r2 << "=" << (((long double)star.m[n]-(long double)stararray0[ilow].m[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].m[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].m[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << "=" << (((long double)star.cm[n]-(long double)stararray0[ilow].cm[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].cm[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].cm[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << endl;
      }else{	//quadratic solution of 0=A*r**2+(C-B-A)*r+B
        r = (-0.5*((long double)C-(long double)B-(long double)A)+sqrt(0.25*((long double)C-(long double)B-(long double)A)*((long double)C-(long double)B-(long double)A)-(long double)B*(long double)A))/(long double)A;
//        if (fmin(fabs(stararray0[iup].m[jup-1]-stararray0[iup].m[jup]),fabs(stararray0[ilow].m[jlow-1]-stararray0[ilow].m[jlow]))>fmin(fabs(stararray0[iup].cm[jup-1]-stararray0[iup].cm[jup]),fabs(stararray0[ilow].cm[jlow-1]-stararray0[ilow].cm[jlow]))){	//use masses
        if (fmin(fmax(stararray0[iup].m[jup-1],stararray0[iup].m[jup])/fmin(stararray0[iup].m[jup-1],stararray0[iup].m[jup]),fmax(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow])/fmin(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow]))>fmin(fmax(stararray0[iup].cm[jup-1],stararray0[iup].cm[jup])/fmin(stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]),fmax(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow])/fmin(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]))){	//use masses
          r1 = getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
          r2 = getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
//          r1 = (((long double)star.m[n]-(long double)stararray0[iup].m[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].m[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].m[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in upper track: {(m-m11*(1-r))*(j22-j21)+r*[m22*(j21-j11)-m21*(j22-j11)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
//          r2 = (((long double)star.m[n]-(long double)stararray0[ilow].m[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].m[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].m[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in lower track: {(m-m21*r)*(j12-j11)+(1-r)*[m12*(j11-j21)-m11*(j12-j21)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
          if (((r1<=0.0)&&(r1+accuracy>0.0))||((r2<=0.0)&&(r2+accuracy>0.0))){	//value not trustable => use core masses
            r1 = getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
            r2 = getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
          }
        }else{	//use core masses
          r1 = getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
          r2 = getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
//          r1 = (((long double)star.cm[n]-(long double)stararray0[iup].cm[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].cm[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].cm[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in upper track: {(c-c11*(1-r))*(j22-j21)+r*[c22*(j21-j11)-c21*(j22-j11)]}/{(c22-c21)*r*(j12-j11)+(c12-c11)*(1-r)*(j22-j21)}
//          r2 = (((long double)star.cm[n]-(long double)stararray0[ilow].cm[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].cm[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].cm[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in lower track: {(c-c21*r)*(j12-j11)+(1-r)*[c12*(j11-j21)-c11*(j12-j21)]}/{(c22-c21)*r*(j12-j11)+(c12-c11)*(1-r)*(j22-j21)}
          if (((r1<=0.0)&&(r1+accuracy>0.0))||((r2<=0.0)&&(r2+accuracy>0.0))){	//value not trustable => use masses
            r1 = getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
            r2 = getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
          }
        }
        if ((0.0<=r)&&(r<=1.0)&&(0.0<=r1)&&(r1<=1.0)&&(0.0<=r2)&&(r2<=1.0)){	//solution found
          if (debug) cerr << "quadratic solution: r=" << r << " r1=" << r1 << " r2=" << r2 << endl;
          mratio = r;	//save ratio
          break;
        }
        if (debug) cerr << "quadratic: r=" << r << " r1=" << r1 << "=" << getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r) << "=" << getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r) << " r2=" << r2 << "=" << getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r) << "=" << getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
//        if (debug) cerr << " variation: r=" << (1.0+accuracy)*r << " r1=" << getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << "=" << getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << " r2=" << getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << "=" << getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << endl;
//        if (debug) cerr << "quadratic: r=" << r << " r1=" << r1 << "=" << (((long double)star.m[n]-(long double)stararray0[iup].m[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].m[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].m[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << "=" << (((long double)star.cm[n]-(long double)stararray0[iup].cm[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].cm[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].cm[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << " r2=" << r2 << "=" << (((long double)star.m[n]-(long double)stararray0[ilow].m[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].m[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].m[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << "=" << (((long double)star.cm[n]-(long double)stararray0[ilow].cm[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].cm[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].cm[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));
        r = (-0.5*((long double)C-(long double)B-(long double)A)-sqrt(0.25*((long double)C-(long double)B-(long double)A)*((long double)C-(long double)B-(long double)A)-(long double)B*(long double)A))/(long double)A;
//        if (fmin(fabs(stararray0[iup].m[jup-1]-stararray0[iup].m[jup]),fabs(stararray0[ilow].m[jlow-1]-stararray0[ilow].m[jlow]))>fmin(fabs(stararray0[iup].cm[jup-1]-stararray0[iup].cm[jup]),fabs(stararray0[ilow].cm[jlow-1]-stararray0[ilow].cm[jlow]))){	//use masses
        if (fmin(fmax(stararray0[iup].m[jup-1],stararray0[iup].m[jup])/fmin(stararray0[iup].m[jup-1],stararray0[iup].m[jup]),fmax(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow])/fmin(stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow]))>fmin(fmax(stararray0[iup].cm[jup-1],stararray0[iup].cm[jup])/fmin(stararray0[iup].cm[jup-1],stararray0[iup].cm[jup]),fmax(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow])/fmin(stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow]))){	//use masses
          r1 = getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
          r2 = getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
//          r1 = (((long double)star.m[n]-(long double)stararray0[iup].m[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].m[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].m[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in upper track: {(m-m11*(1-r))*(j22-j21)+r*[m22*(j21-j11)-m21*(j22-j11)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
//          r2 = (((long double)star.m[n]-(long double)stararray0[ilow].m[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].m[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].m[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in lower track: {(m-m21*r)*(j12-j11)+(1-r)*[m12*(j11-j21)-m11*(j12-j21)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
          if (((r1<=0.0)&&(r1+accuracy>0.0))||((r2<=0.0)&&(r2+accuracy>0.0))){	//value not trustable => use core masses
            r1 = getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
            r2 = getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
          }
        }else{	//use core masses
          r1 = getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
          r2 = getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r);
//          r1 = (((long double)star.cm[n]-(long double)stararray0[iup].cm[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].cm[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].cm[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in upper track: {(c-c11*(1-r))*(j22-j21)+r*[c22*(j21-j11)-c21*(j22-j11)]}/{(c22-c21)*r*(j12-j11)+(c12-c11)*(1-r)*(j22-j21)}
//          r2 = (((long double)star.cm[n]-(long double)stararray0[ilow].cm[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].cm[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].cm[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21));	//get ratio in lower track: {(c-c21*r)*(j12-j11)+(1-r)*[c12*(j11-j21)-c11*(j12-j21)]}/{(c22-c21)*r*(j12-j11)+(c12-c11)*(1-r)*(j22-j21)}
          if (((r1<=0.0)&&(r1+accuracy>0.0))||((r2<=0.0)&&(r2+accuracy>0.0))){	//value not trustable => use masses
            r1 = getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
            r2 = getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r);
          }
        }
        if (debug) cerr << " or r=" << r << " r1=" << r1 << "=" << getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r) << "=" << getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r) << " r2=" << r2 << "=" << getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,r) << "=" << getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,r) << endl;
//        if (debug) cerr << " variation: r=" << (1.0+accuracy)*r << " r1=" << getr1(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << "=" << getr1(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << " r2=" << getr2(star.m[n],stararray0[iup].m[jup-1],stararray0[iup].m[jup],stararray0[ilow].m[jlow-1],stararray0[ilow].m[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << "=" << getr2(star.cm[n],stararray0[iup].cm[jup-1],stararray0[iup].cm[jup],stararray0[ilow].cm[jlow-1],stararray0[ilow].cm[jlow],j11,j12,j21,j22,(1.0+accuracy)*r) << endl;
//        if (debug) cerr << " or r=" << r << " r1=" << r1 << "=" << (((long double)star.m[n]-(long double)stararray0[iup].m[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].m[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].m[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << "=" << (((long double)star.cm[n]-(long double)stararray0[iup].cm[jup-1]*(1.0-(long double)r))*((long double)j22-(long double)j21)+(long double)r*((long double)stararray0[ilow].cm[jlow]*((long double)j21-(long double)j11)-(long double)stararray0[ilow].cm[jlow-1]*((long double)j22-(long double)j11)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << " r2=" << r2 << "=" << (((long double)star.m[n]-(long double)stararray0[ilow].m[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].m[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].m[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].m[jlow]-(long double)stararray0[ilow].m[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].m[jup]-(long double)stararray0[iup].m[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << "=" << (((long double)star.cm[n]-(long double)stararray0[ilow].cm[jlow-1]*(long double)r)*((long double)j12-(long double)j11)+(1.0-(long double)r)*((long double)stararray0[iup].cm[jup]*((long double)j11-(long double)j21)-(long double)stararray0[iup].cm[jup-1]*((long double)j12-(long double)j21)))/(((long double)stararray0[ilow].cm[jlow]-(long double)stararray0[ilow].cm[jlow-1])*(long double)r*((long double)j12-(long double)j11)+((long double)stararray0[iup].cm[jup]-(long double)stararray0[iup].cm[jup-1])*(1.0-(long double)r)*((long double)j22-(long double)j21)) << endl;
        if ((0.0<=r)&&(r<=1.0)&&(0.0<=r1)&&(r1<=1.0)&&(0.0<=r2)&&(r2<=1.0)){	//solution found
          if (debug) cerr << "quadratic solution: r=" << r << " r1=" << r1 << " r2=" << r2 << endl;
          mratio = r;	//save ratio
          break;
        }
      }
      nextpoint:	//goto-marker
      if (j22<j12){	//go to next point at ilow
        if (debug){
          double ratioj = ratio(j22,j12,j11);
          double mup = stararray0[iup].m[jup]-ratioj*(stararray0[iup].m[jup]-stararray0[iup].m[jup-1]);
          double cmup = stararray0[iup].cm[jup]-ratioj*(stararray0[iup].cm[jup]-stararray0[iup].cm[jup-1]);
          double ratiom = ratio(star.m[n],stararray0[ilow].m[jlow],mup);
          cerr << " mlow=" << stararray0[ilow].m[jlow] << " cmlow=" << stararray0[ilow].cm[jlow] << " mup=" << mup << " cmup=" << cmup << " ratioj=" << ratioj;
          cerr << " ratiom=" << ratiom << " c=" << stararray0[ilow].cm[jlow]-ratiom*(stararray0[ilow].cm[jlow]-cmup) << endl;
        }
        jlow++;
        cllowold = cllow;
        cllow = stararray0[ilow].t[jlow];
        j21 = cllowold*totalup;
        j22 = cllow*totalup;
      }else if (j22>j12){	//go to next point at iup
        if (debug){
          double ratioj = ratio(j12,j22,j21);
          double mlow = stararray0[ilow].m[jlow]-ratioj*(stararray0[ilow].m[jlow]-stararray0[ilow].m[jlow-1]);
          double cmlow = stararray0[ilow].cm[jlow]-ratioj*(stararray0[ilow].cm[jlow]-stararray0[ilow].cm[jlow-1]);
          double ratiom = ratio(star.m[n],stararray0[iup].m[jup],mlow);
          cerr << " mup=" << stararray0[iup ].m[jup] << " cmup=" << stararray0[iup ].cm[jup] << " mlow=" << mlow << " cmlow=" << cmlow << " ratioj=" << ratioj;
          cerr << " ratiom=" << ratiom << " c=" << stararray0[iup].cm[jup]-ratiom*(stararray0[iup].cm[jup]-cmlow) << endl;
        }
        jup++;
        clupold = clup;
        clup = stararray0[iup].t[jup];
        j11 = clupold*totallow;
        j12 = clup*totallow;
      }else{	//go to next point at ilow and iup
        jlow++;
        jup++;
        cllowold = cllow;
        cllow = stararray0[ilow].t[jlow];
        clupold = clup;
        clup = stararray0[iup].t[jup];
        j21 = cllowold*totalup;
        j22 = cllow*totalup;
        j11 = clupold*totallow;
        j12 = clup*totallow;
      }
      if ((jup>jmax)&&(r1>1.0)) jmax = nup;
    }
  }

  if (mratio>=0){
    if (SN){
      star.SN.Hecoremass = stararray0[iup].m[nup]-mratio*(stararray0[iup].m[nup]-stararray0[ilow].m[nlow]);	//get He mass of last track point
      if ((star.SN.COcoremass>0.01)||((star.track.cc[star.track.n-1]>0.01)&&(stararray0[ilow].cm[nlow]>0.01))||(star.SN.Hecoremass>2.0*WDmax)) star.SN.COcoremass = stararray0[iup].cm[nup]-mratio*(stararray0[iup].cm[nup]-stararray0[ilow].cm[nlow]);	//get CO core mass of last track point, if star could ignite core helium burning
// (Hecoremass and inimass limits are taken randomly, but conservative)		||(star.SN.Hecoremass>5.0)||(star.inimass>30.0)
    }else{
      getnewtrack(stararray0, ilow, iup, mratio, 0, newtrack, star.rmax);
/*      if (star.cm[n]<newtrack.cm[1]){	//no significant core
        if (star.cm[n]>cmup) cerr << "#Warning: overwrite newtrack.cm[1]" << endl;
        newtrack.cm[1] = star.cm[n];	//overwrite frist track value with current core mass
      }*/
      changetrack(star, n, newtrack);
      if (star.stage[n]==2) star.stage[n] = 2;
      else if (star.last<star.track.TAMS) star.stage[n] = 0;
      else star.stage[n] = 1;
      //free memory
      free(newtrack.m);
      free(newtrack.t);
      free(newtrack.r);
      free(newtrack.cm);
      free(newtrack.llum);
      free(newtrack.lteff);
      free(newtrack.lambda);
      free(newtrack.cc);
      free(newtrack.cf);
      free(newtrack.omega);
    }
  }else{
    if ((star.stage[n]==0)&&(stararray1!=stararray2)){	//Hydrogen burning star -> check for Helium burning star
      star.stage[n] = 1;
      updatestarstrack(star, n, stararray2, nstar2, SN);
    }else if ((jlow>=nlow)&&(jup>=nup)){	//no solution at last track point
      if (fmin(stararray0[iup].m[nup-1],stararray0[iup].m[nup])<star.m[n]){	//check if solution could be possible in next track
        if (iup<nstar0-1){	//while search is in the grid
          iup++;
          ilow++;
        }else{
          star.stage[n] = -3;	//let star outside the gird explode
          screen = true;
          if (debug||true) cerr << "at last track point and end of grid: iup=" << iup << endl;
          goto notrackupdate;
        }
        for (jlow=1;jlow<stararray0[ilow].n;jlow++){
          if (stararray0[ilow].m[jlow]<star.m[n]){
            break;
          }
        }
        for (jup=1;jup<stararray0[iup].n;jup++){
          if (stararray0[iup].m[jup]<star.m[n]){
            break;
          }
        }
        goto nexttrack;
      }
      mratio = ratio(star.m[n],stararray0[iup].m[nup],stararray0[ilow].m[nlow]);	//calculate ratio between last track points
      cmup = stararray0[iup].cm[nup]-mratio*(stararray0[iup].cm[nup]-stararray0[ilow].cm[nlow]);	//calculate largest possible core mass between last track points
      if (debug) cerr << "at last track point: cmup=" << cmup << "Msun mratio=" << mratio << endl;
      if (star.cm[n]>=cmup){
        star.stage[n] = -3;
        for (iup=1;iup<nstar0;iup++){
          nup = stararray0[iup].n-1;
          if (stararray0[iup].cm[nup]>star.cm[n]){	//find track with core mass at the end
            break;
          }
        }
        if (iup>=nstar0){	//too large core mass
          for (iup=1;iup<nstar0-1;iup++){
            if (stararray0[iup].cm[stararray0[iup].n-1]<stararray0[iup-1].cm[stararray0[iup-1].n-1]){	//find maximum core mass at the end
              break;
            }
          }
          star.inimass = stararray0[iup-1].m[0];	//mass with largest core mass
        }else{
          ilow = iup-1;
          nlow = stararray0[ilow].n-1;
          mratio = ratio(star.cm[n],stararray0[iup].cm[nup],stararray0[ilow].cm[nlow]);	//calculate ratio between last track points
          star.inimass = stararray0[iup].m[0]-mratio*(stararray0[iup].m[0]-stararray0[ilow].m[0]);	//calculate initial mass of a star with this core mass
          if (debug) cerr << "new inimass=" << star.inimass << "Msun star.m[" << n << "]=" << star.m[n] << "Msun star.cm[" << n << "]=" << star.cm[n] << "Msun star.cc[" << n << "]=" << star.cc[n] << "Msun" << endl;
        }
      }else goto notrackupdate;
    }else{
      notrackupdate:	//goto-marker
      trackadaption = true;
//        cnotrackupdate++;
//        screen = true;
      if (star.stage[n]!=-3){
        if (cnosol>0) screen = true;
        if (screen) cerr << "#Warning: no track update, because star out of grid: star.m[" << n << "]=" << star.m[n] << "Msun and star.cm[" << n << "]=" << star.cm[n] << "Msun" << endl;
        undophase(star, n, true);
        if (star.rmax<=star.last){	//update rmax
          if (star.stage[n]==0) star.rmax = star.track.TAMS;	//set index of rmax at TAMS
          else star.rmax = star.track.n-1;			//set index of rmax at end track point
          maxr = star.track.r[star.rmax];	//save radius at last track point
          for (i=star.last+1;i<star.track.n;i++){
            if ((maxr<star.track.r[i])&&(star.track.r[i-1]<star.track.r[i])){
              maxr = star.track.r[i];
              if (star.track.r[0]<star.track.r[i]) star.track.r[0] = star.track.r[i];
              if ((i<=star.track.TAMS)||(star.last>=star.track.TAMS)) star.rmax = i;
            }
          }
        }
      }else{
        if (screen) cerr << "#Warning: no track update, because star out of grid: star.m[" << n << "]=" << star.m[n] << "Msun and star.cm[" << n << "]=" << star.cm[n] << "Msun -> star will explode" << endl;
      }
    }
  }
  if (trackadaption){
//    screen = true;
    cnotrackupdate++;	//increase counter
  }
}

/*void updateHetrack(t_star& star, int n){
  int j,jlow,jup;	//index variables
  int ilow=0,iup=0,nlow=0,nup=0;	//index variables: ilow/iup=mass index
  double cmratio,ratiolow,ratioup,mratio;	//core mass ratio, ratios at lower/upper mass track, mass ratio
  double mw,t,r,cm,llum,lteff;	//interpolated values of mass, age, radius, core mass, log10(luminosity), log10(effective temperature)
  double timecorr;	//timecorrection
  t_HRD newtrack;
  bool cmchange=false;	//tells if core mass has to be updated

  //assume a constant core mass while the naked helium phase
  for (iup=1;iup<nhe-1;iup++){	//find mass index
    if(hestararray[iup].cm[0]>star.cm[n]){	//hestararray must be sorted from low (core) mass to high (core) mass (done when table is read in)
      break;
    }
  }
//  if (screen) cout << "iup=" << iup << " n_he=" << nhe << " hestararray[iup].cm[0]=" << hestararray[iup].cm[0] << "Msun" << endl;
  ilow = iup-1;
  //calculate the ratio between the two nearest evolutionary tracks
  cmratio = ratio(star.cm[n],hestararray[iup].cm[0],hestararray[ilow].cm[0]);	//calculate ratio between masses
  if (star.m[n]>hestararray[iup].mw[0]-cmratio*(hestararray[iup].mw[0]-hestararray[ilow].mw[0])){	//no helium star with same mass and core mass possible => set star to beginning of naked helium stage by using its total mass
    cmchange = true;
//    screen = true;
    if (screen) cerr << "#Warning: set star back to beginning of naked helium stage" << endl;
    for (iup=1;iup<nhe-1;iup++){	//find mass index
      if(hestararray[iup].m[0]>star.m[n]){	//hestararray must be sorted from low mass to high mass (done when table is read in)
        break;
      }
    }
//    if (screen) cout << "iup=" << iup << " n_he=" << nhe << " hestararray[iup].m[0]=" << hestararray[iup].m[0] << "Msun" << endl;
    ilow = iup-1;
    //calculate the ratio between the two nearest evolutionary tracks
    cmratio = ratio(star.m[n],hestararray[iup].m[0],hestararray[ilow].m[0]);	//calculate ratio between masses
  }
  newtrack.n = hestararray[ilow].n+hestararray[iup].n-2;
  //reserve memory
  newtrack.m = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.t = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.r = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.cm = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.llum = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.lteff = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.lambda = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.cc = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.cf = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.omega = (double *)malloc(newtrack.n*sizeof(double));
  //check if memory allocation fails
  if (newtrack.m==NULL) cerr << "#Error: memory allocation failed: newtrack.m" << endl;
  if (newtrack.t==NULL) cerr << "#Error: memory allocation failed: newtrack.t" << endl;
  if (newtrack.r==NULL) cerr << "#Error: memory allocation failed: newtrack.r" << endl;
  if (newtrack.cm==NULL) cerr << "#Error: memory allocation failed: newtrack.cm" << endl;
  if (newtrack.llum==NULL) cerr << "#Error: memory allocation failed: newtrack.llum" << endl;
  if (newtrack.lteff==NULL) cerr << "#Error: memory allocation failed: newtrack.lteff" << endl;
  if (newtrack.lambda==NULL) cerr << "#Error: memory allocation failed: newtrack.lambda" << endl;
  if (newtrack.cc==NULL) cerr << "#Error: memory allocation failed: newtrack.cc" << endl;
  if (newtrack.cf==NULL) cerr << "#Error: memory allocation failed: newtrack.cf" << endl;
  if (newtrack.omega==NULL) cerr << "#Error: memory allocation failed: newtrack.omega" << endl;
  //calculate new track
  j = 0;
  newtrack.m[j] = hestararray[iup].mw[0]-cmratio*(hestararray[iup].mw[0]-hestararray[ilow].mw[0]);
  newtrack.t[j] = fmax(hestararray[iup].t[0]-cmratio*(hestararray[iup].t[0]-hestararray[ilow].t[0]),0.0);	//minimal age is 0yr
  newtrack.r[j] = hestararray[iup].r[0]-cmratio*(hestararray[iup].r[0]-hestararray[ilow].r[0]);
  if (newtrack.r[j]<rsmax) newtrack.r[j] = fmax(newtrack.r[j],2.0*G*newtrack.m[j]/(c*c));	//minimal radius is the Schwarzschild radius
  newtrack.cm[j] = hestararray[iup].cm[0]-cmratio*(hestararray[iup].cm[0]-hestararray[ilow].cm[0]);
  if (newtrack.cm[j]>newtrack.m[j]){	//if core mass is larger than total mass -> take the mean for both
    newtrack.cm[j] = 0.5*(newtrack.cm[j]+newtrack.m[j]);
    newtrack.m[j] = newtrack.cm[j];
  }
  newtrack.llum[j] = hestararray[iup].llum[0]-cmratio*(hestararray[iup].llum[0]-hestararray[ilow].llum[0]);
  newtrack.lteff[j] = hestararray[iup].lteff[0]-cmratio*(hestararray[iup].lteff[0]-hestararray[ilow].lteff[0]);
  newtrack.lambda[j] = lambda_const;
  if ((cmratio<0)||(cmratio>1)){	//values are extrapolated
    newtrack.m[j] = fmax(newtrack.m[j],1.0e-10);	//minimal mass of 10^{-10}Msun
    newtrack.cm[j] = fmax(newtrack.cm[j],0.0);	//minimal core mass is 0
  }
  newtrack.cc[j] = newtrack.cm[j];
  newtrack.cf[j] = 0.0;
  newtrack.omega[j] = 0.0;
  if ((star.m[n]>newtrack.m[0])&&(star.stage[n]==0)){	//too less new material to get normal star and too much material the convert the material to helium
//    screen = true;
    star.stage[n]=-3;	//produces black hole - white dwarf systems
    //free memory
    free(newtrack.m);
    free(newtrack.t);
    free(newtrack.r);
    free(newtrack.cm);
    free(newtrack.llum);
    free(newtrack.lteff);
    free(newtrack.lambda);
    free(newtrack.cc);
    free(newtrack.cf);
    free(newtrack.omega);
    return;
  }
  jlow = 1;
  jup = 1;
  nlow = hestararray[ilow].n-1;
  nup = hestararray[iup].n-1;
  for (j=1;j<newtrack.n-1;j++){
    if ((jlow>nlow)||(jup>nup)){
      cerr << "#Error: j=" << j << " jlow=" << jlow << " nlow=" << nlow << " jup=" << jup << " nup=" << nup << endl;
      screen = true;
      break;
    }else if (jlow*nup<jup*nlow){
      //calculate ratio in upper evolutionary track
      ratioup = ratio(jlow*nup,jup*nlow,(jup-1)*nlow);
      //determine values at upper evolutionary track
      mw = hestararray[iup].mw[jup]-ratioup*(hestararray[iup].mw[jup]-hestararray[iup].mw[jup-1]);
      t = hestararray[iup].t[jup]-ratioup*(hestararray[iup].t[jup]-hestararray[iup].t[jup-1]);
      r = hestararray[iup].r[jup]-ratioup*(hestararray[iup].r[jup]-hestararray[iup].r[jup-1]);
      cm = hestararray[iup].cm[jup]-ratioup*(hestararray[iup].cm[jup]-hestararray[iup].cm[jup-1]);
      llum = hestararray[iup].llum[jup]-ratioup*(hestararray[iup].llum[jup]-hestararray[iup].llum[jup-1]);
      lteff = hestararray[iup].lteff[jup]-ratioup*(hestararray[iup].lteff[jup]-hestararray[iup].lteff[jup-1]);
      //interpolate/extrapolate between evolutionary tracks
      newtrack.m[j] = fmax(mw-cmratio*(mw-hestararray[ilow].mw[jlow]),1.0e-10);	//minimal mass of 10^{-10}Msun
      newtrack.t[j] = fmax(t-cmratio*(t-hestararray[ilow].t[jlow]),0.0);	//minimal age is 0yr
      newtrack.r[j] = r-cmratio*(r-hestararray[ilow].r[jlow]);
      newtrack.cm[j] = cm-cmratio*(cm-hestararray[ilow].cm[jlow]);
      if (newtrack.cm[j]>newtrack.m[j]){	//if core mass is larger than total mass -> take the mean for both
        newtrack.cm[j] = 0.5*(newtrack.cm[j]+newtrack.m[j]);
        newtrack.m[j] = newtrack.cm[j];
      }
      newtrack.llum[j] = llum-cmratio*(llum-hestararray[ilow].llum[jlow]);
      newtrack.lteff[j] = lteff-cmratio*(lteff-hestararray[ilow].lteff[jlow]);
      newtrack.lambda[j] = lambda_const;
      newtrack.cc[j] = newtrack.cm[j];
      newtrack.cf[j] = 0.0;
      newtrack.omega[j] = 0.0;
      //next trackpoint
      jlow++;
    }else{
      //calculate ratio in lower evolutionary track
      ratiolow = ratio(jup*nlow,jlow*nup,(jlow-1)*nup);
      //determine values at lower evolutionary track
      mw = hestararray[ilow].mw[jlow]-ratiolow*(hestararray[ilow].mw[jlow]-hestararray[ilow].mw[jlow-1]);
      t = hestararray[ilow].t[jlow]-ratiolow*(hestararray[ilow].t[jlow]-hestararray[ilow].t[jlow-1]);
      r = hestararray[ilow].r[jlow]-ratiolow*(hestararray[ilow].r[jlow]-hestararray[ilow].r[jlow-1]);
      cm = hestararray[ilow].cm[jlow]-ratiolow*(hestararray[ilow].cm[jlow]-hestararray[ilow].cm[jlow-1]);
      llum = hestararray[ilow].llum[jlow]-ratiolow*(hestararray[ilow].llum[jlow]-hestararray[ilow].llum[jlow-1]);
      lteff = hestararray[ilow].lteff[jlow]-ratiolow*(hestararray[ilow].lteff[jlow]-hestararray[ilow].lteff[jlow-1]);
      //interpolate/extrapolate between evolutionary tracks
      newtrack.m[j] = fmax(hestararray[iup].mw[jup]-cmratio*(hestararray[iup].mw[jup]-mw),1.0e-10);	//minimal mass of 10^{-10}Msun
      newtrack.t[j] = fmax(hestararray[iup].t[jup]-cmratio*(hestararray[iup].t[jup]-t),0.0);	//minimal age is 0yr
      newtrack.r[j] = hestararray[iup].r[jup]-cmratio*(hestararray[iup].r[jup]-r);
      newtrack.cm[j] = hestararray[iup].cm[jup]-cmratio*(hestararray[iup].cm[jup]-cm);
      if (newtrack.cm[j]>newtrack.m[j]){	//if core mass is larger than total mass -> take the mean for both
        newtrack.cm[j] = 0.5*(newtrack.cm[j]+newtrack.m[j]);
        newtrack.m[j] = newtrack.cm[j];
      }
      newtrack.llum[j] = hestararray[iup].llum[jup]-cmratio*(hestararray[iup].llum[jup]-llum);
      newtrack.lteff[j] = hestararray[iup].lteff[jup]-cmratio*(hestararray[iup].lteff[jup]-lteff);
      newtrack.lambda[j] = lambda_const;
      newtrack.cc[j] = newtrack.cm[j];
      newtrack.cf[j] = 0.0;
      newtrack.omega[j] = 0.0;
      //next trackpoint
      jup++;
    }
    if ((cmratio<0)||(cmratio>1)){	//values are extrapolated
      newtrack.m[j] = fmin(newtrack.m[j],newtrack.m[j-1]);	//mass can not increase
      newtrack.t[j] = fmax(newtrack.t[j],newtrack.t[j-1]+1.0);	//age can not decrease
      newtrack.cm[j] = fmax(newtrack.cm[j],newtrack.cm[j-1]);	//core mass can not decrease
      newtrack.cc[j] = fmax(newtrack.cc[j],newtrack.cc[j-1]);	//carbon core mass can not decrease
    }
    if (newtrack.r[j]<rsmax) newtrack.r[j] = fmax(newtrack.r[j],2.0*G*newtrack.m[j]/(c*c));	//minimal radius is the Schwarzschild radius
    if ((newtrack.t[j]<newtrack.t[j-1])&&(jlow*nup<jup*nlow)){
      cerr << "#Error: newtrack.t[" << j << "]=" << newtrack.t[j] << "yr newtrack.t[" << j-1 << "]=" << newtrack.t[j-1] << "yr t=" << t << "yr hestararray[" << ilow << "].t[" << jlow << "]=" << hestararray[ilow].t[jlow] << "yr cm_ratio=" << cmratio << endl;
      screen = true;
    }
    if ((newtrack.t[j]<newtrack.t[j-1])&&(jlow*nup>=jup*nlow)){
      cerr << "#Error: newtrack.t[" << j << "]=" << newtrack.t[j] << "yr newtrack.t[" << j-1 << "]=" << newtrack.t[j-1] << "yr hestararray[" << iup << "].t[" << jup << "]=" << hestararray[iup].t[jup] << "yr t=" << t << "yr cm_ratio=" << cmratio << endl;
      screen = true;
    }
  }
  if ((jlow!=nlow)||(jup!=nup)) cerr << "#Warning: jlow=" << jlow << " nlow=" << nlow << " ilow=" << ilow << " jup=" << jup << " nup=" << nup << " iup=" << iup << endl;
  j = newtrack.n-1;
  newtrack.m[j] = fmax(hestararray[iup].mw[nup]-cmratio*(hestararray[iup].mw[nup]-hestararray[ilow].mw[nlow]),1.0e-10);	//minimal mass of 10^{-10}Msun
  newtrack.t[j] = fmax(hestararray[iup].t[nup]-cmratio*(hestararray[iup].t[nup]-hestararray[ilow].t[nlow]),0.0);	//minimal age is 0yr
  newtrack.r[j] = hestararray[iup].r[nup]-cmratio*(hestararray[iup].r[nup]-hestararray[ilow].r[nlow]);
  newtrack.cm[j] = hestararray[iup].cm[nup]-cmratio*(hestararray[iup].cm[nup]-hestararray[ilow].cm[nlow]);
  if (newtrack.cm[j]>newtrack.m[j]){	//if core mass is larger than total mass -> take the mean for both
    newtrack.cm[j] = 0.5*(newtrack.cm[j]+newtrack.m[j]);
    newtrack.m[j] = newtrack.cm[j];
  }
  newtrack.llum[j] = hestararray[iup].llum[nup]-cmratio*(hestararray[iup].llum[nup]-hestararray[ilow].llum[nlow]);
  newtrack.lteff[j] = hestararray[iup].lteff[nup]-cmratio*(hestararray[iup].lteff[nup]-hestararray[ilow].lteff[nlow]);
  newtrack.lambda[j] = lambda_const;
  newtrack.cc[j] = newtrack.cm[j];
  newtrack.cf[j] = 0.0;
  newtrack.omega[j] = 0.0;
  if ((cmratio<0)||(cmratio>1)){	//values are extrapolated
    newtrack.m[j] = fmin(newtrack.m[j],newtrack.m[j-1]);	//mass can not increase
    newtrack.t[j] = fmax(newtrack.t[j],newtrack.t[j-1]+1.0);	//age can not decrease
    newtrack.cm[j] = fmax(newtrack.cm[j],newtrack.cm[j-1]);	//core mass can not decrease
    newtrack.cc[j] = fmax(newtrack.cc[j],newtrack.cc[j-1]);	//carbon core mass can not decrease
  }
  if (newtrack.r[j]<rsmax) newtrack.r[j] = fmax(newtrack.r[j],2.0*G*newtrack.m[j]/(c*c));	//minimal radius is the Schwarzschild radius
  if (newtrack.t[j]<newtrack.t[j-1]){
    cerr << "#Error: newtrack.t[" << j << "]=" << newtrack.t[j] << "yr newtrack.t[" << j-1 << "]=" << newtrack.t[j-1] << "yr hestararray[" << iup << "].t[" << jup << "]=" << hestararray[iup].t[jup] << "yr hestararray[" << ilow << "].t[" << nlow << "]=" << hestararray[ilow].t[nlow] << "yr cm_ratio=" << cmratio << endl;
    screen = true;
  }
  for (jup=1;jup<newtrack.n;jup++){	//find time index in new track
    if(newtrack.m[jup]<star.m[n]){	//newtrack must be sorted from high mass to low mass (done when new track is created)
      break;
    }
  }
  if (jup==newtrack.n){	//beyond last track point
    star.stage[n]=-3;
  }else{	//appley new track
    jlow = jup-1;
    //set new array length
    star.track.n = star.last+1+newtrack.n-jlow;
    //reserve memory of new length
    star.track.m = (double *)realloc(star.track.m,star.track.n*sizeof(double));
    star.track.t = (double *)realloc(star.track.t,star.track.n*sizeof(double));
    star.track.r = (double *)realloc(star.track.r,star.track.n*sizeof(double));
    star.track.cm = (double *)realloc(star.track.cm,star.track.n*sizeof(double));
    star.track.llum = (double *)realloc(star.track.llum,star.track.n*sizeof(double));
    star.track.lteff = (double *)realloc(star.track.lteff,star.track.n*sizeof(double));
    star.track.lambda = (double *)realloc(star.track.lambda,star.track.n*sizeof(double));
    star.track.cc = (double *)realloc(star.track.cc,star.track.n*sizeof(double));
    star.track.cf = (double *)realloc(star.track.cf,star.track.n*sizeof(double));
    star.track.omega = (double *)realloc(star.track.omega,star.track.n*sizeof(double));
    //check if memory allocation fails
    if (star.track.m==NULL) cerr << "#Error: memory reallocation failed: star.track.m" << endl;
    if (star.track.t==NULL) cerr << "#Error: memory reallocation failed: star.track.t" << endl;
    if (star.track.r==NULL) cerr << "#Error: memory reallocation failed: star.track.r" << endl;
    if (star.track.cm==NULL) cerr << "#Error: memory reallocation failed: star.track.cm" << endl;
    if (star.track.llum==NULL) cerr << "#Error: memory reallocation failed: star.track.llum" << endl;
    if (star.track.lteff==NULL) cerr << "#Error: memory reallocation failed: star.track.lteff" << endl;
    if (star.track.lambda==NULL) cerr << "#Error: memory reallocation failed: star.track.lambda" << endl;
    if (star.track.cc==NULL) cerr << "#Error: memory reallocation failed: star.track.cc" << endl;
    if (star.track.cf==NULL) cerr << "#Error: memory reallocation failed: star.track.cf" << endl;
    if (star.track.omega==NULL) cerr << "#Error: memory reallocation failed: star.track.omega" << endl;
    mratio = ratio(star.m[n],newtrack.m[jlow],newtrack.m[jup]);	//calculate ratio between masses
    if ((mratio<0)||(mratio>1)){
      cerr << "#Error: star.m[" << n << "]=" << star.m[n] << "Msun newtrack.m[" << jlow << "]=" << newtrack.m[jlow] << "Msun newtrack.m[" << jup << "]=" << newtrack.m[jup] << "Msun mratio=" << mratio << endl;
      screen = true;
    }
    star.last++;
    //update track for naked helium star
    star.inimass = newtrack.m[0];
    star.track.m[star.last] = newtrack.m[jlow]-mratio*(newtrack.m[jlow]-newtrack.m[jup]);
    star.track.t[star.last] = newtrack.t[jlow]-mratio*(newtrack.t[jlow]-newtrack.t[jup]);
    star.track.r[star.last] = newtrack.r[jlow]-mratio*(newtrack.r[jlow]-newtrack.r[jup]);
    star.track.cm[star.last] = newtrack.cm[jlow]-mratio*(newtrack.cm[jlow]-newtrack.cm[jup]);
    star.track.llum[star.last] = newtrack.llum[jlow]-mratio*(newtrack.llum[jlow]-newtrack.llum[jup]);
    star.track.lteff[star.last] = newtrack.lteff[jlow]-mratio*(newtrack.lteff[jlow]-newtrack.lteff[jup]);
    star.track.lambda[star.last] = newtrack.lambda[jlow]-mratio*(newtrack.lambda[jlow]-newtrack.lambda[jup]);
    star.track.cc[star.last] = newtrack.cc[jlow]-mratio*(newtrack.cc[jlow]-newtrack.cc[jup]);
    star.track.cf[star.last] = newtrack.cf[jlow]-mratio*(newtrack.cf[jlow]-newtrack.cf[jup]);
    star.track.omega[star.last] = newtrack.omega[jlow]-mratio*(newtrack.omega[jlow]-newtrack.omega[jup]);
    for (j=1;j<newtrack.n-jlow;j++){
      star.track.m[star.last+j] = newtrack.m[jlow+j];
      star.track.t[star.last+j] = newtrack.t[jlow+j];
      star.track.r[star.last+j] = newtrack.r[jlow+j];
      star.track.cm[star.last+j] = newtrack.cm[jlow+j];
      star.track.llum[star.last+j] = newtrack.llum[jlow+j];
      star.track.lteff[star.last+j] = newtrack.lteff[jlow+j];
      star.track.lambda[star.last+j] = newtrack.lambda[jlow+j];
      star.track.cc[star.last+j] = newtrack.cc[jlow+j];
      star.track.cf[star.last+j] = newtrack.cf[jlow+j];
      star.track.omega[star.last+j] = newtrack.omega[jlow+j];
    }
    timecorr = star.t[n]-star.track.t[star.last];
    star.rmax = star.last;
    for (j=star.last;j<star.track.n;j++){
      star.track.t[j] += timecorr;
      if (star.track.r[j]>star.track.r[star.rmax]) star.rmax = j;	//update rmax
    }
    if (cmchange) star.cm[n] = star.track.cm[star.last];	//update the core mass
  }
  //free memory
  free(newtrack.m);
  free(newtrack.t);
  free(newtrack.r);
  free(newtrack.cm);
  free(newtrack.llum);
  free(newtrack.lteff);
  free(newtrack.lambda);
  free(newtrack.cc);
  free(newtrack.cf);
  free(newtrack.omega);
}*/

int nextstage(t_star& star, int n){
  int j;	//index variables
//  int ilow=0,iup=0;	//index variables: ilow/iup=mass index
//  int jlow,jup,nlow=0,nup=0;	//index variables
//  double mratio;	//mass ratio
//  double ratiolow,ratioup;	//ratios at lower/upper mass track
//  double mw,t,r,cm,llum,lteff;	//interpolated values of mass, age, radius, core mass, log10(luminosity), log10(effective temperature)
//  t_HRD newtrack;	//new track part

//  if (screen) cout << "nextstage: star.stage=" << star.stage[n] << endl;
  if (star.stage[n]==0){	//end hydrogen burning
/*    //increase array length
    star.track.n++;
    //enlarge reserved memory
    star.track.m = (double *)realloc(star.track.m,star.track.n*sizeof(double));
    star.track.t = (double *)realloc(star.track.t,star.track.n*sizeof(double));
    star.track.r = (double *)realloc(star.track.r,star.track.n*sizeof(double));
    star.track.cm = (double *)realloc(star.track.cm,star.track.n*sizeof(double));
    star.track.llum = (double *)realloc(star.track.llum,star.track.n*sizeof(double));
    star.track.lteff = (double *)realloc(star.track.lteff,star.track.n*sizeof(double));
    star.track.lambda = (double *)realloc(star.track.lambda,star.track.n*sizeof(double));
    star.track.cc = (double *)realloc(star.track.cc,star.track.n*sizeof(double));
    star.track.cf = (double *)realloc(star.track.cf,star.track.n*sizeof(double));
    star.track.omega = (double *)realloc(star.track.omega,star.track.n*sizeof(double));
    //check if memory allocation fails
    if (star.track.m==NULL) cerr << "#Error: memory reallocation failed: star.track.m" << endl;
    if (star.track.t==NULL) cerr << "#Error: memory reallocation failed: star.track.t" << endl;
    if (star.track.r==NULL) cerr << "#Error: memory reallocation failed: star.track.r" << endl;
    if (star.track.cm==NULL) cerr << "#Error: memory reallocation failed: star.track.cm" << endl;
    if (star.track.llum==NULL) cerr << "#Error: memory reallocation failed: star.track.llum" << endl;
    if (star.track.lteff==NULL) cerr << "#Error: memory reallocation failed: star.track.lteff" << endl;
    if (star.track.lambda==NULL) cerr << "#Error: memory reallocation failed: star.track.lambda" << endl;
    if (star.track.cc==NULL) cerr << "#Error: memory reallocation failed: star.track.cc" << endl;
    if (star.track.cf==NULL) cerr << "#Error: memory reallocation failed: star.track.cf" << endl;
    if (star.track.omega==NULL) cerr << "#Error: memory reallocation failed: star.track.omega" << endl;
    //add track for helium burning
//    star.track.m[star.track.n-1] = fmax(star.track.m[star.track.n-2]-(star.track.m[1]-star.track.m[star.track.n-2]-star.addedmass),star.track.cm[star.track.n-2]);	//assume same wind mass loss as for hydrogenburning, but mass can not fall below the core mass
    star.track.m[star.track.n-1] = star.track.m[star.track.n-2]-(star.inimass-star.track.m[star.track.n-2]);	//assume same wind mass loss as for hydrogenburning
    if (star.track.m[star.track.n-1]<0) star.track.m[star.track.n-1] = fmin(10.0,star.track.m[star.track.n-2]);	//only very high mass stars lose more than 50% of its initial mass during hydrogenburning, those will end as black holes
    star.track.t[star.track.n-1] = 1.1*star.track.t[star.track.n-2];	//assume 10% time of hydrogenburning
    star.track.r[star.track.n-1] = star.track.r[star.track.n-2];	//assume constant radius (a little bit increased)
    star.track.cm[star.track.n-1] = star.track.m[star.track.n-1];	//assume core mass = mass at the end of helium burning
    star.track.llum[star.track.n-1] = star.track.llum[star.track.n-2];	//assume constant luminosity
    star.track.lteff[star.track.n-1] = star.track.lteff[star.track.n-2];	//assume constant effective temperature
    star.track.lambda[star.track.n-1] = star.track.lambda[star.track.n-2];	//assume constant lambda
    star.track.cc[star.track.n-1] = star.track.cc[star.track.n-2];	//assume constant carbon core mass
    star.track.cf[star.track.n-1] = star.track.cf[star.track.n-2];	//assume constant concentration factor
    star.track.omega[star.track.n-1] = star.track.omega[star.track.n-2];	//assume constant angular velocity
    if (debug){
      cout << endl << "index\tm in Msun\tt in yr\tr in Rsun\tcm in Msun\tllum\tlteff\tlambda" << endl;
      cout << star.track.n-2 << "\t" << star.track.m[star.track.n-2] << "\t" << star.track.t[star.track.n-2] << "\t" << star.track.r[star.track.n-2] << "\t" << star.track.cm[star.track.n-2] << "\t" << star.track.llum[star.track.n-2] << "\t" << star.track.lteff[star.track.n-2] << "\t" << star.track.lambda[star.track.n-2] << endl;
      cout << star.track.n-1 << "\t" << star.track.m[star.track.n-1] << "\t" << star.track.t[star.track.n-1] << "\t" << star.track.r[star.track.n-1] << "\t" << star.track.cm[star.track.n-1] << "\t" << star.track.llum[star.track.n-1] << "\t" << star.track.lteff[star.track.n-1] << "\t" << star.track.lambda[star.track.n-1] << endl;
    }
    star.rmax = star.track.n-1;*/
/*    star.last = star.track.n-1;*/
    star.last = star.track.TAMS;
    if (star.t[n]<star.track.t[star.last]){
      star.m[n] = star.track.m[star.last];
      star.t[n] = star.track.t[star.last];
      star.r[n] = star.track.r[star.last];
      star.cm[n] = star.track.cm[star.last];
      star.llum[n] = star.track.llum[star.last];
      star.lteff[n] = star.track.lteff[star.last];
      star.lambda[n] = star.track.lambda[star.last];
      star.cc[n] = star.track.cc[star.last];
      star.cf[n] = star.track.cf[star.last];
      star.omega[n] = star.track.omega[star.last];
    }
    //redetermine maximal radius
    star.track.r[0] = star.track.r[star.track.TAMS];
    for (j=star.last;j<star.track.n;j++){
      if (star.track.r[0]<star.track.r[j]){
        star.track.r[0] = star.track.r[j];
        star.rmax = j;
      }
    }
/*    for (iup=0;iup<nstar2;iup++){	//find mass index
      if(stararray2[iup].inimass>star.inimass){	//stararray must be sorted from low mass to high mass (done when table is read in)
        break;
      }
    }
    ilow = iup-1;
    mratio = ratio(star.inimass,stararray2[iup].inimass,stararray2[ilow].inimass);	//calculate ratio between masses
    if (debug) cerr << "star.m[" << n << "]=" << star.m[n] << "Msun stararray2[" << ilow << "].n=" << stararray2[ilow].n << " stararray2[" << iup << "].n=" << stararray2[iup].n << " mratio=" << mratio << endl;
    getnewtrack(stararray2, ilow, iup, mratio, 0, newtrack, star.rmax);
    changetrack(star, n, newtrack);
    //free memory
    free(newtrack.m);
    free(newtrack.t);
    free(newtrack.r);
    free(newtrack.cm);
    free(newtrack.llum);
    free(newtrack.lteff);
    free(newtrack.lambda);
    free(newtrack.cc);
    free(newtrack.cf);
    free(newtrack.omega);*/
    //set stage to helium burning
    star.stage[n] = 1;
  }else if (star.stage[n]==-2){	//get track for naked helium star
/*    for (iup=1;iup<nhe-1;iup++){	//find mass index
      if(hestararray[iup].m[0]>star.m[n]){	//hestararray must be sorted from low mass to high mass (done when table is read in)
        break;
      }
    }
//    if (screen) cout << "iup=" << iup << " n_he=" << nhe << " hestararray[iup].m[0]=" << hestararray[iup].m[0] << "Msun" << endl;
    ilow = iup-1;
    //calculate the ratio between the two nearest evolutionary tracks
    mratio = ratio(star.m[n],hestararray[iup].m[0],hestararray[ilow].m[0]);	//calculate ratio between masses
    //set new array length
    star.track.n = star.last+hestararray[ilow].n+hestararray[iup].n-1;
    //reserve memory of new length
    star.track.m = (double *)realloc(star.track.m,star.track.n*sizeof(double));
    star.track.t = (double *)realloc(star.track.t,star.track.n*sizeof(double));
    star.track.r = (double *)realloc(star.track.r,star.track.n*sizeof(double));
    star.track.cm = (double *)realloc(star.track.cm,star.track.n*sizeof(double));
    star.track.llum = (double *)realloc(star.track.llum,star.track.n*sizeof(double));
    star.track.lteff = (double *)realloc(star.track.lteff,star.track.n*sizeof(double));
    star.track.lambda = (double *)realloc(star.track.lambda,star.track.n*sizeof(double));
    star.track.cc = (double *)realloc(star.track.cc,star.track.n*sizeof(double));
    star.track.cf = (double *)realloc(star.track.cf,star.track.n*sizeof(double));
    star.track.omega = (double *)realloc(star.track.omega,star.track.n*sizeof(double));
    //check if memory allocation fails
    if (star.track.m==NULL) cerr << "#Error: memory reallocation failed: star.track.m" << endl;
    if (star.track.t==NULL) cerr << "#Error: memory reallocation failed: star.track.t" << endl;
    if (star.track.r==NULL) cerr << "#Error: memory reallocation failed: star.track.r" << endl;
    if (star.track.cm==NULL) cerr << "#Error: memory reallocation failed: star.track.cm" << endl;
    if (star.track.llum==NULL) cerr << "#Error: memory reallocation failed: star.track.llum" << endl;
    if (star.track.lteff==NULL) cerr << "#Error: memory reallocation failed: star.track.lteff" << endl;
    if (star.track.lambda==NULL) cerr << "#Error: memory reallocation failed: star.track.lambda" << endl;
    if (star.track.cc==NULL) cerr << "#Error: memory reallocation failed: star.track.cc" << endl;
    if (star.track.cf==NULL) cerr << "#Error: memory reallocation failed: star.track.cf" << endl;
    if (star.track.omega==NULL) cerr << "#Error: memory reallocation failed: star.track.omega" << endl;
    //calculate new track
    j = star.last+1;
    star.track.m[j] = hestararray[iup].mw[0]-mratio*(hestararray[iup].mw[0]-hestararray[ilow].mw[0]);
    if (fabs(star.track.m[j]-star.m[n])>accuracy){
      cerr << "#Error: star.m[" << n << "]=" << star.m[n] << "Msun star.track.m[" << j << "]=" << star.track.m[j] << "Msun=" << hestararray[iup].m[0]-mratio*(hestararray[iup].m[0]-hestararray[ilow].m[0]) << "Msun ilow=" << ilow << " iup=" << iup << endl;
      screen = true;
    }
    star.track.t[j] = fmax(hestararray[iup].t[0]-mratio*(hestararray[iup].t[0]-hestararray[ilow].t[0]),0.0);	//minimal age is 0yr
    star.track.r[j] = hestararray[iup].r[0]-mratio*(hestararray[iup].r[0]-hestararray[ilow].r[0]);
    if (star.track.r[j]<rsmax) star.track.r[j] = fmax(star.track.r[j],2.0*G*star.track.m[j]/(c*c));	//minimal radius is the Schwarzschild radius
    star.track.cm[j] = hestararray[iup].cm[0]-mratio*(hestararray[iup].cm[0]-hestararray[ilow].cm[0]);
    if (star.track.cm[j]>star.track.m[j]){	//if core mass is larger than total mass -> take the mean for both
      star.track.cm[j] = 0.5*(star.track.cm[j]+star.track.m[j]);
      star.track.m[j] = star.track.cm[j];
    }
    star.track.llum[j] = hestararray[iup].llum[0]-mratio*(hestararray[iup].llum[0]-hestararray[ilow].llum[0]);
    star.track.lteff[j] = hestararray[iup].lteff[0]-mratio*(hestararray[iup].lteff[0]-hestararray[ilow].lteff[0]);
    star.track.lambda[j] = lambda_const;
    star.track.cc[j] = star.track.cm[j];
    star.track.cf[j] = 0.0;
    star.track.omega[j] = 0.0;
    star.rmax = j;	//update rmax	//if (star.track.r[j]>star.track.r[star.rmax])
    if ((mratio<0)||(mratio>1)){	//values are extrapolated
      star.track.cm[j] = fmax(star.track.cm[j],0.0);	//minimal core mass is 0
      star.track.cc[j] = fmax(star.track.cc[j],0.0);	//minimal carbon core mass is 0
    }
    jlow = 1;
    jup = 1;
    nlow = hestararray[ilow].n-1;
    nup = hestararray[iup].n-1;
    if (debug){
      cerr << endl << "He-track:" << endl;
    }
    for (j=star.last+2;j<star.track.n-1;j++){
      if ((jlow>nlow)||(jup>nup)){
        cerr << "#Error: j=" << j << " jlow=" << jlow << " nlow=" << nlow << " jup=" << jup << " nup=" << nup << endl;
        screen = true;
        break;
      }else if (jlow*nup<jup*nlow){
        //calculate ratio in upper evolutionary track
        ratioup = ratio(jlow*nup,jup*nlow,(jup-1)*nlow);
        //determine values at upper evolutionary track
        mw = hestararray[iup].mw[jup]-ratioup*(hestararray[iup].mw[jup]-hestararray[iup].mw[jup-1]);
        t = hestararray[iup].t[jup]-ratioup*(hestararray[iup].t[jup]-hestararray[iup].t[jup-1]);
        r = hestararray[iup].r[jup]-ratioup*(hestararray[iup].r[jup]-hestararray[iup].r[jup-1]);
        cm = hestararray[iup].cm[jup]-ratioup*(hestararray[iup].cm[jup]-hestararray[iup].cm[jup-1]);
        llum = hestararray[iup].llum[jup]-ratioup*(hestararray[iup].llum[jup]-hestararray[iup].llum[jup-1]);
        lteff = hestararray[iup].lteff[jup]-ratioup*(hestararray[iup].lteff[jup]-hestararray[iup].lteff[jup-1]);
        //interpolate/extrapolate between evolutionary tracks
        star.track.m[j] = fmax(mw-mratio*(mw-hestararray[ilow].mw[jlow]),1.0e-10);	//minimal mass of 10^{-10}Msun
        star.track.t[j] = fmax(t-mratio*(t-hestararray[ilow].t[jlow]),0.0);	//minimal age is 0yr
        star.track.r[j] = r-mratio*(r-hestararray[ilow].r[jlow]);
        star.track.cm[j] = cm-mratio*(cm-hestararray[ilow].cm[jlow]);
        if (star.track.cm[j]>star.track.m[j]){	//if core mass is larger than total mass -> take the mean for both
          star.track.cm[j] = 0.5*(star.track.cm[j]+star.track.m[j]);
          star.track.m[j] = star.track.cm[j];
        }
        star.track.llum[j] = llum-mratio*(llum-hestararray[ilow].llum[jlow]);
        star.track.lteff[j] = lteff-mratio*(lteff-hestararray[ilow].lteff[jlow]);
        star.track.lambda[j] = lambda_const;
        star.track.cc[j] = star.track.cm[j];
        star.track.cf[j] = 0.0;
        star.track.omega[j] = 0.0;
        //next trackpoint
        jlow++;
      }else{
        //calculate ratio in lower evolutionary track
        ratiolow = ratio(jup*nlow,jlow*nup,(jlow-1)*nup);
        //determine values at lower evolutionary track
        mw = hestararray[ilow].mw[jlow]-ratiolow*(hestararray[ilow].mw[jlow]-hestararray[ilow].mw[jlow-1]);
        t = hestararray[ilow].t[jlow]-ratiolow*(hestararray[ilow].t[jlow]-hestararray[ilow].t[jlow-1]);
        r = hestararray[ilow].r[jlow]-ratiolow*(hestararray[ilow].r[jlow]-hestararray[ilow].r[jlow-1]);
        cm = hestararray[ilow].cm[jlow]-ratiolow*(hestararray[ilow].cm[jlow]-hestararray[ilow].cm[jlow-1]);
        llum = hestararray[ilow].llum[jlow]-ratiolow*(hestararray[ilow].llum[jlow]-hestararray[ilow].llum[jlow-1]);
        lteff = hestararray[ilow].lteff[jlow]-ratiolow*(hestararray[ilow].lteff[jlow]-hestararray[ilow].lteff[jlow-1]);
        //interpolate/extrapolate between evolutionary tracks
        star.track.m[j] = fmax(hestararray[iup].mw[jup]-mratio*(hestararray[iup].mw[jup]-mw),1.0e-10);	//minimal mass of 10^{-10}Msun
        star.track.t[j] = fmax(hestararray[iup].t[jup]-mratio*(hestararray[iup].t[jup]-t),0.0);	//minimal age is 0yr
        star.track.r[j] = hestararray[iup].r[jup]-mratio*(hestararray[iup].r[jup]-r);
        star.track.cm[j] = hestararray[iup].cm[jup]-mratio*(hestararray[iup].cm[jup]-cm);
        if (star.track.cm[j]>star.track.m[j]){	//if core mass is larger than total mass -> take the mean for both
          star.track.cm[j] = 0.5*(star.track.cm[j]+star.track.m[j]);
          star.track.m[j] = star.track.cm[j];
        }
        star.track.llum[j] = hestararray[iup].llum[jup]-mratio*(hestararray[iup].llum[jup]-llum);
        star.track.lteff[j] = hestararray[iup].lteff[jup]-mratio*(hestararray[iup].lteff[jup]-lteff);
        star.track.lambda[j] = lambda_const;
        star.track.cc[j] = star.track.cm[j];
        star.track.cf[j] = 0.0;
        star.track.omega[j] = 0.0;
//        if (jlow*nup==jup*nlow){
//          star.track.m[j] = nextafter(star.track.m[j],-1.0);
//          star.track.t[j] = nextafter(star.track.t[j],-1.0);
//          star.track.r[j] = nextafter(star.track.r[j],-1.0);
//          star.track.cm[j] = nextafter(star.track.cm[j],-1.0);
//          star.track.llum[j] = nextafter(star.track.llum[j],-1000.0);
//          star.track.lteff[j] = nextafter(star.track.lteff[j],-1.0);
//          star.track.lambda[j] = nextafter(star.track.lambda[j],-1000.0);
//          star.track.cc[j] = nextafter(star.track.cc[j],-1.0);
//        }
        //next trackpoint
        jup++;
      }
      if (star.track.r[j]<rsmax) star.track.r[j] = fmax(star.track.r[j],2.0*G*star.track.m[j]/(c*c));	//minimal radius is the Schwarzschild radius
      if (star.track.r[j]>star.track.r[star.rmax]) star.rmax = j;	//update rmax
      if ((mratio<0)||(mratio>1)){	//values are extrapolated
        star.track.m[j] = fmin(star.track.m[j],star.track.m[j-1]);	//mass can not increase
        star.track.t[j] = fmax(star.track.t[j],star.track.t[j-1]+1.0);	//age can not decrease
        star.track.cm[j] = fmax(star.track.cm[j],star.track.cm[j-1]);	//core mass can not decrease
        star.track.cc[j] = fmax(star.track.cc[j],star.track.cc[j-1]);	//carbon core mass can not decrease
      }
      if ((star.track.t[j]<star.track.t[j-1])&&(jlow*nup<jup*nlow)){
        cerr << "#Error: star.track.t[" << j << "]=" << star.track.t[j] << "yr star.track.t[" << j-1 << "]=" << star.track.t[j-1] << "yr t=" << t << "yr hestararray[" << ilow << "].t[" << jlow << "]=" << hestararray[ilow].t[jlow] << "yr m_ratio=" << mratio << endl;
        screen = true;
      }
      if ((star.track.t[j]<star.track.t[j-1])&&(jlow*nup>=jup*nlow)){
        cerr << "#Error: star.track.t[" << j << "]=" << star.track.t[j] << "yr star.track.t[" << j-1 << "]=" << star.track.t[j-1] << "yr hestararray[" << iup << "].t[" << jup << "]=" << hestararray[iup].t[jup] << "yr t=" << t << "yr m_ratio=" << mratio << endl;
        screen = true;
      }
      if (debug){
        cerr << "star.track.m[" << j << "]=" << star.track.m[j] << "Msun star.track.t[" << j << "]=" << star.track.t[j] << "yr star.track.r[" << j << "]=" << star.track.r[j] << "Rsun star.track.cm[" << j << "]=" << star.track.cm[j] << "Msun jlow=" << jlow << " nlow=" << nlow << " jup=" << jup << " nup=" << nup << endl;
      }
    }
    if ((jlow!=nlow)||(jup!=nup)) cerr << "#Warning: jlow=" << jlow << " nlow=" << nlow << " ilow=" << ilow << " jup=" << jup << " nup=" << nup << " iup=" << iup << endl;
    j = star.track.n-1;
    star.track.m[j] = fmax(hestararray[iup].mw[nup]-mratio*(hestararray[iup].mw[nup]-hestararray[ilow].mw[nlow]),1.0e-10);	//minimal mass of 10^{-10}Msun
    star.track.t[j] = fmax(hestararray[iup].t[nup]-mratio*(hestararray[iup].t[nup]-hestararray[ilow].t[nlow]),0.0);	//minimal age is 0yr
    star.track.r[j] = hestararray[iup].r[nup]-mratio*(hestararray[iup].r[nup]-hestararray[ilow].r[nlow]);
    if (star.track.r[j]<rsmax) star.track.r[j] = fmax(star.track.r[j],2.0*G*star.track.m[j]/(c*c));	//minimal radius is the Schwarzschild radius
    star.track.cm[j] = hestararray[iup].cm[nup]-mratio*(hestararray[iup].cm[nup]-hestararray[ilow].cm[nlow]);
    if (star.track.cm[j]>star.track.m[j]){	//if core mass is larger than total mass -> take the mean for both
      star.track.cm[j] = 0.5*(star.track.cm[j]+star.track.m[j]);
      star.track.m[j] = star.track.cm[j];
    }
    star.track.llum[j] = hestararray[iup].llum[nup]-mratio*(hestararray[iup].llum[nup]-hestararray[ilow].llum[nlow]);
    star.track.lteff[j] = hestararray[iup].lteff[nup]-mratio*(hestararray[iup].lteff[nup]-hestararray[ilow].lteff[nlow]);
    star.track.lambda[j] = lambda_const;
    star.track.cc[j] = star.track.cm[j];
    star.track.cf[j] = 0.0;
    star.track.omega[j] = 0.0;
    if (star.track.r[j]>star.track.r[star.rmax]) star.rmax = j;	//update rmax
    if ((mratio<0)||(mratio>1)){	//values are extrapolated
      star.track.m[j] = fmin(star.track.m[j],star.track.m[j-1]);	//mass can not increase
      star.track.t[j] = fmax(star.track.t[j],star.track.t[j-1]+1.0);	//age can not decrease
      star.track.cm[j] = fmax(star.track.cm[j],star.track.cm[j-1]);	//core mass can not decrease
      star.track.cc[j] = fmax(star.track.cc[j],star.track.cc[j-1]);	//carbon core mass can not decrease
    }
    if (star.track.t[j]<star.track.t[j-1]){
      cerr << "#Error: star.track.t[" << j << "]=" << star.track.t[j] << "yr star.track.t[" << j-1 << "]=" << star.track.t[j-1] << "yr hestararray[" << iup << "].t[" << jup << "]=" << hestararray[iup].t[jup] << "yr hestararray[" << ilow << "].t[" << nlow << "]=" << hestararray[ilow].t[nlow] << "yr m_ratio=" << mratio << endl;
      screen = true;
    }
    for (j=star.last+1;j<star.track.n;j++){	//recalculate times
      star.track.t[j] += star.t[n];
    }
    star.inimass = star.track.m[star.last];
    star.last++;*/
/*    for (iup=0;iup<nhe1;iup++){	//find mass index
      if(hestararray1[iup].inimass>star.m[n]){	//hestararray1 must be sorted from low mass to high mass (done when table is read in)
        break;
      }
    }
    if (iup==nhe1){
      cerr << "#Error: He-star mass of " << star.m[n] << "Msun not in table" << endl;
      iup = nhe-1;
    }
    ilow = iup-1;
    mratio = ratio(star.m[n],hestararray1[iup].inimass,hestararray1[ilow].inimass);	//calculate ratio between masses
    if (debug) cerr << "star.m[" << n << "]=" << star.m[n] << "Msun hestararray1[" << ilow << "].inimass=" << hestararray1[ilow].inimass << "Msun hestararray1[" << iup << "].inimass=" << hestararray1[iup].inimass << "Msun hestararray1[" << ilow << "].n=" << hestararray1[ilow].n << " hestararray1[" << iup << "].n=" << hestararray1[iup].n << " mratio=" << mratio << endl;
    getnewtrack(hestararray1, ilow, iup, mratio, -1, newtrack, star.rmax);
    star.last++;	//increase track marker
    star.track.m[star.last] = star.m[n-1];	//write old mass to track
    star.track.t[star.last] = star.t[n-1];	//write old time to track
    star.track.r[star.last] = star.r[n-1];	//write old radius to track
    star.track.cm[star.last] = star.cm[n-1];	//write old core mass to track
    star.track.llum[star.last] = star.llum[n-1];	//write old luminosity to track
    star.track.lteff[star.last] = star.lteff[n-1];	//write old effective temperature to track
    star.track.lambda[star.last] = star.lambda[n-1];	//write old lambda to track
    star.track.cc[star.last] = star.cc[n-1];	//write old carbon core mass to track
    star.track.cf[star.last] = star.cf[n-1];	//write old concentration factor to track
    star.track.omega[star.last] = star.omega[n-1];	//write old angular velocity to track
    star.m[n] = newtrack.m[1];	//write new mass
    star.r[n] = newtrack.r[1];	//write new radius
    star.cm[n] = newtrack.cm[1];	//write new core mass
    star.llum[n] = newtrack.llum[1];	//write new luminosity
    star.lteff[n] = newtrack.lteff[1];	//write new effective temperature
    star.lambda[n] = newtrack.lambda[1];	//write new lambda
    star.cc[n] = newtrack.cc[1];	//write new carbon core mass
    star.cf[n] = newtrack.cf[1];	//write new concentration factor
    star.omega[n] = star.track.omega[star.last];	//write new angular velocity
    changetrack(star, n, newtrack);
    //free memory
    free(newtrack.m);
    free(newtrack.t);
    free(newtrack.r);
    free(newtrack.cm);
    free(newtrack.llum);
    free(newtrack.lteff);
    free(newtrack.lambda);
    free(newtrack.cc);
    free(newtrack.cf);
    free(newtrack.omega);*/
    //set stage to naked helium star
    star.stage[n] = 2;
    star.m[n] = star.cm[n];
    star.cm[n] = star.cc[n];
    updatetrack(star,n);
    star.omega[n] = getOmegaMassloss(star,n);
  }else if(star.stage[n]>2){	//star is WD, NS or BH
    //set new array length
//    star.track.n += 1;
    star.track.n = star.last+3;
    //enlarge reserved memory
    star.track.m = (double *)realloc(star.track.m,star.track.n*sizeof(double));
    star.track.t = (double *)realloc(star.track.t,star.track.n*sizeof(double));
    star.track.r = (double *)realloc(star.track.r,star.track.n*sizeof(double));
    star.track.cm = (double *)realloc(star.track.cm,star.track.n*sizeof(double));
    star.track.llum = (double *)realloc(star.track.llum,star.track.n*sizeof(double));
    star.track.lteff = (double *)realloc(star.track.lteff,star.track.n*sizeof(double));
    star.track.lambda = (double *)realloc(star.track.lambda,star.track.n*sizeof(double));
    star.track.cc = (double *)realloc(star.track.cc,star.track.n*sizeof(double));
    star.track.cf = (double *)realloc(star.track.cf,star.track.n*sizeof(double));
    star.track.omega = (double *)realloc(star.track.omega,star.track.n*sizeof(double));
    //check if memory allocation fails
    if (star.track.m==NULL) cerr << "#Error: memory reallocation failed: star.track.m" << endl;
    if (star.track.t==NULL) cerr << "#Error: memory reallocation failed: star.track.t" << endl;
    if (star.track.r==NULL) cerr << "#Error: memory reallocation failed: star.track.r" << endl;
    if (star.track.cm==NULL) cerr << "#Error: memory reallocation failed: star.track.cm" << endl;
    if (star.track.llum==NULL) cerr << "#Error: memory reallocation failed: star.track.llum" << endl;
    if (star.track.lteff==NULL) cerr << "#Error: memory reallocation failed: star.track.lteff" << endl;
    if (star.track.lambda==NULL) cerr << "#Error: memory reallocation failed: star.track.lambda" << endl;
    if (star.track.cc==NULL) cerr << "#Error: memory reallocation failed: star.track.cc" << endl;
    if (star.track.cf==NULL) cerr << "#Error: memory reallocation failed: star.track.cf" << endl;
    if (star.track.omega==NULL) cerr << "#Error: memory reallocation failed: star.track.omega" << endl;
    //set values
    star.last = star.track.n-2;
    star.rmax = star.last;
    star.inimass = star.m[n];
    //set two new trackpoints, first: current values
    star.track.m[star.last] = star.m[n];
    star.track.t[star.last] = star.t[n];
    star.track.r[star.last] = star.r[n];
    star.track.cm[star.last] = star.cm[n];
    star.track.llum[star.last] = star.llum[n];
    star.track.lteff[star.last] = star.lteff[n];
    star.track.lambda[star.last] = star.lambda[n];
    star.track.cc[star.last] = star.cc[n];
    star.track.cf[star.last] = star.cf[n];
    star.track.omega[star.last] = star.omega[n];
    //second trackpoint: stay unchanged
    star.track.m[star.last+1] = star.m[n];
    star.track.t[star.last+1] = 1.0e+99;	//rests for ever
    star.track.r[star.last+1] = star.r[n];
    star.track.cm[star.last+1] = star.cm[n];
    star.track.llum[star.last+1] = star.llum[n];
    star.track.lteff[star.last+1] = star.lteff[n];
    star.track.lambda[star.last+1] = star.lambda[n];
    star.track.cc[star.last+1] = star.cc[n];
    star.track.cf[star.last+1] = star.cf[n];
    star.track.omega[star.last+1] = star.omega[n];
  }else{	//set stage to end stage
    star.last = star.track.n-1;
    if (star.t[n]<star.track.t[star.last]){
      star.m[n] = star.track.m[star.last];
      star.t[n] = star.track.t[star.last];
      star.r[n] = star.track.r[star.last];
      star.cm[n] = star.track.cm[star.last];
      star.llum[n] = star.track.llum[star.last];
      star.lteff[n] = star.track.lteff[star.last];
      star.lambda[n] = star.track.lambda[star.last];
      star.cc[n] = star.track.cc[star.last];
      star.cf[n] = star.track.cf[star.last];
      star.omega[n] = star.track.omega[star.last];
    }
    star.stage[n] = -3;
  }
  if (star.rmax>=star.track.n){	//recalculate rmax if necessary
    star.rmax = star.last;
    for (j=star.last;j<star.track.n;j++){
      if (star.track.r[j]>star.track.r[star.rmax]) star.rmax = j;
    }
  }
  return 1;
}

double radius_ini(double mass){
  int ilow,iup;	//index variables: ilow/iup=mass index
  double mratio,r;	//mass ratio, radius

  for (iup=1;iup<nstar-1;iup++){
    if (stararray[iup].m[0]>mass){
      break;
    }
  }
  ilow = iup-1;
  mratio = ratio(mass,stararray[iup].m[0],stararray[ilow].m[0]);
  r = stararray[iup].r[0]-mratio*(stararray[iup].r[0]-stararray[ilow].r[0]);
  if (r<rsmax) r = fmax(r,2.0*G*mass/(c*c));	//minimal radius is the Schwarzschild radius
  return r;
}

void getnewtrack(t_HRD* stararray0, int ilow, int iup, double mratio, int jpos, t_HRD& newtrack, int& jrmax){
  int j=1,jlow=0,jup=0;	//index variables: j: newtrack, jlow: stararray0[ilow], jup: stararray0[iup]
  int nlow=stararray0[ilow].n-1,nup=stararray0[iup].n-1;	//number of points in stararray0[ilow], stararray0[iup]
  double ratiolow,ratioup;	//mass ratio, ratios at lower/upper mass track
  double m,t,r,cm,llum,lteff,lambda,cc,cf,omega;	//interpolated values of mass, age, radius, core mass, log10(luminosity), log10(effective temperature), lambda, carbon core mass, concentration factor, angular velocity
  double totalup=0.0, totallow=0.0;	//total length in upper/lower track
  double clup=0.0, cllow=0.0, clupold=0.0, cllowold=0.0;	//current length in upper/lower track
  double EbindperG=0.0;	//envelope binding energy for interpolation
  double TAMStime=0.0;	//time for TAMS

  if (debug) cerr << "ilow=" << ilow << " iup=" << iup << " mratio=" << mratio << " jpos=" << jpos << endl;
  if (debug) cerr << "jlow=" << jlow << " jup=" << jup << " nlow=" << nlow << " nup=" << nup << endl;
  //get initial mass of the new track
  newtrack.inimass = stararray0[iup].inimass-mratio*(stararray0[iup].inimass-stararray0[ilow].inimass);
  //get length for the arrays (+maximal values and +t=0 and +reserve)
  newtrack.n = stararray0[ilow].n+stararray0[iup].n+3;
  //initialize TAMS index
  newtrack.TAMS = newtrack.n-1;
  //reserve memory
  newtrack.m = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.t = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.r = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.cm = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.llum = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.lteff = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.lambda = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.cc = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.cf = (double *)malloc(newtrack.n*sizeof(double));
  newtrack.omega = (double *)malloc(newtrack.n*sizeof(double));
  //check if memory allocation fails
  if (newtrack.m==NULL) cerr << "#Error: memory allocation failed: newtrack.m" << endl;
  if (newtrack.t==NULL) cerr << "#Error: memory allocation failed: newtrack.t" << endl;
  if (newtrack.r==NULL) cerr << "#Error: memory allocation failed: newtrack.r" << endl;
  if (newtrack.cm==NULL) cerr << "#Error: memory allocation failed: newtrack.cm" << endl;
  if (newtrack.llum==NULL) cerr << "#Error: memory allocation failed: newtrack.llum" << endl;
  if (newtrack.lteff==NULL) cerr << "#Error: memory allocation failed: newtrack.lteff" << endl;
  if (newtrack.lambda==NULL) cerr << "#Error: memory allocation failed: newtrack.lambda" << endl;
  if (newtrack.cc==NULL) cerr << "#Error: memory allocation failed: newtrack.cc" << endl;
  if (newtrack.cf==NULL) cerr << "#Error: memory allocation failed: newtrack.cf" << endl;
  if (newtrack.omega==NULL) cerr << "#Error: memory allocation failed: newtrack.omega" << endl;
  //calculate new evolutionary track
  //set initial values for first track point
  newtrack.m[1] = stararray0[iup].m[jup]-mratio*(stararray0[iup].m[jup]-stararray0[ilow].m[jlow]);
  newtrack.t[1] = stararray0[iup].t[jup]-mratio*(stararray0[iup].t[jup]-stararray0[ilow].t[jlow]);
  newtrack.r[1] = stararray0[iup].r[jup]-mratio*(stararray0[iup].r[jup]-stararray0[ilow].r[jlow]);
  newtrack.cm[1] = stararray0[iup].cm[jup]-mratio*(stararray0[iup].cm[jup]-stararray0[ilow].cm[jlow]);
  if (debug) cerr << "newtrack.cm[1]=" << newtrack.cm[1] << "=" << stararray0[iup].cm[jup] << "-" << mratio << "*(" << stararray0[iup].cm[jup] << "-" << stararray0[ilow].cm[jlow] << ")" << "=" << stararray0[iup].cm[jup] << "-" << mratio << "*(" << stararray0[iup].cm[jup]-stararray0[ilow].cm[jlow] << ")" << endl;
  newtrack.llum[1] = stararray0[iup].llum[jup]-mratio*(stararray0[iup].llum[jup]-stararray0[ilow].llum[jlow]);
  newtrack.lteff[1] = stararray0[iup].lteff[jup]-mratio*(stararray0[iup].lteff[jup]-stararray0[ilow].lteff[jlow]);
//  newtrack.lambda[1] = stararray0[iup].lambda[jup]-mratio*(stararray0[iup].lambda[jup]-stararray0[ilow].lambda[jlow]);
  EbindperG = stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
              -mratio*(stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
                      -stararray0[ilow].m[jlow]*(stararray0[ilow].m[jlow]-stararray0[ilow].cm[jlow])/(stararray0[ilow].lambda[jlow]*stararray0[ilow].r[jlow]));	//interpolate current binding energy of the envelope
  if ((newtrack.m[1]-newtrack.cm[1]==0.0)&&(EbindperG==0.0)) newtrack.lambda[1] = ignorev;	//no interpolation
  else newtrack.lambda[1] = newtrack.m[1]*(newtrack.m[1]-newtrack.cm[1])/(EbindperG*newtrack.r[1]);	//interpolate lambda	// E_bind=G*M*M_env / lambda*R
  if(isnan(newtrack.lambda[1])){ screen = true; cerr << endl << "#newtrack.lambda[1]=" << newtrack.lambda[1];}
  newtrack.cc[1] = stararray0[iup].cc[jup]-mratio*(stararray0[iup].cc[jup]-stararray0[ilow].cc[jlow]);
  newtrack.cf[1] = stararray0[iup].cf[jup]-mratio*(stararray0[iup].cf[jup]-stararray0[ilow].cf[jlow]);
  newtrack.omega[1] = stararray0[iup].omega[jup]-mratio*(stararray0[iup].omega[jup]-stararray0[ilow].omega[jlow]);
  //set initial values for maximal values
  jrmax = 1;
  newtrack.m[0] = newtrack.m[1];
  newtrack.t[0] = stararray0[iup].t[nup]-mratio*(stararray0[iup].t[nup]-stararray0[ilow].t[nlow]);
  newtrack.r[0] = newtrack.r[1];
  newtrack.cm[0] = stararray0[iup].cm[nup]-mratio*(stararray0[iup].cm[nup]-stararray0[ilow].cm[nlow]);
  newtrack.llum[0] = stararray0[iup].llum[nup]-mratio*(stararray0[iup].llum[nup]-stararray0[ilow].llum[nlow]);
  newtrack.lteff[0] = newtrack.lteff[1];
  newtrack.lambda[0] = newtrack.lambda[1];
  newtrack.cc[0] = stararray0[iup].cc[nup]-mratio*(stararray0[iup].cc[nup]-stararray0[ilow].cc[nlow]);
  newtrack.cf[0] = newtrack.cf[1];
  newtrack.omega[0] = newtrack.omega[1];
/*  for (jup=1;jup<=nup;jup++){
    totalup += fabs(stararray0[iup].t[jup-1]-stararray0[iup].t[jup]);
  }
  for (jlow=1;jlow<=nlow;jlow++){
    totallow += fabs(stararray0[ilow].t[jlow-1]-stararray0[ilow].t[jlow]);
  }*/
  totalup = stararray0[iup].t[nup];
  totallow = stararray0[ilow].t[nlow];
  if (debug){
    cerr << "totalup=" << totalup << " totallow=" << totallow << endl;
    cerr << "newtrack: newtrack.n=" << newtrack.n << endl << " j jlow jup\tmass(Msun)\tage(yr)\tradius(Rsun)\tcore mass(Msun)\tlg(lum/Lsun)\tlg(T_eff/K)\tlambda\tcarbon core(Msun)\tconc. factor\tomega (1/yr)" << endl;	//write head of track-table
    cerr << " " << j << "    0   0\t" << newtrack.m[j] << "\t" << newtrack.t[j] << "\t" << newtrack.r[j] << "\t" << newtrack.cm[j] << "\t" << newtrack.llum[j] << "\t" << newtrack.lteff[j] << "\t" << newtrack.lambda[j] << "\t" << newtrack.cc[j] << "\t" << newtrack.cf[j] << "\t" << newtrack.omega[j] << "\t" << endl;	//write track-table of the star
  }
  jup = 1;
  jlow = 1;
/*  clup = fabs(stararray0[iup].t[jup-1]-stararray0[iup].t[jup]);
  cllow = fabs(stararray0[ilow].t[jlow-1]-stararray0[ilow].t[jlow]);*/
  clup = stararray0[iup].t[jup];
  cllow = stararray0[ilow].t[jlow];
  //estimate TAMS time
  ratioup = ratio(stararray0[iup].t[stararray0[iup].TAMS],totalup,stararray0[iup].t[1]);
  ratiolow = ratio(stararray0[ilow].t[stararray0[ilow].TAMS],totallow,stararray0[ilow].t[1]);
  TAMStime = ratioup-mratio*(ratioup-ratiolow);
  if (debug) cerr << "ratioup=" << ratioup << " ratiolow=" << ratiolow << " ratioTAMS=" << TAMStime;
  TAMStime = newtrack.t[0]-TAMStime*(newtrack.t[0]-newtrack.t[1]);
  if (debug) cerr << " TAMStime=" << TAMStime << "yr newtrack.t[1]=" << newtrack.t[1] << "yr newtrack.t[0]=" << newtrack.t[0] << "yr" << endl;
  //get values for all times in the two nearest evolutionary tracks
  for (j=2;j<newtrack.n;j++){
    if (newtrack.TAMS>j){
      if ((newtrack.t[j-1]>=TAMStime)||((fabs(log10(stararray0[iup].r[jup])-log10(stararray0[ilow].r[jlow]))>1.0)&&(stararray0[ilow].inimass>0.08))){
//      if ((jup>=stararray0[iup].TAMS)&&(jlow>=stararray0[ilow].TAMS)){
        newtrack.TAMS = j-1;	//set TAMS index
        if (debug) cerr << "newtrack.TAMS=" << newtrack.TAMS << " stararray0[" << iup << "].TAMS=" << stararray0[iup].TAMS << " stararray0[" << ilow << "].TAMS=" << stararray0[ilow].TAMS << " next jlow=" << jlow << " nlow=" << nlow << " next jup=" << jup << " nup=" << nup << endl;
      }
    }
    if ((jlow==stararray0[ilow].n)&&(jup==stararray0[iup].n)){	//both evolutionary tracks reached last time
      //set new array length
      newtrack.n = j;
      //remove not needed array part
      newtrack.m = (double *)realloc(newtrack.m,newtrack.n*sizeof(double));
      newtrack.t = (double *)realloc(newtrack.t,newtrack.n*sizeof(double));
      newtrack.r = (double *)realloc(newtrack.r,newtrack.n*sizeof(double));
      newtrack.cm = (double *)realloc(newtrack.cm,newtrack.n*sizeof(double));
      newtrack.llum = (double *)realloc(newtrack.llum,newtrack.n*sizeof(double));
      newtrack.lteff = (double *)realloc(newtrack.lteff,newtrack.n*sizeof(double));
      newtrack.lambda = (double *)realloc(newtrack.lambda,newtrack.n*sizeof(double));
      newtrack.cc = (double *)realloc(newtrack.cc,newtrack.n*sizeof(double));
      newtrack.cf = (double *)realloc(newtrack.cf,newtrack.n*sizeof(double));
      newtrack.omega = (double *)realloc(newtrack.omega,newtrack.n*sizeof(double));
      //check if memory allocation fails
      if (newtrack.m==NULL) cerr << "#Error: memory reallocation failed: newtrack.m" << endl;
      if (newtrack.t==NULL) cerr << "#Error: memory reallocation failed: newtrack.t" << endl;
      if (newtrack.r==NULL) cerr << "#Error: memory reallocation failed: newtrack.r" << endl;
      if (newtrack.cm==NULL) cerr << "#Error: memory reallocation failed: newtrack.cm" << endl;
      if (newtrack.llum==NULL) cerr << "#Error: memory reallocation failed: newtrack.llum" << endl;
      if (newtrack.lteff==NULL) cerr << "#Error: memory reallocation failed: newtrack.lteff" << endl;
      if (newtrack.lambda==NULL) cerr << "#Error: memory reallocation failed: newtrack.lambda" << endl;
      if (newtrack.cc==NULL) cerr << "#Error: memory reallocation failed: newtrack.cc" << endl;
      if (newtrack.cf==NULL) cerr << "#Error: memory reallocation failed: newtrack.cf" << endl;
      if (newtrack.omega==NULL) cerr << "#Error: memory reallocation failed: newtrack.omega" << endl;
      //last time is maximal time
      newtrack.t[0] = newtrack.t[j-1];
      break;
    }else if (cllow*totalup==clup*totallow){	//same number of steps in the nearest two evolutionary tracks or interpolate at same point
      //interpolate between evolutionary tracks
      newtrack.m[j] = fmax(stararray0[iup].m[jup]-mratio*(stararray0[iup].m[jup]-stararray0[ilow].m[jlow]),1.0e-10);	//new mass: minimal mass of 10^{-10}Msun
//      if (debug) cerr << newtrack.m[j] << "=" << stararray0[iup].m[jup] << "-" << mratio << "*(" << stararray0[iup].m[jup] << "-" << stararray0[ilow].m[jlow] << ")" << endl;
      newtrack.t[j] = fmax(stararray0[iup].t[jup]-mratio*(stararray0[iup].t[jup]-stararray0[ilow].t[jlow]),newtrack.t[j-1]*(1.0+0.01*accuracy));	//new age: minimal age is previous age + 0.01*accuracy*age
      newtrack.r[j] = stararray0[iup].r[jup]-mratio*(stararray0[iup].r[jup]-stararray0[ilow].r[jlow]);	//new radius
      if (newtrack.r[j]<rsmax) newtrack.r[j] = fmax(newtrack.r[j],newtrack.m[j]*rsmax*0.001);	//minimal radius is the Schwarzschild radius
      newtrack.cm[j] = stararray0[iup].cm[jup]-mratio*(stararray0[iup].cm[jup]-stararray0[ilow].cm[jlow]);	//new core mass
//      if (debug) cerr << newtrack.cm[j] << "=" << stararray0[iup].cm[jup] << "-" << mratio << "*(" << stararray0[iup].cm[jup] << "-" << stararray0[ilow].cm[jlow] << ")" << endl;
      if (newtrack.cm[j]>newtrack.m[j]){	//if core mass is larger than total mass -> take the mean for both
        newtrack.cm[j] = 0.5*(newtrack.cm[j]+newtrack.m[j]);
        newtrack.m[j] = newtrack.cm[j];
      }
      newtrack.llum[j] = stararray0[iup].llum[jup]-mratio*(stararray0[iup].llum[jup]-stararray0[ilow].llum[jlow]);	//new luminosity
      newtrack.lteff[j] = stararray0[iup].lteff[jup]-mratio*(stararray0[iup].lteff[jup]-stararray0[ilow].lteff[jlow]);	//new effective temperature
//      newtrack.lambda[j] = stararray0[iup].lambda[jup]-mratio*(stararray0[iup].lambda[jup]-stararray0[ilow].lambda[jlow]);	//new lambda
      EbindperG = stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
                  -mratio*(stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
                          -stararray0[ilow].m[jlow]*(stararray0[ilow].m[jlow]-stararray0[ilow].cm[jlow])/(stararray0[ilow].lambda[jlow]*stararray0[ilow].r[jlow]));	//new binding energy
      if ((newtrack.m[j]-newtrack.cm[j]==0.0)&&(EbindperG==0.0)) newtrack.lambda[j] = ignorev;	//no interpolation
      else newtrack.lambda[j] = newtrack.m[j]*(newtrack.m[j]-newtrack.cm[j])/(EbindperG*newtrack.r[j]);	//new lambda	// E_bind=G*M*M_env / lambda*R
      if(isnan(newtrack.lambda[j])){ screen = true; cerr << endl << "#newtrack.lambda[j]=" << newtrack.lambda[j];}
      newtrack.cc[j] = stararray0[iup].cc[jup]-mratio*(stararray0[iup].cc[jup]-stararray0[ilow].cc[jlow]);	//new core mass
      newtrack.cf[j] = stararray0[iup].cf[jup]-mratio*(stararray0[iup].cf[jup]-stararray0[ilow].cf[jlow]);	//new concentration factor	//###
      newtrack.omega[j] = stararray0[iup].omega[jup]-mratio*(stararray0[iup].omega[jup]-stararray0[ilow].omega[jlow]);	//new angular velocity
      if (debug) cerr << j << "   "<< jlow << "  "<< jup << "\t" << newtrack.m[j] << "\t" << newtrack.t[j] << "\t" << newtrack.r[j] << "\t" << newtrack.cm[j] << "\t" << newtrack.llum[j] << "\t" << newtrack.lteff[j] << "\t" << newtrack.lambda[j] << "\t" << newtrack.cc[j] << "\t" << newtrack.cf[j] << "\t" << newtrack.omega[j] << endl << "newtrack.cf[" << j << "]=" << newtrack.cf[j] << "=" << stararray0[iup].cf[jup] << "-" << mratio << "*(" << stararray0[iup].cf[jup] << "-" << stararray0[ilow].cf[jlow] << ")" << "=stararray0[" << iup << "].cf[" << jup << "]-mratio*(stararray0[" << iup << "].cf[" << jup << "]-stararray0[" << ilow << "].cf[" << jlow << "])" << endl;	//write track-table of the star	// << "\t" << cllowold << "\t" << cllow << "\t" << clupold << "\t" << clup
      //increase index for upper & lower evolutionary track
      jup++;
      jlow++;
      //save old values
      clupold = clup;
      cllowold = cllow;
      //add length
/*      clup += fabs(stararray0[iup].t[jup-1]-stararray0[iup].t[jup]);
      cllow += fabs(stararray0[ilow].t[jlow-1]-stararray0[ilow].t[jlow]);*/
      clup = stararray0[iup].t[jup];
      cllow = stararray0[ilow].t[jlow];
    }else if (cllow*totalup<clup*totallow){
      //calculate ratio in upper evolutionary track
      ratioup = ratio(totalup*cllow/totallow,clup,clupold);
      //determine values at upper evolutionary track
      m = stararray0[iup].m[jup]-ratioup*(stararray0[iup].m[jup]-stararray0[iup].m[jup-1]);	//get mass at upper evolutionary track
      t = stararray0[iup].t[jup]-ratioup*(stararray0[iup].t[jup]-stararray0[iup].t[jup-1]);	//get age at upper evolutionary track
      r = stararray0[iup].r[jup]-ratioup*(stararray0[iup].r[jup]-stararray0[iup].r[jup-1]);	//get radius at upper evolutionary track
      cm = stararray0[iup].cm[jup]-ratioup*(stararray0[iup].cm[jup]-stararray0[iup].cm[jup-1]);	//get core mass at upper evolutionary track
      llum = stararray0[iup].llum[jup]-ratioup*(stararray0[iup].llum[jup]-stararray0[iup].llum[jup-1]);	//get luminosity at upper evolutionary track
      lteff = stararray0[iup].lteff[jup]-ratioup*(stararray0[iup].lteff[jup]-stararray0[iup].lteff[jup-1]);	//get effective temperature at upper evolutionary track
//      lambda = stararray0[iup].lambda[jup]-ratioup*(stararray0[iup].lambda[jup]-stararray0[iup].lambda[jup-1]);	//get lambda at upper evolutionary track
      EbindperG = stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
                  -ratioup*(stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
                           -stararray0[iup].m[jup-1]*(stararray0[iup].m[jup-1]-stararray0[iup].cm[jup-1])/(stararray0[iup].lambda[jup-1]*stararray0[iup].r[jup-1]));	//new binding energy at upper evolutionary track
      if ((m-cm==0.0)&&(EbindperG==0.0)) lambda = ignorev;	//no interpolation
      else lambda = m*(m-cm)/(EbindperG*r);	//get lambda at upper evolutionary track	// E_bind=G*M*M_env / lambda*R
      if(isnan(lambda)||(lambda==0.0)){ screen = true; cerr << endl << "#lambda_upper=" << lambda << " m=" << m << "Msun cm=" << cm << "Msun iup=" << iup << " jup=" << jup << endl << "EbindperG=" << EbindperG << "Msun^2/Rsun=" << stararray0[iup].m[jup] << "*(" << stararray0[iup].m[jup] << "-" << stararray0[iup].cm[jup] << ")/(" << stararray0[iup].lambda[jup] << "*" << stararray0[iup].r[jup] << ")-" << ratioup << "*(" << stararray0[iup].m[jup] << "*(" << stararray0[iup].m[jup] << "-" << stararray0[iup].cm[jup] << ")/(" << stararray0[iup].lambda[jup] << "*" << stararray0[iup].r[jup] << ")-" << stararray0[iup].m[jup-1] << "*(" << stararray0[iup].m[jup-1] << "-" << stararray0[iup].cm[jup-1] << ")/(" << stararray0[iup].lambda[jup-1] << "*" << stararray0[iup].r[jup-1] << "))Msun^2/Rsun";}
      cc = stararray0[iup].cc[jup]-ratioup*(stararray0[iup].cc[jup]-stararray0[iup].cc[jup-1]);	//get carbon core mass at upper evolutionary track
      cf = stararray0[iup].cf[jup]-ratioup*(stararray0[iup].cf[jup]-stararray0[iup].cf[jup-1]);	//get concentration factor at upper evolutionary track
      omega = stararray0[iup].omega[jup]-ratioup*(stararray0[iup].omega[jup]-stararray0[iup].omega[jup-1]);	//get angular velocity at upper evolutionary track
      //interpolate/extrapolate between evolutionary tracks
      newtrack.m[j] = fmax(m-mratio*(m-stararray0[ilow].m[jlow]),1.0e-10);	//new mass: minimal mass of 10^{-10}Msun
      newtrack.t[j] = fmax(t-mratio*(t-stararray0[ilow].t[jlow]),newtrack.t[j-1]*(1.0+0.01*accuracy));	//new age: minimal age is previous age + 0.01*accuracy*age
      newtrack.r[j] = r-mratio*(r-stararray0[ilow].r[jlow]);	//new radius
      if (newtrack.r[j]<rsmax) newtrack.r[j] = fmax(newtrack.r[j],newtrack.m[j]*rsmax*0.001);	//minimal radius is the Schwarzschild radius
      newtrack.cm[j] = cm-mratio*(cm-stararray0[ilow].cm[jlow]);	//new core mass
      if (newtrack.cm[j]>newtrack.m[j]){	//if core mass is larger than total mass -> take the mean for both
        newtrack.cm[j] = 0.5*(newtrack.cm[j]+newtrack.m[j]);
        newtrack.m[j] = newtrack.cm[j];
      }
      newtrack.llum[j] = llum-mratio*(llum-stararray0[ilow].llum[jlow]);	//new luminosity
      newtrack.lteff[j] = lteff-mratio*(lteff-stararray0[ilow].lteff[jlow]);	//new effective temperature
//      newtrack.lambda[j] = lambda-mratio*(lambda-stararray0[ilow].lambda[jlow]);	//new lambda
      EbindperG = m*(m-cm)/(lambda*r)-mratio*(m*(m-cm)/(lambda*r)
                                             -stararray0[ilow].m[jlow]*(stararray0[ilow].m[jlow]-stararray0[ilow].cm[jlow])/(stararray0[ilow].lambda[jlow]*stararray0[ilow].r[jlow]));	//new binding energy
      if ((newtrack.m[j]-newtrack.cm[j]==0.0)) newtrack.lambda[j] = ignorev;	//no interpolation	//&&(EbindperG==0.0)
      else newtrack.lambda[j] = newtrack.m[j]*(newtrack.m[j]-newtrack.cm[j])/(EbindperG*newtrack.r[j]);	//get lambda	// E_bind=G*M*M_env / lambda*R
      if(isnan(newtrack.lambda[j])||(newtrack.lambda[j]==0.0)){ screen = true; cerr << endl << "#newtrack.lambda[j]_upper=" << newtrack.lambda[j] << " newtrack.m[j]=" << newtrack.m[j] << "Msun newtrack.cm[j]=" << newtrack.cm[j] << "Msun ilow=" << ilow << " jlow=" << jlow << endl << "EbindperG=" << EbindperG << "Msun^2/Rsun=" << m << "*(" << m << "-" << cm << ")/(" << lambda << "*" << r << ")-" << mratio << "*(" << m << "*(" << m << "-" << cm << ")/(" << lambda << "*" << r << ")-" << stararray0[ilow].m[jlow] << "*(" << stararray0[ilow].m[jlow] << "-" << stararray0[ilow].cm[jlow] << ")/(" << stararray0[ilow].lambda[jlow] << "*" << stararray0[ilow].r[jlow] << "))Msun^2/Rsun";}
      newtrack.cc[j] = cc-mratio*(cc-stararray0[ilow].cc[jlow]);	//new carbon core mass
      newtrack.cf[j] = cf-mratio*(cf-stararray0[ilow].cf[jlow]);	//new concentration factor
      newtrack.omega[j] = omega-mratio*(omega-stararray0[ilow].omega[jlow]);	//new angular velocity
      if (debug) cerr << j << "   "<< jlow << "  "<< jup-ratioup << "\t" << newtrack.m[j] << "\t" << newtrack.t[j] << "\t" << newtrack.r[j] << "\t" << newtrack.cm[j] << "\t" << newtrack.llum[j] << "\t" << newtrack.lteff[j] << "\t" << newtrack.lambda[j] << "\t" << newtrack.cc[j] << "\t" << newtrack.cf[j] << "\t" << newtrack.omega[j] << endl;	//write track-table of the star	// << "\t" << cllowold << "\t" << cllow << "\t" << clupold << "\t" << clup
      //next trackpoint
      jlow++;
      //save old values
      cllowold = cllow;
      //add length
/*      cllow += fabs(stararray0[ilow].t[jlow-1]-stararray0[ilow].t[jlow]);*/
      cllow = stararray0[ilow].t[jlow];
    }else if (cllow*totalup>clup*totallow){
      //calculate ratio in lower evolutionary track
      ratiolow = ratio(totallow*clup/totalup,cllow,cllowold);
      //determine values at lower evolutionary track
      m = stararray0[ilow].m[jlow]-ratiolow*(stararray0[ilow].m[jlow]-stararray0[ilow].m[jlow-1]);	//get mass at lower evolutionary track
      t = stararray0[ilow].t[jlow]-ratiolow*(stararray0[ilow].t[jlow]-stararray0[ilow].t[jlow-1]);	//get age at lower evolutionary track
      r = stararray0[ilow].r[jlow]-ratiolow*(stararray0[ilow].r[jlow]-stararray0[ilow].r[jlow-1]);	//get radius at lower evolutionary track
      cm = stararray0[ilow].cm[jlow]-ratiolow*(stararray0[ilow].cm[jlow]-stararray0[ilow].cm[jlow-1]);	//get core mass at lower evolutionary track
      llum = stararray0[ilow].llum[jlow]-ratiolow*(stararray0[ilow].llum[jlow]-stararray0[ilow].llum[jlow-1]);	//get luminosity at lower evolutionary track
      lteff = stararray0[ilow].lteff[jlow]-ratiolow*(stararray0[ilow].lteff[jlow]-stararray0[ilow].lteff[jlow-1]);	//get effective temperature at lower evolutionary track
//      lambda = stararray0[ilow].lambda[jlow]-ratiolow*(stararray0[ilow].lambda[jlow]-stararray0[ilow].lambda[jlow-1]);	//get lambda at lower evolutionary track
      EbindperG = stararray0[ilow].m[jlow]*(stararray0[ilow].m[jlow]-stararray0[ilow].cm[jlow])/(stararray0[ilow].lambda[jlow]*stararray0[ilow].r[jlow])
                  -ratiolow*(stararray0[ilow].m[jlow]*(stararray0[ilow].m[jlow]-stararray0[ilow].cm[jlow])/(stararray0[ilow].lambda[jlow]*stararray0[ilow].r[jlow])
                            -stararray0[ilow].m[jlow-1]*(stararray0[ilow].m[jlow-1]-stararray0[ilow].cm[jlow-1])/(stararray0[ilow].lambda[jlow-1]*stararray0[ilow].r[jlow-1]));	//new binding energy at lower evolutionary track
      if ((m-cm==0.0)&&(EbindperG==0.0)) lambda = ignorev;	//no interpolation
      else lambda = m*(m-cm)/(EbindperG*r);	//get lambda at lower evolutionary track	// E_bind=G*M*M_env / lambda*R
      if(isnan(lambda)||(lambda==0.0)){ screen = true; cerr << endl << "#lambda_lower=" << lambda << " m=" << m << "Msun cm=" << cm << "Msun ilow=" << ilow << " jlow=" << jlow << endl << "EbindperG=" << EbindperG << "Msun^2/Rsun=" << stararray0[ilow].m[jlow] << "*(" << stararray0[ilow].m[jlow] << "-" << stararray0[ilow].cm[jlow] << ")/(" << stararray0[ilow].lambda[jlow] << "*" << stararray0[ilow].r[jlow] << ")-" << ratioup << "*(" << stararray0[ilow].m[jlow] << "*(" << stararray0[ilow].m[jlow] << "-" << stararray0[ilow].cm[jlow] << ")/(" << stararray0[ilow].lambda[jlow] << "*" << stararray0[ilow].r[jlow] << ")-" << stararray0[ilow].m[jlow-1] << "*(" << stararray0[ilow].m[jlow-1] << "-" << stararray0[ilow].cm[jlow-1] << ")/(" << stararray0[ilow].lambda[jlow-1] << "*" << stararray0[ilow].r[jlow-1] << "))Msun^2/Rsun";}
      cc = stararray0[ilow].cc[jlow]-ratiolow*(stararray0[ilow].cc[jlow]-stararray0[ilow].cc[jlow-1]);	//get carbon core mass at lower evolutionary track
      cf = stararray0[ilow].cf[jlow]-ratiolow*(stararray0[ilow].cf[jlow]-stararray0[ilow].cf[jlow-1]);	//get concentration factor at lower evolutionary track
      omega = stararray0[ilow].omega[jlow]-ratiolow*(stararray0[ilow].omega[jlow]-stararray0[ilow].omega[jlow-1]);	//get angular velocity at lower evolutionary track
      //interpolate/extrapolate between evolutionary tracks
      newtrack.m[j] = fmax(stararray0[iup].m[jup]-mratio*(stararray0[iup].m[jup]-m),1.0e-10);	//new mass: minimal mass of 10^{-10}Msun
      newtrack.t[j] = fmax(stararray0[iup].t[jup]-mratio*(stararray0[iup].t[jup]-t),newtrack.t[j-1]*(1.0+0.01*accuracy));	//new age: minimal age is previous age + 0.01*accuracy*age
      newtrack.r[j] = stararray0[iup].r[jup]-mratio*(stararray0[iup].r[jup]-r);	//new radius
      if (newtrack.r[j]<rsmax) newtrack.r[j] = fmax(newtrack.r[j],newtrack.m[j]*rsmax*0.001);	//minimal radius is the Schwarzschild radius
      newtrack.cm[j] = stararray0[iup].cm[jup]-mratio*(stararray0[iup].cm[jup]-cm);	//new core mass
      if (newtrack.cm[j]>newtrack.m[j]){	//if core mass is larger than total mass -> take the mean for both
        newtrack.cm[j] = 0.5*(newtrack.cm[j]+newtrack.m[j]);
        newtrack.m[j] = newtrack.cm[j];
      }
      newtrack.llum[j] = stararray0[iup].llum[jup]-mratio*(stararray0[iup].llum[jup]-llum);	//new luminosity
      newtrack.lteff[j] = stararray0[iup].lteff[jup]-mratio*(stararray0[iup].lteff[jup]-lteff);	//new effective temperature
//      newtrack.lambda[j] = stararray0[iup].lambda[jup]-mratio*(stararray0[iup].lambda[jup]-lambda);	//new lambda
      EbindperG = stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
                  -mratio*(stararray0[iup].m[jup]*(stararray0[iup].m[jup]-stararray0[iup].cm[jup])/(stararray0[iup].lambda[jup]*stararray0[iup].r[jup])
                          -m*(m-cm)/(lambda*r));	//new binding energy
      if ((newtrack.m[j]-newtrack.cm[j]==0.0)) newtrack.lambda[j] = ignorev;	//no interpolation	//&&(EbindperG==0.0)
      else newtrack.lambda[j] = newtrack.m[j]*(newtrack.m[j]-newtrack.cm[j])/(EbindperG*newtrack.r[j]);	//get lambda	// E_bind=G*M*M_env / lambda*R
      if(isnan(newtrack.lambda[j])||(newtrack.lambda[j]==0.0)){ screen = true; cerr << endl << "#newtrack.lambda[j]_lower=" << newtrack.lambda[j] << " newtrack.m[j]=" << newtrack.m[j] << "Msun newtrack.cm[j]=" << newtrack.cm[j] << "Msun iup=" << iup << " jup=" << jup << endl << "EbindperG=" << EbindperG << "Msun^2/Rsun=" << stararray0[iup].m[jup] << "*(" << stararray0[iup].m[jup] << "-" << stararray0[iup].cm[jup] << ")/(" << stararray0[iup].lambda[jup] << "*" << stararray0[iup].r[jup] << ")-" << mratio << "*(" << stararray0[iup].m[jup] << "*(" << stararray0[iup].m[jup] << "-" << stararray0[iup].cm[jup] << ")/(" << stararray0[iup].lambda[jup] << "*" << stararray0[iup].r[jup] << ")-" << m << "*(" << m << "-" << cm << ")/(" << lambda << "*" << r << "))Msun^2/Rsun";}
      newtrack.cc[j] = stararray0[iup].cc[jup]-mratio*(stararray0[iup].cc[jup]-cc);	//new carbon core mass
      newtrack.cf[j] = stararray0[iup].cf[jup]-mratio*(stararray0[iup].cf[jup]-cf);	//new concentration factor
      newtrack.omega[j] = stararray0[iup].omega[jup]-mratio*(stararray0[iup].omega[jup]-omega);	//new angular velocity
      if (debug) cerr << j << "   "<< jlow-ratiolow << "  "<< jup << "\t" << newtrack.m[j] << "\t" << newtrack.t[j] << "\t" << newtrack.r[j] << "\t" << newtrack.cm[j] << "\t" << newtrack.llum[j] << "\t" << newtrack.lteff[j] << "\t" << newtrack.lambda[j] << "\t" << newtrack.cc[j] << "\t" << newtrack.cf[j] << "\t" << newtrack.omega[j] << endl;	//write track-table of the star	// << "\t" << cllowold << "\t" << cllow << "\t" << clupold << "\t" << clup
      //next trackpoint
      jup++;
      //save old values
      clupold = clup;
      //add length
/*      clup += fabs(stararray0[iup].t[jup-1]-stararray0[iup].t[jup]);*/
      clup = stararray0[iup].t[jup];
    }else{
      cerr << "#Error: table length: stararray0[" << ilow << "].n=" << stararray0[ilow].n << " stararray0[" << iup << "].n=" << stararray0[iup].n << endl << " Error: j=" << j << " jlow=" << jlow << " nlow=" << nlow << " jup=" << jup << " nup=" << nup << endl;
      screen = true;
    }
    //determine maximal values
    if (newtrack.m[0]<newtrack.m[j]) newtrack.m[0] = newtrack.m[j];
    if (newtrack.r[0]<newtrack.r[j]){
      newtrack.r[0] = newtrack.r[j];
      if (j<=newtrack.TAMS) jrmax = j;
    }
    if (newtrack.cm[0]<newtrack.cm[j]) newtrack.cm[0] = newtrack.cm[j];
    if (newtrack.llum[0]<newtrack.llum[j]) newtrack.llum[0] = newtrack.llum[j];
    if (newtrack.lteff[0]<newtrack.lteff[j]) newtrack.lteff[0] = newtrack.lteff[j];
    if (newtrack.lambda[0]<newtrack.lambda[j]) newtrack.lambda[0] = newtrack.lambda[j];
    if (newtrack.cc[0]<newtrack.cc[j]) newtrack.cc[0] = newtrack.cc[j];
    if (newtrack.cf[0]<newtrack.cf[j]) newtrack.cf[0] = newtrack.cf[j];
    if (newtrack.omega[0]<newtrack.omega[j]) newtrack.omega[0] = newtrack.omega[j];
    //check if maximal time for the star is reached
    if (newtrack.t[0]*(1.0+accuracy)<newtrack.t[j]){	//set end track conditions
      cerr << "#Warning: end track time newtrack.t[0]=" << newtrack.t[0] << "yr=" << stararray0[iup].t[nup]-mratio*(stararray0[iup].t[nup]-stararray0[ilow].t[nlow]) << "yr reached: newtrack.t[" << j << "]=" << newtrack.t[j] << "yr newtrack.t[0]-newtrack.t[" << j << "]=" << newtrack.t[0]-newtrack.t[j] << "yr" << endl;
/*      if (difflow==diffup) cerr << " Warning: newtrack.t[" << j << "]=" << newtrack.t[j] << "yr=" << stararray0[iup].t[jup-1]-mratio*(stararray0[iup].t[jup-1]-stararray0[ilow].t[jlow-1]) << endl;
      else if (difflow<diffup) cerr << " Warning: newtrack.t[" << j << "]=" << newtrack.t[j] << "yr=" << stararray0[iup].t[jup]-ratioup*(stararray0[iup].t[jup]-stararray0[iup].t[jup-1])-mratio*(stararray0[iup].t[jup]-ratiolow*(stararray0[iup].t[jup]-stararray0[iup].t[jup-1])-stararray0[ilow].t[jlow-1]) << endl;
      else if (difflow>diffup) cerr << " Warning: newtrack.t[" << j << "]=" << newtrack.t[j] << "yr=" << stararray0[iup].t[jup-1]-mratio*(stararray0[iup].t[jup-1]-stararray0[ilow].t[jlow]-ratioup*(stararray0[ilow].t[jlow]-stararray0[ilow].t[jlow-1])) << endl;
      else cerr << "#Error: unknown calculation of age" << endl;*/
      debug = true;
      if ((jlow<nlow)||(jup<nup)){
        cerr << "#Error: end track time: jlow=" << jlow << " nlow=" << nlow << " jup=" << jup << " nup=" << nup << endl;
        screen = true;
      }
      jlow = stararray0[ilow].n;
      jup = stararray0[iup].n;
    }else if(newtrack.t[0]<newtrack.t[j]) newtrack.t[0] = newtrack.t[j];	//adjust maximal age with the accuracy
    if (debug && ((jlow>nlow)||(jup>nup))) cerr << "jlow=" << jlow << " nlow=" << nlow << " jup=" << jup << " nup=" << nup << endl;
    if ((newtrack.m[j]<0)||(newtrack.t[j]<0)) cerr << "#Error: newtrack of m_ini=" << newtrack.m[1] << "Msun t_max=" << newtrack.t[0] << "yr: j=" << j << " jlow=" << jlow << " jup=" << jup << " m=" << newtrack.m[j] << "Msun t=" << newtrack.t[j] << "yr r=" << newtrack.r[j] << "Rsun core mass(cm)=" << newtrack.cm[j] << "Msun" << endl;
  }
  if (newtrack.TAMS>=newtrack.n) newtrack.TAMS = newtrack.n-1;	//correct TAMS index
}

void changetrack(t_star& star, int n, t_HRD newtrack){
  int i,j;	//index variables: i: star; j: newtrack
  double mratio=-1.0;	//ratio of current values in newtrack
  double ageshift=0.0;	//age shift between the old and the new track
  double maxr=0.0;	//maximal radius in future
  double EbindperG=0.0;	//envelope binding energy for interpolation
  
  i = star.last;	//get last track position
  if (i==star.track.n-1) j = i-1;	//get previous track position if last track position is at the end
  else j = i+1;	//get next track position
  //check that current values are at the track position
  if (((fmin(star.track.m[i],star.track.m[j])>star.m[n])||(star.m[n]>fmax(star.track.m[i],star.track.m[j])))&&((fmin(star.track.m[i],star.track.m[j])>star.m[n-1])||(star.m[n-1]>fmax(star.track.m[i],star.track.m[j])))){
    cerr << "#Error: wrong track position: star.m[" << n << "]=" << star.m[n] << "Msun star.m[" << n-1 << "]=" << star.m[n-1] << "Msun star.track.m[" << i << "]=" << star.track.m[i] << "Msun star.track.m[" << j << "]=" << star.track.m[j] << "Msun" << endl;
    screen = true;
  }
  if (((fmin(star.track.t[i],star.track.t[j])>star.t[n])||(star.t[n]>fmax(star.track.t[i],star.track.t[j])))&&((fmin(star.track.t[i],star.track.t[j])>star.t[n-1])||(star.t[n-1]>fmax(star.track.t[i],star.track.t[j])))){
    cerr << "#Error: wrong track position: star.t[" << n << "]=" << star.t[n] << "yr star.t[" << n-1 << "]=" << star.t[n-1] << "yr star.track.t[" << i << "]=" << star.track.t[i] << "yr star.track.t[" << j << "]=" << star.track.t[j] << "yr" << endl;
    screen = true;
  }
  if (((fmin(star.track.r[i],star.track.r[j])>star.r[n])||(star.r[n]>fmax(star.track.r[i],star.track.r[j])))&&((fmin(star.track.r[i],star.track.r[j])>star.r[n-1])||(star.r[n-1]>fmax(star.track.r[i],star.track.r[j])))){
    cerr << "#Error: wrong track position: star.r[" << n << "]=" << star.r[n] << "Rsun star.r[" << n-1 << "]=" << star.r[n-1] << "Rsun star.track.r[" << i << "]=" << star.track.r[i] << "Rsun star.track.r[" << j << "]=" << star.track.r[j] << "Rsun" << endl;
    screen = true;
  }
  if (((fmin(star.track.cm[i],star.track.cm[j])>star.cm[n])||(star.cm[n]>fmax(star.track.cm[i],star.track.cm[j])))&&((fmin(star.track.cm[i],star.track.cm[j])>star.cm[n-1])||(star.cm[n-1]>fmax(star.track.cm[i],star.track.cm[j])))){
    cerr << "#Error: wrong track position: star.cm[" << n << "]=" << star.cm[n] << "Msun star.cm[" << n-1 << "]=" << star.cm[n-1] << "Msun star.track.cm[" << i << "]=" << star.track.cm[i] << "Msun star.track.cm[" << j << "]=" << star.track.cm[j] << "Msun" << endl;
    if (debug) cerr << "# star.cm[" << n << "]-star.cm[" << n-1 << "]=" << star.cm[n]-star.cm[n-1] << "Msun star.cm[" << n << "]-star.track.cm[" << i << "]=" << star.cm[n]-star.track.cm[i] << "Msun star.cm[" << n << "]-star.track.cm[" << j << "]=" << star.cm[n]-star.track.cm[j] << "Msun star.cm[" << n-1 << "]-star.track.cm[" << i << "]=" << star.cm[n-1]-star.track.cm[i] << "Msun star.cm[" << n-1 << "]-star.track.cm[" << j << "]=" << star.cm[n-1]-star.track.cm[j] << "Msun star.t[" << n << "]-star.t[" << n-1 << "]=" << star.t[n]-star.t[n-1] << "yr star.t[" << n-1 << "]-star.t[" << n-2 << "]=" << star.t[n-1]-star.t[n-2] << "yr" << endl;
    screen = true;
  }
  for (j=2;j<newtrack.n;j++){
    if ((fmin(newtrack.m[j-1],newtrack.m[j])<=star.m[n]*(1.0+0.1*accuracy))&&(fmax(newtrack.m[j-1],newtrack.m[j])*(1.0+0.1*accuracy)>=star.m[n])){
      if ((fmin(newtrack.cm[j-1],newtrack.cm[j])<=star.cm[n]*(1.0+0.1*accuracy))&&(fmax(newtrack.cm[j-1],newtrack.cm[j])*(1.0+0.1*accuracy)>=star.cm[n])){
        //calculate ratio in newtrack
        if (fabs(1.0-newtrack.m[j]/newtrack.m[j-1])>accuracy) mratio = ratio(star.m[n],newtrack.m[j],newtrack.m[j-1]);	//use mass for the ratio if mass change is large enough
        else mratio = ratio(star.cm[n],newtrack.cm[j],newtrack.cm[j-1]);	//otherwise use the coremass
        break;
      }
    }
  }
  if (j>=newtrack.n) j = newtrack.n-1;	//use last track point
  if ((star.cm[n]<newtrack.cm[1])&&(mratio==-1.0)){	//no significant core
    j = 2;
    mratio = ratio(star.m[n],newtrack.m[j],newtrack.m[j-1]);
  }else if (fabs(mratio-ratio(star.cm[n],newtrack.cm[j],newtrack.cm[j-1]))>1.732051*accuracy*(1.0/fabs(newtrack.cm[j]-newtrack.cm[j-1])+1.0/fabs(newtrack.m[j]-newtrack.m[j-1]))){
    cerr << "#Error: mratio-cmratio=" << fabs(mratio-ratio(star.cm[n],newtrack.cm[j],newtrack.cm[j-1])) << ">" << 1.732051*accuracy*(1.0/fabs(newtrack.cm[j]-newtrack.cm[j-1])+1.0/fabs(newtrack.m[j]-newtrack.m[j-1])) << " mratio=" << mratio << "=" << newtrack.m[j]-star.m[n] << "/" << newtrack.m[j]-newtrack.m[j-1] << "=(" << newtrack.m[j] << "-"  << star.m[n] << ")/(" << newtrack.m[j] << "-" << newtrack.m[j-1] << ") cmratio=" << ratio(star.cm[n],newtrack.cm[j],newtrack.cm[j-1]) << "=" << newtrack.cm[j]-star.cm[n] << "/" << newtrack.cm[j]-newtrack.cm[j-1] << "=(" << newtrack.cm[j] << "-"  << star.cm[n] << ")/(" << newtrack.cm[j] << "-" << newtrack.cm[j-1] << ")" << endl;	//check if mass and core mass give the same position
//    cerr << 1.0-fmin(newtrack.m[j-1],newtrack.m[j])/star.m[n] << " " << 1.0-fmax(newtrack.m[j-1],newtrack.m[j])/star.m[n] << endl;
    screen = true;
  }
  if (debug) cerr << "mratio=" << mratio << endl;
  if ((mratio+1.732051*accuracy/fabs(newtrack.m[j]-newtrack.m[j-1])<0.0)||(mratio>1.0+1.732051*accuracy/fabs(newtrack.m[j]-newtrack.m[j-1]))){
    cerr << "#Error: mratio=" << mratio << " not in [" << -1.732051*accuracy/fabs(newtrack.m[j]-newtrack.m[j-1]) << "," << 1.0+1.732051*accuracy/fabs(newtrack.m[j]-newtrack.m[j-1]) << "] star.cm[" << n << "]=" << star.cm[n] << "Msun star.m[" << n << "]=" << star.m[n] << "Msun newtrack.m[" << j << "]=" << newtrack.m[j] << "Msun newtrack.m[" << j-1 << "]=" << newtrack.m[j-1] << "Msun" << endl;	//check if current mass between two track points
    fprintf(stderr,"# mratio=%.15f=%.15e/%.15e=(%.15e-%.15e)/(%.15e-%.15e)\n",mratio,newtrack.m[j]-star.m[n],newtrack.m[j]-newtrack.m[j-1],newtrack.m[j],star.m[n],newtrack.m[j],newtrack.m[j-1]);
//    cerr << "# mratio=" << mratio << "=" << newtrack.m[j]-star.m[n] << "/" << newtrack.m[j]-newtrack.m[j-1] << " 1-mratio=" << 1.0-mratio << endl;
    screen = true;
  }
  if ((fabs(star.track.m[i]-star.m[n])>accuracy)||(fabs(star.track.t[i]-star.t[n])>accuracy)||(fabs(star.track.r[i]-star.r[n])>accuracy)||(fabs(star.track.cm[i]-star.cm[n])>accuracy)||(fabs(star.track.llum[i]-star.llum[n])>accuracy)||(fabs(star.track.lteff[i]-star.lteff[n])>accuracy)||(fabs(star.track.lambda[i]-star.lambda[n])>accuracy)||(fabs(star.track.cc[i]-star.cc[n])>accuracy)||(fabs(star.track.cf[i]-star.cf[n])>accuracy)||(fabs(star.track.omega[i]-star.omega[n])>accuracy)) star.last++;	//increase trackmarker
  //get length for the arrays (allready used track + newtrack - not needed part of newtrack)
  star.track.n = star.last+newtrack.n-(j-1);
  //reserve memory
  star.track.m = (double *)realloc(star.track.m,star.track.n*sizeof(double));
  star.track.t = (double *)realloc(star.track.t,star.track.n*sizeof(double));
  star.track.r = (double *)realloc(star.track.r,star.track.n*sizeof(double));
  star.track.cm = (double *)realloc(star.track.cm,star.track.n*sizeof(double));
  star.track.llum = (double *)realloc(star.track.llum,star.track.n*sizeof(double));
  star.track.lteff = (double *)realloc(star.track.lteff,star.track.n*sizeof(double));
  star.track.lambda = (double *)realloc(star.track.lambda,star.track.n*sizeof(double));
  star.track.cc = (double *)realloc(star.track.cc,star.track.n*sizeof(double));
  star.track.cf = (double *)realloc(star.track.cf,star.track.n*sizeof(double));
  star.track.omega = (double *)realloc(star.track.omega,star.track.n*sizeof(double));
  //check if memory allocation fails
  if (star.track.m==NULL) cerr << "#Error: memory reallocation failed: star.track.m" << endl;
  if (star.track.t==NULL) cerr << "#Error: memory reallocation failed: star.track.t" << endl;
  if (star.track.r==NULL) cerr << "#Error: memory reallocation failed: star.track.r" << endl;
  if (star.track.cm==NULL) cerr << "#Error: memory reallocation failed: star.track.cm" << endl;
  if (star.track.llum==NULL) cerr << "#Error: memory reallocation failed: star.track.llum" << endl;
  if (star.track.lteff==NULL) cerr << "#Error: memory reallocation failed: star.track.lteff" << endl;
  if (star.track.lambda==NULL) cerr << "#Error: memory reallocation failed: star.track.lambda" << endl;
  if (star.track.cc==NULL) cerr << "#Error: memory reallocation failed: star.track.cc" << endl;
  if (star.track.cf==NULL) cerr << "#Error: memory reallocation failed: star.track.cf" << endl;
  if (star.track.omega==NULL) cerr << "#Error: memory reallocation failed: star.track.omega" << endl;
  //write current values as track point
  star.track.m[star.last] = star.m[n];
  star.track.t[star.last] = star.t[n];
  star.track.r[star.last] = star.r[n];
  star.track.r[star.last] = newtrack.r[j]-mratio*(newtrack.r[j]-newtrack.r[j-1]);
  if (abs(star.stage[n])<3){	// non compact star
    if (star.track.r[star.last]<WDradius(star.track.m[star.last])) star.track.r[star.last] = WDradius(star.track.m[star.last]);	// maximum radius is given for full electron degeneracy (WD radius)
  }
  star.track.cm[star.last] = star.cm[n];
  star.track.llum[star.last] = star.llum[n];
  star.track.llum[star.last] = newtrack.llum[j]-mratio*(newtrack.llum[j]-newtrack.llum[j-1]);
  star.track.lteff[star.last] = star.lteff[n];
  star.track.lteff[star.last] = newtrack.lteff[j]-mratio*(newtrack.lteff[j]-newtrack.lteff[j-1]);
  star.track.lambda[star.last] = star.lambda[n];
//  star.track.lambda[star.last] = newtrack.lambda[j]-mratio*(newtrack.lambda[j]-newtrack.lambda[j-1]);
  EbindperG = newtrack.m[j]*(newtrack.m[j]-newtrack.cm[j])/(newtrack.lambda[j]*newtrack.r[j])
              -mratio*(newtrack.m[j]*(newtrack.m[j]-newtrack.cm[j])/(newtrack.lambda[j]*newtrack.r[j])
                      -newtrack.m[j-1]*(newtrack.m[j-1]-newtrack.cm[j-1])/(newtrack.lambda[j-1]*newtrack.r[j-1]));	//interpolate current binding energy of the envelope
  if ((star.track.m[star.last]-star.track.cm[star.last]==0.0)&&(EbindperG==0.0)) star.track.lambda[star.last] = ignorev;	//no interpolation
  else star.track.lambda[star.last] = star.track.m[star.last]*(star.track.m[star.last]-star.track.cm[star.last])/(EbindperG*star.track.r[star.last]);	//interpolate lambda	// E_bind=G*M*M_env / lambda*R
  if(isnan(star.track.lambda[star.last])){ screen = true; cerr << endl << "#star.track.lambda[star.last]=" << star.track.lambda[star.last] << "=" << star.track.m[star.last] << "*(" << star.track.m[star.last] << "-" << star.track.cm[star.last] << ")/(" << EbindperG << "*" << star.track.r[star.last] << ") j=" << j;}
  star.track.cc[star.last] = star.cc[n];
  star.track.cc[star.last] = newtrack.cc[j]-mratio*(newtrack.cc[j]-newtrack.cc[j-1]);
  star.track.cf[star.last] = star.cf[n];
  star.track.cf[star.last] = newtrack.cf[j]-mratio*(newtrack.cf[j]-newtrack.cf[j-1]);
  star.track.omega[star.last] = star.omega[n];
  star.track.omega[star.last] = newtrack.omega[j]-mratio*(newtrack.omega[j]-newtrack.omega[j-1]);
  if (mratio>1){	//check for unphsical extrapolation
    if (star.track.m[star.last]<=0.0){
      if (screen) cerr << "#star.track.m[" << star.last << "]=" << star.track.m[star.last] << "Msun" << endl;
      star.track.m[star.last] = 1.0e-99;
    }
    if (star.track.r[star.last]<=0.0){
      if (screen) cerr << "#star.track.r[" << star.last << "]=" << star.track.r[star.last] << "Msun" << endl;
      star.track.r[star.last] = 1.0e-99;
    }
    if (star.track.cm[star.last]<=0.0){
      if (screen) cerr << "#star.track.cm[" << star.last << "]=" << star.track.cm[star.last] << "Msun" << endl;
      star.track.cm[star.last] = 1.0e-99;
    }
    if (star.track.cc[star.last]<=0.0){
      if (screen) cerr << "#star.track.cc[" << star.last << "]=" << star.track.cc[star.last] << "Msun" << endl;
      star.track.cc[star.last] = 1.0e-99;
    }
  }
  star.inimass = newtrack.inimass;	//get new initial mass
  if (star.stage[n]==0){	//star is on MS
    star.track.TAMS = star.last+1+newtrack.TAMS-j;	//get new index of TAMS
    if (debug) cerr << "newtrack.TAMS=" << newtrack.TAMS << " star.track.TAMS=" << star.track.TAMS << endl;
  }
  //calculate age shift between the tracks
  ageshift = star.t[n]-(newtrack.t[j]-mratio*(newtrack.t[j]-newtrack.t[j-1]));
  if (debug) cerr << "star.m[" << n << "]=" << star.m[n] << "Msun=" << newtrack.m[j]-mratio*(newtrack.m[j]-newtrack.m[j-1]) << "Msun star.t[" << n << "]=" << star.t[n] << "yr newtrack age=" << newtrack.t[j]-mratio*(newtrack.t[j]-newtrack.t[j-1]) << "yr ageshift=" << ageshift << "yr" << endl;
  if (debug){
    cerr << "index\tm in Msun\tt in yr\tr in Rsun\tcm in Msun\tllum\tlteff\tlambda\tcc in Msun\tcf\tomega in 1/yr" << endl;	//write head of primary track-table
    cerr << star.last << "\t" << star.track.m[star.last] << "\t" << star.track.t[star.last] << "\t" << star.track.r[star.last] << "\t" << star.track.cm[star.last] << "\t" << star.track.llum[star.last] << "\t" << star.track.lteff[star.last] << "\t" << star.track.lambda[star.last] << "\t" << star.track.cc[star.last] << "\t" << star.track.cf[star.last] << "\t" << star.track.omega[star.last] << endl;	//write new track-table of the star
  }
  //prepare maximal values
  star.track.m[0] = star.track.m[1];
  star.track.r[0] = star.track.r[1];
  star.rmax = 1;
  star.track.cm[0] = newtrack.cm[newtrack.n-1];
  star.track.llum[0] = newtrack.llum[newtrack.n-1];
  star.track.lteff[0] = star.track.lteff[1];
  star.track.lambda[0] = star.track.lambda[1];
  star.track.cc[0] = newtrack.cc[newtrack.n-1];
  star.track.cf[0] = star.track.cf[1];
  star.track.omega[0] = star.track.omega[1];
  maxr = star.track.r[star.last];
  for (i=2;i<=star.last;i++){
    //determine maximal values
    if (star.track.m[0]<star.track.m[i]) star.track.m[0] = star.track.m[i];
    if (star.track.r[0]<star.track.r[i]){
      star.track.r[0] = star.track.r[i];
      star.rmax = i;
    }
    if (star.track.cm[0]<star.track.cm[i]) star.track.cm[0] = star.track.cm[i];
    if (star.track.llum[0]<star.track.llum[i]) star.track.llum[0] = star.track.llum[i];
    if (star.track.lteff[0]<star.track.lteff[i]) star.track.lteff[0] = star.track.lteff[i];
    if (star.track.lambda[0]<star.track.lambda[i]) star.track.lambda[0] = star.track.lambda[i];
    if (star.track.cc[0]<star.track.cc[i]) star.track.cc[0] = star.track.cc[i];
    if (star.track.cf[0]<star.track.cf[i]) star.track.cf[0] = star.track.cf[i];
    if (star.track.omega[0]<star.track.omega[i]) star.track.omega[0] = star.track.omega[i];
  }
  //write new track points
  for (i=star.last+1;i<star.track.n;i++){
    if (j>=newtrack.n){
      cerr << "#Error: too many new track points: i=" << i << " star.track.n=" << star.track.n << endl;
      break;
    }
    if (i==star.last+1) star.track.m[i] = fmin(newtrack.m[j],star.track.m[i-1]);	//write mass
    else star.track.m[i] = newtrack.m[j];	//write mass
    star.track.t[i] = fmax(ageshift+newtrack.t[j],star.track.t[i-1]+1.0);	//write age: minimal age is previous age + 1yr
    star.track.r[i] = newtrack.r[j];	//write radius
    if (abs(star.stage[n])<3){	// non compact star
      if (star.track.r[i]<WDradius(star.track.m[i])) star.track.r[i] = WDradius(star.track.m[i]);	// maximum radius is given for full electron degeneracy (WD radius)
    }
    star.track.cm[i] = newtrack.cm[j];	//write core mass
    star.track.llum[i] = newtrack.llum[j];	//write luminosity
    star.track.lteff[i] = newtrack.lteff[j];	//write effective temperature
    star.track.lambda[i] = newtrack.lambda[j];	//write lambda
    star.track.cc[i] = newtrack.cc[j];	//write carbon core mass
    star.track.cf[i] = newtrack.cf[j];	//write concentration factor
    star.track.omega[i] = star.track.omega[i-1];	//write angular velocity
    if (debug) cerr << i << "\t" << star.track.m[i] << "\t" << star.track.t[i] << "\t" << star.track.r[i] << "\t" << star.track.cm[i] << "\t" << star.track.llum[i] << "\t" << star.track.lteff[i] << "\t" << star.track.lambda[i] << "\t" << star.track.cc[i] << "\t" << star.track.cf[i] << "\t" << star.track.omega[i] << endl;	//write new track-table of the star
    j++;
    //determine maximal values
    if (star.track.m[0]<star.track.m[i]) star.track.m[0] = star.track.m[i];
    if (maxr<star.track.r[i]){
      maxr = star.track.r[i];
      if (star.track.r[0]<star.track.r[i]) star.track.r[0] = star.track.r[i];
      if ((i<=star.track.TAMS)||(star.last>=star.track.TAMS)||(abs(star.stage[n])==2)) star.rmax = i;
    }
    if (star.track.cm[0]<star.track.cm[i]) star.track.cm[0] = star.track.cm[i];
    if (star.track.llum[0]<star.track.llum[i]) star.track.llum[0] = star.track.llum[i];
    if (star.track.lteff[0]<star.track.lteff[i]) star.track.lteff[0] = star.track.lteff[i];
    if (star.track.lambda[0]<star.track.lambda[i]) star.track.lambda[0] = star.track.lambda[i];
    if (star.track.cc[0]<star.track.cc[i]) star.track.cc[0] = star.track.cc[i];
    if (star.track.cf[0]<star.track.cf[i]) star.track.cf[0] = star.track.cf[i];
    if (star.track.omega[0]<star.track.omega[i]) star.track.omega[0] = star.track.omega[i];
  }
  if (j<newtrack.n) cerr << "#Error: not all new track points copied: j=" << j << " newtrack.n=" << newtrack.n << endl;
  star.track.t[0] = star.track.t[star.track.n-1];	//write maximal time
}

long double getA(long double j11,long double j12,long double j21,long double j22,long double m11,long double m12,long double m21,long double m22,long double cm11,long double cm12,long double cm21,long double cm22){
  return (j22-j21)*(cm11*m12-cm12*m11) + (j12-j22)*(cm11*m21-cm21*m11) + (j22-j11)*(cm12*m21-cm21*m12)	//=(j22-j21)*(c11*m12-c12*m11)+(j12-j22)*(c11*m21-c21*m11)+(j22-j11)*(c12*m21-c21*m12)
       + (j12-j21)*(cm22*m11-cm11*m22) + (j21-j11)*(cm22*m12-cm12*m22) + (j12-j11)*(cm21*m22-cm22*m21);	//+(j12-j21)*(c22*m11-c11*m22)+(j21-j11)*(c22*m12-c12*m22)+(j12-j11)*(c21*m22-c22*m21)
}

long double getB(long double j21,long double j22,long double m,long double m11,long double m12,long double cm,long double cm11,long double cm12){
  return (j22-j21)*(cm*(m11-m12)+cm11*(m12-m)+cm12*(m-m11));	//=(j22-j21)*[c*(m11-m12)+c11*(m12-m)+c12*(m-m11)]
}

long double getC(long double j11,long double j12,long double m,long double m21,long double m22,long double cm,long double cm21,long double cm22){
  return (j12-j11)*(cm*(m21-m22)+cm21*(m22-m)+cm22*(m-m21));	//=(j12-j11)*[c*(m21-m22)+c21*(m22-m)+c22*(m-m21)]
}

long double getr1(long double m,long double m11,long double m12,long double m21,long double m22,long double j11,long double j12,long double j21,long double j22,long double r){
  return ((m-m11*(1.0-r))*(j22-j21)+r*(m22*(j21-j11)-m21*(j22-j11)))/((m22-m21)*r*(j12-j11)+(m12-m11)*(1.0-r)*(j22-j21));	//get ratio in upper track: {(m-m11*(1-r))*(j22-j21)+r*[m22*(j21-j11)-m21*(j22-j11)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
}

long double getr2(long double m,long double m11,long double m12,long double m21,long double m22,long double j11,long double j12,long double j21,long double j22,long double r){
  return ((m-m21*r)*(j12-j11)+(1.0-r)*(m12*(j11-j21)-m11*(j12-j21)))/((m22-m21)*r*(j12-j11)+(m12-m11)*(1.0-r)*(j22-j21));	//get ratio in lower track: {(m-m21*r)*(j12-j11)+(1-r)*[m12*(j11-j21)-m11*(j12-j21)]}/{(m22-m21)*r*(j12-j11)+(m12-m11)*(1-r)*(j22-j21)}
}

double getCOcore(double& McoreHe){	//calculate the CO core mass from the He and C core masses by taking the CO core of the last track point in the naked Helium star grid
  int i;	//index variable
  double mratio;	//mass ratio of initial He mass

  for (i=0;i<nhe1;i++){	//find mass index
    if(hestararray1[i].inimass>McoreHe){	//hestararray1 must be sorted from low mass to high mass (done when table is read in)
      break;
    }
  }
  if (i>=nhe1) i = nhe1-1;	//out of He star grid --> extrapolate
  mratio = ratio(McoreHe,hestararray1[i].inimass,hestararray1[i-1].inimass);	//get mass ratio of initial He mass
  McoreHe = hestararray1[i].m[hestararray1[i].n-1]-mratio*(hestararray1[i].m[hestararray1[i].n-1]-hestararray1[i-1].m[hestararray1[i-1].n-1]);	//get He mass of last track point
  return hestararray1[i].cm[hestararray1[i].n-1]-mratio*(hestararray1[i].cm[hestararray1[i].n-1]-hestararray1[i-1].cm[hestararray1[i-1].n-1]);	//get CO core mass of last track point
}

void getcoremasses(t_star& star, int n){	//calculate the CO core mass from the He and C core masses by taking the CO core of the last track point in the naked Helium star grid
//  debug=true;
  star.SN.remnantmass = star.m[n];
  star.SN.Hecoremass = star.cm[n];
  star.SN.COcoremass = star.cc[n];
  star.m[n] = star.SN.Hecoremass;
  star.cm[n] = star.SN.COcoremass;
  updatestarstrack(star, n, hestararray1, nhe1, true);	//get final masses of a naked helium star with the current He mass and CO core mass
}

void undophase(t_star& star, int n, bool keept){	//undo the phase n for a star
  star.m[n] = star.m[n-1];
  if (!keept) star.t[n] = star.t[n-1];
  star.r[n] = star.r[n-1];
  star.cm[n] = star.cm[n-1];
  star.llum[n] = star.llum[n-1];
  star.lteff[n] = star.lteff[n-1];
  star.lambda[n] = star.lambda[n-1];
  star.cc[n] = star.cc[n-1];
  star.cf[n] = star.cf[n-1];
  star.omega[n] = star.omega[n-1];
  star.stage[n] = star.stage[n-1];
  while ((star.track.t[star.last]>star.t[n])&&(star.last>0)){	//reset last track point
    star.last--;
  }
  if (debug) cerr << "reset last to " << star.last << endl;
}

void setlasttrack(t_star& star, int n){	//set last track to given phase
  if (debug) cerr << "set last track at " << star.last << endl;
  star.track.m[star.last] = star.m[n];
  star.track.t[star.last] = star.t[n];
  star.track.r[star.last] = star.r[n];
  star.track.cm[star.last] = star.cm[n];
  star.track.llum[star.last] = star.llum[n];
  star.track.lteff[star.last] = star.lteff[n];
  star.track.lambda[star.last] = star.lambda[n];
  star.track.cc[star.last] = star.cc[n];
  star.track.cf[star.last] = star.cf[n];
  star.track.omega[star.last] = star.omega[n];
}
