//      galactic.cpp
//      This subroutine calculates the movement of the system

//#include <iostream>	//in-/output to/from screen
//#include <math.h>	//provides the mathematical functions
#include "ComBinElib.h" // local header file
#include <limits.h>	//provides constants with limits

//using namespace std;

/*global variables*/
integratorfunction integrator;	//contains the used integrator
potentialfunction potential;	//contains the used potential

void movesystem(t_system& system, int m, int n){	//calculates movement of the system
  double y[6];
  double dt=system.prim.t[n]-system.prim.t[m];
  double t=system.prim.t[m];
  int i=0;
  t_star* star;	//star
  t_star* companion;	//companion

  y[0] = system.x[m];	//initial in pc x-position
  y[1] = system.y[m];	//initial in pc y-position
  y[2] = system.z[m];	//initial in pc z-position
  y[3] = system.vx[m]*kms/pc;	//initial x-velocity in pc/yr	//convert velocity from km/s to pc/yr
  y[4] = system.vy[m]*kms/pc;	//initial y-velocity in pc/yr	//convert velocity from km/s to pc/yr
  y[5] = system.vz[m]*kms/pc;	//initial z-velocity in pc/yr	//convert velocity from km/s to pc/yr

  if(snapshots){
    for (i=1;i<15;i++){
      if ((t<i*1.0e+9)&&(system.prim.t[n]>=i*1.0e+9)){
        dt = i*1.0e+9-t;
        break;
      }
    }
  }

  while(dt>0){
    integrator(y,dt);	//integrate over specified time

    if(isnan(y[0])||isnan(y[1])||isnan(y[2])||isnan(y[3])||isnan(y[4])||isnan(y[5])){
      debug = true;
      y[0] = system.x[m];	//initial in pc x-position
      y[1] = system.y[m];	//initial in pc y-position
      y[2] = system.z[m];	//initial in pc z-position
      y[3] = system.vx[m]*kms/pc;	//initial x-velocity in pc/yr	//convert velocity from km/s to pc/yr
      y[4] = system.vy[m]*kms/pc;	//initial y-velocity in pc/yr	//convert velocity from km/s to pc/yr
      y[5] = system.vz[m]*kms/pc;	//initial z-velocity in pc/yr	//convert velocity from km/s to pc/yr
      integrator(y,dt);	//integrate over specified time
    }

    if(snapshots){
      t += dt;
      if(i<14){
        if (system.prim.m[system.n-1]<0.3){
          star = &system.prim;
          companion = &system.sec;
        }else{
          star = &system.sec;
          companion = &system.prim;
        }
        switch(i){
          case 1: fprintf(snapshot1,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 2: fprintf(snapshot2,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 3: fprintf(snapshot3,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 4: fprintf(snapshot4,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 5: fprintf(snapshot5,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 6: fprintf(snapshot6,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 7: fprintf(snapshot7,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 8: fprintf(snapshot8,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 9: fprintf(snapshot9,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 10: fprintf(snapshot10,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 11: fprintf(snapshot11,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 12: fprintf(snapshot12,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
          case 13: fprintf(snapshot13,"%10.3e\t%10.3e %10.3e %10.3e\t%10.3e %10.3e %10.3e\t%10.6f %2d \t%10.6f %2d \n", t, y[0], y[1], y[2], y[3]*pc/kms, y[4]*pc/kms, y[5]*pc/kms, star[0].m[m], star[0].stage[m], companion[0].m[m], companion[0].stage[m]); break;
        }
        fflush(snapshot1);	//update file
        fflush(snapshot2);	//update file
        fflush(snapshot3);	//update file
        fflush(snapshot4);	//update file
        fflush(snapshot5);	//update file
        fflush(snapshot6);	//update file
        fflush(snapshot7);	//update file
        fflush(snapshot8);	//update file
        fflush(snapshot9);	//update file
        fflush(snapshot10);	//update file
        fflush(snapshot11);	//update file
        fflush(snapshot12);	//update file
        fflush(snapshot13);	//update file
      }
      dt = system.prim.t[n]-t;
      for (i=1;i<15;i++){
        if ((t<i*1.0e+9)&&(system.prim.t[n]>=i*1.0e+9)){
          dt = i*1.0e+9-t;
          break;
        }
      }
    }else{
      dt = 0;
    }
  }

  system.x[n]  = y[0];	//save move in x-direction
  system.y[n]  = y[1];	//save move in y-direction
  system.z[n]  = y[2];	//save move in z-direction
  system.vx[n] = y[3]*pc/kms;	//save new x-velocity	//convert velocity back from pc/yr to km/s
  system.vy[n] = y[4]*pc/kms;	//save new y-velocity	//convert velocity back from pc/yr to km/s
  system.vz[n] = y[5]*pc/kms;	//save new z-velocity	//convert velocity back from pc/yr to km/s
}

void rk4(double* y, double time){	//Runge-Kutta 4 integrator: vector y={x,y,z,vx,vy,vz}
  double yini[6]={y[0],y[1],y[2],y[3],y[4],y[5]};	//save initial values
  double ylast[6]={y[0]+y[3]*time,y[1]+y[4]*time,y[2]+y[5]*time,y[3],y[4],y[5]};	//saved results: constant velocity(first guess)
  double yinter[6];	//auxiliary variables: intermediate y
  double k[4][6];	//auxiliary variables: k-factor in rk
  double deltay[6];	//auxiliary variables: intermediate y
  double dt=time;	//time step
  double dt2=0.5*dt;	//half time step
  const double lgmaxi=log10(INT_MAX)+1.0;	//log10 value of maximal integer and add one for step size control
  double err1=2.0,nchange;
  static int i,j;	//loop variables
  int n=1;	//number of steps
  int nlast=n;	//saved number of steps

  if (debug &&(dt!=0.0)) cerr << endl << "time to integrate=" << time << "yr {x,y,z}={" << y[0] << "," << y[1] << "," << y[2] << "}pc {vx,vy,vz}={" << y[3] << "," << y[4] << "," << y[5] << "}pc/yr";
  while (err1>1.0){
    dt  = time/(double)n;	//new time step size
    dt2 = 0.5*dt;	//new half time step size
    for (j=0;j<6;j++) y[j] = yini[j];
    for (i=0;i<n;i++){	//do integration
      if (debug &&(n>16)&&(i%(n/16)==0)) cerr << endl << "n=" << n << ": t=" << i*dt << "yr {x,y,z}={" << y[0] << "," << y[1] << "," << y[2] << "}pc {vx,vy,vz}={" << y[3] << "," << y[4] << "," << y[5] << "}pc/yr";
      for (j=0;j<3;j++) k[0][j] = y[j+3];	//get new k for space
      potential(k[0][3], k[0][4], k[0][5], y[0]*0.001, y[1]*0.001, y[2]*0.001);	//get new k for velocity
      for (j=0;j<6;j++) yinter[j] = y[j] + dt2*k[0][j];	//first intermediate step
      for (j=0;j<3;j++) k[1][j] = yinter[j+3];	//get new k for space
      potential(k[1][3], k[1][4], k[1][5], yinter[0]*0.001, yinter[1]*0.001, yinter[2]*0.001);	//get new k for velocity
      for (j=0;j<6;j++) yinter[j] = y[j] + dt2*k[1][j];	//second intermediate step
      for (j=0;j<3;j++) k[2][j] = yinter[j+3];	//get new k for space
      potential(k[2][3], k[2][4], k[2][5], yinter[0]*0.001, yinter[1]*0.001, yinter[2]*0.001);	//get new k for velocity
      for (j=0;j<6;j++) yinter[j] = y[j] + dt*k[2][j];	//third intermediate step
      for (j=0;j<3;j++) k[3][j] = yinter[j+3];	//get new k for space
      potential(k[3][3], k[3][4], k[3][5], yinter[0]*0.001, yinter[1]*0.001, yinter[2]*0.001);	//get new k for velocity
      for (j=0;j<6;j++) y[j] += dt/3.0*(0.5*(k[0][j]+k[3][j])+k[1][j]+k[2][j]);	//full step
    }
    nchange = (double)n/(double)nlast;	//get ratio of last an current number of time steps
    for (j=0;j<6;j++){
      deltay[j] = (y[j]-ylast[j])/(1.0-nchange*nchange*nchange);	//get truncation error estimate
      if(isnan(deltay[j])) deltay[j] = 0.01;	//=> err1>10/sqrt(3)
    }
    if (debug) cout << endl << "deltay={" << deltay[0] << "," << deltay[1] << "," << deltay[2] << "," << deltay[3] << "," << deltay[4] << "," << deltay[5] << "}";
    err1 = sqrt(1.0/3.0)*sqrt(deltay[0]*deltay[0]+deltay[1]*deltay[1]+deltay[2]*deltay[2])*1.0E+3;	//N=3; scale_x=scale_y=scale_z=1.0E-3 	//space coordinates accurate to 10^-3pc and no constrains on the velocity
    if (debug) cerr << endl << "n=" << n << " dt=" << dt << "yr {x,y,z}={" << y[0] << "," << y[1] << "," << y[2] << "}pc {vx,vy,vz}={" << y[3] << "," << y[4] << "," << y[5] << "}pc/yr";
    if ((dt==dt2)||(err1<1.0)) break;	//stop integration: timestep to small or accuracy reached
    nlast = n;	//save last number of time steps
    n *= fmin(fmax(1.001*pow(err1,0.25),1.0/n+1.0),lgmaxi-log10(n));	//get new number of time steps: use err1, minimal increase of 1 and maximal increase of a factor depending on the magnitude of the time step number
    if (n<0){
      cerr << "#Error: integer overflow n=" << n << endl;
      screen = true;
    }
    for (j=0;j<6;j++) ylast[j] = y[j];	//save last results
  }
}

void MW_potential(double& ax, double& ay, double& az, double x, double y, double z){	//calculates the acceleration in the galactic potential where x,y,z is in kpc and ax,ay,az in pc/yr^2
//Ref: Allen, C., Santillan, A., 1991, RevMexAA, 22, 255
  const double conversion=0.1*kms/pc*kms/pc;	//from 100(km/s)^2 to kpc*pc/yr^2
  const double M1= 606.0;	//in 2.32e+7Msun
  const double b1= 0.3873;	//in kpc
  const double M2= 3690.0;	//in 2.32e+7Msun
  const double a2= 5.3178;	//in kpc
  const double b2= 0.2500;	//in kpc
  const double M3= 4615.0;	//in 2.32e+7Msun
  const double a3= 12.0;	//in kpc
  static double rr = x*x+y*y;	//in kpc^2
  static double RR = rr+z*z;	//in kpc^2
  static double R  = sqrt(RR);	//in kpc
  static double Ra3= R/a3;	//no unit
  static double Ra = pow(Ra3,1.02);	//no unit
  static double Zb = sqrt(z*z+b2*b2);	//in kpc
  static double Zba= a2+Zb;	//in kpc
  static double C1 = M1*pow(RR+b1*b1,-1.5);	//in 2.32e+7Msun/kpc^3=100(km/s)^2/kpc^2
  static double C2 = M2*pow(rr+Zba*Zba,-1.5);	//in 2.32e+7Msun/kpc^3=100(km/s)^2/kpc^2
  static double C3 = M3*Ra/(RR*a3*(1.0+Ra));	//in 2.32e+7Msun/kpc^3=100(km/s)^2/kpc^2

  ax = -conversion*(C1+C2+C3)*x;
  ay = -conversion*(C1+C2+C3)*y;
  az = -conversion*(C1+C2*Zba/Zb+C3)*z;
}
