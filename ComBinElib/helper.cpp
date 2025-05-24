//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
#include <stdio.h>
#include "ComBinElib.h" // local header file

//using namespace std;

double square(double value){	//calculates the square of a value
  if (value==ignorev) return ignorev;
  else return value*value;
}

double cubic(double value){	//calculates the third power of a value
  if (value==ignorev) return ignorev;
  else return value*value*value;
}

double pow10(double value){	//calculates 10 to the power of a value
  if (value==ignorev) return ignorev;
  else return pow(10.0,value);
}

double ignorevToNan(double value){	//converts the ignore value to a quiet nan
  if (value==ignorev) return nan("9");
  else return value;
}

double ratio(double value, double high, double low){	//calculates the ratio and consider zero interval width
  if (high==low){	//if interval width is zero
    if (value<low){	//if value is below the interval
      return 1.0;	//tread as value=low
    }else if (value>high){	//if value is above the interval
      return 0.0;	//tread as value=high
    }else{		//if value is in the interval
      return 0.5;	//tread as value=0.5*(high+low)
    }
  }else{
    return (high-value)/(high-low);	//calculate ratio
  }
}

double interpolate(double ratio, double high, double low){	//calculates the interpolated value, don't allow extrapolation
  if (ratio<=0.0){	//outside the interpolation
    cerr << "#Extrapolation: ratio=" << ratio;
    return high;
  }else if (ratio>=1.0){	//outside the interpolation
    cerr << "#Extrapolation: ratio=" << ratio;
    return low;
  }else{
    return high-ratio*(high-low);
  }
}

double getEnvelopeBindingEnergy(t_system& system, int j){	//calculates the envelope's binding energy
  if (system.prim.r[j]/system.rp[j] >= system.sec.r[j]/system.rs[j]){
    if (system.prim.lambda[j]==ignorev) return ignorev;
    else return G*system.prim.m[j]*(system.prim.m[j]-system.prim.cm[j])/(system.prim.lambda[j]*system.prim.r[j])/cgsEnergy;
  }else{
    if (system.sec.lambda[j]==ignorev) return ignorev;
    else return G*system.sec.m[j]*(system.sec.m[j]-system.sec.cm[j])/(system.sec.lambda[j]*system.sec.r[j])/cgsEnergy;
  }
}
