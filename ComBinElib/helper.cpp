//#include <iostream>	//in-/output to/from screen
//#include <stdlib.h>	//provides standard functions (conversion, memory allocation, ...)
//#include <math.h>	//provides the mathematical functions
#include <stdio.h>
#include "ComBinElib.h" // local header file

//using namespace std;

double square(double value){	//calculates the square of a value
  return value*value;
}

double cubic(double value){	//calculates the third power of a value
  return value*value*value;
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
