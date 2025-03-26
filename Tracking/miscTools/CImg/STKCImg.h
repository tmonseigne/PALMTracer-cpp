




#include "CImg.h"
#include "../../generalParam/NoError.h"
#include "../../miscTools\detectLine/lineDetect.h"
#include "../../array_library_static/allDataLib.h"
#include <stdlib.h>


int display(USrtArray & img, double &gamma, char* fileName=NULL, unsigned short nClick=0, unsigned short *getClick=NULL);
int display(USrtArray & input, double &gamma, char* windowName, contour ***contours, long num_cont, long minNbPoints);