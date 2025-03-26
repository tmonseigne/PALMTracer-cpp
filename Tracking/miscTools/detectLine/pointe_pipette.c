#include "lines.h"
 pointe_pipette(contour ***contours, long num_cont)
{
  contour *cont;
  long    i;
  int *colmin,*rowmin,j,k,*imin,index;
  double dmin,dmin0;


  imin = (int*)calloc(10,sizeof(int));
  colmin = (int*)calloc(10,sizeof(int));
  rowmin = (int*)calloc(10,sizeof(int));
  for (i=0,k=0; i<num_cont; i++) {
    colmin[k]=10000;    
    cont = (*contours)[i];
    /* Suppression des contours courts*/
    if(cont->num>100){
      /* Determination des extremites */
      for(j=0;j<cont->num;j++)
	if((int)cont->col[j]<colmin[k]){
	  colmin[k] = (int)cont->col[j];
	  rowmin[k] =  (int)cont->row[j];
	}
      printf("colmin=%d\t rowmin=%d\n",colmin[k],rowmin[k] );
      imin[k++]=i;
    }
  }

  if(k==2){
    /* Recherche du 3ene segment */
    dmin=10000;
    for (i=0; i<num_cont; i++) {
      cont = (*contours)[i];
      if((int)i!=imin[0] && (int)i!=imin[1]){
	dmin0=10000;
	for(j=0;j<cont->num;j++)
	  if(ABS(colmin[0]-cont->col[j])+ABS(rowmin[0]-cont->row[j])<dmin0)
	    dmin0= ABS(colmin[0]-cont->col[j])+ABS(rowmin[0]-cont->row[j]);
	
      }
      if(dmin0<dmin){index=i;
      dmin = dmin0;}
    }
  }

  printf("index=%d\t dmin=%lf \n",index,dmin);
  for (i=0; i<num_cont; i++) {
    cont = (*contours)[i];
    if(i!=index)cont->num =0;
  }
}

