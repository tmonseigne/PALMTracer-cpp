

#ifndef __contour_def__
#define __contour_def__
typedef enum {
  cont_no_junc,    
  cont_start_junc, 
  cont_end_junc,   
  cont_both_junc,  
  cont_closed      
} contour_class;



typedef struct {
  long  num;                
  float *row;               
  float *col;               
  float *angle;             
  float *response;          
  float *width_l;           
  float *width_r;           
  float *asymmetry;         
  float *contrast;          
  contour_class cont_class; 
} contour;
#endif

extern "C"
{

	int detectLines   (unsigned char *image, long width, long height, contour ***contours, float sigma, float low);
	void deleteContours(long num_cont, contour ***contours);

}