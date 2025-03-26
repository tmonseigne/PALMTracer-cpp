

/**********************************************************************************************************

  To use CImg with stk, to display...


**********************************************************************************************************/


#include "STKCImg.h"




using namespace cimg_library;
// The line below is necessary for compilation with non-standart C++ compilers.
#if !(( defined(_MSC_VER) && _MSC_VER<=1200 ) || defined(__DMC__))
using namespace std;
#endif

int display(USrtArray & input, double &gamma, char* windowName, contour ***contours, long num_cont, long minNbPoints){

	long x,y,c;
	unsigned short countClick=0, XMax, YMax;

	if(input.nbDim==2){
		XMax=input.size(0);
		YMax=input.size(1);
	}else{
		XMax=input.size(1);
		YMax=input.size(2);
	}


	CImg<float> cImg(XMax+30,__max(30*2,YMax),1,3,0), img;
	double mini,maxi, val;

	mini=maxi=pow(*input.data, gamma);
	for(x=0;x<XMax ;x++)
	for(y=0;y<YMax ;y++)
	for(c=0;c<3 ;c++){
		val=pow(input.data[y*XMax+x], gamma);
		if(mini>val) mini=val;
		if(maxi<val) maxi=val;
		cImg(x+30,y,0,c) = (float)val; 
	}
	for(x=0;x<30 ;x++)
	for(y=0;y<60 ;y++)
	for(c=0;c<3 ;c++){
		if (x%30==0 || y%30==0 || x%30==29 || y%30==29)
			cImg(x,y,0,c) = (float)(mini+(maxi-mini)/2);
		else
			cImg(x,y,0,c) = (float)mini;
	}

	for(x=5;x<25 ;x++)
	for(c=0;c<3 ;c++){
		cImg(x,15,0,c)=cImg(x,45,0,c)=cImg(15,x,0,c)=(float)maxi;
	}

	//display the coutours
	long num_points;
	int i, j;
	contour *cont;
	for (i=0; i<num_cont; i++) {
		cont = (*contours)[i];
		num_points = cont->num;
		if(num_points>minNbPoints)
			for (j=0; j<num_points; j++) {
				cImg( (int)cont->col[j]+30,(int)cont->row[j] ,0,0)=(float)maxi;
				cImg( (int)cont->col[j]+30,(int)cont->row[j],0,1)=0;
				cImg( (int)cont->col[j]+30,(int)cont->row[j],0,2)=0;
			}
	}


	double zoom=1;
    
	CImgDisplay disp(cImg,windowName,1,3);
	bool changed=false;
	while (!(disp.closed || disp.key == cimg::keyQ || disp.key==cimg::keyESC || disp.key==cimg::keyENTER)){
		//printf("%d %d \n",disp.mousex, disp.mousey);
		if(disp.button!=0){
			if((disp.mousex>=0)&&(disp.mousey>=0)){
				if(disp.mousex/zoom<30){
					if(disp.mousey/zoom<30){
						gamma*=0.8;
						//printf("click  +  en %d %d gamma %f\n",disp.mousex, disp.mousey, gamma);
						changed=true;
					}else{
						if(disp.mousey/zoom<60){
							gamma*=1.2;
							//printf("click  -  en %d %d gamma %f\n",disp.mousex, disp.mousey, gamma);
							changed=true;
						}
					}
					if(changed){
						mini=maxi=pow(*input.data, gamma);
						for(x=0;x<XMax ;x++)
						for(y=0;y<YMax ;y++){
							val=pow(input.data[y*XMax+x], gamma);
							for(c=0;c<3 ;c++){
								if(mini>val) mini=val;
								if(maxi<val) maxi=val;
								cImg(x+30,y,0,c) = (float)val; 
							}
						}
						for(x=0;x<30 ;x++)
						for(y=0;y<60 ;y++)
						for(c=0;c<3 ;c++){
							if (x%30==0 || y%30==0 || x%30==29 || y%30==29)
								cImg(x,y,0,c) = (float)(mini+(maxi-mini)/2);
							else
								cImg(x,y,0,c) = (float)mini;
						}

						for(x=5;x<25 ;x++)
						for(c=0;c<3 ;c++){
							cImg(x,15,0,c)=cImg(x,45,0,c)=cImg(15,x,0,c)=(float)maxi;
						}


						for (i=0; i<num_cont; i++) {
							cont = (*contours)[i];
							num_points = cont->num;
							if(num_points>minNbPoints)
								for (j=0; j<num_points; j++) {
									cImg( (int)cont->col[j]+30,(int)cont->row[j],0,0)=(float)maxi;
									cImg( (int)cont->col[j]+30,(int)cont->row[j],0,1)=0;
									cImg( (int)cont->col[j]+30,(int)cont->row[j],0,2)=0;
								}
						}
					}
				}
				
			}
			changed=false;
			disp.button=0;
			disp.resize(disp.dimx(),disp.dimy()).display(cImg);
		}
		else{
			double dwx=disp.window_dimx();
			double dx=disp.dimx();
			double dy=disp.dimy();
			if (fabs(dwx-dx)>1e-10){
				zoom=zoom*dwx/dx;
				disp.resize((int)dwx,(int)(dy*dwx/dx)).display(cImg);
			}else
				cImg.wait(25);
		}
	}
	return(RET_OK);
}
int display(USrtArray & input, double &gamma, char* windowName, unsigned short nClick, unsigned short *getClick){


	long x,y,c;
	unsigned short countClick=0, XMax, YMax;

	if(input.nbDim==2){
		XMax=input.size(0);
		YMax=input.size(1);
	}else{
		XMax=input.size(1);
		YMax=input.size(2);
	}


	CImg<float> cImg(XMax+30,__max(30*2,YMax),1,3,0), img;
	double mini,maxi, val;

	mini=maxi=pow(*input.data, gamma);
	for(x=0;x<XMax ;x++)
	for(y=0;y<YMax ;y++)
	for(c=0;c<3 ;c++){
		val=pow(input.data[y*XMax+x], gamma);
		if(mini>val) mini=val;
		if(maxi<val) maxi=val;
		cImg(x+30,y,0,c) = (float)val; 
	}
	for(x=0;x<30 ;x++)
	for(y=0;y<60 ;y++)
	for(c=0;c<3 ;c++){
		if (x%30==0 || y%30==0 || x%30==29 || y%30==29)
			cImg(x,y,0,c) = (float)(mini+(maxi-mini)/2);
		else
			cImg(x,y,0,c) = (float)mini;
	}

	for(x=5;x<25 ;x++)
	for(c=0;c<3 ;c++){
		cImg(x,15,0,c)=cImg(x,45,0,c)=cImg(15,x,0,c)=(float)maxi;
	}

	double zoom=1;
    
	CImgDisplay disp(cImg,windowName,1,3);
	bool changed=false;
	while (!(disp.closed || disp.key == cimg::keyQ || disp.key==cimg::keyESC || disp.key==cimg::keyENTER)){
		//printf("%d %d \n",disp.mousex, disp.mousey);
		if(disp.key == cimg::keyDELETE || disp.key == cimg::keyBACKSPACE)
			if(countClick>0){
				countClick--;
				for(x=-10;x<10 ;x++)
					for(c=0;c<3 ;c++){
						cImg(x+getClick[countClick*2]+30,getClick[countClick*2+1],0,c)=(float)pow(input.data[getClick[countClick*2+1]*XMax+x+getClick[countClick*2]], gamma);
					}
				for(x=-10;x<10 ;x++)
					for(c=0;c<3 ;c++){
						cImg(getClick[countClick*2]+30,x+getClick[countClick*2+1],0,c)=(float)pow(input.data[(x+getClick[countClick*2+1])*XMax+getClick[countClick*2]], gamma);
					}
				disp.resize(disp.dimx(),disp.dimy()).display(cImg);
				disp.key=0;
			}
		if(disp.button!=0){
			if((disp.mousex>=0)&&(disp.mousey>=0)){
				if(disp.mousex/zoom<30){
					if(disp.mousey/zoom<30){
						gamma*=0.8;
						//printf("click  +  en %d %d gamma %f\n",disp.mousex, disp.mousey, gamma);
						changed=true;
					}else{
						if(disp.mousey/zoom<60){
							gamma*=1.2;
							//printf("click  -  en %d %d gamma %f\n",disp.mousex, disp.mousey, gamma);
							changed=true;
						}
					}
					if(changed){
						mini=maxi=pow(*input.data, gamma);
						for(x=0;x<XMax ;x++)
						for(y=0;y<YMax ;y++){
							val=pow(input.data[y*XMax+x], gamma);
							for(c=0;c<3 ;c++){
								if(mini>val) mini=val;
								if(maxi<val) maxi=val;
								cImg(x+30,y,0,c) = (float)val; 
							}
						}
						for(x=0;x<30 ;x++)
						for(y=0;y<60 ;y++)
						for(c=0;c<3 ;c++){
							if (x%30==0 || y%30==0 || x%30==29 || y%30==29)
								cImg(x,y,0,c) = (float)(mini+(maxi-mini)/2);
							else
								cImg(x,y,0,c) = (float)mini;
						}

						for(x=5;x<25 ;x++)
						for(c=0;c<3 ;c++){
							cImg(x,15,0,c)=cImg(x,45,0,c)=cImg(15,x,0,c)=(float)maxi;
						}

						for(int i=0; i<countClick;i++){
						for(x=-10;x<10 ;x++)
							for(c=0;c<3 ;c++){
								cImg(x+getClick[i*2]+30,getClick[i*2+1],0,c)=cImg(getClick[i*2]+30,x+getClick[i*2+1],0,c)=(float)maxi;
							}
						}
					}
				}else{
					if ((nClick>0)&&(countClick<nClick)){
						getClick[countClick*2]=(unsigned short)(disp.mousex/zoom-30);
						getClick[countClick*2+1]=(unsigned short)(disp.mousey/zoom);
						//printf("click -> %d %d\n", getClick[countClick*2], getClick[countClick*2+1]);
						
						for(x=-10;x<10 ;x++)
							for(c=0;c<3 ;c++){
								cImg(x+getClick[countClick*2]+30,getClick[countClick*2+1],0,c)=cImg(getClick[countClick*2]+30,x+getClick[countClick*2+1],0,c)=(float)maxi;
							}
						countClick++;
					}
				}
				
			}
			changed=false;
			disp.button=0;
			disp.resize(disp.dimx(),disp.dimy()).display(cImg);
		}
		else{
			double dwx=disp.window_dimx();
			double dx=disp.dimx();
			double dy=disp.dimy();
			if (fabs(dwx-dx)>1e-10){
				zoom=zoom*dwx/dx;
				disp.resize((int)dwx,(int)(dy*dwx/dx)).display(cImg);
			}else
				cImg.wait(25);
		}
	}
	return(countClick);
}