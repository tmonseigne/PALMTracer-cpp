#include <fstream.h>
#include <stdlib.h>
void read_data(char *,int *, double*,double*);

main(int argc,char* argv[])
{
int nbline;
double *xinit,*yinit;

xinit = new double[100000];
yinit = new double[100000];
read_data(argv[1],&nbline,xinit,yinit);

}

void read_data(char *nomfic, int *nbline, double *xinit, double *yinit)
{
 double a1,a2,a3,a4,a5;
 int i,j,k,npoint;
	char buf[100];
	ifstream infile(nomfic);	
	cout << nomfic << endl;	
	if(!infile){
		cout <<"Pb d'ouverture du fichier " <<nomfic <<endl;
	}
	for(i=0;i<2;i++){
		infile >> buf;
//	infile.read((char*)buf,100*sizeof(unsigned char),'\n');
		cout << buf[0] <<endl;
	}
	infile >> *nbline;

  for(i=0,k=0;i<*nbline;i++){
    infile.get(buf,100,'\n');
    infile >> npoint;
	//cout << npoint <<endl;
    for(j=0;j<npoint;j++,k++)
      infile >> yinit[k] >> xinit[k] >>a1 >> a2 >>a3>>a4>>a5;
  }

  cout << "nblin=" << *nbline <<"\t nbp=" << k<<endl;
   
  infile.close();

}
