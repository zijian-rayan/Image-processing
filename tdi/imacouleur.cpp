#include "imacouleur.h"



/* --------------------------------------------
Conversion RVB -> IST (intensite saturation teinte)
----------------------------------------------*/
void imacouleur:: RVB2IST () {
		int i,j,k;
		double r, v, b, xx, aa, bb, t, s;
		for (i=0; i<nblig; i++)
			for(j=0; j<nbcol; j++) {
				r=(*this)(i,j,0); v=(*this)(i,j,1); b=(*this)(i,j,2);
				xx=r+v+b; 
				(*this)(i,j,0)=xx/3.;
				if (xx>0) s=1.-mini(mini(r,v),b)/xx;
				else s=0.;
				(*this)(i,j,1)=255*s;
				aa=2.*r-v-b;
				bb=sqrt(3.)*(v-b);
				if (aa!=0) t=atan(bb/aa);
				else t=atan(bb/2.);
				if(t>PI/2.) t=t-PI/2.;
				if(t<(-PI/2.)) t=PI/2.-t;
				(*this)(i,j,2)=((t+PI/2.)*(255./PI));
			}
}
/* --------------------------------------------
fonction utile a la conversion RGB->LAB
----------------------------------------------*/
/*float f_lab(float val)
{
	float res;
	if(val >0.008856) res=(float)pow(val, 1./3);
	else res=(float)7.787*val+16./116;
	return(res);
}*/
/* --------------------------------------------
Conversion RVB -> Lab (intensite saturation teinte)
----------------------------------------------*/
/*imacouleur imacouleur::RVB2Lab () {
	int i,j, R,G,B;
	float x,y,z,L,a,b;
		for (i=0; i<nblig; i++)
			for(j=0; j<nbcol; j++) {
				R=(*this)(i,j,0);
				G=(*this)(i,j,1);
				B=(*this)(i,j,2);
				// conversion RGB->XYZ
				x=0.412453*R+  0.357580*G+  0.180423*B;
				y=0.212671*R+  0.715160*G+  0.072169*B;
				z=0.019334*R+  0.119193*G+  0.950227*B;
				// X/Xn, Y/Yn, Z/Zn
				x/=(0.950456*255);
				y/=255;
				z/=(1.088754*255);
				// conversion xyz->Lab
				if(y>0.008856) L= (int)(116*pow(y, 1./3)-16);
				else L= (int)floor(903.3*pow(y, 1./3)) ;
				a=(int)floor(500*(f_lab(x)-f_lab(y))+128);
				b=(int)floor(200*(f_lab(y)-f_lab(z))+128);
				(*this)(i,j,0) =L;
				(*this)(i,j,1) =a;
				(*this)(i,j,2) =b;

			}
		}*/
/* --------------------------------------------
Conversion RVB -> I1I2I3 (Ohta)
----------------------------------------------*/
/*imacouleur imacouleur:: RVB2I1I2I3 (){
	int i,j, R, G, B;
	int I1, I2, I3;
	for (i=0; i<nblig; i++)
		for(j=0; j<nbcol; j++) {
			R=(*this)(i,j,0);
			G=(*this)(i,j,1);
			B=(*this)(i,j,2);
			I1=(R+G+B)/3;
			I2=(R-B);
			I3=(2G-(R+B))/2;
			(*this)(i,j,0)=I1;
			(*this)(i,j,1)=I2;
			(*this)(i,j,2)=I3;
		}

}*/
/* --------------------------------------------
Conversion RVB -> Norme L1 (coordonnées normalisees)
----------------------------------------------*/
/*imacouleur imacouleur:: RVB2NormeL1 (){
	int i,j, R, G, B, I;
	for (i=0; i<nblig; i++)
		for(j=0; j<nbcol; j++) {
			R=(*this)(i,j,0);
			G=(*this)(i,j,1);
			B=(*this)(i,j,2);
			I=R+G+B+1;
			(*this)(i,j,0)=R*255/I;
			(*this)(i,j,1)=G*255/I;
			(*this)(i,j,2)=B*255/I;
		}

}*/