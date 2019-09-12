#ifndef _IMACOULEUR_H
#define _IMACOULEUR_H

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <math.h>
//#include "constantes.h"
//#include "def.h"
//#include "fonctions.h"
//#include "pixel.h"
//#include "imasites.h"
#include "imadata.h"
//#include "statclass.h"


template <class T> class imacouleur : public imadata<T>
{	bool iRVB, iIST, iLAB, iI1I2I3, iL1; // on choisit un espace couleur dans lequel on va travailler
 public:
/* constructeur par defaut */
	imacouleur (int nl=1, int nc=1) : imadata<T>(nl,nc,3) {}
/* constructeur a partir d'une image de donnees a priori RVB */
	imacouleur (imadata<T> &imadon) : imadata<T>(imadon.nlig(),imadon.ncol(),3) {
		int i,j,k;
		int T_canaux[3];
		if (imadon.ncanaux()==3)
			for (k=0; k<3; k++) T_canaux[k]=k;
		else
			for (k=0; k<3; k++) T_canaux[k]=0;
		for (k=0; k<3; k++)
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					(*this)(i,j,k)=imadon(i,j,T_canaux[k]);
		iRVB=1;
		iIST=iLAB=iL1=iI1I2I3=0;
	}
// conversion de l'image de niveaux de gris en fausses couleurs
	imacouleur (imadata<T> &imadon, int icanal) : imadata<T>(imadon.nlig(),imadon.ncol(),3) {
		cout<<" conversion en fausses couleurs\n";
		int ican=icanal,i,j,k,n;
		if (ican<0 || ican>=imadon.ncanaux()) {
			ican=maxi(0,mini(icanal,imadon.ncanaux()-1));
			cout<<" numero du canal traite = "<<ican<<" !!!\n";
		}
		int Tab[256][3];
		for (i=0; i<256; i++) for (j=0; j<3; j++) Tab[i][j]=0;
		int v_h=0, v_b;
		for (n=0; n<256-8; n+=8) {
			v_b=n; v_h=(int)(255*cos(v_h/255.f)); //cout<<n<<" "<<v_b<<" "<<v_h<<"\n";
			Tab[0+n][0]=Tab[0+n][1]=Tab[0+n][2]=v_b;
			Tab[1+n][0]=v_h; Tab[1+n][1]=Tab[1+n][2]=v_b;
			Tab[2+n][0]=v_b; Tab[2+n][1]=v_h; Tab[2+n][2]=v_b;
			Tab[3+n][0]=Tab[3+n][1]=v_b; Tab[3+n][2]=v_h;
			Tab[4+n][0]=Tab[4+n][1]=v_h; Tab[4+n][2]=v_b;
			Tab[5+n][0]=v_h; Tab[5+n][1]=v_b; Tab[5+n][2]=v_h;
			Tab[6+n][0]=v_b; Tab[6+n][1]=Tab[6+n][2]=v_h;
			Tab[7+n][0]=Tab[7+n][1]=Tab[7+n][2]=v_h;
		}
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				n=imadon(i,j,ican);
				for (k=0; k<3; k++) (*this)(i,j,k)=Tab[n][k];
			}
		iRVB=1;
		iIST=iLAB=iL1=iI1I2I3=0;
	}
// conversion RVB vers IST (Intensite, Saturation, Teinte)
	void RVB2IST () {
		int i,j;
		double r, v, b, xx, aa, bb, t, s;
		for (i=0; i<nblig; i++)
			for(j=0; j<nbcol; j++) {
				r=(*this)(i,j,0); v=(*this)(i,j,1); b=(*this)(i,j,2);
				xx=r+v+b;
				(*this)(i,j,0)=(T)(xx/3.);
				if (xx>0) s=1.-mini(mini(r,v),b)/xx;
				else s=0.;
				(*this)(i,j,1)=(T)(255*s);
				aa=2.*r-v-b;
				bb=sqrt(3.)*(v-b);
				if (aa!=0) t=atan(bb/aa);
				else t=atan(bb/2.);
				if(t>PI/2.) t=t-PI/2.;
				if(t<(-PI/2.)) t=PI/2.-t;
				(*this)(i,j,2)=(T)((t+PI/2.)*(255./PI));
				iL1=iRVB=iLAB=iI1I2I3=0;
				iIST=1;
			}
	}
// conversion RVB vers I (Intensite)
	imadata<T> RVB2I () {
		int i,j;
		imadata<T> imares(nblig,nbcol);
		for (i=0; i<nblig; i++)
			for(j=0; j<nbcol; j++) imares(i,j)=(T)(((double)(*this)(i,j,0)+(double)(*this)(i,j,1)+(double)(*this)(i,j,2))/3.);
		return imares;
	}
// conversion RVB vers Lab (Espace uniforme)
	float f_lab(float val){
		float res;
		if (val>0.008856) res=(float)pow(val,1.0f/3);
		else res=(float)(7.787*val+16./116.);
		return res;
	}

	void RVB2Lab () {
	int i,j, R,G,B;
	float x,y,z,L,a,b;
		for (i=0; i<nblig; i++)
			for(j=0; j<nbcol; j++) {
				R=(*this)(i,j,0);
				G=(*this)(i,j,1);
				B=(*this)(i,j,2);
				// conversion RGB->XYZ
				x=0.412453f*R+  0.357580f*G+  0.180423f*B;
				y=0.212671f*R+  0.715160f*G+  0.072169f*B;
				z=0.019334f*R+  0.119193f*G+  0.950227f*B;
				// X/Xn, Y/Yn, Z/Zn
				x/=(float)(0.950456*255);
				y/=(float)255;
				z/=(float)(1.088754*255);
				// conversion xyz->Lab
				if(y>0.008856f) L= (float)(116*pow(y, 1.0f/3)-16);
				else L= (float)floor(903.3*pow(y, 1.0f/3)) ;
				a=(float)floor(500*(f_lab(x)-f_lab(y))+128);
				b=(float)floor(200*(f_lab(y)-f_lab(z))+128);
				(*this)(i,j,0) =(T)L;
				(*this)(i,j,1) =(T)a;
				(*this)(i,j,2) =(T)b;
				iL1=iRVB=iIST=iI1I2I3=0;
				iLAB=1;
			}
		}
// conversion RVB vers I1I2I3 (Ohta, systeme d'axes independants)
	void RVB2I1I2I3 (){
	int i,j, R, G, B;
	int I1, I2, I3;
	for (i=0; i<nblig; i++)
		for(j=0; j<nbcol; j++) {
			R=(*this)(i,j,0);
			G=(*this)(i,j,1);
			B=(*this)(i,j,2);
			I1=(R+G+B)/3;
			I2=(R-B);
			I3=(2*G-(R+B))/2;
			(*this)(i,j,0)=I1;
			(*this)(i,j,1)=I2;
			(*this)(i,j,2)=I3;
			iRVB=iLAB=iIST=iL1=0;
			iI1I2I3=1;
		}

}
// conversion RVB vers norme L1 (espace normalise)
	void RVB2NormeL1 (){
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
			iRVB=iLAB=iIST=iI1I2I3=0;
			iL1=1;
		}
	}
// selection 1 composante parmi 3
	imadata<BYTE> select_composante (int no){
		imadata<BYTE> ima(nblig, nbcol);
		int i,j;
		for (i=0; i<nblig; i++)
			for(j=0; j<nbcol; j++) {
				ima(i,j)=(*this)(i,j,no);
			}
			return(ima);
	}

};

#endif