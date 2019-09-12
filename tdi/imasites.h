#ifndef _IMASITES_H
#define _IMASITES_H

#include <iostream>
#include <iomanip>
using namespace std;
#include "def.h"
#include "fonctions.h"

class imasites
{protected:
	int nblig, nbcol;
	long int nbpix;
	void **Tsites;
public:
	imasites(const int nl=1, const int nc=1) : nblig(nl), nbcol(nc)  {//cout<<" constructeur imasites "<<nblig<<" "<<nbcol<<" \n";
		Tsites=NULL; nbpix=nblig*nbcol;
		Tsites=new void*[nbpix];
		for (long int i=0; i<nbpix; i++) Tsites[i]=NULL;
	}
	~imasites() {//cout<<" destructeur imasites \n";
		if (Tsites!=NULL) {delete[] Tsites; Tsites=NULL;}
	}
  imasites(const imasites &ima_s) {
		Tsites=NULL; Tsites=new void*[nbpix=ima_s.nbpix];
		nblig=ima_s.nblig; nbcol=ima_s.nbcol;
		nbpix=nblig*nbcol; for (long int i=0; i<nbpix; i++) Tsites[i]=NULL;
	}
	imasites& operator=(const imasites &ima_s) {
		if (this != &ima_s) {
			if (Tsites!=NULL) delete[] Tsites;
			nblig=ima_s.nblig; nbcol=ima_s.nbcol; nbpix=nblig*nbcol;
			Tsites=new void*[nbpix=ima_s.nbpix];
			for (long int i=0; i<nbpix; i++) Tsites[i]=ima_s.Tsites[i];
		}
		return *this;
	}
  void* & operator[] (long int i) const {
		if (i<0 || i>=nbpix) {
			cout<<" debordement d''indice "<<i<<"\n";
			if (i<0) i=0;
			if (i>=nbpix) i=nbpix-1;
		}
		return Tsites[i];
	}
	int nlig() const {return nblig;}
	int ncol() const {return nbcol;}
	void affiche() const;
	void sousIma(int,int,int=0,int=0,int=1,int=1);
	void zoomIma(int,int);
};
class imabin;
template <class T> class imadata;

class eltstruct {
protected:
	bool *aT;
	int nl, nc, x0, y0;
public:
	eltstruct (int nlig, int ncol, int yc, int xc) {
		aT=NULL;
		nl=nlig; nc=ncol;
		aT=new bool[nl*nc];
		for (int i=0; i<nl*nc; i++) aT[i]=1;
		x0=xc; y0=yc;
	}
	eltstruct (int nlig, int ncol) {
		aT=NULL;
		nl=nlig; nc=ncol;
		aT=new bool[nl*nc];
		for (int i=0; i<nl*nc; i++) aT[i]=1;
		x0=nc/2; y0=nl/2;
	}
	eltstruct (int n, double thetadeg) {
		aT=NULL;
		nl=n; nc=n;
		aT=new bool[nl*nc];
		int i,j,r=(n-1)/2,inc_j,inc_i;
		for (i=0; i<nl*nc; i++) aT[i]=0;
		x0=nc/2; y0=nl/2;
		aT[y0*nc+x0]=1;
		float dtheta=atan2(1.f,(float)r), thetaI=-(float)PI/2.f, thetaS=thetaI+dtheta, theta_approx, theta=(float)(thetadeg/180.*PI), errtheta, 
					theta2, errtheta2;
		if (theta>PI/2) theta-=(float)PI;     // theta est en radian et \in[-PI/2,PI/2] pour la suite des calculs
		bool trouve=0;                                               
		while (!trouve) {
			if (theta>=thetaI && theta<thetaS) { //cout<<theta*180/PI<<" "<<thetaI*180/PI<<" "<<thetaS*180/PI<<" -> on est dans l'intervalle\n";
				i=0; j=0; theta_approx=0.;
				if (fabs(theta)<mini(fabs((float)(PI/2.)-theta),fabs((float)(PI/2.)+theta))) {//cout<<" -> theta + proche de 0 -> increment index colonnes\n";
					while (i<r) {
						i++; inc_j=+1; theta_approx=atan2((float)(j+inc_j),float(i)); 
						if (theta_approx>PI/2) theta_approx-=(float)PI; //theta_approx \in[-PI/4,PI/4], theta \in[-PI/4,PI/4] 
						errtheta=fabs(theta_approx-theta); //theta_approx-theta \in[-PI/2,PI/2] -> |theta_approx-theta| \in[0,PI/2]
						theta2=atan2((float)j,(float)i); if (theta2>PI/2) theta2-=(float)PI; errtheta2=fabs(theta2-theta);
						if (errtheta2<errtheta) {inc_j=+0; /*theta_approx=atan2((float)(j+inc_j),(float)i); if (theta_approx>PI/2) theta_approx-=(float)PI;*/}
						theta2=atan2((float)(j-1),(float)i); if (theta2>PI/2) theta2-=(float)PI; errtheta2=fabs(theta2-theta);
						if (errtheta2<errtheta) {inc_j=-1; /*theta_approx=atan2((float)(j+inc_j),(float)i); if (theta_approx>PI/2) theta_approx-=(float)PI;*/}
						j+=inc_j; aT[(y0-j)*nc+(x0+i)]=1; aT[(y0+j)*nc+(x0-i)]=1;
					}
				} else {
					if (theta<0) theta+=(float)PI;     // theta \in[+PI/4,3*PI/4]
					while (j<r) {//cout<<" -> theta + proche de +/-PI/2 -> increment index lignes\n";
						j++; inc_i=+1; theta_approx=atan2((float)(j),float(i+inc_i)); 
						if (theta_approx<0) theta_approx+=(float)PI; //theta_approx \in[+PI/4,3*PI/4]
						errtheta=fabs(theta_approx-theta); //theta_approx-theta \in[-PI/2,PI/2] -> |theta_approx-theta| \in[0,PI/2]
						theta2=atan2((float)j,(float)i); if (theta2<0) theta2+=(float)PI; errtheta2=fabs(theta2-theta);
						if (errtheta2<errtheta) {inc_i=+0; /*theta_approx=atan2((float)(j),(float)i+inc_i);*/}
						theta2=atan2((float)j,(float)(i-1)); if (theta2<0) theta2+=(float)PI; errtheta2=fabs(theta2-theta);
						if (errtheta2<errtheta) {inc_i=-1; /*theta_approx=atan2((float)(j),(float)i+inc_i);*/}
						i+=inc_i; aT[(y0-j)*nc+(x0+i)]=1; aT[(y0+j)*nc+(x0-i)]=1;
					}
				}
				trouve=1;
			} else {
				thetaI=thetaS;
				thetaS+=dtheta;
			}
		}
		affiche();
	}
	~eltstruct () {
		if (aT!=NULL) {delete[] aT; aT=NULL;}
	}
	eltstruct (const eltstruct &E) {
		aT=NULL;
		nl=E.nl; nc=E.nc;
		aT=new bool[nl*nc];
		for (int i=0; i<nl*nc; i++) aT[i]=E.aT[i];
		x0=E.x0; y0=E.y0;
	}
	eltstruct& operator= (const eltstruct &E) {
		if (this!=&E) {
			if (aT!=NULL) delete[] aT;
			nl=E.nl; nc=E.nc;
			aT=new bool[nl*nc];
			for (int i=0; i<nl*nc; i++) aT[i]=E.aT[i];
			x0=E.x0; y0=E.y0;
		}
		return *this;
	}
 	bool& operator() (int i, int j) const {
 		return aT[i*nc+j];
 	}
 	eltstruct transpose() {
 		eltstruct Bt(nl,nc,nl-1-y0,nc-1-x0);
 		int i,j,j0,jj0;
 		for (i=0; i<nl; i++) {
 			j0=i*nc;
 			jj0=(nl-1-i)*nc+nc-1;
 			for (j=0; j<nc; j++) Bt.aT[jj0-j]=aT[j0+j];
 		}
 		return Bt;
 	}
	void mise_a_zero() {
		for (int i=0; i<nl; i++)
			for (int j=0; j<nc; j++) aT[i*nc+j]=0;
	}
 	void affiche() {
 		for (int i=0; i<nl; i++) {
 			for (int j=0; j<nc; j++) cout<<setw(2)<<aT[i*nc+j];
 			cout<<"\n";
 		}
 		cout<<" origine en ("<<y0<<", "<<x0<<")\n";
 	}
	friend class imabin;
	friend class imadata<bool>;
	friend class imadata<BYTE>;
	friend class imadata<int>;
	friend class imadata<float>;
	friend class imadata<double>;
};

#endif