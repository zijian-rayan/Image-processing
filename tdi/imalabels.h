#ifndef _IMALABELS_H
#define _IMALABELS_H
#define _CRT_RAND_S

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include <limits.h>

//#include "constantes.h"
#include "def.h"
#include "fonctions.h"
#include "pixel.h"
#include "imasites.h"
#include "imadata.h"
#include "imacouleur.h"
#include "statclass.h"
#include "imaBR.h"
#include "fichiers.h"

class imalabels : public imasites
{	int nblayer;
	BYTE *Tima;
	BYTE valnul;
 public:
/* constructeur par defaut */
	imalabels (int nl=1, int nc=1, int nd=1) : imasites(nl,nc) {
		Tima=NULL;
		long int i; int j;
		nblayer=nd;
		valnul=0;
		Tima=new BYTE[nbpix*nblayer];
		for (i=0; i<nbpix; i++) {
			for (j=0; j<nblayer; j++) Tima[i*nblayer+j]=valnul;
			Tsites[i]=&(Tima[i*nblayer]);
			}
		}

/* constructeur par conversion d'une image en BYTE en une image de labels */
	imalabels (imadata<BYTE> &ima, int k=0)  : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL;
		nblayer=1;
		Tima=new BYTE[nbpix*nblayer];
		valnul=0;
		int i,j; long int ii;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				Tima[ii]=(BYTE)ima(i,j,k);
				Tsites[ii]=&(Tima[ii]);
			}
	}

/* constructeur par simulation d'une realisation MRF */
	imalabels (int, int, int, float, float=1.);

/* constructeur par classifiation non supervisee EM d'une image des donnees */
	imalabels (imadata<float> &, int, char *algo="EMG");

/* constructeur par classifiation supervisee MRF d'une image des donnees */
	imalabels (imadata<float> &, statclasses, float=0., char *algo="ICM", float=0.);

/* constructeur par classifiation supervisee MRF d'une image des donnees optimisation par graph-cut */
//	imalabels (imadata<float> &, float, statclasses);

/* constructeur par classifiation supervisee MRF avec processus lignes d'une image des donnees */
	imalabels (imadata<float> &, statclasses, char, float, float a0=0.9, float a1=1.8, float a2=2.7);

/* constructeur par classifiation supervisee MRF d'une image des donnees */
	imalabels (imadata<float> &, double, double=20., float=0., /*char *algo="ICM", */float=0.);

/* constructeur par classifiation supervisee MRF d'une image des donnees BR */
	imalabels (int, int, imadata<float> &, statclasses, float=0., char *algo="ICM");

/* constructeur correspondant à la classe majoritaire dans le cas d'une image de pixels mixtes */
	imalabels (imaBR<float> & ima) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL;
		nblayer=1;
		Tima=new BYTE[nbpix*nblayer];
		valnul=0;
		int i,j,k,ncl=ima.nclasses(),icl; long int ii;
		double x, a_max;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				x=0;
				for (k=0; k<ncl; k++)
					if ((x=ima.pct_cl(i,j,k))>a_max) {a_max=x; icl=k;}
				ii=i*nbcol+j;
				Tima[ii]=(BYTE)icl;
				Tsites[ii]=&(Tima[ii]);
			}
	}

/* destructeur */
	~imalabels () {
		if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
	}

/* constructeur de recopie */
	imalabels (const imalabels& ima) : imasites(ima) {
		Tima=NULL;
		nblayer=ima.nblayer;
		valnul=ima.valnul;
		Tima=new BYTE[nbpix*nblayer];
		BYTE *adval;
		for (long int i=0; i<nbpix; i++) {
			adval=(BYTE *)ima.Tsites[i];
			for (int j=0; j<nblayer; j++) Tima[i*nblayer+j]=*(adval+j);
			Tsites[i]=&(Tima[i*nblayer]);
		}
	}

/* operateur d'affectation */
	imalabels& operator= (const imalabels& ima) {
		if (this != &ima) {
			imasites *ad1, *ad2;
			ad1=this;
			ad2=(imasites *) &ima;
			*ad1=*ad2;
			if (Tima!=NULL) delete[] Tima;
			nblayer=ima.nblayer;
			valnul=ima.valnul;
			Tima=new BYTE[nbpix*nblayer];
			BYTE * adval;
			for (long int i=0; i<nbpix; i++) {
				adval=(BYTE *)ima.Tsites[i];
				for (int j=0; j<nblayer; j++) Tima[i*nblayer+j]=*(adval+j);
				Tsites[i]=&(Tima[i*nblayer]);
			}
		}
		return *this;
	}

/* # de couches */
	int nlayer() {return nblayer; }

/* operateur d'acces a la valeur (i,j,k) */
	BYTE& operator () (int i, int j, int k=0) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol || k<0 || k>=nblayer) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<","<<k<<")\n";
			if (i<0) i=0; if (i>=nblig) i=nblig-1;
			if (j<0) j=0; if (j>=nbcol) j=nbcol-1;
			if (k<0) k=0; if (k>=nblayer) k=nblayer-1;
		}
		BYTE * adval=(BYTE *)Tsites[i*nbcol+j];
		adval=adval+k;
		return *adval;
	}

/* conversion en imadata<BYTE> */
	imadata<BYTE> conv2imBYTE (int k=0) {
	int i,j,x, ncl=0;
	imadata<BYTE> ima_u(nblig,nbcol);
	BYTE fact=255;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) if ((*this)(i,j,k)>ncl) ncl=(*this)(i,j,k);
	cout<<ncl<<" classes detectees sur l'image des labels\n"; 
	if (ncl>1) fact=255/ncl;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			BYTE *adval=(BYTE*)Tsites[i*nbcol+j];
			x=int(*(adval+k));
			ima_u(i,j)=(BYTE)x*fact;
		}
	return ima_u;
	}

/* conversion en imacouleur */
	imacouleur<BYTE> conv2imFACOL (int k=0) {
	int i,j,x, ncl=0;
	imadata<BYTE> ima_u(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			BYTE *adval=(BYTE*)Tsites[i*nbcol+j];
			x=int(*(adval+k));
			ima_u(i,j)=(BYTE)x;
		}
	imacouleur<BYTE> imacol(ima_u,0);
	return imacol;
	}

/* copie d'une image des labels dans une autre */
	void copiecanal(int k, imalabels &ima, int icanal=0) {
		int i, j, imax=mini(nblig,ima.nblig), jmax=mini(nbcol,ima.nbcol);
		for (i=0; i<imax; i++)
			for (j=0; j<jmax; j++) (*this)(i,j,k)=ima(i,j,icanal);
	}

/* affichage */
	void affiche ();

/* creation de la matrice des probabilites d'appartenance à chaque classes en chaque pixel */
/*  sous forme d'une image des donnees 'nombre_de_classes' canaux */
	imadata<float> imaU0condClas(imadata<float> &, statclasses);

/* comptage des voisins pour chaque label en 4-connexite */
	void cmptvs4connex (int, int, int, int, int*);

/* comptage des voisins pour chaque label en 8-connexite */
	void cmptvs8connex (int, int, int, int, int*);

/* energie locale du processus ligne */
	double U_lineprocess (int,int,int,int,float,float,float);

/* comptage des voisins pour chaque label en 4-connexite avec processus ligne */
	void cmptvs4connexLP (int, int, int, int*);

/* extension en colonnes de la classe icl_m a partir du pixel courant de label courant */
	int HorizLength (int, int, int, int&, int=0);

/* operateur 'difference d'images des labels' */
	imalabels operator - (imalabels&);

/* mise en correspondance des labels de 2 images par minimisation du nombre de pixels de label different */
	void correspondances_labels (imalabels&, bool **T, bool=1);

/* ecriture dans un fichier de sortie */
	void sauve_Ima (char * nomfich="imaLabsauvee.dat", int=0) const;
	void sauve3d_Ima(char * nomfich="imaLab3dsauvee.dat", int=0) const;
	void sauve1d_Ima(char * nomfich="imaLab1dsauvee.dat", int=0) const;
};

#endif