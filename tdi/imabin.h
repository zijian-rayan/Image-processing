#ifndef _IMABIN_H
#define _IMABIN_H

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
//#include <stdlib.h>
#include <math.h>
//#include "constantes.h"
#include "def.h"
#include "fonctions.h"
#include "pixel.h"
#include "liste_pixels.h"
#include "imasites.h"
#include "imadata.h"

class imabin : public imasites
{	bool *Tima;
	void init_zero() { // attention pas de désallocation mémoire car utilisée dans constructeurs
		Tima=NULL; Tima=new bool[nbpix];
		for (long int i=0; i<nbpix; i++) {
			Tima[i]=0;
			Tsites[i]=&(Tima[i]);
		}
	}
	void re_init_zero(const int nl,const int nc) { // avec désallocation mémoire
		if (nblig!=nl || nbcol!=nc) {
			nblig=nl; nbcol=nc;
			if (Tsites!=NULL) {delete[] Tsites; Tsites=NULL;}
			Tsites=new void*[nbpix=(long int)nl*(long int)nc];
			for (long int i=0; i<nbpix; i++) Tsites[i]=NULL;
		}
		if (Tima!=NULL) delete[] Tima;
		init_zero();
	}
 public:
	imabin(const int nl=1, const int nc=1) : imasites(nl,nc) {
		init_zero();
	}
	~imabin() {
		if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
	}
	imabin(const imabin &ima) : imasites(ima) {
		bool *adval;
		init_zero();
		for (long int i=0; i<nbpix; i++) {
			adval=(bool*)ima.Tsites[i];
			Tima[i]=*adval;
		}
	}
	imabin(const imadata<float> &ima, const float seuil, const int k=0) : imasites(ima) {
		int i,j; long int n;
		init_zero();
		for (i=0; i<nblig; i++) {
			n=i*nbcol;
			for (j=0; j<nbcol; j++) {
				if (ima.valpix(i,j,k)>=seuil && n>=0 && n<nbpix) Tima[n]=1;
				n++;
			}
		}
	}
	imabin(const imadata<int> &ima, const float seuil, const int k=0) : imasites(ima) {
		int i,j; long int n;
		init_zero();
		for (i=0; i<nblig; i++) {
			n=i*nbcol;
			for (j=0; j<nbcol; j++) {
				if (ima.valpix(i,j,k)>=seuil && n>=0 && n<nbpix) Tima[n]=1;
				n++;
			}
		}
	}
	imabin(const imadata<BYTE> &ima, const float seuil, const int k=0) : imasites(ima) {
		int i,j; long int n;
		init_zero();
		for (i=0; i<nblig; i++) {
			n=i*nbcol;
			for (j=0; j<nbcol; j++) {
				if (ima.valpix(i,j,k)>=seuil && n>=0 && n<nbpix) Tima[n]=1;
				n++;
			}
		}
	}
	imabin& operator=(const imabin &ima) {
		if (this != &ima) {
			imasites *ad1, *ad2;
			ad1=this;
			ad2=(imasites*) &ima;
			*ad1=*ad2;
			if (Tima!=NULL) delete[] Tima;
			Tima=new bool[nbpix];
			bool *adval;
			for (long int i=0; i<nbpix; i++) {
				adval=(bool*)ima.Tsites[i];
				Tima[i]=*adval;
				Tsites[i]=&(Tima[i]);
			}
		}
		return *this;
	}
	bool& operator () (int i, int j) const {
		if (i<0 || i>=nblig || j<0 || j>=nbcol) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
			if (j>=nbcol) j=nbcol-1;
		}
		bool *adval=(bool*)Tsites[i*nbcol+j];
		return *adval;
	}
	bool& operator () (int i, int j) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
			if (j>=nbcol) j=nbcol-1;
		}
		bool *adval=(bool*)Tsites[i*nbcol+j];
		return *adval;
	}
	operator imadata<float>();
	operator imadata<BYTE>();
	imadata<BYTE> imaunsignedchar (const bool=1);
	void affiche (const int=1) const;
	imabin operator+ (const imabin &);
	imabin operator- (const imabin &);
	imabin operator&& (const imabin &);
	imabin operator|| (const imabin &);
	void mise_a_zero ();
	void mise_a_un ();
	imabin negatif ();
	imabin ajoutbords ();
	imabin retirebords ();
	bool operator== (const imabin &);
	bool operator!= (const imabin &);
	float norm () const;
	void sauve_Ima (char *nomfich="imagesauvee.dat", bool=1) const;
	imadata<int> composantes_connexes (int &, int=4, bool=0);
	imadata<int> boxes_rectangulaires (int &, bool=1);
/* morphologie mathematique binaire */
	imadata<float> Tr_dist (const imabin &, int=1) const;
	imadata<float> Tr_dist (int=1) const;
	imabin dilate (eltstruct);
	imabin dilate (eltstruct, int);
	imabin dilate (eltstruct, eltstruct);
	imabin dilate (const imadata<float> &, float=1.);
	imabin dilate (int=3, float=1.);
	imabin erode (eltstruct);
	imabin erode (eltstruct, int);
	imabin erode (eltstruct, eltstruct);
	imabin erode (const imadata<float> &, float=1.);
	imabin erode (int=3, float=1.);
	imabin ouverture (eltstruct);
	imabin ouverture (eltstruct, int);
	imabin fermeture (eltstruct);
	imabin fermeture (eltstruct, int);
	imabin tophat (eltstruct);
	imabin tophat_c (eltstruct);
	imabin reconstruction_geodesique (imabin &, int=4);
	imabin reconstruction_geodesique (int, int, int=4);
	imabin bouche_trou (int=4);
	imabin bouche_trou2 (int=4);
	imadata<int> CompoConnexes_MM (int &, int=4, bool=1);
//	imabin erode_ultime (eltstruct);               // a ecrire et a valider
	imabin erode_ultime (const imadata<float> &, int=8);
	imabin erode_ultime (int=3);
	imabin transfo_tout_ou_rien (eltstruct, eltstruct, int=1) const;
	imabin zones_influence_geodesique(const imabin &, bool &, int=8);
	imadata<BYTE> detect_coin (int=1, bool=1, bool=0);
	imabin enveloppe_convexe ();
	imabin squelette (int=8, BYTE=1, bool=1);
	imabin squelette (imabin&, int=8, bool=1);
	imabin elagage (int=1, bool=0);
	imabin elagage_saufextrem (int=1, int=8);
	imabin detect_extremites_sq ();  // s'applique à une image de squelette
	imadata<BYTE> calc_orient_sq (bool=0, bool=1); // s'applique à une image de squelette
/* analyse d'image binaire */
	imadata<float> transformee_Fourier();
	imadata<float> transformee_Fourier_Inv();
	imadata<int> transformee_Hough_1pt (bool=1); // spécifique pour espace cumulatif associé à calcul NFA
	imadata<int> transformee_Hough();
	void reconst_hough_transform (imadata<float> &, int, int);
	void reconst_hough_transform (imadata<float> &, int, int, imabin); 
	imadata<int> norme_transformee_Hough(int, int, int=1, int=1);
	imadata<int> norme_transformee_Hough (const imabin &);
	imadata<int> hough_transform();
//	imadata<int> analyse_hough (imadata<int>);
	imadata<int> transformee_Hough_cercles(const int=25);
	imadata<int> transformee_Hough_squares(const int=25);
	imadata<float> coocurrence (int, int, int=1, unsigned short int=0);
/* trace formes geometriques basiques */
	void trace_line (float, float, float=0.f, float=2*PI, bool=0);
	void trace_line (float, float, int, int, int, int);
	void trace_line (int, int, int, int);
};

#endif