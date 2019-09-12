#ifndef _CIBLE_H
#define _CIBLE_H

#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <string>
#include <math.h>
#include <float.h>
using namespace std;


class cible
{	int nocib, x_centre, y_centre, dx_rect, dy_rect, x_ul_rect, y_ul_rect;
	double mod_vit, x_dir_vit, y_dir_vit;
	bool pos_valid, dim_valid, vit_valid;
	static const int valnul;
 public:
	cible(int _nocib=0, bool _pos_valid=0, int _x_centre=0, int _y_centre=0,
		  bool _dim_valid=0, int _dx_rect=0, int _dy_rect=0, int _x_ul_rect=0, int _y_ul_rect=0,
		  bool _vit_valid=0, double _x_dir_vit=0, double _y_dir_vit=0) : pos_valid(_pos_valid),  
		  dim_valid(_dim_valid), vit_valid(_vit_valid) {
		nocib=_nocib;
		if (_pos_valid) {
			x_centre=_x_centre;
			y_centre=_y_centre;
			}
		if (_dim_valid) {
			dx_rect=_dx_rect;
			dy_rect=_dy_rect;
			x_ul_rect=_x_ul_rect;
			y_ul_rect=_y_ul_rect;
			}
		if (_vit_valid) {
			x_dir_vit=_x_dir_vit;
			y_dir_vit=_y_dir_vit;
			mod_vit=pow(x_dir_vit*x_dir_vit+y_dir_vit*y_dir_vit,0.5);
			}
		}
	cible& operator=(const cible &cib) {
		if (this != &cib) {
			nocib=cib.nocib;
			pos_valid=cib.pos_valid; 
			x_centre=cib.x_centre;
			y_centre=cib.y_centre;
			dim_valid=cib.dim_valid;
			dx_rect=cib.dx_rect;
			dy_rect=cib.dy_rect;
			x_ul_rect=cib.x_ul_rect;
			y_ul_rect=cib.y_ul_rect;
			vit_valid=cib.vit_valid;
			mod_vit=cib.mod_vit;
			x_dir_vit=cib.x_dir_vit;
			y_dir_vit=cib.y_dir_vit;
			}
		return *this;
		}
	bool operator==(const cible &cib) {
		bool OK=1;
		if (this != &cib) {
			if ((nocib!=cib.nocib) ||  
				(pos_valid && cib.pos_valid && (x_centre!=cib.x_centre || y_centre!=cib.y_centre) ) ||
				(dim_valid && cib.dim_valid && (dx_rect!=cib.dx_rect || dy_rect!=cib.dy_rect ||
				 x_ul_rect!=cib.x_ul_rect || y_ul_rect!=cib.y_ul_rect) ) ||
				(vit_valid && cib.vit_valid && (x_dir_vit!=cib.x_dir_vit || y_dir_vit!=cib.y_dir_vit) ) )
				OK=0;
			}
		return OK;
		}
	int& no_cible () {return nocib;}
	bool centre_valide () const {return pos_valid;}
	void centre_invalide () {pos_valid=0;}
	int& xcentre () {
		if (!pos_valid) {
			cout<<" centre de la cible inconnu\n";
			x_centre=valnul;
			}
		return x_centre;
		}
	int& ycentre () {
		if (!pos_valid) {
			cout<<" centre de la cible inconnu\n";
			y_centre=valnul;
			}
		return y_centre;
		}
	void centre (const int _x_centre, const int _y_centre) {
		x_centre=_x_centre;
		y_centre=_y_centre;
		if (dim_valid && (x_centre<x_ul_rect || x_centre>x_ul_rect+dx_rect ||
			y_centre<y_ul_rect || y_centre>y_ul_rect+dy_rect) ) {
			cout<<" rectangle englobant ne contient pas le centre de la cible!\n";
			pos_valid=0;
			}
		else pos_valid=1;
		}
	bool rectangle_valide () const {return dim_valid;}
	void rectangle_invalide () {dim_valid=0;}
	int& dx_rectangle () {
		if (!dim_valid) {
			cout<<" rectangle englobant de la cible inconnu\n";
			dx_rect=valnul;
			}
		return dx_rect;
		}
	int& dy_rectangle () {
		if (!dim_valid) {
			cout<<" rectangle englobant de la cible inconnu\n";
			dy_rect=valnul;
			}
		return dy_rect;
		}
	int& x_ul_rectangle () {
		if (!dim_valid) {
			cout<<" rectangle englobant de la cible inconnu\n";
			x_ul_rect=valnul;
			}
		return x_ul_rect;
		}
	int& y_ul_rectangle () {
		if (!dim_valid) {
			cout<<" rectangle englobant de la cible inconnu\n";
			y_ul_rect=valnul;
			}
		return y_ul_rect;
		}
	void rectangle (const int _dx_rect, const int _dy_rect, const int _x_ul_rect, const int _y_ul_rect) {
		dx_rect=_dx_rect;
		dy_rect=_dy_rect;
		x_ul_rect=_x_ul_rect;
		y_ul_rect=_y_ul_rect;
		if (pos_valid && (x_centre<x_ul_rect || x_centre>x_ul_rect+dx_rect ||
			y_centre<y_ul_rect || y_centre>y_ul_rect+dy_rect) ) {
			cout<<" rectangle englobant ne contient pas le centre de la cible!\n";
			dim_valid=0;
			}
		else dim_valid=1;
		}
	bool vitesse_valide () const {return vit_valid;}
	void vitesse_invalide () {vit_valid=0;}
	double& x_vitesse () {
		if (!vit_valid) {
			cout<<" vitesse de la cible inconnue\n";
			x_dir_vit=valnul;
			}
		return x_dir_vit;
		}
	double& y_vitesse () {
		if (!vit_valid) {
			cout<<" vitesse de la cible inconnue\n";
			y_dir_vit=valnul;
			}
		return y_dir_vit;
		}
	double& mod_vitesse () {
		if (!vit_valid) {
			cout<<" vitesse de la cible inconnue\n";
			mod_vit=valnul;
			}
		return mod_vit;
		}
	void vitesse (const double _x_dir_vit, const double _y_dir_vit) {
		x_dir_vit=_x_dir_vit;
		y_dir_vit=_y_dir_vit;
		mod_vit=pow(x_dir_vit*x_dir_vit+y_dir_vit*y_dir_vit,0.5);
		vit_valid=1;
		}
//	void sup (const pixel<T> & pix) {
//		for (int i=0; i<nbcanaux; i++) 
//			if (pix.adval[i]>adval[i]) adval[i]=pix.adval[i];
//		}
	void affiche () const;
//	float distance (pixel<T> &);
	double distance_cibles (cible&);
};

#endif

#ifndef _LISTE_CIBLES_H
#define _LISTE_CIBLES_H

#include <iostream>
using namespace std;
#include "fonctions.h"

struct elt_cible {
	cible cib;
	elt_cible *suivant;
};

class liste_cibles {
	int nbelt;
	elt_cible *debut, *courant, *precedent;
public:
	liste_cibles () {
		debut=NULL;
		courant=debut;
		precedent=debut;
		nbelt=0;
		}
	~liste_cibles () {
		courant=debut;
		while (courant!=NULL) {
			debut=courant->suivant;
			delete courant;
			nbelt--;
			courant=debut;
			}
		}
	liste_cibles (const liste_cibles & L) {
		courant=debut;
		while (courant!=NULL) {
			debut=courant->suivant;
			delete courant;
			nbelt--;
			courant=debut;
		}
		nbelt=0;
		debut=NULL;
		courant=L.debut;
		while (courant!=NULL) {
			elt_cible *nouveau;
			nouveau=new elt_cible;
			nouveau->cib=courant->cib;
			nouveau->suivant=debut;
			debut=nouveau;
			nbelt++;
			}
		if (nbelt!=L.nbelt) cout<<" probleme dans recopie liste\n";
		courant=debut; //L.courant;
		precedent=debut; //L.precedent;
		}
	liste_cibles & operator= (const liste_cibles & L) {
		if (this!=&L) {
			courant=debut;
			while (courant!=NULL) {
				debut=courant->suivant;
				delete courant;
				nbelt--;
				courant=debut;
			}
			nbelt=0;
			debut=NULL;
			courant=L.debut;
			while (courant!=NULL) {
				elt_cible *nouveau;
				nouveau=new elt_cible;
				nouveau->cib=courant->cib;
				nouveau->suivant=debut;
				debut=nouveau;
				nbelt++;
				courant=courant->suivant;
				}
			if (nbelt!=L.nbelt) cout<<" probleme dans recopie liste\n";
			courant=debut; //L.courant;
			precedent=debut; //L.precedent;
			}
		return *this;
		}
	int nb_elts () {
		return nbelt;
		}
	void insere (cible *);
	void insere (cible &);
	cible extrait ();
	void supprime (int);
	bool existe (cible *);
	bool existe (cible &);
	bool existe (int);
	void affiche ();
	int label_max ();
	void liste_cibles::correspondances_cibles (liste_cibles &, bool **, bool=1);
};

#endif