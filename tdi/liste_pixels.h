#ifndef _LISTE_PIXELS_H
#define _LISTE_PIXELS_H

#include "imabin.h"
#include <iostream>
using namespace std;

struct elt_liste {
	int x,y;
	elt_liste *suivant;
};

class liste_pixels {
	int nbelt, xmin, xmax, ymin, ymax;
	bool irect;
	elt_liste *debut, *courant, *precedent, *fin;
public:
	liste_pixels () {
//		cout<<" entree dans le constructeur de liste_pixels pour liste "<<this<<"\n";
		courant=precedent=fin=debut=NULL;
		nbelt=0;
		xmax=ymax=INT_MIN;
		xmin=ymin=INT_MAX;
		irect=0;
	}
	~liste_pixels () {
//		cout<<" entree dans le destructeur de liste_pixels pour liste "<<this<<"\n";
		courant=debut;
		while (courant!=NULL) {
			debut=courant->suivant;
//			cout<<" destruction de l'element d'adresse "<<courant<<"\n";
			delete courant;
			nbelt--;
			courant=debut;
		}
//		cout<<" fin destructeur\n";
	}
	liste_pixels (const liste_pixels &L) {
//		cout<<" entree dans le constructeur de copie de liste_pixels pour liste "<<this<<"\n"; //char aa; cin>>aa;
/*		courant=debut;
		while (courant!=NULL) {
			debut=courant->suivant;
			delete courant;
			nbelt--;
			courant=debut;
		}*/
//		cout<<" "<<L.nbelt; cout<<" "<<L.debut; 
		nbelt=0;
		debut=fin=precedent=NULL; 
		courant=L.debut; //cout<<L.nbelt<<" elements a copier depuis L.debut = @"<<L.debut<<"\n";
		while (courant!=NULL) {
			elt_liste *nouveau=new elt_liste;
			nouveau->x=courant->x;
			nouveau->y=courant->y;
			nouveau->suivant=debut;
			debut=nouveau;
			nbelt++;
			if (nbelt==1) fin=debut;
			courant=courant->suivant;
//			if (nbelt>892 && nbelt<899) cout<<" copie elt "<<nbelt<<" @"<<courant<<"\n";
/*			if (nbelt>L.nbelt) {
				cout<<" attention "<<nbelt<<" depasse le nbre d'elements de L "<<L.nbelt<<"\n";char aa; cin>>aa;
				L.affiche();
				cin>>aa;}*/
		}
//		cout<<" fin recopie "<<nbelt<<" elements liste \n";
		if (nbelt!=L.nbelt) cout<<" probleme dans recopie liste\n";
		xmin=L.xmin; xmax=L.xmax; ymin=L.ymin; ymax=L.ymax; irect=L.irect;
//		cout<<" fin recopie xmin etc. "<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<" "<<irect<<"\n";
		courant=precedent=debut;//cout<<" @deb "<<debut<<" "<<fin<<"\n";
	}
	liste_pixels& operator= (const liste_pixels& L) { //cout<<" operateur d'affectation liste\n";
//		cout<<" entree dans l'operateur d'affectation de liste_pixels pour liste "<<this<<"\n";
		if (this!=&L) {
			courant=debut;
			while (courant!=NULL) {
				debut=courant->suivant;
				delete courant;
				nbelt--;
				courant=debut;
			}
			nbelt=0;
			precedent=debut=fin=NULL;
			courant=L.debut;
			while (courant!=NULL) {
				elt_liste *nouveau;
				nouveau=new elt_liste;
				nouveau->x=courant->x;
				nouveau->y=courant->y;
				nouveau->suivant=debut;
				debut=nouveau;
				nbelt++;
				if (nbelt==1) fin=nouveau;
				courant=courant->suivant;
			}
			if (nbelt!=L.nbelt) cout<<" probleme dans recopie liste\n";
			xmin=L.xmin; xmax=L.xmax; ymin=L.ymin; ymax=L.ymax; irect=L.irect;
			courant=precedent=debut;
		}
		return *this;
	}
	int nb_elts () {return nbelt;}
	int l_min () {return xmin;}
	int l_max () {return xmax;}
	int c_min () {return ymin;}
	int c_max () {return ymax;}
/*	void vide () {//~liste_pixels();}
		courant=debut;
		while (courant!=NULL) {
			debut=courant->suivant;
			delete courant;
			nbelt--;
			courant=debut;
		}
	}*/
	liste_pixels& operator += (liste_pixels&);
	void insere (const elt_liste*);
	void insere (const elt_liste&);
	void insere (const int, const int);
	elt_liste extrait ();
	void supprime (int, int);
	bool existe (elt_liste*);
	bool existe (elt_liste&);
	bool existe (int, int);
	void affiche () const;
	void rectangle_englobant ();
	unsigned int nb_cc (liste_pixels*, unsigned short int=8, bool=0);
	friend class tree;
};

#endif