#ifndef _STATCLASS_H
#define _STATCLASS_H

#include <iostream>
using namespace std;
//#include <stdlib.h>
#include <math.h>
#include "constantes.h"
#include "def.h"
#include "pixel.h"

class stat1class                  // a valider en tant que class et non struct
{	BYTE lab;
	int nbcanaux;
	float *moy, *Cov, *InvCov;
	float detCov;
	float Pap;
	stat1class *suiv;
public:
	stat1class (BYTE l=0, int nd=1) {
//		cout<<" construction de l'objet classe\n";
		moy=Cov=InvCov=NULL; suiv=NULL;
		int i;
		lab=l;
		nbcanaux=nd;
		moy=new float [nbcanaux];
		for (i=0; i<nbcanaux; i++) moy[i]=0.;
		Cov=new float [nbcanaux*nbcanaux];
		for (i=0; i<nbcanaux*nbcanaux; i++) Cov[i]=0.;
		for (i=0; i<nbcanaux; i++) Cov[i*nbcanaux+i]=1.;
		InvCov=new float [nbcanaux*nbcanaux];
		for (i=0; i<nbcanaux*nbcanaux; i++) InvCov[i]=0.;
		for (i=0; i<nbcanaux; i++) InvCov[i*nbcanaux+i]=1.;
		detCov=1.;
		Pap=1.;
	}
	~stat1class () {
//		cout<<" destruction de l'objet classe\n";
		if (moy!=NULL) {delete[] moy; moy=NULL;}
		if (Cov!=NULL) {delete[] Cov; Cov=NULL;}
		if (InvCov!=NULL) {delete[] InvCov; InvCov=NULL;}
	}
	stat1class (const stat1class &cl) {
		moy=Cov=InvCov=NULL; suiv=NULL;
		int i;
		nbcanaux=cl.nbcanaux;
		lab=cl.lab;
		moy=new float[nbcanaux];
		for (i=0; i<nbcanaux; i++) moy[i]=cl.moy[i];
		Cov=new float[nbcanaux*nbcanaux];
		for (i=0; i<nbcanaux*nbcanaux; i++) Cov[i]=cl.Cov[i];
		InvCov=new float[nbcanaux*nbcanaux];
		for (i=0; i<nbcanaux*nbcanaux; i++) InvCov[i]=cl.InvCov[i];
		detCov=cl.detCov;
		suiv=NULL;
		Pap=cl.Pap;
	}
	stat1class& operator= (const stat1class &cl) {
		if (this != &cl) {
			int i;
			nbcanaux=cl.nbcanaux;
			lab=cl.lab;
			if (moy!=NULL) delete[] moy;
			moy=new float [nbcanaux];
			for (i=0; i<nbcanaux; i++) moy[i]=cl.moy[i];
			if (Cov!=NULL) delete[] Cov;
			Cov=new float [nbcanaux*nbcanaux];
			for (i=0; i<nbcanaux*nbcanaux; i++) Cov[i]=cl.Cov[i];
			if (InvCov!=NULL) delete[] InvCov;
			InvCov=new float[nbcanaux*nbcanaux];
			for (i=0; i<nbcanaux*nbcanaux; i++) InvCov[i]=cl.InvCov[i];
			detCov=cl.detCov;
			suiv=NULL;
			Pap=cl.Pap;
			}
		return *this;
	}
	float& mean (int k) {
    	if (k<0 || k>=nbcanaux) {
    		cout<<" debordement d'indice dans stat1class.mean "<<k<<"\n";
    		if (k<0) k=0;
    		if (k>=nbcanaux) k=nbcanaux-1;
    		}
	   	return moy[k];
   	}
	float& cova (int j, int k) {
    	if (j<0 || j>=nbcanaux || k<0 || k>=nbcanaux) {
    		cout<<" debordement d'indice dans stat1class.cova "<<j<<" "<<k<<"\n";
    		if (j<0) j=0;
    		if (j>=nbcanaux) j=nbcanaux-1;
	    	if (k<0) k=0;
    		if (k>=nbcanaux) k=nbcanaux-1;
   			}
	    return Cov[j*nbcanaux+k];
    }
	float& icova (int j, int k) {
    	if (j<0 || j>=nbcanaux || k<0 || k>=nbcanaux) {
    		cout<<" debordement d'indice dans stat1class.icova "<<j<<" "<<k<<"\n";
	    	if (j<0) j=0;
    		if (j>=nbcanaux) j=nbcanaux-1;
    		if (k<0) k=0;
    		if (k>=nbcanaux) k=nbcanaux-1;
    		}
	   	return InvCov[j*nbcanaux+k];
    }
	float& dcova () {
	    return detCov;
    }
	float& proba () {
	    return Pap;
    }
	BYTE& label () {
	    return lab;
    }
	void affiche ();
	friend class statclasses;
};

class statclasses
{	stat1class *debut, *courant;
	int nbclas;
 public:
	statclasses () {
//		cout<<" construction de l'objet statclasses\n";
		nbclas=0;
		debut=NULL;
		courant=debut;
		}
	~statclasses () {
//		cout<<" destruction de l'objet statclasses\n";
		stat1class *cladetruire;
		courant=debut;
		while (courant != NULL) {
			cladetruire=courant;
			courant=courant->suiv;
			delete cladetruire;
			nbclas--;
			}
		}
	void ajoute (const stat1class&);
	statclasses (const statclasses &stcl) {
		nbclas=0;
		debut=NULL;
		courant=stcl.debut;
		while (courant != NULL) {
			ajoute(*courant);
			courant=courant->suiv;
			}
		}
	statclasses& operator= (const statclasses & stcl) {  // a valider
		if (this != &stcl) {
			nbclas=0;
			debut=NULL;
			courant=stcl.debut;
			while (courant != NULL) {
				ajoute(*courant);
				courant=courant->suiv;
				}
			}
		return *this;
		}
	int nclas() {return nbclas;}
	void posdebut() { courant=debut;}
	stat1class cl_courante() {
		stat1class cl(*courant);
		courant=courant->suiv;
		return cl;
		}
	void affiche ();
	bool existe(BYTE);
	stat1class extract(BYTE);
	void statclasses::modif (stat1class &);
	pixel<float> distclasses(const pixel<float> &);
//	pixel<double> distclasses(const pixel<double> &pix);
	pixel<float> U0gaus_classes(const pixel<float> &);
};

#endif