#include "statclass.h"

void stat1class::affiche () {
	int i,j;
	cout<<" classe "<<(int)lab<<"\n";
	cout<<" vecteur moyenne : ";
	for (i=0; i<nbcanaux; i++) cout<<moy[i]<<" ";
	cout<<"\n";
	cout<<" matrice de covariance : \n";
	for (i=0; i<nbcanaux; i++) {
		cout<<"            ";
		for (j=0; j<nbcanaux; j++) cout<<Cov[i*nbcanaux+j]<<" ";
		cout<<"\n";
		}
	cout<<" matrice de covariance inversee : \n";
	for (i=0; i<nbcanaux; i++) {
		cout<<"            ";
		for (j=0; j<nbcanaux; j++) cout<<InvCov[i*nbcanaux+j]<<" ";
		cout<<"\n";
		}
	cout<<" determinant de la matrice de covariance : "<<detCov<<"\n";
	cout<<" probabilite a priori : "<<Pap<<"\n";
   }

void statclasses::ajoute (const stat1class & clas) {
//	cout<<" ajout de l'objet classe\n";
	stat1class * newclas=new stat1class(clas);
	newclas->suiv=debut;
	debut=newclas;
	nbclas++;
	}

bool statclasses::existe (BYTE lab) {
	courant=debut;
	bool trouve=0;
	while (courant != NULL && !trouve) {
		if (courant->lab==lab) trouve=1;
		else courant=courant->suiv;
		}
	return trouve;
	}

stat1class statclasses::extract (BYTE lab) {
	stat1class cl;
	courant=debut;
	bool trouve=0;
	while (courant != NULL && !trouve) {
		if (courant->lab==lab) {
			trouve=1;
			cl=*courant;
			}
		else courant=courant->suiv;
		}
	if (!trouve) {
		cout<<" classe "<<(int)lab<<" non trouvee parmi :\n";
		affiche ();
		}
	return cl;
	}

void statclasses::modif (stat1class & clas) {
	stat1class *cladetruire=NULL, *clprec=NULL, *newclas=new stat1class(clas);
	BYTE lab=newclas->label();
	courant=debut;
	bool trouve=0;
	while (courant != NULL && !trouve) {
		if (courant->lab==lab) {
			trouve=1;
			cladetruire=courant;
			courant=courant->suiv;
			delete cladetruire;
			newclas->suiv=courant;
			if (clprec!=NULL) clprec->suiv=newclas;
			else debut=newclas;
			}
		else {
			clprec=courant;
			courant=courant->suiv;
			}
		}
	if (!trouve) {
		cout<<" classe "<<(int)lab<<" non trouvee parmi :\n";
		affiche ();
		}
//	return cl;
	}

void statclasses::affiche () {
	cout<<nbclas<<" classes a afficher\n";
	courant=debut;
	while (courant != NULL) {
		courant->affiche();
		courant=courant->suiv;
		}
	}

pixel<float> statclasses::distclasses(const pixel<float> &pix) {
	stat1class cl;
	int i,j,nd,lab,k=0;
	float xx, xmin=1.e+9;
	pixel<float> Tdist(nbclas+1);
	posdebut();
	while (courant!=NULL && k<nbclas) {
		k++;
		cl=cl_courante();
		nd=cl.nbcanaux;
		float *X=new float[nd], *Y=new float[nd];
		lab=(int)cl.lab;
		if (lab<=0 || lab>nbclas) {
			cout<<" PB lab<>[1,nbclas] : lab="<<lab;
			cout<<" nbclas = "<<nbclas<<"\n"; cl.affiche();
			lab=k;
			}
		xx=0.;
		for (i=0; i<nd; i++) X[i]=pix[i]-cl.moy[i];
		for (i=0; i<nd; i++) {
			Y[i]=0.f;
			for (j=0; j<nd; j++) Y[i]+=X[j]*cl.InvCov[i*nd+j];
			xx+=X[i]*Y[i];
			}
		Tdist[lab]=xx;
		if (xx<xmin) {
			Tdist[0]=(float)lab;
			xmin=xx;
			}
		if (X!=NULL) delete[] X;
		if (Y!=NULL) delete[] Y;
		}
	return Tdist;
	}
/*
pixel<double> statclasses::distclasses(const pixel<double> &pix) {
	stat1class cl;
	int i,j,nd,lab,k=0;
	double xx, xmin=1.e+9;
	pixel<double> Tdist(nbclas+1);
	posdebut();
	while (courant!=NULL && k<nbclas) {
		k++;
		cl=cl_courante();
		nd=cl.nbcanaux;
		double *X=new double[nd], *Y=new double[nd];
		lab=(int)cl.lab;
		if (lab<=0 || lab>nbclas) {
			cout<<" PB lab<>[1,nbclas] : lab="<<lab;
			cout<<" nbclas = "<<nbclas<<"\n"; cl.affiche();
			lab=k;
			}
		xx=0.;
		for (i=0; i<nd; i++) X[i]=pix[i]-cl.moy[i];
		for (i=0; i<nd; i++) {
			Y[i]=0.f;
			for (j=0; j<nd; j++) Y[i]+=X[j]*cl.InvCov[i*nd+j];
			xx+=X[i]*Y[i];
			}
		Tdist[lab]=xx;
		if (xx<xmin) {
			Tdist[0]=lab;
			xmin=xx;
			}
		if (X!=NULL) delete[] X;
		if (Y!=NULL) delete[] Y;
		}
	return Tdist;
	}
*/
pixel<float> statclasses::U0gaus_classes(const pixel<float> &pix) {
//	const double PI=3.14159;
	stat1class cl;
	int i,j,nd,lab,k=0;
	float xx, xmin=1.e+9;
	pixel<float> Tdist(nbclas+1);
	posdebut();
	while (courant!=NULL && k<nbclas) {
		k++;
		cl=cl_courante();
		nd=cl.nbcanaux;
		float *X=new float[nd], *Y=new float[nd];
		lab=(int)cl.lab;
		if (lab<=0 || lab>nbclas) {
			cout<<" PB lab<>[1,nbclas] : lab="<<lab;
			cout<<" nbclas = "<<nbclas<<"\n";
			lab=k;
			}
		xx=0.;
		for (i=0; i<nd; i++) X[i]=pix[i]-cl.moy[i];
		for (i=0; i<nd; i++) {
			Y[i]=0.f;
			for (j=0; j<nd; j++) Y[i]+=X[j]*cl.InvCov[i*nd+j];
			xx+=X[i]*Y[i];
			}
//		Tdist[lab]=0.5*(xx+log(exp(nd*log(2.*PI))*cl.detCov));
		Tdist[lab]=(float)(0.5*(xx+nd*log(2.*PI*cl.detCov)));
		if (Tdist[lab]<xmin) {
			Tdist[0]=(float)lab;
			xmin=Tdist[lab];
			}
		if (X!=NULL) delete[] X;
		if (Y!=NULL) delete[] Y;
		}
	return Tdist;
	}

