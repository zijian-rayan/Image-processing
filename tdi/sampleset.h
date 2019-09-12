#ifndef _SAMPLESET_H
#define _SAMPLESET_H

#include "def.h"

template <class T> class sampleset;
class statclasses;

template <class T> class sample
{	T val;
	int no;
	BYTE lab;
	sample<T> *suivant;
 public:
	sample() {val=0; no=0; lab=0; suivant=NULL;}
	sample(T t, BYTE l=0, int n=0) {
		val=t; no=n; lab=l; suivant=NULL;
		}
	~sample() {}
	sample(const sample<T> &s) {
		val=s.val;
		no=s.no;
		lab=s.lab;
		suivant=s.suivant;
		}
	sample& operator=(const sample<T> &s) {
		if (this!=&s) {
			val=s.val;
			no=s.no;
			lab=s.lab;
			suivant=s.suivant;
			}
		return *this;
		}
	sample& operator=(const T &t) {
		val=t;
		no=0;
		lab=0;
		suivant=NULL;
		return *this;
		}
	void affiche() {
		cout<<" echantillon ";
		val.affiche(0);
		cout<<" numero = "<<no<<", label = "<<(int)lab<<"\n";
		}
	BYTE label() { return lab;}
	void label(BYTE l) { lab=l;}
	int numero() { return no;}
	void numero(int n) { no=n;}
	friend class sampleset<T>;
};

class imalabels;

template <class T> class sampleset
{   sample<T> * debut, * courant;
    int nbsample;
 public:
	sampleset () {
		nbsample=0;
		debut=NULL;
		courant=debut;
		}
	~sampleset () {
		sample<T> *adetruire;
		courant=debut;
		while (courant != NULL) {
			adetruire=courant;
			courant=courant->suivant;
			delete adetruire;
			nbsample--;
			}
		}
	sampleset (const sampleset<T> & spst) {
		nbsample=0;
		debut=NULL;
		courant=spst.debut;
		while (courant != NULL) {
			sample<T> * s=new sample<T>;
			s->val=courant->val;
			s->no=courant->no;
			s->lab=courant->lab;
			s->suivant=debut;
			debut=s;
			courant=courant->suivant;
			nbsample++;
			}
		}
	sampleset<T>& operator= (const sampleset<T> & spst) {
		if (this != &spst) {
			nbsample=0;
			debut=NULL;
			courant=spst.debut;
			while (courant != NULL) {
				sample<T> * s=new sample<T>;
				s->val=courant->val;
				s->no=courant->no;
				s->lab=courant->lab;
				s->suivant=debut;
				debut=s;
				courant=courant->suivant;
				nbsample++;
				}
			}
			return *this;
		}
	int nsample() {return nbsample;}
	void ajoute (const sample<T> &s) {
		sample<T> *new_s=new sample<T>;
		new_s->val=s.val;
		new_s->no=s.no;
		new_s->lab=s.lab;
		new_s->suivant=debut;
		debut=new_s;
		nbsample++;
	}
	void ajoute (T t, int n=0, BYTE l=0) {
		sample<T> * new_s=new sample<T>;
		new_s->val=t;
		new_s->no=n;
		new_s->lab=l;
		new_s->suivant=debut;
		debut=new_s;
		nbsample++;
	}
	sample<T> extract (int n) {
		courant=debut;
		while (courant != NULL && courant->no != n)
			courant=courant->suivant;
		sample<T> * new_s=new sample<T>;
		if (courant != NULL) {
			new_s->val=courant->val;
			new_s->no=courant->no;
			new_s->lab=courant->lab;
			}
		return *new_s;
		}
	sample<T> remove (int n) {
		courant=debut;
		sample<T> *precedent=NULL, *adetruire=NULL;
		while (courant != NULL && courant->no != n) {
			precedent=courant; courant=precedent->suivant;}
		sample<T> * new_s=new sample<T>;
		if (courant != NULL) {
			new_s->val=courant->val;
			new_s->no=courant->no;
			new_s->lab=courant->lab;
			if (precedent!=NULL) precedent->suivant=courant->suivant;
			else debut=courant->suivant;
			adetruire=courant;
			courant=courant->suivant;
			if (adetruire!=NULL) delete adetruire;
			nbsample--;
			}
		return *new_s;
		}
	BYTE label (int n) {
		courant=debut;
		while (courant != NULL && courant->no != n)
			courant=courant->suivant;
		if (courant != NULL) return courant->lab;
		else return 0;
		}
	void affiche () {
		cout<<" nombre d'echantillons = "<<nbsample<<"\n";
		courant=debut;
		while (courant!=NULL) {
			(*courant).affiche();
			courant=courant->suivant;
			}
		}
	void k_ppv (sampleset<T> &, int=1);
	statclasses k_means (int=2);
};

template <class T> void sampleset<T>::k_ppv (sampleset<T> & classample, int k) {
	const float alpha=100.;
//	const float alpha=1.;
//	const float alpha=0.5;
	sample<T> *cl_courante, *cl_courante2, *T_ppv=new sample<T>[k];
	T val_s, val_cl;
	float dist, dmax_ppv, maxrayon_cl, *T_dist=new float[k];
	int n_ppv, i, i_dmax, nlab, labmax, T_nlab[255];
	BYTE lab;
	if (classample.nsample()<k) k=classample.nsample();
	maxrayon_cl=0;
	cl_courante=classample.debut;
	while (cl_courante!=NULL) {
		i=cl_courante->no;
		lab=cl_courante->lab;
		val_cl=cl_courante->val;
		cl_courante2=classample.debut;
		while (cl_courante2!=NULL) {
			if (cl_courante2->no>=i && cl_courante2->lab==lab) {
				dist=val_cl.distance(cl_courante2->val);
				if (dist>maxrayon_cl) maxrayon_cl=dist;
				}
			cl_courante2=cl_courante2->suivant;
			}
		cl_courante=cl_courante->suivant;
		}
	cout<<" diametre maximal des classes = "<<maxrayon_cl<<"\n";
	courant=debut;
	while (courant!=NULL) {
		val_s=courant->val;
		cl_courante=classample.debut;
		n_ppv=0;
		dmax_ppv=0.;
		while (cl_courante!=NULL) {
			val_cl=cl_courante->val;
			dist=val_s.distance(val_cl);
			if (n_ppv<k && dist<=alpha*maxrayon_cl) {
				T_ppv[n_ppv]=(*cl_courante);
				T_dist[n_ppv]=dist;
				if (dist>dmax_ppv) {
					dmax_ppv=dist;
					i_dmax=n_ppv;
					}
				n_ppv++;
				}
			else {
				if (dist<dmax_ppv) {
					T_ppv[i_dmax]=(*cl_courante);
					T_dist[i_dmax]=dist;
					dmax_ppv=0;
					for (i=0; i<n_ppv; i++) {
						if (T_dist[i]>dmax_ppv) {
							dmax_ppv=T_dist[i];
							i_dmax=i;
							}
						}
					}
				}
			cl_courante=cl_courante->suivant;
			}
//		for (i=0; i<n_ppv; i++) {T_ppv[i].affiche();
//				         cout<<" | "<<T_dist[i]<<"\n";}
		for (i=0; i<255; i++) T_nlab[i]=0;
		if (n_ppv>0) {
			labmax=0;
			for (i=0; i<n_ppv; i++) {
				lab=(T_ppv[i]).label();
				T_nlab[(int)lab]++;
				if ((int)lab>labmax) labmax=(int)lab;
				}
			nlab=0;
			for (i=0; i<=labmax; i++) {
				if (T_nlab[i]>nlab) {
					nlab=T_nlab[i];
					lab=(BYTE)i;
					}				
				}
//			cout<<" => label majoritaire = "<<(int)lab<<"\n";
			}
		else {
			lab=0;
			cout<<" pixel inclassifiable\n";
			val_s.affiche();
			}
		courant->lab=lab;
//		if (courant->no%1000==0) 
//			cout<<" classif. echantillon "<<courant->no<<" val = "<<courant->val[0]<<" : lab = "<<(int)courant->lab<<"\n";
		courant=courant->suivant;
		}
	if (T_ppv!=NULL) {delete[] T_ppv; T_ppv=NULL;}
	if (T_dist!=NULL) {delete[] T_dist; T_dist=NULL;}
	}

template <class T> statclasses sampleset<T>::k_means (int k) {cout<<" entree dans k_means\n";
	sample<T> s0=extract(0);s0.affiche();
	int dim=(s0.val).ncanaux();cout<<" dim = "<<dim<<"\n";
	T moy(dim), var(dim), max=s0.val, min=s0.val, p;
	sample<T> *s;
	BYTE lab;
	int i,j,n;
// initialisation des centres de classes
	s=debut;
	while (s != NULL) {
		p=s->val;
		moy=moy+p;
		var=var+p*p;
		max.sup(p);
		min.inf(p);
		s=s->suivant;
	}
	moy=moy*(1.f/nbsample);
	var=var*(1.f/(nbsample-1))-moy*moy*(1.f*nbsample/(nbsample-1));
	cout<<"  moyenne des echantillons = ";moy.affiche();
	cout<<" variance des echantillons = ";var.affiche();
	cout<<" valeur min des echantillons = ";min.affiche();
	cout<<" valeur max des echantillons = ";max.affiche();
	cout<<" initialisation des "<<k<<" classes\n";
	statclasses Tclas;
	n=1;                   // att.: les numeros de classe commencent a 1
	if (k%2==1) {
		stat1class C((BYTE)n,dim);
		for (i=0; i<dim; i++) C.mean(i)=moy[i];
		Tclas.ajoute(C);
		n++;
	}
	float coef=1.f/(k/2+1);
	for (j=0; j<k/2; j++) {
		lab=(BYTE)n;
		stat1class C1(lab,dim);
		for (i=0; i<dim; i++) C1.mean(i)=min[i]+(moy[i]-min[i])*coef*(j+1);
		Tclas.ajoute(C1); n++;
		lab=(BYTE)n;
		stat1class C2(lab,dim);
		for (i=0; i<dim; i++) C2.mean(i)=max[i]-(max[i]-moy[i])*coef*(j+1);
		Tclas.ajoute(C2); n++;
	}
	Tclas.affiche();
// convergence
	bool fin;
	pixel<float> dist_cl;
	int* T_chg=new int[k+1];
	int nch, iter=0;
	fin=0;
	while (!fin) {
		for (i=0; i<=k; i++) T_chg[i]=0;
		nch=0;
		courant=debut;
		while (courant!=NULL) {             // affectation des labels
			p=courant->val;
			dist_cl=Tclas.distclasses(p);
			i=(int)dist_cl[0];
			j=(int)courant->lab;
			if (i!=j) {
				courant->lab=(BYTE)i;
				T_chg[i]++;
				T_chg[j]++;
			}
			courant=courant->suivant;
		}
		for (i=1; i<=k; i++) nch+=T_chg[i];
		cout<<" iteration "<<setw(3)<<iter;
		cout<<" : # changements = "<<setw(4)<<nch/2<<"\n";
		if (nch>0) {
			statclasses newTcl;
			for (i=1; i<=k; i++)  {       // mise a jour des centres
				lab=(BYTE)i; // des classes ayant change
				if (T_chg[i]==0) {
					stat1class cl1=Tclas.extract(lab);
					newTcl.ajoute(cl1);
					}
				else {
					stat1class cl1(lab,dim);
					pixel<float> centre(dim,0);
					s=debut;
					n=0;
					while (s != NULL) {
						if (s->lab==lab) {
							p=s->val;
							centre=centre+p;
							n++;
							}
						s=s->suivant;
						}
					centre=centre*(1.f/n);
					for (j=0; j<dim; j++) cl1.mean(j)=centre[j];
					newTcl.ajoute(cl1);
					}
				}
//			cout<<" nouvelles classes : \n"; newTcl.affiche();
			Tclas=newTcl;
			iter++;
			}
		else fin=1;
		}
		if (T_chg!=NULL) {delete[] T_chg; T_chg=NULL;}
		return Tclas;
	}

#endif