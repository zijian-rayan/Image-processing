#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include "pixel.h"
#include "sampleset.h"

void sampleset::ajoute (const sample & s) {
	sample * new_s;
	new_s=new sample(s);
	new_s->suivant=debut;
	debut=new_s;
	}

void sampleset::affiche () {
	cout<<" nombre d'echantillons = "<<nbsample<<"\n";
	courant=debut;
	while (courant!=NULL) {
		(*courant).affiche();
		cout<<" label = "<<courant->lab<<"\n";
		courant=courant->suivant;
		}
	}

	void k_ppv (sampleset &, int=1);
void sampleset::k_ppv (sampleset & classample, int k) {
	sample * cl_courante;
	T val_s, val_cl;
	float dist, dmax_ppv;
	int n_ppv, i, i_dmax, nlab;
	sample T_ppv[k];
	float T_dist[k];
	unsigned char lab, labmax;
	int T_nlab[255];
	if (classample.nsample()<k) k=classample.nsample();
	courant=debut;
	while (courant!=NULL) {
		val_s=courant->val;
		cl_courante=classample.debut;
		n_ppv=0;
		dmax_ppv=0.;
		while (cl_courante!=NULL) {
			val_cl=cl_courante->val;
			dist=val_s.distance(val_cl);
			if (n_ppv<k) {
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
					for (i=0; i<k; i++) {
						if (T_dist[i]>dmax_ppv) {
							dmax_ppv=T_dist[i];
							i_dmax=i;
							}
						}
					}
				}
			}
		for (i=0; i<255; i++) T_nlab[i]=0;
		labmax=0;
		for (i=0; i<k; i++) {
			lab=(T_ppv[k]).lab;
			T_nlab[(int)lab]++;
			if (lab>labmax) labmax=lab;
			}
		nlab=0
		for (i=0; i<(int)labmax; i++)
			if (T_nlab[i]>nlab) {
				nlab=T_nlab[i];
				lab=(unsigned char)i;
				}
		courant->lab=lab;
		courant=courant->suivant;
		}
	}

