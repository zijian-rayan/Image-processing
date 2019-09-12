#ifndef _IMAREGIONS_H
#define _IMAREGIONS_H

#include <iostream>
#include <fstream>
//#include <iomanip>
//#include <string>
using namespace std;
#include <math.h>
//#include "float.h"
//#include "constantes.h"
#include "fonctions.h"
#include "imasites.h"
#include "pixel.h"
#include "liste_pixels.h"
#include "imadata.h"
#include "imabin.h"
#include "imalabels.h"
#include "imacouleur.h"

const unsigned int nmaxarreteparsommet=6000;
struct sommet {
	liste_pixels Lpix;
	double val;
	unsigned int n_arretes;
	unsigned int T_arretes[nmaxarreteparsommet];
	bool valid;
};

struct arrete {
	unsigned int pred, succ;
	double val, length;
	bool valid;
};

struct region {
	liste_pixels Lpix;
	double val;
	unsigned int n_arretes;
	unsigned int *T_arretes;
	bool valid;
};

struct regionV {
	liste_pixels Lpix;
	double *val;
	unsigned int n_arretes;
	unsigned int *T_arretes;
	bool valid;
};

class imaregions : public imasites
{	int nblayer;
	long int *Tima;
	long int valnul, nbregions;
	float s_homog;
 public:
/* constructeur par defaut */
	imaregions (int nl=1, int nc=1, int nd=1) : imasites(nl,nc) {
		Tima=NULL;
		long int i; int j;
		nblayer=nd;
		valnul=/*-1*/0; // valnul=0 si les labels de region commencent à 1
		Tima=new long int[nbpix*nblayer];
		for (i=0; i<nbpix; i++) {
			for (j=0; j<nblayer; j++) Tima[i*nblayer+j]=valnul;
			Tsites[i]=&(Tima[i*nblayer]);
			}
		nbregions=0;
		s_homog=-1;
	}

/* constructeur par conversion d'une image en int*4 en une image de regions */
	imaregions (imadata<long int> &ima, int k=0) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL;
		nblayer=1;
		Tima=new long int[nbpix*nblayer];
		valnul=/*-1*/0;
		int i,j, l0=(int)(ima.minI()-1-valnul); long int ii, n;
		nbregions=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				Tima[ii]=n=ima(i,j,k)-l0;
				if (n-valnul>nbregions) nbregions=n-valnul;
				Tsites[ii]=&(Tima[ii]);
			}
		s_homog=-1;
	}

/* constructeur par conversion d'une image en BYTE en une image de regions */
	imaregions (imadata<BYTE> &ima, int k=0) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL;
		nblayer=1;
		Tima=new long int[nbpix*nblayer];
		valnul=/*-1*/0;
		int i,j, l0=(int)(ima.minI()-1-valnul); long int ii, n;
		nbregions=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				Tima[ii]=n=(long int)(ima(i,j,k)-l0);
				if (n-valnul>nbregions) nbregions=n-valnul;
				Tsites[ii]=&(Tima[ii]);
			}
		s_homog=-1;
	}

/* constructeur par conversion d'une image en BYTE en une image de regions */
	imaregions (imadata<BYTE> &ima, char *extent) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL; Tima=new long int[nbpix*nblayer];
		nblayer=1; valnul=0/*-1*/; nbregions=0;
		int i,j,k;
		char *pgm="pgm", *ppm="ppm";
		bool trouve=0;
		int *coul_reg=new int[nbpix], coul, no_reg, ii/*, nb_reg=(int)nbregions*/;
		for (ii=0; ii<nbpix; ii++) coul_reg[ii]=-1;
		if (strcmp(extent,pgm)==0) {
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					coul=(int)ima(i,j);
					trouve=0; ii=0; no_reg=nbregions+1;
					while (!trouve && ii<nbregions && ii<nbpix) {
						if (coul_reg[ii]==coul) {trouve=1; no_reg=ii+1;}
						ii++;
					}
					if (!trouve) coul_reg[nbregions++]=coul;
					Tima[i*nbcol+j]=no_reg;
					Tsites[i*nbcol+j]=&(Tima[i*nbcol+j]);
				}
		} else {
			if (strcmp(extent,ppm)==0) {
				cout<<" conversion de l'image couleur des segments en une 'imaregions'\n";
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						coul=0;	for (k=0; k<3; k++)	coul+=(long int)(pow(256.f,k)*ima(i,j,k)); //cout<<" en ("<<i<<", "<<j<<") couleur="<<coul<<"\n";
						trouve=0; ii=0; 
						while (!trouve && ii<nbregions && ii<nbpix) {
							if (coul_reg[ii]==coul) {trouve=1; no_reg=ii+1;}
							ii++;
						}
						if (!trouve) {coul_reg[nbregions++]=coul; no_reg=nbregions;
/*							coul_reg[nbregions]=coul;
							no_reg=nbregions++;
							cout<<" nouvelle region en ("<<i<<", "<<j<<") couleur="<<coul<<" => #regions ="<<nbregions<<"\n";
							char aa; cin>>aa;*/
						}
						Tima[i*nbcol+j]=no_reg;
						Tsites[i*nbcol+j]=&(Tima[i*nbcol+j]);
					}
			} else
				cout<<" format de donnees non prevu\n";
		}
		if (coul_reg!=NULL) {delete[] coul_reg; coul_reg=NULL;}
		cout<<" nb regions = "<<nbregions<<"\n";
		s_homog=-1;
	}

/* constructeur par etiquettage de regions connexes dans une image binaire */
	imaregions (imabin&, bool=0);

/* constructeur par etiquettage de regions connexes dans une image des labels */
	imaregions (imalabels&, int=0);

/* superposition de 2 segmentations */
	imaregions (imaregions&, imaregions&);

/* verification region 1 composante connexe */
	bool verif_1cc_per_region() {cout<<" nbregions = "<<nbregions<<"\n";
		int i,j,k,ncc,n,nmax,max,ii,jj;
		bool ok=true;
		imabin imab(nblig,nbcol);
		for (k=0; k<=nbregions; k++) {
			imab.mise_a_zero();
			for (i=0; i<nblig; ++i)
				for (j=0; j<nbcol; ++j) if ((*this)(i,j)==k) imab(i,j)=1;
			if (imab.norm()>0) {
				imadata<int> imab_cc=imab.composantes_connexes(ncc,8);
				if (ncc!=1) {
					cout<<" region "<<k<<" a "<<ncc<<" composantes connexes ??? taille respective = \n"; 
					int *tabnpix=new int[ncc+1], *tabvois=new int[nbregions+1]; for (n=0; n<=ncc; n++) tabnpix[n]=0; 
					for (i=0; i<nblig; ++i) for (j=0; j<nbcol; ++j) if (imab(i,j)) tabnpix[imab_cc(i,j)]++;
					nmax=0;
					for (n=1; n<=ncc; n++) {if (tabnpix[nmax]<tabnpix[n]) nmax=n; cout<<" "<<tabnpix[n]; } cout<<"\n";
					for (i=0; i<nblig; ++i) for (j=0; j<nbcol; ++j) if (imab(i,j) && imab_cc(i,j)!=nmax) {
						if (tabnpix[imab_cc(i,j)]==1) {
							for (n=0; n<=nbregions; n++) tabvois[n]=0;
							for (ii=maxi(i-1,0); ii<mini(i+1,nblig); ++ii) 
								for (jj=maxi(j-1,0); jj<mini(j+1,nbcol); ++jj) 
									if ((*this)(ii,jj)<=nbregions) tabvois[(*this)(ii,jj)]++; else cout<<" (*this)(ii,jj) = "<<(*this)(ii,jj)<<"???\n";
							max=0; for (n=0; n<=nbregions; n++) if (n!=(*this)(i,j) && tabvois[max]<tabvois[n]) max=n;
							(*this)(i,j)=max; cout<<" pixel isole remplace par valeur CC "<<max<<"\n";
							tabnpix[imab_cc(i,j)]=0;
						} else {
							if (tabnpix[imab_cc(i,j)]>0 && tabnpix[nmax]>10 && tabnpix[imab_cc(i,j)]<5) {
								n=imab_cc(i,j);
								imabin imabtmp(nblig,nbcol); for (ii=0; ii<nblig; ++ii) for (jj=0; jj<nbcol; ++jj) imabtmp(ii,jj)=(imab_cc(ii,jj)==n?1:0);
								imabin imabtmp2=imabtmp.dilate(eltstruct(3,3))-imabtmp;
								for (n=0; n<=nbregions; n++) tabvois[n]=0;
								for (ii=0; ii<nblig; ++ii) for (jj=0; jj<nbcol; ++jj) if (imabtmp2(ii,jj)) tabvois[(*this)(ii,jj)]++;
								max=0; for (n=0; n<=nbregions; n++) if (tabvois[max]<tabvois[n]) max=n;
								for (ii=0; ii<nblig; ++ii) for (jj=0; jj<nbcol; ++jj) if (imabtmp(ii,jj)) (*this)(i,j)=max; 
								cout<<" petite CC isolee remplacee par valeur CC "<<max<<"\n"; tabnpix[imab_cc(i,j)]=0;
							} 
						}
					} //cout<<" ************\n";
					for (n=1; n<=ncc; n++) 
						if (tabnpix[n]>0 && n!=nmax) {
								cout<<" on duplique CC "<<n<<"....\n"; nbregions++;
								for (i=0; i<nblig; ++i) for (j=0; j<nbcol; ++j) if (imab_cc(i,j)==n) (*this)(i,j)=nbregions; 
							}
					if (tabnpix!=NULL) delete[] tabnpix; if (tabvois!=NULL) delete[] tabvois; //cout<<" !!!!!!!!!!!\n"; 
					char aa; cin>>aa;
					ok=false;
				}
			}
		}
		return ok;
	}

/* critere d'homogeneite d'une region */
	bool verif_homog (float, float);

/* constructeur par croissance de regions */
	imaregions (imadata<float> &, float, char*, bool=0);

/* constructeur par quadtree */
	imaregions (imadata<float> &, float, unsigned short int=2, bool=1);

/* constructeur par fusion de regions dans un graphe, arrêt sur critère d'homogénéité */
	imaregions (imadata<float> &, double);

/* constructeur par fusion de regions dans un graphe, arrêt sur nombre de régions */
	imaregions (imadata<float> &, unsigned int);

/* constructeur par segmentation de Mumford & Sha éventuellement à partir d'1 image de superpixels */
	unsigned long int comptearrete (imadata<int> &, unsigned long int);
	void cree_graph_init (imadata<int> &, imadata<float> &, regionV*, arrete*, unsigned long int &, unsigned long int &);
	imaregions (imadata<float> &, imadata<int> &, unsigned int, char*);
	imaregions (imadata<float> &, unsigned int, char*);

/* constructeur par ligne de partage des eaux */
	imaregions (imadata<float> &, int=4, bool=0);

/* Waterpixels */
	imaregions (imadata<float> &ima, double param_k, long int &nregions, int iconnex=4) : imasites(ima) {
		Tima=NULL; nblayer=1; cout<<" entree dans constructeur Waterpixels\n";
		Tima=new long int[nbpix*nblayer];
		long int n;
		int i,j; 
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {n=i*nbcol+j; Tima[n]=-1; Tsites[n]=&(Tima[n]);}
		valnul=-1;
		const float eps=(float)1.e-6; cout<<" entree dans constructeur waterpixels\n";
/*		imadonneesMS.statbasic(1);
		if (ifiltr) {eltstruct ES(3,3); // ou (sigma^2)/16
			imadonnees=imadonnees.ouverture(ES).fermeture(ES); }*/
		imadata<float> ima_dirGrad, ima_modGrad;
		ima_modGrad=ima.gradient(ima_dirGrad,"Sobel"); ima_modGrad.imaunsignedchar(1).sauve_ImaPGM("ima_modGrad.pgm");
/*		if (iLPE) {
			imadata<float> ima_modGradbis=ima_modGrad+(-2.f);
			ima_modGradbis=ima_modGrad.reconst_geod (ima_modGradbis,eltstruct(3,3));
			imaregions ima_LPE(ima_modGradbis); // LPE
			strncpy(nomfile_out,nomfic_out,lg_nomfic+1); strncpy(nomfile_out+lg_nomfic,"_LPE.pgm",9);
			ima_LPE.conv2imUI().sauve_ImaPGM2(nomfile_out);
		}*/
		int nstep=(int)pow((float)nregions,0.5f), istep=mini(nblig,nbcol)/nstep;
		imabin imaGrid(nblig,nbcol), imaSeed(nblig,nbcol);
		for (i=0; i<nblig; i+=istep) for (j=0; j<nbcol; j++) imaGrid(i,j)=1;
		for (j=0; j<nbcol; j+=istep) for (i=0; i<nblig; i++) imaGrid(i,j)=1;
		float Xmax=(float)ima_modGrad.maxI(), rho=istep/3.f, xmin, xx; // rho est la marge pour ne pas considérer de minima trop proches des bords
		int ii, jj, ncc, nmax, imax;
		imadata<float> ima_cell(istep,istep);
		for (i=0; i<nblig; i+=istep)
			for (j=0; j<nbcol; j+=istep) {
				xmin=FLT_MAX; ima_cell.mise_a_zero(); ima_cell=ima_cell+Xmax;
				for (ii=i+(int)rho; ii<i+istep-(int)rho; ii++)
					for (jj=j+(int)rho; jj<j+istep-(int)rho; jj++)
						if (ii<nblig && jj<nbcol) {
							xx=ima_modGrad(ii,jj);
							ima_cell(ii-i,jj-j)=xx; // ima_cell est une mini image des valeurs du gradient sur la cellule considérée
							if (xx<xmin) xmin=xx;
						}
				imabin ima_cellb(ima_cell,xmin+eps); ima_cellb=ima_cellb.negatif();
				imadata<int> ima_cellb_cc=ima_cellb.composantes_connexes (ncc,iconnex,0);
				if (ncc>1) {
					int *tab=new int[ncc]; for (ii=0; ii<ncc; ii++) tab[ii]=0;
					for (ii=(int)rho; ii<istep-(int)rho; ii++)
						for (jj=(int)rho; jj<istep-(int)rho; jj++)
							if (ima_cellb_cc(ii,jj)>0) tab[ima_cellb_cc(ii,jj)-1]++;
					nmax=0; imax=0;
					for (ii=0; ii<ncc; ii++) if (tab[ii]>nmax) {nmax=tab[ii]; imax=ii+1;}
					if (tab!=NULL) delete[] tab;
				} 
				else imax=1;
				if (imax>0)
					for (ii=i+(int)rho; ii<i+istep-(int)rho; ii++)
						for (jj=j+(int)rho; jj<j+istep-(int)rho; jj++)
							if (ima_cellb_cc(ii-i,jj-j)==imax && ii<nblig && jj<nbcol) imaSeed(ii,jj)=1;
			}
		imaSeed.imaunsignedchar(1).sauve_ImaPGM("ima_Seed.pgm");
		imadata<float> ima_d=imaSeed.Tr_dist(3); ima_d.statbasic(1); 
		cout<<" sauvegarde de ``ima_DstS.pgm''\n"; ima_d.imaunsignedchar(1).sauve_ImaPGM("ima_DstS.pgm");
// image du gradient `regularisé'
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) ima_modGrad(i,j)=(float)(ima_modGrad(i,j)*1.0f+ima_d(i,j)*(2*param_k/istep));
		ima_modGrad.statbasic(1);
		cout<<" sauvegarde de ``ima_GrdRg.pgm''\n"; ima_modGrad.imaunsignedchar(1).sauve_ImaPGM("ima_GrdRg.pgm");
// un seul minimum par superpixel
		float max_modGrad=(float)ima_modGrad.maxI(), min_modGrad=(float)ima_modGrad.minI(); cout<<" max ima_modGrad = "<<max_modGrad<<"\n";
		imadata<float> ima_modGradinv(nblig,nbcol), imaSeed_fl(nblig,nbcol);
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) {
				ima_modGradinv(i,j)=(max_modGrad-ima_modGrad(i,j))/(max_modGrad-min_modGrad)*255; 
				imaSeed_fl(i,j)=imaSeed(i,j)*255.f/*max_modGrad*/;}
		eltstruct ES4(3,3); ES4(0,0)=ES4(2,2)=ES4(0,2)=ES4(2,0)=0;
		ima_modGradinv=ima_modGradinv.reconst_geod(imaSeed_fl,ES4);
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) ima_modGrad(i,j)=max_modGrad-ima_modGradinv(i,j);
		ima_modGrad.statbasic(1);
		cout<<" sauvegarde de ``ima_GrdRgRc.pgm''\n"; ima_modGrad.imaunsignedchar(1).sauve_ImaPGM("ima_GrdRgRc.pgm");
// LPE sur image du gradient `regularisé' à un seul minimum par superpixel
		imaregions ima_WaterPix(ima_modGrad,iconnex,1); // LPE
		ima_WaterPix.sauve_Ima("ima_WatPix.dat"); ima_WaterPix.conv2imUI().sauve_ImaPGM2("ima_WatPix.pgm");
		cout<<" image des waterpixels sauvee sous ``ima_WatPix.pgm''\n";
//		ima_WaterPix.suppression_trop_petites_regions();
//		int nSPix=ima_WaterPix.nregions(); cout<<nSPix<<" superpixels\n";
//		imadata<int> ima_WaterPixI=ima_WaterPix.conv2imUI();
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) Tima[i*nbcol+j]=ima_WaterPix(i,j);
	}

/* superpixels SLIC */
	imaregions (imadata<float> &ima, long int &nregions, double fact) : imasites(ima) {cout<<" entree dans constructeur superpixels SLIC\n";
		Tima=NULL; nblayer=1;
		Tima=new long int[nbpix*nblayer];
		long int n;
		int it=0,i,j,k,i0,j0,ii,jj; 
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {n=i*nbcol+j; Tima[n]=-1; Tsites[n]=&(Tima[n]);}
		valnul=-1;
		imadata<float> ima_di, imagrd=ima.gradient (ima_di);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) ima_di(i,j)=FLT_MAX; 
		const int dim=ima.ncanaux(), itmax=100; 
		const float rr2=(float)nbpix/nregions, rr=pow(rr2,0.5f), drr=rr*1.25f; cout<<" rr2="<<rr2<<" rr="<<rr<<" drr="<<drr<<"\n";
		double **T_Ck=new double*[nregions], **T_Ck_new=new double*[nregions];
		double dr, ds, fact2=pow(fact,2.), dt;
		for (k=0; k<nregions; k++) {T_Ck[k]=new double[dim+3]; T_Ck_new[k]=new double[dim+3];}
		float xi=rr/2.f, xj=rr/2.f; 
		k=0;
		do {//cout<<" k="<<k; cout<<" : i="<<xi<<", j="<<xj;
			i0=around(xi); j0=around(xj);
			for (ii=maxi(0,around(xi)-1); ii<=mini(around(xi)+1,nblig-1); ii++)
				for (jj=maxi(0,around(xj)-1); jj<=mini(around(xj)+1,nbcol-1); jj++) if (imagrd(ii,jj)<imagrd(i0,j0)) {i0=ii; j0=jj;}
			T_Ck[k][dim]=i0; T_Ck[k][dim+1]=j0; //cout<<" -> i0="<<i0<<", j0="<<j0<<"\n";
			for (ii=0; ii<dim; ii++) T_Ck[k][ii]=ima(i0,j0,ii); 
//			for (ii=0; ii<dim+2; ii++) cout<<setw(8)<<setprecision(4)<<T_Ck[k][ii]<<","; cout<<"\n";
			xj+=rr; if (around(xj)>=nbcol) {xj-=nbcol; xi+=rr;} 
			if (around(xi)>=nblig) {cout<<" depassement selon lig et col !!!\n"; //nregions=k+1; 
				xi=(float)(rand()%nblig); xj=(float)(rand()%nbcol);}
		} while (++k<nregions); 
		{imacouleur<BYTE> ima_init(ima.imaunsignedchar());
		for (k=0; k<nregions; k++) {
			i0=(int)T_Ck[k][dim]; j0=(int)T_Ck[k][dim+1];
			for (ii=maxi(0,i0-1); ii<mini(i0+1,nblig-1); ii++)
				for (jj=maxi(0,j0-1); jj<=mini(j0+1,nbcol-1); jj++) {ima_init(ii,jj,0)=255; ima_init(ii,jj,1)=ima_init(ii,jj,2)=0;}
		} ima_init.sauve_ImaPGM("ima_init.ppm"); 
		}
//		cout<<" fin initialisations\n";
		bool fini=0;
		const float err_r=5.f, err_s=2.f, threshold=dim*pow(err_r,2)+2*pow(err_s,2); cout<<" seuil sur les erreurs = "<<threshold<<"\n";
		do {cout<<" iteration "<<it<<", ";
			for (k=0; k<nregions; k++) {
				i0=(int)T_Ck[k][dim]; j0=(int)T_Ck[k][dim+1];
				for (ii=maxi(0,i0-(int)drr); ii<=mini(i0+(int)drr,nblig-1); ii++)
					for (jj=maxi(0,j0-(int)drr); jj<=mini(j0+(int)drr,nbcol-1); jj++) {
						dr=0.; ds=pow(i0-ii,2.)+pow(j0-jj,2.);
						for (i=0; i<dim; i++) dr+=pow(T_Ck[k][i]-ima(ii,jj,i),2.);
						dt=dr+ds/rr2*fact2;
						if (dt<ima_di(ii,jj)) {ima_di(ii,jj)=(float)dt; (*this)(ii,jj)=k+1;}
					}
			} //conv2imBYTE().sauve_ImaPGM("imaSLIC.pgm"); {char aa; cin>>aa;}
			for (k=0; k<nregions; k++) 
				for (ii=0; ii<dim+3; ii++) T_Ck_new[k][ii]=0.;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					k=(*this)(i,j)-1;
					if (k>=0 && k<nregions) {
						for (ii=0; ii<dim; ii++) T_Ck_new[k][ii]+=ima(i,j,ii);
						T_Ck_new[k][dim]+=i; T_Ck_new[k][dim+1]+=j; T_Ck_new[k][dim+2]+=1;
					} else {cout<<" pb dans label superpixel en ("<<i<<","<<j<<") : label = "<<k<<"\n"; char aa; cin>>aa;}
				}
			for (k=0; k<nregions; k++) 
				if (T_Ck_new[k][dim+2]>0) for (ii=0; ii<dim+2; ii++) T_Ck_new[k][ii]/=T_Ck_new[k][dim+2];
//			for (k=0; k<nregions; k++) {cout<<setw(4)<<k<<":"; 
//				for (ii=0; ii<dim+3; ii++) cout<<setw(8)<<setprecision(4)<<T_Ck_new[k][ii]<<","; cout<<"\n"; }
			dt=0.;
			for (k=0; k<nregions; k++) 
				for (ii=0; ii<dim+2; ii++) dt+=pow(T_Ck_new[k][ii]-T_Ck[k][ii],2.);
			cout<<" erreur = "<<dt<<"\n";
			if (dt<threshold || ++it>itmax) fini=1;
			else {
				for (k=0; k<nregions; k++) 
					if (T_Ck_new[k][dim+2]>0) for (ii=0; ii<dim+3; ii++) T_Ck[k][ii]=T_Ck_new[k][ii];
			}
		} while (!fini); conv2imBYTE().sauve_ImaPGM("imaSLIC.pgm"); //{char aa; cin>>aa;}
		for (k=0; k<nregions; k++) {if (T_Ck[k]!=NULL) delete[] T_Ck[k]; if (T_Ck_new[k]!=NULL) delete[] T_Ck_new[k]; }
		if (T_Ck!=NULL) delete[] T_Ck; if (T_Ck_new!=NULL) delete[] T_Ck_new;
		nbregions=nregions;
		s_homog=-1; cout<<" fin constructeur superpixels SLIC\n";
	}

/* constructeur de recopie */
	imaregions (const imaregions &ima) : imasites(ima) {
		Tima=NULL;
		nblayer=ima.nblayer;
		valnul=ima.valnul;
		Tima=new long int[nbpix*nblayer];
		long int *adval;
		for (long int i=0; i<nbpix; i++) {
			adval=(long int*)ima.Tsites[i];
			for (int j=0; j<nblayer; j++) Tima[i*nblayer+j]=*(adval+j);
			Tsites[i]=&(Tima[i*nblayer]);
		}
		nbregions=ima.nbregions;
	}

/* destructeur */
	~imaregions () {
		if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
	}

/* operateur d'affectation */
	imaregions& operator=(const imaregions &ima) {
		if (this != &ima) {
			imasites *ad1, *ad2;
			ad1=this;
			ad2=(imasites *) &ima;
			*ad1=*ad2;
			if (Tima!=NULL) delete[] Tima;
			nblayer=ima.nblayer;
			valnul=ima.valnul;
			Tima=new long int[nbpix*nblayer];
			long int *adval;
			for (long int i=0; i<nbpix; i++) {
				adval=(long int*)ima.Tsites[i];
				for (int j=0; j<nblayer; j++) Tima[i*nblayer+j]=*(adval+j);
				Tsites[i]=&(Tima[i*nblayer]);
			}
			nbregions=ima.nbregions;
		}
		return *this;
	}

/* # de regions */
	long int nregions() {return nbregions;}

/* # de couches */
	int nlayer() {return nblayer;}

/* operateur d'acces a la valeur (i,j,k) */
	long int& operator () (int i, int j, int k=0) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol || k<0 || k>=nblayer) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<","<<k<<")\n";
			if (i<0) i=0; if (i>=nblig) i=nblig-1;
			if (j<0) j=0; if (j>=nbcol) j=nbcol-1;
			if (k<0) k=0; if (k>=nblayer) k=nblayer-1;
		}
		long int *adval=(long int*)Tsites[i*nbcol+j];
		adval=adval+k;
		return *adval;
	}

/* elimination petites regions */
	void suppression_trop_petites_regions (int=5);
	void suppression_plus_petites_regions (int=255);

/* reconnexion des regions non connexes */
	void reconnecte_reg (int npixbig=100, int npixmin=5, int iconnex=8) {cout<<" entree dans reconnecte_reg\n";
		imabin imabk(nblig,nbcol), imab(nblig,nbcol); imab.mise_a_zero();
		int i,j,k,ncc,ii,jj,n,n2; 
		int *tabreg=new int[nbregions+1];
		imadata<int> imacck;
		for (k=0; k<nbregions; k++) {
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) if ((*this)(i,j)==k) imabk(i,j)=1; else imabk(i,j)=0;
			imacck=imabk.composantes_connexes(ncc,iconnex);
			if (ncc>1) {
//				cout<<" region "<<k<<" "<<ncc<<"-composantes connexes\n";
				int *tabcc=new int[ncc+1]; for (ii=0; ii<=ncc; ii++) tabcc[ii]=0;
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) if (imabk(i,j)) tabcc[imacck(i,j)]++;
				jj=-1; n=0;
				for (ii=1; ii<=ncc; ii++) 
					if (tabcc[ii]>n) {n=tabcc[ii]; jj=ii;}
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) if (imabk(i,j) && imacck(i,j)!=jj) imab(i,j)=1;
				if (tabcc!=NULL) delete[] tabcc;
			}
		}
//		for (i=0; i<nblig; i++)
//			for (j=0; j<nbcol; j++) if ((*this)(i,j)==0) imabk(i,j)=1; else imabk(i,j)=0;
		imacck=imab.composantes_connexes(ncc,iconnex);
		cout<<ncc<<" regions non labelisees\n";
		if (nbregions<256) conv2imBYTE().sauve_ImaPGM("ima_0_reg_not_connexes.pgm");
		else if (nbregions<65536) conv2imUI().sauve_ImaPGM2("ima_0_reg_not_connexes.pgm");
		int *tabcc=new int[ncc+1]; for (ii=0; ii<=ncc; ii++) tabcc[ii]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) if (imab(i,j)) tabcc[imacck(i,j)]++;
		for (k=1; k<=ncc; k++)																										 // cas des regions de 1 seul pixel
			if (tabcc[k]==1) { cout<<" composante "<<k<<" a seulement 1 pixel";
				for (ii=0; ii<=nbregions; ii++) tabreg[ii]=0;
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) 
						if (imacck(i,j)==k) {cout<<" en ("<<i<<","<<j<<")";
							for (ii=maxi(0,i-1); ii<=mini(nblig-1,i+1); ii++)
								for (jj=maxi(0,j-1); jj<=mini(nbcol-1,j+1); jj++) 
									if (!imab(ii,jj) && (*this)(ii,jj)>0 && (*this)(ii,jj)<=nbregions) tabreg[(*this)(ii,jj)]++;
							jj=-1; n=n2=0;
							for (ii=1; ii<=nbregions; ii++) if (tabreg[ii]>n) {n=tabreg[ii]; jj=ii;}
							for (ii=1; ii<=nbregions; ii++) if (ii!=jj && tabreg[ii]>n2) {n2=tabreg[ii];}
							if (n2==0) cout<<" exclusivement entoure de label "<<jj<<"\n";
							else cout<<" entoure par "<<n<<" (>"<<n2<<") pixels de label "<<jj<<"\n";
							if (jj>0 && jj<=nbregions) {(*this)(i,j)=jj; tabcc[k]=0; imab(i,j)=0;}
							i=nblig; j=nbcol;
						}
			}
		if (npixmin>1) {																										 // cas des petites regions de quelques pixels
			for (k=1; k<=ncc; k++)
				if (tabcc[k]>0 && tabcc[k]<=npixmin) { cout<<" composante "<<k<<" a seulement "<<tabcc[k]<<" pixel(s) en";
					for (ii=0; ii<=nbregions; ii++) tabreg[ii]=0;
					for (i=0; i<nblig; i++)
						for (j=0; j<nbcol; j++) 
							if (imacck(i,j)==k) {cout<<" ("<<i<<","<<j<<")";
								for (ii=maxi(0,i-1); ii<=mini(nblig-1,i+1); ii++)
									for (jj=maxi(0,j-1); jj<=mini(nbcol-1,j+1); jj++) 
										if (!imab(ii,jj) && (*this)(ii,jj)>0 && (*this)(ii,jj)<=nbregions) tabreg[(*this)(ii,jj)]++;
							}
					jj=-1; n=n2=0;
					for (ii=1; ii<=nbregions; ii++) if (tabreg[ii]>n) {n=tabreg[ii]; jj=ii;}
					for (ii=1; ii<=nbregions; ii++) if (ii!=jj && tabreg[ii]>n2) {n2=tabreg[ii];}
					if (n2==0) cout<<" exclusivement entoure de label "<<jj<<"\n";
					else cout<<" entoure par "<<n<<" (>"<<n2<<") pixels de label "<<jj<<"\n";
					if (jj>0 && jj<=nbregions) {
						for (i=0; i<nblig; i++)
							for (j=0; j<nbcol; j++) if (imacck(i,j)==k) {(*this)(i,j)=jj; imab(i,j)=0;}
						tabcc[k]=0; 
					}
				}
		}
		n=nbregions;
		for (k=1; k<=ncc; k++)																															// cas des grandes regions
			if (tabcc[k]>=npixbig) { cout<<" composante "<<k<<" plus de "<<tabcc[k]<<" pixel(s)";
				for (ii=0; ii<=nbregions; ii++) tabreg[ii]=0;
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) 
						if (imacck(i,j)==k)
							if ((*this)(i,j)>0 && (*this)(i,j)<=nbregions) tabreg[(*this)(i,j)]++;
				for (ii=1; ii<=nbregions; ii++)
					if (tabreg[ii]>npixbig) {
						n++; cout<<" creation d'une nouvelle region de label "<<n<<" pour fractionner region "<<ii<<"\n";
						for (i=0; i<nblig; i++)
							for (j=0; j<nbcol; j++) if (imacck(i,j)==k && (*this)(i,j)==ii) {(*this)(i,j)=n; imab(i,j)=0;}
						tabreg[ii]=0;
					}
			}
		nbregions=n;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) if (imab(i,j)) (*this)(i,j)=0;

		if (tabcc!=NULL) delete[] tabcc;
		if (nbregions<256) conv2imBYTE().sauve_ImaPGM("ima_0_reg_1_pix.pgm");
		else if (nbregions<65536) conv2imUI().sauve_ImaPGM2("ima_0_reg_1_pix.pgm");
		if (tabreg!=NULL) delete[] tabreg;
		cout<<" fin reconnecte_reg\n";
	}

/* conversion en imadata<BYTE> */
	imadata<BYTE> conv2imBYTE (int k=0, bool etal=1) {
		int i,j, nreg=0;
		imadata<BYTE> ima_u(nblig,nbcol);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if ((*this)(i,j,k)>nreg) nreg=(*this)(i,j,k);
		BYTE fact=255;
		if (nreg>0) fact/=nreg;
		if (etal==0) fact=1;
		cout<<" facteur "<<(int)fact<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if ((*this)(i,j,k)>=0) ima_u(i,j)=(BYTE)(*this)(i,j,k)*fact;
				else ima_u(i,j)=0;
		return ima_u;
	}

/* conversion en imadata<unsigned int> */
	imadata<unsigned short int> conv2imUI (int k=0) {cout<<" entree dans conv2imUI\n";
		int i,j, nreg=0;
		imadata<unsigned short int> ima_u(nblig,nbcol);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if ((*this)(i,j,k)>nreg) nreg=(*this)(i,j,k);
//		unsigned int fact=65535;
//		if (nreg>0) fact/=nreg;
//		cout<<" facteur "<<(int)fact<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
//				if ((*this)(i,j,k)>=0) ima_u(i,j)=(unsigned int)(*this)(i,j,k)*fact;
				if ((*this)(i,j,k)>=0) ima_u(i,j)=(unsigned short int)(*this)(i,j,k);
				else ima_u(i,j)=0;
		ima_u.statbasic(1);
		return ima_u;
	}

/* affichage */
	void affiche () const;

/* ecriture dans un fichier de sortie */
	void sauve_Ima(string="imaRegsauvee.dat", int=0) const;	
	void sauve_Ima (char* ="imaRegsauvee.dat", int=0) const;
	void sauve_ImaUc (string="imaReg1dsauvee.dat", int=0) const;
	void sauve_ImaUc (char* ="imaReg1dsauvee.dat", int=0) const;
};

#endif