#ifndef _IMABR_H
#define _IMABR_H

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <math.h>
//#include "constantes.h"
//#include "def.h"
//#include "fonctions.h"
//#include "pixel.h"
//#include "imasites.h"
#include "imadata.h"
//#include "statclass.h"
//#include "imalabels.h"
#include "matrices.h"

class imalabels;

template <class T> class imaBR : public imadata<T>
{	unsigned int nclas; // # de classes tel que les indices de classes vont de 0 à #classes-1
	double **T_A; // matrice decrivant la composition des pixels en termes de classes : #pixelsBR x #classes
	double **T_y; // matrice decrivant les caractéristiques des classes : #classes x #composantes (temporelles, spectrales, ...)
	bool A_valid, y_valid;
 public:
/* constructeur par defaut */
	imaBR (int nl=1, int nc=1, int nd=1, int ncl=1) : imadata<T>(nl,nc,nd), nclas(ncl) {
		int np=nl*nc, i;
		unsigned int k;
		T_A=new double*[np];
		for (i=0; i<np; i++) {
			T_A[i]=new double[nclas];
			T_A[i][0]=1;
			for (k=1; k<nclas; k++) T_A[i][k]=0;
		}
		A_valid=0;
		T_y=NULL; y_valid=0;
	}
/* constructeur a partir d'une image de donnees HR sous hypothèse de mélange linéaire */
	imaBR (imadata<T> imadon, unsigned int rlig=1, unsigned int rcol=1, int ncl=1) : 
		   imadata<T>(imadon.nblig/rlig,imadon.nbcol/rcol,imadon.nbcanaux), nclas(ncl) {
		int np=nl*nc,i,j,k,ii,jj,iimin,iimax,jjmin,jjmax,n;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				val=0;
				iimin=maxi(i*rlig,0); iimax=mini((i+1)*rlig,imadon.nlig());
				jjmin=maxi(j*rcol,0); jjmax=mini((j+1)*rcol,imadon.ncol());
				n=(iimax-iimin)*(jjmax-jjmin);
				for (k=0; k<nbcanaux; k++) {
					if (n>0) {
						for (ii=iimin; ii<iimax; ii++)
							for (jj=jjmin; jj<jjmax; jj++) val+=imadon(ii,jj,k);
						val/=n;
					}
					(*this)(i,j,k)=val;
				}
			}
		T_A=new double*[np];
		for (i=0; i<np; i++) {
			T_A[i]=new double[nclas];
			T_A[i][0]=1;
			for (k=1; k<nclas; k++) T_A[i][k]=0;
		}
		A_valid=0;
		T_y=NULL; y_valid=0;
	}
/* destructeur */
	~imaBR () {
		cout<<" destructeur imaBR\n";
		int np=nbcol*nblig, i;
		unsigned int k;
		if (T_A!=NULL) {
			for (i=0; i<np; i++) 
				if (T_A[i]!=NULL) {delete[] T_A[i]; T_A[i]=NULL;}
			delete[] T_A; T_A=NULL;
		}
		if (T_y!=NULL) {
			for (k=0; k<nclas; k++) 
				if (T_y[k]!=NULL) {delete[] T_y[k]; T_y[k]=NULL;}
			delete[] T_y; T_y=NULL;
		}
	}
/* constructeur de recopie */
	imaBR (const imaBR &ima) : imadata<T>(ima), nclas(ima.nclas) {
		int np=nblig*nbcol, i;
		unsigned int k;
		T_A=NULL; T_y=NULL;
		if (ima.T_A!=NULL) {
			T_A=new double*[np];
			for (i=0; i<np; i++) {
				T_A[i]=new double[nclas];
				for (k=0; k<nclas; k++) T_A[i][k]=ima.T_A[i][k];
			}
		}
		if (ima.T_y!=NULL) {
			T_y=new double*[nclas];
			for (k=0; k<nclas; k++) {
				T_y[k]=new double[nbcanaux];
				for (i=0; i<nbcanaux; i++) T_y[k][i]=ima.T_y[k][i];
			}
		}
		A_valid=ima.A_valid;
		y_valid=ima.y_valid;
		nclas=ima.nclas;
	}
/* operateur d'affectation */
	imaBR& operator= (const imaBR &ima) {
		if (this != &ima) {
			imadata<T> *ad1, *ad2;
			ad1=this;
			ad2=(imadata<T>*)&ima;
			*ad1=*ad2;
			int np=nbcol*nblig, i;
			if (T_A!=NULL) {
				for (i=0; i<np; i++) 
					if (T_A[i]!=NULL) {delete[] T_A[i]; T_A[i]=NULL;}
				delete[] T_A; T_A=NULL;
			}
			if (T_y!=NULL) {
				for (i=0; i<nclas; i++) 
					if (T_y[i]!=NULL) {delete[] T_y[i]; T_y[i]=NULL;}
				delete[] T_y; T_y=NULL;
			} 
			if (ima.T_A!=NULL) {
				T_A=new double*[np];
				for (i=0; i<np; i++) {
					T_A[i]=new double[nclas];
					for (k=0; k<nclas; k++) T_A[i][k]=ima.T_A[i][k];
				}
			}
			if (ima.T_y!=NULL) {
				T_y=new double*[nclas];
				for (i=0; i<nclas; i++) {
					T_y[i]=new double[nbcanaux];
					for (k=0; k<nbcanaux; k++) T_y[i][k]=ima.T_y[i][k];
				}
			}
			A_valid=ima.A_valid;
			y_valid=ima.y_valid;
			nclas=ima.nclas;
		}
		return *this;
	}

/* operateur d'acces au pourcentage de la classe k dans la composition d'un pixel (i,j) */
	double& pct_cl (int, int, int);

/* création de l'image des classes */
	imadata<float> creeimclas ();

/* nombre de classes */
	int nclasses () {return nclas;}

/* calcul de la composition des pixels à partir d'une image des labels HR */
	void labelHR_to_comppixBR (imalabels);

/* calcul de la composition des pixels à partir d'une image des labels HR décalée par rapport à l'image BR */
//	void labelHRd_to_comppixBR (imalabels,int=0,int=0);

/* desaggregation ou estimation des caracteristiques des classes connaissant la composition des pixels */
	void desaggregation (matrice2D<bool> &V=matrice2D<bool>(0,0));

/* erreur quadratique par pixel */
	matrice2D<double> err2_ppix(bool=0);

/* estimation de la composition des pixels connaissant les caracteristiques des classes */
	void decomposition();
};

/* operateur d'acces au pourcentage de la classe k dans la composition d'un pixel (i,j) */
template <class T> double& imaBR<T>::pct_cl (int i, int j, int k) {
	if (i<0 || i>=nblig || j<0 || j>=nbcol || k<0 || k>=(int)nclas) {
		cout<<" debordement d''indice dans ("<<i<<","<<j<<","<<k<<")\n";
		if (i<0) i=0;
		if (i>=nblig) i=nblig-1;
		if (j<0) j=0;
		if (j>=nbcol) j=nbcol-1;
		if (k<0) k=0;
		if (k>=(int)nclas) k=nclas-1;
	}
	double *adval=(double*)T_A[i*nbcol+j];
	adval=adval+k;
	return *adval;
}

/* création de l'image des classes */
template <class T> imadata<float> imaBR<T>::creeimclas () {
	int i,j,k;
	imadata<float> imacl(nblig,nbcol,nclas);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) 
			for (k=0; k<nclas; k++) imacl(i,j,k)=T_A[i*nbcol+j][k];
	return imacl;
}

/* calcul de la composition des pixels à partir d'une image des labels HR */
template <class T> void imaBR<T>::labelHR_to_comppixBR(imalabels imalab) {
	int nl_BR=nblig, nc_BR=nbcol, nl_HR=imalab.nlig(), nc_HR=imalab.ncol(),
			np_BR=nl_BR*nc_BR, np_HR=nl_HR*nc_HR, nHRl=nl_HR/nl_BR, nHRc=nc_HR/nc_BR;
	cout<<" #lig HR/BR = "<<nl_HR<<"/"<<nl_BR<<" = "<<nHRl<<" ?\n";
	cout<<" #col HR/BR = "<<nc_HR<<"/"<<nc_BR<<" = "<<nHRc<<" ?\n";
	unsigned long int *H_clas=NULL;
	int i,j,n,iHR,jHR,nHR;
	unsigned int nmax=256,kmin=nmax,kmax=0,k,k2;
	H_clas=new unsigned long int[nmax];
	for (k=0; k<nmax; k++) H_clas[k]=0;
	for (i=0; i<nl_HR; i++)
		for (j=0; j<nc_HR; j++) {
			k=imalab(i,j);
			if (k<nmax) H_clas[k]++;
			if (k<kmin) kmin=k;
			if (k>kmax) kmax=k;
		}
	BYTE *Teq_cl=new BYTE[nmax];
	for (k=0; k<nmax; k++)
		if (H_clas[k]==0) Teq_cl[k]=0;
		else Teq_cl[k]=k;
	for (k=kmax; k>=0; k--) {
		if (H_clas[k]==0) {
			for (k2=k; k2<kmax; k2++) {
				H_clas[k2]=H_clas[k2+1];
				Teq_cl[k]=k-1;
			}
			kmax--;
		}
	}
	if ((kmax+1)!=nclas) {
		int np=nl_BR*nc_BR;
		for (i=0; i<np; i++) if (T_A[i]!=NULL) delete[] T_A[i];
		nclas=kmax+1;
		for (i=0; i<np; i++) T_A[i]=new double[nclas];
		}
	cout<<" # de classes de l'image des labels HR = "<<nclas<<"\n";
	if (H_clas!=NULL) {delete[] H_clas; H_clas=NULL;}
	n=0;
	for (i=0; i<nl_BR; i++)
		for (j=0; j<nc_BR; j++) { //cout<<n<<" ";//cout<<i<<" "<<j<<" "<<n<<"\n";
			for (k=0; k<nclas; k++) T_A[n][k]=0;
			nHR=0;
			for (iHR=i*nHRl; iHR<mini((i+1)*nHRl,nl_HR); iHR++)
				for (jHR=j*nHRc; jHR<mini((j+1)*nHRc,nc_HR); jHR++) {
					nHR++;
					k=Teq_cl[imalab(iHR,jHR)];
					if (k>=0 && k<nclas) T_A[n][k]+=1;
					else cout<<" Pb dans le label "<<k<<" <0 ou >=nclas="<<nclas<<"\n";
				}
			for (k=0; k<nclas; k++) T_A[n][k]/=nHR;
			n++;
		}
	if (Teq_cl!=NULL) {delete[] Teq_cl; Teq_cl=NULL;}
	A_valid=1;
}

/* calcul de la composition des pixels à partir d'une image des labels HR décalée par rapport à l'image BR */
/*template <class T> void imaBR<T>::labelHRd_to_comppixBR(imalabels imalab, int , int , bool iaf) {
	int nl_BR=nblig, nc_BR=nbcol, nl_HR=imalab.nlig(), nc_HR=imalab.ncol(),
		np_BR=nl_BR*nc_BR, np_HR=nl_HR*nc_HR, nHRl=nl_HR/nl_BR, nHRc=nc_HR/nc_BR;
	if (iaf) {
		cout<<" #lig HR/BR = "<<nl_HR<<"/"<<nl_BR<<" = "<<nHRl<<" ?\n";
		cout<<" #col HR/BR = "<<nc_HR<<"/"<<nc_BR<<" = "<<nHRc<<" ?\n";
	}
		for (dl=-dlm/2; dl<dlm-dlm/2; dl++)
			for (dc=-dcm/2; dc<dcm-dcm/2; dc++) {
				dlc=(dl+dlm/2)*dcm+dc+dcm/2;
				for (i=0; i<nlig_e-nl+1; i+=nl)
					for (j=0; j<ncol_e-nc+1; j+=nc) {
						i2=i/nl; j2=j/nc; k2=i2*ncolBR+j2;
						for (k=0; k<nq; k++) T_M_alpha[dlc](k2,k)=0;
						if (i+dl>=0 && i+dlm+dl<nlig_e && j+dc>=0 && j+dcm+dc<ncol_e) {
							n=0;
							for (ii=0; ii<nl; ii++)
								for (jj=0; jj<nc; jj++) {
									k=imaHRl(i+ii+dl,j+jj+dc);
									T_M_alpha[dlc](k2,k)+=1;
									n++;
								}
							for (k=0; k<nq; k++) T_M_alpha[dlc](k2,k)/=n;
						} 
						else T_M_alpha[dlc](k2,0)=1;
					}
				for (k2=0; k2<npixBR; k2++) {
					x=0.;
					for (k=0; k<nq; k++) x+=fabs(T_M_alpha[dlc](k2,k));
					if (fabs(x-1)>prec) sortie<<" pb ds la matrice des proportions au pixel "<<k2<<" : somme_ak = "<<x<<"\n";
//					else sortie<<" matrice des proportions OK au pixel "<<k2<<" : somme_ak = "<<x<<"\n";
				}
//				sortie<<" matrice des proportions calculee\n";
			}

	n=0;
	for (i=0; i<nl_BR; i++)
		for (j=0; j<nc_BR; j++) { //cout<<n<<" ";//cout<<i<<" "<<j<<" "<<n<<"\n";
			for (k=0; k<nclas; k++) T_A[n][k]=0;
			nHR=0;
			for (iHR=i*nHRl; iHR<mini((i+1)*nHRl,nl_HR); iHR++)
				for (jHR=j*nHRc; jHR<mini((j+1)*nHRc,nc_HR); jHR++) {
					nHR++;
					k=Teq_cl[imalab(iHR,jHR)];
					if (k>=0 && k<nclas) T_A[n][k]+=1;
					else cout<<" Pb dans le label "<<k<<" <0 ou >=nclas="<<nclas<<"\n";
				}
			for (k=0; k<nclas; k++) T_A[n][k]/=nHR;
			n++;
		}
	if (Teq_cl!=NULL) delete[] Teq_cl;
	A_valid=1;
}*/

/* desaggregation ou estimation des caracteristiques des classes connaissant la composition des pixels */
template <class T> void imaBR<T>::desaggregation(matrice2D<bool> &V) {
	int i,j,n=0, nl_BR=nblig, nc_BR=nbcol, dim=nbcanaux, np_BR=0;
	unsigned int k;
	matrice2D<bool> M_V(nl_BR,nc_BR);
	if (V.nlig()==nl_BR && V.ncol()==nc_BR) M_V=V;
	for (i=0; i<nl_BR; i++)
		for (j=0; j<nc_BR; j++) if (!M_V(i,j)) np_BR++;
	matrice2D<float> M_A(np_BR,nclas), M_At(nclas,np_BR), M_X(np_BR,dim), M_y(nclas,dim);
	n=0;
	for (i=0; i<nl_BR; i++)
		for (j=0; j<nc_BR; j++)
			if (!M_V(i,j)) {
				for (k=0; k<(unsigned int)dim; k++) M_X(n,k)=(*this)(i,j,k);
				for (k=0; k<nclas; k++) M_A(n,k)=M_At(k,n)=(float)(T_A[i*nc_BR+j][k]);
				n++;
			}
	matrice2D<float> M_AAt(M_At,M_A,nclas,np_BR,nclas);
	bool OK;
	matrice2D<float> M_AAt_1=M_AAt.inverse(OK);
	cout<<" verification de l'inversion de A.At : inversion ";
	matrice2D<float> Id(nclas,nclas); for (k=0; k<nclas; k++) Id(k,k)=1;
	(matrice2D<float>(M_AAt,M_AAt_1,nclas,nclas,nclas)==Id&&matrice2D<float>(M_AAt_1,M_AAt,nclas,nclas,nclas)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
//	(matrice2D<float>(M_AAt,M_AAt_1,nclas,nclas,nclas)).affiche();
//	(matrice2D<float>(M_AAt_1,M_AAt,nclas,nclas,nclas)).affiche();
	M_y=matrice2D<float>(matrice2D<float>(M_AAt_1,M_At,nclas,nclas,np_BR),M_X,nclas,np_BR,dim);
	cout<<" caracteristiques des classes obtenues connaissant A\n"; M_y.affiche();
	if (T_y==NULL) {
		T_y=new double*[nclas];
		for (k=0; k<nclas; k++) T_y[k]=new double[dim];
	}
	for (k=0; k<nclas; k++)
		for (j=0; j<dim; j++) T_y[k][j]=M_y(k,j);
//	cout<<" matrice de l'erreur par pixel :\n";
//	(matrice2D<float>(matrice2D<float>(M_A,M_y,np_BR,nclas,dim)-M_X)).affiche();
//	matrice2D<float> M_d(M_A,M_y,np_BR,nclas,dim);(M_d-M_X).affiche();
}

/* erreur quadratique par pixel */
template <class T> matrice2D<double> imaBR<T>::err2_ppix(bool isauv) {
	int nl=nblig, nc=nbcol, np=nl*nc, dim=nbcanaux, i,j,k,n;
	matrice2D<double> M_A(np,nclas,T_A), M_y(nclas,dim,T_y), M_Xs(M_A,M_y,np,nclas,dim);
	matrice2D<double> M_d(nl,nc);
	double x;
	for (i=0; i<nl; i++)
		for (j=0; j<nc; j++) {
			n=i*nc+j;
			x=0;
			for (k=0; k<dim; k++) x+=pow(M_Xs(n,k)-(*this)(i,j,k),2);
			M_d(i,j)=pow(x,0.5);
		}
	if (isauv) {
		imadata<double> imatemp(M_d);
		((imadata<BYTE>)imatemp).sauve_ImaPGM("./image_err2ppix.pgm");
	}
	return M_d;
}

/* estimation de la composition des pixels connaissant les caracteristiques des classes */
template <class T> void imaBR<T>::decomposition() {
	int nl_BR=nblig, nc_BR=nbcol, np_BR=nl_BR*nc_BR;
	int i,j,n,n_it=0,init;
	unsigned int k;
//	int nb_lib=nclas+2; // # de variables à optimiser lors de la décomposition = # de classes + paramètre lagrangien (optim. sous contrainte) + 1 pour commencer en 1 (!=0)
	int nb_lib=nclas; // # de variables à optimiser lors de la décomposition = # de classes -1 (contrainte de stochasticité) + 1 pour commencer en 1 (!=0)
	float *a_s=new float[nb_lib]; // composition du pixel s : a_s[0]=s, k ds [1,nclas] a_s[k]=pct de la classe k dans s
	extern int dim1, dim2; dim1=nclas; dim2=nbcanaux; // variables externes définissant la fonction coût à optimiser
	extern float **M; M=new float*[dim1]; for (i=0; i<dim1; i++) M[i]=new float[dim2];
	for (k=0; k<nclas; k++) 
		for (j=0; j<nbcanaux; j++) M[k][j]=(float)(T_y[k][j]); 
	extern float *V; V=new float[dim2];
	extern float *W; W=new float[dim1];
	float **xi, fmin=1.e+9f, eps=1.e-9f, s, t, distmin; 
	xi=new float*[nb_lib]; for (i=0; i<nb_lib; i++) xi[i]=new float[nb_lib];
	for (n=0; n<np_BR; n++) {
		for (j=0; j<nbcanaux; j++) V[j]=(*this)(n/nbcol,n%nbcol,j);
//		a_s[0]=n; for (i=1; i<=nclas; i++) a_s[i]=1/nclas; a_s[nb_lib-1]=0.f;
		a_s[0]=(float)n; for (k=1; k<nclas; k++) a_s[k]=(float)(1/nclas); 
		distmin=0.;                                              // choix d'une 'bonne' initialisation
		for (j=0; j<nbcanaux; j++) {
			t=0.;
			for (k=0; k<nclas; k++) t+=M[k][j];
			t/=nclas;
			distmin+=pow(V[j]-t,2);
		}
		init=nclas;
		for (k=0; k<nclas; k++) {
			t=0.;
			for (j=0; j<nbcanaux; j++) t+=pow(V[j]-M[k][j],2);
			if (t<distmin) {distmin=t; init=i;}
		}
		if (init<nclas) {
			for (k=1; k<nclas; k++) a_s[k]=(float)(1/nclas);
			if (init<nclas-1) a_s[init+1]=1.f;
		}
//		if (n==0) {cout<<" initialisation "<<init<<" : "; for (i=1; i<nclas; i++) cout<<a_s[i]<<" "; cout<<"\n";}
		for (k=0; k<nclas; k++) W[k]=(float)(T_A[n][k]);
		for (i=0; i<nb_lib; i++) 
			for (j=0; j<nb_lib; j++) xi[i][j]=0;
		for (i=1; i<nb_lib; i++) xi[i][i]=1;
		powell(a_s,xi,nb_lib-1,eps,&n_it,&fmin,err2);
		s=0.;
		for (k=1; k<nclas; k++) {T_A[n][k-1]=a_s[k]; s+=a_s[k];} T_A[n][nclas-1]=1.-s;
		cout<<" pix. "<<setw(3)<<(int)a_s[0]<<" min. "<<fixed<<setw(6)<<setprecision(3)<<fmin<<" en "; 
		for (k=0; k<nclas; k++) cout<<setw(6)<<setprecision(3)<<T_A[n][k]<<" "; cout<<" en "<<n_it<<" it.\n";
/*		if (n==0) {
			for (j=0; j<nbcanaux; j++) cout<<V[j]<<" "; cout<<"\n";
			for (i=0; i<nclas; i++) {
				float xx=0.;
				for (j=0; j<nbcanaux; j++) xx+=pow(V[j]-M[i][j],2);
				cout<<" distance a la classe "<<i<<" = "<<xx<<"\n";
			}
		}*/
	}
	for (i=0; i<nb_lib; i++) if (xi[i]!=NULL) delete[] xi[i]; if (xi!=NULL) delete[] xi;
	if (a_s!=NULL) delete[] a_s;
	for (i=0; i<dim1; i++) if (M[i]!=NULL) delete[] M[i]; if (M!=NULL) delete[] M;
	if (V!=NULL) delete[] V;
	if (W!=NULL) delete[] W;
}

#endif