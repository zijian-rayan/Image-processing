#ifndef _FONCTIONS_H
#define _FONCTIONS_H
#define _CRT_RAND_S

#include <iostream>
#include <iomanip>
using namespace std;
#include <limits>
#include "math.h"
#include "constantes.h"

// divers
template <class T> T mini(const T x1, const T x2) {
	if (x1<x2) return x1;
	else return x2;
}
template <class T> T maxi(const T x1, const T x2) {
	if (x1>x2) return x1;
	else return x2;
}
template <class T> int around(const T x) {
	int int_x=(int)x;
	if (x>=0) {
		if ((x-int_x)<=0.5) return int_x;
		else return (int_x+1);
	} else {
		if ((int_x-x)<=0.5) return int_x;
		else return (int_x-1);
	}
}
template <class T> int signe(const T x) {
	if (x>=0) return 1;
	else return -1;
}
template <class T> T modulo2PI(const T x, bool icentre=0) {
	const double PI2=2*PI;
	const int itmax=10;
	double y=x;
	int it=0;
	while (y>=PI2 && it<itmax) {y-=PI2; it++;}
	while (y<0 && it<itmax) {y+=PI2; it++;}
	if (it>=itmax) {cout<<"Pb ds modulo2PI : it = "<<it<<", x = "<<x<<", y = "<<y<<"\n"; char aa; cin>>aa;}
	if (icentre) y-=PI;
	return (T)y;
}
template <class T> bool Porte(const T x) {return (fabs(x)<=0.5);}

template <class T> double logistic(const T x, const T x0, const double a=1., const double b=1.) {return b/(1.+exp(-a*(x-x0)));}

/*int mini(int, int);
unsigned int mini(unsigned int, unsigned int);
float mini(float, float);
double mini(double, double);
int maxi(int, int);
unsigned int maxi(unsigned int, unsigned int);
float maxi(float, float);
double maxi(double, double);
int around(float);
int around(double);*/
bool puissance_de_2(const unsigned long int);
void AffichBinaryNumber (unsigned int);
void AffichBinaryNumber (unsigned int, ostream&);
double val_Cij(const int, const int);
long double facto(const int);

void fct_Tab_permut (const int, int*, const int, int**, int&, int&);
void Tab_permutation (const int, int**);

// géométrie
template<class T> bool intersection_droites_polaires2 (float &fRhoI, float &fThetaI, T _fRho1, T _fTheta1, T _fRho2, T _fTheta2) {
	const double eps=1.e-3;
	bool Iok=0;
	float fRho1=(float)_fRho1, fTheta1=(float)_fTheta1, fRho2=(float)_fRho2, fTheta2=(float)_fTheta2;
	double alpha=fRho1*cos(fTheta2)-fRho2*cos(fTheta1), beta=fRho1*sin(fTheta2)-fRho2*sin(fTheta1), coefnorm=pow(alpha*alpha+beta*beta,0.5);
	if (coefnorm!=0) {
		fThetaI=-(float)atan2(alpha/coefnorm,beta/coefnorm);
		if (fabs(fRho1*cos(fThetaI-fTheta2)-fRho2*cos(fThetaI-fTheta1))>eps) {
			cout<<" Pb dans calcul de l'intersection des droites de coordonnees polaires : fThetaI="<<fThetaI/PI*180;
			cout<<", fRho1*cos(fThetaI-fTheta2) = "<<fRho1*cos(fThetaI-fTheta2)<<", fRho2*cos(fThetaI-fTheta1) = "<<fRho2*cos(fThetaI-fTheta1);
			cout<<" dif = "<<fRho1*cos(fThetaI-fTheta2)-fRho2*cos(fThetaI-fTheta1)<<"\n";
		} else {
			if (cos(fThetaI-fTheta1)!=0) fRhoI=fRho1/cos(fThetaI-fTheta1);
			else fRhoI=fRho2/cos(fThetaI-fTheta2);
			Iok=1;
		}
	}
	else cout<<" alpha = "<<alpha<<", beta = "<<beta<<"\n";
	return Iok;
}

double aire_parallelogramme (double, double, double, double, double, double, double, double, double=1.e-6);
//double aire_parallelogramme (double, double, double, double, double, double);
template<class T> T aire_parallelogramme (T xA, T yA, T xB, T yB, T xC, T yC) {
	return xA*(yB-yC)-xB*(yA-yC)+xC*(yA-yB);
}

// tirages aléatoires, dist histo
double tirage_gauss(const double=1, const double=0);
double tirage_gauss_rapide(const double=1, const double=0);
unsigned int tirage_proba(const double*, const unsigned int);

bool histo_norm(float*, int, bool=0);
double distBhattacharyya_histo(const float*, const float*, int);

template<class T> double distance_euclidienne(const T *v1,const T *v2,const int sz_v) {
	double xx=0.;
	int n;
	for (n=0; n<sz_v; n++) xx+=pow((double)v1[n]-v2[n],2.);
	return xx;
}
template<class T> double distance_Battacharrya(const T *v1,const T *v2,const int sz_v) {
	double xx=0., s_v1=0., s_v2=0.;
	int n;
	for (n=0; n<sz_v; n++) {xx+=(double)v1[n]*v2[n]; s_v1+=(double)v1[n]; s_v2+=(double)v2[n];}
	xx=pow(xx/s_v1/s_v2,0.5);
	return 1.-xx;
}

// fonctions de tri par ordre croissant
template<class T> void tri_par_insertion(T *t, int n) {
	int i,j;
	if (n>1) tri_par_insertion(t,n-1);
	T x=t[n-1];
	i=0;
	while (t[i]<x) i++;
	for (j=n-1; j>i; j--) t[j]=t[j-1];
	t[i]=x;
}
template<class T> void tri_par_selection(T *t, const int n) {
	int i,i_max=n-1;
	T xmax=t[n-1];
	for (i=0; i<n-1; i++)
		if (t[i]>xmax) {
			i_max=i;
			xmax=t[i];
		}
	t[i_max]=t[n-1];
	t[n-1]=xmax;
	if (n>2) tri_par_selection(t,n-1);
}
template<class T> void ordre_decroiss_par_selection(T *t, int *r, const int n) {
	int i,i_min=n-1,rmin=r[i_min];
	T xmin=t[i_min];
	for (i=0; i<n-1; i++)
		if (t[i]<xmin) {
			i_min=i;
			xmin=t[i];
			rmin=r[i];
		}
	t[i_min]=t[n-1];
	t[n-1]=xmin;
	r[i_min]=r[n-1];
	r[n-1]=rmin;
	if (n>2) ordre_decroiss_par_selection(t,r,n-1);
}
template<class T> void tri_bulle(T *t, const int n) {
	bool permut;
	int i;
	T x;
	do {
		permut=0;
		for (i=0; i<n-1; i++)
			if (t[i]>t[i+1]) {
				x=t[i];
				t[i]=t[i+1];
				t[i+1]=x;
				permut=1;
			}
	} while (permut);
}
/*template<class T> int partition(T *t, const int n) {
	int i, n_inf=0, pos_pivot=0; // ou bien : pos_pivot=rand()%n
	T x, pivot=t[pos_pivot];
	for (i=1; i<n; i++)
		if (t[i]<pivot) {
			x=t[i];
			t[i]=t[n_inf];
			t[n_inf]=x;
			if (pos_pivot==n_inf) pos_pivot=i;
			n_inf++;
		}
	if (pos_pivot!=n_inf) {
		x=t[pos_pivot];
		t[pos_pivot]=t[n_inf];
		t[n_inf]=x;
	}
	return n_inf;
}
template<class T> void tri_rapide(T *t, const int n) {
	if (n>1) {
		int i,pos_pivot=partition(t,n);
		if (pos_pivot>1) {
			tri_rapide(t,pos_pivot);
		}
		if (pos_pivot<n-1) {
			int nbis=n-pos_pivot-1;
			T *tbis=new T[nbis];
			for (i=0; i<nbis; i++) tbis[i]=t[i+pos_pivot+1];
			tri_rapide(tbis,nbis);
			for (i=0; i<nbis; i++) t[i+pos_pivot+1]=tbis[i];
			if (tbis!=NULL) delete[] tbis;
		}
	}
}*/
struct val_pos {double val; unsigned int pos;};
/*template<> int partition<val_pos>(val_pos *t, const int n) {
	int i, k, n_inf=0, pos_pivot=0; // ou bien : pos_pivot=rand()%n
	double x, pivot=t[pos_pivot].val;
	for (i=1; i<n; i++)
		if (t[i].val<pivot) {
			x=t[i].val; k=t[i].pos;
			t[i].val=t[n_inf].val; t[i].pos=t[n_inf].pos;
			t[n_inf].val=x; t[n_inf].pos=k;
			if (pos_pivot==n_inf) pos_pivot=i;
			n_inf++;
		}
	if (pos_pivot!=n_inf) {
		x=t[pos_pivot].val; k=t[pos_pivot].pos;
		t[pos_pivot].val=t[n_inf].val; t[pos_pivot].pos=t[n_inf].pos;
		t[n_inf].val=x; t[n_inf].pos=k;
	}
	return n_inf;
}
template<> void tri_rapide<val_pos>(val_pos *t, const int n) {
	if (n>1) {
		int i,pos_pivot=partition(t,n);
		if (pos_pivot>1) {
			tri_rapide(t,pos_pivot);
		}
		if (pos_pivot<n-1) {
			int nbis=n-pos_pivot-1;
			val_pos *tbis=new val_pos[nbis];
			for (i=0; i<nbis; i++) {tbis[i].val=t[i+pos_pivot+1].val; tbis[i].pos=t[i+pos_pivot+1].pos;}
			tri_rapide(tbis,nbis);
			for (i=0; i<nbis; i++) {t[i+pos_pivot+1].val=tbis[i].val; t[i+pos_pivot+1].pos=tbis[i].pos;}
			if (tbis!=NULL) delete[] tbis;
		}
	}
}
*/
/*void tri_par_insertion(double*, const int);
//void tri_par_insertion(int *, int);
void tri_par_selection(double*, const int);
//void tri_par_selection(int *, int);
void ordre_decroiss_par_selection(double*, int*, int);
void tri_bulle(double*, const int);
//void tri_bulle(int*, int);*/
int partition(double*, const int);
void tri_rapide(double*, const int);
//int partition(double*, const int, const int, const bool=0);
//void tri_rapide_partiel(double*, const int, const int, const bool=0);
//void tri_rapide(double*, const int, const bool=0);
int partition(double*, const int, const int, const bool);
void tri_rapide_partiel(double*, const int, const int, const bool);
void tri_rapide(double*, const int, const bool);
//int partition(val_pos*, const int);
int partition(val_pos*, const int, const int, const bool=0);
void tri_rapide_partiel(val_pos*, const int, const int, const bool=0);
void tri_rapide(val_pos*, const int, const bool=0);

// divers++
float regres_lin (float**, const int, float&, float&);
float algoKuhn (float**, const int, bool=0);
float algoKuhn2 (float**, const int, bool=0);
template<class T> void aff_Tab2D (T** Tab, const int d) {
	cout<<"_______________________________\n";
	for (int i=0; i<d; i++) {for (int j=0; j<d; j++) cout<<setw(5)<<Tab[i][j]<<" "; cout<<"\n";}
	cout<<"_______________________________\n";
}
template<class T> void aff_Tabpseudo2D (T* Tab, const int d) {
	cout<<"_______________________________\n";
	for (int i=0; i<d; i++) {for (int j=0; j<d; j++) cout<<setw(5)<<Tab[i*d+j]<<" "; cout<<"\n";}
	cout<<"_______________________________\n";
}
template<class T> void soustrait_min_col_lig(T** Tab, const int d) {
	T x;
	int i,j;
	for (i=0; i<d; i++) {
		x=Tab[i][i]/2;
		for (j=0; j<i; j++)
			if (Tab[i][j]<x || Tab[j][i]<x) x=mini(Tab[i][j],Tab[j][i]);
		for (j=i+1; j<d; j++)
			if (Tab[i][j]<x || Tab[j][i]<x) x=mini(Tab[i][j],Tab[j][i]);
		for (j=0; j<d; j++) {Tab[i][j]-=x; Tab[j][i]-=x;}
	}
}
template<class T> int barrezeros_Kuhn (T** Tab, const int d, bool* zerobarre, bool* zeromarqu, int n, const int ilig, const int jcol) {
	const double epsilon=1.e-10;
	int i,j;
	for (j=0; j<d; j++) 
		if (fabs(Tab[ilig][j])<=epsilon && !zerobarre[ilig*d+j] && !zeromarqu[ilig*d+j]) {
			zerobarre[ilig*d+j]=1; n--; /*cout<<" 0 barre en ("<<ilig+1<<","<<j+1<<")\n";*/}
	for (i=0; i<d; i++) 
		if (fabs(Tab[i][jcol])<=epsilon && !zerobarre[i*d+jcol] && !zeromarqu[i*d+jcol]) {
			zerobarre[i*d+jcol]=1; n--; /*cout<<" 0 barre en ("<<i+1<<","<<jcol+1<<")\n";*/}
	return n;
}
template<class T> int barrezeros_Kuhn2 (T** Tab, const int d, bool* zerobarre, bool* zeromarqu, int n, const int ilig, const int jcol) {
	const double epsilon=1.e-10;
	int i,j;
//	for (j=0; j<d; j++) 
//		if (Tab[ilig][j]==0 && !zerobarre[ilig*d+j] && !zeromarqu[ilig*d+j]) {
//			zerobarre[ilig*d+j]=1; n--; /*cout<<" 0 barre en ("<<ilig+1<<","<<j+1<<")\n";*/}
//	for (i=0; i<d; i++) 
//		if (Tab[i][jcol]==0 && !zerobarre[i*d+jcol] && !zeromarqu[i*d+jcol]) {
//			zerobarre[i*d+jcol]=1; n--; /*cout<<" 0 barre en ("<<i+1<<","<<jcol+1<<")\n";*/}
	for (j=0; j<d; j++) 
		if (fabs(Tab[ilig][j])<=epsilon && !zerobarre[ilig*d+j] && !zeromarqu[ilig*d+j]) {
			zerobarre[ilig*d+j]=1; n--; /*cout<<" 0 barre en ("<<ilig+1<<","<<j+1<<")\n";*/}
	for (i=0; i<d; i++) 
		if (fabs(Tab[i][ilig])<=epsilon && !zerobarre[i*d+ilig] && !zeromarqu[i*d+ilig]) {
			zerobarre[i*d+ilig]=1; n--; /*<<" 0 barre en ("<<i+1<<","<<ilig+1<<")\n";*/}
	for (i=0; i<d; i++) 
		if (fabs(Tab[i][jcol])<=epsilon && !zerobarre[i*d+jcol] && !zeromarqu[i*d+jcol]) {
			zerobarre[i*d+jcol]=1; n--; /*cout<<" 0 barre en ("<<i+1<<","<<jcol+1<<")\n";*/}
	for (j=0; j<d; j++) 
		if (fabs(Tab[jcol][j])<=epsilon && !zerobarre[jcol*d+j] && !zeromarqu[jcol*d+j]) {
			zerobarre[jcol*d+j]=1; n--; /*cout<<" 0 barre en ("<<jcol+1<<","<<j+1<<")\n";*/}
	return n;
}
template<class T> int coord_zerodiag_argmin_nbzeros(T** Tab, const int d) {
	const double epsilon=1.e-10;
	int nzeromin=2*d+2, ilig=0, n, i, j;
	float xmin=Tab[0][0];
	for (i=1; i<d; i++) 
		if (Tab[i][i]<xmin) {xmin=Tab[i][i]; ilig=i;}
//	cout<<" valeur min sur la diagonale = "<<xmin<<"\n";
	for (i=0; i<d; i++)
		if (Tab[i][i]==xmin) {//cout<<" atteinte en ilig = icol = "<<i;
			n=0;
			for (j=0; j<d; j++) {if (Tab[i][j]==0) n++; if (Tab[j][i]==0) n++;} //cout<<" avec "<<n<<" zeros\n";
			if (n<nzeromin) {nzeromin=n; ilig=i; /*cout<<" i.e. moins de zeros que precedemment -> ilig = "<<ilig<<"\n";*/}
		}
	return ilig;
}
template<class T> int coord_zero_argmin_nbzeros_Kuhn(T** Tab, const int d, const bool* zeromarqu, const bool* zerobarre, const int* i_0, const int* j_0, const int nmax) {
	const double epsilon=1.e-10;
	int nzeromin=nmax+2, k=-1, n, ilig, jcol, i, nn;
	for (n=0; n<=nmax; n++) {
		ilig=i_0[n]; jcol=j_0[n]; 
		if (!zeromarqu[ilig*d+jcol] && !zerobarre[ilig*d+jcol]) { //cout<<" zero "<<n<<" "<<ilig+1<<" "<<jcol+1<<" ";
			nn=0;
			for (i=0; i<d; i++) {
				if (fabs(Tab[ilig][i])<=epsilon && !zeromarqu[ilig*d+i] && !zerobarre[ilig*d+i]) nn++;
				if (fabs(Tab[i][jcol])<=epsilon && !zeromarqu[i*d+jcol] && !zerobarre[i*d+jcol]) nn++;
			}
			if (nn>0 && nn<nzeromin) {nzeromin=nn; k=n;} //cout<<" nn = "<<nn<<" nzeromin "<<nzeromin<<" k "<<k<<"\n";
		}
	}
	return k;
}
template<class T> int coord_zero_argmin_nbzeros_Kuhn2(T** Tab, const int d, const bool* zeromarqu, const bool* zerobarre, const int* i_0, const int* j_0, const int nmax) {
	const double epsilon=1.e-10;
	int nzeromin=nmax+2, k=-1, n, ilig, jcol, i, nn;
	for (n=0; n<=nmax; n++) {
		ilig=i_0[n]; jcol=j_0[n]; 
		if (!zeromarqu[ilig*d+jcol] && !zerobarre[ilig*d+jcol]) { //cout<<" zero "<<n<<" "<<ilig+1<<" "<<jcol+1<<" ";
			nn=0;
			for (i=0; i<d; i++) {
				if (fabs(Tab[ilig][i])<=epsilon && !zeromarqu[ilig*d+i] && !zerobarre[ilig*d+i]) nn++;
				if (fabs(Tab[i][jcol])<=epsilon && !zeromarqu[i*d+jcol] && !zerobarre[i*d+jcol]) nn++;
			}
			if (nn>0 && nn<nzeromin) {nzeromin=nn; k=n;} //cout<<" nn = "<<nn<<" nzeromin "<<nzeromin<<" k "<<k<<"\n";
		}
	}
	return k;
}

// fonctions gamma
/*long double gammln(long double);
void gser(long double&, long double, long double, long double&);
void gcf(long double&, long double, long double, long double&);
long double gammq(long double, long double);
long double gammp(long double, long double);
void ln_gcf(long double&, long double, long double, long double&);
long double factln(int n);
void ln_gser(long double&, long double, long double, long double&);
long double ln_gammq(long double, long double);
long double ln_gammp(long double, long double);*/

// fonctions d'optimisation de Numerical recipes
/*void mnbrak(float*, float*, float*, float*, float*, float*, float (*func)(float));
float golden(float, float, float, float (*f) (float), float, float*);
float brent(float, float, float,float (*f)(float), float,float*);
float *vector(int, int);
void free_vector(float*, int, int);
void linmin(float[], float[], int, float*, float (*f)(float));
float f1dim(float);
void powell(float[], float**, int, float, int*, float*, float (*f)(float[]));
float err2(float[]);*/

#endif