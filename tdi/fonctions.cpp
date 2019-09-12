#include "fonctions.h"
/*
int mini(const int n1, const int n2) {
	if (n1<n2) return n1;
	else return n2;
}
unsigned int mini(const unsigned int n1, const unsigned int n2) {
	if (n1<n2) return n1;
	else return n2;
}
float mini(const float x1, const float x2) {
	if (x1<x2) return x1;
	else return x2;
}
double mini(const double x1, const double x2) {
	if (x1<x2) return x1;
	else return x2;
}

int maxi(const int n1, const int n2) {
	if (n1>n2) return n1;
	else return n2;
}
unsigned int maxi(const unsigned int n1, const unsigned int n2) {
	if (n1>n2) return n1;
	else return n2;
}
float maxi(const float x1, const float x2) {
	if (x1>x2) return x1;
	else return x2;
}
double maxi(const double x1, const double x2) {
	if (x1>x2) return x1;
	else return x2;
}

int around(const float x) {
	int int_x=(int)x;
	if (x>=0) {
		if ((x-int_x)<=0.5) return int_x;
		else return (int_x+1);
	} else {
		if ((int_x-x)<=0.5) return int_x;
		else return (int_x-1);
	}
}
int around(const double x) {
	int int_x=(int)x;
	if (x>=0) {
		if ((x-int_x)<=0.5) return int_x;
		else return (int_x+1);
	} else {
		if ((int_x-x)<=0.5) return int_x;
		else return (int_x-1);
	}
}
*/
bool puissance_de_2(const unsigned long int l) {
//	cout<<" puissance_de_2 "<<l<<" "<<(l^(l-1))<<" "<<(l+(l-1))<<"\n";
	return (l!=0 && (l^(l-1))==(l+(l-1)));
}

void AffichBinaryNumber (unsigned int n) {
	unsigned int mask, idec=around(log((float)n)/log(2.f))+2;
	idec=16;
	mask = 1 << idec; // l'algorithme ne marche pas quand n=0 (boucle infinie)
	if (n == 0) cout<<"0";
	else {
		while(!((mask) & n)) mask >>= 1; // on affiche a partir du premier "1"
		do {
			cout<<((n & mask)?"1":"0");
		}	while(mask >>= 1); // on decale vers la droite 'mask'
	}
}

void AffichBinaryNumber (unsigned int n, ostream& f_out) {
	unsigned int mask, idec=around(log((float)n)/log(2.f))+2;
	idec=16;
	mask = 1 << idec; // l'algorithme ne marche pas quand n=0 (boucle infinie)
	if (n == 0) f_out<<"0";
	else {
		while(!((mask) & n)) mask >>= 1; // on affiche a partir du premier "1"
		do {
			f_out<<((n & mask)?"1":"0");
		}	while(mask >>= 1); // on decale vers la droite 'mask'
	}
}
double val_Cij(const int _i, const int _j) {
	int k,i=_i,j=_j;
	if (i>j) {k=i; i=j; j=k;}
	double z=1;
	int n=i;
	if (j-n>n) n=j-n;
	for (k=j; k>n; k--) z*=k;
	for (k=j-n; k>1; k--) z/=k;
	return z;
}

long double facto(const int i) {
	int k;
	long double z=1;
	if (i<0) return -1.;
	else {
		for (k=2; k<=i; k++) z*=k;
		return z;
	}
}

void fct_Tab_permut (const int cardL, int *L, const int n, int **T, int& ligT, int& colT) {
	int i,j;
	if (cardL>0) {
		int *L2=new int[cardL-1];
		for (i=0; i<cardL; i++) {
			if (colT<0) colT=n-cardL; //cout<<" reintialisation des colonnes a "<<colT<<"\n";
			for (j=0; j<facto(cardL-1); j++)
				T[ligT+j][colT]=L[i]; //cout<<" T["<<ligT+j<<"]["<<colT<<"] = "<<T[ligT+j][colT]<<"\n";
			colT++;
			for (j=0; j<cardL; j++) {
				if (j<i) L2[j]=L[j];
				if (j>i) L2[j-1]=L[j];
			}
			fct_Tab_permut (cardL-1,L2,n,T,ligT,colT);
		}
		if (L2!=NULL) delete[] L2;
	} else {
		ligT++;
		colT=-1;
	}
}

void Tab_permutation (const int n, int **T) {
	int ligT=0,colT=-1,i,*L=new int[n];
	for (i=0; i<n; i++) L[i]=i;
	fct_Tab_permut (n,L,n,T,ligT,colT);
	if (L!=NULL) delete[] L;
/*	for (i=0; i<facto(n);i++) {
		for (int j=0; j<n; j++) cout<<" "<<T[i][j]<<"  |   ";
		cout<<"\n";
	}
	char aa; cin>>aa;*/
}

double aire_parallelogramme (double xA, double yA, double xB, double yB, double xC, double yC, double xD, double yD, double prec) {
	if ((fabs(fabs(yB-yA)-fabs(yD-yC))>prec || fabs(fabs(xB-xA)-fabs(xD-xC))>prec) && 
			(fabs(fabs(yC-yA)-fabs(yD-yB))>prec || fabs(fabs(xC-xA)-fabs(xD-xB))>prec)) {
		cout<<" les points ("<<xA<<","<<yA<<"), ("<<xB<<","<<yB<<"), ("<<xC<<","<<yC<<"), ("<<xD<<","<<yD<<") ne forment pas un parallelogramme !!!\n";
		cout<<fabs(fabs(yB-yA)-fabs(yD-yC))<<" "<<fabs(fabs(xB-xA)-fabs(xD-xC))<<" "<<fabs(fabs(yC-yA)-fabs(yD-yB))<<" "<<fabs(fabs(xC-xA)-fabs(xD-xB))<<"\n";
	}
	double aire=-1., aireABC=fabs(xA*(yB-yC)-xB*(yA-yC)+xC*(yA-yB)), aireABD=fabs(xA*(yB-yD)-xB*(yA-yD)+xD*(yA-yB)), 
				 aireACD=fabs(xA*(yC-yD)-xC*(yA-yD)+xD*(yA-yC)), aireBCD=fabs(xB*(yC-yD)-xC*(yB-yD)+xD*(yB-yC)); 
	if (fabs(maxi(maxi(maxi(aireABC,aireABD),aireACD),aireBCD)-mini(mini(mini(aireABC,aireABD),aireACD),aireBCD))>prec) {
		cout<<" Pb les aires des parallelogrammes issus des triplets de points sont differentes : "<<aireABC<<" "<<aireABD<<" "<<aireACD<<" "<<aireBCD;
		cout<<"\n -> les 4 points ne forment pas un parallelogramme ?\n";
	} 
	else { 
		aire=(aireABC+aireABD+aireACD+aireBCD)/4.; 
//		cout<<" aire = "<<aire<<"\n"; 
	} 
	return aire;
}

//double aire_parallelogramme (double xA, double yA, double xB, double yB, double xC, double yC) {
//	return xA*(yB-yC)-xB*(yA-yC)+xC*(yA-yB);
//}

double tirage_gauss(const double v, const double m) { // méthode de Box-Muller (ici v est la variance et m la moyenne)
	double a, b, z;
  unsigned int    number;
  errno_t         err;
	do {
//		a=double(rand())/RAND_MAX;
		err=rand_s(&number);
    if (err!=0) cout<<"The rand_s function failed!\n";
    a=(double)number/(double)UINT_MAX;
	} while (a<=0 || a>1);
	do {
//		b=double(rand())/RAND_MAX;
		err=rand_s(&number);
    if (err!=0) cout<<"The rand_s function failed!\n";
    b=(double)number/(double)UINT_MAX;
	} while (b<0 || b>1);
	z=sqrt(-2.0*v*log(a))*cos(2.0*PI*b); //z=sqrt(-2.0*v*log(a))*sin(2.0*PI*b);
	return z+m;
}
/*double tirage_gauss(const double v, const double m) {
	double a, b, z;
	do {
		a=double(rand())/RAND_MAX;
	} while (a<=0 || a>1);
	do {
		b=double(rand())/RAND_MAX;
	} while (b<0 || b>1);
	z=sqrt(-2.0*v*log(a))*cos(2.0*PI*b);
	return z+m;
}*/

double tirage_gauss_rapide(const double v, const double m) { // méthode de Marsaglia et Bray (ici v est la variance et m la moyenne)
	double a, b, z, w;
  unsigned int    number;
  errno_t         err;
	do {
//		a=double(rand())/RAND_MAX;
		err=rand_s(&number);
    if (err!=0) cout<<"The rand_s function failed!\n";
    a=(double)number/(double)UINT_MAX;
	} while (a<=0 || a>1);
	do {
//		b=double(rand())/RAND_MAX;
		err=rand_s(&number);
    if (err!=0) cout<<"The rand_s function failed!\n";
    b=(double)number/(double)UINT_MAX;
	} while (b<0 || b>1);
	w=a*a+b*b;
	z=a*sqrt(2.*v/w*log(1./w)); //b=a*sqrt(2.*v/w*log(1./w));
	return z+m;
}
unsigned int tirage_proba(const double *P, const unsigned int d) {
	unsigned int i=0;
	double a=double(rand())/RAND_MAX, s=P[0];
	while (a>s && i<d-1) s+=P[++i];
//	if (i<0 || i>=d) cout<<" debordement de tableau "<<i<<">"<<d-1<<"\n";
	return i;
}
bool histo_norm(float *H, int nbin, bool iaf) {
	int i;
	double xx=0.;
	bool iOk=1;
	const float eps=(float)1.e-6;
	for (i=0; i<nbin; i++) {
		if (H[i]<-eps) {cout<<" valeur negative dans l'histogramme : H["<<i<<"]="<<H[i]<<" ???\n"; iOk=0;}
		if (H[i]<eps) H[i]=fabs(H[i]);
		xx+=H[i];
	}
	if (fabs(1-xx)>eps) {
		if (iaf) cout<<" normalisation de l'histogramme car avant normalisation somme = "<<xx<<"\n";
		if (xx>0) {for (i=0; i<nbin; i++) H[i]/=(float)xx;}
		else {
			if (xx==0) {
				cout<<" valeurs toutes a zero dans l'histogramme : somme = "<<xx<<" ???\n"; iOk=0;
				for (i=0; i<nbin; i++) H[i]=1.f/nbin;
			}
		}
	}
	return iOk;
}
double distBhattacharyya_histo(const float *H_1, const float *H_2, int nbin) {
	int i;
	double xx=0.;
	for (i=0; i<nbin; i++) xx+=pow(double(H_1[i]*H_2[i]),0.5);
	xx=1.-xx;
	if (xx>=0) return pow(xx,0.5);
	else {
		cout<<" PB: distance de Bhattacharyya negative ???\n";
		for (i=0; i<nbin; i++) cout<<setw(6)<<setprecision(2)<<H_1[i]; cout<<"\n";
		for (i=0; i<nbin; i++) cout<<setw(6)<<setprecision(2)<<H_2[i]; cout<<"\n"; {char aa; cin>>aa;}
		return -1.;
	}
}

/*
void tri_par_insertion(double *t, int n) {
	int i,j;
	if (n>1) tri_par_insertion(t,n-1);
	double x=t[n-1];
	i=0;
	while (t[i]<x) i++;
	for (j=n-1; j>i; j--) t[j]=t[j-1];
	t[i]=x;
}
void tri_par_selection(double *t, const int n) {
	int i,i_max=n-1;
	double xmax=t[n-1];
	for (i=0; i<n-1; i++)
		if (t[i]>xmax) {
			i_max=i;
			xmax=t[i];
		}
	t[i_max]=t[n-1];
	t[n-1]=xmax;
	if (n>2) tri_par_selection(t,n-1);
}
void ordre_decroiss_par_selection(double *t, int *r, const int n) {
	int i,i_min=n-1,rmin=r[i_min];
	double xmin=t[i_min];
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
void tri_bulle(double *t, const int n) {
	bool permut;
	int i;
	double x;
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
}*/
int partition(double *t, const int n) {
	int i, n_inf=0, pos_pivot=0; // ou bien : pos_pivot=rand()%n
	double x, pivot=t[pos_pivot];
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
void tri_rapide(double *t, const int n) {
	if (n>1) {
		int i,pos_pivot=partition(t,n);
		if (pos_pivot>1) {
			tri_rapide(t,pos_pivot);
		}
		if (pos_pivot<n-1) {
			int nbis=n-pos_pivot-1;
			double *tbis=new double[nbis];
			for (i=0; i<nbis; i++) tbis[i]=t[i+pos_pivot+1];
			tri_rapide(tbis,nbis);
			for (i=0; i<nbis; i++) t[i+pos_pivot+1]=tbis[i];
			if (tbis!=NULL) delete[] tbis;
		}
	}
}
int partition(double *t, const int debut, const int fin, const bool iaf) {
	int i, n_inf=debut, pos_pivot=n_inf; // ou bien : pos_pivot=rand()%n
	double x, pivot=t[n_inf];
	for (i=debut+1; i<=fin; i++)
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
	if (iaf) {
		cout<<"********* fin calcul de partition de "<<debut<<" a "<<fin<<": pivot de valeur "<<pivot<<" en "<<n_inf<<"\n";
		for (int j=debut; j<=fin; j++) cout<<j<<" "<<t[j]<<"\n";
	}
	return n_inf;
}
void tri_rapide_partiel(double *t, const int debut, const int fin, const bool iaf) {
	if (debut<fin) { if (iaf) cout<<" tri_rapide_partiel de "<<debut<<" a "<<fin<<"\n";
		int pos_pivot=partition(t,debut,fin,iaf);
		tri_rapide_partiel(t,debut,pos_pivot-1,iaf);
		int pos_pivot2=pos_pivot;
		while (pos_pivot2<fin && t[pos_pivot2]==t[pos_pivot2+1]) pos_pivot2++;
		tri_rapide_partiel(t,pos_pivot2+1,fin,iaf);
	}
}
void tri_rapide(double *t, const int n, const bool iaf) { // n est la taille du tableau
	int debut=0, fin=n-1, i;                    // fin est le dernier indice du tableau
	double xmin=1.e+39, xmax=-1.e+39;
	for (i=debut; i<=fin; i++) {
		if (t[i]<xmin) xmin=t[i];
		if (t[i]>xmax) xmax=t[i];
	}
	while (debut<fin && t[debut]<=xmin) debut++;
	while (debut<fin && t[fin]>=xmax) fin--;
	if (debut<fin) tri_rapide_partiel(t,debut,fin,iaf);
}
int partition(val_pos *t, const int debut, const int fin, const bool iaf) {
	int i, k, n_inf=debut, pos_pivot=n_inf; // ou bien : pos_pivot=rand()%n
	double x, pivot=t[n_inf].val;
	for (i=debut+1; i<=fin; i++)
		if (t[i].val<pivot) {
			x=t[i].val; k=t[i].pos;
			t[i].val=t[n_inf].val; t[i].pos=t[n_inf].pos;
			t[n_inf].val=x; t[n_inf].pos=k;
			if (pos_pivot==n_inf) pos_pivot=i;
			n_inf++;
/*	if (iaf) {
		cout<<" calcul partition en cours: pivot de valeur "<<pivot<<" en "<<n_inf<<"\n";
		for (int j=debut; j<=fin; j++) cout<<j<<" "<<t[j].val<<"\n";
	}*/
  }
	if (pos_pivot!=n_inf) {
		x=t[pos_pivot].val; k=t[pos_pivot].pos;
		t[pos_pivot].val=t[n_inf].val; t[pos_pivot].pos=t[n_inf].pos;
		t[n_inf].val=x; t[n_inf].pos=k;
	}
	if (iaf) {
		cout<<"********* fin calcul de partition de "<<debut<<" a "<<fin<<": pivot de valeur "<<pivot<<" en "<<n_inf<<"\n";
		for (int j=debut; j<=fin; j++) cout<<j<<" "<<t[j].val<<"\n";
	}
	return n_inf;
}
void tri_rapide_partiel(val_pos *t, const int debut, const int fin, const bool iaf) {
	if (debut<fin) { if (iaf) cout<<" tri_rapide_partiel de "<<debut<<" a "<<fin<<"\n";
		int pos_pivot=partition(t,debut,fin,iaf);
		tri_rapide_partiel(t,debut,pos_pivot-1,iaf);
		int pos_pivot2=pos_pivot;
		while (pos_pivot2<fin && t[pos_pivot2].val==t[pos_pivot2+1].val) pos_pivot2++;
		tri_rapide_partiel(t,pos_pivot2+1,fin,iaf);
	}
}
void tri_rapide(val_pos *t, const int n, const bool iaf) { // n est la taille du tableau
	int debut=0, fin=n-1, i;                    // fin est le dernier indice du tableau
	double xmin=1.e+39, xmax=-1.e+39;
	for (i=debut; i<=fin; i++) {
		if (t[i].val<xmin) xmin=t[i].val;
		if (t[i].val>xmax) xmax=t[i].val;
	}
//	while (debut<fin && t[debut].val==t[debut+1].val) debut++;
	while (debut<fin && t[debut].val<=xmin) debut++;
	while (debut<fin && t[fin].val>=xmax) fin--;
	if (debut<fin) tri_rapide_partiel(t,debut,fin,iaf);
}

float regres_lin (float **Txy, const int n, float &a, float &b) { //cout<<" debut regres_lin\n";
	double S_xi=0., S_yi=0., S_xiyi=0., S_xi2=0., x, y, delta;
	int i;
	for (i=0; i<n; i++) {
		x=Txy[i][0]; y=Txy[i][1];
		S_xi+=x; S_yi+=y; S_xi2+=x*x; S_xiyi+=x*y;
	}
	delta=n*S_xi2-S_xi*S_xi; //cout<<" delta = "<<delta<<"\n";
	float rmse=0.f;
	if (delta!=0) {
		a=(float)((n*S_xiyi-S_xi*S_yi)/delta);
		b=(float)((-S_xi*S_xiyi+S_xi2*S_yi)/delta);
		for (i=0; i<n; i++)
			rmse+=pow(Txy[i][1]-a*Txy[i][0]-b,2);
//		rmse/=n;
		rmse=pow(rmse/(float)n,0.5f);
	}
	else rmse=FLT_MAX;
//	cout<<" fin regres_lin\n";
	return rmse;
}


float algoKuhn (float **Tab, const int dim, bool iaf) { // la matrice doit être carrée
//	if (dim>=25) iaf=1;
	if (iaf) cout<<" algorithme de Kuhn\n";
	const double epsilon=1.e-10;
	const int itermax=100;
	int it=0,i,j,n,nmax,nz,ilig,jcol,k,nzeroindep;
	float x, xmin, xcout=0.f;
	bool fini=0;
	int *i_0=new int[dim*dim], *j_0=new int[dim*dim];
	bool *zerobarre=new bool[dim*dim], *zeromarqu=new bool[dim*dim], *lig_star=new bool[dim], *col_star=new bool[dim]; 
	int **M=new int*[dim]; for (i=0; i<dim; i++) M[i]=new int[dim];
	float *Tabs=new float[dim*dim];
	for (i=0; i<dim*dim; i++) Tabs[i]=Tab[i/dim][i%dim];
	if (iaf) aff_Tab2D(Tab,dim);
/*	if (iaf) {
		for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<setw(5)<<Tab[i][j]<<" "; cout<<"\n";} 
		cout<<"_______________________________\n";
	}*/
	do {
//		if (iaf) cout<<"\n\n\n iteration "<<it++<<"\n\n\n"; 
		if (dim>=25) cout<<"\n\n\n iteration "<<it++<<"\n\n\n"; //if (it>=10) iaf=1;
		for (i=0; i<dim*dim; i++) zerobarre[i]=zeromarqu[i]=0;
		for (i=0; i<dim; i++) {
			x=Tab[i][0];
			for (j=1; j<dim; j++) 
				if (Tab[i][j]<x) x=Tab[i][j];
			for (j=0; j<dim; j++) Tab[i][j]-=x;
		}
		for (i=0; i<dim; i++) {
			x=Tab[0][i];
			for (j=1; j<dim; j++) 
				if (Tab[j][i]<x) x=Tab[j][i];
			for (j=0; j<dim; j++) Tab[j][i]-=x;
		}
		if (iaf) aff_Tab2D(Tab,dim);
/*		if (iaf) {
			for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<setw(5)<<Tab[i][j]<<" "; cout<<"\n";} 
			cout<<"_______________________________\n";
		}*/
		nz=-1;
		for (i=0; i<dim; i++) {
			for (j=0; j<dim; j++) 
				if (Tab[i][j]==0) {nz++; i_0[nz]=i; j_0[nz]=j;}
		}
		nmax=nz;
		if (iaf) {cout<<nz+1<<" zeros dans matrice : "; for (k=0; k<=nz; k++) cout<<"("<<i_0[k]+1<<","<<j_0[k]+1<<") "; cout<<"\n";}
		int it2=0;
		while (nz>=0 && it2<100) {it2++;
			k=coord_zero_argmin_nbzeros_Kuhn(Tab,dim,zeromarqu,zerobarre,i_0,j_0,nmax);
			ilig=i_0[k]; jcol=j_0[k]; 
			zeromarqu[ilig*dim+jcol]=1; nz--; 
			nz=barrezeros_Kuhn (Tab,dim,zerobarre,zeromarqu,nz,ilig,jcol);
//			if (iaf) cout<<" reste "<<nz+1<<" zeros a traiter\n"; 
		}
		if (it2>=100) {cout<<"_______________________________ it2 = "<<it2<<"\n";
		for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<setw(5)<<Tab[i][j]<<" "; cout<<"\n";}
		cout<<"_______________________________\n";char aa; cin>>aa;}
		nz=nmax;
		nzeroindep=0; if (iaf) cout<<" zeros independents en ";
		for (n=0; n<=nz; n++) {
			ilig=i_0[n]; jcol=j_0[n]; 
			if (zeromarqu[ilig*dim+jcol]) {nzeroindep++; if (iaf) cout<<"("<<ilig+1<<","<<jcol+1<<") ";}
		} if (iaf) cout<<" => nombre de zeros independents : "<<nzeroindep<<"\n"; //{char aa; cin>>aa;}
		if (nzeroindep<dim) { 
			for (i=0; i<dim; i++) {lig_star[i]=1; col_star[i]=0;}
			for (n=0; n<=nz; n++) {
				ilig=i_0[n]; jcol=j_0[n]; 
				if (zeromarqu[ilig*dim+jcol]) lig_star[ilig]=0;
			}
			for (n=0; n<=nz; n++) {
				ilig=i_0[n]; jcol=j_0[n]; 
				if (zerobarre[ilig*dim+jcol] && lig_star[ilig]) col_star[jcol]=1;
			}
			for (n=0; n<=nz; n++) {
				ilig=i_0[n]; jcol=j_0[n]; 
				if (zeromarqu[ilig*dim+jcol] && col_star[jcol]) lig_star[ilig]=1;
			}
			if (iaf) {cout<<"  lignes etoiles : "; for (i=0; i<dim; i++) if (lig_star[i]) cout<<i+1<<" "; cout<<"\n";}
			if (iaf) {cout<<" colones etoiles : "; for (i=0; i<dim; i++) if (col_star[i]) cout<<i+1<<" "; cout<<"\n";}
			for (i=0; i<dim; i++) for (j=0; j<dim; j++) M[i][j]=0;
			for (i=0; i<dim; i++) 
				if (!lig_star[i]) for (j=0; j<dim; j++) M[i][j]=1;
			for (i=0; i<dim; i++) 
				if (col_star[i]) for (j=0; j<dim; j++) M[j][i]+=1;
			if (iaf) aff_Tab2D(M,dim);
/*			if (iaf) {
				cout<<"_______________________________\n";
				for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<M[i][j]<<" "; cout<<"\n";}
				cout<<"_______________________________\n";}*/
			xmin=FLT_MAX;
			for (i=0; i<dim; i++) 
				for (j=0; j<dim; j++) 
					if (M[i][j]==0 && xmin>Tab[i][j]) xmin=Tab[i][j];
			if (iaf) cout<<" xmin = "<<xmin<<"\n";
			if (xmin<=0) {
				cout<<"\n\n\n Pb : xmin nul !!!\n\n\n\n";
				for (i=0; i<dim; i++) 
					for (j=0; j<dim; j++) 
						if (M[i][j]==0 && fabs(Tab[i][j])<=epsilon) {if (iaf) cout<<" on etoile la ligne "<<i+1<<"\n";
							for (k=0; k<dim; k++) {M[i][k]+=1;}}
//				for (i=0; i<dim; i++) 
//					for (j=0; j<dim; j++) M[i][j]=M[j][i]=maxi(M[i][j],M[j][i]);
				if (iaf) aff_Tab2D(M,dim);
				xmin=FLT_MAX;
				for (i=0; i<dim; i++) 
					for (j=0; j<dim; j++) 
						if (M[i][j]==0 && xmin>Tab[i][j]) xmin=Tab[i][j];
				if (iaf) cout<<" xmin = "<<xmin<<"\n"; //{char aa; cin>>aa;}
				if (xmin<=0) {fini=1; xcout=999.;}
			}
			for (i=0; i<dim; i++) 
				for (j=0; j<dim; j++) Tab[i][j]+=(M[i][j]-1)*xmin;
			if (iaf) aff_Tab2D(Tab,dim);
/*			if (iaf) {
				cout<<"_______________________________\n";
				for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<setw(5)<<Tab[i][j]<<" "; cout<<"\n";}
				cout<<"_______________________________\n";}*/
		} else {
			for (n=0; n<=nz; n++) {
				ilig=i_0[n]; jcol=j_0[n]; k=ilig*dim+jcol;
				if (zeromarqu[k]) {
					xcout+=Tabs[k]; 
					for (i=0; i<dim; i++) Tab[ilig][i]=0; 
					Tab[ilig][jcol]=1;
				}
			} 
			fini=1;
		}
	} while (!fini && it<itermax);
	if (!fini || it>=itermax) {
		cout<<"_______________________________ it = "<<it<<"\n";
		for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<setw(5)<<Tab[i][j]<<" "; cout<<"\n";}
		cout<<"_______________________________\n";char aa; cin>>aa;
	}
	if (Tabs!=NULL) delete[] Tabs;
	for (i=0; i<dim; i++) if (M[i]!=NULL) delete[] M[i]; if (M!=NULL) delete[] M; 
	if (lig_star!=NULL) delete[] lig_star; if (col_star!=NULL) delete[] col_star;
	if (zerobarre!=NULL) delete[] zerobarre; if (zeromarqu!=NULL) delete[] zeromarqu;
	if (i_0!=NULL) delete[] i_0; if (j_0!=NULL) delete[] j_0; 
	return xcout;
}

float algoKuhn2 (float **Tab, const int dim, bool iaff) { // la matrice doit être carrée
	if (iaff) cout<<" algorithme de Kuhn modifie pour matrice symetrique\n";
	if (dim==1) {float x=Tab[0][0]; Tab[0][0]=1; return x;}
	const double epsilon=1.e-10;
	int it=0,itmax=dim*dim+2,i,j,n,nmax,nz,ilig,jcol,k,nzeroindep,nmaxzeroindep=0,i2,j2;
	float xmin, xcout=0.f, cout_adhoc;
	bool fini=0, i_adhoc=0, ifini2=0;
	int *i_0=new int[dim*dim], *j_0=new int[dim*dim];
	bool *zerobarre=new bool[dim*dim], *zeromarqu=new bool[dim*dim], *lig_star=new bool[dim], *col_star=new bool[dim]; 
	int **M=new int*[dim]; for (i=0; i<dim; i++) M[i]=new int[dim];
	float *Tabs=new float[dim*dim], **Tab_adhoc=new float*[dim];  for (i=0; i<dim; i++) Tab_adhoc[i]=new float[dim];
	for (i=0; i<dim*dim; i++) Tabs[i]=Tab[i/dim][i%dim];
	do {
		++it; if (iaff) cout<<"\n\n\n iteration "<<it<<"\n\n\n"; 
		for (i=0; i<dim*dim; i++) zerobarre[i]=zeromarqu[i]=0;
		soustrait_min_col_lig(Tab,dim);
		if (iaff) aff_Tab2D(Tab,dim);
		nz=-1;
		for (i=0; i<dim; i++) {
			for (j=0; j<dim; j++) 
				if (fabs(Tab[i][j])<=epsilon) {nz++; i_0[nz]=i; j_0[nz]=j;}
		}
		nmax=nz;
		if (iaff) {cout<<nz+1<<" zeros dans matrice : "; for (k=0; k<=nz; k++) cout<<"("<<i_0[k]+1<<","<<j_0[k]+1<<") "; cout<<"\n";}
		if (dim%2==1) {
			jcol=ilig=coord_zerodiag_argmin_nbzeros(Tab,dim);
			zeromarqu[ilig*dim+jcol]=1; if (iaff) cout<<"premier 0 marque en ("<<ilig+1<<","<<jcol+1<<")\n"; 
			if (Tab[ilig][jcol]==0) nz--;
			else {i_0[++nmax]=ilig; j_0[nmax]=jcol;}
			nz=barrezeros_Kuhn2 (Tab,dim,zerobarre,zeromarqu,nz,ilig,jcol);
		}
		while (nz>=0) {
			k=coord_zero_argmin_nbzeros_Kuhn2(Tab,dim,zeromarqu,zerobarre,i_0,j_0,nmax);
			ilig=i_0[k]; jcol=j_0[k]; 
			zeromarqu[ilig*dim+jcol]=1; nz--; 
			if (ilig!=jcol) {zeromarqu[jcol*dim+ilig]=1; nz--; }
			nz=barrezeros_Kuhn2 (Tab,dim,zerobarre,zeromarqu,nz,ilig,jcol);
//			if (iaff) cout<<" reste "<<nz+1<<" zeros a traiter\n"; 
		}
		nz=nmax;
		nzeroindep=0; if (iaff) cout<<" zeros independents en ";
		for (n=0; n<=nz; n++) {
			ilig=i_0[n]; jcol=j_0[n]; 
			if (zeromarqu[ilig*dim+jcol]) {nzeroindep++; if (iaff) cout<<"("<<ilig+1<<","<<jcol+1<<") ";}
		} if (iaff) cout<<" => nombre de zeros independents : "<<nzeroindep<<"\n"; //{char aa; cin>>aa;}
		if (nmaxzeroindep<nzeroindep) nmaxzeroindep=nzeroindep;
		if (nzeroindep==dim) fini=1;
		if (!fini) {
			for (i=0; i<dim; i++) {lig_star[i]=1; col_star[i]=0;}
			for (n=0; n<=nz; n++) {
				ilig=i_0[n]; jcol=j_0[n]; 
				if (zeromarqu[ilig*dim+jcol]) lig_star[ilig]=0;
			}
			if (iaff) {cout<<"  lignes etoiles : "; for (i=0; i<dim; i++) if (lig_star[i]) cout<<i+1<<" "; cout<<"\n";}
			for (n=0; n<=nz; n++) {
				ilig=i_0[n]; jcol=j_0[n]; 
				if (zerobarre[ilig*dim+jcol] && lig_star[ilig]) col_star[jcol]=1;
			}
			if (iaff) {cout<<" colones etoiles : "; for (i=0; i<dim; i++) if (col_star[i]) cout<<i+1<<" "; cout<<"\n";}
			if (i_adhoc && !ifini2) {cout<<" on passe en ad hoc\n";
				int dimred=0;
				dimred=dim-nzeroindep;
				float **Tabred=new float*[dimred];
				for (j=0; j<dimred; j++) Tabred[j]=new float[dimred]; i2=j2=0;
				for (i=0; i<dim; i++) {
					if (lig_star[i]) {j2=0; 
						for (j=0; j<dim; j++) {
							if (!col_star[j] && lig_star[j]) Tabred[i2][j2++]=Tab[i][j]; 
						} i2++;
					}
				}
				if (iaff) {cout<<" matrice reduite de dimension "<<dimred<<"\n"; aff_Tab2D(Tabred,dimred); }
				float coutred=0; Tabred[0][0]=1;
				if (dimred>1) coutred=algoKuhn2(Tabred,dimred);
				if (iaff) cout<<" cout reduit "<<coutred<<"\n"; 
				coutred=0;
				bool *Tab_adhocs=new bool[dim*dim];
				for (i=0; i<dim*dim; i++) Tab_adhocs[i]=0; xcout=coutred;
				for (n=0; n<=nz; n++) {
					ilig=i_0[n]; jcol=j_0[n]; k=ilig*dim+jcol;
					if (zeromarqu[k]) Tab_adhocs[k]=Tab_adhocs[jcol*dim+ilig]=1;
				} i2=0;
				for (i=0; i<dim; i++) {
					if (lig_star[i]) {j2=0; 
						for (j=0; j<dim; j++) {
							if (!col_star[j] && lig_star[j]) {if (Tabred[i2][j2]==1) Tab_adhocs[i*dim+j]=Tab_adhocs[j*dim+i]=1; j2++; }
						} i2++;
					}
				}
				for (i=0; i<dim; i++)
					for (j=0; j<dim; j++)
						if (Tab_adhocs[i*dim+j]==1) coutred+=Tabs[i*dim+j];
				cout<<" solution adhoc trouvee : cout = "<<coutred<<"\n";
				if (iaff) aff_Tabpseudo2D(Tab_adhocs,dim);
				if (coutred<cout_adhoc) { cout<<" on actualise la solution ad-hoc\n";
					cout_adhoc=coutred;
					for (i=0; i<dim; i++)
						for (j=0; j<dim; j++)
							Tab_adhoc[i][j]=Tab_adhocs[i*dim+j];
				} else {cout<<" on n'actualise pas la solution ad-hoc car "<<coutred<<">"<<cout_adhoc<<"\n";}
				for (j=0; j<dimred; j++) if (Tabred[j]!=NULL) delete[] Tabred[j]; 
				if (Tabred!=NULL) delete[] Tabred;
				if (iaff) aff_Tab2D(Tab_adhoc,dim);
				for (i=0; i<dim*dim; i++) Tab_adhocs[i]=0; 
				xmin=coutred; i2=j2=-1;
				int nz2=dim*dim-1;
				for (i=0; i<dim; i++) {
					for (j=0; j<dim; j++) {
						if (fabs(Tab[i][j])>epsilon) {
							nz2--;
							if (Tab[i][j]<xmin) {xmin=Tab[i][j]; i2=i; j2=j;} 
						}
					}
				}
				if (iaff) cout<<" xmin "<<xmin<<" en ("<<i2<<","<<j2<<")\n";
				bool *zerobarre2=new bool[dim*dim], *zeromarqu2=new bool[dim*dim];
				for (i=0; i<dim*dim; i++) zerobarre2[i]=zeromarqu2[i]=0;
				zeromarqu2[i2*dim+j2]=zeromarqu2[j2*dim+i2]=1;
				Tab_adhocs[i2*dim+j2]=Tab_adhocs[j2*dim+i2]=1;
				if (iaff) cout<<" initialement nz2 = "<<nz2+1<<"\n";
				nz2=barrezeros_Kuhn2 (Tab,dim,zerobarre2,zeromarqu2,nz2,i2,j2);
//				if (iaff) cout<<" reste "<<nz2+1<<" zeros a traiter\n";
				while (nz2>=0) {
					k=coord_zero_argmin_nbzeros_Kuhn2(Tab,dim,zeromarqu2,zerobarre2,i_0,j_0,nmax);
					ilig=i_0[k]; jcol=j_0[k]; 
					zeromarqu2[ilig*dim+jcol]=1; nz2--; Tab_adhocs[ilig*dim+jcol]=Tab_adhocs[jcol*dim+ilig]=1; 
					if (ilig!=jcol) {zeromarqu2[jcol*dim+ilig]=1; nz2--;}
					nz2=barrezeros_Kuhn2 (Tab,dim,zerobarre2,zeromarqu2,nz2,ilig,jcol);
//					if (iaff) cout<<" reste "<<nz2+1<<" zeros a traiter\n"; 
				}
				nz2=nmax;
				int nzeroindep2=0;
				for (i=0; i<dim*dim; i++) nzeroindep2+=Tab_adhocs[i];
				if (nzeroindep2==dim) {
					coutred=0;
					for (i=0; i<dim; i++)
						for (j=0; j<dim; j++)
							if (Tab_adhocs[i*dim+j]==1) coutred+=Tabs[i*dim+j];
					if (coutred<cout_adhoc) {
						cout_adhoc=coutred; 
						ifini2=1;
						cout<<" Solution optimale retrouvee : cout = "<<coutred<<"\n"; if (iaff) aff_Tabpseudo2D(Tab_adhocs,dim);
						for (i=0; i<dim*dim; i++) Tab_adhoc[i/dim][i%dim]=Tab_adhocs[i];
					}
				} else {
					cout<<" pas de solution optimale trouvee\n";
				}
				if (zerobarre2!=NULL) delete[] zerobarre2; if (zeromarqu2!=NULL) delete[] zeromarqu2; 
				if (Tab_adhocs!=NULL) delete[] Tab_adhocs;
			}
			if (!ifini2 && !fini && it<itmax) {
				for (n=0; n<=nz; n++) {
					ilig=i_0[n]; jcol=j_0[n]; 
					if (zeromarqu[ilig*dim+jcol] && col_star[jcol] && fabs(Tab[ilig][jcol])<=epsilon) lig_star[ilig]=1;
				}
				if (iaff) {cout<<"  lignes etoiles : "; for (i=0; i<dim; i++) if (lig_star[i]) cout<<i+1<<" "; cout<<"\n";}
				for (i=0; i<dim; i++) for (j=0; j<dim; j++) M[i][j]=0;
				for (i=0; i<dim; i++) 
					if (!lig_star[i]) for (j=0; j<dim; j++) M[i][j]=1;
				for (i=0; i<dim; i++) 
					if (col_star[i]) for (j=0; j<dim; j++) M[j][i]+=1;
				if (iaff) aff_Tab2D(M,dim);
				for (i=0; i<dim; i++) 
					if (!lig_star[i] && !col_star[i]) for (j=0; j<dim; j++) M[j][i]+=1;
				for (i=0; i<dim; i++) 
					if (col_star[i] && lig_star[i]) for (j=0; j<dim; j++) M[i][j]+=1;
				if (iaff) aff_Tab2D(M,dim);
				xmin=FLT_MAX;
				for (i=0; i<dim; i++) 
					for (j=0; j<dim; j++) 
						if (M[i][j]==0 && xmin>Tab[i][j]) xmin=Tab[i][j];
				if (iaff) cout<<" xmin = "<<xmin<<"\n";
				if (xmin<=0) {
					for (i=0; i<dim; i++) 
						for (j=0; j<dim; j++) 
							if (M[i][j]==0 && fabs(Tab[i][j])<=epsilon) {if (iaff) cout<<" on etoile la ligne "<<i+1<<"\n";
								for (k=0; k<dim; k++) {M[i][k]+=1;}}
					for (i=0; i<dim; i++) 
						for (j=0; j<dim; j++) M[i][j]=M[j][i]=maxi(M[i][j],M[j][i]);
					if (iaff) aff_Tab2D(M,dim);
					xmin=FLT_MAX;
					for (i=0; i<dim; i++) 
						for (j=0; j<dim; j++) 
							if (M[i][j]==0 && xmin>Tab[i][j]) xmin=Tab[i][j];
					if (iaff) cout<<" xmin = "<<xmin<<"\n"; //{char aa; cin>>aa;}
					if (xmin<=0) {fini=1; xcout=999.;}
				}
				for (i=0; i<dim; i++) 
					for (j=0; j<dim; j++) Tab[i][j]+=(M[i][j]-1)*xmin;
				if (iaff) aff_Tab2D(Tab,dim);
			}
		} else {
			for (i=0; i<dim; i++) for (j=0; j<dim; j++) Tab[i][j]=0; 
			for (n=0; n<=nz; n++) {
				ilig=i_0[n]; jcol=j_0[n]; k=ilig*dim+jcol;
				if (zeromarqu[k]) {
					xcout+=Tabs[k]; 
					Tab[ilig][jcol]=Tab[jcol][ilig]=1;
				}
			}
		}
		if (it==dim*dim) {
			i_adhoc=1; cout_adhoc=FLT_MAX; cout<<" on n'a pas converge !!!\n"; //iaff=1; 
//			char aa; cin>>aa;
		}
	} while (!fini && it<itmax && !ifini2);
	if (i_adhoc) { cout<<" on recupere une sol. plus ou moins ad hoc\n";
		for (i=0; i<dim; i++) for (j=0; j<dim; j++) Tab[i][j]=Tab_adhoc[i][j];
		xcout=cout_adhoc;
		if (!ifini2) {xcout+=1000;}
//		{char aa; cin>>aa;}
	}
	if (iaff) cout<<"\n\n\n nombre d'iterations = "<<it<<"\n\n\n"; 
	if (Tabs!=NULL) delete[] Tabs;
	for (i=0; i<dim; i++) if (Tab_adhoc[i]!=NULL) delete[] Tab_adhoc[i]; if (Tab_adhoc!=NULL) delete[] Tab_adhoc; 
	for (i=0; i<dim; i++) if (M[i]!=NULL) delete[] M[i]; if (M!=NULL) delete[] M; 
	if (lig_star!=NULL) delete[] lig_star; if (col_star!=NULL) delete[] col_star;
	if (zerobarre!=NULL) delete[] zerobarre; if (zeromarqu!=NULL) delete[] zeromarqu;
	if (i_0!=NULL) delete[] i_0; if (j_0!=NULL) delete[] j_0; 
	return xcout;
}

















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































