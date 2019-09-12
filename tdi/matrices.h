#ifndef _MATRICES_H
#define _MATRICES_H

#include <iostream>
#include <iomanip>
#include <complex>
using namespace std;

/*****************************************************************************************************************/
/*********************************************** classe matrice2D ************************************************/
/*****************************************************************************************************************/

template <class T> class matrice2D
{	int nl, nc;
	T **M;
public:
	matrice2D(int n_l=1, int n_c=1) : nl(n_l), nc(n_c) { // constructeur par defaut
		M=new T*[nl];
		int i, j;
		for (i=0; i<nl; i++) M[i]=new T[nc];
		for (i=0; i<nl; i++) 
			for (j=0; j<nc; j++) M[i][j]=0;
	};
	matrice2D(int n_l, int n_c, T** tab) : nl(n_l), nc(n_c) { // constructeur par copie d'un tableau 2D
		M=new T*[nl];
		int i, j;
		for (i=0; i<nl; i++) M[i]=new T[nc];
		for (i=0; i<nl; i++) 
			for (j=0; j<nc; j++) M[i][j]=tab[i][j];
	};
	matrice2D(matrice2D<T>, matrice2D<T>, int, int); // autres constructeurs
	matrice2D(matrice2D<T>, matrice2D<T>, int, int, int);

	matrice2D(const matrice2D<T> &A) : nl(A.nl), nc(A.nc) { // constructeur de copie
		M=new T*[nl];
		int i, j;
		for (i=0; i<nl; i++) M[i]=new T[nc];
		for (i=0; i<nl; i++) 
			for (j=0; j<nc; j++) M[i][j]=A.M[i][j];
	};
	matrice2D<T>& operator= (const matrice2D<T> &A) { // operateur d'affectation
		if (this!=&A) {
			int i, j;
			for (i=0; i<nl; i++) 
				if (M[i]!=NULL) delete[] M[i];
			if (M!=NULL) delete M;
			nl=A.nl; nc=A.nc;
			M=new T*[nl];
			for (i=0; i<nl; i++) M[i]=new T[nc];
			for (i=0; i<nl; i++) 
				for (j=0; j<nc; j++) M[i][j]=A.M[i][j];
		}
		return (*this);
	};
	~matrice2D() { // destructeur
		for (int i=0; i<nl; i++) 
			if (M[i]!=NULL) {delete[] M[i]; M[i]=NULL;}
		if (M!=NULL) {delete M; M=NULL;}
	};
	int nlig() {return nl;};
	int ncol() {return nc;};
	T determinant2x2 () const {
		T detM_2x2=0;
		detM_2x2=M[0][0]*M[1][1]-M[0][1]*M[1][0];
		return detM_2x2;
	}
	T determinant3x3 () const {
		T detM_3x3=0;
		detM_3x3=M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[1][0]*(M[0][1]*M[2][2]-M[0][2]*M[2][1]);
		detM_3x3+=M[2][0]*(M[1][2]*M[0][1]-M[1][1]*M[0][2]);
		return detM_3x3;
	}
	bool valeurspropres2x2 (double T_vp[2]) const {
		const double eps=1.e-9;
		double sum_vp=M[0][0]+M[1][1], prod_vp=determinant2x2(), det=sum_vp*sum_vp-4*prod_vp;
		if (fabs(det-eps)>0) {
			T_vp[0]=(sum_vp-pow(det,0.5))/2.; T_vp[1]=(sum_vp+pow(det,0.5))/2.;
			return true;
		}
		else return false;
	}
	matrice2D<T> operator + (matrice2D<T> &);
	matrice2D<T> operator - (matrice2D<T> &);
	matrice2D<T> operator + (float);
	matrice2D<T> operator * (float);

	operator matrice2D<double> () ; // conversion de type
	operator matrice2D<float> ();
	operator matrice2D<int> ();

	T& operator () (int i, int j) {
		if (i<0 || i>=nl || j<0 || j>=nc) {
			cout<<"matrice (nblig ="<<nl<<", nbcol ="<<nc<<") : debordement en ("<<i<<","<<j<<")\n";
			if (i<0) i=0; if (i>=nl) i=nl-1;
			if (j<0) j=0; if (j>=nc) j=nc-1;
		}
		T *adval=&M[i][j];
		return *adval;
	};
	void affiche (ostream& os=cout) {
		int i,j;
		for (i=0; i<nl; i++) {
			for (j=0; j<nc; j++) os<<setw(8)<<setprecision(2)<<(*this)(i,j)<<" ";
			os<<"\n";
		}
	};
	bool operator ==(matrice2D<T>); // test d'égalité; cet opérateur affiche la norme sup de la difference
	bool operator !=(matrice2D<T>); // test de non-égalité; cet opérateur n'affiche rien
	matrice2D<complex<double>> transpoC(); // calcul de la matrice hermitienne (transposee conjuguee)
	matrice2D<T> transpo();  // transposition
	void perm_lig(int,int);  // permutation des lignes
	void perm_col(int,int);  // permutation des colonnes
//	matrice2D<double> L_inv (bool&); // inversion matrice triangulaire inferieure
	matrice2D<T> L_inv(bool&); // inversion matrice triangulaire inferieure
//	matrice2D<double> U_inv (bool&); // inversion matrice triangulaire superieure
	matrice2D<T> U_inv(bool&); // inversion matrice triangulaire superieure
	matrice2D<double> LU_decomp(matrice2D<double>&, bool&); // decomposition LU : PA=LU avec L et U triangulaires
	matrice2D<complex<double>> LU_decomp(matrice2D<complex<double>>&, bool&); // decomposition LU : PA=LU avec L et U triangulaires
	matrice2D<double> inverse (bool&);   // inversion matrice
	bool inverse (matrice2D<double>&);   // inversion matrice
	bool inverse(matrice2D<complex<double>>&); // inversion matrice complexe
// analyse 'statistique' de la matrice : application = analyse de texture d'après la matrice de cooccurrence
	double entropy(double=2.);
	double dissimilarity();
	double contrast();
	double homogeneity();
};

template <class T> matrice2D<T>::matrice2D(matrice2D<T> A, matrice2D<T> B, int n_l, int n_c) : nl(n_l), nc(n_c) {
	M=new T*[nl];
	int i,j;
	for (i=0; i<nl; i++) M[i]=new T[nc];
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) M[i][j]=A(i,j)+B(i,j);
}

template <class T> matrice2D<T>::matrice2D(matrice2D<T> A, matrice2D<T> B, int n_l, int n, int n_c) : nl(n_l), nc(n_c) {
	M=new T*[nl];
	int i,j,k;
	T x;
	for (i=0; i<nl; i++) M[i]=new T[nc];
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) {
			x=0;
			for (k=0; k<n; k++) x+=A(i,k)*B(k,j);
			M[i][j]=x;
		}
}

template <class T> matrice2D<T> matrice2D<T>::operator + (matrice2D<T> &A) {
	int nlR=0, ncR=0, i,j;
	if (nl==A.nl && nc==A.nc) {nlR=nl; ncR=nc;}
	else cout<<" Pb : les dim. des matrices ne sont pas les memes\n";
	matrice2D<T> R(nlR,ncR);
	for (i=0; i<nlR; i++) 
		for (j=0; j<ncR; j++) R.M[i][j]=M[i][j]+A.M[i][j];
	return R;
}

template <class T> matrice2D<T> matrice2D<T>::operator - (matrice2D<T> &A) {
	int nlR=0, ncR=0, i,j;
	if ((nl==A.nl) && (nc==A.nc)) {nlR=nl; ncR=nc;}
	else cout<<" Pb : les dim. des matrices ne sont pas les memes\n";
	matrice2D<T> R(nlR,ncR);
	for (i=0; i<nlR; i++) 
		for (j=0; j<ncR; j++) R.M[i][j]=M[i][j]-A.M[i][j];
	return R;
}

template <class T> matrice2D<T> matrice2D<T>::operator * (float u) {
	matrice2D<T> R(nl,nc);
	int i,j;
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) R.M[i][j]=M[i][j]*u;
	return R;
}

template <class T> matrice2D<T> matrice2D<T>::operator + (float u) {
	matrice2D<T> R(nl,nc);
	int i,j;
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) R.M[i][j]=M[i][j]+u;
	return R;
}
template <class T> matrice2D<T>::operator matrice2D<double> () {
	int i,j;
	matrice2D<double> A(nl,nc);
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) A(i,j)=(double)(*this)(i,j);
	return A;
}
template <class T> matrice2D<T>::operator matrice2D<float> () {
	int i,j;
	matrice2D<float> A(nl,nc);
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) A(i,j)=(float)(*this)(i,j);
	return A;
}
template <class T> matrice2D<T>::operator matrice2D<int> () {
	int i,j;
	matrice2D<int> A(nl,nc);
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) A(i,j)=(int)(*this)(i,j);
	return A;
}

template <> void matrice2D<complex<double>>::affiche (ostream&); // specialisation pour le cas complex<double> cf. matrices.cpp

template <> bool matrice2D<complex<double>>::operator == (matrice2D<complex<double>>); // specialisation pour le cas complex<double> cf. matrices.cpp

template <class T> bool matrice2D<T>::operator == (matrice2D<T> A) { // cet operateur affiche la norme sup de la difference
	const double eps=1.e-6;
	bool ieq=1;
	double nsup=0., temp;
	if (nl!=A.nl || nc!=A.nc) ieq=0;
	else {
		int i,j,ii=-1,jj=-1;
		for (i=0; i<nl; i++)
			for (j=0; j<nc; j++) {
				temp=fabs((*this)(i,j)-A(i,j));
				if (temp>nsup) {nsup=temp; ii=i; jj=j;}
			}
		if (nsup>eps) {ieq=0; cout<<" norme sup de la difference = "<<nsup<<" en ("<<ii<<","<<jj<<")\n";}
	}
	return ieq;
}

template <> bool matrice2D<complex<double>>::operator != (matrice2D<complex<double>>); // specialisation pour le cas complex<double> cf. matrices.cpp

template <class T> bool matrice2D<T>::operator != (matrice2D<T> A) { // cet operateur n'affiche rien 
	const double eps=1.e-6;
	bool ieq=1;
	double nsup=0., temp;
	if (nl!=A.nl || nc!=A.nc) ieq=0;
	else {
		int i,j,ii,jj;
		for (i=0; i<nl; i++)
			for (j=0; j<nc; j++) {
				temp=fabs((*this)(i,j)-A(i,j));
//				cout<<(*this)(i,j)<<" "<<A(i,j)<<" "<<temp<<" "<<ieq<<"\n";
				/*if (temp<0) ieq=0;
				else*/
					if (temp>nsup) {nsup=temp; ii=i; jj=j;}
			}
		if (nsup>eps) ieq=0;
//		cout<<" ieq = "<<ieq<<"\n";
	}
	return !ieq;
}

template <class T> matrice2D<T> matrice2D<T>::transpo () {
	int i,j;
	matrice2D<T> R(nc,nl);
	for (i=0; i<nl; i++) 
		for (j=0; j<nc; j++) R.M[j][i]=M[i][j];
	return R;
}

template <class T> void matrice2D<T>::perm_lig (int i1, int i2) {	
	int j;
	T temp;
	if (i1>=0 && i1<nl && i2>=0 && i2<nl) {
		for (j=0; j<nc; j++) {
			temp=M[i1][j];
			M[i1][j]=M[i2][j];
			M[i2][j]=temp;
		}
	}
	else
		cout<<" PB ds matrice2D<T>::perm_lig (int i1, int i2) "<<i1<<" "<<i2<<"\n";
}

template <class T> void matrice2D<T>::perm_col (int j1, int j2) {
	int i;
	T temp;
	if (j1>=0 && j1<nc && j2>=0 && j2<nc) {
		for (i=0; i<nl; i++) {
			temp=M[i][j1];
			M[i][j1]=M[i][j2];
			M[i][j2]=temp;
		}
	}
	else
		cout<<" PB ds matrice2D<T>::perm_col (int j1, int j2) "<<j1<<" "<<j2<<"\n";
}

template <> matrice2D<complex<double>> matrice2D<complex<double>>::L_inv (bool&); // specialisation pour le cas complex<double> cf. matrices.cpp

//template <class T> matrice2D<double> matrice2D<T>::L_inv (bool &OK) {
template <class T> matrice2D<T> matrice2D<T>::L_inv (bool &OK) {
	int n=nl,i,j,k;
	matrice2D<double> L_1(n,n);
	if (nl!=nc) {
		cout<<" inversion d'une matrice triangulaire inferieure non carree ????\n";
		OK=0;
	}
	if (OK) {
//		cout<<"\n inversion de la matrice triangulaire inferieure : \n"; affiche();
		double x;
		for (i=0;i<n;i++) L_1(i,i)=1;
		for (j=0;j<n;j++) {     // boucle sur les colonnes
			for (i=j;i<n;i++) { // boucle sur les lignes au dessous de la diagonale	(L_1 est aussi triang inf)
				x=L_1(i,j); 
				for (k=0; k<i; k++) x-=(*this)(i,k)*L_1(k,j);
				if ((*this)(i,i)!=0) L_1(i,j)=x/(*this)(i,i); 
				else OK=0;
			}
		}
//		cout<<" produit L.L_1 = \n"; (matrice2D<double>((matrice2D<double>)(*this),L_1,n,n,n)).affiche();
//		cout<<" produit L_1.L = \n"; (matrice2D<double>(L_1,(matrice2D<double>)(*this),n,n,n)).affiche();
	}
//	return L_1;
	return matrice2D<T>(L_1);
}

template <> matrice2D<complex<double>> matrice2D<complex<double>>::U_inv (bool&); // specialisation pour le cas complex<double> cf. matrices.cpp

//template <class T> matrice2D<double> matrice2D<T>::U_inv (bool &OK) {
template <class T> matrice2D<T> matrice2D<T>::U_inv (bool &OK) {
	int n=nl,i,j,k;
	matrice2D<double> U_1(n,n);
	if (nl!=nc) {
		cout<<" inversion d'une matrice triangulaire superieure non carree ????\n";
		OK=0;
	}
	if (OK) {
//		cout<<"\n inversion de la matrice triangulaire superieure : \n"; affiche();
		double x;
		for (i=0;i<n;i++) U_1(i,i)=1;
		for (j=n-1;j>=0;j--) {   // boucle sur les colonnes
			for (i=j;i>=0;i--) { // boucle sur les lignes au dessus de la diagonale (U_1 est aussi triang sup)
				x=U_1(i,j); 
				for (k=i+1; k<n; k++) x-=(*this)(i,k)*U_1(k,j); //cout<<k<<" x-som = "<<x<<"\n";}
				if ((*this)(i,i)!=0) U_1(i,j)=x/(*this)(i,i);
				else OK=0;
			}
		}
//		cout<<" produit U.U_1 = \n"; (matrice2D<double>((matrice2D<double>)(*this),U_1,n,n,n)).affiche();
//		cout<<" produit U_1.U = \n"; (matrice2D<double>(U_1,(matrice2D<double>)(*this),n,n,n)).affiche();
	}
//	return U_1;
	return matrice2D<T>(U_1);
}

template <> matrice2D<complex<double>> matrice2D<complex<double>>::LU_decomp(matrice2D<complex<double>>&, bool&); // specialisation pour le cas complex<double> cf. matrices.cpp

template <class T> matrice2D<double> matrice2D<T>::LU_decomp(matrice2D<double> &P, bool &OK) {
	OK=1;
	const double eps=1.e-6;
	matrice2D<double> A(*this), L(nl,nc);
	if (nl!=nc) {
		cout<<" decomposition LU d'une matrice non carree ????\n";
		OK=0;
	}
	if (OK) {
		int n=nl,i,j,lpiv;
		double big,temp;
		matrice2D<double> U((*this)), Id(n,n);
		for (i=0;i<n;i++) Id(i,i)=1.f; // matrice identite
		for (j=0; j<n; j++) {          // boucle sur les colonnes
			big=0.0;
			for (i=0; i<n; i++)        // recherche dans la colonne de l'element le plus grand (pivot)
				if ((temp=fabs((*this)(i,j)))>big) big=temp;
			if (big<eps) {OK=0; /*cout<<"Matrice singuliere\n";*/}
		}
		if (OK) {
			L=Id; P=Id;                // P est la matrice des permutations
			matrice2D<double> Ln(n,n);
			for (j=0;j<n;j++) {        // boucle sur les colonnes
				Ln=Id;
				big=U(j,j); lpiv=j;
				for (i=j+1; i<n; i++)
					if ((temp=fabs(U(i,j)))>big) {big=temp; lpiv=i;}
				if (lpiv!=j) {
//					cout<<" permutation ligne "<<j<<" et "<<lpiv<<"\n";//U.affiche(); 
					U.perm_lig(j,lpiv);  
					P.perm_lig(j,lpiv);  
					L.perm_col(j,lpiv);
				}
				if (n!=nl || n!=nc) cout<<" µµµµµµµµµµµµµµµµµ PB : n!=nl || n!=nc : "<<n<<" "<<nl<<" "<<nc<<"\n";
				if (U(j,j)!=0)
					for (i=j+1; i<n; i++) Ln(i,j)=-U(i,j)/U(j,j); 
				U=matrice2D<double>(Ln,U,n,n,n); 
				L=matrice2D<double>(L,Ln.L_inv(OK),n,n,n);
			}
//			cout<<" matrice L de la decomposition LU :\n"; L.affiche();
//			cout<<" matrice U de la decomposition LU :\n"; U.affiche();
//			cout<<" matrice P de permutation : P.A=L.U\n"; P.affiche();
//			cout<<" produit LU :\n"; (matrice2D<double>(L,U,n,n,n)).affiche();
			if (matrice2D<double>(L,U,n,n,n)!=A) {OK=0; cout<<" Pb ds decomposition LU\n";}
			(*this)=U;
			L=matrice2D<double>(P,L,n,n,n); // sinon L contient P'L
		}
	}
	return L;
}

template <class T> matrice2D<double> matrice2D<T>::inverse(bool &OK) {
	int n=nl, k;
	matrice2D<double> L(nl,nc), U(*this), L_1, U_1, A(*this);
	matrice2D<double> Id(n,n); for (k=0; k<n; k++) Id(k,k)=1;
	OK=1;
	if (nl!=nc) {
		cout<<" inversion d'une matrice non carree ????\n";
		OK=0;
	}
	if (OK) {
		matrice2D<double> P(Id);
		L=U.LU_decomp(P,OK); /*cout<<" matrice L \n"; L.affiche(); cout<<" matrice U \n"; U.affiche();
		cout<<" verification de la decomposition L.U ";
		(matrice2D<double>(L,U,n,n,n)==matrice2D<double>(P,A,n,n,n))?cout<<"OK":cout<<"Not OK";cout<<"\n";
		(matrice2D<double>(L,U,n,n,n)-matrice2D<double>(P,A,n,n,n)).affiche();*/
		if (OK) {
			U_1=U.U_inv(OK);
			/*cout<<" verification de l'inversion de U ";
			(matrice2D<double>(U,U_1,n,n,n)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
			(matrice2D<double>(U,U_1,n,n,n)).affiche();*/
		}
		if (OK) {
			L_1=L.L_inv(OK);
			/*cout<<" verification de l'inversion de L ";
			(matrice2D<double>(L,L_1,n,n,n)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
			(matrice2D<double>(L,L_1,n,n,n)).affiche();*/
		}
		if (OK) {
			L_1=matrice2D<double>(L_1,P,n,n,n);       // PA=LU <=> Ainv=(Uinv)(Linv)P
			/*cout<<" verification de l'inversion de A ";
			(matrice2D<double>(matrice2D<double>(U_1,L_1,n,n,n),A,n,n,n)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
			(matrice2D<double>(matrice2D<double>(U_1,L_1,n,n,n),A,n,n,n)).affiche();*/
		}
	}
/*	(matrice2D<double>(matrice2D<double>(U_1,L_1,n,n,n),A,n,n,n)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
	(matrice2D<double>(matrice2D<double>(U_1,L_1,n,n,n),A,n,n,n)).affiche();*/
	return matrice2D<double>(U_1,L_1,n,n,n);
}

template <class T> bool matrice2D<T>::inverse(matrice2D<double> &Ainv) {
	int n=nl, k;
	matrice2D<double> L(nl,nc), U(*this), L_1, U_1, A(*this);
	matrice2D<double> Id(n,n); for (k=0; k<n; k++) Id(k,k)=1;
	bool inversionOK=1;
	if (nl!=nc) inversionOK=0;
	if (inversionOK) {
		matrice2D<double> P(Id);
		L=U.LU_decomp(P,inversionOK);
		if (inversionOK) {
//			cout<<" matrice L \n"; L.affiche(); cout<<" matrice U \n"; U.affiche(); cout<<" matrice P \n"; P.affiche(); 
			if (matrice2D<double>(L,U,n,n,n)!=matrice2D<double>(P,A,n,n,n)) {
				inversionOK=0;
//				cout<<" L\n"; L.affiche(); cout<<" U\n"; U.affiche();
//				cout<<" P\n"; P.affiche(); cout<<" A\n"; A.affiche();
			}
			else {
//				cout<<" matrice L \n"; L.affiche(); cout<<" matrice U \n"; U.affiche();
/*				cout<<" verification de la decomposition L.U ";
				(matrice2D<double>(L,U,n,n,n)==matrice2D<double>(P,A,n,n,n))?cout<<"OK":cout<<"Not OK";cout<<"\n";
				(matrice2D<double>(L,U,n,n,n)-matrice2D<double>(P,A,n,n,n)).affiche();*/
				U_1=U.U_inv(inversionOK); //cout<<" matrice U_1 OK ? "<<inversionOK<<"\n"; U_1.affiche();  
				if (matrice2D<double>(U,U_1,n,n,n)!=Id) inversionOK=0;
				else {
/*					cout<<" verification de l'inversion de U ";
					(matrice2D<double>(U,U_1,n,n,n)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
					(matrice2D<double>(U,U_1,n,n,n)).affiche();*/
					L_1=L.L_inv(inversionOK); //cout<<" matrice L_1 OK ? "<<inversionOK<<"\n"; L_1.affiche();
					if (matrice2D<double>(L,L_1,n,n,n)!=Id) inversionOK=0;
					else {
/*						cout<<" verification de l'inversion de L ";
						(matrice2D<double>(L,L_1,n,n,n)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
						(matrice2D<double>(L,L_1,n,n,n)).affiche();*/
						L_1=matrice2D<double>(L_1,P,n,n,n);       // PA=LU <=> Ainv=(Uinv)(Linv)P
//						cout<<" matrice L_1.P\n"; L_1.affiche();
/*						cout<<" verification de l'inversion de A ";
						(matrice2D<double>(matrice2D<double>(U_1,L_1,n,n,n),A,n,n,n)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
						(matrice2D<double>(matrice2D<double>(U_1,L_1,n,n,n),A,n,n,n)).affiche();*/
					}
				}
			}
		}
	}
	if (inversionOK) Ainv=matrice2D<double>(U_1,L_1,n,n,n);
	return inversionOK;
}

template <class T> double matrice2D<T>::entropy(double b) {
	double e=0.;
  for (int i=0; i<nl; i++) 
		for (int j=0; j<nc; j++) 
			if (M[i][j]>0) e-=M[i][j]*log(M[i][j]);
  return e/log(b);
}

template <class T> double matrice2D<T>::contrast() {
	double d=0.;
  for (int i=0; i<nl; i++) 
		for (int j=0; j<nc; j++) d+=pow(double(i-j),2)*M[i][j];
  return d;
}

template <class T> double matrice2D<T>::dissimilarity() {
	double d=0.;
  for (int i=0; i<nl; i++) 
		for (int j=0; j<nc; j++) d+=abs(i-j)*M[i][j];
  return d;
}

template <class T> double matrice2D<T>::homogeneity() {
	double h=0.;
  for (int i=0; i<nl; i++) 
		for (int j=0; j<nc; j++) h+=1./(1.+pow((double)i-j,2))*M[i][j];
  return h;
}


/*****************************************************************************************************************/
/*********************************************** classe matrice1D ************************************************/
/*****************************************************************************************************************/


template <class T> class matrice1D
{	int n;
	T *M;
public:
	matrice1D(int _n=1) : n(_n) { // constructeur par defaut
		M=new T[n];
		for (int i=0; i<n; i++) M[i]=0;
	};
	matrice1D(int _n, T* tab) : n(_n) { // constructeur par copie d'un tableau 1D
		M=new T[n];
		for (int i=0; i<n; i++) M[i]=tab[i];
	};
	matrice1D(const matrice1D<T> &A) : n(A.n) { // constructeur de copie
		M=new T[n];
		for (int i=0; i<n; i++) M[i]=A.M[i];
	};
	matrice1D<T>& operator= (const matrice1D<T> &A) { // operateur d'affectation
		if (this!=&A) {
			if (M!=NULL) delete M;
			n=A.n;
			M=new T[n];
			for (int i=0; i<n; i++) M[i]=A.M[i];
		}
		return (*this);
	};
	~matrice1D() { // destructeur
		if (M!=NULL) {delete M; M=NULL;}
	};

	int dim() {return n;};

	matrice1D<T> operator + (matrice1D<T> &);
	matrice1D<T> operator - (matrice1D<T> &);
	matrice1D<T> operator * (float);

	operator matrice1D<double> () ; // conversion de type
	operator matrice1D<float> ();
	operator matrice1D<int> ();

	T& operator () (int i) {
		if (i<0 || i>=n) {
			cout<<" dim ="<<n<<" => debordement d'indice en "<<i<<"\n";
			if (i<0) i=0; if (i>=n) i=n-1;
			}
		T *adval=&M[i];
		return *adval;
	};
	void affiche (ostream &os=cout) {
		for (int i=0; i<n; i++) os<<setw(8)<<setprecision(2)<<(*this)(i)<<" ";
	};

	bool operator == (matrice1D<T> A) {
		const double eps=1.e-6;
		bool ieq=1;
		double nsup=0., temp;
		if (n!=A.n) ieq=0;
		else {
			int i,ii;
			for (i=0; i<n; i++)
				if ((temp=fabs((*this)(i)-A(i)))>nsup) {nsup=temp; ii=i;}
				if (nsup>eps) {ieq=0; cout<<" norme sup de la difference = "<<nsup<<" en ("<<ii<<")\n";}
		}
		return ieq;
	}
// analyse 'statistique' de la matrice : application = analyse de texture à l'ordre 1
	double moment (int=1);	
	double mean () {moment();}
	double contrast (double=255.);
	double entropy (double=2.);
	double energy ();
};

template <class T> matrice1D<T> matrice1D<T>::operator + (matrice1D<T> &A) {
	int nlR=0, i;
	if (n==A.n) {nR=n;}
	else cout<<" Pb : la dim. des matrices n'est pas la meme\n";
	matrice1D<T> R(nR);
	for (i=0; i<nR; i++) R.M[i]=M[i]+A.M[i];
	return R;
}

template <class T> matrice1D<T> matrice1D<T>::operator - (matrice1D<T> &A) {
	int nlR=0, i;
	if (n==A.n) {nR=n;}
	else cout<<" Pb : la dim. des matrices n'est pas la meme\n";
	matrice1D<T> R(nR);
	for (i=0; i<nR; i++) R.M[i]=M[i]-A.M[i];
	return R;
}

template <class T> matrice1D<T> matrice1D<T>::operator * (float u) {
	matrice1D<T> R(n);
	for (int i=0; i<n; i++) R.M[i]=M[i]*u;
	return R;
}

template <class T> double matrice1D<T>::moment(int k) {
	int i;
	double muk=0., mu1=0., s=0.;
	for (i=0; i<n; i++) {mu1+=i*M[i]; s+=M[i];}
	if (s>0) mu1/=s;
	if (k>1) {
		for (i=0; i<n; i++) muk+=pow(i-mu1,k)*M[i];
		if (s>0) muk/=s;
	}
	else muk=mu1;
 	return muk;
}

template <class T> double matrice1D<T>::contrast(double ampl) {
	double maxp=0, minp=DBL_MAX, crst;
    int  i=0;
	while (M[i]==0 && i<n-1) i++;
	minp=i;
	i=n-1;
	while (M[i]==0 && i>0) i--;
	maxp=i;
	crst= ((double)(maxp-minp))/(maxp+minp);
	return crst*ampl;
}

template <class T> double matrice1D<T>::entropy(double b) {
	double e=0.;
    for (int i=0; i<n; i++) 
		if (M[i]>0) e-=M[i]*log(M[i]); 
 	return e/log(b);
            
}

template <class T> double matrice1D<T>::energy() {
	double e=0.;
    for (int i=0; i<n; i++) e+=pow(M[i],2); 
 	return e;
}

// fonctions d'algebre de Numerical Recipes
bool svdcmp(float**, int, int, float*, float**);

#endif