#ifndef _IMADATA_H
#define _IMADATA_H

#include <float.h>
#include "def.h"
#include "constantes.h"
#include "imasites.h"
#include "pixel.h"
#include "matrices.h"
#include "fonctions.h"
#include "fourier.h"

struct eltHisto {
	float nbelts;
	float bornesup;
};

template <class T> class sampleset;

template <class T> class imadata : public imasites
{protected:
	int nbcanaux;
	T *Tima;
	T valnul;
	T *valmin, *valmax, *valmod;
	double *valmoy, *valvar;
	double *CarClas;
	eltHisto *histo;
	int nbinH;
	bool istat, ihist;
 public:
	imadata(int nl=1, int nc=1, int nd=1) : imasites(nl,nc) {
		Tima=valmin=valmax=valmod=NULL; valmoy=valvar=CarClas=NULL; histo=NULL;
		long int i; int j;
		nbcanaux=nd;
		if (typeid(valnul)==typeid(int) || typeid(valnul)==typeid(BYTE) || typeid(valnul)==typeid(bool) ||
				typeid(valnul)==typeid(unsigned int) || typeid(valnul)==typeid(unsigned long int) ) valnul=0;
		else valnul=(T)(-999);
		Tima=new T[nbpix*nbcanaux];
		for (i=0; i<nbpix; i++) {
			for (j=0; j<nbcanaux; j++) Tima[i*nbcanaux+j]=valnul;
			Tsites[i]=&(Tima[i*nbcanaux]);
		}
		valmin=new T[nbcanaux];
		valmax=new T[nbcanaux];
		valmod=new T[nbcanaux];
		valmoy=new double[nbcanaux];
		valvar=new double[nbcanaux];
		for (j=0; j<nbcanaux; j++) {
			valmin[j]=valmax[j]=valmod[j]=valnul; valmoy[j]=valvar[j]=0.;}
		nbinH=256;
		histo=new eltHisto[nbinH*nbcanaux];
		istat=0; ihist=0;
	}
	imadata(matrice2D<T> &A) : imasites(A.nlig(),A.ncol()) {
		Tima=valmin=valmax=valmod=NULL; valmoy=valvar=CarClas=NULL; histo=NULL;
		long int i;
		nbcanaux=1;
		if (typeid(valnul)==typeid(int) || typeid(valnul)==typeid(BYTE) || typeid(valnul)==typeid(bool) ||
				typeid(valnul)==typeid(unsigned int) || typeid(valnul)==typeid(unsigned long int) ) valnul=0;
		else valnul=(T)(-999);
//		cout<<" valnul = "<<valnul<<"\n";
		Tima=new T[nbpix*nbcanaux];
		for (i=0; i<nbpix; i++) {
			Tima[i]=A(i/nbcol,i%nbcol);
			Tsites[i]=&(Tima[i]);
		}
		valmin=new T[nbcanaux]; valmax=new T[nbcanaux]; valmod=new T[nbcanaux]; valmoy=new double[nbcanaux]; valvar=new double[nbcanaux];
		valmin[0]=valmax[0]=valmod[0]=(T)valnul; valmoy[0]=valvar[0]=(double)valnul;
		nbinH=256;
		histo=new eltHisto[nbinH*nbcanaux];
		istat=0; ihist=0;
	}
	~imadata() {
		if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
//		if (istat && valmin!=NULL) {delete[] valmin; valmin=NULL;}
//		if (istat && valmax!=NULL) {delete[] valmax; valmax=NULL;}
//		if (ihist && valmod!=NULL) {delete[] valmod; valmod=NULL;}
//		if (istat && valmoy!=NULL) {delete[] valmoy; valmoy=NULL;}
//		if (istat && valvar!=NULL) {delete[] valvar; valvar=NULL;}
//		if (ihist && histo!=NULL)  {delete[] histo; histo=NULL;}
		if (valmin!=NULL) {delete[] valmin; valmin=NULL;}
		if (valmax!=NULL) {delete[] valmax; valmax=NULL;}
		if (valmod!=NULL) {delete[] valmod; valmod=NULL;}
		if (valmoy!=NULL) {delete[] valmoy; valmoy=NULL;}
		if (valvar!=NULL) {delete[] valvar; valvar=NULL;}
		if (histo!=NULL)  {delete[] histo; histo=NULL;}
	}
	imadata(imabin &imaB) : imasites(imaB.nlig(),imaB.ncol()) {
		Tima=valmin=valmax=valmod=NULL; valmoy=valvar=CarClas=NULL; histo=NULL;
		long int i;
		nbcanaux=1;
		if (typeid(valnul)==typeid(int) || typeid(valnul)==typeid(BYTE) || typeid(valnul)==typeid(bool) ||
				typeid(valnul)==typeid(unsigned int) || typeid(valnul)==typeid(unsigned long int) ) valnul=0;
		else valnul=(T)(-999);
		Tima=new T[nbpix*nbcanaux];
		for (i=0; i<nbpix; i++) {
			Tima[i]=(T)(imaB(i/nbcol,i%nbcol)?255:0);
			Tsites[i]=&(Tima[i]);
		}
		valmin=new T[nbcanaux]; valmax=new T[nbcanaux]; valmod=new T[nbcanaux]; valmoy=new double[nbcanaux]; valvar=new double[nbcanaux];
		valmin[0]=valmax[0]=valmod[0]=(T)valnul; valmoy[0]=valvar[0]=(double)valnul;
		nbinH=256;
		histo=new eltHisto[nbinH*nbcanaux];
		istat=0; ihist=0;
	}
	imadata(const imadata<T> &ima) : imasites(ima) {
		Tima=valmin=valmax=valmod=NULL; valmoy=valvar=CarClas=NULL; histo=NULL;
		nbcanaux=ima.nbcanaux;
		Tima=new T[nbpix*nbcanaux];
		T *adval;
		long int i;
		for (i=0; i<nbpix; i++) {
			adval=(T*)ima.Tsites[i];
			for (int j=0; j<nbcanaux; j++) Tima[i*nbcanaux+j]=*(adval+j);
			Tsites[i]=&(Tima[i*nbcanaux]);
		}
		valnul=ima.valnul;
		valmin=new T[nbcanaux]; valmax=new T[nbcanaux]; valmod=new T[nbcanaux];
		valmoy=new double[nbcanaux]; valvar=new double[nbcanaux];
		for (i=0; i<nbcanaux; i++) {
			valmin[i]=ima.valmin[i]; valmax[i]=ima.valmax[i]; valmod[i]=ima.valmod[i];
			valmoy[i]=ima.valmoy[i]; valvar[i]=ima.valvar[i];
		}
		nbinH=ima.nbinH;
		histo=new eltHisto[nbinH*nbcanaux];
		for (i=0; i<nbinH*nbcanaux; i++) histo[i]=ima.histo[i];
		istat=ima.istat; ihist=ima.ihist;
	}
	imadata& operator=(const imadata<T> &ima) {
		if (this != &ima) {
			imasites *ad1, *ad2;
			ad1=this;
			ad2=(imasites*)&ima;
			*ad1=*ad2;
			if (Tima!=NULL) delete[] Tima;
			nbcanaux=ima.nbcanaux;
			nbpix=ima.nbpix;
			Tima=new T[nbpix*nbcanaux];
			T *adval;
			long int i;
			for (i=0; i<nbpix; i++) {
				adval=(T*)ima.Tsites[i];
				for (int j=0; j<nbcanaux; j++) Tima[i*nbcanaux+j]=*(adval+j);
				Tsites[i]=&(Tima[i*nbcanaux]);
			}
			valnul=ima.valnul;
			if (valmin!=NULL) delete[] valmin; if (valmax!=NULL) delete[] valmax;
			if (valmod!=NULL) delete[] valmod;
			if (valmoy!=NULL) delete[] valmoy; if (valvar!=NULL) delete[] valvar;
			valmin=new T[nbcanaux]; valmax=new T[nbcanaux]; valmod=new T[nbcanaux];
			valmoy=new double[nbcanaux]; valvar=new double[nbcanaux];
			for (i=0; i<nbcanaux; i++) {
				valmin[i]=ima.valmin[i]; valmax[i]=ima.valmax[i]; valmod[i]=ima.valmod[i];
				valmoy[i]=ima.valmoy[i]; valvar[i]=ima.valvar[i];
			}
			if (histo!=NULL) delete[] histo; nbinH=ima.nbinH;
			histo=new eltHisto[nbinH*nbcanaux];
			for (i=0; i<nbinH*nbcanaux; i++) histo[i]=ima.histo[i];
			istat=ima.istat; ihist=ima.ihist;
		}
		return *this;
	}
	int ncanaux() const {
		return nbcanaux;
	}
  T v_nulle() const {
		return valnul;
	}
	T histov(int i, int j=0) const {
		if (histo!=NULL && i>=0 && i<nbinH*nbcanaux) {
			if (j==0) return histo[i].nbelts;
			else return histo[i].bornesup;
		}
		else return valnul;
	}
	T& operator () (int i, int j, int k=0) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol || k<0 || k>=nbcanaux) {
			cout<<" nblig ="<<nblig<<", nbcol ="<<nbcol<<", nbcanaux ="<<nbcanaux<<")\n";
			cout<<" debordement d'indice dans ("<<i<<","<<j<<","<<k<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
 			if (j>=nbcol) j=nbcol-1;
			if (k<0) k=0;
			if (k>=nbcanaux) k=nbcanaux-1;
		}
		T *adval=(T*)Tsites[i*nbcol+j];
		adval=adval+k;
		return *adval;
	}
	T& operator () (int i, int j, int k=0) const {
		if (i<0 || i>=nblig || j<0 || j>=nbcol || k<0 || k>=nbcanaux) {
			cout<<" nblig ="<<nblig<<", nbcol ="<<nbcol<<", nbcanaux ="<<nbcanaux<<")\n";
			cout<<" debordement d'indice dans ("<<i<<","<<j<<","<<k<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
 			if (j>=nbcol) j=nbcol-1;
			if (k<0) k=0;
			if (k>=nbcanaux) k=nbcanaux-1;
		}
		T *adval=(T*)Tsites[i*nbcol+j];
		adval=adval+k;
		return *adval;
	}
	T valpix (int i, int j, int k=0) const {
		if (i<0 || i>=nblig || j<0 || j>=nbcol || k<0 || k>=nbcanaux) {
			cout<<" nblig ="<<nblig<<", nbcol ="<<nbcol<<", nbcanaux ="<<nbcanaux<<")\n";
			cout<<" debordement d'indice dans ("<<i<<","<<j<<","<<k<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
 			if (j>=nbcol) j=nbcol-1;
			if (k<0) k=0;
			if (k>=nbcanaux) k=nbcanaux-1;
		}
//		T x=Tima[(i*nbcol+j)*nbcanaux+k]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
		T x=*((T*)Tsites[i*nbcol+j]+k);
		return x;
	}
	pixel<T> pix (int i, int j) const {
		if (i<0 || i>=nblig || j<0 || j>=nbcol) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
			if (j>=nbcol) j=nbcol-1;
		}
		T *adval=(T*)Tsites[i*nbcol+j];
		pixel<T> pix(nbcanaux);
		for (int k=0; k<nbcanaux; k++) pix[k]=*(adval+k);
		return pix;
	}
	void pix (int i, int j, pixel<T> pix) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
			if (j>=nbcol) j=nbcol-1;
		}
		T *adval=(T*)Tsites[i*nbcol+j];
		for (int k=0; k<nbcanaux; k++) *(adval+k)=pix[k];
	}
	sampleset< pixel<T> > dataset () const {
		sampleset< pixel<T> > S;
		int i,j;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) S.ajoute(pix(i,j),i*nbcol+j);
		cout<<" nombre d'echantillons extraits de l'image = "<<S.nsample()<<"\n";
		return S;
	}
 	sampleset< pixel<T> > dataset (T val_out) const {
 		sampleset< pixel<T> > S;
 		int i,j,ii,k,n=0;
 		bool valid;
 		cout<<" valeur 'nulle' = "<<val_out<<"\n";
 		for (i=0; i<nblig; i++)
 			for (j=0; j<nbcol; j++) {
 				valid=1;
 				ii=i*nbcol+j;
 				for (k=0; k<nbcanaux; k++)
// 					if (Tima[ii*nbcanaux+k]==val_out) valid=0; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
 					if (*((T*)Tsites[ii]+k)==val_out) valid=0;
 				if (valid) S.ajoute(pix(i,j),n++);
 			}
 		cout<<" nombre d'echantillons extraits de l'image = "<<S.nsample()<<"\n";
 		return S;
 	}
 	sampleset< pixel<T> > dataset (imabin imask) const {
 		sampleset< pixel<T> > S;
 		int i,j,n=0;
 		for (i=0; i<nblig; i++)
 			for (j=0; j<nbcol; j++)
 				if (imask(i,j)) S.ajoute(pix(i,j),n++);
 		cout<<" nombre d'echantillons extraits de l'image = "<<S.nsample()<<"\n";
 		return S;
	}
	matrice2D<T> conv2mat2D () const {
		matrice2D<T> M(nblig,nbcol);
		int i, j;
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) M(i,j)=(*this)(i,j);
		return M;
	};
	imadata(const imadata<T> &ima, int icanal) : imasites(ima) {
		Tima=valmin=valmax=valmod=NULL; valmoy=valvar=CarClas=NULL; histo=NULL;
		nbcanaux=1;
		Tima=new T[nbpix];
		T *adval;
		long int i;
		for (i=0; i<nbpix; i++) {
			adval=(T*)ima.Tsites[i];
			Tima[i]=*(adval+icanal);
			Tsites[i]=&(Tima[i]);
		}
		valnul=ima.valnul;
		valmin=new T[nbcanaux]; valmax=new T[nbcanaux]; valmod=new T[nbcanaux];
		valmoy=new double[nbcanaux]; valvar=new double[nbcanaux];
		valmin[0]=ima.valmin[icanal]; valmax[0]=ima.valmax[icanal]; valmod[0]=ima.valmod[icanal];
		valmoy[0]=ima.valmoy[icanal]; valvar[0]=ima.valvar[icanal];
		istat=ima.istat;
		nbinH=256;
		histo=new eltHisto[nbinH*nbcanaux];
		ihist=0;
	}
	void copiecanal(int k, imadata<T> &ima, int icanal=0) {
		int i, j, imax=mini(nblig,ima.nblig), jmax=mini(nbcol,ima.nbcol);
		for (i=0; i<imax; i++)
			for (j=0; j<jmax; j++) (*this)(i,j,k)=ima(i,j,icanal);
		istat=0; ihist=0;
	}
	imadata<T> mult_canal (float, int=0) const;
	imadata<T> add_canal (float, int=0) const;
	operator imadata<float> ();
	operator imadata<int> ();
	operator imadata<BYTE> ();
	imadata<BYTE> imaunsignedchar (bool=0, bool=0);
	imadata<unsigned __int16> imaunsignedshortint (bool=0, bool=0);
	void affiche (int=1) const;
	void inverse_ngr (T,int=0);
	void quantif_ngr (int=8);
  imadata<T> borneval(double,double) const;
	void tracesegmentdroite(const int, const int, const int, const int, const T=255);
	imadata<T> operator + (imadata<T> &) const;
	imadata<T> operator - (imadata<T> &) const;
	imadata<T> operator * (imadata<T> &) const;
	imadata<T> operator / (imadata<T> &) const;
	imadata<T> operator * (float) const;
	imadata<T> operator + (float) const;
	imadata<T> absI () const;
	bool operator == (imadata<T> &);
	bool operator == (const imadata<T> &) const;
	void sauve_Ima (char *nomfich="imagesauvee.dat", int=0) const;
	void sauve_Ima (string nomfich="imagesauvee.dat", int=0) const;
	void sauve_ImaBSQ(char *nomfich="imaBSQsauvee.dat") const;
	void sauve_ImaBSQ(string nomfich="imaBSQsauvee.dat") const;
	void sauve_ImaBIP(char *nomfich="ima3d_sauvee.dat") const;
	void sauve_ImaPGM(string nomfich="ima_sauvee.pgm");
	void sauve_ImaPGM(char *nomfich="image_sauvee.pgm");
	void sauve_ImaPGM2(string nomfich="image_sauvee.pgm");
	void sauve_ImaPGM2(char *nomfich="image_sauvee.pgm");
// statistiques image
	void statbasic (bool=0);
	double minI (int=0);
	double maxI (int=0);
	double moyI (int=0);
	double moyI (int=0) const;
	double varI (int=0);
	double varI (int=0) const;
	double modeI (double &, int=0);
	float histogramme_allc (int=256, float=0.5, bool=0);
	void histogramme (int=256, bool=0, bool=1);
	float seuil_percentile (float, int=0);
	void histogramme (imabin &, int, bool, bool=1);
	float seuil_otsu (imabin &, int=0);
	void supprime_Niveaux_vide (bool=0);
	float seuil_percentile (imabin &, float, int=0);
	void egalisehisto (int=64, bool=0);
	void egalisehisto (imabin &, int, bool);
	imadata<T> ima_diffligne () const;
	imadata<T> ima_diffcolonne () const;
	void mise_a_zero ();
	float norm() const;
	void normalize_all (float=1.);
//	imadata<unsigned int> ima_Integrale_16bin(T=255, T=0);
	void ima_Integrale_16bin(imadata<unsigned int> &, T=255, T=0);
	void ima_Integrale_Nbin(imadata<unsigned int> &, T=255, T=0);
	int hist_from_imaIntegrale_16bin(int *,unsigned short,int,int,int,int,bool=0) const; // attention l'objet appelant doit être une image integrale
	int hist_from_imaIntegrale_Nbin(int *,unsigned short,int,int,int,int,bool=0) const; // attention l'objet appelant doit être une image integrale
	double distHisto_Bhattacharyya (imadata<T> &, int=64);
// bruitage
	void add_gauss_noise (float);
	void add_impul_noise (float);
// filtrages passe-bas
	imadata<T> BR_filtremoyenne (int=3, int=3) const;
	imadata<T> filtremoyenne (int=3, int=3) const;
	imadata<T> filtregaussienne (float=1.41, bool=0) const;
	imadata<T> filtregaussienne2D (float=1.41, bool=0) const;
	imadata<T> filtreNagao (int=3, bool=0, BYTE=7, bool=0) const;
	imadata<T> filtreSymNearNeigh (int=3) const;
	imadata<T> filtremedian (int=3, int=3) const;
	imadata<T> filtremedian1Dadapt (int=5, int=8) const;
	void psnr (const imadata<T> &) const;
	void filtre_sigma_delta(imadata<T> &, T delta, imadata<T> *imasigma=NULL);
	void filtre_sigma_delta(imadata<T> &, imabin, T delta, imadata<T> *imasigma=NULL);
	void filtre_sigma_delta_codebook(imadata<T> &, imadata<T> &, imadata<BYTE> &, T delta);
// transformations spectrales et geometriques
	float boite_englobee(int&, int&, int&, int&) const;
	imadata<float> filtrage_Hamming () const;
	imadata<float> transformee_Fourier(bool=0) const;
	imadata<float> transformee_Fourier_Inv(bool=0) const;
	imadata<float> correlation_phase(const imadata<T> &, int&, int&, float&, bool=0, bool=0) const;
	imadata<float> transformee_FourierMellin(const imadata<T> &, float &, float &, float &, float &, bool=0) const;
	imadata<T> projette_ima_geom_fixe(const double, const double, const double, const double, const int, const int, unsigned char=0) const;
	imadata<T> projette_ima(double=0., double=0., double=1., double=0., bool=0, unsigned char=0) const;
	imadata<T> projette_ima(matrice2D<double>&, int, int, imabin&, unsigned char=1, unsigned char=0) const;
	imadata<T> projette_ima(matrice2D<double>&, unsigned char=1) const;
// morphologie mathematique
//	imadata<T> dilate (const eltstruct);
	imadata<T> dilate (const eltstruct, int=1);
//	imadata<T> erode (const eltstruct);
	imadata<T> erode (const eltstruct, int=1);
	imadata<T> rehausse_contraste (const eltstruct, int=1, float=0.5, float=0.5);
	imadata<T> gradient_m (const eltstruct);
	imadata<T> laplacien_m (const eltstruct);
	imadata<T> ouverture (const eltstruct);
	imadata<T> ouverture (const eltstruct, int);
	imadata<T> fermeture (const eltstruct);
	imadata<T> fermeture (const eltstruct, int);
	imadata<T> tophat (const eltstruct);
	imadata<T> tophat_c (const eltstruct);
	imadata<T> filtrealterne (const eltstruct, int, bool=1);
	imadata<T> reconst_geod (const imadata<T> &, const eltstruct);
	imadata<T> marq_centres_lowerset(const eltstruct,int=0);
	imadata<T> h_max(const eltstruct, float);
	imadata<T> marq_max_regionaux(const eltstruct, float);
// seuillages
	imadata<BYTE> seuil_ima (float) const;
	void seuil_ima (float, imadata<BYTE> &, int=0) const;
	void seuil_hysteresis (float, float, imabin&, int=0);
	imadata<BYTE> seuil_hysteresis (float, float, int=0) const;
	void seuil_hysteresis (float, float, imadata<BYTE> &, int=0);
	imadata<float> init1_segmentationregions (const int, const float, int&);
// image du gradient/laplacien
	imadata<float> gradient (imadata<float> &, char* ="Sobel") const;
	imadata<float> laplacien (char* ="8connex") const;
	imadata<float> gradientOpt (imadata<float> &, float, char* ="Deriche") const;
	imadata<float> steerablefilter (imadata<float> &, float) const;
	imadata<BYTE> chemin_opt4connex (imadata<float> &, int, int, int, int, float=-1.f) const;
	imadata<BYTE> chemin_opt4connex (imadata<float> &, int, int, const imabin &, int &, int &, float=-1.f) const;
// texture
	imadata<float> ima_contrast (int,int=0,bool=0);
	imadata<float> ima_energy (int,int=0,bool=0);
	imadata<float> ima_entropy (int,int=0,bool=0);
	imadata<float> ima_moment (int,int,int=0);
	imadata<float> cooccurence (int,int,int=0,bool=0);
	imadata<float> ima_carac_cooccurence (int,int,int,unsigned char,int=0,bool=0);
	imadata<float> ima1k_LBP (int=0) const;
	imadata<BYTE> ima1k_lenghtLBP(int=0) const;
	imadata<BYTE> ima1k_lenghtLTP(int=0) const;
};

template <class T> void imadata<T>::affiche(int ip) const {
	cout<<" affichage des valeurs de l'image "<<nblig<<" "<<nbcol<<" "<<nbcanaux<<"\n";
	for (int k=0; k<nbcanaux; k++) {
		cout<<" canal "<<k<<"\n";
		for (int i=0; i<nblig; i+=ip) {
			for (int j=0; j<nbcol; j+=ip) cout<<setw(4)<<(*this)(i,j,k);
			cout<<"\n";
		}
		cout<<"\n";
	}
}

template <class T> void imadata<T>::inverse_ngr(T pivot, int k) {
	float x, xmin, xmax;
	if (valnul==0) {
		xmin=0; xmax=255;
		for (int i=0; i<nblig; i++)
			for (int j=0; j<nbcol; j++) {
				x=pivot-(*this)(i,j,k);
				if (x<xmin) x=xmin;
				if (x>xmax) x=xmax;
				(*this)(i,j,k)=x;
			}
	}
	else { // on suppose qu'il n'y aura pas de débordement
		for (int i=0; i<nblig; i++)
			for (int j=0; j<nbcol; j++)
				(*this)(i,j,k)=pivot-(*this)(i,j,k);
	}
	statbasic();
}

template <class T> void imadata<T>::quantif_ngr(int n) {
	int i,j,k;
	float xmin, xmax, xdelta;
	statbasic();
	for (k=0; k<nbcanaux; k++) {
		xmin=valmin[k];
		xmax=valmax[k];
		xdelta=(xmax-xmin)/n;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) (*this)(i,j,k)=(T)(maxi(0,around((float)((*this)(i,j,k)-xmin)/xdelta-0.5f)));
	}
	statbasic();
}

template <class T> imadata<T> imadata<T>::borneval(double xmin,double xmax) const {
	int i,j,k;
	imadata<T> imaRes(*this);
	for (k=0; k<nbcanaux; k++) {
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaRes(i,j,k)=(T)(maxi(xmin,mini(xmax,(double)imaRes(i,j,k))));
	}
	imaRes.statbasic();
	return imaRes;
}

template<class T> void imadata<T>::tracesegmentdroite(const int i1, const int j1, const int i2, const int j2, const T val) {
	int i,j;
	if (i2==i1) {
		for (j=mini(j1,j2); j<=maxi(j1,j2); j++) (*this)(i1,j)=val;
	} else {
		float a=(float)(j2-j1)/(float)(i2-i1), b=j2-a*i2;
		for (i=mini(i1,i2); i<=maxi(i1,i2); i++) if (a*i+b>=0 && a*i+b<nbcol) (*this)(i,(int)(a*i+b))=val;
		if (a!=0) for (j=mini(j1,j2); j<=maxi(j1,j2); j++) if ((j-b)/a>=0 && (j-b)/a<nblig) (*this)((int)((j-b)/a),j)=val;
	}
}

template <class T> imadata<T> imadata<T>::mult_canal (float x, int k) const {
	int i,j;
	imadata<T> imaRes(*this);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (imaRes(i,j,k)!=valnul) imaRes(i,j,k)*=x;
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T> imadata<T>::add_canal (float x, int k) const {
	int i,j;
	imadata<T> imaRes(*this);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (imaRes(i,j,k)!=valnul) imaRes(i,j,k)+=x;
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T>::operator imadata<float>() {
	int i,j,k;
	imadata<float> ima_f(nblig,nbcol,nbcanaux);
	float x;
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			T *adval=(T*)Tsites[i*nbcol+j];
			for (k=0; k<nbcanaux; k++) {
				x=(float)(*(adval+k));
				ima_f(i,j,k)=x;
			}
			ima_f[i*nbcol+j]=&(ima_f(i,j,0));
		}
	}
	return ima_f;
}

template <class T> imadata<T>::operator imadata<int>() {
	int i,j,k,x;
	T val;
	imadata<int> ima_u(nblig,nbcol,nbcanaux);
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			T *adval=(T*)Tsites[i*nbcol+j];
			for (k=0; k<nbcanaux; k++) {
				val=*(adval+k);
				if (val==valnul) x=0;
				else x=int(val);
				ima_u(i,j,k)=x;
			}
		}
	}
	return ima_u;
}

template <class T> imadata<T>::operator imadata<BYTE>() { //cout<<" operator imadata<BYTE>\n"; 
	int i,j,k,x;
	T val;
	imadata<BYTE> ima_u(nblig,nbcol,nbcanaux);
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			T *adval=(T*)Tsites[i*nbcol+j];
			for (k=0; k<nbcanaux; k++) {
				val=*(adval+k);
				if (val==valnul) x=0;
				else x=int(val);
				x=mini(maxi(0,x),255);
				ima_u(i,j,k)=(BYTE)x;
			}
//			ima_u[i*nbcol+j]=&(ima_u(i,j,0));
		}
	}
	return ima_u;
}

template <class T> imadata<BYTE> imadata<T>::imaunsignedchar(bool etaldyn, bool aff) {
	int i,j,k,x;
	T val;
	imadata<BYTE> ima_u(nblig,nbcol,nbcanaux);
	if (etaldyn || !istat) statbasic(aff);
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			T *adval=(T*)Tsites[i*nbcol+j];
			for (k=0; k<nbcanaux; k++) {
				val=*(adval+k);
				if (val==valnul) x=0;
				else {
					if (etaldyn) x=int((val-valmin[k])*255./(valmax[k]-valmin[k]));
					else x=around((float)val);
				}
				ima_u(i,j,k)=(BYTE)(mini(maxi(x,0),255));
			}
		}
	}
	return ima_u;
}

template <class T> imadata<unsigned __int16> imadata<T>::imaunsignedshortint(bool etaldyn, bool aff) {
	int i,j,k,x;
	T val;
	imadata<unsigned __int16> ima_u(nblig,nbcol,nbcanaux);
	if (etaldyn || !istat) statbasic(aff);
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			T *adval=(T*)Tsites[i*nbcol+j];
			for (k=0; k<nbcanaux; k++) {
				val=*(adval+k);
				if (val==valnul) x=0;
				else {
					if (etaldyn) x=int((val-valmin[k])*65535./(valmax[k]-valmin[k]));
					else x=around((float)val);
				}
				ima_u(i,j,k)=(unsigned __int16)(mini(maxi(x,0),65535));
			}
		}
	}
	return ima_u;
}

template <class T> imadata<T> imadata<T>::operator + (imadata<T> &ima2) const {
	int nlig=mini(nblig,ima2.nlig()), ncol=mini(nbcol,ima2.ncol()),
	    ncanaux=mini(nbcanaux,ima2.ncanaux()), i, j, k;
	imadata<T> imaRes(nlig,ncol,ncanaux);
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++)
			for (k=0; k<ncanaux; k++)
				imaRes(i,j,k)=(*this)(i,j,k)+ima2(i,j,k);
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T> imadata<T>::operator - (imadata<T> &ima2) const {
	int nlig=mini(nblig,ima2.nlig()), ncol=mini(nbcol,ima2.ncol()),
	    ncanaux=mini(nbcanaux,ima2.ncanaux()), i, j, k;
	imadata<T> imaRes(nlig,ncol,ncanaux);
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++)
			for (k=0; k<ncanaux; k++)
				imaRes(i,j,k)=(*this)(i,j,k)-ima2(i,j,k);
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T> imadata<T>::operator * (imadata<T> &ima2) const {
	cout<<" multiplication pixel a pixel\n";
	int nlig=mini(nblig,ima2.nlig()), ncol=mini(nbcol,ima2.ncol()),
	    ncanaux=mini(nbcanaux,ima2.ncanaux()), i, j, k;
	imadata<T> imaRes(nlig,ncol,ncanaux);
	T valeur_nulle=imaRes.v_nulle(); cout<<" valeur_nulle = "<<valeur_nulle<<"\n";
	const float eps=(float)1.e-6;
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++)
			for (k=0; k<ncanaux; k++)
				if (ima2(i,j,k)!=valeur_nulle && ima2(i,j,k)!=valeur_nulle) imaRes(i,j,k)=(*this)(i,j,k)*ima2(i,j,k);
				else imaRes(i,j,k)=valeur_nulle;
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T> imadata<T>::operator / (imadata<T> &ima2) const {
	cout<<" division pixel a pixel\n";
	int nlig=mini(nblig,ima2.nlig()), ncol=mini(nbcol,ima2.ncol()),
	    ncanaux=mini(nbcanaux,ima2.ncanaux()), i, j, k;
	imadata<T> imaRes(nlig,ncol,ncanaux);
	T valeur_nulle=imaRes.v_nulle(); cout<<" valeur_nulle = "<<valeur_nulle<<"\n";
	const float eps=1.e-6f;
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++)
			for (k=0; k<ncanaux; k++)
				if (ima2(i,j,k)!=0 && ima2(i,j,k)!=valeur_nulle) imaRes(i,j,k)=(*this)(i,j,k)/ima2(i,j,k);
				else imaRes(i,j,k)=valeur_nulle;
/*				else
					if (abs((*this)(i,j,k))<eps) imaRes(i,j,k)=0;
					else imaRes(i,j,k)=valeur_nulle;*/
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T> imadata<T>::operator * (float x) const {
	int i, j, k;
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++)
				if ((*this)(i,j,k)!=valnul) imaRes(i,j,k)=(T)((*this)(i,j,k)*x);
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T> imadata<T>::operator + (float x) const {
	int i, j, k;
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++)
				imaRes(i,j,k)=(*this)(i,j,k)+x;
	imaRes.statbasic();
	return imaRes;
}

template <class T> imadata<T> imadata<T>::absI () const {
	int i, j, k;
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++)
				imaRes(i,j,k)=abs((*this)(i,j,k));
	imaRes.statbasic();
	return imaRes;
}

template <class T> bool imadata<T>::operator == (imadata<T> &ima2) {
	bool b=1;
	int i=0,j=0,k=0;
	while (b && i<nblig && j<nbcol && k<nbcanaux) {
		if ((*this)(i,j,k)!=ima2(i,j,k)) b=0;
		else {
			if (j<nbcol-1) j++;
			else {
				if (i<nblig-1) {i++; j=0;}
				else {k++; i=0; j=0;}
			}
		}
	}
	return b;
}

template <class T> bool imadata<T>::operator == (const imadata<T> &ima2) const {
	bool b=1;
	int i=0,j=0,k=0;
	while (b && i<nblig && j<nbcol && k<nbcanaux) {
		if ((*this)(i,j,k)!=ima2(i,j,k)) b=0;
		else {
			if (j<nbcol-1) j++;
			else {
				if (i<nblig-1) {i++; j=0;}
				else {k++; i=0; j=0;}
			}
		}
	}
	return b;
}

template <class T> void imadata<T>::sauve_Ima(string nomfich, int icanal) const {
	char *nomfic=(char *)nomfich.c_str();
	sauve_Ima(nomfic,icanal);
}

template <class T> void imadata<T>::sauve_Ima(char *nomfich, int icanal) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);
	}
	T *adval=(T*)Tsites[0];
	int itype=0, sizeval=0;
	if (typeid(valnul)==typeid(float)) {itype=3; sizeval=sizeof(float); }
	else
		if (typeid(valnul)==typeid(BYTE)) {itype=0; sizeval=sizeof(BYTE);}
		else {cout<<" type inconnu \n"; sizeval=0;}
//	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col., de type = "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (long int i=0; i<nbpix; i++) {
		adval=(T*)Tsites[i]+icanal;
		sortie.write((char *)adval,sizeval);
	}
	sortie.close();
}

template <class T> void imadata<T>::sauve_ImaBSQ(char *nomfich) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {cout<<" ouverture de "<<nomfich<<" impossible\n"; exit(-1);}
	T *adval=(T*)Tsites[0];
	int itype=0, sizeval=0;
	if (typeid(valnul)==typeid(float)) {itype=3; sizeval=sizeof(float); }
	else
		if (typeid(valnul)==typeid(BYTE)) {itype=0; sizeval=sizeof(BYTE);}
		else {
			if (typeid(valnul)==typeid(int)) {itype=0; sizeval=sizeof(int);}
			else {cout<<" type inconnu \n"; sizeval=0;}
		}
	cout<<" sauvegarde de image "<<nblig<<" lig.& "<<nbcol<<" col.& "<<nbcanaux<<" canaux, de type "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (int k=0; k<nbcanaux; k++)
		for (long int i=0; i<nbpix; i++) {adval=(T*)Tsites[i]+k; sortie.write((char *)adval,sizeval);}
	sortie.close();
}

template <class T> void imadata<T>::sauve_ImaBSQ(string nomfichier) const {
	ofstream sortie;
	char *nomfich=new char[nomfichier.size()+1];
	for (BYTE j=0; j<nomfichier.size(); j++) nomfich[j]=nomfichier[j];
	nomfich[nomfichier.size()]='\0';
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);
	}
	T *adval=(T*)Tsites[0];
	int itype=0, sizeval=0;
	if (typeid(valnul)==typeid(float)) {itype=3; sizeval=sizeof(float); }
	else
		if (typeid(valnul)==typeid(BYTE)) {itype=0; sizeval=sizeof(BYTE);}
		else {
			if (typeid(valnul)==typeid(int)) {itype=0; sizeval=sizeof(int);}
			else {cout<<" type inconnu \n"; sizeval=0;}
		}
//	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col.& "<<nbcanaux<<" canaux, de type "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (int k=0; k<nbcanaux; k++)
		for (long int i=0; i<nbpix; i++) {
			adval=(T*)Tsites[i]+k;
			sortie.write((char *)adval,sizeval);
		}
	sortie.close();
	if (nomfich!=NULL) {delete[] nomfich; nomfich=NULL;}
}

template <class T> void imadata<T>::sauve_ImaBIP(char *nomfich) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);
	}
	int itype=0, sizeval=0;
	T *adval=(T*)Tsites[0];
	if (typeid(valnul)==typeid(float)) {itype=3; sizeval=sizeof(float); }
	else
		if (typeid(valnul)==typeid(BYTE)) {itype=0; sizeval=sizeof(BYTE);}
		else {
			if (typeid(valnul)==typeid(int)) {itype=0; sizeval=sizeof(int);}
			else {cout<<" type inconnu \n"; sizeval=0;}
		}
//	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col.& "<<nbcanaux<<" canaux, de type "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (long int i=0; i<nbpix; i++)
		for (int k=0; k<nbcanaux; k++) {
			adval=(T*)Tsites[i]+k;
			sortie.write((char *)adval,sizeval);
		}
	sortie.close();
}

template <class T> void imadata<T>::sauve_ImaPGM(string nomfich) {
	char *nomfic=(char *)nomfich.c_str();
	sauve_ImaPGM(nomfic);
}

template <class T> void imadata<T>::sauve_ImaPGM(char *nomfich) {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {cout<<" ouverture de "<<nomfich<<" impossible\n";exit (-1);}
	double vmin=minI(), vmax=maxI(), x;
	for (int k=0; k<nbcanaux; k++) {
		x=(double)minI(k);
		if (x<vmin) vmin=x;
		x=(double)maxI(k);
		if (x<vmax) vmax=x;
	}
	bool ichgdyn=0;
	if (vmin<0 || vmax>255) {
		ichgdyn=1;
		cout<<" valmin = "<<vmin<<" valmax = "<<vmax;
		cout<<" => changement de la dynamique pour etre entre 0 et 255\n";
	}
	imadata<BYTE> ima_ui=imaunsignedchar(ichgdyn);
	char aa='P'; sortie.write((char *)&aa,sizeof(char));
	aa='5'; if (ncanaux()==3) aa='6'; sortie.write((char *)&aa,sizeof(char));
	aa='\n'; sortie.write((char *)&aa,sizeof(char));
	int n=nbcol, p=(int)floor(log10((float)n)), ii;
	while (p>=0) {
		ii=(int)floor(n/pow(10.f,p)); aa=ii+48;
		sortie.write((char *)&aa,sizeof(char));//cout<<aa;
		n-=ii*(int)pow(10.f,p); p--;
	}
//	cout<<"\n";
	aa=' '; sortie.write((char *)&aa,sizeof(char));
	n=nblig; p=(int)floor(log10((float)n));
	while (p>=0) {
		ii=(int)floor(n/pow(10.f,p)); aa=ii+48;
		sortie.write((char *)&aa,sizeof(char));//cout<<aa;
		n-=ii*(int)pow(10.f,p); p--;
	}
//	cout<<"\n";
	aa='\n'; sortie.write((char *)&aa,sizeof(char));
	n=255; p=(int)floor(log10((float)n));
	while (n>0) {
		ii=(int)floor(n/pow(10.f,p)); aa=ii+48;
		sortie.write((char *)&aa,sizeof(char));
		n-=ii*(int)pow(10.f,p); p=(int)floor(log10((float)n));
	}
	aa='\n'; sortie.write((char *)&aa,sizeof(char));
	BYTE *adval;
	for (int i=0; i<nblig; i++)
		for (int j=0; j<nbcol; j++)
			for (int k=0; k<nbcanaux; k++) {
				adval=&ima_ui(i,j,k);
				sortie.write((char *)adval,sizeof(char));
			}
	sortie.close();
}

template <class T> void imadata<T>::sauve_ImaPGM2(string nomfich) {
	char *nomfic=(char *)nomfich.c_str();
	sauve_ImaPGM2(nomfic);
}

template <class T> void imadata<T>::sauve_ImaPGM2(char *nomfich) {cout<<" entree dans sauve_ImaPGM2\n";
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {cout<<" ouverture de "<<nomfich<<" impossible\n";exit (-1);}
	double vmin=minI(), vmax=maxI(), x;
	for (int k=0; k<nbcanaux; k++) {
		x=(double)minI(k); if (x<vmin) vmin=x;
		x=(double)maxI(k); if (x<vmax) vmax=x;
	}
	bool ichgdyn=0;
	if (vmin<0 || vmax>65535) {
		ichgdyn=1; cout<<" valmin = "<<vmin<<" valmax = "<<vmax<<" => changement de la dynamique pour etre entre 0 et 65535\n";
	}
	imadata<unsigned __int16> ima_ui=imaunsignedshortint(ichgdyn);
	char aa='P'; sortie.write((char *)&aa,sizeof(char));
	aa='5'; if (ncanaux()==3) aa='6'; sortie.write((char *)&aa,sizeof(char));
	aa='\n'; sortie.write((char *)&aa,sizeof(char));
	int n=nbcol, p=(int)floor(log10((float)n)), ii;
	while (p>=0) {
		ii=(int)floor(n/pow(10.f,p)); aa=ii+48;
		sortie.write((char *)&aa,sizeof(char));//cout<<aa;
		n-=ii*(int)pow(10.f,p); p--;
	}
//	cout<<"\n";
	aa=' '; sortie.write((char *)&aa,sizeof(char));
	n=nblig; p=(int)floor(log10((float)n));
	while (p>=0) {
		ii=(int)floor(n/pow(10.f,p)); aa=ii+48;
		sortie.write((char *)&aa,sizeof(char));//cout<<aa;
		n-=ii*(int)pow(10.f,p); p--;
	}
//	cout<<"\n";
	aa='\n'; sortie.write((char *)&aa,sizeof(char));
	aa='6'; sortie.write((char *)&aa,sizeof(char));
	aa='5'; sortie.write((char *)&aa,sizeof(char));
	aa='5'; sortie.write((char *)&aa,sizeof(char));
	aa='3'; sortie.write((char *)&aa,sizeof(char));
	aa='5'; sortie.write((char *)&aa,sizeof(char));
	aa='\n'; sortie.write((char *)&aa,sizeof(char));
	unsigned __int16 *adval;
	for (int i=0; i<nblig; i++)
		for (int j=0; j<nbcol; j++)
			for (int k=0; k<nbcanaux; k++) { //cout<<ima_ui(i,j,k)<<" ";
				adval=&ima_ui(i,j,k);
				sortie.write((char *)adval+1,sizeof(BYTE));
				sortie.write((char *)adval,sizeof(BYTE));
//				sortie.write((char *)adval,sizeof(unsigned __int16));
			}
	sortie.close();
}

template <class T> void imadata<T>::statbasic(bool affich) {
	long int i,ii; int j;
	double x, *vmin=new double[nbcanaux], *vmax=new double[nbcanaux];
	for (j=0; j<nbcanaux; j++) {
		vmin[j]=+DBL_MAX; vmax[j]=-DBL_MAX;
		valmoy[j]=valvar[j]=0.;
		unsigned long int n=0;
		for (i=0; i<nbpix; i++) {
			ii=i*nbcanaux;
//			x=(double)Tima[ii+j]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
			x=(double)(*((T*)Tsites[i]+j));
			if (x!=(double)valnul) {
				if (x<vmin[j]) vmin[j]=x;
				if (x>vmax[j]) vmax[j]=x;
				valmoy[j]+=x;
				valvar[j]+=pow(x,2);
				n++;
			}
		}
		if (n>0) valmoy[j]=valmoy[j]/n;
		if (n>1) valvar[j]=(valvar[j]-n*pow(valmoy[j],2))/(n-1);
		else valvar[j]=-999;
	}
	for (j=0; j<nbcanaux; j++) {
		valmin[j]=(T)vmin[j];
		valmax[j]=(T)vmax[j];
	}
	if (vmin!=NULL) {delete[] vmin; vmin=NULL;}
	if (vmax!=NULL) {delete[] vmax; vmax=NULL;}
	istat=1;
	if (affich) {
		cout<<" canal |  minimum  |  maximum  |  moyenne  | variance\n";
		for (j=0; j<nbcanaux; j++) {
			cout<<setw(6)<<j<<"   "<<setw(8)<<(float)valmin[j]<<"   "<<setw(9)<<(float)valmax[j];
			cout<<"   "<<setw(9)<<(float)valmoy[j]<<"   "<<setw(9)<<(float)valvar[j]<<"\n";
		}
	}
}

template <class T> double imadata<T>::modeI (double &prec, int i) {
	if (!ihist) {
/*		unsigned int n=0;
		for (long int j=0; j<nbpix; j++)
			if ((*((T*)Tsites[j]+i))!=valnul) n++;
		if (n>0) histogramme(n/10+2);*/
		histogramme(256,0,0);
	}
	const float nsigma=6.;
	prec=nsigma*sqrt(valvar[i])/(nbinH-2);
	return valmod[i];
}

template <class T> double imadata<T>::minI (int i) {
	if (!istat) statbasic();
	return valmin[i];
	}

template <class T> double imadata<T>::maxI (int i) {
	if (!istat) statbasic();
	return valmax[i];
	}

template <class T> double imadata<T>::moyI (int i) {
	if (!istat) statbasic();
	return valmoy[i];
}

template <class T> double imadata<T>::moyI (int i) const {
	if (!istat) {cout<<" stat non sauvegardees car image const...\n";
		double x=0.; 
		for (int i=0; i<nblig; ++i)
			for (int j=0; j<nbcol; ++j) x+=(*this)(i,j);
		return x/(nblig*nbcol);
	}
	else return valmoy[i];
}

template <class T> double imadata<T>::varI (int i) {
	if (!istat) statbasic();
	return valvar[i];
}

template <class T> double imadata<T>::varI (int i) const {
	if (!istat) {cout<<" stat non sauvegardees car image const...\n";
		double x=0., xx=0.; 
		for (int i=0; i<nblig; ++i)
			for (int j=0; j<nbcol; ++j) {x+=(*this)(i,j); xx+=(*this)(i,j)*(*this)(i,j);}
		int n=nblig*nbcol;
		return (xx-x*x/n)/n;
	}
	else return valvar[i];
}

template <class T> float imadata<T>::histogramme_allc (int nbin, float pc, bool affich) { 
	if (affich) cout<<" calcul de l'histogramme sur tous les canaux\n";
	const float nsigma=6.f ;
	if (nbin!=nbinH) {
		nbinH=nbin;
		if (histo!=NULL) delete[] histo;
		histo=new eltHisto[nbinH];
	}
	int j,k,n; T x;
	float xmin, xmax, pas;
	statbasic(1);
	xmin=valmin[0]; xmax=valmax[0];
	for (j=0; j<nbcanaux; j++) {xmin=mini((float)valmin[j],xmin); xmax=maxi((float)valmax[j],xmax);}
	pas=(xmax-xmin)/(float)nbin;
	if (affich) cout<<" pas de l'histo = "<<pas<<", sup 1er bin = "<<xmin<<"\n";
	for (j=0; j<nbinH; j++) {
		histo[j].nbelts=0.f;
		histo[j].bornesup=xmin+j*pas;
	}
	histo[nbinH-1].bornesup=maxi((float)(xmax+pas/10.),(float)(histo[nbinH-1].bornesup));
	for (k=0; k<nbcanaux; k++)
		for (j=0; j<nbpix; j++) {
			x=*((T*)Tsites[j]+k);
			if (x!=valnul) {
				if (x<=xmin) (histo[0].nbelts)+=1.;
				else {
					n=(int)mini((floor)((x-xmin)/pas),(float)nbinH-1); // cas d'un histogramme ?pas constant
					(histo[n].nbelts)+=1.;
				}
			}
		}
	float x_nelt=0.f, n_pc;	
	for (j=0; j<nbinH; j++) x_nelt+=histo[j].nbelts;
	if (affich) {
		cout<<" histogramme global tous canaux\n";
		for (j=0; j<nbinH; j++) cout<<setw(3)<<j<<"   "<<setw(8)<<histo[j].bornesup<<"   "<<setw(8)<<histo[j].nbelts<<"\n";
		cout<<" en tout "<<x_nelt<<" pixels\n";
	}
	ihist=1;
	n_pc=pc*x_nelt; //cout<<" n_pc "<<n_pc<<"\n";
	x_nelt=0.f; j=0;
	while (x_nelt<n_pc && j<nbinH) {x_nelt+=histo[j++].nbelts; /*cout<<x_nelt<<" ";*/} //cout<<"\n";
	if (j>0 && x_nelt-n_pc>n_pc-(x_nelt-histo[j-1].nbelts)) j--;
	j=mini((int)j,nbinH-1);
	return (float)(histo[j].bornesup);
}

template <class T> void imadata<T>::histogramme (int nbin, bool affich, bool ioptim) { 
	if (affich) cout<<" calcul de l'histogramme canal par canal\n";
	const float nsigma=6.f ;
	if (nbin!=nbinH) {
		nbinH=nbin;
		if (histo!=NULL) delete[] histo;
		histo=new eltHisto[nbinH*nbcanaux];
	}
	int i,ii,n; long int j; T x;
	float x_nelt, *xmin=new float[nbcanaux], *pas=new float[nbcanaux];
//	if (!istat) statbasic (affich);
	statbasic(1);
	for (j=0; j<nbcanaux; j++) {
		pas[j]=mini((float)(valmax[j]-valmin[j])/nbin,(float)(nsigma*sqrt(valvar[j])/(nbinH-2)));
		xmin[j]=maxi((float)valmin[j],(float)(valmoy[j]-nsigma/2.*sqrt(valvar[j])));
		if (valmax[j]<xmin[j]+(nbinH-1)*pas[j] && valmin[j]<xmin[j]) xmin[j]=valmax[j]-(nbinH-1)*pas[j];
		if (valmin[j]>xmin[j] && valmax[j]>xmin[j]+(nbinH-1)*pas[j]) xmin[j]=valmin[j]+pas[j];
		if (!ioptim) {xmin[j]=(float)valmin[j]; pas[j]=(float)(valmax[j]-valmin[j])/nbin;}
		if (affich) cout<<" pas de l'histo = "<<pas[j]<<", sup 1er bin = "<<xmin[j]<<"\n";
	}
	for (i=0; i<nbcanaux; i++) {
		ii=i*nbinH;
		for (j=0; j<nbinH; j++) {
			histo[ii+j].nbelts=0.f;
			histo[ii+j].bornesup=xmin[i]+j*pas[i];
		}
		histo[ii+nbinH-1].bornesup=maxi((float)(valmax[i]+pas[i]/10.),(float)(histo[ii+nbinH-1].bornesup));
	}
	for (i=0; i<nbcanaux; i++) {
		ii=i*nbinH;
		for (j=0; j<nbpix; j++) {
//			x=Tima[j*nbcanaux+i]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
			x=*((T*)Tsites[j]+i);
			if (x!=valnul) {
				if (x<=xmin[i]) (histo[ii].nbelts)+=1.;
				else {
					n=0;
					do n++; while (n<nbinH && x>histo[ii+n].bornesup);
					(histo[ii+n].nbelts)+=1.;
				}
			}
		}
	}
	for (i=0; i<nbcanaux; i++) {
		ii=i*nbinH;
		x_nelt=0.f;
		int jmod=0;
		for (j=0; j<nbinH; j++)
			if (histo[ii+j].nbelts>x_nelt) {
				x_nelt=histo[ii+j].nbelts;
				jmod=j;
			}
		valmod[i]=(T)(histo[ii+jmod].bornesup-pas[i]/2.f);
	}
	if (affich)
		for (i=0; i<nbcanaux; i++) {
			ii=i*nbinH;
			cout<<" canal "<<i<<"\n";
			x_nelt=0.f;
			for (j=0; j<nbinH; j++) {
				cout<<setw(3)<<j<<"   ";
				cout<<setw(8)<<histo[ii+j].bornesup<<"   ";
				cout<<setw(8)<<histo[ii+j].nbelts<<"\n";
				x_nelt+=histo[ii+j].nbelts;
			}
			cout<<" en tout "<<x_nelt<<" pixels\n";
		}
	ihist=1;
	if (pas!=NULL) {delete[] pas; pas=NULL;}
	if (xmin!=NULL) {delete[] xmin; xmin=NULL;}	
}

template <class T> float imadata<T>::seuil_otsu (imabin &imab, int k) {
	statbasic(); histogramme();
	double xx=histo[k*nbinH+0].bornesup, mu1=0, n1=histo[k*nbinH+0].nbelts, mu2=valmoy[k], n2=nblig*nbcol, D_inter, D_inter_max, D_inter_max2=0.;
	int i,j,jmax=0,jmax2=-1;
	for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) if ((*this)(i,j,k)<=xx) mu1+=(*this)(i,j,k);
	mu2=mu2*n2-mu1;
	n2-=n1; 
	if (n2>0) mu2/=n2; else {cout<<" 0- n2 nul ? "<<n2<<"\n"; /*char aa; cin>>aa;*/}
	if (n1>0) mu1/=n1; else {cout<<" 0- n1 nul ? "<<n1<<"\n"; /*char aa; cin>>aa;*/}
	D_inter_max=n1*n2*pow(mu1-mu2,2);
	for (j=1; j<nbinH-1; j++) {
		xx=(histo[k*nbinH+j].bornesup+histo[k*nbinH+j-1].bornesup)/2.*histo[k*nbinH+j].nbelts;
		mu1=mu1*n1+xx;
		n1+=histo[k*nbinH+j].nbelts; if (n1>0) mu1/=n1; else {cout<<j<<"- n1 nul ? "<<n1<<"\n"; /*char aa; cin>>aa;*/}
		mu2=mu2*n2-xx;
		n2-=histo[k*nbinH+j].nbelts; if (n2>0) mu2/=n2; else {cout<<j<<"- n2 nul ? "<<n2<<"\n"; /*char aa; cin>>aa;*/}
		D_inter=n1*n2*pow(mu1-mu2,2);
//		cout<<j<<" "<<mu1<<" "<<mu2<<" "<<n1<<" "<<n2<<" "<<D_inter<<" "<<D_inter_max<<"\n";
		if (D_inter>D_inter_max) {D_inter_max2=D_inter_max; jmax2=jmax; D_inter_max=D_inter; jmax=j;}
	}
	double s1=histo[k*nbinH+jmax].bornesup, s2=histo[k*nbinH+jmax2].bornesup;
	if (jmax>0) s1=(s1+histo[k*nbinH+jmax-1].bornesup)/2.;
	if (jmax2>0) s2=(s2+histo[k*nbinH+jmax2-1].bornesup)/2.;
	cout<<" 1er seuil optimal trouve = "<<s1<<", 2nd seuil optimal trouve = "<<s2<<"\n";
	imab=imabin(*this,(float)s1);
//	return imabin(*this,s1);
	return (float)s1;
}

template <class T> void imadata<T>::supprime_Niveaux_vide (bool affich) {
	if (!istat) statbasic (affich);
	const int nbin=255;
	int T_corr[nbin];
	int i,j,k;
	bool trou=0;
	for (k=0; k<nbin; k++) T_corr[k]=k;
	histogramme(nbin,affich);
	T vmin, vmax, v;
	vmin=valmin[0];
	for (k=0; k<nbin; k++) {
		if (histo[k].nbelts>0) {
			if (trou) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						v=(*this)(i,j);
						if (v>vmin) (*this)(i,j)=v-(vmax-vmin);
					}
			}
			vmin=(T)(histo[k].bornesup);
		} else {
			trou=1;
			vmax=(T)(histo[k].bornesup);
		}
	}
}

template <class T> float imadata<T>::seuil_percentile (float pc, int k) {
	if (!ihist) histogramme((int)256,(bool)0);
	int ii=k*nbinH, j_pc=0, j=0;
	float x_nelt=0.f, n_pc=pc*nbpix;
	while (x_nelt<n_pc && j<nbinH) {x_nelt+=histo[ii+j].nbelts; j++;};
	j_pc=j;
	if (j_pc>0 && x_nelt-n_pc>n_pc-(x_nelt-histo[ii+j_pc-1].nbelts)) j_pc--;
	return histo[ii+j_pc].bornesup;
}

template <class T> void imadata<T>::histogramme (imabin &uMasq, int nbin, bool affich, bool ioptim) { cout<<" calcul de l'histo avec masque\n";
	const float nsigma=6.f;
	if (nbin!=nbinH) {
		nbinH=nbin;
		if (histo!=NULL) delete[] histo;
		histo=new eltHisto[nbinH*nbcanaux];
	}
	int i,ii,n; long int j;
	T x;
	float *xmin=new float[nbcanaux], *xmax=new float[nbcanaux], x_nelt;
	statbasic (affich);
	float *pas=new float[nbcanaux];
	for (j=0; j<nbcanaux; j++) {
		xmin[j]=maxi((float)(valmoy[j]-nsigma/2.*sqrt(valvar[j])),(float)valmin[j]);
		xmax[j]=mini((float)(valmoy[j]+nsigma/2.*sqrt(valvar[j])),(float)valmax[j]);
		if (!ioptim) {xmin[j]=(float)valmin[j]; xmax[j]=(float)valmax[j];}
		pas[j]=(xmax[j]-xmin[j])/nbin; xmin[j]+=pas[j];
/*		pas[j]=nsigma*(float)sqrt(valvar[j])/(nbinH-2);
		xmin[j]=(float)(valmoy[j]-nsigma/2.*sqrt(valvar[j]));
		if (valmax[j]<xmin[j]+(nbinH-1)*pas[j] && valmin[j]<xmin[j])
			xmin[j]=valmax[j]-(nbinH-1)*pas[j];
		if (valmin[j]>xmin[j] && valmax[j]>xmin[j]+(nbinH-1)*pas[j])
			xmin[j]=valmin[j]+pas[j];*/
		if (affich) cout<<" pas = "<<pas[j]<<", sup 1er bin = "<<xmin[j]<<"\n";
	}
	for (i=0; i<nbcanaux; i++) {
		ii=i*nbinH;
		histo[ii].nbelts=0.f; histo[ii].bornesup=xmin[i];
		for (j=1; j<nbinH-1; j++) {
			histo[ii+j].nbelts=0.f;
			histo[ii+j].bornesup=xmin[i]+j*pas[i];
		}
		histo[ii+nbinH-1].nbelts=0.f;
		histo[ii+nbinH-1].bornesup=(float)(valmax[i]+pas[i]/10.f);
	}
	for (i=0; i<nbcanaux; i++) {
		ii=i*nbinH;
		for (j=0; j<nbpix; j++) {
//			x=Tima[j*nbcanaux+i]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
			x=*((T*)Tsites[j]+i); //cout<<j<<" "<<x<<" "<<histo[ii].nbelts<<"\n";
			if (x!=valnul && uMasq(j/nbcol,j%nbcol)==1) {
				if (x<=xmin[i]) (histo[ii].nbelts)+=1.;
				else {
					n=0;
					do n++;
					while (n<nbinH && x>histo[ii+n].bornesup);
					(histo[ii+n].nbelts)+=1.;
				}
			}
		}
	}
	for (i=0; i<nbcanaux; i++) {
		ii=i*nbinH;
		x_nelt=0.f;
		int jmod=0;
		for (j=0; j<nbinH; j++)
			if (histo[ii+j].nbelts>x_nelt) {
				x_nelt=histo[ii+j].nbelts;
				jmod=j; //cout<<jmod<<" "<<x_nelt<<" "<<histo[ii+jmod].bornesup<<"\n";
			}
		valmod[i]=(T)(histo[ii+jmod].bornesup-pas[i]/2.f);
	}
	if (affich)
		for (i=0; i<nbcanaux; i++) {
			ii=i*nbinH;
			cout<<" canal "<<i<<"\n";
			x_nelt=0.f;
			for (j=0; j<nbinH; j++) {
				cout<<setw(3)<<j<<"   ";
				cout<<setw(8)<<histo[ii+j].bornesup<<"   ";
				cout<<setw(8)<<histo[ii+j].nbelts<<"\n";
				x_nelt+=histo[ii+j].nbelts;
			}
			cout<<" en tout "<<x_nelt<<" pixels, mode en "<<valmod[i]<<"\n";
		}
	ihist=1;
	if (pas!=NULL) {delete[] pas; pas=NULL;}
	if (xmin!=NULL) {delete[] xmin; xmin=NULL;}
	if (xmax!=NULL) {delete[] xmax; xmax=NULL;}
}

template <class T> float imadata<T>::seuil_percentile (imabin &uMasq, float pc, int k) {
	histogramme(uMasq,256,0,0);
	int ii=k*nbinH, j_pc=0, j=0;
	float x_nelt=0, n_pc=pc*uMasq.norm();
//	cout<<" pc="<<pc<<" uMasq.norm()="<<uMasq.norm()<<" nlig="<<nblig<<" ncol="<<nbcol<<" n_pc = "<<n_pc<<"\n";
	while (x_nelt<n_pc && j<nbinH) {x_nelt+=histo[ii+j].nbelts; j++;};
	j_pc=j; //cout<<" j_pc = "<<j_pc<<" => seuil = "<<histo[ii+j_pc].bornesup<<"\n";
	if (j_pc>0 && x_nelt-n_pc>n_pc-(x_nelt-histo[ii+j_pc-1].nbelts)) j_pc--;
	return histo[ii+j_pc].bornesup;
}

template <class T> void imadata<T>::egalisehisto (int nbin, bool affich) {
	int npix_bin,j,ii2,jj,n; long int i,ii;
	const int nbin2=512;
	statbasic(affich);
	ihist=0; cout<<" calcul de l'histo fin\n"; (*this).histogramme(nbin2,affich); 
	float x_nelt=0.f, **x_nelt_max=new float*[nbcanaux], *pas=new float[nbcanaux];
	for (i=0; i<nbcanaux; i++) {
		x_nelt_max[i]=new float[3];
		for (n=0; n<3; n++) x_nelt_max[i][n]=0;
		ii2=i*nbinH;                     // nbinH=nbin2 après appel 'histogramme(nbin2,affich)'
		for (j=0; j<nbinH; j++) {
			if (histo[ii2+j].nbelts>x_nelt_max[i][0]) x_nelt_max[i][0]=histo[ii2+j].nbelts;
			else {
				if (histo[ii2+j].nbelts>x_nelt_max[i][1]) x_nelt_max[i][1]=histo[ii2+j].nbelts;
				else 
					if (histo[ii2+j].nbelts>x_nelt_max[i][2]) x_nelt_max[i][2]=histo[ii2+j].nbelts;
			}
		}
		if (affich) cout<<" maxima de l'histo fin : "<<x_nelt_max[i][0]<<" "<<x_nelt_max[i][1]<<" "<<x_nelt_max[i][2]<<"\n";
	}
	eltHisto *histoE=NULL; histoE=new eltHisto[nbin*nbcanaux];
	T x,y,tab[nbin2];
	for (i=0; i<nbcanaux; i++) {
		ii2=i*nbinH;
		x_nelt=0.f;
		for (j=0; j<nbinH; j++) x_nelt+=histo[ii2+j].nbelts;
		npix_bin=(int)maxi(x_nelt/nbin,x_nelt_max[i][2]); cout<<" npix par bin "<<npix_bin<<"\n";
		ii=i*nbin;
		x_nelt=0.f; jj=-1;
		for (j=0; j<nbinH; j++) {         // on balaie l'histogramme 'fin' de l'image initiale
			x_nelt+=histo[ii2+j].nbelts;  // on cumule jusqu'?depasser le nombre de pixels par bin sur l'histo égalis?
			if (x_nelt>=npix_bin) {
				if ((x_nelt-npix_bin)>(npix_bin-x_nelt+histo[ii2+j].nbelts)) {x_nelt-=histo[ii2+j].nbelts; j--;}
				jj++;
				if (jj>=nbin) {
					cout<<" Pb dim histoE jj = "<<jj<<" >= "<<nbin<<"\n";
					(histoE[ii+nbin-1].nbelts)+=x_nelt;
				}
				else { 
					histoE[ii+jj].nbelts=x_nelt;
					histoE[ii+jj].bornesup=histo[ii2+j].bornesup;
				}
				x_nelt-=npix_bin; 
			}
/*			if (x_nelt>=npix_bin) {
				while (x_nelt>=npix_bin && jj<nbin) {
					jj++;
					if (jj>=nbin) {
						cout<<" Pb dim histoE jj = "<<jj<<" >= "<<nbin<<"\n";
						(histoE[ii+nbin-1].nbelts)+=x_nelt;
					}
					else {
						histoE[ii+jj].nbelts=x_nelt;
						x_nelt-=npix_bin;
						histoE[ii+jj].bornesup=histo[ii2+j].bornesup;
					}
				}
			}*/
		}
		if (x_nelt>0) {
			jj++;
			if (jj>=nbin) {
				cout<<" Pb dim histoE jj = "<<jj<<" >= "<<nbin<<" n = "<<x_nelt<<"\n";
				(histoE[ii+nbin-1].nbelts)+=x_nelt;
			}
			else {
				histoE[ii+jj].nbelts=x_nelt;
				histoE[ii+jj].bornesup=histo[ii2+nbinH-1].bornesup;
			}
		}
		for (j=jj+1; j<nbin; j++) {
			histoE[ii+j].nbelts=0.f;
			histoE[ii+j].bornesup=0.f;
		}
	}
	if (affich) {
		for (i=0; i<nbcanaux; i++) {
			ii=i*nbin;
			cout<<" canal "<<i<<"\n";
			x_nelt=0;
			for (j=0; j<nbin; j++) {
				if (histoE[ii+j].nbelts>0) cout<<ii+j<<"   "<<histoE[ii+j].bornesup<<"   "<<histoE[ii+j].nbelts<<"\n";
				x_nelt+=histoE[ii+j].nbelts;
			}
			cout<<" en tout "<<x_nelt<<" pixels\n";
		}
	}
	for (j=0; j<nbcanaux; j++) { 
		pas[j]=(valmax[j]-valmin[j])/nbin;
		for (i=0; i<nbin2; i++) {
			x=histo[j*nbin2+i].bornesup;
			n=-1;
			do n++;
			while ((n<nbin-1) && (histoE[j*nbin+n].bornesup<x));
			if (n==0) {x=valmin[j]; y=histo[j*nbin2].bornesup;}
			else {x=valmin[j]+n*pas[j]; y=x+pas[j];}
			tab[i]=x+(y-x)/2; //if (affich) cout<<i<<" -> "<<tab[i]<<"\n";
		}
		pas[j]=histo[j*nbin2+1].bornesup-histo[j*nbin2].bornesup;
		for (i=0; i<nbpix; i++) {
			ii=i*nbcanaux;
//			n=mini(maxi((int)floor((Tima[ii+j]-histo[j*nbin2].bornesup)/pas[j])+1,0),nbin2-1); // non valable dans le cas o?on a sous-echantillonn?l'image !!!
			n=mini(maxi((int)floor((*((T*)Tsites[i]+j)-histo[j*nbin2].bornesup)/pas[j])+1,0),nbin2-1);
//			Tima[ii+j]=tab[n]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
			*((T*)Tsites[i]+j)=tab[n];
		}
	}
/*	for (i=0; i<nbpix; i++) {
		ii=i*nbcanaux;
		for (j=0; j<nbcanaux; j++) {
//			x=Tima[ii+j]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
			x=*((T*)Tsites[i]+j);
			n=-1;
			do n++;
			while ((n<nbin-1) && (histoE[j*nbin+n].bornesup<x));
			if (n==0) x=valmin[j];
			else x=valmin[j]+n*pas[j];
			y=x+pas[j];
//			Tima[ii+j]=x+(y-x)/2; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
			*((T*)Tsites[i]+j)=x+(y-x)/2;
		}
	}*/
	statbasic();
	if (histoE!=NULL) {delete[] histoE; histoE=NULL;}
	if (pas!=NULL) {delete[] pas; pas=NULL;}
	for (i=0; i<nbcanaux; i++) if (x_nelt_max[i]!=NULL) delete[] x_nelt_max[i];
	if (x_nelt_max!=NULL) delete[] x_nelt_max;
    }

template <class T> void imadata<T>::egalisehisto (imabin &uMasq, int nbin, bool affich) {
	int npix_bin,j,ii2,jj; long int i,ii;
	const int nbin2=512;
	float x_elt;
	float *vmin=NULL; vmin=new float[nbcanaux];
	float *vmax=NULL; vmax=new float[nbcanaux];
	float *pas=NULL; pas=new float[nbcanaux];
	eltHisto *histoE=NULL; histoE=new eltHisto[nbin*nbcanaux];
	histogramme(uMasq,nbin2,affich);
	T x,y;
	for (i=0; i<nbcanaux; i++) {
		ii2=i*nbinH;                     // nbinH=nbin2 apres appel 'histogramme(nbin2,affich)'
		x_elt=0;
		for (j=0; j<nbinH; j++) x_elt+=histo[ii2+j].nbelts;
		npix_bin=x_elt/nbin;
//		cout<<" n = "<<n<<" npix_bin = "<<npix_bin<<"\n";
		ii=i*nbin;
		x_elt=0.f;
		jj=-1;
		for (j=0; j<nbinH; j++) {         // on balaie l'histogramme 'fin' de l'image initiale
			x_elt+=histo[ii2+j].nbelts;
			if (x_elt>=npix_bin) {
				jj++;
				if (jj>=nbin) {
					cout<<" Pb dim histoE jj = "<<jj<<" >= "<<nbin<<"\n";
					(histoE[ii+nbin-1].nbelts)+=x_elt;
					}
				else {
					histoE[ii+jj].nbelts=x_elt;
					x_elt-=(histoE[ii+jj].nbelts);
					histoE[ii+jj].bornesup=histo[ii2+j].bornesup;
					}
				}
			}
		if (x_elt>0) {
			jj++;
			if (jj>=nbin) {
				cout<<" Pb dim histoE jj = "<<jj<<" >= "<<nbin<<" n = "<<x_elt<<"\n";
				(histoE[ii+nbin-1].nbelts)+=x_elt;
				}
			else {
				histoE[ii+jj].nbelts=x_elt;
//				if (jj>0) x_elt-=(histoE[ii+jj-1].nbelts);
				histoE[ii+jj].bornesup=histo[ii2+nbinH-1].bornesup;
				}
			}
		for (j=jj+1; j<nbin; j++) {
			histoE[ii+j].nbelts=0.f;
			histoE[ii+j].bornesup=0.f;
			}
		}
	if (affich) {
		for (i=0; i<nbcanaux; i++) {
			ii=i*nbin;
			cout<<" canal "<<i<<"\n";
			x_elt=0;
			for (j=0; j<nbin; j++) {
				cout<<ii+j<<"   "<<histoE[ii+j].bornesup<<"   ";
				cout<<histoE[ii+j].nbelts<<"\n";
				x_elt+=histoE[ii+j].nbelts;
				}
			cout<<" en tout "<<x_elt<<" pixels\n";
			}
		}
	for (j=0; j<nbcanaux; j++) {
		vmin[j]=valmax[j];
		vmax[j]=valmin[j];
		}
	for (i=0; i<nbpix; i++) {
		ii=i*nbcanaux;
		if (uMasq(i/nbcol,i%nbcol)==1) {
			for (j=0; j<nbcanaux; j++) {
//				x=Tima[ii+j]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
				x=*((T*)Tsites[i]+j);
				if (x<vmin[j]) vmin[j]=x;
				if (x>vmax[j]) vmax[j]=x;
				}
			}
		}
	for (j=0; j<nbcanaux; j++) pas[j]=(vmax[j]-vmin[j])/nbin;
	int n;
	for (i=0; i<nbpix; i++) {
		ii=i*nbcanaux;
		if (uMasq(i/nbcol,i%nbcol)==1) {
			for (j=0; j<nbcanaux; j++) {
//				x=Tima[ii+j]; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
				x=*((T*)Tsites[i]+j);
				n=-1;
				do n++;
				while ((n<nbin-1) && (histoE[j*nbin+n].bornesup<x));
				if (n==0) x=vmin[j];
				else x=vmin[j]+n*pas[j];
				y=x+pas[j];
//				Tima[ii+j]=x+(y-x)/2; // non valable dans le cas o?on a sous-echantillonn?l'image !!!
				*((T*)Tsites[i]+j)=x+(y-x)/2;
				}
			}
		}
	if (histoE!=NULL) {delete[] histoE; histoE=NULL;}
	if (pas!=NULL) {delete[] pas; pas=NULL;}
	if (vmin!=NULL) {delete[] vmin; vmin=NULL;}
	if (vmax!=NULL) {delete[] vmax; vmax=NULL;}
    }

template <class T> imadata<T> imadata<T>::ima_diffligne () const {
	int i,j,k;
	imadata<T> imaRes(nblig-1,nbcol,nbcanaux);
	for (i=0; i<nblig-1; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++)
				imaRes(i,j,k)=(*this)(i+1,j,k)-(*this)(i,j,k);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::ima_diffcolonne () const {
	int i,j,k;
	imadata<T> imaRes(nblig,nbcol-1,nbcanaux);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol-1; j++)
			for (k=0; k<nbcanaux; k++)
				imaRes(i,j,k)=(*this)(i,j+1,k)-(*this)(i,j,k);
	return imaRes;
}

template <class T> void imadata<T>::mise_a_zero () {
	int i,j,k;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++) (*this)(i,j,k)=(T)0;
}

template <class T> float imadata<T>::norm() const {
	int i,j,k;
	float n=0.f;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++)
				n+=(*this)(i,j,k);
	return n;
}

template <class T> void imadata<T>::normalize_all(float val_max) {
	int i,j,k;
	float xmin,xmax,fact;
	for (k=0; k<nbcanaux; k++) {
		xmin=valmin[k]; xmax=valmax[k];
		if (xmax-xmin>0) {
			fact=val_max/(xmax-xmin);
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) (*this)(i,j,k)=((*this)(i,j,k)-xmin)*fact;
		}
	}
}

/*template <class T> imadata<unsigned int> imadata<T>::ima_Integrale_16bin(T valmax, T valmin) {
	const int nbin=16;
	float fact=(float)nbin/(float)(valmax-valmin);
	imadata<unsigned int> imaIntegrale(nblig,nbcol,nbin);
	imaIntegrale.mise_a_zero();
	int i=0, j=0, i0, j0, ii, jj, k;
	imaIntegrale(i,j,int(((*this)(i,j)-valmin)*fact))++;
	do {
		i0=i; j0=j; i=mini(i+1,nblig-1); j=mini(j+1,nbcol-1);
		if (j>j0) {
			for (k=0; k<nbin; k++) imaIntegrale(0,j,k)=imaIntegrale(0,j-1,k); imaIntegrale(0,j,int(((*this)(0,j)-valmin)*fact))++;
			for (ii=1; ii<i; ii++) {
				for (k=0; k<nbin; k++) imaIntegrale(ii,j,k)=imaIntegrale(ii,j-1,k)+imaIntegrale(ii-1,j,k)-imaIntegrale(ii-1,j-1,k); 
				imaIntegrale(ii,j,int(((*this)(ii,j)-valmin)*fact))++;
			}
		}
		if (i>i0) {
			for (k=0; k<nbin; k++) imaIntegrale(i,0,k)=imaIntegrale(i-1,0,k); imaIntegrale(i,0,int(((*this)(i,0)-valmin)*fact))++;
			for (jj=1; jj<j; jj++) {
				for (k=0; k<nbin; k++) imaIntegrale(i,jj,k)=imaIntegrale(i,jj-1,k)+imaIntegrale(i-1,jj,k)-imaIntegrale(i-1,jj-1,k); 
				imaIntegrale(i,jj,int(((*this)(i,jj)-valmin)*fact))++;
			}
		}
		for (k=0; k<nbin; k++) imaIntegrale(i,j,k)=imaIntegrale(i,j-1,k)+imaIntegrale(i-1,j,k)-imaIntegrale(i-1,j-1,k); 
		imaIntegrale(i,j,int(((*this)(i,j)-valmin)*fact))++;
	} while (i<nblig-1 || j<nbcol-1);
	return imaIntegrale;
}*/

template <class T> void imadata<T>::ima_Integrale_16bin(imadata<unsigned int> &imaIntegrale, T valmax, T valmin) {
	const int nbin=16;
	float fact=(float)nbin/(float)(valmax+1-valmin);
	imaIntegrale.mise_a_zero();
	int i=0, j=0, i0, j0, ii, jj, k;
	imaIntegrale(i,j,int(((*this)(i,j)-valmin)*fact))++;
	do {
		i0=i; j0=j; i=mini(i+1,nblig-1); j=mini(j+1,nbcol-1);
		if (j>j0) {
			for (k=0; k<nbin; k++) imaIntegrale(0,j,k)=imaIntegrale(0,j-1,k); imaIntegrale(0,j,int(((*this)(0,j)-valmin)*fact))++;
			for (ii=1; ii<i; ii++) {
				for (k=0; k<nbin; k++) imaIntegrale(ii,j,k)=imaIntegrale(ii,j-1,k)+imaIntegrale(ii-1,j,k)-imaIntegrale(ii-1,j-1,k); 
				imaIntegrale(ii,j,int(((*this)(ii,j)-valmin)*fact))++;
			}
		}
		if (i>i0) {
			for (k=0; k<nbin; k++) imaIntegrale(i,0,k)=imaIntegrale(i-1,0,k); imaIntegrale(i,0,int(((*this)(i,0)-valmin)*fact))++;
			for (jj=1; jj<j; jj++) {
				for (k=0; k<nbin; k++) imaIntegrale(i,jj,k)=imaIntegrale(i,jj-1,k)+imaIntegrale(i-1,jj,k)-imaIntegrale(i-1,jj-1,k); 
				imaIntegrale(i,jj,int(((*this)(i,jj)-valmin)*fact))++;
			}
		}
		for (k=0; k<nbin; k++) imaIntegrale(i,j,k)=imaIntegrale(i,j-1,k)+imaIntegrale(i-1,j,k)-imaIntegrale(i-1,j-1,k); 
		imaIntegrale(i,j,int(((*this)(i,j)-valmin)*fact))++;
	} while (i<nblig-1 || j<nbcol-1);
}

template <class T> void imadata<T>::ima_Integrale_Nbin(imadata<unsigned int> &imaIntegrale, T valmax, T valmin) {
	int nbin=imaIntegrale.ncanaux();
	float fact=(float)nbin/(float)(valmax+1-valmin);
	imaIntegrale.mise_a_zero();
	int i=0, j=0, i0, j0, ii, jj, k;
	imaIntegrale(i,j,mini(nbin-1,maxi(0,(int)floor(((*this)(i,j)-valmin)*fact))))++;
	do {
		i0=i; j0=j; i=mini(i+1,nblig-1); j=mini(j+1,nbcol-1);
		if (j>j0) {
			for (k=0; k<nbin; k++) imaIntegrale(0,j,k)=imaIntegrale(0,j-1,k); imaIntegrale(0,j,mini(nbin-1,maxi(0,(int)floor(((*this)(0,j)-valmin)*fact))))++;
			for (ii=1; ii<i; ii++) {
				for (k=0; k<nbin; k++) imaIntegrale(ii,j,k)=imaIntegrale(ii,j-1,k)+imaIntegrale(ii-1,j,k)-imaIntegrale(ii-1,j-1,k); 
				imaIntegrale(ii,j,mini(nbin-1,maxi(0,(int)floor(((*this)(ii,j)-valmin)*fact))))++;
			}
		}
		if (i>i0) {
			for (k=0; k<nbin; k++) imaIntegrale(i,0,k)=imaIntegrale(i-1,0,k); imaIntegrale(i,0,mini(nbin-1,maxi(0,(int)floor(((*this)(i,0)-valmin)*fact))))++;
			for (jj=1; jj<j; jj++) {
				for (k=0; k<nbin; k++) imaIntegrale(i,jj,k)=imaIntegrale(i,jj-1,k)+imaIntegrale(i-1,jj,k)-imaIntegrale(i-1,jj-1,k); 
				imaIntegrale(i,jj,mini(nbin-1,maxi(0,(int)floor(((*this)(i,jj)-valmin)*fact))))++;
			}
		}
		for (k=0; k<nbin; k++) imaIntegrale(i,j,k)=imaIntegrale(i,j-1,k)+imaIntegrale(i-1,j,k)-imaIntegrale(i-1,j-1,k); 
		imaIntegrale(i,j,mini(nbin-1,maxi(0,(int)floor(((*this)(i,j)-valmin)*fact))))++;
	} while (i<nblig-1 || j<nbcol-1);
	imaIntegrale.statbasic();
}

template <class T> int imadata<T>::hist_from_imaIntegrale_16bin(int *Histo16, unsigned short nbin, int i0, int j0, int i2, int j2, bool icumul) const { // attention l'image appellante doit être une image integrale
	if (nbin!=16) {cout<<" attention l'histogramme doit avoir 16 bins\n"; return -1;}
	else {//cout<<" entree dans imadata<T>::hist_from_imaIntegrale_16bin\n";
		int k, nn=0;
		if (!icumul) {for (k=0; k<nbin; k++) Histo16[k]=0;}
		if (i0>0 && j0>0) {
			i0--; j0--;
			for (k=0; k<nbin; k++) Histo16[k]+=(*this)(i2,j2,k)-(*this)(i0,j2,k)-(*this)(i2,j0,k)+(*this)(i0,j0,k);
			nn=(i2-i0)*(j2-j0);
		} else {
			if (i0==0 && j0>0) {
				j0--;
				for (k=0; k<nbin; k++) Histo16[k]+=(*this)(i2,j2,k)-(*this)(i2,j0,k);
				nn=(i2+1)*(j2-j0);
			} else {
				if (i0>0 && j0==0) {
					i0--;
					for (k=0; k<nbin; k++) Histo16[k]+=(*this)(i2,j2,k)-(*this)(i0,j2,k);
					nn=(i2-i0)*(j2+1);
				} else {
					if (i0==0 && j0==0) {
						for (k=0; k<nbin; k++) Histo16[k]+=(*this)(i2,j2,k);
						nn=(i2+1)*(j2+1);
					}
				}
			}
		}
		return nn;
	}
}
template <class T> int imadata<T>::hist_from_imaIntegrale_Nbin(int *HistoI, unsigned short nbin, int i0, int j0, int i2, int j2, bool icumul) const { // attention l'image appellante doit être une image integrale
	if (nbin!=nbcanaux) {cout<<" attention l'histogramme doit avoir le meme nombre de bins qu'il y a de canaux dans l'image integrale\n"; return -1;}
	else {//cout<<" entree dans imadata<T>::hist_from_imaIntegrale_16bin\n";
		int k, nn=0;
		if (!icumul) {for (k=0; k<nbin; k++) HistoI[k]=0;}
		if (i0>0 && j0>0) {
			i0--; j0--;
			for (k=0; k<nbin; k++) HistoI[k]+=(*this)(i2,j2,k)-(*this)(i0,j2,k)-(*this)(i2,j0,k)+(*this)(i0,j0,k);
			nn=(i2-i0)*(j2-j0);
		} else {
			if (i0==0 && j0>0) {
				j0--;
				for (k=0; k<nbin; k++) HistoI[k]+=(*this)(i2,j2,k)-(*this)(i2,j0,k);
				nn=(i2+1)*(j2-j0);
			} else {
				if (i0>0 && j0==0) {
					i0--;
					for (k=0; k<nbin; k++) HistoI[k]+=(*this)(i2,j2,k)-(*this)(i0,j2,k);
					nn=(i2-i0)*(j2+1);
				} else {
					if (i0==0 && j0==0) {
						for (k=0; k<nbin; k++) HistoI[k]+=(*this)(i2,j2,k);
						nn=(i2+1)*(j2+1);
					}
				}
			}
		}
		return nn;
	}
}

template <class T> double imadata<T>::distHisto_Bhattacharyya (imadata<T> &ima2, int nbin) {
	if (ihist && ima2.ihist) 
		if (nbinH!=ima2.nbinH) {if (nbinH!=nbin) ihist=0; if (ima2.nbinH!=nbin) ima2.ihist=0;}
	if (!ihist) histogramme(nbin,0,0);
	if (!ima2.ihist) ima2.histogramme(nbin,0,0);
	int j;
	float n1=0.f, n2=0.f;
	double xx=0.;
	for (j=0; j<nbinH; j++) {
		n1+=histo[j].nbelts; n2+=ima2.histo[j].nbelts;
		xx+=pow((double)histo[j].nbelts*ima2.histo[j].nbelts,0.5);
	}
	xx/=pow((double)n1*n2,0.5);
	return xx;
}

template <class T> void imadata<T>::add_gauss_noise (float v) {
	int i,j,k;
	T x;
	double a,b,z;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++) {
				do {
					a=double(rand())/RAND_MAX;
//					a=drand48();
				} while (a<=0 || a>1);
				do {
					b=double(rand())/RAND_MAX;
//					b=drand48();
				} while (b<0 || b>1);
				z=sqrt(-2.0*(double)v*log(a))*cos(2.0*PI*b); //cout<<z<<" ";
				x=(*this)(i,j,k)+(T)z;
				(*this)(i,j,k)=x;
			}
	statbasic (1);
	}

template <class T> void imadata<T>::add_impul_noise (float p) {
	int i,j,k;
	if (!istat) statbasic ();
	double a,b;
	if (p>1) cout<<" attention p>1 => bruit impulsif > 100%\n";
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			for (k=0; k<nbcanaux; k++) {
//				a=drand48();
				a=double(rand())/RAND_MAX;
				if (a<p) {
//					b=drand48();
					b=double(rand())/RAND_MAX;
//					if (b<0.5) (*this)(i,j,k)=valmin[k];
//					else (*this)(i,j,k)=valmax[k];
					if (b<0.5) (*this)(i,j,k)=0;
					else (*this)(i,j,k)=255;
					}
				}
	statbasic (1);
	}

template <class T> imadata<T> imadata<T>::BR_filtremoyenne (int nl, int nc) const {
	int i,j,k,i0,i2,j0,j2,nl2,nc2,ii,jj,nsum;
	double sum;
	nl2=nl/2;
	nc2=nc/2;
	int ilig=0, icol=0;
	for (i=nl2; i<nblig; i+=nl) ilig++;
	for (j=nc2; j<nbcol; j+=nc) icol++;
	imadata<T> imaRes(ilig,icol,nbcanaux);
	cout<<" image BR : # lignes = "<<ilig<<" # colonnes = "<<icol<<"\n";
	for (k=0; k<nbcanaux; k++) {
		ilig=0;
		for (i=nl2; i<nblig; i+=nl) {
			i0=maxi(0,i-nl2);
//			i2=mini(i+nl2,nblig-1);
			i2=mini(i+nl-nl2-1,nblig-1);
			icol=0;
			for (j=nc2; j<nbcol; j+=nc) {
				j0=maxi(0,j-nc2);
//				j2=mini(j+nc2,nbcol-1);
				j2=mini(j+nc-nc2-1,nbcol-1);
				sum=0.;
				nsum=0;
				for (jj=j0; jj<=j2; jj++)
					for (ii=i0; ii<=i2; ii++) {
						sum+=(*this)(ii,jj,k);
						nsum++;
					}
				imaRes(ilig,icol,k)=(T)(sum/nsum);
				icol++;
			}
			ilig++;
		}
	}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::filtremoyenne(int nl, int nc) const {
	int i,j,k,i0,i2,j0,j2,nl2,nc2,ii,jj,nsum;
	T sum;
	nl2=nl/2;
	nc2=nc/2;
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-nl2);
			i2=mini(i+nl2,nblig-1);
			sum=0;
			nsum=0;
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-nc2);
				j2=mini(j+nc2,nbcol-1);
				if (j>nc2) {                // quand on deplace la fenetre de moyenne
				                                               // on pert une colonne
					for (ii=i0; ii<=i2; ii++) sum-=(*this)(ii,j-nc2-1,k);
					nsum-=i2-i0+1;
				}
				if (j<nbcol-nc2) {                 // on n'est pas au bout de l'image
					if (j==0) {                         // initialisation de la ligne
						for (jj=0; jj<mini(nc2+1,nbcol); jj++) {
						    for (ii=i0; ii<=i2; ii++) sum+=(*this)(ii,jj,k);
						    nsum+=i2-i0+1;
						}
					}
					else {                  // on acquiert une colonne supplémentaire
						for (ii=i0; ii<=i2; ii++) sum+=(*this)(ii,j+nc2,k);
						nsum+=i2-i0+1;
					}
				}
				imaRes(i,j,k)=sum/nsum;
			}
		}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::filtregaussienne(float sigma, bool affich) const {
	int i,j=0,k,i0,i2,j0,j2,ii,jj,dim,dim2;
	float sigma2, fact, norm;
	double sum, coef , center;
	double *noyau=NULL;
	sigma2=2*sigma*sigma;
	fact=1/sqrt(2*(float)PI)/sigma;
	dim2=maxi(int(2*sigma),1);
	dim=2*dim2+1;
	noyau=new double[dim];
	imadata<T> imaRes(*this);
  center = dim/2;
  sum = 0;

/************************** debut a programmer ****************************/
for(ii=0;ii<dim;ii++)      //dim
	{
			noyau[ii]=fact*exp(-(ii-center)*(ii-center)/(sigma2));
			sum = sum + noyau[ii];
	}


	for(ii=0;ii<dim;ii++)   //output
	{
			noyau[ii]/=sum;
			cout<<noyau[ii]<<"  ";
		cout<<endl<<endl;
	}


	for (k=0; k<nbcanaux; k++) // convolution par ligne
		for (i=0;i<nblig; i++)
		{
		if(j<=dim2)
		{
			for (j = 0;j<=dim2; j++)
				  {
					double *tab = NULL;
					double res_conv = 0;
					tab = new double[dim];
					coef = 0;
					for (jj=0; jj<=dim2-j;jj++)
						tab[jj] = 0;
					for (jj = dim-dim2+j;jj<dim; jj++)
					{
						coef = coef + noyau[jj];
						tab[jj] = ((*this)(i,j-dim2+1,k));
						res_conv += noyau[jj]*tab[jj]; 

					}   
					imaRes(i,j,k) = (int) res_conv/sum;
					delete[] tab;
				  }
		}
		else if((j>dim2) && (j < nbcol - dim2 +1))
		{ 
			  for (j = dim2;j<nbcol-dim2-1; j++)
			  {
				double *tab = NULL;
				double res_conv = 0;
				tab = new double[dim];
				coef = 0;

				for (jj = 0;jj<dim; jj++)
				{
					coef = coef + noyau[jj];
					tab[jj] = ((*this)(i,j-dim2+1,k));
					res_conv += noyau[jj]*tab[jj]; 

				}   
				imaRes(i,j,k) = (int) res_conv/sum;
				delete[] tab;
			  }
		}
		else if (j>= nbcol-dim2-1)
		{
		for (j = nbcol-dim2-1;j<dim; j++)
				  {
					double *tab = NULL;
					double res_conv = 0;
					tab = new double[dim];
					coef = 0;
					for (jj=dim-(j+dim2-nbcol); jj<=dim;jj++)
						tab[jj] = 0;
					for (jj = 0;jj<dim-(j+dim2-nbcol); jj++)
					{
						coef = coef + noyau[jj];
						tab[jj] = ((*this)(i,j-dim2+1,k));
						res_conv += noyau[jj]*tab[jj]; 

					}   
					imaRes(i,j,k) = (int) res_conv/sum;
					delete[] tab;
				  }
		}
	}








/***************************  fin a programmer ****************************/

	if (noyau!=NULL) {delete[] noyau; noyau=NULL;}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::filtregaussienne2D (float sigma, bool affich) const {
//template <class T> imadata<T> imadata<T>::filtregaussienne (float sigma, bool affich) {
	int i,j,k,i0,i2,j0,j2,ii,jj,dim,dim2;
	float sigma2, fact, norm;
	double sum, coef;
	double *noyau=NULL;
	sigma2=2*sigma*sigma;
	fact=1/sqrt(2*(float)PI)/sigma;
	dim2=int(2*sigma);
	dim=2*dim2+1;
	noyau=new double[dim*dim];
	imadata<T> imaRes;

/************************** debut a programmer ****************************/




/***************************  fin a programmer ****************************/

	if (noyau!=NULL) {delete[] noyau; noyau=NULL;}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::filtreNagao (int nl, bool imed, BYTE nnoy, bool iaf) const {
	const BYTE nnoymax=11; nnoy=maxi(mini(nnoy,nnoymax),(BYTE)1);
	int nl2=(nl-1)/2, k_nl=around((nl*nl-nl2*nl2)/(float)(2*nl2-1)), dim=maxi(2*nl-1,(nl2+k_nl)*2-1), ic=(dim-1)/2; 
	cout<<" nl = "<<nl<<" nl2 = "<<nl2<<" k_nl = "<<k_nl<<" dim = "<<dim<<" ic = "<<ic<<"\n";
	cout<<" nl*nl = "<<nl*nl<<"; (nl2+1)*(2*nl-1) = "<<(nl2+1)*(2*nl-1)<<"; nl2*nl2+k_nl*(2*nl2-1) = "<<nl2*nl2+k_nl*(2*nl2-1)<<"\n";
	int i,j,k;
// calcul des noyaux
	bool*** T_noyau=new bool**[nnoy];
	for (i=0; i<nnoy; i++) {T_noyau[i]=new bool*[dim]; for (j=0; j<dim; j++) T_noyau[i][j]=new bool[dim]; }
	for (k=0; k<nnoy; k++) for (i=0; i<dim; i++) for (j=0; j<dim; j++) T_noyau[k][i][j]=0; 
	int i0,j0,i2,j2,l; cout<<" nnoy = "<<(int)nnoy<<"\n";
	k=0; i0=j0=ic-nl2; i2=j2=ic+nl2; for (i=i0; i<=i2; i++) for (j=j0; j<=j2; j++) T_noyau[k][i][j]=1; // noyau isotrope
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}
	if (nnoy>++k) {i0=ic-(nl-1); i2=i0+(2*nl-1); j0=ic; j2=j0+nl2; for (i=i0; i<i2; i++) for (j=j0; j<=j2; j++) T_noyau[k][i][j]=1;} // noyau bord ?droite
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {i0=ic-(nl-1); i2=i0+(2*nl-1); j0=ic-nl2; j2=j0+nl2; for (i=i0; i<i2; i++) for (j=j0; j<=j2; j++) T_noyau[k][i][j]=1;} // noyau bord ?gauche
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {j0=ic-(nl-1); j2=j0+(2*nl-1); i0=ic; i2=i0+nl2; for (i=i0; i<=i2; i++) for (j=j0; j<j2; j++) T_noyau[k][i][j]=1;} // noyau bord au-dessus
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {j0=ic-(nl-1); j2=j0+(2*nl-1); i0=ic-nl2; i2=i0+nl2; for (i=i0; i<=i2; i++) for (j=j0; j<j2; j++) T_noyau[k][i][j]=1;} // noyau bord au-dessous
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {i0=ic-(nl2-1)/2; i2=i0+nl2-1; j0=ic-(nl2-1)/2; j2=j0+nl2-1; for (l=-k_nl/2; l<=k_nl/2; l++) for (i=i0-l; i<=i2-l; i++) for (j=j0+l; j<=j2+l; j++) T_noyau[k][i][j]=1;} // noyau 1ère bissectrice up
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {i0=ic-(nl2-1)/2; i2=i0+nl2-1; j0=ic-(nl2-1)/2; j2=j0+nl2-1; for (l=-k_nl/2; l<=k_nl/2; l++) for (i=i0-l; i<=i2-l; i++) for (j=j0-l; j<=j2-l; j++) T_noyau[k][i][j]=1;} // noyau 1ère bissectrice up
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	cout<<(int)nnoy<<" "<<k<<" "<<(int)((int)nnoy>=k)<<"\n";
	if (nnoy>++k) {i0=ic-(nl2-1); i2=i0+nl2-1; j0=ic; j2=j0+nl2-1; for (l=0; l<=k_nl; l++) for (i=i0-l; i<=i2-l; i++) for (j=j0+l; j<=j2+l; j++) T_noyau[k][i][j]=1;} // noyau 1ère bissectrice up
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {i0=ic; i2=i0+nl2-1; j0=ic-(nl2-1); j2=j0+nl2-1; for (l=0; l<=k_nl; l++) for (i=i0+l; i<=i2+l; i++) for (j=j0-l; j<=j2-l; j++) T_noyau[k][i][j]=1;} // noyau 1ère bissectrice down
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {i0=ic-(nl2-1); i2=i0+nl2-1; j0=ic-(nl2-1); j2=j0+nl2-1; for (l=0; l<=k_nl; l++) for (i=i0-l; i<=i2-l; i++) for (j=j0-l; j<=j2-l; j++) T_noyau[k][i][j]=1;} // noyau 2ème bissectrice up
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (nnoy>++k) {i0=ic; i2=i0+nl2-1; j0=ic; j2=j0+nl2-1; for (l=0; l<=k_nl; l++) for (i=i0+l; i<=i2+l; i++) for (j=j0+l; j<=j2+l; j++) T_noyau[k][i][j]=1;} // noyau 2ème bissectrice up
//	cout<<" noyau "<<k<<" ok\n"; for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}}
	if (iaf)
		for (k=0; k<nnoy; k++) {cout<<" ***************** k = "<<k<<" *****************\n";
			for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<"  "<<(int)T_noyau[k][i][j];	cout<<"\n";	}
		}
	imadata<T> imaRes(*this);
	imadata<BYTE> imaNoyau(nblig,nbcol,nbcanaux);
	int dim2=(dim-1)/2, ii, jj, n, m;
	T x, sum, sum2, y;
	if (!imed) {
		double varmin, var, moy;
		for (l=0; l<nbcanaux; l++)
			for (i=0; i<nblig; i++) {
				i0=maxi(i-dim2,0); i2=mini(i+dim2,nblig-1);
				for (j=0; j<nbcol; j++) {
					j0=maxi(j-dim2,0); j2=mini(j+dim2,nbcol-1);
					varmin=DBL_MAX; y=(*this)(i,j,l); m=0;
					for (k=0; k<nnoy; k++) {
						n=0; sum=sum2=(T)0; var=0.;
						for (ii=i0; ii<=i2; ii++)
							for (jj=j0; jj<=j2; jj++)
								if (T_noyau[k][ii-i+dim2][jj-j+dim2]) {n++; x=(*this)(ii,jj,l); sum+=x; sum2+=x*x; }
//								if (T_noyau[k][ii-i+dim2][jj-j+dim2]) {n++; x=(*this)(ii,jj,l); sum+=x; var+=(double)(x-(*this)(i,j,l))*(x-(*this)(i,j,l));}
						if (n>1){ 
							moy=sum/n; var=sum2/(n-1)-moy*moy;
							if (var<varmin) {varmin=var; y=(T)moy; m=k;}
						}
					}
					imaRes(i,j,l)=y; imaNoyau(i,j,l)=m+1;
				}
			}
	}
	if (imed) {
		int nmax=maxi(maxi(nl*nl,(nl2+1)*(2*nl-1)),nl2*nl2+k_nl*(2*nl2-1)); cout<<" nmax = "<<nmax<<"\n";
		double distQmin, distQ;
		double *listeval=new double[nmax+1];
		for (l=0; l<nbcanaux; l++)
			for (i=0; i<nblig; i++) {
				i0=maxi(i-dim2,0); i2=mini(i+dim2,nblig-1);
				for (j=0; j<nbcol; j++) {
					j0=maxi(j-dim2,0); j2=mini(j+dim2,nbcol-1);
					distQmin=DBL_MAX; y=(*this)(i,j,l); m=0;
					for (k=0; k<nnoy; k++) {
						n=0;
						for (ii=i0; ii<=i2; ii++)
							for (jj=j0; jj<=j2; jj++)
								if (T_noyau[k][ii-i+dim2][jj-j+dim2]) {listeval[n++]=(double)(*this)(ii,jj,l); 
																											if (n>nmax) cout<<" *************** PB !!!!!!!!!!!!!!! "<<n<<"\n";}
						if (n>1){ 
							tri_rapide(listeval,n);
//							distQ=abs(listeval[3*n/4]-listeval[n/4]);
							distQ=abs(listeval[n/2]-(*this)(i,j,l));
							if (distQ<distQmin) {distQmin=distQ; y=(T)listeval[n/2]; m=k;} 
						}
					}
					imaRes(i,j,l)=y; imaNoyau(i,j,l)=m+1;
				}
			}
		if (listeval!=NULL) delete[] listeval;
	}
	imaNoyau.imaunsignedchar(1,1).sauve_ImaPGM("imaNoyau.pgm");
	for (i=0; i<nnoy; i++) for (j=0; j<dim; j++) if (T_noyau[i][j]!=NULL) delete[] T_noyau[i][j];
	for (i=0; i<nnoy; i++) if (T_noyau[i]!=NULL) delete[] T_noyau[i];
	if (T_noyau!=NULL) delete[] T_noyau;
	return imaRes;
}

/*template <class T> imadata<T> imadata<T>::filtreNagao () const {
	int dim=5, dim2=2;
	int i,j,k,i0,i2,di0,di2,j0,j2,dj0,dj2,ii,jj,nsum;
	T sum, moy, x;
	
	T** fenetre=new T*[dim];
	for (i=0; i<dim; i++) fenetre[i]=new T[dim];
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	double var, varmin;
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-dim2); i2=mini(i+dim2,nblig-1);
			di0=i0-(i-dim2);
			di2=i+dim2-i2;
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-dim2); j2=mini(j+dim2,nbcol-1);
				dj0=j0-(j-dim2);
				dj2=j+dim2-j2;
				if (j==0) {
					for (ii=0; ii<dim; ii++)
						for (jj=0; jj<dim; jj++) fenetre[ii][jj]=valnul;
					for (ii=i0; ii<=i2; ii++)
						for (jj=j0; jj<=j2; jj++)
							fenetre[ii-i+dim2][jj+dim2]=(*this)(ii,jj,k);
					}
				else {
					for (jj=0; jj<dim-1; jj++)
						for (ii=i0-i+dim2; ii<=i2-i+dim2; ii++)
							fenetre[ii][jj]=fenetre[ii][jj+1];
					if (j<nbcol-dim2)
						for (ii=i0; ii<=i2; ii++)
							fenetre[ii-i+dim2][dim-1]=(*this)(ii,j+dim2,k);
					}
				varmin=1.e+9;
				nsum=1; sum=fenetre[dim2][dim2]; var=(double)sum*(double)sum;
				for (ii=dim2+1; ii<mini(dim,dim-di2); ii++)                        // D3
					for (jj=maxi(dim2-1,dj0); jj<mini(dim-1,dim-dj2); jj++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (nsum>1) { sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} }
				nsum=1; sum=fenetre[dim2][dim2]; var=(double)sum*(double)sum;
				for (ii=maxi(0,di0); ii<dim2; ii++)                                // D2
					for (jj=maxi(dim2-1,dj0); jj<mini(dim-1,dim-dj2); jj++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (nsum>1) { sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} }
				nsum=1; sum=fenetre[dim2][dim2]; var=(double)sum*(double)sum;
				for (jj=maxi(0,dj0); jj<dim2; jj++)                                // D4
					for (ii=maxi(dim2-1,di0); ii<mini(dim-1,dim-di2); ii++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (nsum>1) { sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} }
				nsum=1; sum=fenetre[dim2][dim2]; var=(double)sum*(double)sum;
				for (jj=dim2+1; jj<mini(dim,dim-dj2); jj++)                        // D1
					for (ii=maxi(dim2-1,di0); ii<mini(dim-1,dim-di2); ii++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (nsum>1) { sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} }
				nsum=0; sum=0; var=0.;
				for (jj=maxi(1,dj0); jj<mini(dim-1,dim-dj2); jj++)                // D5
					for (ii=maxi(1,di0); ii<mini(dim-1,dim-di2); ii++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
//				if (nsum>1) {
					sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} //}
				nsum=0; sum=0; var=0.;
				for (jj=dim2; jj<mini(dim,dim-dj2); jj++)                         // D7
					for (ii=dim2; ii<mini(dim,dim-di2); ii++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (dj2==0) {x=fenetre[dim2][dim-1]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
				if (di2==0) {x=fenetre[dim-1][dim2]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
//				if (nsum>1) {
					sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} //}
				nsum=0; sum=0; var=0.;
				for (jj=maxi(0,dj0); jj<=dim2; jj++)                              // D8
					for (ii=dim2; ii<mini(dim,dim-di2); ii++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (dj0==0) {x=fenetre[dim2][0]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
				if (di2==0) {x=fenetre[dim-1][dim2]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
//				if (nsum>1) {
					sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} //}
				nsum=0; sum=0; var=0.;
				for (jj=maxi(0,dj0); jj<=dim2; jj++)                              // D9
					for (ii=maxi(0,di0); ii<=dim2; ii++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (dj0==0) {x=fenetre[dim2][0]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
				if (di0==0) {x=fenetre[0][dim2]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
//				if (nsum>1) {
					sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} //}
				nsum=0; sum=0; var=0.;
				for (jj=dim2; jj<mini(dim,dim-dj2); jj++)                         // D6
					for (ii=maxi(0,di0); ii<=dim2; ii++) {
						x=fenetre[ii][jj]; nsum++;
						sum+=x; var+=(double)x*(double)x;}
				if (dj2==0) {x=fenetre[dim2][dim-1]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
				if (di0==0) {x=fenetre[0][dim2]; nsum--;
					sum-=x; var-=(double)x*(double)x;}
//				if (nsum>1) {
					sum=sum/nsum; var=var/nsum-sum*sum;
					if (var<varmin) {varmin=var; moy=sum;} //}
				imaRes(i,j,k)=moy;
				}
			}
	for (i=0; i<dim; i++) if (fenetre[i]!=NULL) {delete[] fenetre[i]; fenetre[i]=NULL;}
	if (fenetre!=NULL) {delete[] fenetre; fenetre=NULL;}
	return imaRes;
}*/

template <class T> imadata<T> imadata<T>::filtreSymNearNeigh (int nl) const {
	int dim=nl, dim2=(nl-1)/2;
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	int i,j,k,i0,i2,j0,j2,ii,jj,iis,jjs,nsum;
	T x0, sum, x, xOK;
	double dist;
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++) {
			i0=i-dim2; i2=i+dim2;
			for (j=0; j<nbcol; j++) {
				j0=j-dim2; j2=j+dim2;
				x0=(*this)(i,j,k); sum=0; nsum=0;
				for (ii=i0; ii<=i; ii++) {
					iis=i+(i-ii);
					for (jj=j0; jj<=j2; jj++) {
						if ((ii<i) || (jj<j)) {
							jjs=j+(j-jj); dist=DBL_MAX; xOK=valnul;
							if ((ii>=0) && (jj>=0) && (ii<nblig) && (jj<nbcol)) {xOK=(*this)(ii,jj,k); dist=fabs((double)x0-xOK); }
							if ((iis>=0) && (jjs>=0) && (iis<nblig) && (jjs<nbcol)) {x=(*this)(iis,jjs,k); if (fabs((double)x0-x)<dist) xOK=x; }
							if (xOK!=valnul) {sum+=xOK; nsum++;	}
						}
					}
				}
				if (nsum>0) imaRes(i,j,k)=sum/nsum;
			}
		}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::filtremedian (int nl, int nc) const { 
//	statbasic(1);
	int i,j,k,i0,i2,j0,j2,nl2,nc2,n,ii,jj,nn,nlist;
	T *listeval=NULL, x;
	if (nl%2==0) {nl++; cout<<" forcage nb lig fenetre impair -> "<<nl<<"\n";}
	if (nc%2==0) {nc++; cout<<" forcage nb col fenetre impair -> "<<nc<<"\n";}
	nl2=nl/2; nc2=nc/2; nlist=nl*nc;
	listeval=new T[nlist];
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++) {
      i0=maxi(0,i-nl2); i2=mini(i+nl2,nblig-1);
			for (n=0; n<nlist; n++) listeval[n]=0;
			nlist=0;
			for (j=0; j<nbcol; j++) {
        j0=maxi(0,j-nc2); j2=mini(j+nc2,nbcol-1);
				if (j>nc2) {
					for (ii=i0; ii<=i2; ii++) {
						x=(*this)(ii,j-nc2-1,k);
						n=0;
						while (listeval[n]!=x && n<nlist) n++;
						nlist--;
						for (nn=n; nn<nlist; nn++) listeval[nn]=listeval[nn+1];
					}
				}
				if (j<nbcol-nc2) {
					if (j==0) {
						nlist=0;
						for (jj=0; jj<mini(nc2+1,nbcol); jj++)
							for (ii=i0; ii<=i2; ii++) {
								x=(*this)(ii,jj,k);
								n=0;
								while (n<nlist && listeval[n]<x) n++;
								for (nn=nlist-1; nn>=n; nn--) listeval[nn+1]=listeval[nn];
								listeval[n]=x;
								nlist++;
							}
					}
					else
						for (ii=i0; ii<=i2; ii++) {
							x=(*this)(ii,j+nc2,k);
							n=0;
							while (listeval[n]<x && n<nlist) n++;
              for (nn=nlist-1; nn>=n; nn--) listeval[nn+1]=listeval[nn];
							listeval[n]=x;
							nlist++;
						}
				}
				imaRes(i,j,k)=listeval[nlist/2];
			}
		}
	if (listeval!=NULL) {delete[] listeval; listeval=NULL;}
	imaRes.statbasic(1);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::filtremedian1Dadapt (int nl, int ndir=4) const { 
//	statbasic(1);
	int i,j,k,l,nl2,n,ii,jj,nlist,idir;//i0,i2,j0,j2,nc2,nn;
	if (nl%2==0) {nl++; cout<<" forcage nb lig fenetre impair -> "<<nl<<"\n";}
	if (nl<4) ndir=4; cout<<" nbre de directions considerees = "<<ndir<<"\n";
	nl2=nl/2; nlist=nl; 
	T *listeval=new T[nlist], x;
	int **coord_i=new int*[ndir], **coord_j=new int*[ndir];
	for (l=0; l<ndir; l++) { //cout<<"  * dir *  "<<l<<"\n";
		coord_i[l]=new int[nlist]; coord_j[l]=new int[nlist]; 
		switch (l) {
			case 0: for (i=0; i<nlist; i++) coord_i[l][i]=0; for (j=0; j<nlist; j++) coord_j[l][j]=j-nl2; break;
			case 1: for (i=0; i<nlist; i++) coord_i[l][i]=i-nl2; for (j=0; j<nlist; j++) coord_j[l][j]=0; break;
			case 2: for (i=0; i<nlist; i++) coord_i[l][i]=i-nl2; for (j=0; j<nlist; j++) coord_j[l][j]=j-nl2; break;
			case 3: for (i=0; i<nlist; i++) coord_i[l][i]=i-nl2; for (j=0; j<nlist; j++) coord_j[l][j]=nl2-j; break;
			case 4: for (i=0; i<nlist; i++) coord_i[l][i]=(i-nl2)/2; for (j=0; j<nlist; j++) coord_j[l][j]=j-nl2; break;
			case 5: for (i=0; i<nlist; i++) coord_i[l][i]=(i-nl2)/2; for (j=0; j<nlist; j++) coord_j[l][j]=nl2-j; break;
			case 6: for (i=0; i<nlist; i++) coord_i[l][i]=i-nl2; for (j=0; j<nlist; j++) coord_j[l][j]=(j-nl2)/2; break;
			case 7: for (i=0; i<nlist; i++) coord_i[l][i]=i-nl2; for (j=0; j<nlist; j++) coord_j[l][j]=(nl2-j)/2; break;
			default: for (i=0; i<nlist; i++) coord_i[l][i]=0; for (j=0; j<nlist; j++) coord_j[l][j]=0; break;
		}
	}
	for (l=0; l<ndir; l++) {
		cout<<" direction "<<l<<" : "; for (i=0; i<nlist; i++) cout<<"("<<coord_i[l][i]<<","<<coord_j[l][i]<<") "; 
		cout<<"\n";	}
	imadata<T> imaRes(*this);
	float xmoy, xvar, xn, xvarmin;
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) { //cout<<"   *   "<<i<<","<<j;//<<"\n";
				xvarmin=FLT_MAX; idir=-1;
				for (l=0; l<ndir; l++) {
					xmoy=xvar=xn=0.f; 
					for (n=0; n<nl; n++) {
						ii=mini(maxi(0,i+coord_i[l][n]),nblig-1); 
						jj=mini(maxi(0,j+coord_j[l][n]),nbcol-1);
						x=(*this)(ii,jj,k);
						xmoy+=x; xn+=1; xvar+=x*x;
					}
					if (xn>1) {xmoy/=xn; xvar=(xvar-xn*xmoy*xmoy)/(xn-1);}
					if (xvar<xvarmin) {xvarmin=xvar; idir=l;}
				} //cout<<" xmoy = "<<xmoy<<", xvar = "<<xvar<<"\n";
				if (idir>=0 && idir<ndir) {
					for (n=0; n<nlist; n++) {
						ii=mini(maxi(0,i+coord_i[idir][n]),nblig-1); 
						jj=mini(maxi(0,j+coord_j[idir][n]),nbcol-1);
						listeval[n]=(*this)(ii,jj,k);
					}
					//cout<<" avant tri "; 
					tri_par_insertion(listeval,nlist); //cout<<" fin tri\n";
					imaRes(i,j)=listeval[nl2];
				} else cout<<" Pb : pas de direction privilegiee trouvee : idir = "<<idir<<"\n";
			}
	for (i=0; i<ndir; i++) {
		if (coord_i[i]!=NULL) delete[] coord_i[i]; if (coord_j[i]!=NULL) delete[] coord_j[i];}
	if (coord_i!=NULL) delete[] coord_i; if (coord_j!=NULL) delete[] coord_j;
	if (listeval!=NULL) delete[] listeval;
	imaRes.statbasic(1);
	return imaRes;
}

template <class T> void imadata<T>::psnr(const imadata<T> &ima) const {
	double x, sx;
	int i,j,k,n;
	for (k=0; k<nbcanaux; k++) {
		sx=0.;
		n=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				x=(*this)(i,j,k);
				if (x!=valnul) {
					sx+=pow((x-ima(i,j,k)),2);
					n++;
				}
			}
		if (n>0) sx/=n;
		cout<<" canal "<<k<<" mse = "<<sx<<" => psnr = "<<20*log(255.0)-10*log(sx);
	}
}

template <class T> void imadata<T>::filtre_sigma_delta(imadata<T> &ima_t, T delta, imadata<T> *imasigma) {
	int i,j,k,imax=mini(nblig,ima_t.nlig()),jmax=mini(nbcol,ima_t.ncol()),kmax=mini(nbcanaux,ima_t.ncanaux());
	for (k=0; k<kmax; k++)
		for (i=0; i<imax; i++)
			for (j=0; j<jmax; j++) {
//				float x=(*this)(i,j,k);
				if ((*this)(i,j,k)>ima_t(i,j,k)+delta) (*this)(i,j,k)-=delta;
				else if ((*this)(i,j,k)<ima_t(i,j,k)-delta) (*this)(i,j,k)+=delta;
			}
	if (imasigma!=NULL && imasigma->ncanaux()>=kmax && imasigma->nlig()>=imax && imasigma->ncol()>=jmax) {
		for (k=0; k<kmax; k++)
			for (i=0; i<imax; i++)
				for (j=0; j<jmax; j++) {
					if (fabs((*this)(i,j,k)-ima_t(i,j,k))>(*imasigma)(i,j,k)) (*imasigma)(i,j,k)+=1;
					else if (fabs((*this)(i,j,k)-ima_t(i,j,k))<(*imasigma)(i,j,k)) (*imasigma)(i,j,k)-=1;
				}
	}
}

template <class T> void imadata<T>::filtre_sigma_delta(imadata<T> &ima_t, imabin imasq, T delta, imadata<T> *imasigma) {
	int i,j,k,imax=mini(nblig,ima_t.nlig()),jmax=mini(nbcol,ima_t.ncol()),kmax=mini(nbcanaux,ima_t.ncanaux());
	for (i=0; i<imax; i++)
		for (j=0; j<jmax; j++)
			if (imasq(i,j)) {
				for (k=0; k<kmax; k++) {
					if ((*this)(i,j,k)>ima_t(i,j,k)+delta) (*this)(i,j,k)-=delta;
					else if ((*this)(i,j,k)<ima_t(i,j,k)-delta) (*this)(i,j,k)+=delta;
				}
			}
	if (imasigma!=NULL && imasigma->ncanaux()>=kmax && imasigma->nlig()>=imax && imasigma->ncol()>=jmax) {
		for (i=0; i<imax; i++)
			for (j=0; j<jmax; j++)
				if (imasq(i,j)) {
					for (k=0; k<kmax; k++) {
						if (fabs((*this)(i,j,k)-ima_t(i,j,k))>(*imasigma)(i,j,k)) (*imasigma)(i,j,k)+=1;
						else if (fabs((*this)(i,j,k)-ima_t(i,j,k))<(*imasigma)(i,j,k)) (*imasigma)(i,j,k)-=1;
					}
				}
	}
}

template <class T> void imadata<T>::filtre_sigma_delta_codebook(imadata<T> &ima_t, imadata<T> &imasigma, imadata<BYTE> &imalast_t, T delta) { // l'image ?t est scalaire
	int i,j,k,imax=mini(nblig,ima_t.nlig()),jmax=mini(nbcol,ima_t.ncol()),kdict=nbcanaux-1,k_argmin;
	const int gamma=2, delai_min=5/*10*/;
	float d_min, d_k, xx;
	for (i=0; i<imax; i++)
		for (j=0; j<jmax; j++) {
			xx=ima_t(i,j); k_argmin=-1; d_min=FLT_MAX;
			for (k=0; k<kdict; k++) {                 // recherche de la valeur la plus proche parmi le dictionaire
				if ((*this)(i,j,k)!=valnul) {
					d_k=fabs(xx-(*this)(i,j,k));
					if (d_k<d_min) {k_argmin=k; d_min=d_k;}
				}
			}
			for (k=0; k<kdict; k++) imalast_t(i,j,k)++;
//			if (d_min<gamma*(float)imasigma(i,j,k_argmin)) { // il existe 1 element du dictionaire qui correspond
			if (d_min<gamma*maxi((float)imasigma(i,j,k_argmin),(float)delai_min)) { // il existe 1 element du dictionaire qui correspond
				imalast_t(i,j,k_argmin)=0;
				if ((*this)(i,j,k_argmin)>xx+delta) (*this)(i,j,k_argmin)-=delta;
				else if ((*this)(i,j,k_argmin)<xx-delta) (*this)(i,j,k_argmin)+=delta;
				if (d_min>imasigma(i,j,k_argmin)) imasigma(i,j,k_argmin)+=1;
				else if (d_min<imasigma(i,j,k_argmin)) imasigma(i,j,k_argmin)-=1;
			} else {                                  // il n'existe pas d'element du dictionaire qui corresponde
				if (imalast_t(i,j,kdict)==0) {(*this)(i,j,kdict)=xx;}
				else {
					d_k=fabs(xx-(*this)(i,j,kdict));
//					if (d_k<gamma*(float)imasigma(i,j,kdict)) {    // actualisation de l'element d'apprentissage en cours
					if (d_k<gamma*maxi((float)imasigma(i,j,kdict),(float)delai_min)) {    // actualisation de l'element d'apprentissage en cours
						if ((*this)(i,j,kdict)>xx+delta) (*this)(i,j,kdict)-=delta;
						else if ((*this)(i,j,kdict)<xx-delta) (*this)(i,j,kdict)+=delta;
						if (d_k>imasigma(i,j,kdict)) imasigma(i,j,kdict)+=1;
						else if (d_k<imasigma(i,j,kdict)) imasigma(i,j,kdict)-=1;
					} else {                               // reinitialisation de l'element d'apprentissage
						imalast_t(i,j,kdict)=0; (*this)(i,j,kdict)=xx;
					}
				}
				imalast_t(i,j,kdict)++;
			}
			if (imalast_t(i,j,kdict)>delai_min) { //cout<<" mise a jour du dictionaire en "<<i<<" "<<j<<" "<<xx<<" "<<(*this)(i,j,0)<<" "<<(*this)(i,j,1)<<" "<<(*this)(i,j,kdict)<<"\n";// mise ?jour du dictionnaire
				k_argmin=0; d_k=imalast_t(i,j,0);
				for (k=1; k<kdict; k++)            // recherche de l'element non observe le + longtemps
					if (imalast_t(i,j,k)>d_k) {k_argmin=k; d_k=imalast_t(i,j,k);}
				imalast_t(i,j,k_argmin)=imalast_t(i,j,kdict)=0;
				(*this)(i,j,k_argmin)=(*this)(i,j,kdict);
				imasigma(i,j,k_argmin)=imasigma(i,j,kdict); 
				(*this)(i,j,kdict)=imasigma(i,j,kdict)=0;
			}
			for (k=1; k<kdict; k++)            // elimination des elements non observes
				if ((*this)(i,j,k)!=valnul && imalast_t(i,j,k)>delai_min*10) {(*this)(i,j,k)=valnul; imasigma(i,j,k)=0;}

		}
}

template <class T> imadata<float> imadata<T>::transformee_Fourier(bool iaf) const {
	int		i,j,u,v,nbligTF=nblig,nbcolTF=nbcol;
	const double PI_2=2*acos(-1.), xnli=PI_2/nbligTF, xnco=PI_2/nbcolTF, fact_TFlig=sqrt((double)nblig), fact_TFcol=sqrt((double)nbcol);
	const int yc=nbcolTF/2, xc=nbligTF/2;
	double z,xx,yy,cosz,sinz;
	if (nblig==nbcol && puissance_de_2(nblig)) {
		fourier objetfourier;
		matrice2D<double> mat2D=conv2mat2D(), TFre(nblig,nbcol), TFim(nblig,nbcol);
//	objetfourier.tfd2d (mat2D,TFre,TFim,(unsigned int)nbcol,(unsigned int)nblig);
		objetfourier.fft2d(mat2D,TFre,TFim,(unsigned int)nbcol,(unsigned int)nblig);
		imadata<double> imaFourier(nblig,nbcol,3); 
		imaFourier.copiecanal(0,imadata<double>(TFre)); imaFourier.copiecanal(1,imadata<double>(TFim));
		for (i=0; i<imaFourier.nlig(); i++)
			for (j=0; j<imaFourier.ncol(); j++) imaFourier(i,j,2)=pow(pow(imaFourier(i,j,0),2)+pow(imaFourier(i,j,1),2),0.5);
		imaFourier.statbasic(iaf);
		(imadata<float>(imaFourier)).imaunsignedchar().sauve_ImaPGM("ImaFFT2.ppm");
		imadata<float> logTF_2(nblig,nbcol); logTF_2.mise_a_zero();
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) if (imaFourier(i,j,2)>0) logTF_2(i,j)=(float)log(imaFourier(i,j,2)); logTF_2.imaunsignedchar(1).sauve_ImaPGM("./logFFT2.pgm"); 
		return imaFourier;
	} 
	else {
		imadata<float> imaFourier(nbligTF,nbcolTF,3);
		if (iaf) cout<<" image de la transformee de Fourier : "<<imaFourier.nlig()<<" lignes et "<<imaFourier.ncol()<<" colonnes\n";
		imaFourier.mise_a_zero ();
		imadata<double> imaF_1D(nbligTF,nbcolTF,2);
		for (j=0; j<nbcolTF; j++) {
			for (i=0; i<nbligTF; i++) {
				xx=yy=0;
				for (u=0; u<nblig; u++) {
					z=-u*i*xnli;
					xx+=(*this)(u,j)*cos(z);
					yy+=(*this)(u,j)*sin(z);
				}
				imaF_1D(i,j,0)=xx/fact_TFlig; imaF_1D(i,j,1)=yy/fact_TFlig;
			}
		}
		for (i=0; i<nbligTF; i++) {
			for (j=0; j<nbcolTF; j++) {
				xx=yy=0;
				for (v=0; v<nbcol; v++) {
					z=-v*j*xnco; cosz=cos(z); sinz=sin(z);
					xx+=imaF_1D(i,v,0)*cosz-imaF_1D(i,v,1)*sinz;
					yy+=imaF_1D(i,v,0)*sinz+imaF_1D(i,v,1)*cosz;
				}
				imaFourier(i,j,0)=(float)(xx/fact_TFcol); imaFourier(i,j,1)=(float)(yy/fact_TFcol);
			}
		}
		for (i=0; i<nbligTF/2; i++)
			for (j=0; j<nbcolTF/2; j++) {
				xx=imaFourier(i,j,0); yy=imaFourier(i,j,1);
				imaFourier(i,j,0)=imaFourier(i+xc,j+yc,0); imaFourier(i,j,1)=imaFourier(i+xc,j+yc,1);
				imaFourier(i+xc,j+yc,0)=(float)xx; imaFourier(i+xc,j+yc,1)=(float)yy;
				xx=imaFourier(i+xc,j,0); yy=imaFourier(i+xc,j,1);
				imaFourier(i+xc,j,0)=imaFourier(i,j+yc,0); imaFourier(i+xc,j,1)=imaFourier(i,j+yc,1);
				imaFourier(i,j+yc,0)=(float)xx; imaFourier(i,j+yc,1)=(float)yy;
			}
		const double eps=1.e-9;
		for (i=1; i<nbligTF/2; i++)
			for (j=1; j<nbcolTF/2; j++)
				if (abs(imaFourier(i,j,0)-imaFourier(nbligTF-i,nbcolTF-j,0))>eps || (imaFourier(i,j,1)+imaFourier(nbligTF-i,nbcolTF-j,1))>eps) {
//					cout<<" TF non symetrique / centre : en "<<i<<" "<<j<<" partie reelle = "<<imaFourier(i,j,0)<<" "<<imaFourier(nbligTF-i,nbcolTF-j,0)<<", partie imaginaire = "<<imaFourier(i,j,1)<<" "<<imaFourier(nbligTF-i,nbcolTF-j,1)<<"\n";
					//char aa; cin>>aa;
				}
		for (i=0; i<imaFourier.nlig(); i++)
			for (j=0; j<imaFourier.ncol(); j++) imaFourier(i,j,2)=pow(pow(imaFourier(i,j,0),2)+pow(imaFourier(i,j,1),2),0.5f);
		imaFourier.statbasic(iaf); imaFourier.imaunsignedchar().sauve_ImaPGM("ImaTF2D.ppm");
		imadata<float> logTF_1(nbligTF,nbcolTF); logTF_1.mise_a_zero();
		for (i=0; i<nbligTF; i++) for (j=0; j<nbcolTF; j++) if (imaFourier(i,j,2)>0) logTF_1(i,j)=log(imaFourier(i,j,2)); logTF_1.imaunsignedchar(1).sauve_ImaPGM("./logTF2D.pgm"); 
		return imaFourier;
	}
}

template <class T> imadata<float> imadata<T>::transformee_Fourier_Inv(bool iaf) const {
	const int yc=nbcol/2, xc=nblig/2;
	const double PI_2=2.*acos(-1.), xnli=PI_2/nblig, xnco=PI_2/nbcol;
  int i,j,u,v;
	double z,xx,yy,cosz,sinz;
	imadata<float> imaFourier(*this), imaFourierInv(nblig,nbcol,3);
	if (iaf) cout<<" image de la transformee de Fourier inverse : "<<imaFourierInv.nlig()<<" lignes et "<<imaFourierInv.ncol()<<" colonnes\n";
	imaFourierInv.mise_a_zero ();
	for (i=0; i<nblig/2; i++)
		for (j=0; j<nbcol/2; j++) {
			xx=imaFourier(i,j,0); yy=imaFourier(i,j,1);
			imaFourier(i,j,0)=imaFourier(i+xc,j+yc,0); imaFourier(i,j,1)=imaFourier(i+xc,j+yc,1);
			imaFourier(i+xc,j+yc,0)=(float)xx; imaFourier(i+xc,j+yc,1)=(float)yy;
			xx=imaFourier(i+xc,j,0); yy=imaFourier(i+xc,j,1);
			imaFourier(i+xc,j,0)=imaFourier(i,j+yc,0); imaFourier(i+xc,j,1)=imaFourier(i,j+yc,1);
			imaFourier(i,j+yc,0)=(float)xx; imaFourier(i,j+yc,1)=(float)yy;
		}
	imadata<double> imaF_1D(nblig,nbcol,2);
	for (j=0; j<nbcol; j++) {
		for (i=0; i<nblig; i++) {
			xx=0; yy=0;
			for (u=0; u<nblig; u++) {
				z=u*i*xnli; cosz=cos(z); sinz=sin(z);
				xx+=imaFourier(u,j,0)*cosz-imaFourier(u,j,1)*sinz;
				yy+=imaFourier(u,j,0)*sinz+imaFourier(u,j,1)*cosz;
			}
			imaF_1D(i,j,0)=xx;
			imaF_1D(i,j,1)=yy;
		}
	}
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			xx=0; yy=0;
			for (v=0; v<nbcol; v++) {
				z=v*j*xnco; cosz=cos(z); sinz=sin(z);
				xx+=imaF_1D(i,v,0)*cosz-imaF_1D(i,v,1)*sinz;
				yy+=imaF_1D(i,v,0)*sinz+imaF_1D(i,v,1)*cosz;
			}
			imaFourierInv(i,j,0)=(float)xx;
			imaFourierInv(i,j,1)=(float)yy;
		}
	}
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			imaFourierInv(i,j,0)=imaFourierInv(i,j,0)/nblig/nbcol;
			imaFourierInv(i,j,1)=imaFourierInv(i,j,1)/nblig/nbcol;
		}
	for (i=0; i<imaFourierInv.nlig(); i++)
		for (j=0; j<imaFourierInv.ncol(); j++) imaFourierInv(i,j,2)=pow(pow(imaFourierInv(i,j,0),2)+pow(imaFourierInv(i,j,1),2),0.5f);
	imaFourierInv.statbasic(iaf);
	return imaFourierInv;
}

template <class T> float imadata<T>::boite_englobee(int &lmin1, int &lmax1, int &cmin1, int &cmax1) const {
	int i,j,lmin,lmax,cmin,cmax,ii,jj;
	const float eps=(float)1.e-6;
	imadata<float> ima1=filtremedian(3,3); //ima1.imaunsignedchar().sauve_ImaPGM("ima1.pgm"); {char aa; cin>>aa;}
	imabin ima1b(ima1,eps);
	bool ligneOk;
	float a_max=0.f;
	lmin1=-1; lmax1=nblig; cmin1=-1; cmax1=nbcol;
	for (i=0; i<nblig; i++) { //cout<<i<<" "<<lmin1<<" "<<cmin1<<" "<<lmax1<<" "<<cmax1<<"\n";
		lmin=lmax=i; cmin=-1; cmax=nbcol;
		for (j=0; j<nbcol; j++) {
			if (cmin<0 && ima1b(i,j)) cmin=j;
			if (cmin>=0 && !ima1b(i,j)) {cmax=j; j=nbcol;}
		}
		if (cmin>=0) {
			ii=mini(i+1,nblig-1); ligneOk=1;
			do {
				for (jj=cmin; jj<cmax; jj++) if (!ima1b(ii,jj)) ligneOk=0;
				if (ligneOk) lmax=ii++;
			} while (ligneOk && ii<nblig);
			if ((cmax-cmin)*(lmax-lmin)>a_max) {a_max=(float)(cmax-cmin)*(lmax-lmin); lmin1=lmin; lmax1=lmax+1; cmin1=cmin; cmax1=cmax; }
		}
	}
	return a_max;
}

template <class T> imadata<float> imadata<T>::filtrage_Hamming () const {
	int i,j,k/*,k2*/;
	double coef;
	int nl_I=nblig, nc_I=nbcol, nk_I=nbcanaux, i0_1=0, j0_1=0, i0_2=0, j0_2=0;
	imadata<float> im_1(nl_I,nc_I,nk_I);
	matrice2D<float> M_FH_l(nl_I,1), M_FH_c(1,nc_I); 
	coef=2.*PI/(nl_I-1); for (i=0; i<nl_I; i++) M_FH_l(i,0)=(float)(0.54-0.46*cos(i*coef)); coef=2.*PI/(nc_I-1); for (j=0; j<nc_I; j++) M_FH_c(0,j)=(float)(0.54-0.46*cos(j*coef));
	for (i=0; i<nl_I; i++) 
		for (j=0; j<nc_I; j++) 
			for (k=0; k<nk_I; k++) im_1(i,j,k)=(*this)(i+i0_1,j+j0_1,k)*M_FH_l(i,0)*M_FH_c(0,j);
	return im_1;
}

template <class T> imadata<float> imadata<T>::correlation_phase(const imadata<T> &ima2, int &imax, int &jmax, float &corrmax, bool if_Hamming, bool ielim_noir) const {
	int i,j,lmin1=-1,lmax1=nblig,cmin1=-1,cmax1=nbcol,lmin2=-1,lmax2=ima2.nlig(),cmin2=-1,cmax2=ima2.ncol();
	if (ielim_noir) {
		float a_max1=boite_englobee(lmin1,lmax1,cmin1,cmax1), a_max2=ima2.boite_englobee(lmin2,lmax2,cmin2,cmax2);
		cout<<" boite englobee dans image 1 = ("<<lmin1<<","<<cmin1<<") -> ("<<lmax1<<","<<cmax1<<") => aire = "<<a_max1<<" pixels\n";
		cout<<" boite englobee dans image 2 = ("<<lmin2<<","<<cmin2<<") -> ("<<lmax2<<","<<cmax2<<") => aire = "<<a_max2<<" pixels\n";
	} else {
		lmin1=cmin1=lmin2=cmin2=0; lmax1=nblig; cmax1=nbcol; lmax2=ima2.nlig(); cmax2=ima2.ncol();
	}
	imadata<float> ima_1(*this), ima_2(ima2);
	if (if_Hamming) {ima_1=ima_1.filtrage_Hamming(); ima_1.imaunsignedchar(1).sauve_ImaPGM("./ima1_PH.pgm"); ima_2=ima_2.filtrage_Hamming(); ima_2.imaunsignedchar(1).sauve_ImaPGM("./ima2_PH.pgm");}
	imadata<float> TF_1, TF_2;
	if (lmin1==0 && cmin1==0 && lmax1==nblig && cmax1==nbcol) TF_1=ima_1.transformee_Fourier();
	else {
		imadata<float> ima1bis(lmax1-lmin1,cmax1-cmin1);
		for (i=lmin1; i<lmax1; i++) for (j=cmin1; j<cmax1; j++) ima1bis(i-lmin1,j-cmin1)=ima_1(i,j); ima1bis.imaunsignedchar().sauve_ImaPGM("ima1bis.pgm");
		TF_1=ima1bis.transformee_Fourier();
	}
	if (lmin2==0 && cmin2==0 && lmax2==ima2.nlig() && cmax2==ima2.ncol()) TF_2=ima_2.transformee_Fourier();
	else {
		imadata<float> ima2bis(lmax2-lmin2,cmax2-cmin2);
		for (i=lmin2; i<lmax2; i++) for (j=cmin2; j<cmax2; j++) ima2bis(i-lmin2,j-cmin2)=ima_2(i,j); ima2bis.imaunsignedchar().sauve_ImaPGM("ima2bis.pgm");
		TF_2=ima2bis.transformee_Fourier();
	}
//	imadata<float> TF_1=transformee_Fourier(), TF_2=ima2.transformee_Fourier(); //cout<<" fin TF dans correlation de phase "; {char aa; cin>>aa;}
	int nl_TF=mini(TF_1.nlig(),TF_2.nlig()), nc_TF=mini(TF_1.ncol(),TF_2.ncol());
	imadata<float> ima_corrphase(nl_TF,nc_TF,3); ima_corrphase.mise_a_zero();
	float z;
	for (i=0; i<nl_TF; i++)
		for (j=0; j<nc_TF; j++) {
			ima_corrphase(i,j,0)=(TF_1(i,j,0)*TF_2(i,j,0)+TF_1(i,j,1)*TF_2(i,j,1));
			ima_corrphase(i,j,1)=(TF_1(i,j,1)*TF_2(i,j,0)-TF_1(i,j,0)*TF_2(i,j,1));
			z=TF_1(i,j,2)*TF_2(i,j,2);
			ima_corrphase(i,j,2)=pow(pow(ima_corrphase(i,j,0),2)+pow(ima_corrphase(i,j,1),2),0.5f)/z;
			ima_corrphase(i,j,0)/=z;
			ima_corrphase(i,j,1)/=z;
		}
  ima_corrphase=ima_corrphase.transformee_Fourier_Inv(); 
	cout<<" apres correlation de pharse image de la TF inv de dim "<<ima_corrphase.nlig()<<" par "<<ima_corrphase.ncol()<<"\n"; ima_corrphase.sauve_ImaBSQ("ima_corrphase.dat");
	float zmax=ima_corrphase(0,0,2), zmax2=-999.;
	int imax2=-1,jmax2=-1;
	imax=jmax=0;
	for (i=0; i<nl_TF; i++) {
		for (j=0; j<nc_TF; j++) {
			z=ima_corrphase(i,j,2);
			if (z>zmax) {
				zmax2=zmax; imax2=imax; jmax2=jmax;
				zmax=z; imax=i; jmax=j;
			}
			else if (z>zmax2) {zmax2=z; imax2=i; jmax2=j;}
		}
	}
	cout<<" 1er max = "<<zmax<<" en ("<<imax<<","<<jmax<<"), 2eme max = "<<zmax2<<" en ("<<imax2<<","<<jmax2<<")\n";
//	imax-=nl_TF/2; imax2-=nl_TF/2; jmax-=nc_TF/2; jmax2-=nc_TF/2;                 // inutile si la TF inverse n'est pas centree
/*	if (imax>nl_TF/2) imax=nl_TF-imax; if (imax2>nl_TF/2) imax2=nl_TF-imax2;
	if (jmax>nc_TF/2) jmax=nc_TF-jmax; if (jmax2>nc_TF/2) jmax2=nc_TF-jmax2;*/
	imax=-imax; jmax=-jmax; imax2=-imax2; jmax2=-jmax2;
	if (imax<-nl_TF/2) imax+=nl_TF; if (imax2<-nl_TF/2) imax2+=nl_TF;
	if (jmax<-nc_TF/2) jmax+=nc_TF; if (jmax2<-nc_TF/2) jmax2+=nc_TF;
	corrmax=zmax;
	cout<<" ***** MAXIMUM = "<<corrmax<<" en ("<<imax<<","<<jmax<<") ***** \n"; ima_corrphase.sauve_Ima("ima_corrphase.dat");
	cout<<" ***** En tenant compte du fait que TF calculees sur des sous-image(s), MAXIMUM en ("<<imax+lmin2-lmin1<<","<<jmax+cmin2-cmin1<<") ***** \n";
	imax=imax-lmin2+lmin1; jmax=jmax-cmin2+cmin1;
/*	const int npos=5;
	int k,k2;
	val_pos Tzmax[npos]; for (k=0; k<npos; k++) {Tzmax[k].val=-999; Tzmax[k].pos=0;}
	bool maxloc, trouve; cout<<" nl_TF = "<<nl_TF<<", nc_TF = "<<nc_TF<<"\n";
	for (i=0; i<nl_TF; i++)
		for (j=0; j<nc_TF; j++) {
			z=ima_corrphase(i,j); 
			maxloc=trouve=0; k=0;
			while (!maxloc && k<npos) {
				if (z>Tzmax[k].val) {
					for (k2=npos-1; k2>k; k2--) {Tzmax[k2].val=Tzmax[k2-1].val; Tzmax[k2].pos=Tzmax[k2-1].pos; }
					Tzmax[k].val=z; Tzmax[k].pos=i*nc_TF+j; maxloc=1;
				}
				k++;
			}
		}
	z=(float)Tzmax[0].val; cout<<" max a "<<z<<" -> seuil de significativite des max a "<<z-log(2.f)<<"\n"; z-=log(2.f);
	for (k=0; k<npos; k++) if (Tzmax[k].val<z) Tzmax[k].val=-999.;
	imadata<float> im_2_to_1;
	float corrmaxN=0.,corr,var1,var2;
	int nblig2,nbcol2,n,dx=imax,dy=jmax,ii,jj;
	for (k=0; k<npos; k++) 
		if (Tzmax[k].val>-999.) {
			dx=-(int)(Tzmax[k].pos/nc_TF); dy=-(int)(Tzmax[k].pos%nc_TF);
			if (dx<-nl_TF/2) dx+=nl_TF; if (dy<-nc_TF/2) dy+=nc_TF;
			cout<<"\n "<<k+1<<"eme max = "<<Tzmax[k].val<<" en ("<<Tzmax[k].pos/nc_TF<<","<<Tzmax[k].pos%nc_TF<<") -> dx = "<<dx<<" & dy = "<<dy<<"\n";
			im_2_to_1=ima_2.projette_ima(-dx,-dy,1,0);
			corr=var1=var2=0.; n=0;
			nblig2=im_2_to_1.nlig(); nbcol2=im_2_to_1.ncol();
			for (i=0; i<mini(ima_1.nlig(),nblig2); i++)
				for (j=0; j<mini(ima_1.ncol(),nbcol2); j++) 
					if (im_2_to_1(i,j)>0 && ima_1(i,j)>0) {corr+=im_2_to_1(i,j)*ima_1(i,j); var1+=ima_1(i,j)*ima_1(i,j); var2+=im_2_to_1(i,j)*im_2_to_1(i,j); n++;}
			if (n>0) {
				var1/=n; var2/=n; corr/=n; corr/=pow(var1*var2,0.5f);
				cout<<" Translation testee entre les 2 images : ("<<dx<<","<<dy<<") => correlation = "<<corr<<" >? corr_max"<<corrmaxN<<"\n";
				if (corr>corrmaxN) {corrmaxN=corr; imax=-dx; jmax=-dy;}
			}
		}*/
	return ima_corrphase;
}

template <class T> imadata<float> imadata<T>::transformee_FourierMellin (const imadata<T> &ima2, float &_dx, float &_dy, float &_rapport, float &_angle, bool i_fft) const {
	int i,j,k,k2;
	double coef;
	int nl_I=mini(nblig,ima2.nlig()), nc_I=mini(nbcol,ima2.ncol()), nk_I=mini(nbcanaux,ima2.ncanaux()), i0_1=0, j0_1=0, i0_2=0, j0_2=0;
//	int nl_I=maxi(nblig,ima2.nlig()), nc_I=maxi(nbcol,ima2.ncol()), nk_I=maxi(nbcanaux,ima2.ncanaux()), i0_1=0, j0_1=0, i0_2=0, j0_2=0;
	if (i_fft) {
		if (!puissance_de_2(nl_I) || !puissance_de_2(nc_I) || nl_I!=nc_I) {
			nl_I=nc_I=(int)pow(2,floor(log((double)mini(nl_I,nc_I))/log(2.))); 
			i0_1=(nblig-nl_I)/2; j0_1=(nbcol-nc_I)/2; i0_2=(ima2.nlig()-nl_I)/2; j0_2=(ima2.ncol()-nc_I)/2;
			cout<<" initialement "<<mini(nblig,ima2.nlig())<<" lig. & "<<mini(nbcol,ima2.ncol())<<" col. => "<<nl_I<<" lig. & "<<nc_I<<" col. pour avoir image carree en puissance de 2\n";
			cout<<" offset sur l'image 1 = ("<<i0_1<<","<<j0_1<<") & sur l'image 2 = ("<<i0_2<<","<<j0_2<<")\n";
		}
	}
	imadata<float> im_1(nl_I,nc_I,nk_I), im_2(nl_I,nc_I,nk_I); im_1.mise_a_zero(); im_2.mise_a_zero();
	matrice2D<float> M_FH_l(nl_I,1), M_FH_c(1,nc_I); 
	coef=2.*PI/(nl_I-1); for (i=0; i<nl_I; i++) M_FH_l(i,0)=(float)(0.54-0.46*cos(i*coef)); coef=2.*PI/(nc_I-1); for (j=0; j<nc_I; j++) M_FH_c(0,j)=(float)(0.54-0.46*cos(j*coef));
	for (i=0; i<nl_I; i++) 
		for (j=0; j<nc_I; j++) 
			for (k=0; k<nk_I; k++) {
				im_1(i,j,k)=(*this)((i+i0_1)%nblig,(j+j0_1)%nbcol,k)*M_FH_l(i,0)*M_FH_c(0,j); 
				im_2(i,j,k)=ima2((i+i0_2)%ima2.nlig(),(j+j0_2)%ima2.ncol(),k)*M_FH_l(i,0)*M_FH_c(0,j);
			}
	im_1.imaunsignedchar(1).sauve_ImaPGM("./ima1_PH.pgm"); im_2.imaunsignedchar(1).sauve_ImaPGM("./ima2_PH.pgm");

	imadata<float> TF_1=im_1.transformee_Fourier(); imadata<float> TF_2=im_2.transformee_Fourier(); 
	TF_1.sauve_ImaPGM("./TF_1.pgm"); TF_2.sauve_ImaPGM("./TF_2.pgm");
	int nl_TF=mini(TF_1.nlig(),TF_2.nlig()), nc_TF=mini(TF_1.ncol(),TF_2.ncol());
	imadata<float> logTF_1(nl_TF,nc_TF), logTF_2(nl_TF,nc_TF); logTF_1.mise_a_zero(); logTF_2.mise_a_zero();
	for (i=0; i<nl_TF; i++) for (j=0; j<nc_TF; j++) if (TF_1(i,j,2)>0) logTF_1(i,j)=pow(TF_1(i,j,2),2.f); //logTF_1.imaunsignedchar(1).sauve_ImaPGM("./logTF_1.pgm"); 
	for (i=0; i<nl_TF; i++) for (j=0; j<nc_TF; j++) if (TF_2(i,j,2)>0) logTF_2(i,j)=pow(TF_2(i,j,2),2.f); //logTF_2.imaunsignedchar(1).sauve_ImaPGM("./logTF_2.pgm"); 
	matrice2D<float> Mcos_l(nl_TF,1), Mcos_c(1,nc_TF); 
	coef=PI/(nl_TF-1); for (i=0; i<nl_TF; i++) Mcos_l(i,0)=(float)cos(-PI/2.+i*coef); coef=PI/(nc_TF-1); for (j=0; j<nc_TF; j++) Mcos_c(0,j)=(float)cos(-PI/2.+j*coef);
	matrice2D<float> M_X(Mcos_l,Mcos_c,nl_TF,1,nc_TF), M_PH(nl_TF,nc_TF);
	for (i=0; i<nl_TF; i++) for (j=0; j<nc_TF; j++) M_PH(i,j)=(1.f-M_X(i,j))*(2.f-M_X(i,j)); imadata<float>(M_PH).imaunsignedchar(1).sauve_ImaPGM("./filtrePH.pgm");
//	for (i=0; i<nl_TF; i++) for (j=0; j<nc_TF; j++) M_PH(i,j)=(1.f-Mcos_l(i,0)*Mcos_c(0,j))*(2.f-Mcos_l(i,0)*Mcos_c(0,j)); imadata<float>(M_PH).imaunsignedchar(1).sauve_ImaPGM("./filtrePH.pgm");
	for (i=0; i<nl_TF; i++) for (j=0; j<nc_TF; j++) {logTF_1(i,j)*=M_PH(i,j); logTF_2(i,j)*=M_PH(i,j);} logTF_1.imaunsignedchar(1).sauve_ImaPGM("./logTF_1ph.pgm"); logTF_2.imaunsignedchar(1).sauve_ImaPGM("./logTF_2ph.pgm");
	int nlnRho=200,nTheta=360,iu,jv,imax=0,jmax=0;
//	int nlnRho=256,nTheta=256,iu,jv,imax=0,jmax=0;
	imadata<float> TF_1_logpol(nlnRho,nTheta), TF_2_logpol(nlnRho,nTheta); TF_1_logpol.mise_a_zero(); TF_2_logpol.mise_a_zero();
	float rhomax=(float)mini(nl_TF/2,nc_TF/2),rhomin=12.f,coef_lnRho=log((float)rhomax/(float)rhomin)/(nlnRho-1), 
				thetamax=359.f,thetamin=0.f,pas_theta=(thetamax-thetamin)/(nTheta-1);
	double theta, rho, x, y;
	bool inter_pol=1;
	if (inter_pol) {
		double d1, d2, d3, d4;
		for (i=0; i<nlnRho; i++) {
			rho=rhomin*exp(coef_lnRho*i);
			for (j=0; j<nTheta; j++) {
//				theta=(double)j*360./nTheta*PI/180.; 
				theta=(double)(j*pas_theta+thetamin)*PI/180.; 
				y=rho*cos(theta)+nl_TF/2.; x=rho*sin(theta)+nc_TF/2.;
				iu=(int)floor(y); jv=(int)floor(x);
				if (iu>=0 && iu<nl_TF-1 && jv>=0 && jv<nc_TF-1) {
					d1=(1-abs(y-iu))*(1-abs(x-jv)); d2=(1-abs(iu+1-y))*(1-abs(x-jv));
					d3=(1-abs(y-iu))*(1-abs(jv+1-x)); d4=(1-abs(iu+1-y))*(1-abs(jv+1-x));
					TF_1_logpol(i,j)=float((logTF_1(iu,jv)*d1+logTF_1(iu+1,jv)*d2+logTF_1(iu,jv+1)*d3+logTF_1(iu+1,jv+1)*d4)/(d1+d2+d3+d4));
					TF_2_logpol(i,j)=float((logTF_2(iu,jv)*d1+logTF_2(iu+1,jv)*d2+logTF_2(iu,jv+1)*d3+logTF_2(iu+1,jv+1)*d4)/(d1+d2+d3+d4));
				}
			}
		}
	}
	else {
		imabin imbTF(nl_TF,nc_TF); imbTF.mise_a_zero();
		for (i=0; i<nlnRho; i++) {
			rho=rhomin*exp(coef_lnRho*i);
			for (j=0; j<nTheta; j++) {
//				theta=(double)j*360./nTheta*PI/180.; 
				theta=(double)(j*pas_theta+thetamin)*PI/180.; 
				y=rho*cos(theta)+nl_TF/2.; x=rho*sin(theta)+nc_TF/2.;
				iu=around(rho*cos(theta))+nl_TF/2; jv=around(rho*sin(theta))+nc_TF/2;
				if (iu>=0 && iu<nl_TF && jv>=0 && jv<nc_TF) {TF_1_logpol(i,j)=logTF_1(iu,jv); TF_2_logpol(i,j)=logTF_2(iu,jv); imbTF(iu,jv)=1; }
			}
		}
		imbTF.imaunsignedchar().sauve_ImaPGM("./imb_logpol_TF.pgm");
	}
	float xfact, xshift, xx, xmax;
	cout<<" stat TF_1 en logpolaire :\n"; TF_1_logpol.statbasic(1); xshift=float(-TF_1_logpol.minI()); xfact=100.f/(float)(TF_1_logpol.maxI()-TF_1_logpol.minI()); cout<<" xshift = "<<xshift<<" xfact = "<<xfact<<"\n";
	TF_1_logpol=(TF_1_logpol+xshift)*xfact; TF_1_logpol.statbasic(1); TF_1_logpol.sauve_ImaPGM("./TF_1_logpol.pgm");
	cout<<" stat TF_2 en logpolaire :\n"; TF_2_logpol.statbasic(1); xshift=float(-TF_2_logpol.minI()); xfact=100.f/(float)(TF_2_logpol.maxI()-TF_2_logpol.minI()); cout<<" xshift = "<<xshift<<" xfact = "<<xfact<<"\n";
	TF_2_logpol=(TF_2_logpol+xshift)*xfact; TF_2_logpol.statbasic(1); TF_2_logpol.sauve_ImaPGM("./TF_2_logpol.pgm");
//	imadata<float> ima_corrphase_TF=TF_1_logpol.correlation_phase(TF_2_logpol,imax,jmax,xmax);
	imadata<float> ima_corrphase_TF=TF_1_logpol.correlation_phase(TF_2_logpol,imax,jmax,xmax,1);
	cout<<" par correlation de phase sur les TF des TF en coord. logpolaires, on trouve rotation d'angle "<<imax<<" et homothetie de parametre "<<jmax<<"\n";
	imadata<float> log_corrphase_TF(nlnRho,nTheta); log_corrphase_TF.mise_a_zero();
	for (i=0; i<nlnRho; i++) 
		for (j=0; j<nTheta; j++) if (ima_corrphase_TF(i,j,2)>0) log_corrphase_TF(i,j)=log(ima_corrphase_TF(i,j,2)); 
	log_corrphase_TF.sauve_ImaPGM("./log_corrphase_TF.pgm");
	const int npos=5;
	val_pos Tzmax[npos]; for (k=0; k<npos; k++) {Tzmax[k].val=-999; Tzmax[k].pos=0;}
	bool maxloc, trouve; float z; cout<<" nl_TF = "<<nl_TF<<", nc_TF = "<<nc_TF<<"\n";
	for (i=0; i<nlnRho; i++)
		for (j=0; j<nTheta; j++) {
			z=log_corrphase_TF(i,j); 
			maxloc=trouve=0; k=0;
			while (!maxloc && k<npos) {
				if (z>Tzmax[k].val) {
					for (k2=k-1; k2>=0; k2--) {
						xx=abs((int)(Tzmax[k2].pos%nTheta)-j)*pas_theta;
						if (!trouve && Tzmax[k2].pos/nTheta==i && (abs(xx-180)<pas_theta || abs(xx-360)<pas_theta || abs(xx-540)<pas_theta)) {
//							cout<<" trouve ! max en "<<k2<<" "<<j*pas_theta<<"deg., "<<i<<" val "<<z<<" coincide avec max "<<k2<<" en "<<(Tzmax[k2].pos%nTheta)*pas_theta<<"deg., "<<Tzmax[k2].pos/nTheta<<" val "<<Tzmax[k2].val<<"\n"; 
//							char aa; cin>>aa; 
							trouve=1;}
					}
					if (!trouve) {
						for (k2=npos-1; k2>k; k2--) {Tzmax[k2].val=Tzmax[k2-1].val; Tzmax[k2].pos=Tzmax[k2-1].pos; }
						Tzmax[k].val=z; Tzmax[k].pos=i*nTheta+j; maxloc=1;
					}
				}
				k++;
			}
		}
	z=(float)Tzmax[0].val; cout<<" max a "<<z<<" -> seuil de significativite des max a "<<z-log(2.f)<<"\n"; z-=log(2.f);
	for (k=0; k<npos; k++) if (Tzmax[k].val<z) Tzmax[k].val=-999.;
	float phi, alpha, angle=0, rapport=1, dx=0, dy=0, maxcorrph=0.;
	imadata<float> im_2_to_1, im_corr_phase;
	bool resize=0; unsigned char fill=0/*2*/;
	for (k=0; k<npos; k++) 
		if (Tzmax[k].val>-999.) {
			imax=-(int)(Tzmax[k].pos/nTheta); jmax=-(int)(Tzmax[k].pos%nTheta);
			if (imax<-nlnRho/2) imax+=nlnRho; if (jmax<-nTheta/2) jmax+=nTheta;
			phi=(float)(jmax*pas_theta+thetamin); alpha=exp(coef_lnRho*imax);
			cout<<"\n "<<k+1<<"eme max = "<<Tzmax[k].val<<" en ("<<Tzmax[k].pos/nTheta<<","<<Tzmax[k].pos%nTheta<<") -> phi = "<<phi<<"deg. & alpha = "<<alpha<<"\n";
			im_2_to_1=im_2.projette_ima(0.,0.,alpha,-phi,resize,fill); im_2.projette_ima(0.,0.,alpha,-phi,resize,fill).imaunsignedchar(1).sauve_ImaPGM("./im2_corr.pgm");
/*			double xxcorr=0., varI1=0., varI2=0., x1, x2;
			for (i=0; i<nl_I; i++)
				for (j=0; j<nc_I; j++) 
					for (k2=0; k2<nk_I; k2++)
						if (im_1(i,j,k2)>im_1.v_nulle() && im_2_to_1(i,j,k2)>im_1.v_nulle()) {x1=im_1(i,j,k2); x2=im_2_to_1(i,j,k2); varI1+=x1*x1; varI2+=x2*x2; xxcorr+=x1*x2;}
			cout<<" apres projection : correlation = "<<xxcorr/pow(varI1*varI2,0.5)<<"\n";*/
			im_corr_phase=im_1.correlation_phase(im_2_to_1,imax,jmax,xmax);
			if (xmax>maxcorrph) {maxcorrph=xmax; dx=float(-imax); dy=float(-jmax); rapport=alpha; angle=-phi; cout<<" max de correlation de phase en "<<dx<<" "<<dy<<" de valeur "<<maxcorrph<<"\n";}
		}
	cout<<"\n\n finalement : rapport="<<rapport<<", rotation = "<<angle<<" translation = ("<<dx/cos(angle*PI/180.)*rapport<<","<<dy/sin(angle*PI/180.)*rapport<<")\n\n";
	ima2.projette_ima(0,0,rapport,angle).imaunsignedchar(1).sauve_ImaPGM("./ima2_corr_phi_alpha.pgm");
	ima2.projette_ima(0,0,rapport,angle).projette_ima(dx,dy,1,0).imaunsignedchar(1).sauve_ImaPGM("./ima2_corr_phi_alpha_dx_dy.pgm");
//	ima2.projette_ima(dx/rapport/cos(angle),dy/rapport/sin(angle),rapport,angle).imaunsignedchar(1).sauve_ImaPGM("./ima2_corr.pgm");
//	ima2.projette_ima(dx/rapport/sin(angle),dy/rapport/cos(angle),rapport,angle).imaunsignedchar(1).sauve_ImaPGM("./ima2_corr.pgm");
	float cosangle=(float)cos(-angle*PI/180.), sinangle=(float)sin(-angle*PI/180.), dx2=(dx*cosangle-dy*sinangle)/rapport, dy2=(dx*sinangle+dy*cosangle)/rapport;
	cout<<" transformation trouvee : translation = ("<<dx2<<","<<dy2<<"), angle = "<<angle<<", rapport = "<<rapport<<"\n";
	ima2.projette_ima(dx2,dy2,rapport,angle).imaunsignedchar(1).sauve_ImaPGM("./ima2_corr.pgm");
	cout<<angle<<" "<<cosangle<<" "<<sinangle<<" "<<dx2<<" "<<dy2<<"\n";
	_dx=dx2; _dy=dy2; _rapport=rapport; _angle=angle;
	return log_corrphase_TF;
}

template <class T> imadata<T> imadata<T>::projette_ima_geom_fixe(const double dx, const double dy, const double hfact, const double phi, const int nblig2, const int nbcol2, unsigned char fill) const { //cout<<" projette_ima_0 ";
	cout<<" dx = "<<dx<<", dy = "<<dy<<", hfact = "<<hfact<<", phi = "<<phi<<"\n";
	const int i0=nblig/2, j0=nbcol/2;
	const double cosphi=cos(phi/180.*PI), sinphi=sin(phi/180.*PI);
	matrice2D<double> A(3,2);
	A(0,1)=-dx+i0-nblig2/2*cosphi/hfact-nbcol2/2*sinphi/hfact; A(0,0)=-dy+j0+nblig2/2*sinphi/hfact-nbcol2/2*cosphi/hfact; 
	A(2,1)=cosphi/hfact; A(2,0)=-sinphi/hfact; A(1,1)=sinphi/hfact; A(1,0)=cosphi/hfact; 
//	cout<<nblig<<" "<<nblig2<<" "<<nbcol<<" "<<nbcol2<<" "<<i0<<" "<<j0<<" A = \n"; A.affiche();
	imabin imab;
	imadata<T> imaproj=projette_ima(A,nblig2,nbcol2,imab,1,fill); imaproj.statbasic(1); 
	return imaproj;
}

template <class T> imadata<T> imadata<T>::projette_ima(const double dx, const double dy, const double hfact, const double phi, bool resize, unsigned char fill) const { //cout<<" projette_ima_1 ";
	cout<<" dx = "<<dx<<", dy = "<<dy<<", hfact = "<<hfact<<", phi = "<<phi<<"\n";
	const int i0=nblig/2, j0=nbcol/2;
	const double cosphi=cos(phi/180.*PI), sinphi=sin(phi/180.*PI);
	double x,y,/*x0,y0,*/xmin=DBL_MAX,xmax=DBL_MIN,ymin=DBL_MAX,ymax=DBL_MIN;
	int i,j/*,k*/;
	int nblig2=nblig,nbcol2=nbcol/*,i2,j2*/;
	if (resize) {
		for (i=0; i<nblig; i+=nblig-1)
			for (j=0; j<nbcol; j+=nbcol-1) {
				x=dx+((i-i0)*cosphi-(j-j0)*sinphi)*hfact; y=dy+((i-i0)*sinphi+(j-j0)*cosphi)*hfact;
				if (x<xmin) xmin=x; if (x>xmax) xmax=x; if (y<ymin) ymin=y; if (y>ymax) ymax=y;
			}
		nblig2=(int)(xmax-xmin+2); nbcol2=(int)(ymax-ymin+2);
	}
	matrice2D<double> A(3,2);
	A(0,1)=-dx+i0-nblig2/2*cosphi/hfact-nbcol2/2*sinphi/hfact; A(0,0)=-dy+j0+nblig2/2*sinphi/hfact-nbcol2/2*cosphi/hfact; 
//	A(0,1)=-dx+i0-nblig2/2*cosphi/hfact-nbcol2/2*sinphi/hfact; A(0,0)=-dy+j0+nblig2/2*sinphi/hfact-nbcol2/2*cosphi/hfact; 
	A(2,1)=cosphi/hfact; A(2,0)=-sinphi/hfact; A(1,1)=sinphi/hfact; A(1,0)=cosphi/hfact; 
	imabin imab;
	imadata<T> imaproj=projette_ima(A,nblig2,nbcol2,imab,1,fill); imaproj.statbasic(1);
	return imaproj;
}

template <class T> imadata<T> imadata<T>::projette_ima(matrice2D<double> &A, int nblig2, int nbcol2, imabin &imab, unsigned char i_interpol, unsigned char fill) const { //cout<<" projette_ima_2 : fill = "<<(int)fill<<"\n";
	double x,y,xmin=DBL_MAX,xmax=DBL_MIN,ymin=DBL_MAX,ymax=DBL_MIN,d1,d2,d3,d4,z;
	int ii,jj,i,j,k,d=A.nlig();
	imadata<T> imaproj(nblig2,nbcol2,nbcanaux);
	imab=imabin(nblig2,nbcol2); imab.mise_a_un();
	for (i=0; i<nblig2; i++) // boucle sur les lignes de l'image sur laquelle on recale
		for (j=0; j<nbcol2; j++) { // boucle sur les colonnes de l'image sur laquelle on recale
			x=A(0,0)+A(1,0)*j+A(2,0)*i; if (d>5) x+=A(3,0)*pow((float)j,2)+A(4,0)*pow((float)i,2)+A(5,0)*j*i; // estimation de la coordonnee colonne dans l'image ?recaler
			y=A(0,1)+A(1,1)*j+A(2,1)*i; if (d>5) y+=A(3,1)*pow((float)j,2)+A(4,1)*pow((float)i,2)+A(5,1)*j*i; // estimation de la coordonnee ligne dans l'image ?recaler
			jj=around(x); ii=around(y);
			if (i_interpol==0) { // interpolation au plus proche voisin
				jj=around(x); ii=around(y);
				if (ii>=0 && ii<nblig && jj>=0 && jj<nbcol)
					for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=(*this)(ii,jj,k);
				else {
					imab(i,j)=0;
					if (fill>0) { // 3 façons de remplir (hormis 0) : (1):duplication des bords, (2): periodisation de l'image, (3): symetrie de l'image 
						switch ((int)fill) {
							case 1: cout<<" cas 1 pour remplissage\n"; ii=mini(maxi(ii,0),nblig-1); jj=mini(maxi(jj,0),nbcol-1); break;
							case 2: cout<<" cas 2 pour remplissage\n"; if (ii>=nblig) ii=ii%nblig; else if (ii<0) ii=(ii+(1-ii/nblig)*nblig)%nblig; 
											if (jj>=nbcol) jj=jj%nbcol; else if (jj<0) jj=(jj+(1-jj/nbcol)*nbcol)%nbcol; break;
							case 3: cout<<" cas 3 pour remplissage\n"; if (ii>=nblig) ii=nblig-1-ii%nblig; else if (ii<0) ii=nblig-1-(ii+(1-ii/nblig)*nblig)%nblig; 
											if (jj>=nbcol) jj=nbcol-1-jj%nbcol; else if (jj<0) jj=nbcol-1-(jj+(1-jj/nbcol)*nbcol)%nbcol; break;
							default: ii=mini(maxi(ii,0),nblig-1); jj=mini(maxi(jj,0),nbcol-1); 
						}
						for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=(*this)(ii,jj,k);
					}
					else for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=0;
				}
			} 
			else {						 // interpolation a l'ordre 1
				ii=(int)floor(y); jj=(int)floor(x);                                  
				if (ii>=0 && ii<nblig-1 && jj>=0 && jj<nbcol-1) {
					d1=(1-abs(y-ii))*(1-abs(x-jj)); d2=(1-abs(ii+1-y))*(1-abs(x-jj));
					d3=(1-abs(y-ii))*(1-abs(jj+1-x)); d4=(1-abs(ii+1-y))*(1-abs(jj+1-x));
//				d1=1-maxi(abs(y-ii),abs(x-jj)); d2=1-maxi(abs(ii+1-y),abs(x-jj));
//				d3=1-maxi(abs(y-ii),abs(jj+1-x)); d4=1-maxi(abs(ii+1-y),abs(jj+1-x));
					for (k=0; k<nbcanaux; k++) {
						z=((*this)(ii,jj,k)*d1+(*this)(ii+1,jj,k)*d2+(*this)(ii,jj+1,k)*d3+(*this)(ii+1,jj+1,k)*d4)/(d1+d2+d3+d4);
						imaproj(i,j,k)=(T)z;
					}
				}
				else {
					if (ii>=0 && ii<nblig && jj>=0 && jj<nbcol) for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=(*this)(ii,jj,k);
					else {
						imab(i,j)=0;
						if (fill>0) { // 3 façons de remplir : (1):duplication des bords, (2): periodisation de l'image, (3): symetrie de l'image 
							switch ((int)fill) {
								case 1: ii=mini(maxi(ii,0),nblig-1); jj=mini(maxi(jj,0),nbcol-1); break;
								case 2: if (ii>=nblig) ii=ii%nblig; else if (ii<0) ii=(ii+(1-ii/nblig)*nblig)%nblig; 
												if (jj>=nbcol) jj=jj%nbcol; else if (jj<0) jj=(jj+(1-jj/nbcol)*nbcol)%nbcol; break;
								case 3: if (ii>=nblig) ii=nblig-1-ii%nblig; else if (ii<0) ii=nblig-1-(ii+(1-ii/nblig)*nblig)%nblig; 
												if (jj>=nbcol) jj=nbcol-1-jj%nbcol; else if (jj<0) jj=nbcol-1-(jj+(1-jj/nbcol)*nbcol)%nbcol; break;
								default: ii=mini(maxi(ii,0),nblig-1); jj=mini(maxi(jj,0),nbcol-1); 
							}
							for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=(*this)(ii,jj,k);
						}
						else {for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=0;}
					}
				}
			}
		}
	imaproj.statbasic(1);
	return imaproj;
}

template <class T> imadata<T> imadata<T>::projette_ima(matrice2D<double> &A, unsigned char i_interpol) const {
	double x,y,xmin=DBL_MAX,xmax=DBL_MIN,ymin=DBL_MAX,ymax=DBL_MIN,d1,d2,d3,d4,z;
	int ii,jj,i,j,k,d=A.nlig(),nblig2=nblig,nbcol2=nbcol;
	imadata<T> imaproj(nblig2,nbcol2,nbcanaux);
	for (i=0; i<nblig2; i++) // boucle sur les lignes de l'image sur laquelle on recale
		for (j=0; j<nbcol2; j++) { // boucle sur les colonnes de l'image sur laquelle on recale
			x=A(0,0)+A(1,0)*j+A(2,0)*i; if (d>5) x+=A(3,0)*pow((float)j,2)+A(4,0)*pow((float)i,2)+A(5,0)*j*i; // estimation de la coordonnee colonne dans l'image ?recaler
			y=A(0,1)+A(1,1)*j+A(2,1)*i; if (d>5) y+=A(3,1)*pow((float)j,2)+A(4,1)*pow((float)i,2)+A(5,1)*j*i; // estimation de la coordonnee ligne dans l'image ?recaler
			if (i_interpol==0) { // interpolation au plus proche voisin
				jj=around(x); ii=around(y);
				if (ii>=0 && ii<nblig && jj>=0 && jj<nbcol) 
					for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=(*this)(ii,jj);
				else 
					for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=0;
			} else {						 // interpolation a l'ordre 1
				ii=(int)floor(y); jj=(int)floor(x);                                  
				if (ii>=0 && ii<nblig-1 && jj>=0 && jj<nbcol-1) {
					d1=(1-abs(y-ii))*(1-abs(x-jj)); d2=(1-abs(ii+1-y))*(1-abs(x-jj));
					d3=(1-abs(y-ii))*(1-abs(jj+1-x)); d4=(1-abs(ii+1-y))*(1-abs(jj+1-x));
//				d1=1-maxi(abs(y-ii),abs(x-jj)); d2=1-maxi(abs(ii+1-y),abs(x-jj));
//				d3=1-maxi(abs(y-ii),abs(jj+1-x)); d4=1-maxi(abs(ii+1-y),abs(jj+1-x));
					for (k=0; k<nbcanaux; k++) {
						z=((*this)(ii,jj,k)*d1+(*this)(ii+1,jj,k)*d2+(*this)(ii,jj+1,k)*d3+(*this)(ii+1,jj+1,k)*d4)/(d1+d2+d3+d4);
						imaproj(i,j,k)=(T)z;
					}
				} 
				else {
					if (ii>=0 && ii<nblig && jj>=0 && jj<nbcol) for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=(*this)(ii,jj);
					else for (k=0; k<nbcanaux; k++) imaproj(i,j,k)=0;
				}
			}
		}
	return imaproj;
}

/*template <class T> imadata<T> imadata<T>::dilate(const eltstruct B) {
	int i,j,k,i0,i2,j0,j2,ii,jj,iB;
	int nl=B.nl, nc=B.nc, x0=B.x0, y0=B.y0;
	T z, *z0=new T[nbcanaux];
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	if (!istat) statbasic();
	for (i=0; i<nblig; i++) {
		i0=maxi(i-y0,0); i2=mini(i+nl-1-y0,nblig-1);
		for (j=0; j<nbcol; j++) {
			j0=maxi(j-x0,0); j2=mini(j+nc-1-x0,nbcol-1);
			for (k=0; k<nbcanaux; k++) z0[k]=mini(valnul,valmin[k]);
			for (ii=i0; ii<=i2; ii++) {
				iB=(ii-i+y0)*nc;
				for (jj=j0; jj<=j2; jj++)
					if (B.aT[iB+jj-j+x0]) {
						for (k=0; k<nbcanaux; k++) {
							z=(*this)(ii,jj,k);
							if (z>z0[k]) z0[k]=z;
						}
				}
			}
			for (k=0; k<nbcanaux; k++) imaRes(i,j,k)=z0[k];
		}
	}
	if (z0!=NULL) {delete[] z0; z0=NULL;}
	return imaRes;
}*/

template <class T> imadata<T> imadata<T>::dilate(const eltstruct B, int n) {
	if (n<=0) return *this;
	else {
		int i,j,k,i0,i2,j0,j2,ii,jj,iB,nn;
		int nl=B.nl, nc=B.nc, x0=B.x0, y0=B.y0;
		T z, *z0=new T[nbcanaux];
		imadata<T> imaRes(nblig,nbcol,nbcanaux), imaC(nblig,nbcol,nbcanaux);
		if (!istat) statbasic();
		imaC=(*this);
		for (nn=0; nn<n; nn++) {
			for (i=0; i<nblig; i++) {
				i0=maxi(i-y0,0); i2=mini(i+nl-1-y0,nblig-1);
				for (j=0; j<nbcol; j++) {
					j0=maxi(j-x0,0); j2=mini(j+nc-1-x0,nbcol-1);
					for (k=0; k<nbcanaux; k++) z0[k]=mini(valnul,valmin[k]);
					for (ii=i0; ii<=i2; ii++) {
						iB=(ii-i+y0)*nc;
						for (jj=j0; jj<=j2; jj++)
							if (B.aT[iB+jj-j+x0]) {
								for (k=0; k<nbcanaux; k++) {
									z=imaC(ii,jj,k);
									if (z>z0[k]) z0[k]=z;
                }
              }
          }
          for (k=0; k<nbcanaux; k++) imaRes(i,j,k)=z0[k];
        }
			}
			imaC=imaRes;
		}
		if (z0!=NULL) {delete[] z0; z0=NULL;}
		return imaRes;
  }
}

/*template <class T> imadata<T> imadata<T>::erode(const eltstruct B) {
	int i,j,k,i0,i2,j0,j2,ii,jj,iB;
	int nl=B.nl, nc=B.nc, x0=B.x0, y0=B.y0;
	T z, *z0=new T[nbcanaux];
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	if (!istat) statbasic();
	for (i=0; i<nblig; i++) {
		i0=maxi(i-y0,0); i2=mini(i+nl-1-y0,nblig-1);
		for (j=0; j<nbcol; j++) {
			j0=maxi(j-x0,0); j2=mini(j+nc-1-x0,nbcol-1);
			for (k=0; k<nbcanaux; k++) z0[k]=maxi(valnul,valmax[k]);
			for (ii=i0; ii<=i2; ii++) {
				iB=(ii-i+y0)*nc;
				for (jj=j0; jj<=j2; jj++)
					if (B.aT[iB+jj-j+x0]) {
						for (k=0; k<nbcanaux; k++) {
							z=(*this)(ii,jj,k);
							if (z<z0[k]) z0[k]=z;
						}
				}
			}
			for (k=0; k<nbcanaux; k++) imaRes(i,j,k)=z0[k];
		}
	}
	if (z0!=NULL) {delete[] z0; z0=NULL;}
	return imaRes;
}*/

template <class T> imadata<T> imadata<T>::erode(const eltstruct B, int n) {
	if (n<=0) return *this;
	else {
		int i,j,k,i0,i2,j0,j2,ii,jj,iB,nn;
		int nl=B.nl,nc=B.nc,x0=B.x0,y0=B.y0;
		T z, *z0=new T[nbcanaux];
		imadata<T> imaRes(nblig,nbcol,nbcanaux), imaC(nblig,nbcol,nbcanaux);
		if (!istat) statbasic();
		imaC=(*this);
		for (nn=0; nn<n; nn++) {
			for (i=0; i<nblig; i++) {
				i0=maxi(i-y0,0); i2=mini(i+nl-1-y0,nblig-1);
				for (j=0; j<nbcol; j++) {
					j0=maxi(j-x0,0); j2=mini(j+nc-1-x0,nbcol-1);
					for (k=0; k<nbcanaux; k++) z0[k]=maxi(valnul,valmax[k]);
					for (ii=i0; ii<=i2; ii++) {
						iB=(ii-i+y0)*nc;
						for (jj=j0; jj<=j2; jj++)
							if (B.aT[iB+jj-j+x0]) {
								for (k=0; k<nbcanaux; k++) {
									z=imaC(ii,jj,k);
									if (z<z0[k]) z0[k]=z;
								}
							}
					}
					for (k=0; k<nbcanaux; k++) imaRes(i,j,k)=z0[k];
				}
			}
			imaC=imaRes;
		}
		if (z0!=NULL) {delete[] z0; z0=NULL;}
		return imaRes;
	}
}

template <class T> imadata<T> imadata<T>::rehausse_contraste(const eltstruct B, int n, float alpha, float beta) {
	imadata<T> imaRes(nblig,nbcol,nbcanaux), imaE(nblig,nbcol,nbcanaux),
	           imaD(nblig,nbcol,nbcanaux);
	T x, xe, xd, delta;
	int i,j,k;
	imaE=erode(B,n);
	imaD=dilate(B,n);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				x=(*this)(i,j,k);
				xe=imaE(i,j,k);
				xd=imaD(i,j,k);
				delta=xd-xe;
				if (x<xe+alpha*delta) imaRes(i,j,k)=xe;
				else {
					if (x<=xd-beta*delta) imaRes(i,j,k)=x;
					else imaRes(i,j,k)=xd;
				}
			}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::gradient_m (const eltstruct B) {
	imadata<T> imaRes(nblig,nbcol,nbcanaux), imaE=erode(B), imaD=dilate(B);
	int i,j,k;
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				imaRes(i,j,k)=imaD(i,j,k)-imaE(i,j,k);
	return imaRes;
	}

template <class T> imadata<T> imadata<T>::laplacien_m (const eltstruct B) {
	imadata<T> imaRes(nblig,nbcol,nbcanaux), imaE=erode(B), imaD=dilate(B);
	int i,j,k;
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				imaRes(i,j,k)=imaD(i,j,k)+imaE(i,j,k)-2*(*this)(i,j,k);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::ouverture(const eltstruct B) {
	int i,j,k;
	int nl1=B.nl-1,dl0=nl1-B.y0,dl2=B.y0,nc1=B.nc-1,dc0=nc1-B.x0,dc2=B.x0;
	imadata<T> imaC(nblig+nl1,nbcol+nc1,nbcanaux);
	T x;
	for (k=0; k<nbcanaux; k++) {
		for (i=0; i<dl0; i++) {
			x=(*this)(0,0,k); for (j=0; j<dc0; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(0,j,k);
		}
		for (j=nbcol+dc0; j<nbcol+nc1; j++) {
			x=(*this)(0,nbcol-1,k); for (i=0; i<dl0; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,nbcol-1,k);
		}
		for (i=nblig+dl0; i<nblig+nl1; i++) {
			x=(*this)(nblig-1,nbcol-1,k); for (j=nbcol+dc0; j<nbcol+nc1; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(nblig-1,j,k);
		}
		for (j=0; j<dc0; j++) {
			x=(*this)(nblig-1,0,k); for (i=nblig+dl0; i<nblig+nl1; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,0,k);
		}
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaC(i+dl0,j+dc0,k)=(*this)(i,j,k);
	}
	eltstruct Bt(B); Bt=Bt.transpose();
	imadata<T> imaR=imaC.erode(B).dilate(Bt);
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaRes(i,j,k)=imaR(i+dl0,j+dc0,k);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::ouverture(const eltstruct B, int n) {
	int i,j,k;
	int nl1=(B.nl-1)*n,dl0=(B.nl-1-B.y0)*n,dl2=(B.y0)*n;
	int nc1=(B.nc-1)*n,dc0=(B.nc-1-B.x0)*n,dc2=(B.x0)*n;
	imadata<T> imaC(nblig+nl1,nbcol+nc1,nbcanaux);
	T x;
	for (k=0; k<nbcanaux; k++) {
		for (i=0; i<dl0; i++) {
			x=(*this)(0,0,k); for (j=0; j<dc0; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(0,j,k);
		}
		for (j=nbcol+dc0; j<nbcol+nc1; j++) {
			x=(*this)(0,nbcol-1,k); for (i=0; i<dl0; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,nbcol-1,k);
		}
		for (i=nblig+dl0; i<nblig+nl1; i++) {
			x=(*this)(nblig-1,nbcol-1,k); for (j=nbcol+dc0; j<nbcol+nc1; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(nblig-1,j,k);
		}
		for (j=0; j<dc0; j++) {
			x=(*this)(nblig-1,0,k); for (i=nblig+dl0; i<nblig+nl1; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,0,k);
		}
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaC(i+dl0,j+dc0,k)=(*this)(i,j,k);
	}
	eltstruct Bt(B); Bt=Bt.transpose();
	imadata<T> imaR=imaC.erode(B,n).dilate(Bt,n);
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaRes(i,j,k)=imaR(i+dl0,j+dc0,k);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::fermeture(const eltstruct B) {
	int i,j,k;
	int nl1=B.nl-1,dl0=nl1-B.y0,dl2=B.y0,nc1=B.nc-1,dc0=nc1-B.x0,dc2=B.x0;
	imadata<T> imaC(nblig+nl1,nbcol+nc1,nbcanaux);
	T x;
	for (k=0; k<nbcanaux; k++) {
		for (i=0; i<dl0; i++) {
			x=(*this)(0,0,k); for (j=0; j<dc0; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(0,j,k);
		}
		for (j=nbcol+dc0; j<nbcol+nc1; j++) {
			x=(*this)(0,nbcol-1,k); for (i=0; i<dl0; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,nbcol-1,k);
		}
		for (i=nblig+dl0; i<nblig+nl1; i++) {
			x=(*this)(nblig-1,nbcol-1,k); for (j=nbcol+dc0; j<nbcol+nc1; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(nblig-1,j,k);
		}
		for (j=0; j<dc0; j++) {
			x=(*this)(nblig-1,0,k); for (i=nblig+dl0; i<nblig+nl1; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,0,k);
		}
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaC(i+dl0,j+dc0,k)=(*this)(i,j,k);
	}
	eltstruct Bt(B); Bt=Bt.transpose();
	imadata<T> imaR=imaC.dilate(B).erode(Bt), imaRes(nblig,nbcol,nbcanaux);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaRes(i,j,k)=imaR(i+dl0,j+dc0,k);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::fermeture(const eltstruct B, int n) {
	int i,j,k;
	int nl1=(B.nl-1)*n,dl0=(B.nl-1-B.y0)*n,dl2=(B.y0)*n;
	int nc1=(B.nc-1)*n,dc0=(B.nc-1-B.x0)*n,dc2=(B.x0)*n;
	imadata<T> imaC(nblig+nl1,nbcol+nc1,nbcanaux);
	T x;
	for (k=0; k<nbcanaux; k++) {
		for (i=0; i<dl0; i++) {
			x=(*this)(0,0,k); for (j=0; j<dc0; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(0,j,k);
		}
		for (j=nbcol+dc0; j<nbcol+nc1; j++) {
			x=(*this)(0,nbcol-1,k); for (i=0; i<dl0; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,nbcol-1,k);
		}
		for (i=nblig+dl0; i<nblig+nl1; i++) {
			x=(*this)(nblig-1,nbcol-1,k); for (j=nbcol+dc0; j<nbcol+nc1; j++) imaC(i,j,k)=x;
			for (j=0; j<nbcol; j++) imaC(i,j+dc0,k)=(*this)(nblig-1,j,k);
		}
		for (j=0; j<dc0; j++) {
			x=(*this)(nblig-1,0,k); for (i=nblig+dl0; i<nblig+nl1; i++) imaC(i,j,k)=x;
			for (i=0; i<nblig; i++) imaC(i+dl0,j,k)=(*this)(i,0,k);
		}
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaC(i+dl0,j+dc0,k)=(*this)(i,j,k);
	}
	eltstruct Bt(B); Bt=Bt.transpose();
	imadata<T> imaR=imaC.dilate(B,n).erode(Bt,n);
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) imaRes(i,j,k)=imaR(i+dl0,j+dc0,k);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::tophat(const eltstruct B) {
	imadata<T> imaRes(nblig,nbcol,nbcanaux), imaC=ouverture(B);
	imaRes=(*this)-imaC;
	return imaRes;
}

template <class T> imadata<T> imadata<T>::tophat_c(const eltstruct B) {
	imadata<T> imaRes(nblig,nbcol,nbcanaux);
	imaRes=fermeture(B)-(*this);
	return imaRes;
}

template <class T> imadata<T> imadata<T>::filtrealterne(const eltstruct B, int n, bool sens) {
	imadata<T> imaRes(nblig,nbcol,nbcanaux), imaC;
	imaRes=(*this);
	if (n>0) {
		int i;
		if (!sens) {
			for (i=1; i<=n; i++) {
				imaC=imaRes;
				imaRes=imaC.ouverture(B,i);
				imaC=imaRes;
				imaRes=imaC.fermeture(B,i);
			}
		}
		else {
			for (i=1; i<=n; i++) {
				imaC=imaRes;
				imaRes=imaC.fermeture(B,i);
				imaC=imaRes;
				imaRes=imaC.ouverture(B,i);
			}
		}
	}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::reconst_geod(const imadata<T> &imaMq, const eltstruct B) {
	imadata<T> imaRes(nblig,nbcol,nbcanaux), imaI;
	imaRes=imaMq;
	int n=1, i, j;
	while (n>0) {
		imaI=imaRes.dilate(B);
		n=0;
		for (i=0; i<imaI.nblig; i++)
			for (j=0; j<imaI.nbcol; j++) {
				if (imaI(i,j)>(*this)(i,j)) imaI(i,j)=(*this)(i,j);
				if (imaI(i,j)!=imaRes(i,j)) n++;
			}
		if (n>0) imaRes=imaI;
	}
	return imaRes;
}

template <class T> imadata<T> imadata<T>::marq_centres_lowerset(const eltstruct B, int nm) {
	const int nmax=(nm>0?nm:100);
	T xmin, xmax, dx, x, x0;
	int i,j,k,n,ncc,icc,ncc_iccO,jcc;
	imabin imab_LS(nblig,nbcol), imab_L0(nblig,nbcol), imab, imabO;
	imadata<int> imacc, imacc_iccO;
	imadata<float> imadist;
	imadata<T> imares(*this);
	bool iO,iD;
	for (k=0; k<nbcanaux; k++) {
		xmin=valmin[k], xmax=valmax[k];
		dx=(xmax-xmin)/nmax;
		for (n=0; n<nmax; n++) {
			x=xmin+dx*(n+1); x0=x-dx;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) imab_LS(i,j)=((*this)(i,j,k)<x?1:0);
			imacc=imab_LS.composantes_connexes(ncc);
			for (icc=1; icc<=ncc; icc++) {
				iO=0;
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						imab(i,j)=(imacc(i,j)==icc?1:0);
						if (!iO && imab_L0(i,j)) iO=1;
					}
				if (iO) {
					imabO=imab.ouverture(B);
					imacc_iccO=imab.composantes_connexes(ncc_iccO);
					if (ncc_iccO>1) {
						for (jcc=1; jcc<=ncc_iccO; jcc++) {
							iD=1;
							for (i=0; i<nblig; i++)
								for (j=0; j<nbcol; j++)
									if (imacc_iccO(i,j)==jcc && iD && imab_L0(i,j)) iD=0;
							if (iD) {
								for (i=0; i<nblig; i++)
									for (j=0; j<nbcol; j++)
										imab(i,j)=(imacc_iccO(i,j)==jcc?0:1);
								imadist=imab.Tr_dist();
								dmax=imadist.maxI();
								for (i=0; i<nblig; i++)
									for (j=0; j<nbcol; j++)
										if (imadist(i,j)==dmax) imares(i,j)-=dx;
							}
						}
					}
				}
			}
		}
	}
	return imares;
}

template <class T> imadata<T> imadata<T>::h_max(const eltstruct B, float h) {
	imadata<float> imares(*this);
//	h=-h;
	imares=imares.reconst_geod(imares+(-h),B);
	return imares;
}

template <class T> imadata<T> imadata<T>::marq_max_regionaux(const eltstruct B, float h) {
	imadata<float> imares(*this);
	imabin ima_mask(*this,h);
	imares=imares.reconst_geod(imares+(-h),B);
	imares=imares-imares.reconst_geod(imares+(-1),B);
	imares=imares*imadata<float>(ima_mask);
	return imares;
}
/*
template <class T> imadata<T> imadata<T>::marq_max_regionaux(const eltstruct B, float h) {
	imadata<float> imares(*this);
	h=-h;
	imares=imares.reconst_geod(imares+h,B);
	imares=imares-imares.reconst_geod(imares+(-1),B);
	return imares;
}
*/

template <class T> imadata<BYTE> imadata<T>::seuil_ima (float s) const {
	imadata<BYTE> imaR(nblig,nbcol);
	int i=0,j=0;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j)>s) imaR(i,j)=255;
			else imaR(i,j)=0;
	return imaR;
}

template <class T> void imadata<T>::seuil_ima (float s, imadata<BYTE> &imaR, int k) const {
	int i=0,j=0;
	cout<<" seuillage du canal "<<k<<" de l'image avec seuil = "<<s<<"\n";
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			if ((*this)(i,j,k)>s) imaR(i,j,k)=1;
			else imaR(i,j,k)=0;
		}
}

template <class T> void imadata<T>::seuil_hysteresis (float sh, float sb, imabin &imab, int k) {
/*	int i=0,j=0;
	cout<<" seuillage par hysteresis de l'image avec seuil haut = "<<sh<<" et seuil bas = "<<sb<<"\n";
	liste_pixels Lpix;
	elt_liste E;
	imabin ima_deja(nblig,nbcol), ima_lpix(nblig,nbcol);
	imab=ima_deja;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j,k)>=sh) {Lpix.insere(i,j); ima_lpix(i,j)=1; }
	while (Lpix.nb_elts()>0) {
		E=Lpix.extrait(); ima_lpix(E.x,E.y)=0;
		imab(E.x,E.y)=1;
		ima_deja(E.x,E.y)=1;
		for (i=maxi(E.x-1,0); i<=mini(E.x+1,nblig-1); i++)
			for (j=maxi(E.y-1,0); j<=mini(E.y+1,nbcol-1); j++)
				if (!ima_deja(i,j) && !ima_lpix(i,j) && (*this)(i,j,k)>=sb) {
					Lpix.insere(i,j); ima_lpix(i,j)=1;}
				else ima_deja(E.x,E.y)=1;
	}*/
	cout<<" seuillage par hysteresis (MM) de l'image avec seuil haut = "<<sh<<" et seuil bas = "<<sb<<"\n";
	imabin ima_sh(*this,sh,k), ima_sb(*this,sb,k);
	cout<<" images binaires superieur seuil haut et superieur seuil bas calculees\n";
	int i,j,ncc,n,iconnex=8;
  imadata<int> imacc=ima_sb.composantes_connexes(ncc,iconnex);
	cout<<" nombre de composantes connexes pour seuil bas = "<<ncc<<"\n";
	bool iccOK;
	if (imab.nlig()<nblig || imab.ncol()<nbcol) imab=imabin(nblig,nbcol); imab.mise_a_zero();
	for (n=1; n<=ncc; n++) {
		iccOK=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (!iccOK && imacc(i,j)==n && ima_sh(i,j)) iccOK=1;
		if (iccOK) {
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) if (imacc(i,j)==n) imab(i,j)=1;;
		}
	}
}

template <class T> imadata<BYTE> imadata<T>::seuil_hysteresis (float sh, float sb, int k) const {
	int i=0,j=0;
	cout<<" seuillage par hysteresis de l'image avec seuil haut = "<<sh<<" et seuil bas = "<<sb<<"\n";
	liste_pixels Lpix;
	elt_liste E;
	imabin ima_deja(nblig,nbcol), ima_lpix(nblig,nbcol);
  imadata<BYTE> imaR(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j,k)>=sh) {Lpix.insere(i,j); ima_lpix(i,j)=1; }
	while (Lpix.nb_elts()>0) {
		E=Lpix.extrait(); ima_lpix(E.x,E.y)=0;
		imaR(E.x,E.y)=255;
		ima_deja(E.x,E.y)=1;
		for (i=maxi(E.x-1,0); i<=mini(E.x+1,nblig-1); i++)
			for (j=maxi(E.y-1,0); j<=mini(E.y+1,nbcol-1); j++)
//				if (!ima_deja(i,j) && !Lpix.existe(i,j) && (*this)(i,j,k)>sb) Lpix.insere(i,j);
				if (!ima_deja(i,j) && !ima_lpix(i,j) && (*this)(i,j,k)>=sb) {
					Lpix.insere(i,j); ima_lpix(i,j)=1;}
				else ima_deja(E.x,E.y)=1;
	}
  return imaR;
}

template <class T> void imadata<T>::seuil_hysteresis (float sh, float sb, imadata<BYTE> &imaR, int k) {
/*	int i=0,j=0;
	cout<<" seuillage par hysteresis de l'image avec seuil haut = "<<sh<<" & seuil bas = "<<sb<<"\n";
	liste_pixels Lpix;
	elt_liste E;
	imabin ima_deja(nblig,nbcol), ima_lpix(nblig,nbcol);
	if (imaR.nlig()<nblig || imaR.ncol()<nbcol) {imaR=imadata<BYTE>(nblig,nbcol,nbcanaux); imaR.mise_a_zero();}
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j,k)>=sh) {Lpix.insere(i,j); ima_lpix(i,j)=1; }
	while (Lpix.nb_elts()>0) {
		E=Lpix.extrait(); ima_lpix(E.x,E.y)=0;
		imaR(E.x,E.y,k)=1;
		ima_deja(E.x,E.y)=1;
		for (i=maxi(E.x-1,0); i<=mini(E.x+1,nblig-1); i++)
			for (j=maxi(E.y-1,0); j<=mini(E.y+1,nbcol-1); j++)
//				if (!ima_deja(i,j) && !Lpix.existe(i,j) && (*this)(i,j,k)>sb) Lpix.insere(i,j);
				if (!ima_deja(i,j) && !ima_lpix(i,j) && (*this)(i,j,k)>=sb) {
					Lpix.insere(i,j); ima_lpix(i,j)=1;}
				else ima_deja(E.x,E.y)=1;
	}*/
	cout<<" seuillage par hysteresis (MM) de l'image avec seuil haut = "<<sh<<" et seuil bas = "<<sb<<"\n";
	imabin ima_sh(*this,sh,k), ima_sb(*this,sb,k);
	cout<<" images binaires superieur seuil haut et superieur seuil bas calculees : respectivement normes = "<<ima_sh.norm()<<" et "<<ima_sb.norm()<<"\n";
	int i,j,ncc,n,iconnex=8;
  imadata<int> imacc=ima_sb.composantes_connexes (ncc,iconnex);
	cout<<" nombre de composantes connexes pour seuil bas = "<<ncc<<"\n";
	bool iccOK;
	if (imaR.nlig()<nblig || imaR.ncol()<nbcol) imaR=imadata<BYTE>(nblig,nbcol,nbcanaux); imaR.mise_a_zero();
	for (n=1; n<=ncc; n++) {
		iccOK=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (!iccOK && imacc(i,j)==n && ima_sh(i,j)) iccOK=1;
		if (iccOK) {
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) if (imacc(i,j)==n) imaR(i,j)=1;;
		}
	}
}

template <class T> imadata<float> imadata<T>::init1_segmentationregions (const int nmax_reg_init, const float delta_pas_quantif, int& nreg_init) {
	const float epsilon=(float)1.e-6;
	float pas_quantif=1.f, xx, x_it;
	int i, j, k, l;
	bool fini_init=0;
	nreg_init=nblig*nbcol;
	imadata<float> ima2(*this), imaHomoge(nblig,nbcol); for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) imaHomoge(i,j)=(float)i*nbcol+j;
	int *tab_np=new int[nblig*nbcol];
	imabin ima1reg(nblig,nbcol);
	do {
		x_it=0;
		do {
			l=0; x_it+=1;
			for (i=0; i<nblig; i++) 
				for (j=0; j<nbcol-1; j++)
					if (fabs(ima2(i,j)-ima2(i,j+1))<=epsilon && imaHomoge(i,j+1)!=imaHomoge(i,j)) {xx=mini(imaHomoge(i,j+1),imaHomoge(i,j)); imaHomoge(i,j+1)=imaHomoge(i,j)=xx; l++;}
			for (i=0; i<nblig-1; i++) 
				for (j=0; j<nbcol; j++) 
					if (fabs(ima2(i,j)-ima2(i+1,j))<=epsilon && imaHomoge(i+1,j)!=imaHomoge(i,j)) {xx=mini(imaHomoge(i+1,j),imaHomoge(i,j)); imaHomoge(i+1,j)=imaHomoge(i,j)=xx; l++;}
			if (((int)x_it)%100==0) cout<<" it "<<(int)x_it<<" l = "<<l<<"\n";
		} while (l>0 && x_it<10000);
		for (l=0; l<nblig*nbcol; l++) tab_np[l]=0;
		for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) tab_np[(int)imaHomoge(i,j)]++;
		k=0; for (l=0; l<nblig*nbcol; l++) if (tab_np[l]>0) k++;
		cout<<" a l'initialisation "<<k<<" regions de contraste <= "<<pas_quantif<<"\n";
//			if (k>nmax_reg_init) {
		if (k>nmax_reg_init*(1-0.25)) { // en fait calcul précédent seulement approximatif -> on prend 25% de marge d'erreur sur l'estimation du nombre de regions
			pas_quantif+=delta_pas_quantif; if (pas_quantif<1.f) cout<<" attention pas quatification < 1 ???????????\n";
			for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) {ima2(i,j)=floor((*this)(i,j)/pas_quantif);}
		} else {
			fini_init=1; nreg_init=0; // toujours parce que calcul précédent seulement approximatif -> recalcul exact du nombre de regions
			imaHomoge.mise_a_zero(); ima2.statbasic(1);
			int ng_min=(int)ima2.minI(), ng_max=(int)ima2.maxI(), ng; 
			imadata<int> imaCC; int ncc_ng;
			for (ng=ng_min; ng<=ng_max; ng++) {
				ima1reg.mise_a_zero();
				for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) if (fabs(ima2(i,j)-ng)<delta_pas_quantif) ima1reg(i,j)=1;
				if (ima1reg.norm()>0) {
					imaCC=ima1reg.composantes_connexes(ncc_ng,8);
					for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) if (imaCC(i,j)>0) imaHomoge(i,j)=(float)(imaCC(i,j)+nreg_init);
					nreg_init+=ncc_ng;
				}
			}
		}
	} while (!fini_init);
	if (tab_np!=NULL) delete[] tab_np;
	imaHomoge.sauve_Ima("imaHomoge.dat"); ima2.sauve_ImaPGM("ima_quantif.pgm"); //char aa; cin>>aa;
	return imaHomoge;
}

template <class T> imadata<float> imadata<T>::gradient (imadata<float> &ima_dir, char *masq) const {
	int taille_masq=3; if (strcmp(masq,"MDIF")==0) taille_masq=5;
	const int dim=taille_masq, dim2=dim/2;
	double *noyau=new double[dim*dim], coef=taille_masq, sum;
	imadata<float> ima_gdVH(nblig,nbcol,2), ima_grd(nblig,nbcol,nbcanaux);
	ima_dir=ima_grd;
	int dir,i,i0,i2,j,j0,j2,ii,jj,k;
	for (k=0; k<nbcanaux; k++) {
		for (dir=0; dir<2; dir++) {
			if (dir==0) {                                   //gradient horizontal
				for (i=0; i<dim*dim; i++) noyau[i]=0;
				if (strcmp(masq,"Diff1")==0) {i0=dim2*dim; noyau[i0]=1; noyau[i0+dim-1]=-1; coef=1.;}
				if (strcmp(masq,"Prewitt")==0 || strcmp(masq,"Sobel")==0) {
					for (i=0; i<dim; i++) {i0=i*dim; noyau[i0]=1; noyau[i0+dim-1]=-1;}
					if (strcmp(masq,"Sobel")==0) {i0=dim2*dim; noyau[i0]=2; noyau[i0+dim-1]=-2; coef=4.;}
        }
				if (strcmp(masq,"MDIF")==0)
					for (i=0; i<=dim2; i++) {i0=i*dim; i2=(dim-1-i)*dim; noyau[i0+1]=noyau[i2+1]=i+1; noyau[i0+dim-2]=noyau[i2+dim-2]=-i-1;}
			}
			else {                                          //gradient vertical
				for (i=0; i<dim*dim; i++) noyau[i]=0;
				i0=(dim-1)*dim;
				if (strcmp(masq,"Diff1")==0) {noyau[dim2]=1; noyau[i0+dim2]=-1; coef=1.;}
				if (strcmp(masq,"Prewitt")==0 || strcmp(masq,"Sobel")==0) {
					for (i=0; i<dim; i++) {noyau[i]=1; noyau[i0+i]=-1;}
					if (strcmp(masq,"Sobel")==0) {noyau[dim2]=2; noyau[i0+dim2]=-2; coef=4.;}
				}
				if (strcmp(masq,"MDIF")==0) {
					i0=dim; i2=(dim-2)*dim;
					for (i=0; i<=dim2; i++) {noyau[i0+i]=noyau[i0+dim-1-i]=i+1; noyau[i2+i]=noyau[i2+dim-1-i]=-i-1;}
				}
			}
//			for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<noyau[i*dim+j]<<" "; cout<<"\n";}
//			sum=0.; for (i=0; i<dim*dim; i++) sum+=abs(noyau[i]); sum/=2.;
//			for (i=0; i<dim*dim; i++) noyau[i]/=sum;
//			for (i=0; i<dim; i++) {for (j=0; j<dim; j++) cout<<noyau[i*dim+j]<<" "; cout<<"\n";}
			for (i=0; i<nblig; i++) {
				i0=maxi(0,i-dim2); i2=mini(i+dim2,nblig-1); ii=mini(i-i0,i2-i); i0=i-ii; i2=i+ii;
				for (j=0; j<nbcol; j++) {
					j0=maxi(0,j-dim2); j2=mini(j+dim2,nbcol-1); jj=mini(j-j0,j2-j); j0=j-jj; j2=j+jj;
					sum=0.;
					for (ii=i0; ii<=i2; ii++)
						for (jj=j0; jj<=j2; jj++) sum+=(*this)(ii,jj,k)*noyau[(ii-i+dim2)*dim+jj-j+dim2];
					ima_gdVH(i,j,dir)=(float)(sum/coef);
				}
			}
		}
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {                  //norme et orientation du gradient
				ima_grd(i,j,k)=pow(pow(ima_gdVH(i,j,0),2.f)+pow(ima_gdVH(i,j,1),2.f),0.5f);
//				if (ima_gdVH(i,j,0)==0.) ima_dir(i,j,k)=(float)(PI/2.);
//				else ima_dir(i,j,k)=atan(ima_gdVH(i,j,1)/ima_gdVH(i,j,0));
				ima_dir(i,j,k)=atan2(-ima_gdVH(i,j,1),ima_gdVH(i,j,0));
//				if (ima_grd(i,j,k)>0.001) cout<<i<<" "<<j<<" "<<ima_gdVH(i,j,0)<<" "<<ima_gdVH(i,j,1)<<" "<<ima_grd(i,j,k)<<" "<<ima_dir(i,j,k)<<"\n";
			}
	}
	cout<<" stat ima_gdVH\n"; ima_gdVH.statbasic(1);
	((imadata<float>(ima_gdVH,0)+4)*20).imaunsignedchar().sauve_ImaPGM("./ima_gdH.pgm");
	((imadata<float>(ima_gdVH,1)+4)*20).imaunsignedchar().sauve_ImaPGM("./ima_gdV.pgm");
	if (noyau!=NULL) {delete[] noyau; noyau=NULL;}
	ima_grd.statbasic(1); ima_dir.statbasic(1);
	return ima_grd;
}

template <class T> imadata<float> imadata<T>::laplacien(char *masqL) const {
	imadata<float> ima_Lpc(nblig,nbcol,nbcanaux);
// calcul de l'image du Laplacien
	double *noyauL=NULL, sum, coef;
	const int taille_masqL=3, dimL=taille_masqL, dimL2=dimL/2;
	noyauL=new double[dimL*dimL];
	int i,j,k,ii,jj,i0,i2,j0,j2;
	for (i=0; i<dimL*dimL; i++) noyauL[i]=-1;
	noyauL[dimL2*dimL+dimL2]=8;
	if (strcmp(masqL,"4connex")==0) {
		noyauL[0]=noyauL[dimL-1]=noyauL[(dimL-1)*dimL]=noyauL[dimL*dimL-1]=0;
		noyauL[dimL2*dimL+dimL2]=4;
	}
	for (k=0; k<nbcanaux; k++) {
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-dimL2); i2=mini(i+dimL2,nblig-1);
			ii=mini(i-i0,i2-i); i0=i-ii; i2=i+ii;
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-dimL2); j2=mini(j+dimL2,nbcol-1);
				jj=mini(j-j0,j2-j); j0=j-jj; j2=j+jj;
				sum=0;
				for (ii=i0; ii<=i2; ii++)
					for (jj=j0; jj<=j2; jj++) {
						coef=noyauL[(ii-i+dimL2)*dimL+jj-j+dimL2];
						sum+=(*this)(ii,jj,k)*coef;
					}
				ima_Lpc(i,j,k)=(float)sum;
			}
		}
	}
	char nomfich[80]="imaLpc_"; strcat_s(nomfich,masqL); strcat_s(nomfich,".pgm");
	ima_Lpc.sauve_ImaPGM(nomfich);
	cout<<" stat ima_laplacien\n"; ima_Lpc.statbasic(1);
	if (noyauL!=NULL) delete[] noyauL;
	return ima_Lpc;
}

template <class T> imadata<float> imadata<T>::gradientOpt (imadata<float> &ima_dir, float alpha, char *filtre) const {
	imadata<float> ima_grd(nblig,nbcol,nbcanaux); ima_dir=ima_grd;
  imadata<double> ima_gdVH(nblig,nbcol,2);
	double *vect_B1=new double[maxi(nblig,nbcol)];
	double *vect_B2=new double[maxi(nblig,nbcol)];
	double e_alpha=exp(-alpha), /*coef_s=1-e_alpha, */coef_c=-pow(1-e_alpha,2)/e_alpha, 
		   coef_b=pow(1-e_alpha,2)/(1+2*alpha*e_alpha-pow(e_alpha,2));
	double cflA=coef_b, cflAm1=coef_b*e_alpha*(alpha-1), cflBm1=2*e_alpha, cflBm2=-pow(e_alpha,2),
		   cflAp1=coef_b*e_alpha*(alpha+1), cflAp2=-coef_b*pow(e_alpha,2), cflBp1=cflBm1, 
		   cflBp2=cflBm2, cfdAm1=coef_c*e_alpha, cfdBm1=cflBm1, cfdBm2=cflBm2, cfdAp1=-cfdAm1,
		   cfdBp1=cflBm1, cfdBp2=cflBm2;
	int i,j,k;
	for (k=0; k<nbcanaux; k++) {
		if (strcmp(filtre,"Deriche")==0 || strcmp(filtre,"Canny-Deriche")==0) {
// gradient horizontal = lissage vertical puis derivation horizontale
			for (j=0; j<nbcol; j++) { // lissage vertical
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=cflA*(*this)(i,j,k)+cflAm1*(*this)(i-1,j,k)+cflBm1*vect_B1[i-1];
					if (i>1) vect_B1[i]+=cflBm2*vect_B1[i-2];
				}
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=cflAp1*(*this)(i+1,j,k)+cflBp1*vect_B2[i+1];
					if (i<nblig-2) vect_B2[i]+=cflAp2*(*this)(i+2,j,k)+cflBp2*vect_B2[i+2];
				}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,0)=vect_B1[i]+vect_B2[i];
			}
			for (i=0; i<nblig; i++) { // derivation horizontale
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=cfdAm1*ima_gdVH(i,j-1,0)+cfdBm1*vect_B1[j-1];
					if (j>1) vect_B1[j]+=cfdBm2*vect_B1[j-2];
				}
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=cfdAp1*ima_gdVH(i,j+1,0)+cfdBp1*vect_B2[j+1];
					if (j<nbcol-2) vect_B2[j]+=cfdBp2*vect_B2[j+2];
				}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,0)=vect_B1[j]+vect_B2[j];
			}
// gradient vertical = lissage horizontal puis derivation verticale
			for (i=0; i<nblig; i++) { // lissage horizontal
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=cflA*(*this)(i,j,k)+cflAm1*(*this)(i,j-1,k)+cflBm1*vect_B1[j-1];
					if (j>1) vect_B1[j]+=cflBm2*vect_B1[j-2];
				}
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=cflAp1*(*this)(i,j+1,k)+cflBp1*vect_B2[j+1];
					if (j<nbcol-2) vect_B2[j]+=cflAp2*(*this)(i,j+2,k)+cflBp2*vect_B2[j+2];
				}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,1)=vect_B1[j]+vect_B2[j];
			}
			for (j=0; j<nbcol; j++) { // derivation verticale
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=cfdAm1*ima_gdVH(i-1,j,1)+cfdBm1*vect_B1[i-1];
					if (i>1) vect_B1[i]+=cfdBm2*vect_B1[i-2];
				}
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=cfdAp1*ima_gdVH(i+1,j,1)+cfdBp1*vect_B2[i+1];
					if (i<nblig-2) vect_B2[i]+=cfdBp2*vect_B2[i+2];
				}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,1)=vect_B1[i]+vect_B2[i];
			}
		}
		/*else*/ if (strcmp(filtre,"Shen")==0 || strcmp(filtre,"Shen-Castan")==0) {
// gradient horizontal
			for (i=0; i<nblig; i++) { // derivation horizontale
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=e_alpha*(vect_B1[j-1]-(*this)(i,j,k))+(*this)(i,j,k);
				}
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=e_alpha*(vect_B2[j+1]-(*this)(i,j,k))+(*this)(i,j,k);
				}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,0)=vect_B1[j]-vect_B2[j];
			}
			for (j=0; j<nbcol; j++) { // lissage vertical
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=e_alpha*(vect_B1[i-1]-ima_gdVH(i,j,0))+ima_gdVH(i,j,0);
				}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,0)=vect_B1[i];
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=e_alpha*(vect_B2[i+1]-ima_gdVH(i,j,0))+ima_gdVH(i,j,0);
				}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,0)=vect_B2[i];
			}
// gradient vertical
			for (j=0; j<nbcol; j++) { // derivation verticale
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=e_alpha*(vect_B1[i-1]-(*this)(i,j,k))+(*this)(i,j,k);
				}
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=e_alpha*(vect_B2[i+1]-(*this)(i,j,k))+(*this)(i,j,k);
				}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,1)=vect_B1[i]-vect_B2[i];
			}
			for (i=0; i<nblig; i++) { // lissage horizontal
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=e_alpha*(vect_B1[j-1]-ima_gdVH(i,j,1))+ima_gdVH(i,j,1);
				}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,1)=vect_B1[j];
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=e_alpha*(vect_B2[j+1]-ima_gdVH(i,j,1))+ima_gdVH(i,j,1);
				}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,1)=vect_B2[j];
			}
		}
// norme et direction du gradient
//	cout<<" calcul norme et orientation gradient\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ima_grd(i,j,k)=(float)(pow(pow(ima_gdVH(i,j,0),2.)+pow(ima_gdVH(i,j,1),2.),0.5));
				if (ima_gdVH(i,j,0)==0.) ima_dir(i,j,k)=(float)(PI/2.);
				else ima_dir(i,j,k)=(float)atan(ima_gdVH(i,j,1)/ima_gdVH(i,j,0));
			}
		cout<<" stat ima_gdVH\n"; ima_gdVH.statbasic(1);
		((imadata<float>(ima_gdVH,0)+4)*20).imaunsignedchar().sauve_ImaPGM("./ima_gdH.pgm");
		((imadata<float>(ima_gdVH,1)+4)*20).imaunsignedchar().sauve_ImaPGM("./ima_gdV.pgm");
	}
	if (vect_B1!=NULL) delete[] vect_B1;
	if (vect_B2!=NULL) delete[] vect_B2;
	ima_grd.statbasic(1); ima_dir.statbasic(1);
	return ima_grd;
}

template <class T> imadata<float> imadata<T>::steerablefilter (imadata<float> &ima_grd_theta, float sigma) const {
//	imadata<float> ima_filtreGauss=(*this).filtregaussienne(sigma);
//	ima_filtreGauss.imaunsignedchar().sauve_ImaPGM("ima_filtreGauss.pgm");
	int i,j,k,i0,i2,j0,j2,ii,jj,dim2=int(2*sigma),dim=2*dim2+1;
	float sigma2=2*sigma*sigma, fact=2.f/sigma2;
	double **noyauX2=new double*[dim], **noyauXY=new double*[dim], sumX2, sumY2, sumXY, xxx;
	for (j=0; j<dim; j++) {noyauX2[j]=new double[dim]; noyauXY[j]=new double[dim];}
	for (i=0; i<dim; i++) 
		for (j=0; j<dim; j++) {
			xxx=exp(-fact*((j-dim2)*(j-dim2)+(i-dim2)*(i-dim2))/2.)*fact/2/PI;
			noyauX2[i][j]=xxx*fact*((j-dim2)*(j-dim2)*fact-1.); 
			noyauXY[i][j]=xxx*fact*fact*(j-dim2)*(i-dim2);
		}
	imadata<T> imaRes(nblig,nbcol,nbcanaux*3);
	for (k=0; k<nbcanaux; k++)
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-dim2); i2=mini(i+dim2,nblig-1);
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-dim2); j2=mini(j+dim2,nbcol-1);
				sumX2=sumY2=sumXY=0.;
				for (ii=i0; ii<=i2; ii++)
					for (jj=j0; jj<=j2; jj++) {
						sumY2+=(*this)(ii,jj,k)*noyauX2[ii-i+dim2][jj-j+dim2]; 
						sumX2+=(*this)(ii,jj,k)*noyauX2[jj-j+dim2][ii-i+dim2]; 
						sumXY+=(*this)(ii,jj,k)*noyauXY[ii-i+dim2][jj-j+dim2];
					}
				imaRes(i,j,k)=(float)sumX2;
				imaRes(i,j,nbcanaux+k)=(float)sumY2;
				imaRes(i,j,nbcanaux*2+k)=(float)sumXY;
			}
		}
	imaRes.sauve_ImaBSQ("ima_steerablefilterbasis.dat"); 
	int n_dir=ima_grd_theta.ncanaux()/nbcanaux, n; //cout<<" n_dir = "<<n_dir<<" ima_grd_theta.ncanaux() = "<<ima_grd_theta.ncanaux()<<" nbcanaux = "<<nbcanaux<<"\n";
	float coef[3], theta, xx;
	for (n=0; n<n_dir; n++) {
		switch (n) {
			case 0: theta=0.f; break; 
			case 1: theta=(float)PI/2.f; break; 
			case 2: theta=(float)PI/4.; break; 
			case 3: theta=3*(float)PI/4.f; break; 
			case 4: theta=(float)PI/8.f; break; 
			case 5: theta=3*(float)PI/8.f; break; 
			case 6: theta=5*(float)PI/8.f; break; 
			case 7: theta=7*(float)PI/8.f; break; 
			default: theta=0.f; 
		} cout<<" + direction "<<theta/PI*180<<"\n";
		coef[0]=cos(theta)*cos(theta); coef[1]=-sin(2*theta); coef[2]=sin(theta)*sin(theta);
		for (k=0; k<nbcanaux; k++)
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					xx=coef[0]*imaRes(i,j,k)+coef[2]*imaRes(i,j,nbcanaux+k)+coef[1]*imaRes(i,j,2*nbcanaux+k);
					ima_grd_theta(i,j,nbcanaux*n+k)=xx;
				}
	}
/*	for (j=0; j<dim; j++) noyau[j]=fact*exp(-(j-dim2)*(j-dim2)/sigma2);
	float fact=1/sqrt(2*(float)PI)/sigma;
	double sum, coef;
	if (affich) {
		cout<<" noyau de dimension "<<dim<<"\n";
		coef=noyau[dim2]; norm=0.;
		for (j=0; j<dim; j++) {
			cout<<setw(5)<<1/1000.*around(1000*noyau[j]/coef)<<"  ";
			norm+=(float)noyau[j];
		}
		cout<<" norme = "<<norm<<", norme/max = "<<norm/coef<<"\n";
		for (j=0; j<dim; j++) cout<<setw(5)<<around(noyau[j]*coef/noyau[0]*10)<<"  ";
		cout<<" norme = "<<norm*coef/noyau[0]*10<<"\n";
	}*/
	for (j=0; j<dim; j++) {if (noyauX2[j]!=NULL) delete[] noyauX2[j]; if (noyauXY[j]!=NULL) delete[] noyauXY[j];}
	if (noyauX2!=NULL) delete[] noyauX2; if (noyauXY!=NULL) delete[] noyauXY;
	return imaRes;
}

template <class T> imadata<BYTE> imadata<T>::chemin_opt4connex (imadata<float> &ima_cost, int y0, int x0, int y2, int x2, float cost_min) const {
	cout<<" debut recherche chemin_opt4connex entre ("<<y0<<","<<x0<<") et ("<<y2<<","<<x2<<")\n";
	imadata<unsigned long int> ima_prev(nblig,nbcol);
	int i,j;
	for (i=0; i<nblig; ++i)
		for (j=0; j<nbcol; ++j) {ima_prev(i,j)=nblig*nbcol; ima_cost(i,j)=FLT_MAX;}
	ima_cost(y0,x0)=0;
	const int nmax=10000, kconnex=4, Tvois[kconnex][2]={{-1,0},{0,-1},{+1,0},{0,+1}}, itermax=nmax;
	int Tcoord[nmax][2], n=-1, iter=0, x=x0, y=y0, dx, dy, xs, ys, k, xx,yy, next; 
	float cost, lcost;
	bool fini=false;
// initialisation avec chemin 4-connexe de coût min parmi les chemins de longueur minimale
	Tcoord[++n][0]=y0; Tcoord[n][1]=x0; cout<<" n = "<<n<<" "<<Tcoord[n][0]<<" "<<Tcoord[n][1]<<"\n";
	if (cost_min<0) {
		cost_min=ima_cost(y2,x2);
		do {
//			if (abs(y2-y)>=abs(x2-x)) {dx=0; dy=(y2-y)/abs(y2-y);} else {dy=0; dx=(x2-x)/abs(x2-x);}
			dx=dy=0;
			if (abs(y2-y)>0) dy=(y2-y)/abs(y2-y); if (abs(x2-x)>0) dx=(x2-x)/abs(x2-x);
			if (dy!=0 && dx!=0) if ((*this)(y+dy,x)<=(*this)(y,x+dx)) dx=0; else dy=0;
			ys=y+dy; xs=x+dx; 
//			cout<<" y2-y = "<<y2-y<<" x2-x = "<<x2-x<<" dy = "<<dy<<" dx = "<<dx<<" y = "<<y<<" x = "<<x<<" ys = "<<ys<<" xs = "<<xs<<"\n";
			cost=ima_cost(y,x)+(*this)(ys,xs); 
			if (ys>=0 && ys<nblig && xs>=0 && xs<nbcol && cost<cost_min && cost<ima_cost(ys,xs)) {
				ima_cost(ys,xs)=cost; ima_prev(ys,xs)=y*nbcol+x;
				Tcoord[++n][0]=ys; Tcoord[n][1]=xs;
				if (xs==x2 && ys==y2) {cout<<" extremite finale atteinte\n";
					cost_min=cost; fini=true;
					for (j=0; j<=n; ++j) {
						yy=Tcoord[j][0]; xx=Tcoord[j][1];
						if (yy>=0 && yy<nblig && xx>=0 && xx<nbcol && ima_cost(yy,xx)>cost_min) Tcoord[j][0]=Tcoord[j][1]=-1;
					}
				}
			}
			y=Tcoord[n][0]; x=Tcoord[n][1]; 
		} while (!fini);
		cout<<" par chemin direct de longueur "<<n<<" : cout = "<<cost_min<<"\n";
		ima_cost.sauve_Ima("ima_cost_init.dat");
	}
	next=0;																													// on repart de l'origine
	do {
		y=Tcoord[next][0]; x=Tcoord[next][1]; 
		Tcoord[next][0]=Tcoord[next][1]=-1;
		for (k=0; k<kconnex; ++k) {
			ys=y+Tvois[k][0]; xs=x+Tvois[k][1];
			if (ys>=0 && ys<nblig && xs>=0 && xs<nbcol) {
				cost=ima_cost(y,x)+(*this)(ys,xs); 
				if (cost<cost_min && cost<ima_cost(ys,xs)) {
					ima_cost(ys,xs)=cost; ima_prev(ys,xs)=y*nbcol+x;
					Tcoord[++n][0]=ys; Tcoord[n][1]=xs;
					if (xs==x2 && ys==y2) {cout<<" extremite finale atteinte\n";
						cost_min=cost;
						for (j=0; j<=n; ++j) {
							yy=Tcoord[j][0]; xx=Tcoord[j][1];
							if (yy>=0 && yy<nblig && xx>=0 && xx<nbcol && ima_cost(yy,xx)>cost_min) Tcoord[j][0]=Tcoord[j][1]=-1;
						}
					}
				}
			}
		}
		if (n>nmax-kconnex-1 || iter%(nmax/10)==0) {//cout<<" iteration "<<iter<<" nettoyage du tableau Tcoord\n";
			for (k=0; k<=n; ++k)
				if (Tcoord[k][0]==-1 || Tcoord[k][1]==-1) {
					for (j=k; j<n; ++j) {Tcoord[j][0]=Tcoord[j+1][0]; Tcoord[j][1]=Tcoord[j+1][1];}
					--n;
				}
		}
		fini=true;
/*		while (n>=0 && fini) {y=Tcoord[n][0]; x=Tcoord[n][1]; if (y>=0 && y<nblig && x>=0 && x<nbcol) fini=false; else --n; }*/
		next=0; lcost=cost_min;
		for (j=0; j<=n; ++j) {
			y=Tcoord[j][0]; x=Tcoord[j][1]; 
			if (y>=0 && y<nblig && x>=0 && x<nbcol && ima_cost(y,x)<lcost) {next=j; lcost=ima_cost(y,x); }
		}
		if (lcost<cost_min) fini=false; //cout<<" lcost "<<lcost<<" y = "<<Tcoord[next][0]<<" x = "<<Tcoord[next][1]<<"\n";
		if (iter%1000==0) {cout<<" iter "<<++iter<<" n = "<<n<<"\n"; ima_cost.sauve_Ima("ima_cost.dat"); /*char aa; cin>>aa;*/}
	} while (!fini);
	cout<<" fin recherche de chemin\n"; ima_cost.sauve_Ima("ima_cost.dat");
	imadata<BYTE> ima_chemin(nblig,nbcol);  for (j=0; j<nblig; j++) for (k=0; k<nbcol; k++) ima_chemin(j,k)=1; //ima_chemin.mise_a_zero(); 
	ima_chemin(y0,x0)=ima_chemin(y2,x2)=255;
	y=y2; x=x2; int lg=0;
	do {//cout<<" ("<<y<<","<<x<<") .... ";
		if (y>=0 && y<nblig && x>=0 && x<nbcol) {
			ima_chemin(y,x)=255; lg++;
			yy=ima_prev(y,x)/nbcol; x=ima_prev(y,x)%nbcol; y=yy;
		}
	} while (x>=0 && x<nbcol && y>=0 && y<nblig); 
	ima_chemin.statbasic(1);
//	ima_chemin.sauve_ImaPGM("ima_chemin_tmp.pgm");
	cout<<" fin chemin_opt4connex de cout "<<cost_min<<" et de longueur "<<lg<<"\n";
	return ima_chemin;
}

template <class T> imadata<BYTE> imadata<T>::chemin_opt4connex (imadata<float> &ima_cost, int y0, int x0, const imabin &Iobj, int &y2, int &x2, float cost_min) const {
	cout<<" debut recherche chemin_opt4connex entre ("<<y0<<","<<x0<<") et un ensemble de pixels\n";
	imadata<unsigned long int> ima_prev(nblig,nbcol);
	int i,j;
	for (i=0; i<nblig; ++i)
		for (j=0; j<nbcol; ++j) {ima_prev(i,j)=nblig*nbcol; ima_cost(i,j)=FLT_MAX;}
	ima_cost(y0,x0)=0;
	const int nmax=10000, kconnex=4, Tvois[kconnex][2]={{-1,0},{0,-1},{+1,0},{0,+1}}, itermax=nmax;
	int Tcoord[nmax][2], n=-1, iter=0, x=x0, y=y0, dx, dy, xs, ys, k, xx,yy, next; 
	float cost, lcost, d2, d2min;
	bool fini=false;
// initialisation avec chemin 4-connexe de coût min parmi les chemins de longueur minimale
	Tcoord[++n][0]=y0; Tcoord[n][1]=x0; cout<<" n = "<<n<<" "<<Tcoord[n][0]<<" "<<Tcoord[n][1]<<"\n";
	if (cost_min<0) {
		imadata<float> imadist_Iobj=Iobj.Tr_dist();
		do {
			d2=imadist_Iobj(y,x); 
			if (x>0 && y>0 && x<nbcol-1 && y<nblig-1) d2min=mini(imadist_Iobj(y-1,x),mini(imadist_Iobj(y+1,x),mini(imadist_Iobj(y,x-1),imadist_Iobj(y,x+1))));
			else {
				d2min=FLT_MAX;
				if (y>0) d2min=mini(d2min,imadist_Iobj(y-1,x)); if (y<nblig-1) d2min=mini(d2min,imadist_Iobj(y+1,x)); 
				if (x>0) d2min=mini(d2min,imadist_Iobj(y,x-1)); if (x<nbcol-1) d2min=mini(d2min,imadist_Iobj(y,x+1)); 		
			}
			if (d2min<=d2) {
				if (y>0 && d2min==imadist_Iobj(y-1,x)) y-=1; 
				else 
					if (y<nblig-1 && d2min==imadist_Iobj(y+1,x)) y+=1; 
					else 
						if (x>0 && d2min==imadist_Iobj(y,x-1)) x-=1; 
						else if (x<nbcol-1 && d2min==imadist_Iobj(y,x+1)) x+=1; 
			} else {cout<<" Pb dans calcul de point le + proche sur l'ensemble des points objectifs\n";}
		} while (imadist_Iobj(y,x)>0);
		y2=y; x2=x; y=y0; x=x0; // (y2,x2) est le point le + proche (en distance) sur la CC
		cost_min=ima_cost(y2,x2);
		do {
			dx=dy=0;
			if (abs(y2-y)>0) dy=(y2-y)/abs(y2-y); if (abs(x2-x)>0) dx=(x2-x)/abs(x2-x);
			if (dy!=0 && dx!=0) if ((*this)(y+dy,x)<=(*this)(y,x+dx)) dx=0; else dy=0;
			ys=y+dy; xs=x+dx; 
			cost=ima_cost(y,x)+(*this)(ys,xs); 
			if (ys>=0 && ys<nblig && xs>=0 && xs<nbcol && cost<cost_min && cost<ima_cost(ys,xs)) {
				ima_cost(ys,xs)=cost; ima_prev(ys,xs)=y*nbcol+x;
				Tcoord[++n][0]=ys; Tcoord[n][1]=xs;
				if (xs==x2 && ys==y2) {//cout<<" extremite finale atteinte\n";
					cost_min=cost; fini=true;
					for (j=0; j<=n; ++j) {
						yy=Tcoord[j][0]; xx=Tcoord[j][1];
						if (yy>=0 && yy<nblig && xx>=0 && xx<nbcol && ima_cost(yy,xx)>cost_min) Tcoord[j][0]=Tcoord[j][1]=-1;
					}
				}
			}
			y=Tcoord[n][0]; x=Tcoord[n][1]; 
		} while (!fini);
		cout<<" par chemin direct de longueur "<<n<<" : cout = "<<cost_min<<"\n";
		ima_cost.sauve_Ima("ima_cost_init.dat");
	}
	next=0;																													// on repart de l'origine
	do {
		y=Tcoord[next][0]; x=Tcoord[next][1]; 
		Tcoord[next][0]=Tcoord[next][1]=-1;
		for (k=0; k<kconnex; ++k) {
			ys=y+Tvois[k][0]; xs=x+Tvois[k][1];
			if (ys>=0 && ys<nblig && xs>=0 && xs<nbcol) {
				cost=ima_cost(y,x)+(*this)(ys,xs); 
				if (cost<cost_min && cost<ima_cost(ys,xs)) {
					ima_cost(ys,xs)=cost; ima_prev(ys,xs)=y*nbcol+x;
					Tcoord[++n][0]=ys; Tcoord[n][1]=xs;
//					if (xs==x2 && ys==y2) {cout<<" extremite finale atteinte\n";
					if (Iobj(ys,xs)) {cout<<" objectif final atteint\n";
						cost_min=cost; x2=xs; y2=ys; // (y2,x2) est le point sur la CC atteint par chemin de coût min
						for (j=0; j<=n; ++j) {
							yy=Tcoord[j][0]; xx=Tcoord[j][1];
							if (yy>=0 && yy<nblig && xx>=0 && xx<nbcol && ima_cost(yy,xx)>cost_min) Tcoord[j][0]=Tcoord[j][1]=-1;
						}
					}
				}
			}
		}
		if (n>nmax-kconnex-1 || iter%(nmax/10)==0) {//cout<<" iteration "<<iter<<" nettoyage du tableau Tcoord\n";
			for (k=0; k<=n; ++k)
				if (Tcoord[k][0]==-1 || Tcoord[k][1]==-1) {
					for (j=k; j<n; ++j) {Tcoord[j][0]=Tcoord[j+1][0]; Tcoord[j][1]=Tcoord[j+1][1];}
					--n;
				}
		}
		fini=true;
/*		while (n>=0 && fini) {y=Tcoord[n][0]; x=Tcoord[n][1]; if (y>=0 && y<nblig && x>=0 && x<nbcol) fini=false; else --n; }*/
		next=0; lcost=cost_min;
		for (j=0; j<=n; ++j) {
			y=Tcoord[j][0]; x=Tcoord[j][1]; 
			if (y>=0 && y<nblig && x>=0 && x<nbcol && ima_cost(y,x)<lcost) {next=j; lcost=ima_cost(y,x); }
		}
		if (lcost<cost_min) fini=false; //cout<<" lcost "<<lcost<<" y = "<<Tcoord[next][0]<<" x = "<<Tcoord[next][1]<<"\n";
		if (iter%1000==0) {cout<<" iter "<<++iter<<" n = "<<n<<"\n"; ima_cost.sauve_Ima("ima_cost.dat"); /*char aa; cin>>aa;*/}
	} while (!fini);
	cout<<" fin recherche de chemin\n"; ima_cost.sauve_Ima("ima_cost.dat");
	imadata<BYTE> ima_chemin(nblig,nbcol);  for (j=0; j<nblig; j++) for (k=0; k<nbcol; k++) ima_chemin(j,k)=1; //ima_chemin.mise_a_zero(); 
	ima_chemin(y0,x0)=ima_chemin(y2,x2)=255;
	y=y2; x=x2; int lg=0;
	do {//cout<<" ("<<y<<","<<x<<") .... ";
		if (y>=0 && y<nblig && x>=0 && x<nbcol) {
			ima_chemin(y,x)=255; lg++;
			yy=ima_prev(y,x)/nbcol; x=ima_prev(y,x)%nbcol; y=yy;
		}
	} while (x>=0 && x<nbcol && y>=0 && y<nblig); 
	ima_chemin.statbasic(1);
//	ima_chemin.sauve_ImaPGM("ima_chemin_tmp.pgm");
	cout<<" fin chemin_opt4connex de cout "<<cost_min<<" et de longueur "<<lg<<"\n";
	return ima_chemin;
}

/* -----------------------------------------------------------
Calcul des images de parametres de texture calcules dans des fenetres glissantes
K : demi-longeur de la fenetre
N: nb bins dans l'histogramme
  ----------------------------------------------------------*/

/* -----------------------------------------------------------
Contraste
 ----------------------------------------------------------*/
template <class T> imadata<float> imadata<T>::ima_contrast (int K, int N, bool iaff) {
	imadata<float> ima(nblig,nbcol,nbcanaux);
	if (N<=0) N=(2*K+1)/*K*K/4*/;
	matrice1D<float> h(N);
	int i,j,k,ii,jj,n,i0,i2,j0,j2;
	double maxval;
	const float eps=1.e-6f, increment=1.f/pow(2*K+1.0f,2);
	for (k=0; k<nbcanaux; k++) {
		maxval=maxI(k);
		cout<<" canal "<<k<<" maxval = "<<maxval<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (n=0; n<N; n++) h(n)=0;
				i0=maxi(0,i-K); i2=mini(nblig-1,i+K);
				j0=maxi(0,j-K); j2=mini(nbcol-1,j+K);
				for (ii=i0; ii<=i2; ii++)
					for(jj=j0; jj<=j2; jj++) {
						h((int)((*this)(ii,jj,k)/(maxval+eps)*N))+=increment;}
				ima(i,j,k)=(float)h.contrast();
			}
	}
	return ima;
}
/* -----------------------------------------------------------
Entropie
 ----------------------------------------------------------*/
template <class T> imadata<float> imadata<T>::ima_entropy (int K, int N, bool iaff) {
	imadata<float> ima(nblig,nbcol,nbcanaux); 
	if (N<=0) N=(2*K+1)/*K*K/4*/;
	cout<<" N = "<<N<<" -> log(N) = "<<log((float)N)/log(2.f)<<"\n";
	matrice1D<float> h(N);
	int i,j,k,ii,jj,n,i0,i2,j0,j2;
	double maxval; //statbasic();
	const double eps=1.e-6, increment=1./pow(2*K+1.,2);
	for (k=0; k<nbcanaux; k++) { 
		maxval=maxI(k); 
		for (i=0; i<nblig; i++) {
			for (j=0; j<nbcol; j++) { 
				for (n=0; n<N; n++) h(n)=0;
				i0=maxi(0,i-K); i2=mini(nblig-1,i+K);
				j0=maxi(0,j-K); j2=mini(nbcol-1,j+K);
				for (ii=i0; ii<=i2; ii++)
					for(jj=j0; jj<=j2; jj++)
						h((int)((*this)(ii,jj,k)/(maxval+eps)*N))+=(float)increment;
				ima(i,j,k)=(float)h.entropy();
			}
		}
	}
	const double increment_g=1./(nblig*nbcol);
	matrice1D<float> H256(256), H128(128), H64(64);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) { 
			H256((int)((*this)(i,j)/256.001f*256))+=(float)increment_g;
			H128((int)((*this)(i,j)/256.001f*128))+=(float)increment_g;
			H64((int)((*this)(i,j)/256.001f*64))+=(float)increment_g;
		}
	cout<<" Entropie globale calculee sur histo 256 bins = "<<(float)H256.entropy()<<"\n";
	cout<<" Entropie globale calculee sur histo 128 bins = "<<(float)H128.entropy()<<"\n";
	cout<<" Entropie globale calculee sur histo  64 bins = "<<(float)H64.entropy()<<"\n";
	return ima;
}
/* -----------------------------------------------------------
Energie
 ----------------------------------------------------------*/
template <class T> imadata<float> imadata<T>::ima_energy (int K, int N, bool iaff) {
	imadata<float> ima(nblig,nbcol,nbcanaux);
	if (N<=0) N=(2*K+1)/*K*K/4*/;
	matrice1D<double> h(N);
	int i,j, k, ii, jj, n, i0, i2, j0, j2;
	double maxval;
	const double eps=1.e-6, increment=1./pow(2*K+1.0f,2);
	for (k=0; k<nbcanaux; k++) {
		maxval=maxI(k);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (n=0; n<N; n++) h(n)=0;
				i0=maxi(0,i-K); i2=mini(nblig-1,i+K);
				j0=maxi(0,j-K); j2=mini(nbcol-1,j+K);
				for (ii=i0; ii<=i2; ii++)
					for(jj=j0; jj<=j2; jj++)
						h((int)((*this)(ii,jj,k)/(maxval+eps)*N))+=increment;
				ima(i,j,k)=(float)h.energy();
			}
	}
//	ima.statbasic(1);
	return ima;
}
/* -----------------------------------------------------------
Moment centre d'ordre k
 ----------------------------------------------------------*/
template <class T> imadata<float> imadata<T>::ima_moment (int m, int K, int N) {
	imadata<float> ima(nblig,nbcol,nbcanaux);
	if (N<=0) N=(2*K+1)/*K*K/4*/;
	matrice1D<float> h(N);
	int i,j,k,ii,jj,n,i0,i2,j0,j2;
	double maxval;
	const float eps=1.e-6f, increment=1.f/pow(2*K+1.f,2);
	for (k=0; k<nbcanaux; k++) {
		maxval=maxI(k);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (n=0; n<N; n++) h(n)=0;
				i0=maxi(0,i-K); i2=mini(nblig-1,i+K);
				j0=maxi(0,j-K); j2=mini(nbcol-1,j+K);
				for (ii=i0; ii<=i2; ii++)
					for(jj=j0; jj<=j2; jj++)
						h((int)((*this)(ii,jj,k)/(maxval+eps)*N))+=increment;
				ima(i,j,k)=(float)h.moment(m);
			}
	}
	return(ima);
}

template <class T> imadata<float> imadata<T>::cooccurence (int delta_x, int delta_y, int N, bool iaff) {
	if (N<=0) N=32; cout<<" matrice de coocurrence de taille "<<N<<"x"<<N<<"\n"; 
	matrice2D<float> coocc(N,N);
	int i,j,k,i0,i2,j0,j2,ii,jj;
	double maxval;
	for (k=0; k<nbcanaux; k++) {
		maxval=256/*maxI(k)*/;
//		cout<<" canal "<<k<<" maxval = "<<maxval<<"\n";
		for (ii=0; ii<N; ii++) for (jj=0; jj<N; jj++) coocc(ii,jj)=0;
		i0=maxi(0,-delta_x); i2=mini(nblig-1,nblig-1-delta_x);
		j0=maxi(0,-delta_y); j2=mini(nbcol-1,nbcol-1-delta_y);
		const float eps=1.e-6f, increment=/*1*/255.f/(i2-i0+1.f)/(j2-j0+1.f); 
		for (i=i0; i<i2; i++)
			for (j=j0; j<j2; j++)
					coocc((int)((*this)(i,j,k)/(maxval+eps)*N),(int)((*this)(i+delta_x,j+delta_y,k)/(maxval+eps)*N))+=increment;
	} 
	imadata<float> ima(coocc);
	cout<<" matrice de cooccurence :\n"; ima.statbasic(1);
	return ima;
}

template <class T> imadata<float> imadata<T>::ima_carac_cooccurence (int K, int delta_x, int delta_y, unsigned char icarac, int N, bool iaff) {
	const float eps=1.e-6f, increment=1.f/pow(2*K+1.0f,2);
	imadata<float> ima(nblig,nbcol,nbcanaux);
	if (N<=0) N=(2*K+1)/*K*K/4*/; cout<<" matrice de coocurrence de taille "<<N<<"x"<<N<<"\n"; cout<<" caracteristique calculee "<<(int)icarac<<"\n";
	matrice2D<float> coocc(N,N);
	int i,j,k,i0,i2,j0,j2,ii,jj;
	double maxval;
	for (k=0; k<nbcanaux; k++) {
		maxval=maxI(k);
		cout<<" canal "<<k<<" maxval = "<<maxval<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (ii=0; ii<N; ii++) for (jj=0; jj<N; jj++) coocc(ii,jj)=0;
				i0=maxi(maxi(0,-delta_x),i-K); i2=mini(mini(nblig-1,nblig-1-delta_x),i+K);
				j0=maxi(maxi(0,-delta_y),j-K); j2=mini(mini(nbcol-1,nbcol-1-delta_y),j+K);
				for (ii=i0; ii<=i2; ii++)
					for(jj=j0; jj<=j2; jj++) {
						coocc((int)((*this)(ii,jj,k)/(maxval+eps)*N),(int)((*this)(ii+delta_x,jj+delta_y,k)/(maxval+eps)*N))+=increment;}
				switch (icarac) {
					case 0: ima(i,j,k)=(float)coocc.entropy(); break;
					case 1: ima(i,j,k)=(float)coocc.contrast(); break;
					case 2: ima(i,j,k)=(float)coocc.dissimilarity(); break;
					case 3: ima(i,j,k)=(float)coocc.homogeneity(); break;
					default: break;
				}
			}
	}
	cout<<" image de la caracteristique choisie calculee a partir de la matrice de cooccurence :\n"; ima.statbasic(1);
	return ima;
}

template <class T> imadata<float> imadata<T>::ima1k_LBP(int k) const {
	imadata<float> ima(nblig,nbcol);
	int i,j;
	T xx;
	for (i=1; i<nblig-1; i++)
		for (j=1; j<nbcol-1; j++) {
			xx=(*this)(i,j,k);
			ima(i,j)=pow(2.,0)*((*this)(i-1,j-1,k)>=xx?+1:-1)+pow(2.,1)*((*this)(i-1,j,k)>=xx?+1:-1)+pow(2.,2)*((*this)(i-1,j+1,k)>=xx?+1:-1)
							+pow(2.,3)*((*this)(i,j+1,k)>=xx?+1:-1)+pow(2.,4)*((*this)(i+1,j+1,k)>=xx?+1:-1)+pow(2.,5)*((*this)(i+1,j,k)>=xx?+1:-1)
							+pow(2.,6)*((*this)(i+1,j-1,k)>=xx?+1:-1)+pow(2.,7)*((*this)(i,j-1,k)>=xx?+1:-1);
		}
	return(ima);
}

template <class T> imadata<BYTE> imadata<T>::ima1k_lenghtLBP(int k) const {
	imadata<BYTE> ima(nblig,nbcol);
	int i,j,n,l,lmax;
	T xx;
	const int nneighbors=8;
	bool Tnsup[nneighbors];
	for (i=1; i<nblig-1; i++)
		for (j=1; j<nbcol-1; j++) {
			xx=(*this)(i,j,k);
			for (n=0; n<3; n++) Tnsup[n]=((*this)(i-1,j-1+n,k)>=xx?true:false);
			Tnsup[3]=((*this)(i,j+1,k)>=xx?true:false);
			for (n=4; n<7; n++) Tnsup[n]=((*this)(i+1,j+1-(n-4),k)>=xx?true:false);
			Tnsup[7]=((*this)(i,j-1,k)>=xx?true:false);
			l=lmax=0;
			for (n=0; n<2*nneighbors; n++) 
				if (Tnsup[(n+1)%nneighbors]==Tnsup[n%nneighbors]) {l++; if (l>lmax) lmax=l;}
				else l=0;
			ima(i,j)=mini(lmax,nneighbors);
		}
	return(ima);
}

template <class T> imadata<BYTE> imadata<T>::ima1k_lenghtLTP(int k) const { cout<<" calcul de ima1k_lenghtLTP\n";
	imadata<BYTE> ima(nblig,nbcol);
	int i,j,n,l,lmax;
	T xx, yy;
	const int nneighbors=8;
	BYTE Tnsup[nneighbors], zz;
	const float delta=1.f;
	for (i=1; i<nblig-1; i++)
		for (j=1; j<nbcol-1; j++) {
			xx=(*this)(i,j,k);
			for (n=0; n<8; n++) Tnsup[n]=1;
			for (n=0; n<3; n++) {
				yy=(*this)(i-1,j-1+n,k);
				if (yy>xx+delta) Tnsup[n]=2; else if (yy<xx-delta) Tnsup[n]=0;}
			n=3; yy=(*this)(i,j+1,k); if (yy>xx+delta) Tnsup[n]=2; else if (yy<xx-delta) Tnsup[n]=0;
			for (n=4; n<7; n++) {
				yy=(*this)(i+1,j+1-(n-4),k);
				if (yy>xx+delta) Tnsup[n]=2; else if (yy<xx-delta) Tnsup[n]=0;}
			n=7; yy=(*this)(i,j-1,k); if (yy>xx+delta) Tnsup[n]=2; else if (yy<xx-delta) Tnsup[n]=0;
			l=lmax=0; 
			n=-1; do {n++;} while (n<nneighbors && Tnsup[n%nneighbors]==1); 
			if (n<nneighbors) { 
				zz=Tnsup[n%nneighbors]; 
				do { 
					n++;
//					if (Tnsup[n%nneighbors]==zz) {l++; if (l>lmax) lmax=l; }
					if (Tnsup[n%nneighbors]==zz || Tnsup[n%nneighbors]==1) {l++; if (l>lmax) lmax=l; }
					else if (Tnsup[n%nneighbors]!=1) {l=1; zz=Tnsup[n%nneighbors]; }
				} while (n<2*nneighbors && lmax<=nneighbors);
			}
			ima(i,j)=mini(lmax,nneighbors);
		}
	return(ima);
}

#endif