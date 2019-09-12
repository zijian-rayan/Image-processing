#ifndef _PIXEL_H
#define _PIXEL_H

template <class T> class pixel
{	int nbcanaux;
	T *adval;
//	void * Tadvois[8];
 public:
	pixel(int n=1, int init=0) {
		adval=NULL;
		adval=new T[nbcanaux=n];
		for (int i=0; i<nbcanaux; i++) adval[i]=(T)init;
		}
	~pixel() {
		if (adval!=NULL) {delete adval; adval=NULL;}
		}
	pixel(const pixel<T> & pix) {
		adval=NULL;
		adval=new T[nbcanaux=pix.nbcanaux];
		for (int i=0; i<nbcanaux; i++) adval[i]=pix.adval[i];
//		for (int i=0; i<8; i++) Tadvois[i]=pix.Tadvois[i];
		}
	pixel<T>& operator=(const pixel<T> & pix) {
		if (this != &pix) {
			if (adval!=NULL) delete adval;
			adval=new T[nbcanaux=pix.nbcanaux];
			for (int i=0; i<nbcanaux; i++) adval[i]=pix.adval[i];
//			for (int i=0; i<8; i++) Tadvois[i]=pix.Tadvois[i];
			}
		return *this;
		}
	T& operator [] (int i) {
		if (i<0 || i>=nbcanaux) {
			cout<<" class pixel : debordement d'indice "<<i;
			cout<<" #canaux "<<nbcanaux<<"\n";
			if (i<0) i=0;
			if (i>=nbcanaux) i=nbcanaux-1;
			}
		return adval[i];
		}
	T operator [] (int i) const {
		if (i<0 || i>=nbcanaux) {
			cout<<" class pixel : debordement d'indice "<<i;
			cout<<" #canaux "<<nbcanaux<<"\n";
			if (i<0) i=0;
			if (i>=nbcanaux) i=nbcanaux-1;
			}
		return adval[i];
		}
	int ncanaux () { return nbcanaux;}
	pixel<T> operator+(const pixel<T> & pix) {
		pixel<T> p(nbcanaux);
		for (int i=0; i<nbcanaux; i++) p[i]=adval[i]+pix.adval[i];
		return p;
		}
	pixel<T> operator-(const pixel<T> & pix) {
		pixel<T> p(nbcanaux);
		for (int i=0; i<nbcanaux; i++) p[i]=adval[i]-pix.adval[i];
		return p;
		}
	pixel<T> operator*(const pixel<T> & pix) {
		pixel<T> p(nbcanaux);
		for (int i=0; i<nbcanaux; i++) p[i]=adval[i]*pix.adval[i];
		return p;
		}
	pixel<T> operator*(float a) {
		pixel<T> p(nbcanaux);
		for (int i=0; i<nbcanaux; i++) p[i]=adval[i]*a;
		return p;
		}
	void sup (const pixel<T> & pix) {
		for (int i=0; i<nbcanaux; i++) 
			if (pix.adval[i]>adval[i]) adval[i]=pix.adval[i];
		}
	void inf (const pixel<T> & pix) {
		for (int i=0; i<nbcanaux; i++) 
			if (pix.adval[i]<adval[i]) adval[i]=pix.adval[i];
		}
	void affiche (bool=1);
	float distance (pixel<T> &);
};

template <class T> void pixel<T>::affiche(bool iback) {
	cout<<" (";
	for (int i=0; i<nbcanaux-1; i++) cout<<setw(3)<<adval[i]<<",";
	cout<<setw(3)<<adval[nbcanaux-1]<<")";
	if (iback) cout<<"\n";
	}

template <class T> float pixel<T>::distance(pixel<T> & pix) {
	float dist=0, d;
	for (int i=0; i<nbcanaux; i++) {
		d=adval[i]-pix.adval[i];
//		cout<<i<<" "<<adval[i]<<" "<<pix.adval[i]<<" "<<d<<" "<<dist<<"\n";
		dist+=d*d;
		}
	return dist;
	}

#endif