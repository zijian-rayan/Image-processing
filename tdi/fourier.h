#ifndef _FOURIER_H
#define _FOURIER_H

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
//#include <stdlib.h>
#include <math.h>
#include "constantes.h"
#include "def.h"
#include "matrices.h"
//#include "fonctions.h"
//#include "pixel.h"
//#include "liste_pixels.h"
//#include "imasites.h"
//#include "imadata.h"

class fourier
{
public:
	fourier(void) {}
	~fourier(void) {}

	void fft (matrice1D<double> &x_int_re, matrice1D<double> &x_int_im, unsigned int taille) {
		unsigned int size_2=(taille>>1),tmp1,i;
		double tmp,tmpcos,tmpsin,base=2*PI/taille;
		matrice1D<double> pair_re(size_2), pair_im(size_2), impair_re(size_2), impair_im(size_2);
		for (i=0; i<size_2; i++){
			tmp1=(i<<1);
			pair_re(i)=x_int_re(tmp1);
			pair_im(i)=x_int_im(tmp1);
			impair_re(i)=x_int_re(tmp1+1);
			impair_im(i)=x_int_im(tmp1+1);
		}
		if (taille>2) {
			fft (pair_re,pair_im,size_2);
			fft (impair_re,impair_im,size_2);
		}
		for (i=0; i<size_2; i++) {
			tmp=base*i; tmpcos=cos(tmp); tmpsin=sin(tmp);
			x_int_re(i)=(pair_re(i)+impair_re(i)*tmpcos+impair_im(i)*tmpsin)/SQ_2;
			x_int_im(i)=(pair_im(i)+impair_im(i)*tmpcos-impair_re(i)*tmpsin)/SQ_2;
			x_int_re(i+size_2)=(pair_re(i)-impair_re(i)*tmpcos-impair_im(i)*tmpsin)/SQ_2;
			x_int_im(i+size_2)=(pair_im(i)-impair_im(i)*tmpcos+impair_re(i)*tmpsin)/SQ_2;
		}
	}
	void fft_inv(matrice1D<double> &x_int_re, matrice1D<double> &x_int_im, unsigned int taille) {
		unsigned int size_2=(taille>>1),tmp1,i;	
		double tmp,tmpcos,tmpsin,base=2*PI/taille;
		matrice1D<double> pair_re(size_2), pair_im(size_2), impair_re(size_2), impair_im(size_2);
		for (i=0; i<size_2; i++) {
			tmp1=(i<<1);
			pair_re(i)=x_int_re(tmp1);
			pair_im(i)=x_int_im(tmp1);
			impair_re(i)=x_int_re(tmp1+1);
			impair_im(i)=x_int_im(tmp1+1);
		}
		if (taille>2){
			fft_inv (pair_re,pair_im,size_2);
			fft_inv (impair_re,impair_im,size_2);
		}
		for (i=0; i<size_2; i++){
			tmp=base*i; tmpcos=cos(tmp); tmpsin=sin(tmp);
			x_int_re(i)=(pair_re(i)+impair_re(i)*tmpcos-impair_im(i)*tmpsin)/SQ_2;
			x_int_im(i)=(pair_im(i)+impair_im(i)*tmpcos+impair_re(i)*tmpsin)/SQ_2;
			x_int_re(i+size_2)=(pair_re(i)-impair_re(i)*tmpcos+impair_im(i)*tmpsin)/SQ_2;
			x_int_im(i+size_2)=(pair_im(i)-impair_im(i)*tmpcos-impair_re(i)*tmpsin)/SQ_2;
		}
	}
	void tfd2d (matrice2D<double> &x, matrice2D<double> &xtre, matrice2D<double> &xtim, unsigned int width, unsigned int height) {
		unsigned int i,j,k;
		const double SQ_WIDTH=sqrt((double)width), SQ_HEIGHT=sqrt((double)height);
		double tmp,tmpcos,tmpsin,v_re=0,v_im=0,base1=2*PI/width,base2=2*PI/height;
		matrice2D<double> xint_im(height,width), xint_re(height,width);
		for (i=0; i<height; i++) {
			for (k=0; k<width; k++) {
				for (j=0; j<width; j++) {
					tmp=base1*k*j;
					v_re+=x(i,j)*cos(tmp);
					v_im-=x(i,j)*sin(tmp);
				}
				xint_im(i,k)=v_im/SQ_WIDTH;
				xint_re(i,k)=v_re/SQ_WIDTH;
				v_re=0;
				v_im=0;
			}
		}
		for (j=0; j<width; j++) {
			for (k=0; k<height; k++) {
				for (i=0; i<height; i++) {
					tmp=base2*k*i; tmpcos=cos(tmp); tmpsin=sin(tmp);
					v_re+=xint_re(i,j)*tmpcos+xint_im(i,j)*tmpsin;
					v_im+=xint_im(i,j)*tmpcos-xint_re(i,j)*tmpsin;
				}
				xtre((k+(height>>1))%height,(j+(width>>1))%width)=v_re/SQ_HEIGHT;
				xtim((k+(height>>1))%height,(j+(width>>1))%width)=v_im/SQ_HEIGHT;
				v_re=0;
				v_im=0;
			}
		}
	}
	void tfd2d_inv (matrice2D<double> &xtre, matrice2D<double> &xtim, matrice2D<double> &xrec_re, matrice2D<double> &xrec_im, unsigned int width,unsigned int height) {
		unsigned int  i,j,k ;
		double tmp,tmpcos,tmpsin,u_re=0,u_im=0,base1=2*PI/width,base2=2*PI/height;
		const double SQ_WIDTH=sqrt((double)width), SQ_HEIGHT=sqrt((double)height);
		matrice2D<double> xint_re(height,width), xint_im(height,width);
		for(i=0; i<height; i++) {
			for(k=0; k<width; k++) {
				for(j=0; j<width; j++) {
					tmp=base1*k*j; tmpcos=cos(tmp); tmpsin=sin(tmp);
					u_re+=xtre((i+(height>>1))%height,(j+(width>>1))%width)*tmpcos- xtim((i+(height>>1))%height,(j+(width>>1))%width)*tmpsin;
					u_im+=xtim((i+(height>>1))%height,(j+(width>>1))%width)*tmpcos+ xtre((i+(height>>1))%height,(j+(width>>1))%width)*tmpsin;
				}
				xint_re(i,k)=u_re/SQ_WIDTH;
				xint_im(i,k)=u_im/SQ_WIDTH;
				u_re=0;
				u_im=0;
			}
		}
		for(j=0; j<width; j++) {			
			for(k=0; k<height; k++) {
				for(i=0; i<height; i++) {
					tmp=base2*k*i; tmpcos=cos(tmp); tmpsin=sin(tmp);
					u_re+=xint_re(i,j)*tmpcos-xint_im(i,j)*tmpsin;
					u_im+=xint_im(i,j)*tmpcos+xint_re(i,j)*tmpsin;
				}
				xrec_re(k,j)=u_re/SQ_HEIGHT;
				xrec_im(k,j)=u_im/SQ_HEIGHT;
				u_re=0;
				u_im=0;
			}
		}
	}
	void fft2d(matrice2D<double> &x, matrice2D<double> &xtre, matrice2D<double> &xtim, unsigned int width, unsigned int height) {
		unsigned int i,j,k ;
		matrice2D<double> xint_re(height,width), xint_im(height,width);
		matrice1D<double> x_int_l(width), x_int2_l(width), x_int_c(height), x_int2_c(height);
		for (k=0; k<height; k++) {
			for (j=0; j<width; j++) {
				x_int_l(j)=x(k,j);	
				x_int2_l(j)=0;
			}
			fft(x_int_l,x_int2_l,width);
			for (j=0; j<width; j++) {
				xint_re(k,j)=x_int_l(j);
				xint_im(k,j)=x_int2_l(j);
			}
		}
		for (k=0; k<width; k++) {
			for (i=0; i<height; i++) {
				x_int_c(i)=xint_re(i,k);
				x_int2_c(i)=xint_im(i,k);
			}
			fft(x_int_c,x_int2_c,height) ;
			for (i=0; i<height; i++) {
				xtre((i+(height>>1))%height,(k+(width>>1))%width)=x_int_c(i);
				xtim((i+(height>>1))%height,(k+(width>>1))%width)=x_int2_c(i);
			}
		}
	}
	void fft2d_inv(matrice2D<double> &xtre, matrice2D<double> &xtim, matrice2D<double> &xrec_re, matrice2D<double> &xrec_im, unsigned int width, unsigned int height) {
		unsigned int i,j,k;
		matrice2D<double> xint_re(height,width), xint_im(height,width);
		matrice1D<double> x_int_l(width), x_int2_l(width), x_int_c(height), x_int2_c(height);
		for (k=0; k<height; k++) {
			for (j=0; j<width; j++) {
				x_int_l(j)=xtre((k-(height>>1))%height,(j+(width>>1))%width);	
				x_int2_l(j)=xtim((k-(height>>1))%height,(j+(width>>1))%width);
			}
			fft_inv(x_int_l,x_int2_l,width);
			for (j=0; j<width; j++) {
				xint_re(k,j)=x_int_l(j);
				xint_im(k,j)=x_int2_l(j);
			}
		}
		for (k=0; k<width; k++) {
			for (i=0; i<height; i++) {
				x_int_c(i)=xint_re(i,k);
				x_int2_c(i)=xint_im(i,k);
			}
			fft_inv(x_int_c,x_int2_c,height);
			for (i=0; i<height; i++) {
				xrec_re(i,k)=x_int_c(i);
				xrec_im(i,k)=x_int2_c(i);
			}
		}
	}
};
#endif