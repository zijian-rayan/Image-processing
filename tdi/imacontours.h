#ifndef _IMACONTOURS_H
#define _IMACONTOURS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;
#include <math.h>
#include "def.h"
#include "constantes.h"
#include "fonctions.h"
#include "imasites.h"
#include "pixel.h"
#include "liste_pixels.h"
#include "imadata.h"
#include "imabin.h"

class keypoints
{
	int n_key_points_max, n_key_points, lg_descr;
	float **Tab_key_points, **Tdescr_key_points;

	bool approx_Taylor (imadata<float> &imaD, int nblig, int nbcol, int nb_ech, int i, int j, int k, const float th_discard, matrice2D<double> &x_hat, float &val_x_corr) {
		const double epsilon=1.e-12;
		int i0,i2,j0,j2,k0,k2;
		matrice2D<double> MatD2(3,3), MatD2_inv(3,3), vect_D(3,1), x_hat2(3,1);
		double det_MatD2;
		i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1); k0=maxi(k-1,0); k2=mini(j+1,nb_ech-1);
		MatD2(0,0)=(imaD(i,j2,k)+imaD(i,j0,k)-2*imaD(i,j,k))/4.;
		MatD2(1,1)=(imaD(i2,j,k)+imaD(i0,j,k)-2*imaD(i,j,k))/4.;
		MatD2(0,1)=MatD2(1,0)=(imaD(i2,j2,k)+imaD(i0,j0,k)-imaD(i2,j0,k)-imaD(i0,j2,k))/4.;
		MatD2(0,2)=MatD2(2,0)=(imaD(i,j2,k2)+imaD(i,j0,k0)-imaD(i,j0,k2)-imaD(i,j2,k0))/4.;
		MatD2(1,2)=MatD2(2,1)=(imaD(i2,j,k2)+imaD(i0,j,k0)-imaD(i0,j,k2)-imaD(i2,j,k0))/4.;
		MatD2(2,2)=(imaD(i,j,k2)+imaD(i,j,k0)-2*imaD(i,j,k))/4.;
		det_MatD2=MatD2.determinant3x3();
		if (fabs(det_MatD2)<epsilon) {
			//cout<<" Pb determinant nul ("<<det_MatD2<<") pour le point ("<<i<<","<<j<<","<<k<<") !!!\n"; 
			return false;}
		else {
			vect_D(0,0)=(imaD(i,j2,k)-imaD(i,j0,k))/2.; vect_D(1,0)=(imaD(i2,j,k)-imaD(i0,j,k))/2.; vect_D(2,0)=(imaD(i,j,k2)-imaD(i,j,k0))/2.;
			matrice2D<double> Matbid(MatD2);
			Matbid(0,0)=vect_D(0,0); Matbid(1,0)=vect_D(1,0); Matbid(2,0)=vect_D(2,0); 
			x_hat(0,0)=-Matbid.determinant3x3()/det_MatD2;
			Matbid=MatD2; Matbid(0,1)=vect_D(0,0); Matbid(1,1)=vect_D(1,0); Matbid(2,1)=vect_D(2,0); 
			x_hat(1,0)=-Matbid.determinant3x3()/det_MatD2;
			Matbid=MatD2; Matbid(0,2)=vect_D(0,0); Matbid(1,2)=vect_D(1,0); Matbid(2,2)=vect_D(2,0); 
			x_hat(2,0)=-Matbid.determinant3x3()/det_MatD2;
			if (fabs(x_hat(0,0))>0.5 || fabs(x_hat(1,0))>0.5 || fabs(x_hat(2,0))>0.5) {
//				cout<<" point cle non stable en ("<<i<<","<<j<<","<<k<<"):\n"; 
				if (fabs(x_hat(0,0))>th_discard || fabs(x_hat(1,0))>th_discard) return false;
				else val_x_corr=(float)(fabs(imaD(i,j,k)+0.5*(matrice2D<double>(vect_D.transpo(),x_hat,1,3,1)(0,0)))); return true;
			}
			else {//cout<<" point cle stable\n"; 
				val_x_corr=(float)(fabs(imaD(i,j,k)+0.5*(matrice2D<double>(vect_D.transpo(),x_hat,1,3,1)(0,0)))); return true;
//			cout<<" -> x-\hat{x} = ("<<j-x_hat(0,0)<<","<<i-x_hat(1,0)<<","<<k-x_hat(2,0)<<")\n";
//			cout<<j<<" "<<i<<" "<<k<<" -> \hat{x} = ("<<x_hat(0,0)<<" "<<x_hat(1,0)<<" "<<x_hat(2,0)<<")\n"; char aa; cin>>aa; 
			}
		}
/*		if (!MatD2.inverse(MatD2_inv)) {cout<<" Pb dans inversion matrice pour le point ("<<i<<","<<j<<","<<k<<") !!! determinant = "<<det_MatD2<<"\n"; ima_keypoints(i,j)=0; }
		else {
//cout<<" Inversion matrice OK"; //(matrice2D<double>(MatD2,MatD2_inv,3,3,3)).affiche(); (matrice2D<double>(MatD2_inv,MatD2,3,3,3)).affiche();
			vect_D(0,0)=(imaD(i,j2,k)-imaD(i,j0,k))/2.; vect_D(1,0)=(imaD(i2,j,k)-imaD(i0,j,k))/2.; vect_D(2,0)=(imaD(i,j,k2)-imaD(i,j,k0))/2.;
			x_hat2=matrice2D<double>(MatD2_inv,vect_D,3,3,1)*(-1.);
			xx=fabs(imaD(i,j,k)+0.5*(matrice2D<double>(vect_D.transpo(),x_hat2,1,3,1)(0,0)));
			imaDxhat(i,j)=xx;
			if (xx<=0.1*255) ima_keypoints(i,j)=0; // dans Lowe, 0.03 pour image entre 0 et 1
			if (fabs(x_hat2(0,0))>0.5 || fabs(x_hat2(1,0))>0.5 || fabs(x_hat2(2,0))>0.5) //{
			cout<<" point cle non stable en ("<<i<<","<<j<<","<<k<<"):\n"; //MatD2_inv.affiche(); vect_D.transpo().affiche(); x_hat2.transpo().affiche(); char aa; cin>>aa;}
			else cout<<" point cle stable\n";
//			cout<<" -> x-\hat{x} = ("<<j-x_hat2(0,0)<<","<<i-x_hat2(1,0)<<","<<k-x_hat2(2,0)<<")\n";
			cout<<" *** "<<j<<" "<<i<<" "<<k<<" "<<x_hat2(0,0)<<" "<<x_hat2(1,0)<<" "<<x_hat2(2,0)<<" "<<"\n";
		}*/
	}

public:
	keypoints (imadata<float> &, float, float=2000, bool=1, bool=0/*, char *detecteur*/); // détecteur de Harris (éventuellement à compléter -> Harris-Laplace etc

	keypoints (imadata<float> &ima, const float sigma, const float param_k, const int nb_ech, const float seuilDoG=0.03f, bool iaffich=0) {
		const bool i_interpol=0, i_Harris=1, i_Hessian=0;
		int nblig=ima.nlig(), nbcol=ima.ncol(), i,j,k;
//		nb_ech=log(pow((double)nblig*nblig+(double)nbcol*nbcol,0.5)/sigma);
// Construction de l'image des DoG (Difference of Gaussians)
		imadata<float> imaL(nblig,nbcol,nb_ech+1);
		imaL.copiecanal(0,ima);
		for (k=1; k<nb_ech+1; k++)
			imaL.copiecanal(k,ima.filtregaussienne(pow(param_k,(float)k)*sigma));
		imaL.sauve_ImaBSQ("imaL.dat"); // avec les paramètres de Lowe: sigmaL=1.6 (k=0), sigmaL=2.01 (k=1), sigmaL=2.54 (k=2), sigmaL=3.2 (k=3), sigmaL=4.03 (k=4), sigmaL=5.08 (k=5)
		imadata<float> imaD(nblig,nbcol,nb_ech);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) 
				for (k=0; k<nb_ech; k++) imaD(i,j,k)=imaL(i,j,k+1)-imaL(i,j,k);
		imaD.statbasic(1);
		imadata<float> imaDbis(imaD); //		imaD.normalize_all(); 
		float fact_norm;
		for (k=0; k<nb_ech; k++) {
			fact_norm=(float)maxi(fabs(imaD.maxI(k)),fabs(imaD.minI(k)));
			if (fact_norm>0) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) imaD(i,j,k)/=fact_norm;
			}
		}
		imaD.sauve_ImaBSQ("imaDoG.dat");
		float Lx, Ly, xx;
// Detection des extremas dans l'espace des échelles
		imadata<BYTE> ima_keypoints(nblig,nbcol);
		ima_keypoints.mise_a_zero();
		bool extremI, extremS, fini;
		int lg_f,i0,i2,ii,j0,j2,jj,k0,k2,kk,nb_keypoints=0;
		for (k=0; k<nb_ech; k++) {
			k0=maxi(0,k-1); k2=mini(nb_ech-1,k+1);
			lg_f=(int)floor(pow(param_k,(float)k)*sigma);
			for (i=0; i<nblig; i++) {
				i0=maxi(0,i-lg_f); i2=mini(nblig-1,i+lg_f); 
				for (j=0; j<nbcol; j++) {
					j0=maxi(0,j-lg_f); j2=mini(nbcol-1,j+lg_f); 
					xx=imaD(i,j,k);
					extremI=extremS=1; fini=0;
					ii=i0; jj=j0; kk=k0;
					while ((extremI || extremS) && !fini) {
						if (imaD(ii,jj,kk)>xx && (ii!=i || jj!=j || kk!=k)) extremI=0; // > ou >= ??? sur forum SIFT on garde point même si extrema au sens large
						if (imaD(ii,jj,kk)<xx && (ii!=i || jj!=j || kk!=k)) extremS=0;
						kk++;
						if (kk>k2) {jj++; kk=k0;}
						if (jj>j2) {ii++; jj=j0;}
						if (ii>i2) fini=1;
					}
					if ((extremI || extremS) && fabs(xx)>seuilDoG) {
						ima_keypoints(i,j)=k+1; nb_keypoints++;
//						if (iaffich) cout<<" point-cle trouve en ("<<i<<","<<j<<","<<k<<") avec DoG = "<<xx<<"\n";
					}
				}
			}
		}
		if (nb_keypoints>0.01*nblig*nbcol*nb_ech) {
			cout<<" avec seuilDoG par defaut : "<<nb_keypoints<<" points-cle -> trop nombreux -> modification du seuil pour rester en dessous de 1% de points cles\n";
			imadata<float> imaD2(imaD.absI()); 		imaD2.sauve_ImaBSQ("imaDoG_2.dat");
			const int nbinH_D2=256;
			float seuilDoG_2=imaD2.histogramme_allc(nbinH_D2,1.f-0.01f); cout<<" nouveau seuil = "<<seuilDoG_2<<"\n";
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if (ima_keypoints(i,j)>0 && fabs(imaD(i,j,ima_keypoints(i,j)-1))<seuilDoG_2) ima_keypoints(i,j)=0;
		}
		(ima_keypoints*(255.f/nb_ech)).sauve_ImaPGM("ima_keypoints_all.pgm");
// Amelioration de la précision par interpolation des coordonnées
		if (i_interpol) {
//			float Dx, Dy, Dxx, Dyy, Dxy, Dxs, Dys, Dss;
			matrice2D<double> /*MatD2(3,3), MatD2_inv(3,3), vect_D(3,1),*/ x_hat(3,1)/*, x_hat2(3,1)*/;
//			double det_MatD2; 
			imadata<float> imaDxhat(nblig,nbcol); imaDxhat.mise_a_zero(); 
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) 
					if (ima_keypoints(i,j)>0) {
						k=ima_keypoints(i,j)-1;
						if (!approx_Taylor(imaD,nblig,nbcol,nb_ech,i,j,k,pow(param_k,(float)k)*sigma,x_hat,xx)) ima_keypoints(i,j)=0;
						else {
							imaDxhat(i,j)=xx;
							if (xx<0.03) ima_keypoints(i,j)=0; // dans Lowe, 0.03 pour image entre 0 et 1
						}
					}
			imaDxhat.sauve_Ima("imaDxhat.dat");
			(ima_keypoints*(255.f/nb_ech)).sauve_ImaPGM("ima_keypoints_interp.pgm");
		} //else {}
// Elimination des points-clés sur les contours
		if (i_Harris) { cout<<" Elimination des points-cles sur les contours selon critere de Harris\n";
			const float kappa=0.04f/*0.15f*/; 
			imadata<float> imaIxIy(nblig,nbcol,3*nb_ech), imaMc(nblig,nbcol,nb_ech); 
			imaIxIy.mise_a_zero(); imaMc.mise_a_zero();
			for (k=1; k<nb_ech+1; k++) {		// Detection de Harris nécessite préfiltrage de l'image à traiter
				k0=3*(k-1);
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						if (j<nbcol-1 && j>0) Lx=fabs(imaL(i,j+1,k)-imaL(i,j-1,k))/2; 
						else {if (j>0) Lx=fabs(imaL(i,j,k)-imaL(i,j-1,k))/2; else Lx=fabs(imaL(i,j+1,k)-imaL(i,j,k))/2;}
						if (i<nblig-1 && i>0) Ly=fabs(imaL(i+1,j,k)-imaL(i-1,j,k))/2; 
						else {if (i>0) Ly=fabs(imaL(i,j,k)-imaL(i-1,j,k))/2; else Ly=fabs(imaL(i+1,j,k)-imaL(i,j,k))/2;}
						imaIxIy(i,j,0+k0)=Lx*Lx; imaIxIy(i,j,1+k0)=Ly*Ly; imaIxIy(i,j,2+k0)=Lx*Ly;
					} 
			} imaIxIy.sauve_ImaBSQ("imaIxIy.dat");
/*			imadata<float> ima_dirgrdL(nblig,nbcol,nb_ech+1), ima_normgrdL=imaL.gradient(ima_dirgrdL);
			ima_normgrdL.sauve_ImaBSQ("ima_normgrdL.dat"); (ima_dirgrdL*(180./PI)).sauve_ImaBSQ("ima_dirgrdL.dat");
			for (k=1; k<nb_ech+1; k++) {
				k0=3*(k-1);
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						Lx=fabs(ima_normgrdL(i,j,k)*cos(ima_dirgrdL(i,j,k))); Ly=fabs(ima_normgrdL(i,j,k)*sin(ima_dirgrdL(i,j,k)));
						imaIxIy(i,j,0+k0)=Lx*Lx; imaIxIy(i,j,1+k0)=Ly*Ly; imaIxIy(i,j,2+k0)=Lx*Ly;}
			}*/
			imaIxIy=imaIxIy.filtregaussienne(sigma); // OBLIGATOIRE sinon le déterminant de la matrice est nul !!!
			imaIxIy.sauve_ImaBSQ("imaIxIy_filtree.dat"); 
//			matrice2D<double> M_Harris(2,2); double T_vp[2];
			for (k=0; k<nb_ech; k++) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
//						M_Harris(0,0)=imaIxIy(i,j,0+3*k); M_Harris(1,1)=imaIxIy(i,j,1+3*k); M_Harris(0,1)=M_Harris(1,0)=imaIxIy(i,j,2+3*k);
//						if (M_Harris.valeurspropres2x2(T_vp)) imaMc(i,j,k)=mini(T_vp[0],T_vp[1]);
//						imaMc(i,j,k)=2*(imaIxIy(i,j,0+3*k)*imaIxIy(i,j,1+3*k)-pow(imaIxIy(i,j,2+3*k),2))/(imaIxIy(i,j,0+3*k)+imaIxIy(i,j,1+3*k)+1.e-6);
						imaMc(i,j,k)=imaIxIy(i,j,0+3*k)*imaIxIy(i,j,1+3*k)-pow(imaIxIy(i,j,2+3*k),2)-kappa*pow(imaIxIy(i,j,0+3*k)+imaIxIy(i,j,1+3*k),2);
/*						xx=(imaIxIy(i,j,0+3*k)*imaIxIy(i,j,1+3*k)-pow(imaIxIy(i,j,2+3*k),2));
						if (xx>0) imaMc(i,j,k)=pow(imaIxIy(i,j,0+3*k)+imaIxIy(i,j,1+3*k),2)/(imaIxIy(i,j,0+3*k)*imaIxIy(i,j,1+3*k)-pow(imaIxIy(i,j,2+3*k),2));
						else {
							if (imaIxIy(i,j,0+3*k)+imaIxIy(i,j,1+3*k)==0) imaMc(i,j,k)=0;
							else imaMc(i,j,k)=-1;
						}*/
					}
			}
			imaMc.sauve_ImaBSQ("imaMc.dat"); //imaMc.histogramme(256,1);
//			imadata<float> imaMc_bis(imaMc); 
//			for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) imaMc_bis(i,j)=mini(maxi(-1.f,imaMc_bis(i,j)),2500.f);
//			imaMc_bis.sauve_Ima("imaMc_bis.dat");
			double prec, seuil_Mc=imaMc.modeI(prec); cout<<" seuil sur Mc = "<<seuil_Mc<<"\n";
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					if (ima_keypoints(i,j)>0) {
						k=ima_keypoints(i,j)-1;
//						if (iaffich) cout<<" point-cle en ("<<i<<","<<j<<","<<k<<") s.t. DoG = "<<imaD(i,j,ima_keypoints(i,j)-1)<<" -> critere Harris Mc = "<<imaMc(i,j)<<"\n";
//						if (imaMc(i,j,k)<=seuil_Mc) ima_keypoints(i,j)=0; 
						xx=imaMc(i,j,k); int i_Mcmax=i, j_Mcmax=j;
						lg_f=(int)floor(pow(param_k,(float)k)*sigma);
						i0=maxi(i-lg_f,0); i2=mini(i+lg_f,nblig-1); j0=maxi(j-lg_f,0); j2=mini(j+lg_f,nbcol-1); 
						for (ii=i0; ii<=i2; ii++)
							for (jj=j0; jj<=j2; jj++) 
								if (imaMc(ii,jj,k)>xx) {xx=imaMc(ii,jj,k); i_Mcmax=ii; j_Mcmax=jj;}
//						if ((i==37&&j==83)||(i==37&&j==89)||(i==77&&j==127)||(i==37&&j==167)) cout<<k<<" "<<lg_f<<" "<<i<<" "<<j<<" max Mc sur vois = "<<xx<<"\n";
						if (xx<=seuil_Mc) ima_keypoints(i,j)=0;
						else
							if (i!=i_Mcmax || j!=j_Mcmax) {
//								cout<<" point-cle en ("<<i<<","<<j<<","<<k<<") peu precis en localisation : d'apres critere Harris plutot en ("<<i_Mcmax<<","<<j_Mcmax<<")\n";
								ima_keypoints(i_Mcmax,j_Mcmax)=k+1; ima_keypoints(i,j)=0;
							}
// si on veut faire test de maximum local (cf. lignes commentées juste en dessous), il faut être précis en localisation
/*						extremI=1; fini=0;
						i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1); ii=i0; jj=j0; 
						while (extremI && !fini) {
							if (imaMc(ii,jj,k)>=xx && (ii!=i || jj!=j)) extremI=0;
							jj++;	if (jj>j2) {ii++; jj=j0;}	if (ii>i2) fini=1;
						}
//						if (!extremI || imaMc(i,j)<=1500) ima_keypoints(i,j)=0;
						if (!extremI || imaMc(i,j,k)<=seuil_Mc) ima_keypoints(i,j)=0;*/
					}
					if (ima_keypoints(i,j)>0 && iaffich) cout<<" point-cle en ("<<i<<","<<j<<","<<ima_keypoints(i,j)-1<<") s.t. DoG = "<<imaD(i,j,ima_keypoints(i,j)-1)<<" -> critere Harris Mc = "<<imaMc(i,j)<<"\n";
				}
		}
		if (i_Hessian) { cout<<" Elimination des points-cles sur les contours selon critere base sur la matrice hessienne\n";
			const double r_th_H=1.; // 10. dans article de Lowe
			double Hxx, Hyy, Hxy, Tr_H, Det_H;
			imadata<float> ima_rH(nblig,nbcol); 
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) 
					if (ima_keypoints(i,j)>=0) {
						k=maxi(ima_keypoints(i,j)-1,0);
						i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
						Hxx=(imaD(i,j2,k)+imaD(i,j0,k)-2*imaD(i,j,k))/4.; Hyy=(imaD(i2,j,k)+imaD(i0,j,k)-2*imaD(i,j,k))/4.; 
						Hxy=(imaD(i2,j2,k)+imaD(i0,j0,k)-imaD(i2,j0,k)-imaD(i0,j2,k))/4.; 
						Tr_H=Hxx+Hyy; Det_H=Hxx*Hyy-(Hxy*Hxy);
						if (Tr_H*Tr_H/Det_H>=pow(r_th_H+1,2)/r_th_H) ima_keypoints(i,j)=0;
						ima_rH(i,j)=(float)mini(maxi(-128.,Tr_H*Tr_H/Det_H),128.);
					}
			ima_rH.sauve_Ima("ima_rH.dat");
		}
		(ima_keypoints*(255.f/nb_ech)).sauve_ImaPGM("ima_keypoints_Harris.pgm");
// Estimation d'orientation et sauvegarde des key-points dans un tableau
// axe y même sens et orientation que l'axe des lignes (vers le bas) et axe x même sens et orientation que l'axe des colonnes (vers la gauche)
		int n; n_key_points=0;
		for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) if (ima_keypoints(i,j)>0) n_key_points++;
		n_key_points_max=3*n_key_points;
		Tab_key_points=new float*[n_key_points_max];
		for (n=0; n<n_key_points_max; n++) Tab_key_points[n]=new float[4];
		n_key_points=0;
		const int nb_bin_Horient=36;
		const float fact_HO=(float)(nb_bin_Horient/(2.*PI));
		float Horient[nb_bin_Horient], max0, max1, max2;
		int theta, theta0, theta1, theta2;
		imadata<float> ima_dirgrdL(nblig,nbcol,nb_ech+1), ima_normgrdL=imaL.gradient(ima_dirgrdL);
		ima_normgrdL.sauve_ImaBSQ("ima_normgrdL.dat"); (ima_dirgrdL*(float)(180./PI)).sauve_ImaBSQ("ima_dirgrdL.dat");

		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (ima_keypoints(i,j)>0) { 
					k=ima_keypoints(i,j)-1;
					for (theta=0; theta<nb_bin_Horient; theta++) Horient[theta]=0; 
					i0=maxi(0,i-7); i2=mini(nblig-1,i+8); j0=maxi(0,j-7); j2=mini(nbcol-1,j+8);
					for (ii=i0; ii<=i2; ii++) 
						for (jj=j0; jj<=j2; jj++) {
//							if (jj<nbcol-1 && jj>0) Lx=(imaL(ii,jj+1,k)-imaL(ii,jj-1,k)); 
//							else {if (jj>0) Lx=(imaL(ii,jj,k)-imaL(ii,jj-1,k)); else Lx=(imaL(ii,jj+1,k)-imaL(ii,jj,k));}
//							if (ii<nblig-1 && ii>0) Ly=(imaL(ii+1,jj,k)-imaL(ii-1,jj,k)); 
//							else {if (ii>0) Ly=(imaL(ii,jj,k)-imaL(ii-1,jj,k)); else Ly=(imaL(ii+1,jj,k)-imaL(ii,jj,k));}
							Lx=ima_normgrdL(ii,jj,k)*cos(ima_dirgrdL(ii,jj,k)); Ly=ima_normgrdL(ii,jj,k)*sin(ima_dirgrdL(ii,jj,k));
							xx=(float)pow(Lx*Lx+Ly*Ly,0.5f)*exp(-0.5f*((ii-i)*(ii-i)+(jj-j)*(jj-j))/pow(1.5f*(k+1)*sigma,2));
							theta=(int)floor(atan2(Ly,Lx)*fact_HO); 
							while (theta<0) theta+=nb_bin_Horient;
//							if (theta<0) theta+=(int)(2*PI*fact_HO); theta=theta%nb_bin_Horient;
							Horient[theta]+=xx;
						}
/*					if ((i>=105&&i<=115)&&(j>=14&&j<=24)) {
						cout<<" coin triangle U image 2 "<<i<<" "<<j<<" \n"; 
						for (theta=0; theta<nb_bin_Horient; theta++) cout<<" "<<fixed<<setw(4)<<setprecision(0)<<Horient[theta]; cout<<"\n";
						for (ii=i0; ii<=i2; ii++) {
							for (jj=j0; jj<=j2; jj++) {
//								if (jj<nbcol-1 && jj>0) Lx=(imaL(ii,jj+1,k)-imaL(ii,jj-1,k)); 
//								else {if (jj>0) Lx=(imaL(ii,jj,k)-imaL(ii,jj-1,k)); else Lx=(imaL(ii,jj+1,k)-imaL(ii,jj,k));}
//								if (ii<nblig-1 && ii>0) Ly=(imaL(ii+1,jj,k)-imaL(ii-1,jj,k)); 
//								else {if (ii>0) Ly=(imaL(ii,jj,k)-imaL(ii-1,jj,k)); else Ly=(imaL(ii+1,jj,k)-imaL(ii,jj,k));}
								Lx=ima_normgrdL(ii,jj,k)*cos(ima_dirgrdL(ii,jj,k)); Ly=ima_normgrdL(ii,jj,k)*sin(ima_dirgrdL(ii,jj,k));
								cout<<" ("<<fixed<<setw(4)<<setprecision(0)<<Lx<<","<<fixed<<setw(4)<<setprecision(0)<<Ly<<")";
							}
							cout<<"\n";
						}
					}
					if ((i>= 33&&i<= 43)&&(j>=60&&j<=70)) {
						cout<<" coin triangle U image 1 "<<i<<" "<<j<<" \n"; 
						for (theta=0; theta<nb_bin_Horient; theta++) cout<<" "<<fixed<<setw(4)<<setprecision(0)<<Horient[theta]; cout<<"\n";
						for (ii=i0; ii<=i2; ii++) {
							for (jj=j0; jj<=j2; jj++) {
//								if (jj<nbcol-1 && jj>0) Lx=(imaL(ii,jj+1,k)-imaL(ii,jj-1,k)); 
//								else {if (jj>0) Lx=(imaL(ii,jj,k)-imaL(ii,jj-1,k)); else Lx=(imaL(ii,jj+1,k)-imaL(ii,jj,k));}
//								if (ii<nblig-1 && ii>0) Ly=(imaL(ii+1,jj,k)-imaL(ii-1,jj,k)); 
//								else {if (ii>0) Ly=(imaL(ii,jj,k)-imaL(ii-1,jj,k)); else Ly=(imaL(ii+1,jj,k)-imaL(ii,jj,k));}
								Lx=ima_normgrdL(ii,jj,k)*cos(ima_dirgrdL(ii,jj,k)); Ly=ima_normgrdL(ii,jj,k)*sin(ima_dirgrdL(ii,jj,k));
								cout<<" ("<<fixed<<setw(4)<<setprecision(0)<<Lx<<","<<fixed<<setw(4)<<setprecision(0)<<Ly<<")";
							}
							cout<<"\n";
						}
					}*/
					max0=max1=max2=0.f; theta0=theta1=theta2=-1;
					for (theta=0; theta<nb_bin_Horient; theta++) if (Horient[theta]>max0) {max0=Horient[theta]; theta0=theta;}					
					for (theta=0; theta<nb_bin_Horient; theta++) if (theta!=theta0 && Horient[theta]>max1) {max1=Horient[theta]; theta1=theta;}
					if (max1<=0.8*max0) theta1=-1; // pic secondaires à au moins 80% du pic principal d'après article de Lowe
					else {
						for (theta=0; theta<nb_bin_Horient; theta++) if (theta!=theta0 && theta!=theta1 && Horient[theta]>max2) {max2=Horient[theta]; theta2=theta;}
						if (max2<=0.8*max0) theta2=-1;
					}
//					cout<<" point-cle en ("<<i<<","<<j<<","<<ima_keypoints(i,j)-1<<") orientation(s) = "<<theta0*360.f/nb_bin_Horient;
					Tab_key_points[n_key_points][0]=(float)i; Tab_key_points[n_key_points][1]=(float)j; 
					Tab_key_points[n_key_points][2]=(float)(ima_keypoints(i,j)-1); 
					Tab_key_points[n_key_points++][3]=(float)theta0; 
					if (n_key_points>=n_key_points_max) cout<<" Pb dans dimensionnement du tableau Tab_key_points "<<n_key_points<<" >= "<<n_key_points_max<<"\n";
					if (theta1>=0) {//cout<<" et "<<theta1*360.f/nb_bin_Horient; //}
						for (n=0; n<3; n++) Tab_key_points[n_key_points][n]=Tab_key_points[n_key_points-1][n]; Tab_key_points[n_key_points++][3]=(float)theta1;}
					if (n_key_points>=n_key_points_max) cout<<" Pb dans dimensionnement du tableau Tab_key_points "<<n_key_points<<" >= "<<n_key_points_max<<"\n";
					if (theta2>=0) {//cout<<" et "<<theta2*360.f/nb_bin_Horient; //}
						for (n=0; n<3; n++) Tab_key_points[n_key_points][n]=Tab_key_points[n_key_points-1][n]; Tab_key_points[n_key_points++][3]=(float)theta2;}
					if (n_key_points>=n_key_points_max) cout<<" Pb dans dimensionnement du tableau Tab_key_points "<<n_key_points<<" >= "<<n_key_points_max<<"\n";
//					cout<<"\n";
				}
// Tracé des résultats dans une image
		float xtheta;
		const int l_o=5;
		for (k=nb_ech; k>0; k--) {
			for (n=0; n<n_key_points; n++)
				if (Tab_key_points[n][2]+1==k) {
					i=(int)Tab_key_points[n][0]; j=(int)Tab_key_points[n][1]; k2=k*k;
					i0=maxi(0,i-k); i2=mini(nblig-1,i+k); j0=maxi(0,j-k); j2=mini(nbcol-1,j+k);
					for (ii=i0; ii<=i2; ii++) for (jj=j0; jj<=j2; jj++) if ((ii-i)*(ii-i)+(jj-j)*(jj-j)<=k2) ima_keypoints(ii,jj)=k;
					xtheta=((Tab_key_points[n][3]+0.5f)*360.f/nb_bin_Horient); // xtheta encore en degres
					if (xtheta<=90 || xtheta>=270) {j0=j; j2=mini(nbcol-1,j+l_o*k);} else {j0=maxi(0,j-l_o*k); j2=j;}
					if (xtheta<=180) {i0=i; i2=mini(nblig-1,i+l_o*k);} else {i0=maxi(0,i-l_o*k); i2=i;} 
					xtheta=xtheta/(float)(180.*PI); // xtheta encore en radians
					if (fabs(xtheta-PI/2)>=0.01 && fabs(xtheta-3*PI/2)>=0.01) 
						for (jj=j0; jj<=j2; jj++) {ii=(int)((jj-j)*tan(xtheta))+i; if (ii>=i0 && ii<=i2) ima_keypoints(mini(maxi(ii,0),nblig-1),jj)=k;} 
					if (fabs(xtheta-0)>=0.01 && fabs(xtheta-PI)>=0.01) 
						for (ii=i0; ii<=i2; ii++) {jj=(int)((ii-i)*cos(xtheta)/sin(xtheta))+j; if (jj>=j0 && jj<=j2) ima_keypoints(ii,mini(maxi(jj,0),nbcol-1))=k;}
				}
		}
		(ima_keypoints*(255.f/nb_ech)).sauve_ImaPGM("ima_keypoints_orient.pgm");
// Calcul du descripteur en chaque point-clé
		const int sz_calc_descr=16, nb_bin_Hdesc=8;
		const float fact_HD=(float)(nb_bin_Hdesc/(2.*PI));
		Tdescr_key_points=new float*[n_key_points];
		lg_descr=sz_calc_descr/4*sz_calc_descr/4*nb_bin_Hdesc;
		for (n=0; n<n_key_points; n++) {
			Tdescr_key_points[n]=new float[lg_descr];
			for (theta=0; theta<lg_descr; theta++) Tdescr_key_points[n][theta]=0;
		}
		imadata<float> imagette(2*sz_calc_descr+1,2*sz_calc_descr+1);
		int i_H, j_H;
		float d_bin;
		for (n=0; n<n_key_points; n++) { //cout<<" calcul du descripteur du point-cle "<<n<<"\n";
			i=(int)Tab_key_points[n][0]; j=(int)Tab_key_points[n][1]; k=(int)Tab_key_points[n][2]; xtheta=((Tab_key_points[n][3]+0.5f)*360.f/nb_bin_Horient);
			imagette.mise_a_zero(); 
			i0=maxi(0,i-sz_calc_descr); i2=mini(nblig-1,i+sz_calc_descr); j0=maxi(0,j-sz_calc_descr); j2=mini(nbcol-1,j+sz_calc_descr);
			for (ii=i0; ii<=i2; ii++) for (jj=j0; jj<=j2; jj++) imagette(ii-i+sz_calc_descr,jj-j+sz_calc_descr)=imaL(ii,jj,k);
// soit on fait rotation de l'imagette, soit on décale de -xtheta l'histogramme des orientations du vecteur caractéristique
			imadata<float> imagette_rot;
			const bool irotat_effective=1;
			if (irotat_effective) imagette_rot=imagette.projette_ima(0,0,1,-xtheta);
			else imagette_rot=imagette;
//			if ((i>=72&&i<=84)&&(j>=121&&j<=135)) {cout<<" coin carre LL image 2 "<<i<<" "<<j<<" rotation "<<-xtheta<<"\n"; imagette_rot.imaunsignedchar(1).sauve_ImaPGM("imagette2_rot.pgm"); /*char aa; cin>>aa;*/}
//			if ((i>=75&&i<=96)&&(j>=155&&j<=182)) {cout<<" coin carre LL image 1 "<<i<<" "<<j<<" rotation "<<-xtheta<<"\n"; imagette_rot.imaunsignedchar(1).sauve_ImaPGM("imagette1_rot.pgm"); /*char aa; cin>>aa;*/}
//			if ((i>=105&&i<=115)&&(j>=14&&j<=24)) {cout<<" coin triangle U image 2 "<<i<<" "<<j<<" rotation "<<-xtheta<<"\n"; imagette_rot.imaunsignedchar(1).sauve_ImaPGM("imagette2_rot.pgm"); /*char aa; cin>>aa;*/}
//			if ((i>= 33&&i<= 43)&&(j>=60&&j<=70)) {cout<<" coin triangle U image 1 "<<i<<" "<<j<<" rotation "<<-xtheta<<"\n"; imagette_rot.imaunsignedchar(1).sauve_ImaPGM("imagette1_rot.pgm"); /*char aa; cin>>aa;*/}
			i0=j0=sz_calc_descr/2; i2=j2=3*sz_calc_descr/2;
			for (ii=i0; ii<i2; ii++) {
				i_H=(ii-i0)/(sz_calc_descr/4)*(sz_calc_descr/4);	// i_H=(ii-i0)/(sz_calc_descr/4) mais tjs utilisé *sz_calc_descr/4
				for (jj=j0; jj<j2; jj++) {
					j_H=(jj-j0)/(sz_calc_descr/4);									// j_H=(jj-j0)/(sz_calc_descr/4) ET PAS tjs utilisé *sz_calc_descr/4
					Lx=imagette_rot(ii,jj+1)-imagette_rot(ii,jj-1); Ly=imagette_rot(ii+1,jj)-imagette_rot(ii-1,jj); 
					xx=pow(Lx*Lx+Ly*Ly,0.5f);
//					xx*=exp(-0.5*(pow((float)ii-sz_calc_descr,2)+pow((float)jj-sz_calc_descr,2))/pow(1.5*(k+1)*sigma,2)); // article Lowe pondération gaussienne centrée sur le point-clé
//				  xx*=exp(-0.5*(pow((float)ii-(i_H+sz_calc_descr/8),2)+pow((float)jj-(j_H*(sz_calc_descr/4)+sz_calc_descr/8),2))/pow(1.5*(k+1)*sigma,2));
					theta=(int)floor((atan2(Ly,Lx)+(1.f-irotat_effective)*(-xtheta/180.*PI))*fact_HD);
					d_bin=(float)((atan2(Ly,Lx)+(1.f-irotat_effective)*(-xtheta/180.*PI))*fact_HD-theta);
/*					if ((i>=33&&i<=43)&&(j>=60&&j<=70) || (i>=105&&i<=115)&&(j>=14&&j<=24)) {
						cout<<" xx "<<xx<<" d_bin "<<d_bin<<" corr xtheta "<<(1.f-irotat_effective)*(-xtheta/180.*PI)<<" atan2 ";
						cout<<atan2(Ly,Lx)<<" theta/PI*8 "<<(atan2(Ly,Lx)+(1.f-irotat_effective)*(-xtheta/180.*PI))*fact_HD<<" ";
						cout<<" theta/PI*8 floor "<<theta<<" theta/PI*8 positif "<<(theta+nb_bin_Hdesc)%nb_bin_Hdesc<<" ";} */
//					if (theta<0) theta+=(int)(2.*PI*fact_HD); theta=theta%nb_bin_Hdesc; // pb d'arrondi ?????????????
					while (theta<0) theta+=nb_bin_Hdesc; theta=theta%nb_bin_Hdesc;
/*					if ((i>=33&&i<=43)&&(j>=60&&j<=70) || (i>=105&&i<=115)&&(j>=14&&j<=24)) cout<<theta<<"\n";*/
//					Tdescr_key_points[n][i_H*nb_bin_Hdesc+j_H*nb_bin_Hdesc+theta]+=1;
					Tdescr_key_points[n][i_H*nb_bin_Hdesc+j_H*nb_bin_Hdesc+theta]+=xx*(1.f-d_bin);
					Tdescr_key_points[n][i_H*nb_bin_Hdesc+j_H*nb_bin_Hdesc+theta+1]+=xx*d_bin;
				}
			}
/*			if ((i>=33&&i<=43)&&(j>=60&&j<=70) || (i>=105&&i<=115)&&(j>=14&&j<=24)) {
//			if ((i>=72&&i<=84)&&(j>=121&&j<=135) || (i>=75&&i<=96)&&(j>=155&&j<=182)) {
				for (theta=0; theta<lg_descr; theta++) {cout<<" "<<fixed<<setw(6)<<setprecision(0)<<Tdescr_key_points[n][theta]; if ((theta+1)%nb_bin_Hdesc==0) cout<<" |\n";}
				cout<<"\n **********\n"; char aa; cin>>aa;}*/
		}
	}

/* destructeur */
	~keypoints () {
		int n;
		for (n=0; n<n_key_points; n++) if (Tdescr_key_points[n]!=NULL) delete[] Tdescr_key_points[n];
		if (Tdescr_key_points!=NULL) delete[] Tdescr_key_points;
		for (n=0; n<n_key_points_max; n++) if (Tab_key_points[n]!=NULL) delete[] Tab_key_points[n];
		if (Tab_key_points!=NULL) delete[] Tab_key_points;
	}

/* constructeur de copie */
	keypoints (const keypoints &k_pts) {
		int n,i;
		n_key_points_max=k_pts.n_key_points_max;
		Tab_key_points=NULL; Tab_key_points=new float*[n_key_points_max];
		for (n=0; n<n_key_points_max; n++) {
			Tab_key_points[n]=new float[4];
			for (i=0; i<4; i++) Tab_key_points[n][i]=k_pts.Tab_key_points[n][i];
		}
		n_key_points=k_pts.n_key_points;
		Tdescr_key_points=NULL; Tdescr_key_points=new float*[n_key_points];
		lg_descr=k_pts.lg_descr;
		for (n=0; n<n_key_points; n++) {
			Tdescr_key_points[n]=new float[lg_descr];
			for (i=0; i<lg_descr; i++) Tdescr_key_points[n][i]=k_pts.Tdescr_key_points[n][i];
		}
	}

/* operateur d'affectation */
	keypoints& operator=(const keypoints &k_pts) {
		if (this != &k_pts) {
			int n,i;
			for (n=0; n<n_key_points_max; n++) if (Tab_key_points[n]!=NULL) delete[] Tab_key_points[n];
			if (Tab_key_points!=NULL) delete[] Tab_key_points;
			n_key_points_max=k_pts.n_key_points_max;
			Tab_key_points=new float*[n_key_points_max];
			for (n=0; n<n_key_points_max; n++) {
				Tab_key_points[n]=new float[4];
				for (i=0; i<4; i++) Tab_key_points[n][i]=k_pts.Tab_key_points[n][i];
			}
			for (n=0; n<n_key_points; n++) if (Tdescr_key_points[n]!=NULL) delete[] Tdescr_key_points[n];
			if (Tdescr_key_points!=NULL) delete[] Tdescr_key_points;
			n_key_points=k_pts.n_key_points;
			Tdescr_key_points=new float*[n_key_points];
			lg_descr=k_pts.lg_descr;
			for (n=0; n<n_key_points; n++) {
				Tdescr_key_points[n]=new float[lg_descr];
				for (i=0; i<lg_descr; i++) Tdescr_key_points[n][i]=k_pts.Tdescr_key_points[n][i];
			}		
		}
		return *this;
	}

/* accesseur n_key_points */
	int nb_keypoints() const {return n_key_points;}

/* operateur () */
	float operator() (int n, int i=0) const {
		if (i<0 || i>3 || n<0 || n>=n_key_points_max) {
			cout<<" debordement d'indice dans operator() de keypoints\n"; i=mini(maxi(0,i),3); n=mini(maxi(0,n),n_key_points_max-1);}
		return Tab_key_points[n][i];
	}

	void conv2imBYTE(imadata<BYTE> &ima) {
		ima.mise_a_zero();
		int nblig=ima.nlig(),nbcol=ima.ncol(),n,k,k2,i,j,i0,i2,j0,j2,ii,jj;
		for (n=0; n<n_key_points; n++)
			if (Tab_key_points[n][2]>=0) { 
				k=(int)Tab_key_points[n][2]+1; i=(int)Tab_key_points[n][0]; j=(int)Tab_key_points[n][1]; k2=k*k;
				i0=maxi(0,i-k); i2=mini(nblig-1,i+k); j0=maxi(0,j-k); j2=mini(nbcol-1,j+k);
				for (ii=i0; ii<=i2; ii++) for (jj=j0; jj<=j2; jj++) if ((ii-i)*(ii-i)+(jj-j)*(jj-j)<=k2) ima(ii,jj)=k;
/*				xtheta=((Tab_key_points[n][3]+0.5)*360.f/nb_bin_Horient); // xtheta encore en degres
				if (xtheta<=90 || xtheta>=270) {j0=j; j2=mini(nbcol-1,j+l_o*k);} else {j0=maxi(0,j-l_o*k); j2=j;}
				if (xtheta<=180) {i0=i; i2=mini(nblig-1,i+l_o*k);} else {i0=maxi(0,i-l_o*k); i2=i;} 
				xtheta=xtheta/180.*PI; // xtheta encore en radians
				if (fabs(xtheta-PI/2)>=0.01 && fabs(xtheta-3*PI/2)>=0.01) 
					for (jj=j0; jj<=j2; jj++) {ii=(jj-j)*tan(xtheta)+i; if (ii>=i0 && ii<=i2) ima_keypoints(mini(maxi(ii,0),nblig-1),jj)=k;} 
				if (fabs(xtheta-0)>=0.01 && fabs(xtheta-PI)>=0.01) 
					for (ii=i0; ii<=i2; ii++) {jj=(ii-i)*cos(xtheta)/sin(xtheta)+j; if (jj>=j0 && jj<=j2) ima_keypoints(ii,mini(maxi(jj,0),nbcol-1))=k;}*/
			} //cout<<" conv2imBYTE Ok\n";
	}

/* appariement de deux ensembles de points-clés sans descripteurs associés */ 
	int search_keypoints_corresp (const keypoints &k_pts, int *T_corr, int dim_T_corr, bool iaff=1) const {
		const int d=3;
		int n_key_points1=n_key_points, n_key_points2=k_pts.n_key_points, k1, k2, imax=-1, jmax=-1;
		for (k1=0; k1<n_key_points1; k1++) {
			if (Tab_key_points[k1][0]<imax) imax=(int)Tab_key_points[k1][0];	if (Tab_key_points[k1][1]<jmax) jmax=(int)Tab_key_points[k1][1];}
		for (k2=0; k2<n_key_points2; k2++) {
			if (k_pts.Tab_key_points[k2][0]<imax) imax=(int)k_pts.Tab_key_points[k2][0];	if (k_pts.Tab_key_points[k2][1]<jmax) jmax=(int)k_pts.Tab_key_points[k2][1];}
		int nblig=imax, nbcol=jmax;
		imadata<BYTE> imakeypts1(nblig,nbcol), imakeypts2(nblig,nbcol);
		imakeypts1.mise_a_zero(); 
		for (k1=0; k1<n_key_points1; k1++) {imakeypts1((int)Tab_key_points[k1][0],(int)Tab_key_points[k1][1])=255;}
		imakeypts2.mise_a_zero(); 
		for (k2=0; k2<n_key_points2; k2++) {imakeypts2((int)k_pts.Tab_key_points[k2][0],(int)k_pts.Tab_key_points[k2][1])=255;}
		float dx,dy,rapport,angle;
		imadata<float> imaFM=imakeypts1.transformee_FourierMellin(imakeypts1,dx,dy,rapport,angle,1);
		return 1;
/*		matrice2D<double> X(d,d), Y(d,2), Xt(d,d), XXt(d,d), XXt_1(d,d), Id(d,d);
		for (i=0; i<d; i++) Id(i,i)=1;
		do {
			for (k=0; k<d; k++) {
				k1=rand()%n_key_points1; k2=rand()%n_key_points2;
				X(k,0)=1; X(k,1)=Tab_key_points[k1][0]; X(k,2)=Tab_key_points[k1][1];
				if (d>5) {X(k,3)=pow(X(k,1),2); X(k,4)=pow(X(k,2),2); X(k,5)=X(k,1)*X(k,2);}
				Y(k,0)=k_pts.Tab_key_points[k2][0]; Y(k,1)=k_pts.Tab_key_points[k2][1];
			}
			Xt=X.transpo();
			XXt=matrice2D<double>(Xt,X,d,n,d); 
			XXt_1=XXt.inverse(invOK);
			cout<<" verification de l'inversion de X.Xt : inversion ";
			(matrice2D<double>(XXt,XXt_1,d,d,d)==Id && matrice2D<double>(XXt_1,XXt,d,d,d)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
			A=matrice2D<double>(matrice2D<double>(XXt_1,Xt,d,d,n),Y,d,n,2);
			matrice2D<double> D=Y-matrice2D<double>(X,A,n,d,2); 
			cout<<" erreur d'estimation par point d'amer\n"; D.affiche();
		} while ();
		return nb_pts_corr;*/
	}

/* appariement de deux ensembles de points-clés avec descripteurs associés */ 
	int search_descriptors_corresp (const keypoints &k_pts, int *T_corr, int dim_T_corr, bool iaff=1) const { 
		const float dist_non_assoc=5000.f/*FLT_MAX*/;
		if (lg_descr!=k_pts.lg_descr) {cout<<" descripteurs de tailles differentes -> pas de comparaison !!!\n"; return 0;}
		int n_key_points2=k_pts.n_key_points, nb_pts_corr=0, i, j, j_distmin, j_distmin2;
		double *Tdist=new double[n_key_points2], distmin, distmin2, dist;
		if (iaff) {
			for (i=0; i<n_key_points; i++)
				cout<<" im.1 point-cle "<<i<<" en ("<<Tab_key_points[i][0]<<","<<Tab_key_points[i][1]<<","<<Tab_key_points[i][2]<<") orientation(s) = "<<Tab_key_points[i][3]<<"\n";
			for (i=0; i<k_pts.n_key_points; i++)
				cout<<" im.2 point-cle "<<i<<" en ("<<k_pts.Tab_key_points[i][0]<<","<<k_pts.Tab_key_points[i][1]<<","<<k_pts.Tab_key_points[i][2]<<") orientation(s) = "<<k_pts.Tab_key_points[i][3]<<"\n";
		}
//		float **M2D_dist=new float*[n_key_points+n_key_points2];
//		for (i=0; i<n_key_points+n_key_points2; i++) {
//			M2D_dist[i]=new float[n_key_points+n_key_points2];
//			for (j=0; j<n_key_points+n_key_points2; j++) M2D_dist[i][j]=dist_non_assoc;
//		}
		for (i=0; i<dim_T_corr; i++) T_corr[i]=-1;
		int nmax_corr_ppoint=dim_T_corr/n_key_points;
		cout<<" nmax_corr_ppoint = "<<nmax_corr_ppoint<<"\n";
		for (i=0; i<n_key_points; i++) {
			distmin=distmin2=DBL_MAX; j_distmin=j_distmin2=-1;
			for (j=0; j<n_key_points2; j++) {
				Tdist[j]=DBL_MAX; 
				if (fabs(Tab_key_points[i][2]-k_pts.Tab_key_points[j][2])<=1) {
//					dist=distance_euclidienne<float>(Tdescr_key_points[i],k_pts.Tdescr_key_points[j],lg_descr);
					dist=distance_Battacharrya<float>(Tdescr_key_points[i],k_pts.Tdescr_key_points[j],lg_descr);
					Tdist[j]=dist; if (dist<distmin) {distmin=dist; j_distmin=j;}
//					M2D_dist[i][j]=(float)dist;
				}
			}
			if (nmax_corr_ppoint>2) {
				val_pos *T_dist_pos=new val_pos[n_key_points2];
				for (j=0; j<n_key_points2; j++) {T_dist_pos[j].pos=j; T_dist_pos[j].val=Tdist[j];}
//				cout<<" avant tri rapide : "; for (j=0; j<n_key_points2; j++) cout<<fixed<<setprecision(3)<<Tdist[j]<<" "; cout<<" \n";
				tri_rapide(T_dist_pos,n_key_points2);
//				cout<<" apres tri rapide : "; for (j=0; j<n_key_points2; j++) cout<<fixed<<setprecision(3)<<Tdist[j]<<" "; cout<<" \n";
				for (j=0; j<nmax_corr_ppoint; j++) {
//					T_corr[i*nmax_corr_ppoint+j]=T_dist_pos[n_key_points2-1-j].pos; nb_pts_corr++;
					T_corr[i*nmax_corr_ppoint+j]=T_dist_pos[j].pos; nb_pts_corr++;
				}
				if (iaff) {
					cout<<" dist min entre point ("<<fixed<<setprecision(0)<<Tab_key_points[i][0]<<","<<Tab_key_points[i][1]<<") de im.1 et celui (";
					cout<<k_pts.Tab_key_points[T_corr[i*nmax_corr_ppoint]][0]<<","<<k_pts.Tab_key_points[T_corr[i*nmax_corr_ppoint]][1]<<") de im.2 = ";
					cout<<setprecision(8)<<scientific<<T_dist_pos[0].val<<" et dist. min2 = "<<setprecision(8)<<scientific<<T_dist_pos[1].val<<" obtenue pour point (";
					cout<<fixed<<setprecision(0)<<k_pts.Tab_key_points[T_corr[i*nmax_corr_ppoint+1]][0]<<","<<k_pts.Tab_key_points[T_corr[i*nmax_corr_ppoint+1]][1]<<")\n";
				}
				if (T_dist_pos!=NULL) delete[] T_dist_pos;
			}
			else {
				T_corr[i*nmax_corr_ppoint]=j_distmin; nb_pts_corr++;
				for (j=0; j<n_key_points2; j++)
					if (j!=j_distmin && Tdist[j]<distmin2) {distmin2=Tdist[j]; j_distmin2=j;}
				if (nmax_corr_ppoint==2) {T_corr[i*nmax_corr_ppoint+1]=j_distmin2; nb_pts_corr++;}
				if (iaff) {
					cout<<" dist min entre point ("<<fixed<<setprecision(0)<<Tab_key_points[i][0]<<","<<Tab_key_points[i][1]<<") de im.1 et celui (";
					cout<<k_pts.Tab_key_points[j_distmin][0]<<","<<k_pts.Tab_key_points[j_distmin][1]<<") de im.2 = "<<setprecision(8)<<scientific<<distmin;
					cout<<" et dist. min2 = "<<setprecision(8)<<scientific<<distmin2<<" obtenue pour point (";
					cout<<fixed<<setprecision(0)<<k_pts.Tab_key_points[j_distmin2][0]<<","<<k_pts.Tab_key_points[j_distmin2][1]<<")\n";
				}
			}
/*			if (distmin/distmin2<0.8) {
				if (iaff) {
					cout<<" corresp. entre point ("<<fixed<<setprecision(0)<<Tab_key_points[i][0]<<","<<Tab_key_points[i][1]<<") de im.1 et celui (";
					cout<<k_pts.Tab_key_points[j_distmin][0]<<","<<k_pts.Tab_key_points[j_distmin][1]<<") de im.2 : dist = "<<setprecision(8)<<scientific<<distmin;
					cout<<" et ratio des dist. min et min2 = "<<setprecision(8)<<scientific<<distmin/distmin2<<"\n";
				}
				T_corr[i]=j_distmin; nb_pts_corr++;
			} else {
//				cout<<" pas de correspondance non ambigue pour le point-cle ("<<Tab_key_points[i][0]<<","<<Tab_key_points[i][1];
//				cout<<") de l'image 1 : ratio des distances min et min2 = "<<distmin/distmin2<<"\n";
			}*/
		}
//		for (i=0; i<n_key_points; i++) {
//			for (j=0; j<n_key_points2; j++) cout<<setw(6)<<setprecision(3)<<M2D_dist[i][j]/1000; cout<<"\n";}
//		cout<<" cout_min="<<algoKuhn (M2D_dist,n_key_points+n_key_points2)<<"\n";
//		for (i=0; i<n_key_points; i++)
//			for (j=0; j<n_key_points2; j++) 
//				if (M2D_dist[i][j]==0) {
//					cout<<" point-cle "<<i<<" ("<<Tab_key_points[i][0]<<","<<Tab_key_points[i][1]<<") en correspondance ";
//					cout<<" avec point-cle "<<j<<" ("<<Tab_key_points[j][0]<<","<<Tab_key_points[j][1]<<")\n";
//				}

		if (Tdist!=NULL) delete[] Tdist;
//		for (i=0; i<n_key_points+n_key_points2; i++) if (M2D_dist[i]!=NULL) delete[] M2D_dist[i];
//		if (M2D_dist!=NULL) delete[] M2D_dist;
		return nb_pts_corr;
	}
};

class imacontours : public imasites
{	int nblayer;
	BYTE *Tima;
	BYTE valnul;
 public:
/* constructeur par defaut */
	imacontours (int nl=1, int nc=1, int nd=1) : imasites(nl,nc) {
		Tima=NULL;
		nblayer=nd;
		valnul=0;
		Tima=new BYTE[nbpix*nblayer];
		for (long int i=0; i<nbpix; i++) {
			for (int j=0; j<nblayer; j++) Tima[i*nblayer+j]=valnul;
			Tsites[i]=&(Tima[i*nblayer]);
		}
	}

/* constructeur par conversion d'une image binaire en une image de contours */
	imacontours (imabin& ima) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL;
		nblayer=1;
		Tima=new BYTE[nbpix*nblayer];
		valnul=0;
		int i,j,ii;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				Tima[ii]=(BYTE)ima(i,j);
				Tsites[ii]=&(Tima[ii]);
			}
	}

/* constructeur par conversion d'une image en BYTE en une image de contours */
	imacontours (imadata<BYTE>& ima, int k=0) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL;
		nblayer=1;
		Tima=new BYTE[nbpix*nblayer];
		valnul=0;
		int i,j,ii;
		bool icont;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				icont=(ima(i,j,k)>0?1:0);
				Tima[ii]=(BYTE)icont;
				Tsites[ii]=&(Tima[ii]);
			}
	}

/* constructeur à partir d'une estimation du gradient par Prewitt, Sobel, ou MDIF */
	imacontours (imadata<float> &, char* ="MDIF", unsigned short int=0, bool=0);

/* constructeur à partir d'une estimation du Laplacien et du gradient */
	imacontours (imadata<float> &, char* ="8connex", char* ="Sobel", bool=0);

/* constructeur par filtrage optimal : Deriche ou Shen-Castan */
	imacontours (imadata<float> &, float=1., char* ="Deriche", bool=0);

/* constructeur par Huertas Medioni */ 
	imacontours (imadata<float> &, float=2., float=4232., bool=0);

/* constructeur par detection de points d'intérêt selon Harris */ 
//	imacontours (imadata<float> &, float=2., bool=0, char* ="Harris");

/* destructeur */
	~imacontours () {
		if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
	}

/* constructeur de copie */
	imacontours (const imacontours & ima) : imasites(ima) {
		Tima=NULL;
		nblayer=ima.nblayer;
		valnul=ima.valnul;
		Tima=new BYTE[nbpix*nblayer];
		BYTE * adval;
		for (long int i=0; i<nbpix; i++) {
			adval=(BYTE *)ima.Tsites[i];
			for (int j=0; j<nblayer; j++) Tima[i*nblayer+j]=*(adval+j);
			Tsites[i]=&(Tima[i*nblayer]);
		}
	}

/* operateur d'affectation */
	imacontours& operator=(const imacontours & ima) {
		if (this != &ima) {
			imasites *ad1, *ad2;
			ad1=this;
			ad2=(imasites *) &ima;
			*ad1=*ad2;
			if (Tima!=NULL) delete[] Tima;
			nblayer=ima.nblayer;
			valnul=ima.valnul;
			Tima=new BYTE[nbpix*nblayer];
			BYTE * adval;
			for (long int i=0; i<nbpix; i++) {
				adval=(BYTE *)ima.Tsites[i];
				for (int j=0; j<nblayer; j++) Tima[i*nblayer+j]=*(adval+j);
				Tsites[i]=&(Tima[i*nblayer]);
			}
		}
		return *this;
	}

/* # de couches */
	int nlayer() {
		return nblayer;
	}

/* operateur d'acces a la valeur (i,j,k) */
	BYTE& operator () (int i, int j, int k=0) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol || k<0 || k>=nblayer) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<","<<k<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
			if (j>=nbcol) j=nbcol-1;
			if (k<0) k=0;
			if (k>=nblayer) k=nblayer-1;
		}
		BYTE * adval=(BYTE *)Tsites[i*nbcol+j];
		adval=adval+k;
		return *adval;
	}

/* conversion en imadata<BYTE> */
	imadata<BYTE> conv2im1dBYTE (int k=0) {
		int i,j,x;
		imadata<BYTE> ima_u(nblig,nbcol);
		BYTE fact=255; 
		k=mini(maxi(k,0),nblayer-1);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				BYTE *adval=(BYTE*)Tsites[i*nbcol+j];
				x=int(*(adval+k));
				ima_u(i,j)=(BYTE)x*fact;
			}
		return ima_u;
	}

/* conversion en imadata<BYTE> */
	imadata<BYTE> conv2imBYTE () {
		int i,j,k,x;
		imadata<BYTE> ima_u(nblig,nbcol,nblayer);
		BYTE fact=255; 
		for (k=0; k<nblayer; k++)
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					BYTE *adval=(BYTE*)Tsites[i*nbcol+j];
					x=int(*(adval+k));
					ima_u(i,j,k)=(BYTE)x*fact;
				}
		return ima_u;
	}

/* calcul de l'image des maxima locaux de la norme du gradient dans la direction du gradient */
//	void imagrd_maxloc (imadata<float> &, imadata<float> &, imadata<BYTE> &);
	void imagrd_maxloc (imadata<float> &, imadata<float> &);

/* prolongation des contours en vue de leur fermeture */
	void prolongecontours (imadata<float> &, unsigned short int=5, float=-1.);
	void edgedrawing(imadata<float> &, unsigned short int=8, unsigned short int=4, float=-1.f, bool=0);

/* ecriture dans un fichier de sortie */
	void sauve_Ima (char* ="imaContours.dat", int=0) const;
};

struct noeud {
	int x,y,p;
	BYTE code; 
	unsigned int prof;
	double cout;
	bool stop;
};

#endif