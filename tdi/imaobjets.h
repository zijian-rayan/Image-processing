#ifndef _IMAOBJETS_H
#define _IMAOBJETS_H

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include "constantes.h"
#include "imabin.h"
#include "imaregions.h"

class objet
{	bool i_valid, i_boite_englobante, i_barycentre, i_perimetre, i_surface, i_cc, i_trou, i_nbEuler, i_dir_axis,
		 i_ratio, i_lab;
	unsigned int l_min, l_max, c_min, c_max;
	float l_baryc, c_baryc, perim, surf, dir_axis, delta_dir, lg_MinA, lg_MajA, ratio_MajMinA, err_ratio;
	int nb_cc[2], nb_trou[2], nbEuler[2];
	BYTE lab;
public:
	objet() {
		i_valid=i_boite_englobante=i_barycentre=i_nbEuler=i_perimetre=i_surface=i_dir_axis=i_ratio=i_cc=i_trou=i_lab=0;
		l_min=l_max=c_min=c_max=0; lab=0;
		l_baryc=c_baryc=perim=surf=dir_axis=0.f; delta_dir=360.f;
		lg_MinA=lg_MajA=ratio_MajMinA=0.f; err_ratio=FLT_MAX;
		for (int i=0; i<2; i++) nb_cc[i]=nb_trou[i]=nbEuler[i]=0;
	}
	objet(const objet& obj) {
		i_valid=obj.i_valid; i_boite_englobante=obj.i_boite_englobante; i_barycentre=obj.i_barycentre; 
		i_perimetre=obj.i_perimetre; i_surface=obj.i_surface; i_cc=obj.i_cc; i_trou=obj.i_trou; 
		i_nbEuler=obj.i_nbEuler; i_dir_axis=obj.i_dir_axis; i_ratio=obj.i_ratio; i_lab=obj.i_lab;
		lab=obj.lab; l_min=obj.l_min; l_max=obj.l_max; c_min=obj.c_min; c_max=obj.c_max;
		l_baryc=obj.l_baryc; c_baryc=obj.c_baryc; perim=obj.perim; surf=obj.surf;
		dir_axis=obj.dir_axis; delta_dir=obj.delta_dir; 
		lg_MinA=obj.lg_MinA; lg_MajA=obj.lg_MajA; ratio_MajMinA=obj.ratio_MajMinA; err_ratio=obj.err_ratio;
		for (int i=0; i<2; i++) {
			nb_cc[i]=obj.nb_cc[i]; nb_trou[i]=obj.nb_trou[i]; nbEuler[i]=obj.nbEuler[i];
		}
	}
	objet& operator= (const objet& obj) {
		if (this!=&obj) {
			i_valid=obj.i_valid; i_boite_englobante=obj.i_boite_englobante; i_barycentre=obj.i_barycentre; 
			i_perimetre=obj.i_perimetre; i_surface=obj.i_surface; i_cc=obj.i_cc; i_trou=obj.i_trou; 
			i_nbEuler=obj.i_nbEuler; i_dir_axis=obj.i_dir_axis; i_ratio=obj.i_ratio; i_lab=obj.i_lab;
			lab=obj.lab; l_min=obj.l_min; l_max=obj.l_max; c_min=obj.c_min; c_max=obj.c_max;
			l_baryc=obj.l_baryc; c_baryc=obj.c_baryc; perim=obj.perim; surf=obj.surf;
			dir_axis=obj.dir_axis; delta_dir=obj.delta_dir; 
			lg_MinA=obj.lg_MinA; lg_MajA=obj.lg_MajA; ratio_MajMinA=obj.ratio_MajMinA; err_ratio=obj.err_ratio;
			for (int i=0; i<2; i++) {
				nb_cc[i]=obj.nb_cc[i]; nb_trou[i]=obj.nb_trou[i]; nbEuler[i]=obj.nbEuler[i];
			}
		}
		return *this;
	}
	bool& valide() {return i_valid;}
	void affiche () {
		if (i_valid) {
			if (i_boite_englobante) cout<<" boite englobante : UL=(lig.="<<l_min<<",col.="<<c_min<<"), LR=(lig.="<<l_max<<",col.="<<c_max<<")\n";
			else cout<<" boite englobante non calculee\n";
			if (i_barycentre) cout<<" barycentre en (lig.="<<l_baryc<<",col.="<<c_baryc<<")\n";
			else cout<<" barycentre non calcule\n";
			if (i_perimetre) cout<<" perimetre = "<<perim<<" pixels\n";
			else cout<<" perimetre non calcule\n";
			if (i_surface) cout<<" surface = "<<surf<<" pixels\n";
			else cout<<" surface non calculee\n";
			if (i_dir_axis) cout<<" direction axe principal = "<<dir_axis<<" deg. a +/- "<<delta_dir<<" pres\n";
			else cout<<" direction axe principal non calculee\n";
			if (i_ratio) {
				cout<<" longueur petit_axe = "<<lg_MinA<<" & longueur grand_axe = "<<lg_MajA<<"\n";
				cout<<" ratio petit_axe / grand_axe = "<<ratio_MajMinA<<" a +/- "<<err_ratio<<" pres\n";
			}
			else cout<<" longueurs des axes et ratio petit_axe / grand_axe non calcules\n";
			if (i_cc) cout<<" nombre de composantes 4-connexes = "<<nb_cc[0]<<" et 8_connexes = "<<nb_cc[1]<<"\n";
			else cout<<" nombre de composantes connexes non calcule\n";
			if (i_trou) cout<<" nombre de 'trous' 8-connexes = "<<nb_trou[0]<<" et 4_connexes = "<<nb_trou[1]<<"\n";
			else cout<<" nombre de trous non calcule\n";
			if (i_nbEuler) cout<<" nombre d'Euler en 4-connexite = "<<nbEuler[0]<<" et 8-connexite = "<<nbEuler[1]<<"\n";
			else cout<<" nombre d'Euler non calcule\n";
			if (i_lab) cout<<" label = "<<(int)lab<<"\n";
			else cout<<" label non calcule\n";
		}
	}
	bool boite_englobante(int _lmin, int _lmax, int _cmin, int _cmax) {
		bool iOK=0;
		if (i_valid) {
			iOK=1;
			if (_lmin<=_lmax && _cmin<=_cmax) {
				l_min=_lmin; l_max=_lmax; 
				c_min=_cmin; c_max=_cmax;
				i_boite_englobante=1;
			}
			else {
				cout<<" coordonnees de la boite englobante ??? "<<_lmin<<" "<<_lmax<<" "<<_cmin<<" "<<_cmax<<"\n";
				iOK=0;
			}
			i_boite_englobante=iOK;
		}
		return iOK;
	}
	unsigned int l_min_obj() {if (i_boite_englobante) return l_min; else return UINT_MAX;}
	unsigned int l_max_obj() {if (i_boite_englobante) return l_max; else return UINT_MAX;}
	unsigned int c_min_obj() {if (i_boite_englobante) return c_min; else return UINT_MAX;}
	unsigned int c_max_obj() {if (i_boite_englobante) return c_max; else return UINT_MAX;}
	float l_baryc_obj() {if (i_barycentre) return l_baryc; else return -1;}
	float c_baryc_obj() {if (i_barycentre) return c_baryc; else return -1;}
	float perim_obj() {if (i_perimetre) return perim; else return -1;}
	float surf_obj() {if (i_surface) return surf; else return -1;}
	float dir_axis_obj() {if (i_dir_axis) return dir_axis; else return FLT_MAX;}
	float delta_dir_obj() {if (i_dir_axis) return delta_dir; else return FLT_MAX;}
	int nb_cc_obj(int c=4) {if (i_cc) {if (c==4) return nb_cc[0]; else return nb_cc[1];} else return -1;}
	int nb_trou_obj(int c=4) {if (i_trou) {if (c==4) return nb_trou[0]; else return nb_trou[1];} else return -1;}
	int nbEuler_obj(int c=4) {if (i_nbEuler) {if (c==4) return nbEuler[0]; else return nbEuler[1];} else return INT_MAX;}
	void barycentre(float _lbaryc, float _cbaryc) {
		if (i_valid) {
			l_baryc=_lbaryc; 
			c_baryc=_cbaryc;
			i_barycentre=1;
		}
	}
	void surface(float _surf) {
		if (i_valid) {
			surf=_surf; 
			i_surface=1;
		}
	}
	void nb_compconnexes(int _nb_c4c, int _nb_c8c=0) {
		if (i_valid) {
			nb_cc[0]=_nb_c4c;
			nb_cc[1]=_nb_c8c;
			i_cc=1;
		}
	}
	void nb_Euler(int _n_E4c, int _n_E8c=0) {
		if (i_valid) {
			nbEuler[0]=_n_E4c;
			nbEuler[1]=_n_E8c;
			i_nbEuler=1;
		}
	}
	void nb_trous(int _nb_t4c, int _nb_t8c=0) {
		if (i_valid) {
			nb_trou[0]=_nb_t4c;
			nb_trou[1]=_nb_t8c;
			i_trou=1;
		}
	}
	friend class imaobjets;
};

class imaobjets : public imasites
{	int nobj;
	BYTE *Tima;
	BYTE valnul;
	objet *T_objets;
 	bool i_basic_param, i_basic_topo, i_second_param, i_classif;
public:
/* constructeur par defaut */
	imaobjets (int nl=1, int nc=1, int no=0) : imasites(nl,nc) {
		Tima=NULL; valnul=0;
		Tima=new BYTE[nbpix];
		for (long int i=0; i<nbpix; i++) {
			Tima[i]=valnul;
			Tsites[i]=&(Tima[i]);
		}
		nobj=no;
		if (nobj<=0) T_objets=NULL;
		else T_objets=new objet[nobj];
		i_basic_param=i_basic_topo=i_second_param=0;
	}
/* constructeur par conversion d'une image en BYTE en une image d'objets */
	imaobjets (imadata<BYTE> &ima, int k=0, unsigned int obj0=0) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL; Tima=new BYTE[nbpix];
		valnul=0;
		int i,j,l,objN=(int)obj0; long int ii;
		const int val_max=255;
		bool T_noObj[val_max+1];
		for (l=0; l<=val_max; l++) T_noObj[l]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				Tima[ii]=(BYTE)ima(i,j,k);
				nobj=(int)Tima[ii];
				if (!T_noObj[nobj]) {
					T_noObj[nobj]=1; 
					cout<<" objet de label "<<nobj<<" trouve\n";
					if (nobj>objN) objN=nobj; 
				}
				Tsites[ii]=&(Tima[ii]);
			}
		nobj=0;
		for (l=0; l<=objN; l++)
			if (T_noObj[l]) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) 
						if ((*this)(i,j)==l) (*this)(i,j)=nobj;
				nobj++;
			}
		cout<<nobj<<" objets dans l'image\n";
		T_objets=new objet[nobj];
//		cout<<" calcul des parametres geometriques de base : \n"; 
		basic_param ();
//		cout<<" calcul des parametres topologiques de base : \n"; 
		basic_topologie (obj0);
//		cout<<" calcul des parametres secondaires : \n"; 
		second_param (obj0);
	}

/* constructeur par conversion d'une image de regions en une image d'objets */
	imaobjets (imaregions &ima, int k=0, unsigned int obj0=0) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL; Tima=new BYTE[nbpix];
		valnul=0;
		int i,j,l,objN=(int)obj0; long int ii;
		const int val_max=ima.nregions(); cout<<" # regions = "<<val_max<<"\n";
		bool *T_noObj=new bool[val_max+1];
		for (l=0; l<=val_max; l++) T_noObj[l]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				nobj=maxi((int)(ima(i,j,k)),0);
				Tima[ii]=nobj;
				if (!T_noObj[nobj]) {
					T_noObj[nobj]=1; 
//					cout<<" objet de label "<<nobj<<" trouve en ("<<i<<","<<j<<")\n";
					if (nobj>objN) objN=nobj; 
				}
				Tsites[ii]=&(Tima[ii]);
			}
		nobj=0;
		for (l=0; l<=objN; l++)
			if (T_noObj[l]) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) 
						if ((*this)(i,j)==l) (*this)(i,j)=nobj;
				nobj++;
			}
		cout<<nobj<<" objets dans l'image\n";
		if (T_noObj!=NULL) {delete[] T_noObj; T_noObj=NULL;}
		T_objets=new objet[nobj];
		basic_param ();
		basic_topologie (obj0);
		second_param (obj0);
	}
/* destructeur */
	~imaobjets () {
		if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
		if (T_objets!=NULL) {delete[] T_objets; T_objets=NULL;}
	}
/* constructeur de recopie */
	imaobjets (const imaobjets &ima) : imasites(ima) {
		Tima=NULL; Tima=new BYTE[nbpix];
		valnul=ima.valnul;
		BYTE *adval;
		for (long int i=0; i<nbpix; i++) {
			adval=(BYTE *)ima.Tsites[i];
			Tima[i]=*(adval);
			Tsites[i]=&(Tima[i]);
		}
		nobj=ima.nobj;
		T_objets=new objet[nobj];
		for (int j=0; j<nobj; j++)
			T_objets[j]=ima.T_objets[j];
//		basic_param ();
//		basic_topologie ();
	}
/* operateur d'affectation */
	imaobjets& operator= (const imaobjets &ima) {
		if (this != &ima) {
			imasites *ad1, *ad2;
			ad1=this;
			ad2=(imasites *) &ima;
			*ad1=*ad2;
			if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
			Tima=new BYTE[nbpix];
			valnul=ima.valnul;
			BYTE * adval;
			for (long int i=0; i<nbpix; i++) {
				adval=(BYTE *)ima.Tsites[i];
				Tima[i]=*(adval);
				Tsites[i]=&(Tima[i]);
			}
			if (T_objets!=NULL) delete[] T_objets;
			nobj=ima.nobj;
			T_objets=new objet[nobj];
			for (int j=0; j<nobj; j++)
				T_objets[j]=ima.T_objets[j];
//			basic_param ();
//			basic_topologie ();
		}
		return *this;
	}
/* operateur d'acces a la valeur (i,j) */
	BYTE& operator () (int i, int j) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
			if (j>=nbcol) j=nbcol-1;
		}
		BYTE * adval=(BYTE *)Tsites[i*nbcol+j];
		adval=adval;
		return *adval;
	}
/* conversion en imadata<BYTE> */
	imadata<BYTE> conv2imBYTE() const {
		int i,j,x;
		imadata<BYTE> ima_u(nblig,nbcol);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				BYTE *adval=(BYTE*)Tsites[i*nbcol+j];
				x=int(*(adval));
				if (x>=0) ima_u(i,j)=(BYTE)x;
			}
		return ima_u;
	}
/* # objets */
	int nb_obj () {return nobj;}
/* accesseur objet */
	objet objet_i(int i) {return T_objets[i];}
/* affichage */
	void affiche (bool=0);
/* difference entre deux images de type imaobjets */
	imaobjets operator - (imaobjets &);
/* ecriture dans un fichier de sortie */
	void sauve_Ima  (char *nomfich="imaObjsauvee.dat") const;
	void sauve1d_Ima(char *nomfich="imaObj1dsauvee.dat") const;
/* parametres basiques (barycentre, surface, boites englobante) des objets */
	void basic_param (unsigned int=UINT_MAX);
/* nombre de composantes connexes */
	void basic_topologie (unsigned int=UINT_MAX);
/* parametres du 2ème ordre (direction axe principal) des objets */
	void second_param (unsigned int=UINT_MAX);
/* detection objets de forme disque */
	imabin detection_disques (float=0.9, float=1.1);
/* regroupement d'objets alignés */
	void regroupe_alignements (unsigned int=0);
	void regroupe_2objets (unsigned int, unsigned int, unsigned int**);
	void regroupe_alignements (bool**, unsigned int**, float, float=0, float=0, unsigned int=0);
/* classification d'objets de type marquage route */
	imalabels clas_marquages_route (unsigned int=0, float=0.5f, bool=0);
};

class objetsgeom
{protected:
  int dx, dy; //dimensions de la boite englobante
	float xO, yO; //coordonnées du centroide relativement à la boite englobante
	bool *imOG;
public:
  objetsgeom(float _dx=0, float _dy=0): dx((int)ceil(_dx)), dy((int)ceil(_dy)) {
//		xO=around(dx/2.); yO=around(dy/2.);
		xO=dx/2.f; yO=dy/2.f;
		imOG=NULL; imOG=new bool[dx*dy];
		for (int i=0; i<dx*dy; i++) imOG[i]=false;
	}
	virtual ~objetsgeom() {cout<<" entree dans destructeur ~objetsgeom() :"; 
		if (imOG!=NULL) {delete[] imOG; cout<<" imOG detruit\n"; imOG=NULL; } else cout<<" inutile... \n";
	}
  objetsgeom(const objetsgeom &O) {
		dx=O.dx; dy=O.dy; xO=O.xO; yO=O.yO; 
		imOG=NULL; imOG=new bool[dx*dy];
		for (int i=0; i<dx*dy; i++) imOG[i]=O.imOG[i];
	}
	objetsgeom& operator= (const objetsgeom &O) {
		if (this!=&O) {
			dx=O.dx; dy=O.dy; xO=O.xO; yO=O.yO; 
			if (imOG!=NULL) {delete[] imOG; imOG=NULL;}
			imOG=new bool[dx*dy];
			for (int i=0; i<dx*dy; i++) imOG[i]=O.imOG[i];
		}
		return *this;
	}
	virtual void Affiche() const {
		int n, n_p=0, n_s=0, n8c, ii, jj, ii0, ii2, jj0, jj2;
		imabin imab(dy+2,dx+2); // on rajoute un bord nul autour de l'objet pour pouvoir faie erosion dilatation correctement
		for (int i=0; i<dy; i++) {
			ii0=maxi(i-1,0); ii2=mini(i+1,dy-1);
			for (int j=0; j<dx; j++) {
				n=i*dx+j;
				cout<<" "<<(int)(imOG[n]);
				if (imOG[n]) {
					imab(i+1,j+1)=1;
					n_s++;
					n8c=0; jj0=maxi(j-1,0); jj2=mini(j+1,dx-1);
					for (ii=ii0; ii<=ii2; ii++)
						for (jj=jj0; jj<=jj2; jj++) if (imOG[ii*dx+jj]) n8c++;
					n_p+=(int)(0.5*(9.-n8c));
/*					if ((j==0 || j>0 && !imOG[i*dx+j-1]) || (j==dx-1 || j<dx-1 && !imOG[i*dx+j+1]) || 
						(i==0 || i>0 && !imOG[(i-1)*dx+j]) || (i==dy-1 || i<dy-1 && !imOG[(i+1)*dx+j])) n_p++;*/
/*					if (j==0 || j>0 && !imOG[i*dx+j-1]) n_p++;
					if (j==dx-1 || j<dx-1 && !imOG[i*dx+j+1]) n_p++;
					if (i==0 || i>0 && !imOG[(i-1)*dx+j]) n_p++;
					if (i==dy-1 || i<dy-1 && !imOG[(i+1)*dx+j]) n_p++;*/
				}
			}
			cout<<"\n";
		}
		cout<<"  Surface  = "<<Surface()<<" approx. "<<n_s<<"\n";
		cout<<" Perimetre = "<<Perimetre()<<" approx. "<<n_p<<"\n";
		eltstruct S3(3,3);
//		(imab-imab.erode(S3)).affiche(); (imab.dilate(S3)-imab).affiche();
//		float p8_min=(imab-imab.erode(S3)).norm(), p8_max=(imab.dilate(S3)-imab).norm();
//		cout<<" Perimetre bord 8connexite approx. in ["<<p8_min<<","<<p8_max<<"]\n";
		S3(0,0)=S3(0,2)=S3(2,0)=S3(2,2)=0;
//		(imab-imab.erode(S3)).affiche(); (imab.dilate(S3)-imab).affiche();
		float p4_min=(imab-imab.erode(S3)).norm(), p4_max=(imab.dilate(S3)-imab).norm();
		cout<<" Perimetre bord 4connexite approx. in ["<<p4_min<<","<<p4_max<<"]\n";
	}
	void add_2_ImaNg(imadata<BYTE> &imaObNg, double mu, double sigma, int yUL, int xUL) {
		int nl=imaObNg.nlig(), nc=imaObNg.ncol(), i,j,ii,jj;
		double var=sigma*sigma, xx;
		for (ii=0; ii<dy; ii++) {
			i=ii+yUL;
			if (i>=0 && i<nl)
				for (jj=0; jj<dx; jj++) {
					j=jj+xUL;
					if (j>=0 && j<nc) {
						if (imOG[ii*dx+jj]==true) {
							xx=tirage_gauss(var,mu);
							xx=mini(maxi(xx,0.),255.);
							imaObNg(i,j)=mini(maxi((BYTE)xx,imaObNg(i,j)),(BYTE)255); // suppose des objets clairs sur fond sombre
						}
					}
				}
		}
	}
  virtual float Perimetre() const = 0;
  virtual float Surface() const = 0;
	int dim_x() {return dx;}
	int dim_y() {return dy;}
//	void Infos();
	friend class imaobjetsgeom;
};

class rectangle : public objetsgeom 
{ float largeur, longueur, angle_dx; // angle_dx est l'angle entre le grand axe du rectangle (qui donne la longueur) et l'axe horizontal
	double xA, yA, xB, yB, xC, yC, xD, yD;
public:
	rectangle (float _largeur, float _longueur, float _angle_dx) : largeur(_largeur), longueur(_longueur), angle_dx(_angle_dx/180.f*acos(-1.f)) {
		double demilong=(double)longueur/2., demilarg=(double)largeur/2., costheta=cos(_angle_dx/180.*acos(-1.)), sintheta=sin(_angle_dx/180.*acos(-1.)), x, y;
		xA=demilong*costheta-demilarg*sintheta; yA=demilong*sintheta+demilarg*costheta; 
		xB=-demilong*costheta-demilarg*sintheta; yB=-demilong*sintheta+demilarg*costheta; 
		xC=-demilong*costheta+demilarg*sintheta; yC=-demilong*sintheta-demilarg*costheta; 
		xD=demilong*costheta+demilarg*sintheta; yD=demilong*sintheta-demilarg*costheta;
		dx=(int)(maxi(xA,maxi(xB,maxi(xC,xD)))-mini(xA,mini(xB,mini(xC,xD)))+1);
		dy=(int)(maxi(yA,maxi(yB,maxi(yC,yD)))-mini(yA,mini(yB,mini(yC,yD)))+1);
		xO=-around(2*mini(xA,mini(xB,mini(xC,xD))))/2.f-0.5f; yO=-around(2*mini(yA,mini(yB,mini(yC,yD))))/2.f-0.5f; // decalage de 0.5 pour mettre les valeurs 'dans' les pixels et non 'entre' les pixels
//		xO=(dx-2)/2.f; yO=(dy-2)/2.f;
		imOG=NULL; imOG=new bool[dx*dy];
		int i,j;
		for (i=0; i<dx*dy; i++) imOG[i]=false;
		for (i=0; i<dx; i++) 
			for (j=0; j<dy; j++) {
				x=(i-xO)*costheta+(j-yO)*sintheta;
				y=-(i-xO)*sintheta+(j-yO)*costheta;
				if (x>=-demilong && x<demilong && y>=-demilarg && y<demilarg) imOG[j*dx+i]=true;
			}
	}
	void Affiche() const {
    cout<<" Vignette de taille "<<dy<<"x"<<dx<<" representant le rectangle de taille "<<longueur<<"x"<<largeur<<",\n";
		cout<<"centre en ("<<yO<<","<<xO<<"), d'orientation "<<angle_dx*180/PI<<", et de sommets ("<<xA<<","<<yA<<"),\n";
		cout<<"("<<xB<<","<<yB<<"), ("<<xC<<","<<yC<<"), ("<<xD<<","<<yD<<")\n";
		objetsgeom::Affiche();
	}
	float Perimetre() const {return 2*(largeur+longueur);}
	float Surface() const {return largeur*longueur;}
};

class strip : public objetsgeom
{ float largeur, angle_dx, width; // angle_dx est l'angle entre le grand axe du rectangle (qui donne la longueur) et l'axe horizontal,
																	// width est la largeur de la vignette representant la bande 
	double xA, yA, xB, yB, xC, yC, xD, yD;
public:
	strip (float _largeur, float _angle_dx, float _width) : largeur(_largeur), angle_dx(_angle_dx/180.f*acos(-1.f)), width(_width) {
		const double eps=1.e-6;
		double costheta=eps*around(cos(_angle_dx/180.*acos(-1.))/eps), sintheta=eps*around(sin(_angle_dx/180.*acos(-1.))/eps);
		cout<<" ++++++++++"<<1./pow(2.,0.5)<<" "<<fabs(costheta)<<" "<<sintheta<<"\n";
		double longueur=width/maxi(fabs(costheta),fabs(sintheta)); cout<<" longueur = "<<longueur<<"\n";
		double demilong=longueur/2., demilarg=largeur/2., x, y;
		xA=demilong*costheta-demilarg*sintheta; yA=demilong*sintheta+demilarg*costheta; 
		xB=-demilong*costheta-demilarg*sintheta; yB=-demilong*sintheta+demilarg*costheta; 
		xC=-demilong*costheta+demilarg*sintheta; yC=-demilong*sintheta-demilarg*costheta; 
		xD=demilong*costheta+demilarg*sintheta; yD=demilong*sintheta-demilarg*costheta;
		double dec_x=mini(xA,mini(xB,mini(xC,xD))), dec_y=mini(yA,mini(yB,mini(yC,yD)));
		xA-=dec_x; xB-=dec_x; xC-=dec_x; xD-=dec_x; yA-=dec_y; yB-=dec_y; yC-=dec_y; yD-=dec_y; 
		dx=(int)ceil(maxi(xA,maxi(xB,maxi(xC,xD))));
		dy=(int)ceil(maxi(yA,maxi(yB,maxi(yC,yD)))); cout<<" dx = "<<dx<<" dy = "<<dy<<"\n";
//		dx=maxi(xA,maxi(xB,maxi(xC,xD)))-mini(xA,mini(xB,mini(xC,xD)));
//		dy=maxi(yA,maxi(yB,maxi(yC,yD)))-mini(yA,mini(yB,mini(yC,yD)));
		xO=around(dx)/2.f-0.5f; yO=around(dy)/2.f-0.5f; // decalage de 0.5 pour mettre les valeurs 'dans' les pixels et non 'entre' les pixels
//		xO=-around(2*mini(xA,mini(xB,mini(xC,xD))))/2.f-0.5f; yO=-around(2*mini(yA,mini(yB,mini(yC,yD))))/2.f-0.5f; // decalage de 0.5 pour mettre les valeurs 'dans' les pixels et non 'entre' les pixels
		cout<<" Vignette de taille "<<dy<<"x"<<dx<<" representant la bande de largeur "<<largeur<<" et longueur "<<longueur<<",\n";
		cout<<"centre en ("<<yO<<","<<xO<<"), d'orientation "<<angle_dx*180/PI<<", et de sommets ("<<yA<<","<<xA<<"),\n";
		cout<<"("<<yB<<","<<xB<<"), ("<<yC<<","<<xC<<"), ("<<yD<<","<<xD<<")\n";
		double xM=(demilong-0.5)*costheta+xO, yM=(demilong-0.5)*sintheta+yO, xN=-(demilong-0.5)*costheta+xO, yN=-(demilong-0.5)*sintheta+yO,
			     xP=-(demilarg-0.5)*sintheta+xO, yP=(demilarg-0.5)*costheta+yO, xQ=(demilarg-0.5)*sintheta+xO, yQ=-(demilarg-0.5)*costheta+yO;
		cout<<" centres des segments en ("<<yM<<","<<xM<<"), ("<<yN<<","<<xN<<"), ("<<yP<<","<<xP<<"), ("<<yQ<<","<<xQ<<")\n";
		int x_0=around(mini(xM,mini(xN,mini(xP,xQ)))), y_0=around(mini(yM,mini(yN,mini(yP,yQ)))), 
			  x_2=around(maxi(xM,maxi(xN,maxi(xP,xQ)))), y_2=around(maxi(yM,maxi(yN,maxi(yP,yQ))));
		cout<<" y_0 = "<<y_0<<" x_0 = "<<x_0<<" y_2 = "<<y_2<<" x_2 = "<<x_2<<"\n";
		dx=x_2-x_0+1; dy=y_2-y_0+1; xO-=x_0; yO-=y_0; cout<<" dx = "<<dx<<" dy = "<<dy<<"\n";
		if (imOG!=NULL) {cout<<" destruction de imOG avant reallocation\n"; delete[] imOG; imOG=NULL;} imOG=new bool[dx*dy];
		int i,j;
		for (i=0; i<dx*dy; i++) imOG[i]=false;
		imOG[around((yM-y_0))*dx+around(xM-x_0)]=imOG[around((yN-y_0))*dx+around(xN-x_0)]=imOG[around((yP-y_0))*dx+around(xP-x_0)]=imOG[around((yQ-y_0))*dx+around(xQ-x_0)]=true;
		for (i=0; i<dx; i++) 
			for (j=0; j<dy; j++) {
				x=(i-xO)*costheta+(j-yO)*sintheta;
				y=-(i-xO)*sintheta+(j-yO)*costheta;
				if (x>=-demilong && x<demilong && y>=-demilarg && y<demilarg) imOG[j*dx+i]=true;
			}
	}
	void Affiche() const {
		objetsgeom::Affiche();
	}
	float Perimetre() const {
		double costheta=cos(angle_dx), sintheta=sin(angle_dx);
		return fabs(costheta)>1./pow(2.,0.5)?(float)(2*dx/fabs(costheta)):(float)(2*dy/fabs(sintheta));}
	float Surface() const {
		double costheta=cos(angle_dx), sintheta=sin(angle_dx);
		return fabs(costheta)>1./pow(2.,0.5)?(float)(largeur*dx/fabs(costheta)):(float)(largeur*dy/fabs(sintheta));}
};

class imaobjetsgeom : public imasites
{	int nobj;
	BYTE *Tima;
	BYTE valnul;
//	objet *T_objets;
// 	bool i_basic_param, i_basic_topo, i_second_param, i_classif;
public:
/* constructeur par defaut */
	imaobjetsgeom (int nl=1, int nc=1, int no=0) : imasites(nl,nc) {
		valnul=0;
		Tima=NULL; Tima=new BYTE[nbpix];
		for (long int i=0; i<nbpix; i++) {
			Tima[i]=valnul;
			Tsites[i]=&(Tima[i]);
		}
		nobj=no;
//		if (nobj<=0) T_objets=NULL;
//		else T_objets=new objet[nobj];
//		i_basic_param=i_basic_topo=i_second_param=0;
	}
/* constructeur par conversion d'une image binaire en une image d'objets */
	imaobjetsgeom (imabin &ima) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL; Tima=new BYTE[nbpix];
		valnul=0;
		int i,j; long int ii;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				Tima[ii]=(BYTE)ima(i,j);
				Tsites[ii]=&(Tima[ii]);
			}
		nobj=1;
	}
/* constructeur par conversion d'une image en BYTE en une image d'objets */
	imaobjetsgeom (imadata<BYTE> &ima, int k=0, unsigned int obj0=0) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL; Tima=new BYTE[nbpix];
		valnul=0;
		int i,j,l,objN=(int)obj0; long int ii;
		const int val_max=255;
		bool T_noObj[val_max+1];
		for (l=0; l<=val_max; l++) T_noObj[l]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				Tima[ii]=(BYTE)ima(i,j,k);
				nobj=(int)Tima[ii];
				if (!T_noObj[nobj]) {
					T_noObj[nobj]=1; cout<<" objet de label "<<nobj<<" trouve\n";
					if (nobj>objN) objN=nobj; 
				}
				Tsites[ii]=&(Tima[ii]);
			}
		nobj=0;
		for (l=0; l<=objN; l++)
			if (T_noObj[l]) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) 
						if ((*this)(i,j)==l) (*this)(i,j)=nobj;
				nobj++;
			}
		cout<<nobj<<" objets dans l'image\n";
//		T_objets=new objet[nobj];
//		cout<<" calcul des parametres geometriques de base : \n"; 
//		basic_param ();
//		cout<<" calcul des parametres topologiques de base : \n"; 
//		basic_topologie (obj0);
//		cout<<" calcul des parametres secondaires : \n"; 
//		second_param (obj0);
	}

/* constructeur par conversion d'une image de regions en une image d'objets */
	imaobjetsgeom (imaregions &ima, int k=0, unsigned int obj0=0) : imasites(ima.nlig(),ima.ncol()) {
		Tima=NULL; Tima=new BYTE[nbpix];
		valnul=0;
		int i,j,l,objN=(int)obj0; long int ii;
		const int val_max=ima.nregions(); cout<<" # regions = "<<val_max<<"\n";
		bool *T_noObj=new bool[val_max+1];
		for (l=0; l<=val_max; l++) T_noObj[l]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii=i*nbcol+j;
				nobj=maxi((int)(ima(i,j,k)),0);
				Tima[ii]=nobj;
				if (!T_noObj[nobj]) {
					T_noObj[nobj]=1; // cout<<" objet de label "<<nobj<<" trouve en ("<<i<<","<<j<<")\n";
					if (nobj>objN) objN=nobj; 
				}
				Tsites[ii]=&(Tima[ii]);
			}
		nobj=0;
		for (l=0; l<=objN; l++)
			if (T_noObj[l]) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) 
						if ((*this)(i,j)==l) (*this)(i,j)=nobj;
				nobj++;
			}
		cout<<nobj<<" objets dans l'image\n";
//		T_noObj=NULL; T_objets=new objet[nobj];
//		basic_param ();
//		basic_topologie (obj0);
//		second_param (obj0);
	}
/* destructeur */
	~imaobjetsgeom () {
		if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
//		if (T_objets!=NULL) {delete[] T_objets; T_objets=NULL;}
	}
/* constructeur de recopie */
	imaobjetsgeom (const imaobjetsgeom &ima) : imasites(ima) {
		Tima=NULL; Tima=new BYTE[nbpix];
		valnul=ima.valnul;
		BYTE *adval;
		for (long int i=0; i<nbpix; i++) {
			adval=(BYTE *)ima.Tsites[i];
			Tima[i]=*(adval);
			Tsites[i]=&(Tima[i]);
		}
		nobj=ima.nobj;
//		T_objets=new objet[nobj];
//		for (int j=0; j<nobj; j++) T_objets[j]=ima.T_objets[j];
//		basic_param ();
//		basic_topologie ();
	}
/* operateur d'affectation */
	imaobjetsgeom& operator= (const imaobjetsgeom &ima) {
		if (this != &ima) {
			imasites *ad1, *ad2;
			ad1=this;
			ad2=(imasites *) &ima;
			*ad1=*ad2;
			if (Tima!=NULL) {delete[] Tima; Tima=NULL;}
			Tima=new BYTE[nbpix];
			valnul=ima.valnul;
			BYTE * adval;
			for (long int i=0; i<nbpix; i++) {
				adval=(BYTE *)ima.Tsites[i];
				Tima[i]=*(adval);
				Tsites[i]=&(Tima[i]);
			}
//			if (T_objets!=NULL) delete[] T_objets;
			nobj=ima.nobj;
//			T_objets=new objet[nobj];
//			for (int j=0; j<nobj; j++) T_objets[j]=ima.T_objets[j];
//			basic_param ();
//			basic_topologie ();
		}
		return *this;
	}
/* operateur d'acces a la valeur (i,j) */
	BYTE& operator () (int i, int j) {
		if (i<0 || i>=nblig || j<0 || j>=nbcol) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<")\n";
			if (i<0) i=0;
			if (i>=nblig) i=nblig-1;
			if (j<0) j=0;
			if (j>=nbcol) j=nbcol-1;
		}
		return *((BYTE *)Tsites[i*nbcol+j]);
	}
	BYTE operator () (const int i, const int j) const {
		int ii=i, jj=j;
		if (i<0 || i>=nblig || j<0 || j>=nbcol) {
			cout<<" debordement d''indice dans ("<<i<<","<<j<<")\n";
			if (i<0) ii=0;
			if (i>=nblig) ii=nblig-1;
			if (j<0) jj=0;
			if (j>=nbcol) jj=nbcol-1;
		}
		return *((BYTE *)Tsites[ii*nbcol+jj]);
	}
/* conversion en imadata<BYTE> */
	imadata<BYTE> conv2imBYTE() const {
		int i,j,x;
		imadata<BYTE> ima_u(nblig,nbcol);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				BYTE *adval=(BYTE*)Tsites[i*nbcol+j];
				x=int(*(adval));
				if (x>=0) ima_u(i,j)=(BYTE)x;
			}
		return ima_u;
	}
/* # objets */
	int nb_obj () {return nobj;}
/* initialisation a zero */
	void mise_a_zero() {
		int i,j;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) (*this)(i,j)=0;
	}
/* ajout d'un objet à l'image des objets */
	void add_obj_ima(objetsgeom&,int=0,int=0,float=1);
/* accesseur objet */
//	objet objet_i(int i) {return T_objets[i];}
/* affichage */
//	void affiche (bool=0);
/* difference entre deux images de type imaobjets */
//	imaobjets operator - (imaobjets &);
/* ecriture dans un fichier de sortie */
//	void sauve_Ima  (char *nomfich="imaObjsauvee.dat") const;
//	void sauve1d_Ima(char *nomfich="imaObj1dsauvee.dat") const;
/* parametres basiques (barycentre, surface, boites englobante) des objets */
//	void basic_param (unsigned int=UINT_MAX);
/* nombre de composantes connexes */
//	void basic_topologie (unsigned int=UINT_MAX);
/* parametres du 2ème ordre (direction axe principal) des objets */
//	void second_param (unsigned int=UINT_MAX);
/* detection objets de forme disque */
//	imabin detection_disques (float=0.9, float=1.1);
/* regroupement d'objets alignés */
//	void regroupe_alignements (unsigned int=0);
//	void regroupe_2objets (unsigned int, unsigned int, unsigned int**);
//	void regroupe_alignements (bool**, unsigned int**, float, float=0, float=0, unsigned int=0);
/* classification d'objets de type marquage route */
//	imalabels clas_marquages_route (unsigned int=0, float=0.5f, bool=0);
	void trace_ellipse (double, double, double, double, double, BYTE=255, BYTE=128);
};

#endif