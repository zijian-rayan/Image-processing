#ifndef _PARAMS2_H
#define _PARAMS2_H

#include "stdafx.h"
#include "tdi/tdi.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <windows.h>
#include <string.h>

using namespace std;
#include "math.h"
static bool strip_valide;

class type_strips{
public:
	int nl;//nb lignes de l'image où se trouve les bandes
	int nc;//nb colonnes
	int nb_res;//nb des bandes
	float *angle;//tableau des angles des bandes
	int *large;//tableau des largueur des bandes
	int *decalage;//tableau des decalages des bandes
	double *res;//tableau des resultats NFA des bandes
	int numero_max;//le numero du resultat NFA le plus grand
	double res_max;//le resultat NFA le plus grand
	
	type_strips(int nb_resultat)//constructeur classe
	{
		angle=new float[nb_resultat];
		decalage=new int[nb_resultat];
		res=new double[nb_resultat];
		large=new int[nb_resultat];
		nb_res=nb_resultat;
	}
	type_strips(type_strips& st)//constructeur copie
	{
		nb_res=st.nb_res;
		angle=new float[nb_res];
		decalage=new int[nb_res];
		res=new double[nb_res];
		large=new int[nb_res];
		nl=st.nl;
		nc=st.nc;
		numero_max=st.numero_max;
		res_max=st.res_max;
		int i;
		for(i=0;i<nb_res;i++)
		{
			angle[i]=st.angle[i];
			decalage[i]=st.decalage[i];
			res[i]=st.res[i];
			large[i]=st.large[i];
		}
	}
	
	~type_strips()//destructeur
	{
		delete[] angle;
		delete[] decalage;
		delete[] res;
		delete[] large;
	}
	//methode pour afficher les bandes
	ostream& Show(ostream& os) const
	{
		os<<"Le nombre des strips interessants dans cette image de "<<nl<<"x"<<nc<<": "<<nb_res<<"\n";
		for(int i=0;i<nb_res;i++)
		{
			os<<"largeur = "<<large[i]<<" ; decalage = "<<decalage[i]<<" ; angle = "<<angle[i]<<" ; NF = "<<res[i];
			if(i==numero_max)	os<<"		Valeur MAX";
			os<<"\n";
		}
		os<<"NFA Max: "<<res_max<<"  strip numero "<<numero_max<<"\n\n";
		return os;
	}
	//addition de 2 classe de type_strips 
	type_strips add(const type_strips st) const
	{
		int nb_total=nb_res+st.nb_res;
		type_strips resultat(nb_total);
		
		if((nl!=st.nl)||(nc!=st.nc))	cout<<"ATTENTION!!! Les tailles d'image ne sont pas pareilles!"; 
		resultat.nl=nl;
		resultat.nc=nc;
		
		int i;
		for(i=0;i<nb_res;i++)
		{
			resultat.angle[i]=angle[i];
			resultat.decalage[i]=decalage[i];
			resultat.res[i]=res[i];
			resultat.large[i]=large[i];
		}
		
		for(i=0;i<st.nb_res;i++)
		{
			resultat.angle[i+nb_res]=st.angle[i];
			resultat.decalage[i+nb_res]=st.decalage[i];
			resultat.res[i+nb_res]=st.res[i];
			resultat.large[i+nb_res]=st.large[i];
		}
		
		if(res_max>=st.res_max)
		{
			resultat.numero_max=numero_max;
			resultat.res_max=res_max;
		}
		else
		{
			resultat.numero_max=st.numero_max+nb_res;
			resultat.res_max=st.res_max;
		}
		return resultat;
	}
	//ajouter une bande dans la classe
	type_strips add(float _angle,int _large,int _decalage,double _res)
	{
		type_strips resultat(nb_res+1);
		
		resultat.nl=nl;
		resultat.nc=nc;
		
		int i;
		for(i=0;i<nb_res;i++)
		{
			resultat.angle[i]=angle[i];
			resultat.decalage[i]=decalage[i];
			resultat.res[i]=res[i];
			resultat.large[i]=large[i];
		}
		
		resultat.angle[nb_res]=_angle;
		resultat.decalage[nb_res]=_decalage;
		resultat.res[nb_res]=_res;
		resultat.large[nb_res]=_large;
		
		if(_res>res_max)
		{
			resultat.numero_max=nb_res;
			resultat.res_max=_res;
		}
		else
		{
			resultat.numero_max=numero_max;
			resultat.res_max=res_max;
		}
		return resultat;
	}
	//verifie si une bande est dans la classe
	bool belong(float _angle,int _large,int _decalage)
	{
		for(int i=0;i<nb_res;i++)
		{
			if((angle[i]==_angle)&&(decalage[i]==_decalage)&&(large[i]==_large))	return true;
		}
		return false;
	}
	//rappel de la fonction de creation de bande --------------------------------------------------------???????????????????????????????????????
	static imabin gestalts_creation_strip(int nl,int nc,float angle,int decalage,int large)
	{
		imabin imb(nl,nc);
		int i,k;
		float y1=0,y2=0,x1=0,x2=0;
		if((angle<=45)||((angle>=135)&&(angle<=225))||((angle>=315)&&(angle<360)))
		{
			for(i=0;i<nc;i++)
			{			  
				y1=-decalage/cos(angle*(float)PI/180)-tan(angle*(float)PI/180)*(i-nc/2)+nl/2;
				y2=-(decalage+large)/cos(angle*(float)PI/180)-tan(angle*(float)PI/180)*(i-nc/2)+nl/2;
				
				if(y1>y2)
				{	
					for(k=(int)y2;k<y1;k++)
						if((k>=0)&&(k<nl))	imb(k,i)=1;
				}
				else
				{
					for(k=(int)y2;k>y1;k--)
						if((k>=0)&&(k<nl))	imb(k,i)=1;
				}
			}
		}
		else
		{
			for(i=0;i<nl;i++)
			{			  
				x1=-decalage/sin(angle*(float)PI/180)+(nl/2-i)/tan(angle*(float)PI/180)+nc/2;
				x2=-(decalage+large)/sin(angle*(float)PI/180)+(nl/2-i)/tan(angle*(float)PI/180)+nc/2;
				
				if(x1>x2)
				{	
					for(k=(int)x2;k<x1;k++)
						if((k>=0)&&(k<nc))	imb(i,k)=1;
				}
				else
				{
					for(k=(int)x2;k>x1;k--)
						if((k>=0)&&(k<nc))	imb(i,k)=1;
				}
			}
			//  }
		}
		return imb;
	}
	//tracer tous les bandes dans la classe par une image binaire
	imabin trace()
	{
		imabin ima(nl,nc);
		int i;
		for(i=0;i<nb_res;i++)
		{
			imabin ima_strip=gestalts_creation_strip(nl,nc,angle[i],decalage[i],large[i]);
			ima=ima||ima_strip;
		}
		return ima;
	}
	//operateur egal
	type_strips& operator =(const type_strips st)
	{
		nb_res=st.nb_res;
		angle=new float[nb_res];
		decalage=new int[nb_res];
		res=new double[nb_res];
		large=new int[nb_res];
		nl=st.nl;
		nc=st.nc;
		numero_max=st.numero_max;
		res_max=st.res_max;
		int i;
		for(i=0;i<nb_res;i++)
		{
			angle[i]=st.angle[i];
			decalage[i]=st.decalage[i];
			res[i]=st.res[i];
			large[i]=st.large[i];
		}
		return *this;
	}
	
}; 

//ostream& operator<<(ostream& os, const type_strips& f) {return f.Show_strips(os);};

class type_fenetres{
public:
	int nl;//nb ligne de l'image ou se trouve les bandes
	int nc;//nb colonne
	int nb_res;//nb des fenetres
	int taille;//taille de fenetres
	int *compte_horizontal;//tableau des positions horizontales des fenetres(compter en decalage de la taille de la fenetre)
	int *compte_vertical;//tableau des positions verticales des fenetres
	double *res;//tableau des resultats NFA des bandes
	int numero_max;//le numero du resultat NFA le plus grand
	double res_max;//le resultat NFA le plus grand
	
	type_fenetres(int nb_resultat)//constructeur classe
	{
		compte_horizontal=new int[nb_resultat];
		compte_vertical=new int[nb_resultat];
		res=new double[nb_resultat];
		nb_res=nb_resultat;
	}

	type_fenetres(type_fenetres& st)//constructeur copie
	{
		nb_res=st.nb_res;
		compte_horizontal=new int[nb_res];
		compte_vertical=new int[nb_res];
		res=new double[nb_res];
		nl=st.nl;
		nc=st.nc;
		taille=st.taille;
		numero_max=st.numero_max;
		res_max=st.res_max;
		int i;
		for(i=0;i<nb_res;i++)
		{
			compte_horizontal[i]=st.compte_horizontal[i];
			compte_vertical[i]=st.compte_vertical[i];
			res[i]=st.res[i];
		}
	}
	
	~type_fenetres()//destructeur
	{
		delete[] compte_horizontal;
		delete[] compte_vertical;
		delete[] res;
	}
	
	//methode pour afficher les bandes
	ostream& Show(ostream& os) const
	{
		os<<"Le nombre des fenetres interessants dans cette image de "<<nl<<"x"<<nc<<": "<<nb_res<<"  la taille de fenetre: "<<taille<<"\n";
		for(int i=0;i<nb_res;i++)
		{
			os<<"position horizontale = "<<compte_horizontal[i]<<" ; position verticale = "<<compte_vertical[i]<<" ; NF = "<<res[i];
			if(i==numero_max)	os<<"		Valeur MAX";
			os<<"\n";
		}
		os<<"NFA Max: "<<res_max<<"  fenetre numero "<<numero_max<<"\n\n";
		return os;
	}

	
	//addition de 2 classe de type_strips 
	type_fenetres add(const type_fenetres st) const
	{
		int nb_total=nb_res+st.nb_res;
		type_fenetres resultat(nb_total);
		
		if((nl!=st.nl)||(nc!=st.nc))	cout<<"ATTENTION!!! Les tailles d'image ne sont pas pareilles!"; 
		if(taille!=st.taille)	cout<<"ATTENTION!!! La tailles de fenetre ne sont pas pareilles!"; 
		resultat.nl=nl;
		resultat.nc=nc;
		resultat.taille=taille;
		
		int i;
		for(i=0;i<nb_res;i++)
		{
			resultat.compte_horizontal[i]=compte_horizontal[i];
			resultat.compte_vertical[i]=compte_vertical[i];
			resultat.res[i]=res[i];
		}
		
		for(i=0;i<st.nb_res;i++)
		{
			resultat.compte_horizontal[i]=st.compte_horizontal[i];
			resultat.compte_vertical[i]=st.compte_vertical[i];
			resultat.res[i+nb_res]=st.res[i];
		}
		
		if(res_max>=st.res_max)
		{
			resultat.numero_max=numero_max;
			resultat.res_max=res_max;
		}
		else
		{
			resultat.numero_max=st.numero_max+nb_res;
			resultat.res_max=st.res_max;
		}
		return resultat;
	}
	
	//ajouter une bande dans la classe
	type_fenetres add(int cpt_hor,int cpt_ver,double _res)
	{
		type_fenetres resultat(nb_res+1);
		
		resultat.nl=nl;
		resultat.nc=nc;
		resultat.taille=taille;
		
		int i;
		for(i=0;i<nb_res;i++)
		{
			resultat.compte_horizontal[i]=compte_horizontal[i];
			resultat.compte_vertical[i]=compte_vertical[i];
			resultat.res[i]=res[i];
		}
		
		resultat.compte_horizontal[nb_res]=cpt_hor;
		resultat.compte_vertical[nb_res]=cpt_ver;
		resultat.res[nb_res]=_res;
		
		if(_res>res_max)
		{
			resultat.numero_max=nb_res;
			resultat.res_max=_res;
		}
		else
		{
			resultat.numero_max=numero_max;
			resultat.res_max=res_max;
		}
		return resultat;
	}
	
	//vertifier si une bande est dans la classe
	bool belong(int cpt_hor,int cpt_ver)
	{
		for(int i=0;i<nb_res;i++)
			if((compte_horizontal[i]==cpt_hor)&&(compte_vertical[i]==cpt_ver))	return true;

		return false;
	}
	
	//tracer tous les bandes dans la classe par une image binaire
	imabin trace()
	{
		imabin ima(nl,nc);
		int i,m,n;
		for(i=0;i<nb_res;i++)
			for(m=compte_horizontal[i]*taille;m<taille*(compte_horizontal[i]+1);m++)
				for(n=compte_vertical[i]*taille;n<taille*(compte_vertical[i]+1);n++)
					if((m>=0)&&(m<nc)&&(n>=0)&&(n<nl))	ima(n,m)=1;
		return ima;
	}
	
	//operateur egal
	type_fenetres& operator =(const type_fenetres st)
	{
		nb_res=st.nb_res;
		compte_horizontal=new int[nb_res];
		compte_vertical=new int[nb_res];
		res=new double[nb_res];
		nl=st.nl;
		nc=st.nc;
		taille=st.taille;
		numero_max=st.numero_max;
		res_max=st.res_max;
		int i;
		for(i=0;i<nb_res;i++)
		{
			compte_horizontal[i]=st.compte_horizontal[i];
			compte_vertical[i]=st.compte_vertical[i];
			res[i]=st.res[i];
		}
		return *this;
	}
	
}; 







  /**************************************************************************************/

#endif // _PARAMS2_H