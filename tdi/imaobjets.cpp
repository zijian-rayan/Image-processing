#include "constantes.h"
#include "imaobjets.h"

void imaobjets::affiche(bool iaf_ima) {
	int i,j;
	if (iaf_ima) {
		BYTE *adval, no_obj;
		for (i=0; i<nblig; i++) {
			for (j=0; j<nbcol; j++) {
				adval=(BYTE *)Tsites[i*nbcol+j];
				no_obj=*(adval);
				cout<<(int)no_obj<<" ";
			}
			cout<<"\n";
		}
		cout<<"----------------------------------------------------------\n";
	}
	for (i=0; i<nobj; i++) {
		if (T_objets[i].valide()) {
			cout<<" objet no."<<i<<" : \n";
			T_objets[i].affiche();
		}
	}
}

imaobjets imaobjets::operator - (imaobjets &ima2) {
	int nlig=mini(nblig,ima2.nblig), ncol=mini(nbcol,ima2.nbcol), i, j;
	imaobjets imaRes(nlig,ncol);
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++) {
			if ((*this)(i,j)==ima2(i,j)) imaRes(i,j)=1;
			else imaRes(i,j)=0;
		}
	return imaRes;
}

void imaobjets::sauve_Ima(char *nomfich) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);
	}
	BYTE* adval=(BYTE*)Tsites[0];
	int itype, sizeval;
	itype=0;
	sizeval=sizeof(BYTE);
	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col., de type ";
	cout<<" bytes => itype = "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (long int i=0; i<nbpix; i++) {
		adval=(BYTE*)Tsites[i];
		sortie.write((char *)adval,sizeval);
	}
	sortie.close();
}

void imaobjets::sauve1d_Ima(char *nomfich) const {
	imadata<BYTE> I1c=conv2imBYTE();
	I1c.sauve_ImaBSQ(nomfich);
}

void imaobjets::basic_param(unsigned int obj0) {
	int no, n, lmin, lmax, cmin, cmax, i, j;
	float lbary, cbary;
	bool iOK;
	for (no=0; no<nobj; no++) {
		if (no!=obj0) {
			lbary=cbary=0.f; n=0;
			lmin=nblig; cmin=nbcol; lmax=0; cmax=0;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if ((*this)(i,j)==no) {
						n++;
						lbary+=i; cbary+=j;
						if (i<lmin) lmin=i; if (i>lmax) lmax=i;
						if (j<cmin) cmin=j; if (j>cmax) cmax=j;
					}
			if (n>0) {
				T_objets[no].valide()=1;
				lbary/=n; cbary/=n;
				T_objets[no].barycentre(lbary,cbary);
//				cout<<" barycentre en ("<<lbary<<","<<cbary<<"), i_barycentre = "<<T_objets[no].i_barycentre<<"\n";
				T_objets[no].surface((float)n);
//				cout<<" surface = "<<n<<", i_surface = "<<T_objets[no].i_surface<<"\n";
				iOK=T_objets[no].boite_englobante(lmin,lmax,cmin,cmax);
//				cout<<" boite englobante de ("<<lmin<<","<<cmin<<") a ("<<lmax<<","<<cmax<<"), i_boite_englobante = "<<T_objets[no].i_boite_englobante<<"\n";
			}
		}
	}
	i_basic_param=1;
}

void imaobjets::second_param(unsigned int obj0) {
	const double eps=1.e-6;
	int no, i, j, n, imin1, imax1, jmin1, jmax1, imin2, imax2, jmin2, jmax2;
	if (!i_basic_param) basic_param();
	float i_b,j_b;
	double alpha_m, alpha_v, alpha, d, dmax, cos2_a, sin2_a, cos_asin_a, c_x1, c_y1, 
		   c_x2, c_y2, xp, yp;
	imabin imab(nblig,nbcol), imabsq(nblig,nbcol);
	imadata<BYTE> imabp(nblig,nbcol), imasq(nblig,nbcol);
	for (no=0; no<nobj; no++) { //cout<<" objet "<<no<<"\n";
		if (no!=obj0) {
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					if ((*this)(i,j)==no) imab(i,j)=1;
					else imab(i,j)=0;
				}
			imabsq=imab.squelette(8,1,0);
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) if (imabsq(i,j)) imasq(i,j)=1;
			alpha_m=alpha_v=0.0f; n=0; dmax=0.;
			i_b=T_objets[no].l_baryc; j_b=T_objets[no].c_baryc;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if (imabsq(i,j) && (fabs(i-i_b)+fabs(j-j_b)>1)) { 
						alpha=atan2(-(double)(i-i_b),(double)(j-j_b))/PI*180.; if (fabs(alpha+180)<1.e-3) alpha+=360;
						if (n>0 && fabs(alpha+180-alpha_m/n)<fabs(alpha-alpha_m/n)) alpha+=180.;
						else
							if (n>0 && fabs(alpha-180-alpha_m/n)<fabs(alpha-alpha_m/n)) alpha-=180.;
						alpha_m+=alpha; alpha_v+=pow(alpha,2);
						d=pow(pow((double)(i-i_b),2)+pow((double)(j-j_b),2),0.5);
						if (d>dmax) dmax=d;
						n++;
					}
			if (n>1) {
				alpha_v=alpha_v/(n-1.)-pow(alpha_m/n,2)*n/(n-1.);
				alpha_m/=n; if (alpha_m>=180) alpha_m-=180; if (fabs(alpha_m-180)<1.e-3) alpha_m=0;
				T_objets[no].dir_axis=(float)alpha_m;
				T_objets[no].delta_dir=(float)maxi(pow(alpha_v,0.5),atan(2.f/dmax)*180/PI);
				T_objets[no].i_dir_axis=1;
			}
			if (T_objets[no].i_dir_axis) {
				imin1=jmin1=imin2=jmin2=INT_MAX;
				imax1=jmax1=imax2=jmax2=INT_MIN;
				cos2_a=pow(cos(alpha_m/180*PI),2); sin2_a=pow(sin(alpha_m/180*PI),2); cos_asin_a=cos(alpha_m/180*PI)*sin(alpha_m/180*PI);
				c_x1=sin2_a*j_b+cos_asin_a*i_b; c_y1=cos_asin_a*j_b+cos2_a*i_b;
				c_x2=cos2_a*j_b-cos_asin_a*i_b; c_y2=-cos_asin_a*j_b+sin2_a*i_b;
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						if (imab(i,j)) {
							xp=j*cos2_a-i*cos_asin_a+c_x1;
							yp=-j*cos_asin_a+i*sin2_a+c_y1;
							if (xp>=0 && xp<nbcol-1 && yp>=0 && yp<nblig-1) {
								imabp((int)yp,(int)xp)=no;
								if (yp<imin1) imin1=(int)yp; if (yp>imax1) imax1=(int)yp;
								if (xp<jmin1) jmin1=(int)xp; if (xp>jmax1) jmax1=(int)xp;
							}
							xp=j*sin2_a+i*cos_asin_a+c_x2;
							yp=j*cos_asin_a+i*cos2_a+c_y2;
							if (xp>=0 && xp<nbcol-1 && yp>=0 && yp<nblig-1) {
								imabp((int)yp,(int)xp)=no+100;
								if (yp<imin2) imin2=(int)yp; if (yp>imax2) imax2=(int)yp;
								if (xp<jmin2) jmin2=(int)xp; if (xp>jmax2) jmax2=(int)xp;
							}
						}
/*						for (int ii=T_objets[no].l_min_obj(); ii<T_objets[no].l_max_obj(); ii++)
							for (int jj=T_objets[no].c_min_obj(); jj<T_objets[no].c_max_obj(); jj++) {
								if (fabs((cos(alpha_m/180*PI)*jj+sin(alpha_m/180*PI)*(-ii))-(cos(alpha_m/180*PI)*j+sin(alpha_m/180*PI)*(-i)))<1 && 
									fabs((-sin(alpha_m/180*PI)*jj+cos(alpha_m/180*PI)*(-ii))-(-sin(alpha_m/180*PI)*j_b+cos(alpha_m/180*PI)*(-i_b)))<1) {
									imabp(ii,jj)=no;
								}
								if (fabs(-(sin(alpha_m/180*PI)*jj+cos(alpha_m/180*PI)*ii)+(sin(alpha_m/180*PI)*j+cos(alpha_m/180*PI)*i))<1 && 
									fabs(-(cos(alpha_m/180*PI)*jj-sin(alpha_m/180*PI)*ii)+(cos(alpha_m/180*PI)*j_b-sin(alpha_m/180*PI)*i_b))<1) {
									imabp(ii,jj)=100+no;
								}
							}*/
						if (imabsq(i,j)) imabp(i,j)=255;
					}
				T_objets[no].lg_MinA=pow(pow(jmax2-jmin2+1.f,2)+pow(imax2-imin2+1.f,2),0.5f);
				T_objets[no].lg_MajA=pow(pow(jmax1-jmin1+1.f,2)+pow(imax1-imin1+1.f,2),0.5f);
				T_objets[no].ratio_MajMinA=T_objets[no].lg_MinA/T_objets[no].lg_MajA;
//				T_objets[no].ratio_MajMinA=pow(pow(jmax2-jmin2+1,2)+pow(imax2-imin2+1,2),0.5)/pow(pow(jmax1-jmin1+1,2)+pow(imax1-imin1+1,2),0.5);
				T_objets[no].err_ratio=pow(pow(jmax2-jmin2+2.f,2)+pow(imax2-imin2+2.f,2),0.5f)/maxi(1.f,pow(pow((float)(jmax1-jmin1),2)+pow((float)(imax1-imin1),2),0.5f))-
                               pow(pow((float)(jmax2-jmin2),2)+pow((float)(imax2-imin2),2),0.5f)/pow(pow(jmax1-jmin1+2.f,2)+pow(imax1-imin1+2.f,2),0.5f);
				T_objets[no].i_ratio=1;
			}
		}
	}
	cout<<" sauvegarde de l'image du squelette\n"; imasq.sauve_ImaPGM("ima_squelette.pgm"); 
	cout<<" sauvegarde de l'image des axes des objets\n"; imabp.sauve_ImaPGM("ima_axeobjets.pgm");
//	affiche(); char aa; cin>>aa;
/*	double *thetaU=new double[nobj], *thetaB=new double[nobj], *thetaL=new double[nobj], *thetaR=new double[nobj];
	for (no=0; no<nobj; no++) 
		thetaU[no]=thetaB[no]=thetaL[no]=thetaR[no]=-999;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			no=(*this)(i,j);
			if (no!=obj0) {
				if (i==T_objets[no].l_min && thetaU[no]==-999) {
					thetaU[no]=atan2(-(double)(i-T_objets[no].l_baryc),(double)(j-T_objets[no].c_baryc))/PI*180.;
				}
				if (j==T_objets[no].c_min && thetaL[no]==-999) {
					thetaL[no]=atan2(-(double)(i-T_objets[no].l_baryc),(double)(j-T_objets[no].c_baryc))/PI*180.;
				}
			}
		}
	for (i=nblig-1; i>=0; i--)
		for (j=nbcol-1; j>=0; j--) {
			no=(*this)(i,j);
			if (no!=obj0) {
				if (i==T_objets[no].l_max && thetaB[no]==-999) {
					thetaB[no]=atan2(-(double)(i-T_objets[no].l_baryc),(double)(j-T_objets[no].c_baryc))/PI*180.;
				}
				if (j==T_objets[no].c_max && thetaR[no]==-999) {
					thetaR[no]=atan2(-(double)(i-T_objets[no].l_baryc),(double)(j-T_objets[no].c_baryc))/PI*180.;
				}
			}
		}
	for (no=0; no<nobj; no++) {
		if (no!=obj0) {
			if (thetaU[no]<0) thetaU[no]+=180.;
			if (thetaL[no]<0) thetaL[no]+=180.;
			if (thetaB[no]<0) thetaB[no]+=180.;
			if (thetaR[no]<0) thetaR[no]+=180.;
			cout<<no<<" "<<T_objets[no].dir_axis<<" "<<thetaU[no]<<" "<<thetaB[no]<<" "<<thetaL[no]<<" "<<thetaR[no]<<"\n";
		}
	}
	if (thetaU!=NULL) delete[] thetaU; if (thetaB!=NULL) delete[] thetaB;
	if (thetaL!=NULL) delete[] thetaL; if (thetaR!=NULL) delete[] thetaR;*/
	i_second_param=1;
}

void imaobjets::basic_topologie(unsigned int obj0) {
	if (!i_basic_param) basic_param();
	int no, i, j, n_c4c, n_c8c, s, a, d, t, q, n_E4c, n_E8c;
	for (no=0; no<nobj; no++)
		if (no!=obj0 && T_objets[no].valide()==1) { //cout<<" objet "<<no<<"\n";
			imabin imab(nblig,nbcol);
			imadata<int> imacc(nblig,nbcol);
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if ((*this)(i,j)==no) imab(i,j)=1; 
//			cout<<" # composantes_connexes 4-connex = ";
			imacc=imab.composantes_connexes (n_c4c,4,0); //cout<<n_c4c<<"\n";
//			cout<<" # composantes_connexes 8-connex = ";
			imacc=imab.composantes_connexes (n_c8c,8,0); //cout<<n_c8c<<"\n";
			T_objets[no].nb_compconnexes(n_c4c,n_c8c); 
			s=a=d=t=q=0;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if (imab(i,j)==1) {
						s++;
						if (i>0 && imab(i-1,j)==1) a++;
						if (j>0 && imab(i,j-1)==1) a++;
						if (i>0 && j>0 && imab(i-1,j-1)==1) d++;
						if (i>0 && j<nbcol-1 && imab(i-1,j+1)==1) d++;
						if (i>0 && j>0 && imab(i,j-1)==1 && imab(i-1,j-1)==1) t++;
						if (i>0 && j<nbcol-1 && imab(i-1,j)==1 && imab(i-1,j+1)==1) t++;
						if (i<nblig-1 && j<nbcol-1 && imab(i,j+1)==1 && imab(i+1,j+1)==1) t++;
						if (i<nblig-1 && j>0 && imab(i+1,j)==1 && imab(i+1,j-1)==1) t++;
						if (i>0 && j>0 && imab(i,j-1)==1 && imab(i-1,j-1)==1 && imab(i-1,j)==1) q++;
					}
			n_E4c=s-a+q;
			n_E8c=s-a-d+t-q;
			T_objets[no].nb_Euler(n_E4c,n_E8c);
			T_objets[no].nb_trous(n_c4c-n_E4c,n_c8c-n_E8c);
		}
	i_basic_topo=1;
}

void imaobjets::regroupe_alignements (unsigned int obj0) {
	const float ddmax=45.f;
	const int rES=3; // precision de la direction pour rES=1 : 22.5°, rES=2 : 13.3°, rES=3 : 9.2°, rES=4 : 7.0°, rES=5 : 5.7°, rES=6 : 4.7°, rES=7 : 4.1°
	const int itmax=(int)pow(nblig*nblig+nbcol*nbcol,0.5);
	if (!i_second_param) second_param(obj0);
	int no,i,j,it,no2,nobj2=nobj;
	double d, dd;
	imabin imab_no(nblig,nbcol), imab_nn(nblig,nbcol), imab(nblig,nbcol);
	BYTE *Teq=new BYTE[nobj];
	for (no=0; no<nobj; no++) {
		Teq[no]=no;
		if (no!=obj0) {
			cout<<" objet "<<no<<" : ";
			if (T_objets[no].i_dir_axis) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						imab_no(i,j)=((*this)(i,j)==no);
						imab_nn(i,j)=((*this)(i,j)!=no && (*this)(i,j)!=obj0);
					}
				imab=imab_no&&imab_nn;
				d=T_objets[no].dir_axis;
				dd=T_objets[no].delta_dir;
				cout<<" direction axe principal "<<d;
				if (dd<=ddmax) {
					cout<<", precision "<<dd<<"\n";
					eltstruct ES(2*rES+1,d);
					it=0;
					while (imab.norm()==0 && it<itmax) {
						imab_no=imab_no.dilate(ES); it++;
						imab=imab_no&&imab_nn;
					}
					if (imab.norm()!=0) {
						for (i=0; i<nblig; i++)
							for (j=0; j<nbcol; j++)
								if (imab(i,j)) {
									no2=(*this)(i,j);
									cout<<" objet "<<(int)(*this)(i,j)<<" touche en "<<it<<" dilatations\n";
								}
						Teq[no]=Teq[no2]=mini(Teq[no],Teq[no2]);
						nobj2--;
					} else {
						cout<<" pas d'objet touche en "<<it<<" dilatations\n";
					}
				}
				else {
					cout<<" mal definie\n";
				}
			}
			else {
				cout<<" direction axe principal non definie\n";
			}
		}
	}
	if (nobj2<nobj) {
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				(*this)(i,j)=Teq[(*this)(i,j)];
		nobj=nobj2;
		cout<<" apres fusion des objets alignes : "<<nobj<<" objets dans l'image\n";
		if (T_objets!=NULL) delete[] T_objets;
		T_objets=new objet[nobj];
		basic_param ();
		basic_topologie (obj0);
		second_param (obj0);
	}
	if (Teq!=NULL) delete[] Teq;
}

void imaobjets::regroupe_2objets (unsigned int no1, unsigned int no2, unsigned int **Tdist) {
	int k, i_no1, i_no2;
	unsigned int d;
	for (k=0; k<nobj; k++) {
		d=mini(Tdist[k][no1],Tdist[k][no2]);
		Tdist[k][no2]=Tdist[no2][k]=Tdist[k][no1]=Tdist[no1][k]=d;
	}
	const float y=T_objets[no1].l_baryc-T_objets[no2].l_baryc, x=T_objets[no1].c_baryc-T_objets[no2].c_baryc; 
	const float s=T_objets[no1].surf+T_objets[no2].surf; 
	float z, a=atan2(-y,x)/(float)PI*180.f; if (a<0) a+=180.f; if (fabs(a-180)<1.e-3) a=0;
	if (T_objets[no1].i_dir_axis) {
		if (!T_objets[no2].i_dir_axis || T_objets[no1].delta_dir<T_objets[no2].delta_dir) {i_no1=no1; i_no2=no2;}
		else {i_no1=no2; i_no2=no1;}
	} else {
		if (T_objets[no2].i_dir_axis) {i_no1=no2; i_no2=no1;}
		else {i_no1=i_no2=-1;}
	}
	if (i_no1>=0 && i_no2>=0) {
		if (fabs(a+180-T_objets[i_no1].dir_axis)<fabs(a-T_objets[i_no1].dir_axis)) a+=180.;
		else 
			if (fabs(a-180-T_objets[i_no1].dir_axis)<fabs(a-T_objets[i_no1].dir_axis)) a-=180.;
		z=(T_objets[i_no1].dir_axis*T_objets[i_no1].surf+a*T_objets[i_no2].surf)/s;
	} else z=a; 
	T_objets[no1].dir_axis=T_objets[no2].dir_axis=z;
	T_objets[no1].i_dir_axis=T_objets[no2].i_dir_axis=1;
/*	if (T_objets[no1].l_baryc<T_objets[no2].l_baryc) {
		i=mini(T_objets[no1].l_baryc+1,nblig-1)-maxi(T_objets[no2].l_baryc-1,0);
		if (T_objets[no1].c_baryc>T_objets[no2].c_baryc)
			j=mini(T_objets[no1].c_baryc-1,0)-mini(T_objets[no2].l_baryc+1,nbcol-1);
		else
			j=mini(T_objets[no1].c_baryc+1,nbcol-1)-mini(T_objets[no2].l_baryc-1,0);
	} else {
		i=maxi(T_objets[no1].l_baryc-1,0)-mini(T_objets[no2].l_baryc+1,nblig-1);
		if (T_objets[no1].c_baryc>T_objets[no2].c_baryc)
			j=mini(T_objets[no1].c_baryc-1,0)-mini(T_objets[no2].l_baryc+1,nbcol-1);
		else
			j=mini(T_objets[no1].c_baryc+1,nbcol-1)-mini(T_objets[no2].l_baryc-1,0);
	}
	xx=atan2(-(double)i,T_objets[no1].c_baryc-T_objets[no2].c_baryc)/PI*180.; 
	yy=atan2(-T_objets[no1].l_baryc+T_objets[no2].l_baryc,(double)j)/PI*180.;
	T_objets[no1].delta_dir=T_objets[no2].delta_dir=mini(mini(fabs(x-y),fabs(x+180.-y)),fabs(x-y-180.));*/
	T_objets[no1].delta_dir=T_objets[no2].delta_dir=atan(2.f/pow(x*x+y*y,0.5f))*180/(float)PI;
//	cout<<" recalcul des parametres de dir. : "<<T_objets[no1].dir_axis<<" "<<T_objets[no1].delta_dir<<"\n";
	T_objets[no1].l_min=T_objets[no2].l_min=mini(T_objets[no1].l_min,T_objets[no2].l_min);
	T_objets[no1].l_max=T_objets[no2].l_max=maxi(T_objets[no1].l_max,T_objets[no2].l_max);
	T_objets[no1].c_min=T_objets[no2].c_min=mini(T_objets[no1].c_min,T_objets[no2].c_min);
	T_objets[no1].c_max=T_objets[no2].c_max=maxi(T_objets[no1].c_max,T_objets[no2].c_max);
	z=(T_objets[no1].l_baryc*T_objets[no1].surf+T_objets[no2].l_baryc*T_objets[no2].surf)/s;
	T_objets[no1].l_baryc=T_objets[no2].l_baryc=z;
	z=(T_objets[no1].c_baryc*T_objets[no1].surf+T_objets[no2].c_baryc*T_objets[no2].surf)/s;
	T_objets[no1].c_baryc=T_objets[no2].c_baryc=z;
	T_objets[no1].perim=T_objets[no2].perim=T_objets[no1].perim+T_objets[no2].perim;
	T_objets[no1].surf=T_objets[no2].surf=s;
	for (k=0; k<2; k++) {
		T_objets[no1].nb_cc[k]=T_objets[no2].nb_cc[k]=T_objets[no1].nb_cc[k]+T_objets[no2].nb_cc[k];
		T_objets[no1].nb_trou[k]=T_objets[no2].nb_trou[k]=T_objets[no1].nb_trou[k]+T_objets[no2].nb_trou[k];
		T_objets[no1].nbEuler[k]=T_objets[no2].nbEuler[k]=T_objets[no1].nbEuler[k]+T_objets[no2].nbEuler[k];
	}
	if (T_objets[no1].ratio_MajMinA>=T_objets[no2].ratio_MajMinA) {
		T_objets[no2].ratio_MajMinA=T_objets[no1].ratio_MajMinA;
		T_objets[no2].err_ratio=T_objets[no1].err_ratio;
	} else {
		T_objets[no1].ratio_MajMinA=T_objets[no2].ratio_MajMinA;
		T_objets[no1].err_ratio=T_objets[no2].err_ratio;
	}
} 

void imaobjets::regroupe_alignements (bool **Talign, unsigned int **Tdist, float distVmx, float distVmn, float fact_aplati, 
									  unsigned int obj0) {
	const float a_dist=(distVmx-distVmn)/nblig/(1.f-fact_aplati), b_dist=(distVmn-fact_aplati*distVmx)/(1.f-fact_aplati),
				err_a_min=7.0f, delta_dir_max=30.f, rapSurfMin=0.1f;
	const unsigned int itmax=10;
	bool permut, icontinue;
	int i,j,k,l,m,n,no1,no2,no3,nfus,nobj2=nobj,it=0;
	float distV, a_k, da_k, a_n, da_n, x, y, a/*, b*/, r_s;
	if (nobj>255) cout<<"\n Pb : trop d'objets : "<<nobj<<" pour cette methode !!!!!!!!!!!!!!\n\n";
	else {
		BYTE *Teq=new BYTE[nobj], *Tord=new BYTE[nobj], no;
		bool **Tregroup=new bool*[nobj], **T_actua=new bool*[nobj], ichg; 
		for (k=0; k<nobj; k++) {
			Tregroup[k]=new bool[nobj]; 
			T_actua[k]=new bool[nobj]; 
			for (n=0; n<nobj; n++) T_actua[k][n]=0;
			Teq[k]=k; Tord[k]=k;
		}
		no=Tord[nobj-1];
		Tord[nobj-1]=Tord[obj0];
		Tord[obj0]=no;
		do {
			permut=0;
			for (i=0; i<nobj-2; i++)
				if (T_objets[(int)Tord[i]].surf<T_objets[(int)Tord[i+1]].surf) {
					no=(int)Tord[i];
					Tord[i]=Tord[i+1];
					Tord[i+1]=(BYTE)no;
					permut=1;
				}
		} while (permut);
		sauve_Ima("./imaobjets_ini.dat");

		for (k=0; k<nobj; k++)
			for (n=0; n<nobj; n++) Talign[k][n]=0;
		do {
			nfus=0;
			for (k=0; k<nobj; k++) {
				no1=(int)Tord[k];
				for (n=0; n<nobj; n++) Tregroup[no1][n]=0;
				if (no1!=obj0 && (int)Teq[no1]==no1) {
					distV=maxi(a_dist*T_objets[no1].l_baryc+b_dist,distVmn); 
//					cout<<"\n distance considere autour de l'objet "<<no1<<" = "<<distV<<"\n";
//					if (T_objets[no1].i_dir_axis && T_objets[no1].delta_dir<=45) {
					if (T_objets[no1].i_dir_axis && T_objets[no1].delta_dir<=delta_dir_max) {
						cout<<"\n recherche des objets alignes sur l'objet "<<no1<<" de dir."<<T_objets[no1].dir_axis<<" a "<<T_objets[no1].delta_dir<<" deg. pres\n";
						a_k=T_objets[no1].dir_axis; if (a_k<0) {cout<<" pb : dir axis de "<<no1<<" non dans [0,PI[ : "<<a_k<<"\n"; T_objets[no1].affiche(); a_k+=180.;}
						da_k=maxi(err_a_min,T_objets[no1].delta_dir/2.f);
						for (n=0; n<nobj; n++) {
							no2=(int)Tord[n];
							if (no2!=obj0 && no2!=no1 && (int)Teq[no2]==no2) {
								y=T_objets[no1].l_baryc-T_objets[no2].l_baryc;
								x=T_objets[no1].c_baryc-T_objets[no2].c_baryc;
								a=atan2(-y,x)/(float)PI*180.f; if (a<0) a+=180.f; if (fabs(a-180)<1.e-3) a=0;
								a_n=T_objets[no2].dir_axis; if (a_n<0) {cout<<" pb : dir axis de "<<no2<<" non dans [0,PI[ : "<<a_n<<"\n"; T_objets[no2].affiche(); a_n+=180.;}
								da_n=maxi(err_a_min,T_objets[no2].delta_dir/2.f);
								if ((fabs(a_k-a)<da_k || fabs(a_k+180-a)<da_k || fabs(a_k-180-a)<da_k) && 
									(fabs(a_n-a)<da_n || fabs(a_n+180-a)<da_n || fabs(a_n-180-a)<da_n)) {
//								if (a>=a_k-da_k && a<=a_k+da_k && a>=a_n-da_n && a<=a_n+da_n) {
									Talign[no1][no2]=Talign[no2][no1]=1; 
//									if ((da_n<=30. && Tdist[no1][no2]<5*distV) || Tdist[no1][no2]<distV) {
									if ((da_n<=delta_dir_max && Tdist[no1][no2]<5*distV) || Tdist[no1][no2]<distV) {
										Tregroup[no1][no2]=Tregroup[no2][no1]=1; 
										cout<<" a = "<<a<<" objet "<<no2<<" aligne de dir."<<a_n<<" a "<<da_n<<" deg. pres et a dist "<<Tdist[no1][no2]<<"\n";
										if ((int)Teq[no1]!=no2) Teq[no2]=(BYTE)no1;
//										nobj2--; nfus++;
//										regroupe_2objets (no1,no2,Tdist);
									} //else cout<<" objet "<<no2<<" aligne mais non regroupe car ? "<<da_n<<" "<<Tdist[no1][no2]<<" "<<distV<<"\n";
								} //else cout<<" objet "<<no2<<" non aligne : a = "<<a<<" a_n = "<<a_n<<" deg. a "<<da_n<<" pres\n";
							} //else cout<<" objet "<<no2<<" non examine car no2 = obj0 ("<<obj0<<") ou no1 ("<<no1<<") ou !="<<(int)Teq[no2]<<"\n";
						}
					} else {
						for (n=0; n<nobj; n++) {
							no2=(int)Tord[n];
							if (no2!=obj0 && no2!=no1 && (int)Teq[no2]==no2 && Tdist[no1][no2]<distV && !Tregroup[no1][no2]) {
								y=T_objets[no1].l_baryc-T_objets[no2].l_baryc;
								x=T_objets[no1].c_baryc-T_objets[no2].c_baryc;
								a=atan2(-y,x)/(float)PI*180.f; if (a<0) a+=180.f; if (fabs(a-180)<1.e-3) a=0;
								for (m=0; m<nobj; m++) {
									no3=(int)Tord[m];
//									if (no3!=obj0 && no3!=no1 && no3!=no2 && (int)Teq[no3]==no3) {
									if (no3!=obj0 && no3!=no1 && no3!=no2 && (int)Teq[no3]==no3 && mini(Tdist[no1][no3],Tdist[no2][no3])<distV) {
										y=T_objets[no1].l_baryc-T_objets[no3].l_baryc;
										x=T_objets[no1].c_baryc-T_objets[no3].c_baryc;
										a_k=atan2(-y,x)/(float)PI*180.f; if (a_k<0) a_k+=180.f; if (fabs(a_k-180)<1.e-3) a_k=0;
										y=T_objets[no2].l_baryc-T_objets[no3].l_baryc;
										x=T_objets[no2].c_baryc-T_objets[no3].c_baryc;
										a_n=atan2(-y,x)/(float)PI*180.f; if (a_n<0) a_n+=180.f; if (fabs(a_n-180)<1.e-3) a_n=0;
										if ((fabs(a_k-a)<err_a_min || fabs(a_k+180-a)<err_a_min || fabs(a_k-180-a)<err_a_min) && 
											(fabs(a_n-a)<err_a_min || fabs(a_n+180-a)<err_a_min || fabs(a_n-180-a)<err_a_min) && 
											(fabs(a_n-a_k)<err_a_min || fabs(a_n+180-a_k)<err_a_min || fabs(a_n-180-a_k)<err_a_min)) {
											cout<<" objets "<<no1<<", "<<no2<<" et "<<no3<<" proches : dist="<<mini(Tdist[no1][no3],Tdist[no2][no3]);
											cout<<", et alignes "<<" de dir."<<a<<" = "<<a_k<<" = "<<a_n<<" deg.\n";
											r_s=mini(mini(T_objets[no1].surf,T_objets[no2].surf)/maxi(T_objets[no1].surf,T_objets[no2].surf),
													 mini(mini(T_objets[no2].surf,T_objets[no3].surf)/maxi(T_objets[no2].surf,T_objets[no3].surf),
													 mini(T_objets[no1].surf,T_objets[no3].surf)/maxi(T_objets[no1].surf,T_objets[no3].surf)));
											if (r_s>rapSurfMin) {
												Talign[no1][no2]=Talign[no2][no1]=1; 
												Talign[no1][no3]=Talign[no3][no1]=1; 
												Talign[no2][no3]=Talign[no3][no2]=1; 
												Tregroup[no1][no2]=Tregroup[no2][no1]=1; 
												Tregroup[no1][no3]=Tregroup[no3][no1]=1; 
												Tregroup[no2][no3]=Tregroup[no3][no2]=1; 
												if ((int)Teq[no1]!=no2) Teq[no2]=(BYTE)no1;
												if ((int)Teq[no1]!=no3) Teq[no3]=(BYTE)no1;
											} else cout<<" mais rapport des surface trop grand : "<<r_s<<" < "<<rapSurfMin<<"\n";
										}
									}
								}

							}
						}						
					}
				}
				for (n=0; n<nobj; n++) {
					no2=(int)Tord[n];
					if (Tregroup[no1][no2]) {
						regroupe_2objets (no1,no2,Tdist);
						nobj2--; nfus++; //cout<<" regroupe objets "<<no1<<" et "<<no2<<" => maintenant "<<nobj2<<" objets\n";
//						if (no1>=30 || no2>=30) {cout<<"\n objets "<<no1<<" et "<<no2<<"\n"; T_objets[no1].affiche(); T_objets[no2].affiche(); char aa; cin>>aa;} 
						for (i=0; i<nobj; i++) {
//							cout<<" "<<T_actua[no1][i];
							if (T_actua[no1][i]) {
								T_objets[i]=T_objets[no1]; 
//								cout<<" actualisation objet "<<i<<" par objet "<<no1<<"\n"; /*T_objets[i].affiche();*/
								for (j=0; j<nobj; j++) {
									if (T_actua[i][j]) {T_objets[j]=T_objets[i];
//									cout<<" actualisation objet "<<j<<" par objet "<<i<<" = objet "<<no1<<"\n"; /*T_objets[j].affiche();*/
									T_actua[no1][j]=1;}
								}
							}
							if (T_actua[no2][i]) {
								T_objets[i]=T_objets[no2]; 
//								cout<<" actualisation objet "<<i<<" par objet "<<no2<<"\n"; /*T_objets[i].affiche();*/
								for (j=0; j<nobj; j++) {
									if (T_actua[i][j]) {T_objets[j]=T_objets[i];
//									cout<<" actualisation objet "<<j<<" par objet "<<i<<" = objet "<<no2<<"\n"; /*T_objets[j].affiche();*/
									T_actua[no2][j]=1;}
								}
							}
						}
//						cout<<"\n";
						T_actua[no1][no2]=1; //cout<<" objet "<<no2<<" sera a actualiser si autre fusion\n";
						for (i=0; i<nobj; i++)
							if ( (Tregroup[no1][i] || Tregroup[i][no1]) && (Tregroup[no2][i] || Tregroup[i][no2]) )
								Tregroup[no2][i]=Tregroup[i][no2]=0;
					}
				}
			}
			if (nfus>0) {
				cout<<" iteration "<<++it<<" : apres fusion des objets alignes : "<<nobj2<<" objets reellement dans l'image\n";
				icontinue=1; 
			} else 
				icontinue=0;
		} while (icontinue && it<itmax);
		cout<<"\n"; for (k=0; k<nobj; k++) cout<<k<<" "<<(int)Teq[k]<<"\n";
		bool *T_ricochet=new bool[nobj];
		do {
			ichg=0;
			for (i=0; i<nobj; i++) {
				j=(int)Teq[i];
				n=(int)Teq[j];
				if (j!=0 && n!=j) {
					m=mini(j,n);
					cout<<" ricochet dans la table d'equivalence : Teq["<<i<<"] = "<<j<<" et Teq["<<j<<"] = "<<n<<"\n";
					k=(int)Teq[j];
					if (k!=j) {
						for (l=0; l<nobj; l++) T_ricochet[l]=0;
						T_ricochet[i]=T_ricochet[j]=T_ricochet[k]=1;
						while (k!=(int)Teq[k] && T_ricochet[(int)Teq[k]]!=1) {
							cout<<" ricochet se poursuit : Teq["<<k<<"] = "<<k<<"\n";
							k=(int)Teq[k]; T_ricochet[k]=1;
						}
					}
					if (k!=j) {
						cout<<" chg en "<<k<<" : "<<(int)Teq[k]<<" devient "<<m<<"\n"; 
						Teq[k]=(BYTE)m;
					}
					if (m==j) {
						cout<<" chg en "<<n<<" : "<<(int)Teq[n]<<" devient "<<m<<"\n"; 
						Teq[n]=BYTE(m); ichg=1;
						for (k=0; k<nobj; k++) 
							if ((int)Teq[k]==n) {
								cout<<" chg en "<<k<<" : "<<(int)Teq[k]<<" devient "<<m<<"\n"; 
								Teq[k]=BYTE(m); ichg=1; }
					} else 
						for (k=0; k<nobj; k++) 
							if ((int)Teq[k]==j) {
								cout<<" chg en "<<k<<" : "<<(int)Teq[k]<<" devient "<<m<<"\n"; 
								Teq[k]=BYTE(m); ichg=1; }
				}
			}
		} while (ichg);
		if (T_ricochet!=NULL) delete[] T_ricochet;
		cout<<"\n"; for (k=0; k<nobj; k++) cout<<k<<" "<<(int)Teq[k]<<"\n";
		nobj2=0;
		for (k=0; k<nobj; k++)
			if ((int)Teq[k]>=nobj2) {
				if ((int)Teq[k]>nobj2) {
					no=Teq[k];
					for (i=0; i<nobj; i++)
						if (Teq[i]==no) Teq[i]=(BYTE)nobj2;
				}
				nobj2++;
			}
		cout<<" apres re-arrangement de la table d'equivalence : "<<nobj2<<" objets\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				(*this)(i,j)=Teq[(*this)(i,j)];
		for (k=0; k<nobj; k++)
			if ((int)Teq[k]!=k)
				T_objets[(int)Teq[k]]=T_objets[k];
		for (k=0; k<nobj; k++)
			for (n=k; n<nobj; n++) {
				Tdist[(int)Teq[k]][(int)Teq[n]]=Tdist[k][n];
				Talign[(int)Teq[k]][(int)Teq[n]]=Talign[k][n];
			}
		for (k=0; k<nobj; k++)
			for (n=0; n<k; n++) {
				Tdist[k][n]=Tdist[n][k];
				Talign[k][n]=Talign[n][k];
			}
		for (k=0; k<nobj; k++) {
			if (Tregroup[k]!=NULL) delete[] Tregroup[k];
			if (T_actua[k]!=NULL) delete[] T_actua[k];
		}
		if (Tregroup!=NULL) delete[] Tregroup; 
		if (T_actua!=NULL) delete[] T_actua;
		nobj=nobj2;
		cout<<" fin a l'iteration "<<it<<" : apres fusion des objets alignes : "<<nobj<<" objets reellement dans l'image\n";
		if (Teq!=NULL) delete[] Teq;
		if (Tord!=NULL) delete[] Tord;
//		affiche();
/*		if (T_objets!=NULL) delete[] T_objets;
		T_objets=new objet[nobj]; 
		basic_param ();
		basic_topologie (obj0);
		second_param (obj0);
		cout<<"\n"; affiche();*/
		sauve_Ima("./imaobjets_fus.dat");

/*		imadata<unsigned int> ima_align(nblig,nbcol);
		imabin imab_align(nblig,nbcol);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) ima_align(i,j)=0;
		for (k=1; k<nobj; k++)
			for (n=1; n<nobj; n++) 
				if (k!=n && Talign[k][n]) {
					cout<<"alignement objets "<<k<<" et "<<n<<" poids "<<(T_objets[k].surf+T_objets[n].surf)/2<<"\n";
					for (i=0; i<nblig; i++)
						for (j=0; j<nbcol; j++) imab_align(i,j)=0;
					y=T_objets[k].l_baryc-T_objets[n].l_baryc;
					x=T_objets[k].c_baryc-T_objets[n].c_baryc;
					if (x!=0) {
						a=-y/x;
						b=(nblig-1-T_objets[k].l_baryc)-a*T_objets[k].c_baryc;
						cout<<" dir "<<T_objets[k].dir_axis<<" "<<T_objets[n].dir_axis<<" Eq.: y="<<a<<"*x+"<<b<<"\n";
						for (j=0; j<nbcol; j++) {
							i=(int)(nblig-1-(a*j+b));
							if (i>=0 && i<nblig && !imab_align(i,j)) {ima_align(i,j)+=T_objets[k].surf/2; imab_align(i,j)=1;}
						}
						if (a!=0) 
							for (i=0; i<nblig; i++) {
								j=(int)((nblig-1-i-b)/a);
								if (j>=0 && j<nbcol && !imab_align(i,j)) {ima_align(i,j)+=T_objets[k].surf/2; imab_align(i,j)=1;}
							}
					} else {
						j=(int)(T_objets[k].c_baryc);
						for (i=0; i<nblig; i++) ima_align(i,j)+=T_objets[k].surf/2;
					}
				}
		ima_align.sauve_ImaPGM("ima_align.pgm");
*/
	}
/*	imadata<BYTE> imaseg(nblig,nbcol);
	for (k=0; k<nobj; k++) 
		if (k!=obj0)
			if (T_objets[k].i_dir_axis  && T_objets[k].delta_dir<=45) {
				cout<<" recherche des objets alignes sur l'objet "<<k<<"\n";
				for (m=-d_a; m<=d_a; m+=d_a)
					if ((T_objets[k].dir_axis+m)!=90) {
						a=tan((T_objets[k].dir_axis+m)/180.*PI);
						for (i=-1; i<=1; i++)
							if ((*this)((int)T_objets[k].l_baryc,(int)T_objets[k].c_baryc+i)==k) {
								b=(nblig-1-T_objets[k].l_baryc)-a*(T_objets[k].c_baryc+i);
								for (j=0; j<nbcol; j++) {
									y=nblig-1-(a*j+b);
//									if (k==35 && y>115 && y<140) cout<<"("<<(int)y<<","<<j<<") ";
									if (y>=0 && y<nblig) {
										imaseg((int)y,j)=k;
										n=(*this)((int)y,j);
										if (n!=obj0 && !Talign[k][n]) {Talign[k][n]=Talign[n][k]=1; cout<<" objet "<<n<<" aligne\n";}
									}
								}
								if (a!=0)
									for (j=0; j<nblig; j++) {
										x=((nblig-1-j)-b)/a;
//										if (k==35 && j>115 && j<140) cout<<"* "<<m<<" "<<i<<"("<<j<<","<<(int)x<<") ";
										if (x>=0 && x<nbcol) {
											imaseg(j,(int)x)=k;
											n=(*this)(j,(int)x);
											if (n!=obj0 && !Talign[k][n]) {Talign[k][n]=Talign[n][k]=1; cout<<" objet "<<n<<" aligne\n";}
										}
									}
							}
					}
			}
	for (k=1; k<nobj; k++)
		for (n=1; n<nobj; n++) 
			if (k!=n && Talign[k][n]) {
				x=T_objets[k].dir_axis-T_objets[k].delta_dir;
				y=T_objets[n].dir_axis+T_objets[n].delta_dir;
				if (x>y) Talign[k][n]=Talign[n][k]=0;
				else {
					x=T_objets[k].dir_axis+T_objets[k].delta_dir;
					y=T_objets[n].dir_axis-T_objets[n].delta_dir;
					if (x<y) Talign[k][n]=Talign[n][k]=0;
					else cout<<" alignement objets "<<k<<" et "<<n<<"\n";
				}
			}
	imaseg.sauve_ImaPGM("imaprolongeseg.pgm");*/
}

/* ----------------------------------------------
Detection de disques a partir d'une image d'objets
renvoit une image binaire contenant les disques detectes

 seuil1, seuil2: peuvent etre utilisés dans la decision
 pour le seuillage des parametres de forme
-----------------------------------------------*/
imabin imaobjets::detection_disques(float seuil1, float seuil2) {
	int i,j,k;
	imabin out(nblig,nbcol);

	cout<<"DETECTION DISQUES dans une image a "<<nobj<<" objets\n";
	affiche();
//	aire[k]=T_objets[k].surf;
//	ic[k]=T_objets[k].l_baryc;
//	jc[k]=T_objets[k].c_baryc;
//	perim[k]=T_objets[k].perim;

//-----         METHOD 1 		-------
// calcul du perimetre et de la signature des contours
// le calcul du perimetre se fait a partir du gradient morphologique binaire
//--------------------------------------
	imabin objet(nblig,nbcol);				//image binaire de l'objet
	imabin contour(nblig,nbcol);			//image de contour

	double *r_min=new double[nobj];			//tab de distance min par rapport au centre
	double *r_max=new double[nobj];			//tab de distance max par rapport au centre
	double r;
	int diag=(int)pow(pow((float)nbcol,2)+pow((float)nblig,2),0.5f);
	cout<<" diag = "<<diag<<"\n";
	for (k=0; k<nobj; k++) {
		T_objets[k].perim=0;
		r_min[k]=(float)diag;
		r_max[k]=0.f;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if ((*this)(i,j)==k) objet(i,j)=1;
				else objet(i,j)=0;
		contour=objet-objet.erode(4,1);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (contour(i,j)) {
					T_objets[k].perim+=1;
					r=pow(pow((double)i-T_objets[k].l_baryc,2)+pow((double)j-T_objets[k].c_baryc,2),0.5);
					if (r<r_min[k]) r_min[k]=r;
					if (r>r_max[k]) r_max[k]=r;
				}
		cout<<" objet "<<k<<" : l_baryc = "<<T_objets[k].l_baryc<<", c_baryc = "<<T_objets[k].c_baryc<<"\n";
		cout<<" objet "<<k<<" : rmin = "<<r_min[k]<<", rmax = "<<r_max[k]<<"\n";
	}
//-----         METHOD 2 		-----
//------------ calcul du rapport des moment m20/m02
//--------------------------------------
	float *ratio=new float[nobj];		//tab de rapport m20/m02
	float *m20=new float[nobj];
	float *m02=new float[nobj];
	for (k=0; k<nobj; k++) m20[k]=m02[k]=0.0f;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			k=(*this)(i,j);
			if (k>0) {
				m20[k]+=(pow((i-T_objets[k].l_baryc),2)*pow((j-T_objets[k].c_baryc),0));
				m02[k]+=(pow((i-T_objets[k].l_baryc),0)*pow((j-T_objets[k].c_baryc),2));
			}
		}
	for (k=0; k<nobj; k++) {
		if (m02[k]!=0) ratio[k]=m20[k]/m02[k];
		else ratio[k]=m20[k];
	}
	if (m20!=NULL) delete[] m20;
	if (m02!=NULL) delete[] m02;

//--------------Affichage ------------------------
	for (k=0; k<nobj; k++) {
		cout<<" objet "<<k<<" : aire = "<<T_objets[k].surf<<", barycentre en ("<<T_objets[k].l_baryc<<", "<<T_objets[k].c_baryc<<")\n";
		cout<<"\tSignature : "<<r_min[k]/r_max[k]<<"\tRapport : "<<ratio[k]<<  "\n";
	}
	if (r_min!=NULL) delete[] r_min;
	if (r_max!=NULL) delete[] r_max;

//--------------Decision ------------------------
	int *disque=new int[nobj]; // tab mis a 1 si l'objet k est un disque
	cout << " Disques extraits par le rapport de moments : \n";
	for (k=0; k<nobj; k++) disque[k]=0;
	for (k=0; k<nobj; k++)
		if((ratio[k]<seuil2) && (ratio[k]>seuil1))	{
			cout<< "Objet "<<k<< "\n";
			disque[k]=1;
		}
	if (ratio!=NULL) delete[] ratio;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++){
			if (disque[(*this)(i,j)]==0) out(i,j)=0;
			else out(i,j)=1;
		}
	if (disque!=NULL) delete[] disque;
	return out;
}

imalabels imaobjets::clas_marquages_route (unsigned int obj0, float _fact_aplati, bool iaf) {
	const int itermax=11, nbcl=7; /* classes :	1='passage pieton'
												2='ligne continue de separation de voie'
												3='ligne discontinue de separation de voie'
												4='ligne gauche de bord de voie'
												5='ligne droite de bord de voie'
												6='ligne centrale de separation de voie'												 
												7='autre marquage'												 
												0='rejet' */
	const int n_lab_cl0=5, n_lab_lig=5, n_lab_piet=1, n_lab_ldis=4;
	BYTE lab_lignes[n_lab_lig]={2,3,4,5,6}, lab_pieton[n_lab_piet]={1}, lab_notlig[nbcl+1-n_lab_lig]={0,1,7},
		 lab_notpiet[nbcl+1-n_lab_piet]={0,2,3,4,5,6,7}, lab_ligdis[n_lab_ldis]={3,4,5,6}, lab_cl0[n_lab_cl0]={0,1,2,3,7};
	cout<<" debut clas_marquages_route\n";
	float fact_aplati=_fact_aplati; 
//	const float ratioligne=1./0.15, fact_aplati=1., alpha1=10./(1./0.1);
	const float ratioligne=3.33f, epsil=0.1f*ratioligne, alpha1=1.f;
//	const float ratioUzebra=1.f/0.5f, ratioLzebra=1.f/0.25f, alpha2=2.0f/(ratioUzebra-ratioLzebra);
	const float ratioUzebra=8.0f, ratioLzebra=3.0f, alpha2=1.0f, coefdistpp=2.0f, coefrsurpp=0.5f, coefdenspp=0.33f;
	const float betaV=2.f, distVmx=10.f, distVmn=1.f, rmse_max=1000., betaL=betaV/5;
//	const float d_a=1.f;
	const int np_min=5;
	int i,j,k,l,/*m,*/n,nch,iter=0;
	float distV, a_dist, b_dist;
/* calculs preliminaires : table des distances entre objets et regroupement d'objets */
	bool i_lig, iLR;
	unsigned int **Tdist=new unsigned int*[nobj], *d_min=new unsigned int[nobj], d;
	for (k=0; k<nobj; k++) Tdist[k]=new unsigned int[nobj];
	imabin imab(nblig,nbcol);
	imadata<float> imadist;
	for (k=0; k<nobj; k++) {
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if ((*this)(i,j)==k) imab(i,j)=1;
				else imab(i,j)=0;
		imadist=imab.Tr_dist();
		for (n=0; n<nobj; n++) d_min[n]=UINT_MAX;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				n=(*this)(i,j); d=(unsigned int)imadist(i,j);
				if (d<d_min[n]) d_min[n]=d;
			}
		for (n=0; n<nobj; n++) Tdist[k][n]=Tdist[n][k]=d_min[n];
	}
	bool **Talign=new bool*[nobj]; 
	for (k=0; k<nobj; k++)
		Talign[k]=new bool[nobj];
	cout<<"\n avant regroupement des objets alignes : "<<nobj<<" objets\n"; affiche();
	regroupe_alignements (Talign,Tdist,5*distVmx,5*distVmn,fact_aplati,obj0);
	cout<<"\n apres regroupement des objets alignes : "<<nobj<<" objets\n"; affiche();
//	regroupe_alignements (Talign,Tdist,5*distVmx,5*distVmn,fact_aplati,obj0);
//	cout<<"\n apres 2eme regroupement des objets alignes : "<<nobj<<" objets\n"; affiche();
/* densite des lignes en particulier non continues */
	int lmin, lmax, cmin, cmax;
	bool i_marq;
	float *T_dens=new float[nobj];
	for (k=0; k<nobj; k++) { 
		if (k!=obj0) {
			if (T_objets[k].i_cc && T_objets[k].nb_cc[0]==1) T_dens[k]=1.f;
			else {
				T_dens[k]=0.f;
				if (T_objets[k].dir_axis>=45 && T_objets[k].dir_axis<=135) {
					lmin=T_objets[k].l_min; lmax=T_objets[k].l_max; cmin=T_objets[k].c_min; cmax=T_objets[k].c_max;
					for (i=lmin; i<=lmax; i++) {
						i_marq=0; j=cmin;
						do {
							if ((*this)(i,j)==k) i_marq=1;
							j++;
						} while (!i_marq && j<=cmax);
						T_dens[k]+=(float)i_marq;
					}
					T_dens[k]/=(float)(lmax-lmin+1);
				} else {
					lmin=T_objets[k].l_min; lmax=T_objets[k].l_max; cmin=T_objets[k].c_min; cmax=T_objets[k].c_max;
					for (j=cmin; j<=cmax; j++) {
						i_marq=0; i=lmin;
						do {
							if ((*this)(i,j)==k) i_marq=1;
							i++;
						} while (!i_marq && i<=lmax);
						T_dens[k]+=(float)i_marq;
					}
					T_dens[k]/=(float)(cmax-cmin+1);
				}
				cout<<" objet "<<k<<" #composantes connexes = "<<T_objets[k].nb_cc[0]<<" densite des pointilles = "<<T_dens[k]<<"\n";
			}
		}
	}
/* rmse sur la fct distance entre les deux objets */
	float **dvsl=new float*[nobj], **Txy=new float*[nblig], a_dvsl, b_dvsl, rmse, **T_rmse=new float*[nobj];
	int **nvsl=new int*[nobj], **T_nobjLR=new int*[nobj];
	for (k=0; k<nobj; k++) {dvsl[k]=new float[nblig]; nvsl[k]=new int[nblig]; 
							T_rmse[k]=new float[nobj]; T_nobjLR[k]=new int[2];}
	for (i=0; i<nblig; i++) Txy[i]=new float[2];
	for (k=0; k<nobj; k++) { 
		if (k!=obj0) {
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if ((*this)(i,j)==k) imab(i,j)=1;
					else imab(i,j)=0;
			imadist=imab.Tr_dist();
			for (i=0; i<nblig; i++) {
//				for (n=0; n<nobj; n++) dvsl[n][i]=-1;
				for (n=0; n<nobj; n++) {dvsl[n][i]=0.f; nvsl[n][i]=0;}
				i_lig=0;
				for (j=0; j<nbcol; j++)
					if (!i_lig && imab(i,j)) i_lig=1;
				if (i_lig) {
					for (j=0; j<nbcol; j++) {
						n=(*this)(i,j);
//						if (n!=obj0 && n!=k && (imadist(i,j)<dvsl[n][i] || dvsl[n][i]<0) ) dvsl[n][i]=imadist(i,j);
						if (n!=obj0 && n!=k) {dvsl[n][i]+=imadist(i,j); nvsl[n][i]++;}
					}
				}
				for (n=0; n<nobj; n++) 
					if (nvsl[n][i]>0) dvsl[n][i]/=nvsl[n][i];
					else dvsl[n][i]=-1.f;
			}
			T_nobjLR[k][0]=T_nobjLR[k][1]=0;
			for (n=0; n<nobj; n++) {
				j=0;
				for (i=0; i<nblig; i++)
					if (dvsl[n][i]>=0) {Txy[j][1]=dvsl[n][i]; Txy[j][0]=(float)i; j++;}
				rmse=FLT_MAX;
				if (j>np_min) {
					rmse=regres_lin(Txy,j,a_dvsl,b_dvsl);
//					if (rmse<rmse_max) {
						cout<<" objet "<<k<<" mis en relation avec objet "<<n<<" => a = "<<a_dvsl<<" & rmse = "<<rmse<<"\n";
						if (T_objets[k].c_baryc>T_objets[n].c_baryc) T_nobjLR[k][0]++;
						else T_nobjLR[k][1]++;
//					}
				}
				T_rmse[k][n]=rmse;
			}
		}
	}
	for (k=0; k<nobj; k++) {
		if (dvsl[k]!=NULL) delete[] dvsl[k];
		if (nvsl[k]!=NULL) delete[] nvsl[k];
	}
	for (i=0; i<nblig; i++) if (Txy[i]!=NULL) delete[] Txy[i];
	if (dvsl!=NULL) delete[] dvsl;
	if (nvsl!=NULL) delete[] nvsl;
	if (Txy!=NULL) delete[] Txy;
/* calcul du terme d'attache au donnees : parametre de forme = ratio petit axe / grand axe */
	double **U0=new double*[nobj], **U=new double*[nobj], xx, yy, zz, a_r, b_r, r, a_debordHoriz, Umin, 
		   a_rZu, b_rZu, a_rZl, b_rZl, rZl=ratioLzebra, rZu=ratioUzebra;
	for (k=0; k<nobj; k++) {U0[k]=new double[nbcl+1]; U[k]=new double[nbcl+1];}
	a_r=(1.-1./(ratioligne-epsil))/nblig/(1.-fact_aplati);
//	b_r=1.-(1.-1./(ratioligne-epsil))/(1.-fact_aplati);
	b_r=(-fact_aplati+1./(ratioligne-epsil))/(1.-fact_aplati);
	a_debordHoriz=10./nblig/fact_aplati;
	a_rZu=(1.-1./(ratioUzebra-epsil))/nblig/(1.-fact_aplati); b_rZu=(-fact_aplati+1./(ratioUzebra-epsil))/(1.-fact_aplati);
	a_rZl=(1.-1./(ratioLzebra-epsil))/nblig/(1.-fact_aplati); b_rZl=(-fact_aplati+1./(ratioLzebra-epsil))/(1.-fact_aplati);
	for (k=0; k<nobj; k++) {
		for (i=0; i<=nbcl; i++) U0[k][i]=0.;
		if (T_objets[k].i_ratio && T_objets[k].ratio_MajMinA!=0.) {
			cout<<" \n objet "<<k<<" initialement : U = "; for (i=0; i<=nbcl; i++) cout<<U0[k][i]<<" "; cout<<"\n";
			xx=1./T_objets[k].ratio_MajMinA;
			r=maxi(ratioligne*(a_r*T_objets[k].l_baryc+b_r),(double)ratioligne/3);
			yy=log(1.+exp(-alpha1*(xx-r))); 
			for (j=0; j<n_lab_lig; j++) U0[k][(int)lab_lignes[j]]+=yy; 
			yy=log(1.+exp(+alpha1*(xx-r)));
			for (j=0; j<=nbcl-n_lab_lig; j++) U0[k][(int)lab_notlig[j]]+=yy; 
			cout<<" apres clas ligne : U = "; for (i=0; i<=nbcl; i++) cout<<U0[k][i]<<" "; cout<<"\n";
			yy=-log(1./(1.+exp(-alpha2*(xx-ratioLzebra)))+1./(1.+exp(alpha2*(xx-ratioUzebra)))-1.);
//			rZu=maxi(ratioUzebra*(a_rZu*T_objets[k].l_baryc+b_rZu),(double)ratioUzebra/3);
//			rZl=maxi(ratioLzebra*(a_rZl*T_objets[k].l_baryc+b_rZl),(double)ratioLzebra/3);
//			yy=-log(1./(1.+exp(-alpha2*(xx-rZl)))+1./(1.+exp(alpha2*(xx-rZu)))-1.);
//			cout<<" $$$$$ "<<ratioUzebra<<" "<<ratioLzebra<<" devient "<<rZu<<" "<<rZl<<" a la ligne "<<T_objets[k].l_baryc<<" "<<yy<<"\n";
			cout<<" obj."<<k<<" xx "<<xx<<" $$$$$ energie passage pieton "<<yy;
			for (j=0; j<n_lab_piet; j++) U0[k][(int)lab_pieton[j]]+=yy; 
			yy=-log(2-1./(1.+exp(-alpha2*(xx-ratioLzebra)))-1./(1.+exp(alpha2*(xx-ratioUzebra))));
//			yy=-log(2-1./(1.+exp(-alpha2*(xx-rZl)))-1./(1.+exp(alpha2*(xx-rZu))));
			cout<<" & not piet "<<yy<<" $$$$$ "<<exp(-yy)+1./(1.+exp(-alpha2*(xx-rZl)))+1./(1.+exp(alpha2*(xx-rZu)))-1.<<"\n";
			for (j=0; j<=nbcl-n_lab_piet; j++) U0[k][(int)lab_notpiet[j]]+=yy; 
			cout<<" apres clas pass piet : U = "; for (i=0; i<=nbcl; i++) cout<<U0[k][i]<<" "; cout<<"\n";
		} else {
			for (i=0; i<=nbcl; i++) U0[k][i]=EXP_1;
		}
		xx=a_debordHoriz*maxi((double)fact_aplati*nblig-T_objets[k].l_baryc,(double)0);
		yy=maxi((double)0,EXP_1*(1./(1.+exp(-xx))-0.5));
		U0[k][0]-=yy;
		cout<<" avant terme alignement : U = "; for (i=0; i<=nbcl; i++) cout<<U0[k][i]<<" "; cout<<"\n";
//		{char aa; cin>>aa;}
		if (k!=obj0) {
			xx=T_dens[k];
			cout<<" objet "<<k<<" #composantes = "<<T_objets[k].nb_cc[0]<<" densite des pointilles = "<<T_dens[k]<<"\n";
			if (xx<0.75) {U0[k][3]-=betaL; U0[k][6]-=betaL; cout<<" attache aux donnees clas. centre diminuee de "<<betaL<<"\n";}
			else {U0[k][2]-=betaL; U0[k][4]-=betaL; U0[k][5]-=betaL; cout<<" attache aux donnees clas. bords route diminuee de "<<betaL<<"\n";}
//			yy=0.5*log(1.+exp(-alpha1*(xx-0.75))); U0[k][4]-=yy; U0[k][5]-=yy;
//			yy=0.5*log(1.+exp(+alpha1*(xx-0.75))); U0[k][3]-=yy; U0[k][6]-=yy;
		}
		for (l=0; l<n_lab_ldis; l++) U0[k][(int)lab_ligdis[l]]-=(T_objets[k].nb_cc[1]-1)*betaV;
		cout<<" apres terme alignement : U = "; for (i=0; i<=nbcl; i++) cout<<U0[k][i]<<" "; cout<<"\n";
		iLR=0;
		for (n=0; n<nobj; n++) 
			if (n!=obj0 && n!=k && T_rmse[k][n]<rmse_max) {
				if (T_nobjLR[k][0]>0) {
					if (T_nobjLR[k][1]==0) {U0[k][5]-=betaL*2;}
					else iLR=1;
				} else
					if (T_nobjLR[k][0]==0) {U0[k][4]-=betaL*2;}
			}
		if (iLR) U0[k][6]-=betaL*2;
		cout<<" apres mise en correspondance (rmse) : U = "; for (i=0; i<=nbcl; i++) cout<<U0[k][i]<<" "; cout<<"\n";
	}
	for (k=0; k<nobj; k++) {
		l=0; Umin=U0[k][l];
		for (i=1; i<n_lab_cl0; i++)
			if (U0[k][lab_cl0[i]]<Umin) {l=i; Umin=U0[k][i];}
		T_objets[k].lab=l; 
		T_objets[k].i_lab=1;
	}
	if (iaf)
		for (k=0; k<nobj; k++) {
			cout<<"obj."<<setw(2)<<k<<" "; 
			for (i=0; i<=nbcl; i++) cout<<setw(6)<<setprecision(3)<<U0[k][i]<<" "; 
			cout<<" => l="<<(int)T_objets[k].lab<<"\n";
		}
	imalabels imalab(nblig,nbcol,2);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			imalab(i,j,0)=imalab(i,j,1)=T_objets[(*this)(i,j)].lab;

/*	for (k=0; k<nobj; k++) {
		for (n=0; n<nobj; n++) cout<<setprecision(3)<<setw(4)<<Tdist[k][n];
		cout<<"\n";
	}*/
	a_dist=(distVmx-distVmn)/(float)nblig/(1.f-fact_aplati);
	b_dist=(distVmn-fact_aplati*distVmx)/(1.f-fact_aplati);
//	cout<<distVmx<<" "<<distVmn<<" "<<nblig<<" a_dist = "<<a_dist<<" b_dist = "<<b_dist<<"\n";
	do { char aa; cin>>aa;
		for (k=0; k<nobj; k++)
			if (k!=obj0) {
				U[k][0]=U0[k][0];
				for (l=1; l<=nbcl; l++) U[k][l]=U0[k][l]+betaV/2;
			}
		for (k=0; k<nobj; k++)
			if (k!=obj0) {
				distV=maxi(a_dist*T_objets[k].l_baryc+b_dist,distVmn);
				cout<<"obj."<<k<<" avant vois. U = "; for (i=0; i<=nbcl; i++) cout<<U[k][i]<<" "; cout<<"\n";
				for (n=0; n<nobj; n++) 
					if (n!=k && Tdist[k][n]<distV) {
//					if (n!=obj0 && n!=k && Tdist[k][n]<distV) {
						l=T_objets[n].lab;
						if (l>0) U[k][l]-=betaV;
					}
				cout<<" apres vois. standard : U = "; for (i=0; i<=nbcl; i++) cout<<U[k][i]<<" "; cout<<"\n";
//				cout<<"obj."<<k<<" distV = "<<distV<<" U = "; for (i=0; i<=nbcl; i++) cout<<U[k][i]<<" "; cout<<"\n";

				if (T_objets[k].i_ratio && T_objets[k].i_surface && T_objets[k].i_boite_englobante &&
					(double)T_objets[k].surf/(T_objets[k].l_max-T_objets[k].l_min)/(T_objets[k].c_max-T_objets[k].c_min)>coefdenspp) {
					xx=coefdistpp*T_objets[k].lg_MinA;
					for (n=0; n<nobj; n++) 
						if (n!=obj0 && n!=k && Tdist[k][n]<xx && T_objets[n].i_surface && T_objets[n].i_boite_englobante) {
							yy=mini((double)T_objets[k].surf,(double)T_objets[n].surf)/maxi((double)T_objets[k].surf,(double)T_objets[n].surf);
							zz=(double)T_objets[n].surf/(T_objets[n].l_max-T_objets[n].l_min)/(T_objets[n].c_max-T_objets[n].c_min);
							if (yy>coefrsurpp && zz>coefdenspp) {
//							if (yy>coefrsurpp) {
								cout<<" lab pieton a priori OUI "<<k<<" "<<n<<" "<<yy<<" "<<coefrsurpp<<" "<<zz<<" "<<coefdenspp<<"\n";
								for (j=0; j<n_lab_piet; j++) {
									U[k][(int)lab_pieton[j]]-=betaV/1.;
									U[n][(int)lab_pieton[j]]-=betaV/1.;
								}
								for (j=0; j<=nbcl-n_lab_piet; j++) {
									U[k][(int)lab_notpiet[j]]+=betaV/1.;
									U[n][(int)lab_notpiet[j]]+=betaV/1.;
								}
							} else {
								cout<<" lab pieton a priori NON "<<k<<" "<<n<<" "<<yy<<" "<<coefrsurpp<<" "<<zz<<" "<<coefdenspp<<"\n";
								T_objets[k].affiche(); T_objets[n].affiche();
//								char aa; cin>>aa;
							}
						}
				}
				cout<<" apres vois. pass piet : U = "; for (i=0; i<=nbcl; i++) cout<<U[k][i]<<" "; cout<<"\n";
			}
		nch=0;
		for (k=0; k<nobj; k++)
			if (k!=obj0) {
				l=0; Umin=U[k][l];
				for (i=1; i<=nbcl; i++)
					if (U[k][i]<Umin) {l=i; Umin=U[k][i];}
				if (l!=T_objets[k].lab) {
					nch++; //cout<<" changement de "<<(int)T_objets[k].lab<<" vers "<<l<<"\n";
					T_objets[k].lab=l;
				}
			}
		if (iaf)
			for (k=0; k<nobj; k++) {
				cout<<"obj."<<k<<" "; 
				for (i=0; i<=nbcl; i++) cout<<setw(6)<<setprecision(3)<<U[k][i]<<" "; 
				cout<<" => l="<<(int)T_objets[k].lab<<"\n";
			}
		cout<<" iteration "<<++iter<<" : # changements = "<<nch<<"\n"; //char aa; cin>>aa;
	} while (nch>0 && iter<itermax);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) imalab(i,j,1)=T_objets[(*this)(i,j)].lab;
	for (k=0; k<nobj; k++) {
		if (U0[k]!=NULL) delete[] U0[k];
		if (U[k]!=NULL) delete[] U[k];
		if (Tdist[k]!=NULL) delete[] Tdist[k];
		if (Talign[k]!=NULL) delete[] Talign[k];
		if (T_rmse[k]!=NULL) delete[] T_rmse[k];
		if (T_nobjLR[k]!=NULL) delete[] T_nobjLR[k];
	}
	if (U0!=NULL) delete[] U0;
	if (U!=NULL) delete[] U;
	if (Tdist!=NULL) delete[] Tdist;
	if (Talign!=NULL) delete[] Talign;
	if (d_min!=NULL) delete[] d_min;
	if (T_dens!=NULL) delete[] T_dens;
	if (T_rmse!=NULL) delete[] T_rmse;
	if (T_nobjLR!=NULL) delete[] T_nobjLR;
	return imalab;

/*	const float ddmax=45.f;
	const int rES=3; // precision de la direction pour rES=1 : 22.5°, rES=2 : 13.3°, rES=3 : 9.2°, rES=4 : 7.0°, rES=5 : 5.7°, rES=6 : 4.7°, rES=7 : 4.1°
	const int itmax=(int)pow(nblig*nblig+nbcol*nbcol,0.5);
	if (!i_second_param) second_param(obj0);
	int it,no2,nobj2=nobj;
	double d, dd;
	imabin imab_no(nblig,nbcol), imab_nn(nblig,nbcol), imab(nblig,nbcol);
	BYTE *Teq=new BYTE[nobj];
	for (no=0; no<nobj; no++) {
		Teq[no]=no;
		if (no!=obj0) {
			cout<<" objet "<<no<<" : ";
			if (T_objets[no].i_dir_axis) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) {
						imab_no(i,j)=((*this)(i,j)==no);
						imab_nn(i,j)=((*this)(i,j)!=no && (*this)(i,j)!=obj0);
					}
				imab=imab_no&&imab_nn;
				d=T_objets[no].dir_axis;
				dd=T_objets[no].delta_dir;
				cout<<" direction axe principal "<<d;
				if (dd<=ddmax) {
					cout<<", precision "<<dd<<"\n";
					eltstruct ES(2*rES+1,d);
					it=0;
					while (imab.norm()==0 && it<itmax) {
						imab_no=imab_no.dilate(ES); it++;
						imab=imab_no&&imab_nn;
					}
					if (imab.norm()!=0) {
						for (i=0; i<nblig; i++)
							for (j=0; j<nbcol; j++)
								if (imab(i,j)) {
									no2=(*this)(i,j);
									cout<<" objet "<<(int)(*this)(i,j)<<" touche en "<<it<<" dilatations\n";
								}
						Teq[no]=Teq[no2]=mini(Teq[no],Teq[no2]);
						nobj2--;
					} else {
						cout<<" pas d'objet touche en "<<it<<" dilatations\n";
					}
				}
				else {
					cout<<" mal definie\n";
				}
			}
			else {
				cout<<" direction axe principal non definie\n";
			}
		}
	}
	if (nobj2<nobj) {
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				(*this)(i,j)=Teq[(*this)(i,j)];
		nobj=nobj2;
		cout<<" apres fusion des objets alignes : "<<nobj<<" objets dans l'image\n";
		if (T_objets!=NULL) delete[] T_objets;
		T_objets=new objet[nobj];
		basic_param ();
		basic_topologie (obj0);
		second_param (obj0);
	}
	if (Teq!=NULL) delete[] Teq;*/
}

void imaobjetsgeom::add_obj_ima(objetsgeom &O,int yUL, int xUL, float pct) {
	int nl=O.dy, nc=O.dx, i,j,ii,jj;
	for (ii=0; ii<nl; ii++) {
		i=ii+yUL;
		if (i>=0 && i<nblig)
			for (jj=0; jj<nc; jj++) {
				j=jj+xUL;
				if (j>=0 && j<nbcol) {
					if (O.imOG[ii*nc+jj]==true && (pct>=1.f || rand()<=(int)(RAND_MAX*pct)))
						Tima[i*nbcol+j]=1;
				}
			}
	}
}

void imaobjetsgeom::trace_ellipse (double i, double j, double a_u, double b_u, double o_u, BYTE val1, BYTE val2) {
	double x,c_u=pow(pow(a_u,2.)-pow(b_u,2.),0.5), cos_o_u=cos(o_u), sin_o_u=sin(o_u), i_F=i+c_u*sin_o_u, j_F=j+c_u*cos_o_u, i_G=i-c_u*sin_o_u, j_G=j-c_u*cos_o_u; 
//	cout<<" ellipse foyers en ("<<i_F<<","<<j_F<<") et ("<<i_G<<","<<j_G<<"), parametre c="<<c_u<<"\n";
	int ii,jj,nl=nblig,nc=nbcol;
	for (ii=maxi(0,(int)(i-a_u)); ii<=mini((int)(i+a_u),nl-1); ii++)
		for (jj=maxi(0,(int)(j-a_u)); jj<=mini((int)(j+a_u),nc-1); jj++) {
			x=pow(pow(ii-i_F,2.)+pow(jj-j_F,2.),0.5)+pow(pow(ii-i_G,2.)+pow(jj-j_G,2.),0.5);
			if (x<2*a_u+1 /*&& x>2*a_u-0.5*/) (*this)(ii,jj)=val2;
		}
	if (abs(cos_o_u)>1./pow(2.,0.5)) {
		for (jj=maxi(0,(int)(j-a_u*abs(cos_o_u))); jj<=mini((int)(j+a_u*abs(cos_o_u)),nc-1); jj++) {
			ii=(int)(i+sin_o_u/cos_o_u*(jj-j)); (*this)(mini(maxi(ii,0),nl-1),mini(maxi(jj,0),nc-1))=val1;}
		for (ii=maxi(0,(int)(i-b_u*abs(cos_o_u))); ii<=mini((int)(i+b_u*abs(cos_o_u)),nl-1); ii++) {
			jj=(int)(j-sin_o_u/cos_o_u*(ii-i)); (*this)(mini(maxi(ii,0),nl-1),mini(maxi(jj,0),nc-1))=val1;}
	} else {
		for (ii=maxi(0,(int)(i-a_u*abs(sin_o_u))); ii<=mini((int)(i+a_u*abs(sin_o_u)),nl-1); ii++) {
			jj=(int)(j+cos_o_u/sin_o_u*(ii-i)); (*this)(mini(maxi(ii,0),nl-1),mini(maxi(jj,0),nc-1))=val1;}
		for (jj=maxi(0,(int)(j-b_u*abs(sin_o_u))); jj<=mini((int)(j+b_u*abs(sin_o_u)),nc-1); jj++) {
			ii=(int)(i-cos_o_u/sin_o_u*(jj-j)); (*this)(mini(maxi(ii,0),nl-1),mini(maxi(jj,0),nc-1))=val1;}
	}
}
