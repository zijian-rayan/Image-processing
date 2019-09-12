#include "tree.h"

void tree::Arbre_Min(imadata<int> ima, bool iaf) {
	if (iaf) {
		cout<<"\n *********************************************************************** \n";
		cout<<"\n ********************* Construction de l'Arbre Min ********************* \n";
		cout<<"\n *********************************************************************** \n";
		ima.statbasic(1);
	}
	ima.supprime_Niveaux_vide (); if (iaf) {ima.statbasic(1); ima.affiche();}
	const int nl=ima.nlig(),nc=ima.ncol(),val_min=(int)ima.minI(),val_max=(int)ima.maxI();
	int i,i0,i2,ii,j,j0,j2,jj,imloc,jmloc,k,l,val,val_N,i_val_N,j_val_N,n=0;
	bool minloc;
	for (i=0; i<nl; i++) {
		i0=maxi(i-1,0); i2=mini(i+1,nl-1);
		for (j=0; j<nc; j++) {
			j0=maxi(j-1,0); j2=mini(j+1,nc-1);
			minloc=1; val=ima(i,j);
			for (ii=i0; ii<=i2; ii++)
				for (jj=j0; jj<=j2; jj++)
					if (minloc && ima(ii,jj)>val) minloc=0;
			n+=(int)minloc;
		}
	}
	if (iaf) cout<<" nombre de minima locaux sur l'image = "<<n<<"\n";
	tree *T_tree=new tree[n];
	int **lien_T_tree=new int*[n];
	for (k=0; k<n; k++) {
		lien_T_tree[k]=new int[3];
		lien_T_tree[k][0]=0; lien_T_tree[k][1]=lien_T_tree[k][2]=-1;
	}
	n=0;
	imadata<int> imaCC(nl,nc), ima_noTree(nl,nc);
	imabin imaCCb(nl,nc);
	int nfils=0,nfrere=0,type_lien=0;
	unsigned int nb_ccN,nb_trous,nb_tree=1,i_tree=0;
	bool init, fini, fincc, insere_arbreP=0;
	liste_pixels L_CC, L_A, L_B, L_N, L_vide, L_ccN;
	elt_tree E;
	elt_liste p;

	for (i=0; i<nl; i++) {
		i0=maxi(i-1,0); i2=mini(i+1,nl-1);
		for (j=0; j<nc; j++) {
			val=ima(i,j); 
			if (imaCC(i,j)<=0 && (val==val_min || n>0)) {
				j0=maxi(j-1,0); j2=mini(j+1,nc-1);
				minloc=1; ii=i0; jj=j0; imloc=i; jmloc=j;
				while (minloc && ii<=i2 && jj<=j2) {
					if (ima(ii,jj)<val) {minloc=0; imloc=-1; jmloc=-1;}
					else {
						jj++;
						if (jj>j2) {jj=j0; ii++;}
					}
				}
				if (minloc) {
					if (iaf) cout<<" minimum local trouve en ("<<imloc<<","<<jmloc<<")\n";
					n++; L_A.insere(imloc,jmloc); fini=0; init=1;
				} else fini=1;
			}
			while (!fini) {
				if (init) {
					if (iaf) {
						cout<<" a l'initialisation : val = "<<val<<" n = "<<n<<"\n liste A = \n"; L_A.affiche();
						cout<<"\n liste N = \n"; L_N.affiche(); cout<<"\n liste CC = \n"; L_CC.affiche(); cout<<"\n";
					}
					E.val=(float)val; E.ident=n; E.fils=E.frere=NULL; 
					imaCCb.mise_a_zero(); 
					L_B=L_CC;  
					while (L_B.nb_elts()>0) {
						p=L_B.extrait(); imaCCb(p.x,p.y)=1;
					}
					val_N=val_max; i_val_N=j_val_N=0;
					L_B=L_N;  
					while (L_B.nb_elts()>0) {
						p=L_B.extrait(); l=p.x; k=p.y;
						if (ima(l,k)<val_N) {val_N=ima(l,k); i_val_N=l; j_val_N=k;}
					}
					init=0; fincc=0;
				}
				L_B=L_A;
				while (L_B.nb_elts()>0) {
					p=L_B.extrait(); l=p.x; k=p.y;
					imaCCb(l,k)=1; 
					i0=maxi(l-1,0); i2=mini(l+1,nl-1); j0=maxi(k-1,0); j2=mini(k+1,nc-1);
					for (ii=i0; ii<=i2; ii++)
						for (jj=j0; jj<=j2; jj++)
							if (imaCCb(ii,jj)==0 && imaCC(ii,jj)<=0 && !L_N.existe(ii,jj) && !L_B.existe(ii,jj)) {
								L_N.insere(ii,jj);
								if (ima(ii,jj)<val_N) {val_N=ima(ii,jj); i_val_N=ii; j_val_N=jj;}
								if (ima_noTree(ii,jj)>0 && ima_noTree(ii,jj)!=i_tree+1) {
									if (iaf) cout<<" rencontre d'un arbre partiel deja construit en ("<<ii<<","<<jj<<")\n";
									i_val_N=ii; j_val_N=jj;
									insere_arbreP=1;
								}
							}
				}
				if (insere_arbreP) {
					type_lien=2; nfrere=lien_T_tree[i_tree][2];
					if (iaf) cout<<" insertion dans l'arbre "<<lien_T_tree[i_tree][1]<<" d'une composante frere du numero "<<nfrere<<"\n";
					T_tree[lien_T_tree[i_tree][1]].insere(T_tree[i_tree].racine,nfrere,type_lien);
					for (ii=0; ii<nl; ii++)
						for (jj=0; jj<nc; jj++) if (ima_noTree(ii,jj)==i_tree+1) ima_noTree(ii,jj)=lien_T_tree[i_tree][1];
					lien_T_tree[i_tree][0]=0; //lien_T_tree[i_tree][1]=-1; lien_T_tree[i_tree][2]=-1;
					insere_arbreP=0;
					i_tree=lien_T_tree[i_tree][1];
					type_lien=1; nfils=nfrere; 
					if (iaf) {cout<<"\n etat courant de l'arbre\n"; T_tree[i_tree].affiche(); }
				}
				if (val_N>val) {
					if (iaf) {cout<<" la composante parent est atteinte\n"; imaCCb.affiche();}
					fincc=1;
				} else {
					if (val_N==val) {
						if (iaf) cout<<" la composante n'est pas finie\n";
						while (L_A.nb_elts()>0) {
							p=L_A.extrait(); l=p.x; k=p.y;
							L_CC.insere(l,k); imaCC(l,k)=n; ima_noTree(l,k)=i_tree+1;
						}
						L_B=L_N; L_N=L_vide; val_N=val_max; 
						while (L_B.nb_elts()>0) {
							p=L_B.extrait(); l=p.x; k=p.y;
							if (ima(l,k)==val) L_A.insere(l,k);
							else {L_N.insere(l,k); if (ima(l,k)<val_N) val_N=ima(l,k); }
						}
						if (L_A.nb_elts()==0 && L_N.nb_elts()==0) {
							fincc=1; fini=1; if (iaf) cout<<" en fait si\n";
						}
					} else {
						if (iaf) cout<<" la composante est interrompue pour cause de rencontre d'un autre fils en ("<<i_val_N<<","<<j_val_N<<")\n";
						if (T_tree[i_tree].nb_elts==0) {
							if (iaf) cout<<" en fait l'arbre courant etait encore vide\n";

						} else {
							lien_T_tree[nb_tree][0]=1; lien_T_tree[nb_tree][1]=i_tree; lien_T_tree[nb_tree][2]=nfils;
							if (iaf) {cout<<" demarrage d'un nouvel arbre partiel qui sera relie a l'arbre "<<lien_T_tree[nb_tree][1]<<" par le noeud "<<lien_T_tree[nb_tree][2]<<"\n";}
							i_tree=nb_tree; nb_tree++;
							nfils=nfrere=type_lien=0;
							insere_arbreP=0;
						}
						while (L_CC.nb_elts()>0) {
							p=L_CC.extrait(); l=p.x; k=p.y;
							ima(l,k)=val; imaCC(l,k)=0;
						}
//						if (iaf) ima.affiche();
//						while (L_N.nb_elts()>0) {
//							p=L_N.extrait(); l=p.x; k=p.y;
//							if (ima(l,k)<val) {imloc=l; jmloc=k; val=ima(l,k);}
//						}
						while (L_N.nb_elts()>0) p=L_N.extrait(); //l=p.x; k=p.y;
						imloc=i_val_N; jmloc=j_val_N; val=ima(imloc,jmloc);
						if (iaf) cout<<" minimum local trouve en ("<<imloc<<","<<jmloc<<")\n";
						L_A=L_vide; L_A.insere(imloc,jmloc); 
						val=ima(imloc,jmloc); init=1; 
					}
				}
				if (fincc) {
					if (iaf) {cout<<" image de la composante connexe :\n"; imaCCb.affiche();}
					while (L_A.nb_elts()>0) {
						p=L_A.extrait(); l=p.x; k=p.y; //cout<<" on ajoute pixel de coord. "<<l<<" "<<k<<"\n";
						L_CC.insere(l,k); imaCC(l,k)=n; ima_noTree(l,k)=i_tree+1;
					}
					if (iaf) {imaCC.affiche(); ima_noTree.affiche();}
					E.L=L_CC;
					L_B=L_N; //if (iaf) {cout<<" liste L_B : "; L_B.affiche(); cout<<"\n";}
					for (i=0; i<nl; i++) {
						if (imaCCb(i,0)==1) L_B.insere(i,-1);
						if (imaCCb(i,nc-1)==1) L_B.insere(i,nc);
					}
					for (j=0; j<nc; j++) {
						if (imaCCb(0,j)==1) L_B.insere(-1,j);
						if (imaCCb(nl-1,j)==1) L_B.insere(nl,j);
					}
//					if (iaf) {cout<<" liste L_B : "; L_B.affiche(); cout<<"\n";}
					nb_ccN=L_B.nb_cc(&L_ccN);
					if (nb_ccN>1) {
						nb_trous=nb_ccN-1;
						p=L_ccN.extrait();
						E.L_ini_hole=L_ccN;
					} else E.L_ini_hole=L_vide;
					if (type_lien==0) {
						if (iaf) cout<<" insertion dans l'arbre de la premiere composante\n";
						T_tree[i_tree].insere(&E,0,0);
						nfils=n; 
						if (iaf) {cout<<"\n etat courant de l'arbre\n"; T_tree[i_tree].affiche(); char aa; cin>>aa;}
					} else {
						if (type_lien==1) {
							if (iaf) cout<<" insertion dans l'arbre d'une composante pere\n";
							T_tree[i_tree].insere_pere(&E,nfils);
							nfils=n; nfrere=0; 
							if (iaf) {cout<<"\n etat courant de l'arbre\n"; T_tree[i_tree].affiche(); char aa; cin>>aa;}
						} else {
//							if (type_lien==2) {
//								if (iaf) cout<<" insertion dans l'arbre d'une composante frere\n";
//								T_tree[i_tree].insere(&E,nfrere,type_lien);
//								nfrereaine=n; if (iaf) cout<<" num. frere aine = "<<nfrereaine<<"\n";
//							}
						}
					}
					if (!fini) {
						type_lien=1;
						L_B=L_N; L_N=L_vide;  
						while (L_B.nb_elts()>0) {
							p=L_B.extrait(); l=p.x; k=p.y;
							if (ima(l,k)==val_N) L_A.insere(l,k);
							else {L_N.insere(l,k); if (ima(l,k)<val_N) val_N=ima(l,k); }
						}
						n++; val=val_N; init=1;
					}
				}
//				char aa; cin>>aa;
			}
		}
	}
	if (iaf) {
		cout<<"\n *************************\n Arbre Min obtenu : \n";
		T_tree[i_tree].affiche(); //cout<<"nb etages = "<<Min_Tree.calc_etages()<<"\n"; 
		cout<<"\n *************************\n";
	}
	(*this)=T_tree[i_tree];
	for (k=0; k<n; k++) if (lien_T_tree[k]!=NULL) delete[] lien_T_tree[k];
	if (lien_T_tree!=NULL) delete[] lien_T_tree;
	if (T_tree!=NULL) delete[] T_tree;
}

void tree::Arbre_Max(imadata<int> ima, bool iaf) {
	if (iaf) {
		cout<<"\n *********************************************************************** \n";
		cout<<"\n ********************* Construction de l'Arbre Max ********************* \n";
		cout<<"\n *********************************************************************** \n";
		ima.statbasic(1);
	}
	ima.supprime_Niveaux_vide (); if (iaf) {ima.statbasic(1); ima.affiche();}
	const int nl=ima.nlig(),nc=ima.ncol(),val_min=(int)ima.minI(),val_max=(int)ima.maxI();
	int i,i0,i2,ii,j,j0,j2,jj,imloc,jmloc,k,l,val,val_N,i_val_N,j_val_N,n=0;
	bool maxloc;
	for (i=0; i<nl; i++) {
		i0=maxi(i-1,0); i2=mini(i+1,nl-1);
		for (j=0; j<nc; j++) {
			j0=maxi(j-1,0); j2=mini(j+1,nc-1);
			maxloc=1; val=ima(i,j);
			for (ii=i0; ii<=i2; ii++)
				for (jj=j0; jj<=j2; jj++)
					if (maxloc && ima(ii,jj)>val) maxloc=0;
			n+=(int)maxloc;
		}
	}
	if (iaf) cout<<" nombre de maxima locaux sur l'image = "<<n<<"\n";
	tree *T_tree=new tree[n];
	int **lien_T_tree=new int*[n];
	for (k=0; k<n; k++) {
		lien_T_tree[k]=new int[3];
		lien_T_tree[k][0]=0; lien_T_tree[k][1]=lien_T_tree[k][2]=-1;
	}
	n=0;
	imadata<int> imaCC(nl,nc), ima_noTree(nl,nc);
	imabin imaCCb(nl,nc);
	int nfils=0,nfrere=0,type_lien=0;
	unsigned int nb_ccN,nb_trous,nb_tree=1,i_tree=0;
	bool init, fini, fincc, insere_arbreP=0;
	liste_pixels L_CC, L_A, L_B, L_N, L_vide, L_ccN;
	elt_tree E;
	elt_liste p;

	for (i=0; i<nl; i++) {
		i0=maxi(i-1,0); i2=mini(i+1,nl-1);
		for (j=0; j<nc; j++) {
			val=ima(i,j); 
			if (imaCC(i,j)<=0 && (val==val_max || n>0)) {
				j0=maxi(j-1,0); j2=mini(j+1,nc-1);
				maxloc=1; ii=i0; jj=j0; imloc=i; jmloc=j;
				while (maxloc && ii<=i2 && jj<=j2) {
					if (ima(ii,jj)>val) {maxloc=0; imloc=-1; jmloc=-1;}
					else {
						jj++;
						if (jj>j2) {jj=j0; ii++;}
					}
				}
				if (maxloc) {
					if (iaf) cout<<" maximum local trouve en ("<<imloc<<","<<jmloc<<")\n";
					n++; L_A.insere(imloc,jmloc); fini=0; init=1;
				} else fini=1;
			}
			while (!fini) {
				if (init) {
					if (iaf) {
						cout<<" a l'initialisation : val = "<<val<<" n = "<<n<<"\n liste A = \n"; L_A.affiche();
						cout<<"\n liste N = \n"; L_N.affiche(); cout<<"\n liste CC = \n"; L_CC.affiche(); cout<<"\n";
					}
					E.val=(float)val; E.ident=n; E.fils=E.frere=NULL; 
					imaCCb.mise_a_zero(); 
					L_B=L_CC;  
					while (L_B.nb_elts()>0) {
						p=L_B.extrait(); imaCCb(p.x,p.y)=1;
					}
					val_N=val_min; i_val_N=j_val_N=0;
					L_B=L_N;  
					while (L_B.nb_elts()>0) {
						p=L_B.extrait(); l=p.x; k=p.y;
						if (ima(l,k)>val_N) {val_N=ima(l,k); i_val_N=l; j_val_N=k;}
					}
					init=0; fincc=0;
				}
				L_B=L_A;
				while (L_B.nb_elts()>0) {
					p=L_B.extrait(); l=p.x; k=p.y;
					imaCCb(l,k)=1; 
					i0=maxi(l-1,0); i2=mini(l+1,nl-1); j0=maxi(k-1,0); j2=mini(k+1,nc-1);
					for (ii=i0; ii<=i2; ii++)
						for (jj=j0; jj<=j2; jj++)
							if (imaCCb(ii,jj)==0 && imaCC(ii,jj)<=0 && !L_N.existe(ii,jj) && !L_B.existe(ii,jj)) {
								L_N.insere(ii,jj);
								if (ima(ii,jj)>val_N) {val_N=ima(ii,jj); i_val_N=ii; j_val_N=jj;}
								if (ima_noTree(ii,jj)>0 && ima_noTree(ii,jj)!=i_tree+1) {
									if (iaf) cout<<" rencontre d'un arbre partiel deja construit en ("<<ii<<","<<jj<<")\n";
									i_val_N=ii; j_val_N=jj;
									insere_arbreP=1;
								}
							}
				}
				if (insere_arbreP) {
					type_lien=2; nfrere=lien_T_tree[i_tree][2];
					if (iaf) cout<<" insertion dans l'arbre "<<lien_T_tree[i_tree][1]<<" d'une composante frere du numero "<<nfrere<<"\n";
					T_tree[lien_T_tree[i_tree][1]].insere(T_tree[i_tree].racine,nfrere,type_lien);
					for (ii=0; ii<nl; ii++)
						for (jj=0; jj<nc; jj++) if (ima_noTree(ii,jj)==i_tree+1) ima_noTree(ii,jj)=lien_T_tree[i_tree][1];
					lien_T_tree[i_tree][0]=0;
					i_tree=lien_T_tree[i_tree][1];
					insere_arbreP=0; type_lien=1; nfils=nfrere; 
					if (iaf) {cout<<"\n etat courant de l'arbre\n"; T_tree[i_tree].affiche(); }
				}
				if (val_N<val) {
					if (iaf) {cout<<" la composante parent est atteinte\n"; imaCCb.affiche();}
					fincc=1;
				} else {
					if (val_N==val) {
						if (iaf) cout<<" la composante n'est pas finie\n";
						while (L_A.nb_elts()>0) {
							p=L_A.extrait(); l=p.x; k=p.y;
							L_CC.insere(l,k); imaCC(l,k)=n; ima_noTree(l,k)=i_tree+1;
						}
						L_B=L_N; L_N=L_vide; val_N=val_min; 
						while (L_B.nb_elts()>0) {
							p=L_B.extrait(); l=p.x; k=p.y;
							if (ima(l,k)==val) L_A.insere(l,k);
							else {L_N.insere(l,k); if (ima(l,k)>val_N) val_N=ima(l,k);}
						}
						if (L_A.nb_elts()==0 && L_N.nb_elts()==0) {
							fincc=1; fini=1; if (iaf) cout<<" en fait si\n";
						}
					} else {
						if (iaf) cout<<" la composante est interrompue pour cause de rencontre d'un autre fils en ("<<i_val_N<<","<<j_val_N<<")\n";
						if (T_tree[i_tree].nb_elts==0) {
							if (iaf) cout<<" en fait l'arbre courant etait encore vide\n";
						} else {
							lien_T_tree[nb_tree][0]=1; lien_T_tree[nb_tree][1]=i_tree; lien_T_tree[nb_tree][2]=nfils;
							if (iaf) {cout<<" demarrage d'un nouvel arbre partiel qui sera relie a l'arbre "<<lien_T_tree[nb_tree][1]<<" par le noeud "<<lien_T_tree[nb_tree][2]<<"\n"; char aa; cin>>aa;}
							i_tree=nb_tree; nb_tree++;
							nfils=nfrere=type_lien=0;
							insere_arbreP=0;
						}
						while (L_CC.nb_elts()>0) {
							p=L_CC.extrait(); l=p.x; k=p.y;
							ima(l,k)=val; imaCC(l,k)=0;
						}
						while (L_N.nb_elts()>0) p=L_N.extrait();
						imloc=i_val_N; jmloc=j_val_N; val=ima(imloc,jmloc);
						if (iaf) cout<<" maximum local trouve en ("<<imloc<<","<<jmloc<<")\n";
						L_A=L_vide; L_A.insere(imloc,jmloc); 
						val=ima(imloc,jmloc); init=1; 
					}
				}
				if (fincc) {
					if (iaf) {cout<<" image de la composante connexe :\n"; imaCCb.affiche();}
					while (L_A.nb_elts()>0) {
						p=L_A.extrait(); l=p.x; k=p.y;
						L_CC.insere(l,k); imaCC(l,k)=n; ima_noTree(l,k)=i_tree+1;
					}
					if (iaf) {imaCC.affiche(); ima_noTree.affiche();}
					E.L=L_CC;
					L_B=L_N;
					for (i=0; i<nl; i++) {
						if (imaCCb(i,0)==1) L_B.insere(i,-1);
						if (imaCCb(i,nc-1)==1) L_B.insere(i,nc);
					}
					for (j=0; j<nc; j++) {
						if (imaCCb(0,j)==1) L_B.insere(-1,j);
						if (imaCCb(nl-1,j)==1) L_B.insere(nl,j);
					}
					nb_ccN=L_B.nb_cc(&L_ccN);
					if (nb_ccN>1) {
						nb_trous=nb_ccN-1;
						p=L_ccN.extrait();
						E.L_ini_hole=L_ccN;
					} else E.L_ini_hole=L_vide;
					if (type_lien==0) {
						if (iaf) cout<<" insertion dans l'arbre de la premiere composante\n";
						T_tree[i_tree].insere(&E,0,0);
						nfils=n; 
						if (iaf) {cout<<"\n etat courant de l'arbre\n"; T_tree[i_tree].affiche(); char aa; cin>>aa;}
					} else {
						if (type_lien==1) {
							if (iaf) cout<<" insertion dans l'arbre d'une composante pere\n";
							T_tree[i_tree].insere_pere(&E,nfils);
							nfils=n; nfrere=0; 
							if (iaf) {cout<<"\n etat courant de l'arbre\n"; T_tree[i_tree].affiche(); char aa; cin>>aa;}
						}
					}
					if (!fini) {
						type_lien=1;
						L_B=L_N; L_N=L_vide;  
						while (L_B.nb_elts()>0) {
							p=L_B.extrait(); l=p.x; k=p.y;
							if (ima(l,k)==val_N) L_A.insere(l,k);
							else {L_N.insere(l,k); if (ima(l,k)>val_N) val_N=ima(l,k); }
						}
						n++; val=val_N; init=1;
					}
				}
//				char aa; cin>>aa;
			}
		}
	}
	if (iaf) {
		cout<<"\n *************************\n Arbre Max obtenu : \n";
		T_tree[i_tree].affiche();
		cout<<"\n *************************\n";
	}
	(*this)=T_tree[i_tree];
	for (k=0; k<n; k++) if (lien_T_tree[k]!=NULL) delete[] lien_T_tree[k];
	if (lien_T_tree!=NULL) delete[] lien_T_tree;
	if (T_tree!=NULL) delete[] T_tree;
}

int tree::calc_etages (bool iaf) {
	int n;
	nb_etages=0;
	liste_ad_elt_tree liste_adres(nb_elts);
	liste_adres.L[(liste_adres.n_adres)++].adres_elt=racine;
	courant=racine;
	courant->etage=0; nb_etages++;
	if (iaf) {cout<<" ############## element \n"; affiche(courant); cout<<" a l'etage "<<(int)(courant->etage)<<"\n";}
	while (liste_adres.n_adres>0 && courant!=NULL) {
		(liste_adres.n_adres)--;
		if (courant->frere!=NULL) { if (iaf) cout<<" trouve frere\n";
			n=(liste_adres.n_adres)++;
			liste_adres.L[n].adres_elt=courant->frere;
			liste_adres.L[n].ad_origine=courant;
			liste_adres.L[n].type_origin=2;
		}
		if (courant->fils!=NULL) { if (iaf) cout<<" trouve fils\n";
			n=(liste_adres.n_adres)++;
			liste_adres.L[n].adres_elt=courant->fils;
			liste_adres.L[n].ad_origine=courant;
			liste_adres.L[n].type_origin=1;
		}
		n=liste_adres.n_adres-1;
		courant=liste_adres.L[n].adres_elt;
		if (courant!=NULL && n>=0) {
			courant->etage=(liste_adres.L[n].ad_origine)->etage+(2-liste_adres.L[n].type_origin);
			if (iaf) {cout<<" ############## element \n"; affiche(courant); cout<<" a l'etage "<<(int)(courant->etage)<<"\n";}
			n=(int)(courant->etage)+1;
			if (n>nb_etages) nb_etages=n;
		}
	}
	courant=racine;
	return nb_etages;
}

imadata<int> tree::image_Tree (int nl, int nc, bool iaf) {
	int n_et=calc_etages(iaf); cout<<" nombre d'etages = "<<n_et<<"\n";
	imadata<int> ima(nl,nc,nb_etages);
	int i, n, val, n_adres;
	elt_tree **liste_adres=new elt_tree*[nb_elts];
	for (i=0; i<nb_elts; i++) liste_adres[i]=NULL;
	elt_liste *courantL;
	n_adres=0;
	liste_adres[n_adres++]=racine;
	courant=liste_adres[n_adres-1];
	while (n_adres>0 && courant!=NULL) {
		n_adres--;
		if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere;
		if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils;
		courantL=(courant->L).debut;
		n=courant->etage;
		val=courant->ident; //val=courant->val;
		while (courantL!=NULL) {
			ima(courantL->x,courantL->y,n)=val;
			courantL=courantL->suivant;
		}
		courant=liste_adres[n_adres-1];
	}
	if (nb_elts>0 && liste_adres!=NULL) delete[] liste_adres;
	courant=racine;
	ima.sauve_ImaBSQ("Ima_Arbre.dat");
	return ima;
}
