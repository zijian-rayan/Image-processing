#include "imasites.h"

void imasites::affiche() const {
	for (long int i=0; i<nbpix; i++) cout<<Tsites[i]<<" ";
	cout<<"\n";
}

void imasites::sousIma(int nblig2, int nbcol2, int lig0, int col0, int l_ech, int c_ech) {
	nblig2=mini(nblig2,nblig); nbcol2=mini(nbcol2,nbcol);
	lig0=maxi(lig0,0); col0=maxi(col0,0);
	l_ech=maxi(l_ech,1); c_ech=maxi(c_ech,1);
	int i, j, i2=-1;
	long int nbpix2=((nblig2-lig0)/l_ech+1)*((nbcol2-col0)/c_ech+1), n2=-1;
	void **Tsites2=new void*[nbpix2];
	for (i=lig0; i<nblig2; i+=l_ech) {
		i2++;
		for (j=col0; j<nbcol2; j+=c_ech) Tsites2[++n2]=Tsites[(long int)i*nbcol+j];
	}
	if (Tsites!=NULL) delete[] Tsites;
	Tsites=Tsites2;
	Tsites2=NULL;
  nblig=i2+1; nbpix=n2+1; nbcol=nbpix/nblig;
}

void imasites::zoomIma(int l_zoom, int c_zoom) {
	l_zoom=maxi(l_zoom,1);
	c_zoom=maxi(c_zoom,1);
	int nblig2=l_zoom*nblig, nbcol2=c_zoom*nbcol, i, j, i2, j2;
	long int nbpix2=nblig2*nbcol2, n2=-1;
	void **Tsites2=new void*[nbpix2];
	for (i=0; i<nblig; i++)
		for (i2=0; i2<l_zoom; i2++)
			for (j=0; j<nbcol; j++)
				for (j2=0; j2<c_zoom; j2++)
					Tsites2[++n2]=Tsites[(long int)i*nbcol+j];
	if (Tsites!=NULL) delete[] Tsites;
	Tsites=Tsites2;
	Tsites2=NULL;
	nblig=nblig2;
	nbcol=nbcol2;
	nbpix=nbpix2;
}

