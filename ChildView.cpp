// ChildView.cpp : implementation of the CChildView class
#include "stdafx.h"
#include "TDI_GUI.h"
#include "ChildView.h"
#include "MainFrm.h"
#include <float.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void FillPalette(RGBQUAD* palette, BYTE type_pal=0)
{
  int i;
  switch(type_pal)
  {
  case 0: // niveaux de gris
    for( i = 0; i < 256; i++ )
    {
      palette[i].rgbBlue = palette[i].rgbGreen = palette[i].rgbRed = (BYTE)i;
      palette[i].rgbReserved = 0;
    }
    break;
  case 1: // niveaux de gris
    for( i = 0; i < 256; i++ )
    {
      palette[255-i].rgbBlue = palette[255-i].rgbGreen = palette[255-i].rgbRed = (BYTE)i;
      palette[255-i].rgbReserved = 0;
    }
    break;
  case 2: // jet bleu->rouge
    for( i = 0; i < 256; i++ )
    {
      int ig=(i<128)?i:(255-i), ib=255-i;
      palette[i].rgbGreen = (ig<32)?0:((ig<3*32)?(ig-32)*4:255);
      palette[i].rgbBlue =(ib<3*32)?0:((ib<5*32)?(ib-3*32)*4:((ib<7*32)?(255):(255-(ib-7*32)*4)));  
      palette[i].rgbRed = (i<3*32)?0:((i<5*32)?(i-3*32)*4:((i<7*32)?(255):(255-(i-7*32)*4)));
      palette[i].rgbReserved = 0;
    }
    break;
  case 3: // jet rouge->bleu
    for( i = 0; i < 256; i++ )
    {
      int ig=(i<128)?i:(255-i), ib=255-i;
      palette[255-i].rgbGreen = (ig<32)?0:((ig<3*32)?(ig-32)*4:255);
      palette[255-i].rgbBlue =(ib<3*32)?0:((ib<5*32)?(ib-3*32)*4:((ib<7*32)?(255):(255-(ib-7*32)*4)));  
      palette[255-i].rgbRed = (i<3*32)?0:((i<5*32)?(i-3*32)*4:((i<7*32)?(255):(255-(i-7*32)*4)));
      palette[255-i].rgbReserved = 0;
    }
    break;
  case 4: // hot noir->jaune->rouge
    for( i = 0; i < 128; i++ )
    {
      palette[i].rgbRed = palette[i].rgbGreen = i*2;
      palette[i].rgbReserved = palette[i].rgbBlue =0;  
    }
    for( i = 0; i < 128; i++ )
    {
      palette[i+128].rgbRed = 255; palette[i+128].rgbGreen = 255-i*2;
      palette[i+128].rgbReserved = palette[i+128].rgbBlue =0;  
    }
    break;
  case 5: // hot inv. rouge->jaune->noir
    for( i = 0; i < 128; i++ )
    {
      palette[255-i].rgbRed = palette[255-i].rgbGreen = i*2;
      palette[255-i].rgbReserved = palette[255-i].rgbBlue =0;  
    }
    for( i = 0; i < 128; i++ )
    {
      palette[127-i].rgbRed = 255; palette[127-i].rgbGreen = 255-i*2;
      palette[127-i].rgbReserved = palette[127-i].rgbBlue =0;  
    }
    break;
  case 6: // hash palette
    palette[0].rgbBlue=0;
    palette[0].rgbRed=0;
    palette[0].rgbGreen=0;
    palette[0].rgbReserved=0;
    for( i = 0; i < 255; i++ )
    {
      int ii=((i&1)<<7)|((i&2)<<5)|((i&4)<<3)|((i&8)<<1)|((i&16)>>1)|((i&32)>>3)|((i&64)>>5)|((i&128)>>7);
      palette[ii+1].rgbBlue=i;
      palette[ii+1].rgbRed=255-i;
      palette[ii+1].rgbGreen=(i<128)?(i<<1):((256-i)<<1);
      palette[ii+1].rgbReserved=0;
    }
    break;
  }
}

// seulement bpp = 8, 24, 32
void  FillBitmapInfo(BITMAPINFO* bmi, DWORD width, DWORD height, BYTE bpp, bool XisDown=true, BYTE type_pal=0)
{
  BITMAPINFOHEADER* bmih = &(bmi->bmiHeader);    
  memset( bmih, 0, sizeof(*bmih));
  bmih->biSize = sizeof(BITMAPINFOHEADER);
  bmih->biWidth = width;
  bmih->biHeight = height;
  if(XisDown) bmih->biHeight= -bmih->biHeight;
  bmih->biPlanes = 1;
  bmih->biBitCount = (unsigned short)bpp;
  bmih->biCompression = BI_RGB;
  if( bpp == 8 ) FillPalette(bmi->bmiColors,type_pal);
}

// CChildView
BEGIN_MESSAGE_MAP(CChildView,CWnd )
	//{{AFX_MSG_MAP(CChildView)
	ON_WM_PAINT()
	ON_WM_CREATE()
	ON_COMMAND(ID_PAL_GRAY, OnPalGray)
	ON_UPDATE_COMMAND_UI(ID_PAL_GRAY, OnUpdatePalGray)
	ON_COMMAND(ID_PAL_BLUE_RED, OnPalBlueRed)
	ON_UPDATE_COMMAND_UI(ID_PAL_BLUE_RED, OnUpdatePalBlueRed)
	ON_COMMAND(ID_PAL_RED_BLUE, OnPalRedBlue)
	ON_UPDATE_COMMAND_UI(ID_PAL_RED_BLUE, OnUpdatePalRedBlue)
	ON_COMMAND(ID_ZOOM_IN, OnZoomIn)
	ON_UPDATE_COMMAND_UI(ID_ZOOM_IN, OnUpdateZoomIn)
	ON_COMMAND(ID_ZOOM_OUT, OnZoomOut)
	ON_UPDATE_COMMAND_UI(ID_ZOOM_OUT, OnUpdateZoomOut)
	ON_WM_MOUSEMOVE()
	ON_COMMAND(ID_PAL_GRAY_INV, OnPalGrayInv)
	ON_UPDATE_COMMAND_UI(ID_PAL_GRAY_INV, OnUpdatePalGrayInv)
	ON_COMMAND(ID_PAL_HOT, OnPalHot)
	ON_COMMAND(ID_PAL_HOT_INV, OnPalHotInv)
	ON_COMMAND(ID_PAL_HASH, OnPalHash)
	ON_UPDATE_COMMAND_UI(ID_PAL_HASH, OnUpdatePalHash)
	ON_UPDATE_COMMAND_UI(ID_PAL_HOT, OnUpdatePalHot)
	ON_UPDATE_COMMAND_UI(ID_PAL_HOT_INV, OnUpdatePalHotInv)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

CChildView::CChildView(): type_pal(0), dimx(256), dimy(256), comp_dc(0), org(0,0), zoom_in(2), is_gray(true)
{
}

bool CChildView::LoadPGMPPM(const char* _fname, BYTE _type_pal)
{
  if(!_fname || *_fname==0) return false;
  FILE* pf=fopen(_fname,"rb");
  if(!pf) return false;
  bool res=LoadPGMPPM(pf,_type_pal);
  fclose(pf);
  return res;
}

#define BUF_SIZE 1024
bool CChildView::LoadPGMPPM(FILE* pf, BYTE _type_pal)
{
  char buf[BUF_SIZE];
  if(!fgets(buf,BUF_SIZE,pf) || (strcmp(buf,"P5\n")!=0 && strcmp(buf,"P6\n")!=0))
    return false; // on ne lit que les sous-formats P5 et P6
  is_gray= (buf[1]=='5');
  do
  {
    fgets(buf,BUF_SIZE,pf);
  } while(buf[0]=='#' && buf[strlen(buf)-1]=='\n');
  if(sscanf(buf,"%lu%lu",&dimx,&dimy)!=2) return false;
  DWORD nblev;
  if((!fgets(buf,BUF_SIZE,pf)) || (sscanf(buf,"%lu",&nblev)!=1)) return false;
  if(is_gray && nblev!=255 && nblev!=65535) return false; // on ne lit que du 8 bits/pixel
  bpp= is_gray ? 8 : 24; is_byte=1;
	if (nblev==65535) {bpp*=2; is_byte=0;}
	printf(" image en niveau de gris ? %d ou couleur; nbre de niveaux = %d -> bpp = %d\n",(int)is_gray,nblev,bpp);
  CleanDC();
  comp_dc=::CreateCompatibleDC(0);
	if (is_byte) {
	  BYTE binfo[sizeof(BITMAPINFO) + 255*sizeof(RGBQUAD)];
		FillBitmapInfo( (BITMAPINFO*)binfo, dimx, dimy, bpp, true, type_pal=_type_pal);
		HBITMAP hbmp=CreateDIBSection(comp_dc,(BITMAPINFO*)binfo,DIB_RGB_COLORS,(void**)&pimg,0,0);
		himg=SelectObject(comp_dc,hbmp);
		dimx_bmp=((dimx*bpp>>3)+3)&0xFFFFFFFC; // la ligne BMP est un multiple de 4 octets 
		for(DWORD y=0; y<dimy; y++) {
	    fread(pimg+dimx_bmp*y,(dimx*bpp>>3),1,pf);
		  if(!is_gray) {
	      BYTE* p=pimg+dimx_bmp*y;
		    for(DWORD x=dimx; x; x--, p+=3) {BYTE v=p[0]; p[0]=p[2]; p[2]=v;}
			}
		}
	} else {
	  BYTE binfo[sizeof(BITMAPINFO) + 255*sizeof(RGBQUAD)];
		FillBitmapInfo( (BITMAPINFO*)binfo, dimx, dimy, bpp, true, type_pal=_type_pal);
		HBITMAP hbmp=CreateDIBSection(comp_dc,(BITMAPINFO*)binfo,DIB_RGB_COLORS,(void**)&pimg,0,0);
		himg=SelectObject(comp_dc,hbmp);
		dimx_bmp=((dimx*bpp>>3)+3)&0xFFFFFFFC; // la ligne BMP est un multiple de 4 octets 
		for(DWORD y=0; y<dimy; y++) {
	    fread(pimg+dimx_bmp*y,(dimx*bpp>>3),1,pf);
//			BYTE* p=pimg+dimx_bmp*y;
//			for(DWORD x=dimx; x; x--, p+=2) {printf("%d %d |",(int)p[0],(int)p[1]);/*BYTE v=p[0]; v=0; p[0]=v;*/}
/*		  if(!is_gray) {
	      BYTE* p=pimg+dimx_bmp*y;
		    for(DWORD x=dimx; x; x--, p+=3) {BYTE v=p[0]; p[0]=p[2]; p[2]=v;}
			}*/
		}
	}
  zoom_in=512/max(dimx,dimy);
  if(zoom_in<1) zoom_in=1;
  if(zoom_in>ZOOM_MAX) zoom_in=ZOOM_MAX;
  return true;
}

void CChildView::ChangePalette(BYTE _type_pal)
{
  RGBQUAD pal[256];
  FillPalette(pal,type_pal=_type_pal);
  SetDIBColorTable(comp_dc,0,256,pal);
  Invalidate();
  ((CMainFrame*)GetParent())->m_wndToolBar.SendMessage(WM_IDLEUPDATECMDUI,1,0);
}

void CChildView::CleanDC()
{
  if(!comp_dc) return;
  if(himg) DeleteObject(SelectObject(comp_dc,himg));
  ::DeleteDC(comp_dc);
}

CChildView::~CChildView()
{
  CleanDC();
}

BOOL CChildView::PreCreateWindow(CREATESTRUCT& cs) 
{
	if (!CWnd::PreCreateWindow(cs)) return FALSE;
	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(NULL, IDC_CROSS), 0/*HBRUSH(COLOR_WINDOW+1)*/, NULL);
	return TRUE;
}

void CChildView::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
  dc.SetMapMode(MM_ISOTROPIC);
  dc.SetWindowOrg(0,0); dc.SetWindowExt(1,1);
  dc.SetViewportOrg(org); dc.SetViewportExt(zoom_in,zoom_in);
  if(comp_dc) ::BitBlt(dc.m_hDC, 0, 0, dimx, dimy, comp_dc, 0, 0, SRCCOPY);
  else
  {
    CRect r; GetClientRect(r);
    dc.FillSolidRect(r,0xA0A0A0);
    dc.TextOut(0,0,"Erreur lors du chargement d'image !");
  }
}

int CChildView::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CWnd ::OnCreate(lpCreateStruct) == -1) return -1;
	return 0;
}

void CChildView::OnPalGray() 
{
  ChangePalette(0);	
}

void CChildView::OnPalGrayInv() 
{
  ChangePalette(1);	
}

void CChildView::OnPalBlueRed() 
{
  ChangePalette(2);	
}

void CChildView::OnPalRedBlue() 
{
  ChangePalette(3);	
}

void CChildView::OnPalHot() 
{
  ChangePalette(4);	
}

void CChildView::OnPalHotInv() 
{
  ChangePalette(5);	
}

void CChildView::OnPalHash()
{
  ChangePalette(6);	
}

void CChildView::OnUpdatePalGray(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(is_gray);
	pCmdUI->SetCheck(type_pal==0);
}

void CChildView::OnUpdatePalGrayInv(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(is_gray);
	pCmdUI->SetCheck(type_pal==1);
}

void CChildView::OnUpdatePalBlueRed(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(is_gray);
	pCmdUI->SetCheck(type_pal==2);
}

void CChildView::OnUpdatePalRedBlue(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(is_gray);
	pCmdUI->SetCheck(type_pal==3);
}

void CChildView::OnUpdatePalHot(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(is_gray);
	pCmdUI->SetCheck(type_pal==4);
}

void CChildView::OnUpdatePalHotInv(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(is_gray);
	pCmdUI->SetCheck(type_pal==5);
}

void CChildView::OnUpdatePalHash(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(is_gray);
	pCmdUI->SetCheck(type_pal==6);
}

void CChildView::OnZoomIn() 
{
  if(zoom_in<ZOOM_MAX) {zoom_in++; ((CMainFrame*)GetParent())->AdjustSize();}
}

void CChildView::OnUpdateZoomIn(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(zoom_in<ZOOM_MAX);
}

void CChildView::OnZoomOut() 
{
  if(zoom_in>1) {zoom_in--; ((CMainFrame*)GetParent())->AdjustSize();}
}

void CChildView::OnUpdateZoomOut(CCmdUI* pCmdUI) 
{
  pCmdUI->Enable(zoom_in>1);
}

void CChildView::OnMouseMove(UINT nFlags, CPoint point) 
{
  if(!comp_dc) return;
	point-=org;
  point.x/=zoom_in; point.y/=zoom_in;
  char buf[128];
  if(is_gray)
    if (is_byte) wsprintf(buf,"%d x %d : %d",point.x,point.y,pimg[dimx_bmp*point.y+point.x]);
		else wsprintf(buf,"%d x %d : %d",point.x,point.y,pimg[dimx_bmp*point.y+2*point.x]*256+pimg[dimx_bmp*point.y+2*point.x+1]);
  else
    wsprintf(buf,"%d x %d : (%d,%d,%d)",point.x,point.y,pimg[dimx_bmp*point.y+3*point.x+2],
      pimg[dimx_bmp*point.y+3*point.x+1],pimg[dimx_bmp*point.y+3*point.x]);    
	((CMainFrame*)GetParent())->SetMessageText(buf);
	CWnd ::OnMouseMove(nFlags, point);
}

