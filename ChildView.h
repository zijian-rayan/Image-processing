// ChildView.h : interface of the CChildView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_CHILDVIEW_H__B22483AE_1E94_4B2E_A17F_739D218433A5__INCLUDED_)
#define AFX_CHILDVIEW_H__B22483AE_1E94_4B2E_A17F_739D218433A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// CChildView window
#define ZOOM_MAX 6
class CChildView : public CWnd
{
  HDC comp_dc;
  HGDIOBJ himg;
  BYTE* pimg, bpp;
  DWORD dimx, dimy, dimx_bmp;
  BYTE type_pal;
  CPoint org;
  int zoom_in;
  bool is_gray; bool is_byte;
public:
	CChildView();
	virtual ~CChildView();

  bool LoadPGMPPM(const char* fname, BYTE _type_pal=0);
  bool LoadPGMPPM(FILE* pf, BYTE _type_pal=0);
  void GetViewSize(CPoint& pv) {pv.x=dimx*zoom_in; pv.y=dimy*zoom_in;}
  void ChangePalette(BYTE _pal_type);
  void CleanDC();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CChildView)
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	//}}AFX_VIRTUAL

// Implementation
public:

	// Generated message map functions
protected:
	//{{AFX_MSG(CChildView)
	afx_msg void OnPaint();
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnPalGray();
	afx_msg void OnUpdatePalGray(CCmdUI* pCmdUI);
	afx_msg void OnPalBlueRed();
	afx_msg void OnUpdatePalBlueRed(CCmdUI* pCmdUI);
	afx_msg void OnPalRedBlue();
	afx_msg void OnUpdatePalRedBlue(CCmdUI* pCmdUI);
	afx_msg void OnZoomIn();
	afx_msg void OnUpdateZoomIn(CCmdUI* pCmdUI);
	afx_msg void OnZoomOut();
	afx_msg void OnUpdateZoomOut(CCmdUI* pCmdUI);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnPalGrayInv();
	afx_msg void OnUpdatePalGrayInv(CCmdUI* pCmdUI);
	afx_msg void OnPalHot();
	afx_msg void OnPalHotInv();
	afx_msg void OnPalHash();
	afx_msg void OnUpdatePalHash(CCmdUI* pCmdUI);
	afx_msg void OnUpdatePalHot(CCmdUI* pCmdUI);
	afx_msg void OnUpdatePalHotInv(CCmdUI* pCmdUI);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_CHILDVIEW_H__B22483AE_1E94_4B2E_A17F_739D218433A5__INCLUDED_)
