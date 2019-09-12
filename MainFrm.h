// MainFrm.h : interface of the CMainFrame class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_MAINFRM_H__51A2BF3D_54FF_40EE_888A_12861C981183__INCLUDED_)
#define AFX_MAINFRM_H__51A2BF3D_54FF_40EE_888A_12861C981183__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "ChildView.h"
#define WM_IDLEUPDATECMDUI  0x0363

class CMainFrame : public CFrameWnd
{
public:
	CString fname;
  CPoint wsize;
	CMainFrame(CString& _fname, BYTE type_pal=0);
protected: 
	DECLARE_DYNAMIC(CMainFrame)

// Attributes
public:

// Operations
public:
  void AdjustSize()
  {
    CRect r1, r2;
    CPoint pv;
    m_wndView.GetViewSize(pv);
    m_wndView.GetClientRect(r1);
    GetWindowRect(r2);
    wsize.x=r2.Width()+pv.x-r1.Width(); wsize.y=r2.Height()+pv.y-r1.Height();
    SetWindowPos(0,0,0,wsize.x,wsize.y,SWP_NOMOVE|SWP_NOZORDER);
    m_wndToolBar.SendMessage(WM_IDLEUPDATECMDUI,1,0);
  }

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMainFrame)
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual BOOL OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CMainFrame();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	CStatusBar  m_wndStatusBar;
	CToolBar    m_wndToolBar;
	CChildView  m_wndView;

// Generated message map functions
protected:
	//{{AFX_MSG(CMainFrame)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnSetFocus(CWnd *pOldWnd);
	afx_msg void OnGetMinMaxInfo(MINMAXINFO FAR* lpMMI);
	afx_msg void OnActivate(UINT nState, CWnd* pWndOther, BOOL bMinimized);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MAINFRM_H__51A2BF3D_54FF_40EE_888A_12861C981183__INCLUDED_)
