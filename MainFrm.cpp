// MainFrm.cpp : implementation of the CMainFrame class
#include "stdafx.h"
#include "TDI_GUI.h"
#include "MainFrm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// CMainFrame
IMPLEMENT_DYNAMIC(CMainFrame, CFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CFrameWnd)
	//{{AFX_MSG_MAP(CMainFrame)
	ON_WM_CREATE()
	ON_WM_SETFOCUS()
	ON_WM_GETMINMAXINFO()
	ON_WM_ACTIVATE()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

static UINT indicators[] =
{
	ID_SEPARATOR,           // status line indicator
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

CMainFrame::CMainFrame(CString& _fname, BYTE type_pal):fname(_fname), wsize(0,0)
{
	LoadFrame(IDR_MAINFRAME, WS_OVERLAPPEDWINDOW, 0, 0);
//  FloatControlBar(&m_wndToolBar,CPoint(0,0));
  SetMenu(0);
	ShowWindow(SW_SHOW);
  m_wndView.LoadPGMPPM(fname,type_pal);
  AdjustSize();
//  CRect r; GetWindowRect(r);
//  m_wndToolBar.GetToolBarCtrl().SetRows(4,TRUE,r);
//  m_wndToolBar.GetToolBarCtrl().AutoSize();
//  m_wndToolBar.SetWindowPos(0,r.right-20,r.top,0,0,SWP_NOSIZE |SWP_NOZORDER);
//  m_wndToolBar.MoveWindow(r.right-20,r.top,r.right+20,r.top+60);
}

CMainFrame::~CMainFrame()
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CFrameWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	// create a view to occupy the client area of the frame
	if (!m_wndView.Create(NULL, NULL, AFX_WS_DEFAULT_VIEW,
		CRect(0, 0, 0, 0), this, AFX_IDW_PANE_FIRST, NULL))
	{
		TRACE0("Failed to create view window\n");
		return -1;
	}
	
	if (!m_wndToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP | CBRS_GRIPPER | CBRS_TOOLTIPS
    | CBRS_FLYBY | CBRS_SIZE_DYNAMIC) || !m_wndToolBar.LoadToolBar(IDR_MAINFRAME))
	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}

	if (!m_wndStatusBar.Create(this) || !m_wndStatusBar.SetIndicators(indicators,sizeof(indicators)/sizeof(UINT)))
	{
		TRACE0("Failed to create status bar\n");
		return -1;      // fail to create
	}

	m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
	EnableDocking(CBRS_ALIGN_ANY);
	DockControlBar(&m_wndToolBar);
	return 0;
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	if( !CFrameWnd::PreCreateWindow(cs) ) return FALSE;
	cs.dwExStyle &= ~WS_EX_CLIENTEDGE;
  cs.dwExStyle |= WS_EX_TOOLWINDOW;
	cs.lpszClass = AfxRegisterWndClass(0);
  m_strTitle=fname;
	return TRUE;
}

void CMainFrame::OnSetFocus(CWnd* pOldWnd)
{
  m_wndView.SetFocus();
}

BOOL CMainFrame::OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
{
	if (m_wndView.OnCmdMsg(nID, nCode, pExtra, pHandlerInfo)) return TRUE;
	return CFrameWnd::OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
}

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CFrameWnd::Dump(dc);
}
#endif //_DEBUG


void CMainFrame::OnGetMinMaxInfo(MINMAXINFO FAR* lpMMI) 
{
  if(wsize.x)
  {
    lpMMI->ptMinTrackSize=lpMMI->ptMaxTrackSize=wsize;
  }
	CFrameWnd::OnGetMinMaxInfo(lpMMI);
}

void CMainFrame::OnActivate(UINT nState, CWnd* pWndOther, BOOL bMinimized) 
{
	CFrameWnd::OnActivate(nState, pWndOther, bMinimized);
  if(m_wndToolBar.m_hWnd==0) return;
#if 0
  if(nState==WA_INACTIVE) ShowControlBar(&m_wndToolBar,FALSE,FALSE);
  else
  {
    CRect r; GetWindowRect(r);
    m_wndToolBar.SetWindowPos(0,r.right-20,r.top,0,0,SWP_NOSIZE |SWP_NOZORDER);
    ShowControlBar(&m_wndToolBar,TRUE,FALSE);
  }
#endif
}
