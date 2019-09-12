// TDI_GUI.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "TDI_GUI.h"
#include "TDI_GUIDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CTDI_GUIApp theApp;

// CTDI_GUIApp
BEGIN_MESSAGE_MAP(CTDI_GUIApp, CWinApp)
	//{{AFX_MSG_MAP(CTDI_GUIApp)
	//}}AFX_MSG
	ON_COMMAND(ID_HELP, CWinApp::OnHelp)
END_MESSAGE_MAP()

CTDI_GUIApp::CTDI_GUIApp()
{
}

BOOL CTDI_GUIApp::InitInstance()
{
#ifdef _AFXDLL
	Enable3dControls();			// Call this when using MFC in a shared DLL
#else
	Enable3dControlsStatic();	// Call this when linking to MFC statically
#endif
	CTDI_GUIDlg dlg;
	m_pMainWnd = &dlg;
	dlg.DoModal();
	return FALSE;
}
