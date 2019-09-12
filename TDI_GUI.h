// TDI_GUI.h : main header file for the TDI_GUI application
//

#if !defined(AFX_TDI_GUI_H__AA4036A7_3EAD_4B2C_BE00_E8FFA49151FB__INCLUDED_)
#define AFX_TDI_GUI_H__AA4036A7_3EAD_4B2C_BE00_E8FFA49151FB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CTDI_GUIApp:
// See TDI_GUI.cpp for the implementation of this class
//

class CTDI_GUIApp : public CWinApp
{
public:
	CTDI_GUIApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTDI_GUIApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CTDI_GUIApp)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_TDI_GUI_H__AA4036A7_3EAD_4B2C_BE00_E8FFA49151FB__INCLUDED_)
