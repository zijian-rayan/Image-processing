// TDI_GUIDlg.h : header file
#if !defined(AFX_TDI_GUIDLG_H__15F8CDB0_2C52_46C9_87F4_59F2D07DD1C3__INCLUDED_)
#define AFX_TDI_GUIDLG_H__15F8CDB0_2C52_46C9_87F4_59F2D07DD1C3__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/*************************** CTDI_GUIDlg dialog ******************************/
#define MAX_CTRLS     10
#define OFS_ID_EDIT   20
#define OFS_ID_BUTTON 40
#define OFS_ID_COMBO  60
#define OFS_ID_VIEW   80
class CTDI_GUIDlg : public CDialog
{
  CEdit *pe[MAX_CTRLS];
  CStatic *ps[MAX_CTRLS];
  CButton *pb[MAX_CTRLS];
  CComboBox *pc[MAX_CTRLS];
  DWORD nb_e, nb_s, nb_b, nb_c, *idx_ctrl;
  void CleanCtrls()
  {
    if(idx_ctrl) {delete[] idx_ctrl; idx_ctrl=0;}
    while(nb_e) { --nb_e; if(pe[nb_e]) {delete pe[nb_e]; pe[nb_e]=0;}}
    while(nb_s) { --nb_s; if(ps[nb_s]) {delete ps[nb_s]; ps[nb_s]=0;}}
    while(nb_b) { --nb_b; if(pb[nb_b]) {delete pb[nb_b]; pb[nb_b]=0;}}
    while(nb_c) { --nb_c; if(pc[nb_c]) {delete pc[nb_c]; pc[nb_c]=0;}}
  }
  CFont font;
public:
	CTDI_GUIDlg(CWnd* pParent = NULL);
  ~CTDI_GUIDlg() {CleanCtrls();}

// Dialog Data
	//{{AFX_DATA(CTDI_GUIDlg)
	enum { IDD = IDD_TDI_GUI_DIALOG };
	CComboBox	m_cat;
	CStatic	m_desc;
	CComboBox	m_cmd;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTDI_GUIDlg)
	public:
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CTDI_GUIDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnRun();
	afx_msg void OnSelChangeCommand();
	afx_msg void OnSelChangeCat();
	//}}AFX_MSG
	afx_msg void OnClickOpenSave(UINT nID);
	afx_msg void OnClickView(UINT nID);
	DECLARE_MESSAGE_MAP()
};

/***************************** CEditDrop window ********************************/

class CEditDrop : public CEdit
{
// Construction
public:
	CEditDrop();

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CEditDrop)
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CEditDrop();

	// Generated message map functions
protected:
	//{{AFX_MSG(CEditDrop)
	afx_msg void OnDropFiles(HDROP hDropInfo);
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_TDI_GUIDLG_H__15F8CDB0_2C52_46C9_87F4_59F2D07DD1C3__INCLUDED_)
