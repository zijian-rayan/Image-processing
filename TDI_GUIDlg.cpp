// TDI_GUIDlg.cpp : implementation file
//

#include "stdafx.h"
#include "TDI_GUI.h"
#include "TDI_GUIDlg.h"
#include "params.h"
#include "MainFrm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/******************* CAboutDlg dialog used for App About ******************/
class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

/****************************** CEditDrop *********************************/
BEGIN_MESSAGE_MAP(CEditDrop, CEdit)
	//{{AFX_MSG_MAP(CEditDrop)
	ON_WM_DROPFILES()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

CEditDrop::CEditDrop() {}

CEditDrop::~CEditDrop() {}

void CEditDrop::OnDropFiles(HDROP hDI) 
{
  char fname[1024];
  ::DragQueryFile(hDI,0,fname,1024);
  SetWindowText(fname);
  ::DragFinish(hDI);
	CEdit::OnDropFiles(hDI);
}

/************************* CTDI_GUIDlg dialog *******************************/
BEGIN_MESSAGE_MAP(CTDI_GUIDlg, CDialog)
	//{{AFX_MSG_MAP(CTDI_GUIDlg)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_RUN, OnRun)
	ON_CBN_SELCHANGE(IDC_COMMAND, OnSelChangeCommand)
	ON_CBN_SELCHANGE(IDC_CAT, OnSelChangeCat)
	//}}AFX_MSG_MAP
  ON_CONTROL_RANGE(BN_CLICKED, OFS_ID_BUTTON, OFS_ID_BUTTON+MAX_CTRLS-1, OnClickOpenSave)
  ON_CONTROL_RANGE(BN_CLICKED, OFS_ID_VIEW, OFS_ID_VIEW+MAX_CTRLS-1, OnClickView)
END_MESSAGE_MAP()

CTDI_GUIDlg::CTDI_GUIDlg(CWnd* pParent):CDialog(CTDI_GUIDlg::IDD, pParent),
nb_e(0), nb_s(0), nb_b(0), nb_c(0), idx_ctrl(0)
{
  memset(pe,0,sizeof(CEdit*)*MAX_CTRLS);
  memset(ps,0,sizeof(CStatic*)*MAX_CTRLS);
  memset(pb,0,sizeof(CButton*)*MAX_CTRLS);
  memset(pc,0,sizeof(CComboBox*)*MAX_CTRLS);
	//{{AFX_DATA_INIT(CTDI_GUIDlg)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
  font.CreatePointFont(80,"MS Sans Serif");
}

void CTDI_GUIDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CTDI_GUIDlg)
	DDX_Control(pDX, IDC_CAT, m_cat);
	DDX_Control(pDX, IDC_DESC, m_desc);
	DDX_Control(pDX, IDC_COMMAND, m_cmd);
	//}}AFX_DATA_MAP
}

BOOL CTDI_GUIDlg::OnInitDialog()
{
	CDialog::OnInitDialog();
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	
  for(int i=0; i<Func_Type::GetNbCat(); i++)
    m_cat.AddString(Func_Type::CatNames[i]);
  m_cat.AddString("Toutes");
  m_cat.SetCurSel(Func_Type::GetNbCat());
  OnSelChangeCat();

  return TRUE;  // return TRUE  unless you set the focus to a control
}

void CTDI_GUIDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

void CTDI_GUIDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

HCURSOR CTDI_GUIDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CTDI_GUIDlg::OnRun() 
{
  CString para_txt;
  char buf_msg[256];
  int nsel=m_cmd.GetCurSel(), plen;
  Func_Type& fct=Func_Type::func[nsel];
  for(DWORD i=0; i<fct.nb_param; i++)
  {
    Param_Type& par=fct.params[i];
    switch(par.type)
    {
    case 0:
    case 1:
      plen=pe[idx_ctrl[i]]->GetWindowTextLength();
      if(plen==0)
      {
        wsprintf(buf_msg,"Il manque le nom du fichier\npour \"%s\"",par.desc);
        MessageBox(buf_msg,"Erreur",MB_OK |MB_ICONERROR);
        return;
      }
      if(par.p_fname.fname) delete[] par.p_fname.fname;
      par.p_fname.fname=new char[plen+1];
      pe[idx_ctrl[i]]->GetWindowText(par.p_fname.fname,plen+1);
      break;
    case 2:
      pe[idx_ctrl[i]]->GetWindowText(para_txt);
      if(sscanf(para_txt,"%d",&par.p_int.val)!=1 || par.p_int.val<par.p_int.val_min || par.p_int.val>par.p_int.val_max)
      {
        sprintf(buf_msg,"Valeur de \"%s\" invalide !\nIl faut un entier entre %d et %d",par.desc,par.p_int.val_min,par.p_int.val_max);
        MessageBox(buf_msg,"Erreur",MB_OK |MB_ICONERROR);
        par.p_int.val=par.p_int.val_def;
        return;
      }
      break;
    case 3:
      pe[idx_ctrl[i]]->GetWindowText(para_txt);
      if(sscanf(para_txt,"%f",&par.p_float.val)!=1 || par.p_float.val<par.p_float.val_min || par.p_float.val>par.p_float.val_max)
      {
        sprintf(buf_msg,"Valeur de \"%s\" invalide !\nIl faut un réel entre %.2f et %.2f",par.desc,par.p_float.val_min,par.p_float.val_max);
        MessageBox(buf_msg,"Erreur",MB_OK |MB_ICONERROR);
        par.p_float.val=par.p_float.val_def;
        return;
      }
      break;
    case 4:
      par.p_sel.sel_cur=(BYTE)pc[idx_ctrl[i]]->GetCurSel();
      break;
    }
  }
  HANDLE hth=CreateThread(0,0,Func_Type::Launch,&fct,0,0);
}

void CTDI_GUIDlg::OnSelChangeCat() 
{
	int nsel=m_cat.GetCurSel();
  m_cmd.ResetContent();
  if(nsel==Func_Type::GetNbCat())
  {
    for(DWORD i=0; i<Func_Type::GetNbFunc(); i++)
      m_cmd.AddString(Func_Type::func[i].name);
  }
  else
  {
    for(DWORD i=0; i<Func_Type::GetNbFunc(); i++)
      if(Func_Type::func[i].cat==nsel) m_cmd.AddString(Func_Type::func[i].name);
  }
  m_cmd.SetCurSel(0); OnSelChangeCommand();
}

#define GUI_YMIN   70
#define GUI_YMAX  320
#define GUI_Y0MIN  (GUI_YMIN+15)
#define GUI_YSTEP  22
#define GUI_DY     18
#define GUI_XMIN   10
#define GUI_XMAX  470
#define GUI_X0MIN  15
#define GUI_X0MAX 170
#define GUI_X1MIN 175
#define GUI_X1MAX 365
#define GUI_X3MAX 260
#define GUI_X4MIN 370
#define GUI_X4MAX 415
#define GUI_X5MIN 420
#define GUI_X5MAX 465
#define GUI_DX     50
#define GUI_DY_EDIT 3
void CTDI_GUIDlg::OnSelChangeCommand() 
{
	CleanCtrls();
  int nsel=Func_Type::GetRealIdx(m_cat.GetCurSel(),m_cmd.GetCurSel());
  int ypos=GUI_Y0MIN, ypos_inp=0, ypos_out=0;
  Func_Type& fct=Func_Type::func[nsel];
  m_desc.SetWindowText(fct.desc);
  idx_ctrl=new DWORD[fct.nb_param];
  for(DWORD i=0; i<fct.nb_param; i++)
  {
    Param_Type& par=fct.params[i];
    if(par.type==1 && ypos_out==0)
    {
      ypos_inp=ypos;
      ypos_out=ypos_inp+5;
      ypos=ypos_inp+GUI_YSTEP;
    }
    ps[nb_s]=new CStatic;
    ps[nb_s]->Create(par.desc,WS_CHILD|WS_VISIBLE,CRect(GUI_X0MIN,ypos,GUI_X0MAX,ypos+GUI_DY),this);
    ps[nb_s++]->SetFont(&font);
    switch(par.type)
    {
    case 0:
      idx_ctrl[i]=nb_e;
      pe[nb_e]=new CEditDrop;
      pe[nb_e]->Create(ES_AUTOHSCROLL|ES_READONLY|WS_CHILD|WS_VISIBLE|ES_RIGHT,CRect(GUI_X1MIN,ypos-GUI_DY_EDIT,GUI_X1MAX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_EDIT);
      ::SetWindowLong(pe[nb_e]->m_hWnd,GWL_EXSTYLE,WS_EX_CLIENTEDGE|WS_EX_ACCEPTFILES|WS_EX_RIGHT);
      pe[nb_e]->SetWindowPos(0,0,0,0,0,SWP_FRAMECHANGED|SWP_NOMOVE|SWP_NOSIZE|SWP_NOREPOSITION);
      pe[nb_e++]->SetFont(&font);
      pb[nb_b]=new CButton;
      pb[nb_b]->Create("Ouvrir",WS_CHILD|WS_VISIBLE,CRect(GUI_X4MIN,ypos-GUI_DY_EDIT,GUI_X4MAX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_BUTTON);
      pb[nb_b++]->SetFont(&font);
      pb[nb_b]=new CButton;
      pb[nb_b]->Create("Voir",WS_CHILD|WS_VISIBLE,CRect(GUI_X5MIN,ypos-GUI_DY_EDIT,GUI_X5MAX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_VIEW);
      pb[nb_b++]->SetFont(&font);
      break;
    case 1:
      idx_ctrl[i]=nb_e;
      pe[nb_e]=new CEditDrop;
      pe[nb_e]->Create(ES_AUTOHSCROLL|WS_CHILD|WS_VISIBLE|ES_RIGHT,CRect(GUI_X1MIN,ypos-GUI_DY_EDIT,GUI_X1MAX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_EDIT);
      ::SetWindowLong(pe[nb_e]->m_hWnd,GWL_EXSTYLE,WS_EX_CLIENTEDGE|WS_EX_ACCEPTFILES|WS_EX_RIGHT);
      pe[nb_e]->SetWindowPos(0,0,0,0,0,SWP_FRAMECHANGED|SWP_NOMOVE|SWP_NOSIZE|SWP_NOREPOSITION);
      pe[nb_e++]->SetFont(&font);
      pb[nb_b]=new CButton;
      pb[nb_b]->Create("Sauver",WS_CHILD|WS_VISIBLE,CRect(GUI_X4MIN,ypos-GUI_DY_EDIT,GUI_X4MAX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_BUTTON);
      pb[nb_b++]->SetFont(&font);
      pb[nb_b]=new CButton;
      pb[nb_b]->Create("Voir",WS_CHILD|WS_VISIBLE,CRect(GUI_X5MIN,ypos-GUI_DY_EDIT,GUI_X5MAX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_VIEW);
      pb[nb_b++]->SetFont(&font);
      break;
    case 2:
      idx_ctrl[i]=nb_e;
      pe[nb_e]=new CEdit;
      pe[nb_e]->Create(ES_RIGHT|WS_CHILD|WS_VISIBLE,CRect(GUI_X1MIN,ypos-GUI_DY_EDIT,GUI_X1MIN+GUI_DX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_EDIT);
      ::SetWindowLong(pe[nb_e]->m_hWnd,GWL_EXSTYLE,WS_EX_CLIENTEDGE);
      pe[nb_e]->SetWindowPos(0,0,0,0,0,SWP_FRAMECHANGED|SWP_NOMOVE|SWP_NOSIZE|SWP_NOREPOSITION);
      {char buf[64]; sprintf(buf,"%d",par.p_int.val); pe[nb_e]->SetWindowText(buf);}
      pe[nb_e++]->SetFont(&font);
      break;
    case 3:
      idx_ctrl[i]=nb_e;
      pe[nb_e]=new CEdit;
      pe[nb_e]->Create(ES_RIGHT|WS_CHILD|WS_VISIBLE,CRect(GUI_X1MIN,ypos-GUI_DY_EDIT,GUI_X1MIN+GUI_DX,ypos+GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_EDIT);
      ::SetWindowLong(pe[nb_e]->m_hWnd,GWL_EXSTYLE,WS_EX_CLIENTEDGE);
      pe[nb_e]->SetWindowPos(0,0,0,0,0,SWP_FRAMECHANGED|SWP_NOMOVE|SWP_NOSIZE|SWP_NOREPOSITION);
      {char buf[64]; sprintf(buf,"%.3f",par.p_float.val); pe[nb_e]->SetWindowText(buf);}
      pe[nb_e++]->SetFont(&font);
      break;
    case 4:
      idx_ctrl[i]=nb_c;
      pc[nb_c]=new CComboBox;
      pc[nb_c]->Create(WS_CHILD|WS_VISIBLE|WS_VSCROLL|CBS_DROPDOWNLIST,CRect(GUI_X1MIN,ypos-GUI_DY_EDIT,GUI_X3MAX,ypos+5*GUI_DY-GUI_DY_EDIT),this,i+OFS_ID_COMBO);
      {
        char buf[64];
				char const *plim;
        const char *p1=par.p_sel.sel_str;
        while((plim=strchr(p1,'|'))!=0)
        {
          memcpy(buf,p1,plim-p1); buf[plim-p1]=0;
          pc[nb_c]->AddString(buf);
          p1=plim+1;
        }
        pc[nb_c]->SetCurSel(par.p_sel.sel_cur);
      }
      pc[nb_c++]->SetFont(&font);
      break;
    }
    ypos+=GUI_YSTEP;
  }
  pb[nb_b]=new CButton;
  pb[nb_b]->Create("Entrée",WS_CHILD|WS_VISIBLE|BS_GROUPBOX,CRect(GUI_XMIN,GUI_YMIN,GUI_XMAX,ypos_inp),this,0xffff);
  pb[nb_b++]->SetFont(&font);
  if(ypos_out)
  {
    pb[nb_b]=new CButton;
    pb[nb_b]->Create("Sortie",WS_CHILD|WS_VISIBLE|BS_GROUPBOX,CRect(GUI_XMIN,ypos_out,GUI_XMAX,ypos),this,0xffff);
    pb[nb_b++]->SetFont(&font);
  }
}

void CTDI_GUIDlg::OnClickOpenSave(UINT nID)
{
  DWORD idx=nID-OFS_ID_BUTTON;
  Param_Type& par=Func_Type::func[m_cmd.GetCurSel()].params[idx];
  if(par.type==0)
  {
    CFileDialog cf(TRUE,"*.*",0,OFN_FILEMUSTEXIST|OFN_PATHMUSTEXIST|OFN_EXTENSIONDIFFERENT,"Fichiers PGM|*.pgm|Fichiers PPM|*.ppm|Fichiers TXT|*.txt|");
    cf.m_ofn.Flags &= ~OFN_EXPLORER; // pour eviter le bug Windows/Explorer
    if(cf.DoModal()==IDOK)
    {
      pe[idx_ctrl[idx]]->SetWindowText(cf.GetPathName());
    }
//    else { DWORD res=::CommDlgExtendedError(); }
  }
  else if(par.type==1)
  {
    CFileDialog cf(FALSE,"*.*",0,OFN_OVERWRITEPROMPT|OFN_PATHMUSTEXIST|OFN_EXTENSIONDIFFERENT,"Fichiers PGM|*.pgm|Fichiers PPM|*.ppm|Fichiers TXT|*.txt|");
    cf.m_ofn.Flags &= ~OFN_EXPLORER; // pour eviter le bug Windows/Explorer
    if(cf.DoModal()==IDOK)
    {
      pe[idx_ctrl[idx]]->SetWindowText(cf.GetPathName());
    }
  }
}

void CTDI_GUIDlg::OnClickView(UINT nID)
{
  DWORD idx=nID-OFS_ID_VIEW;
  Param_Type& par=Func_Type::func[m_cmd.GetCurSel()].params[idx];
  CString fname;
  pe[idx_ctrl[idx]]->GetWindowText(fname);

  CMainFrame* pFrame = new CMainFrame(fname);
}


