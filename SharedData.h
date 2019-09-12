#include <windows.h>

template <class CData> class CSharedData
{
  HANDLE hMapFile;
  CData* ptr;

public:
  CSharedData(const TCHAR* name, bool server=true): ptr(0), hMapFile(0)
  {
    if(server) hMapFile = CreateFileMapping(INVALID_HANDLE_VALUE, 0, PAGE_READWRITE, 0, sizeof(CData), name);
    else hMapFile = OpenFileMapping(FILE_MAP_ALL_ACCESS, FALSE, name);
    if(hMapFile == INVALID_HANDLE_VALUE) {hMapFile=0; return;}
    ptr=(CData*)MapViewOfFile(hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, sizeof(CData));
  }
  ~CSharedData()
  {
    if(ptr) UnmapViewOfFile(ptr);
    if(hMapFile) CloseHandle(hMapFile);
  }
  CData* operator->() {return ptr;}
  bool IsOK() { return ptr!=0; }
};

#define SH_NAME "MyShData"

class RobotState
{
public:
  static const int DIST_SENS_NB=8;
  static const int FLOOR_SENS_NB=3;
  static const int LIGHT_SENS_NB=8;
  static const int LED_NB=10;
  static const int IMG_SZ=3*52*39;
public:
  // up_stream (client->server)
  DWORD read_cam:1;
  DWORD read_wheel:1;
  DWORD read_dsensor:1;
  DWORD read_fsensor:1;
  DWORD read_lsensor:1;
  DWORD write_led:1;
  DWORD write_speed:8;
  double speed_l, speed_r;
  int led[LED_NB];
  // down_stream (server->client)
  double dist_l, dist_r;
  double dsens[DIST_SENS_NB];
  double fsens[FLOOR_SENS_NB];
  double lsens[LIGHT_SENS_NB];
  BYTE img[IMG_SZ];
  void Init() { memset(this,0,sizeof(RobotState)); }
  void ReadAll() { read_cam=read_wheel=read_dsensor=read_fsensor=read_lsensor=1; }
};

