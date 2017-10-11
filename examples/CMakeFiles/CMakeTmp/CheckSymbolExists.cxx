/* */
#include <QtCore/qglobal.h>

int main(int argc, char** argv)
{
  (void)argv;
#ifndef Q_WS_WIN
  return ((int*)(&Q_WS_WIN))[argc];
#else
  (void)argc;
  return 0;
#endif
}

