project(PNG)
cmake_minimum_required(VERSION 3.9)

add_library(
  mypng
  png.c
  pngerror.c
  pngget.c
  pngmem.c
  pngpread.c
  pngread.c
  pngrio.c
  pngrtran.c
  pngrutil.c
  pngset.c
  pngwio.c
  pngwrite.c
  pngwtran.c
  pngwutil.c
  pngtrans.c)
