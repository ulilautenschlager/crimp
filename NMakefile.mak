getopt.obj: src\getopt.c
  cl /TC /MT /W4 /c /O2 /EHsc /Fe"getopt.obj" /I"src" src\getopt.c

crimp.obj: src\crimp.c
  cl /TC /MT /W4 /c /O2 /EHsc /Fe"crimp.obj" /I"src" src\crimp.c

crimp.exe: getopt.obj crimp.obj
  link /out:crimp.exe getopt.obj crimp.obj
  del getopt.obj
  del crimp.obj
   
all: crimp.exe
