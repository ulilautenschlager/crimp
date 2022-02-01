getopt.obj: getopt_port\getopt.c
  cl /TC /MT /W4 /c /O2 /EHsc /Fe"getopt.obj" /I"getopt_port" getopt_port\getopt.c

crimp.obj: crimp.c
  cl /TC /MT /W4 /c /O2 /EHsc /Fe"crimp.obj" /I"getopt_port" crimp.c

crimp.exe: getopt.obj crimp.obj
  link /out:crimp.exe getopt.obj crimp.obj
  del getopt.obj
  del crimp.obj
   
all: crimp.exe
