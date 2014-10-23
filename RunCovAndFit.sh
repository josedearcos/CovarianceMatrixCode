#!/bin/bash


g++ main.C -pthread  -pthread -stdlib=libc++ -m64 -I/Users/royal/root/include -L/Users/royal/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,/Users/royal/root/lib -stdlib=libc++ -lm -ldl
./a.out;
g++ mainFitter.C  -pthread -stdlib=libc++ -m64 -I/Users/royal/root/include -L/Users/royal/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,/Users/royal/root/lib -stdlib=libc++ -lm -ldl
./a.out;
#cd Plotter/
root plot_covariance_matrix.C;
root plot_flux_fractions.C;
root plot_systematics.C;


