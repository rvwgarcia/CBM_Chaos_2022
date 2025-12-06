#Opens Mandelbrot set data (kappa,Ps) pairs and executes cxmodelChaosICs.cpp with the parameter pair as inputs.
#Input arguments are Ps NT kmin kmax dk nSamples X0[0] X0[1]
#!/bin/bash

while read Ps kappa; do
./cxmodelChaosICs << EOF
$Ps 1e6 $kappa $kappa 1 1 0.31 0.1211
EOF
done < "LyapArgs.dat"