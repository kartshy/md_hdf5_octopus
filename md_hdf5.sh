rm file.h5
cd ref
mpiexec -n 4 ./miniMD_cray
cd ../ext
mpiexec -n 4 ./miniAnalyser_cray

