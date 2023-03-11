# g++ `root-config --glibs --cflags` ./main.cxx -o main.exe

# ./main.exe ./test.root 1000 0 30000

for c in 100 200 400 800 1600;do
    hadd -f -j 8 ./fetch/raw${c}.root ./toyMC/raw_ctag${c}*.root
done