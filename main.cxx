#include "ToyMC.h"

#include <TFile.h>
#include <string>

using namespace std;

int main(int argc,char** argv) {
    string name_outf = string(argv[1]);
    int ctag = atoi(argv[2]);
    int seed = atoi(argv[3]);
    int Nfills = atoi(argv[4]);

    TFile *fout = new TFile(name_outf.c_str(),"RECREATE");

    ToyMC generator(ctag,seed);

    generator.BookFileHists(fout);

    generator.Loop(Nfills);    

    generator.Finalize();

    return 1;    
}