#ifndef __MYVERTEX_HH__
#define __MYVERTEX_HH__
#include "TROOT.h"
#include "TObject.h"
using namespace std;
#include <vector>

class myVertex : public TObject {
public:

    myVertex() {
        ;
    }

    ~myVertex() {
        ;
    }

    //General



    float px, py, pz, z, mass, dz_Ver_match, Energy, mt, jetMass, eta_SC, phi_SC;
    float position_Rho;

     int tracksSize;
   float normalizedChi2;
    float ndof;
    unsigned int Num_Vertex;

    ClassDef(myVertex, 1)
};
#endif
