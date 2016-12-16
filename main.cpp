#define SRCDATADIR "/usr/local/libint/2.3.0-beta.1/share/libint/2.3.0-beta.1/basis"
#define LINALGWRAP_HAVE_ARMADILLO 1
#include <iostream>
#include <libint2.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <linalgwrap/SmallMatrix.hh>

using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using libint2::BasisSet;
using libint2::Atom;
using namespace std;
using libint2::read_dotxyz;
using libint2::make_point_charges;
using namespace linalgwrap;

SmallMatrix<double> computeOneBodyIntegrals(Operator op, BasisSet basisSet, const vector<Atom> atoms) {
    int nbf = basisSet.nbf();
    int rows = nbf;
    int cols = nbf;
    SmallMatrix<double> S(rows, cols);

    Engine s_engine(op,  // will compute overlap ints
                    basisSet.max_nprim(),    // max # of primitives in shells this engine will accept
                    basisSet.max_l()         // max angular momentum of shells this engine will accept
    );
    s_engine.set_params(make_point_charges(atoms));

    const auto& res = s_engine.results();
    auto shell2bf = basisSet.shell2bf();
    double** overlapMatrix = new double*[7];
    for (int m = 0; m < 7; ++m) {
        overlapMatrix[m] = new double[7];
    }

    for (int i = 0; i != basisSet.size(); ++i) {
        for (int k = 0; k != basisSet.size() ; ++k) {
            s_engine.compute(basisSet[i], basisSet[k]);
            auto shellSizeI = basisSet[i].size();
            auto shellSizeK = basisSet[k].size();
            auto bf1 = shell2bf[i];
            auto bf2 = shell2bf[k];
            auto integral_shell = res[0];
            for (auto j = 0; j < shellSizeI; ++j) {
                for (auto l = 0; l < shellSizeK; ++l) {
//                    cout << j+bf1 << " " << l+bf2 << " " << integral_shell[j*shellSizeK+l] << endl;
                    S(j+bf1,l+bf2) = integral_shell[j*shellSizeK+l];
                }
            }
        }
    }
    return S;
}

void overlap(BasisSet basisSet, const vector<Atom> atoms) {
    Engine s_engine(Operator::overlap,  // will compute overlap ints
                    basisSet.max_nprim(),    // max # of primitives in shells this engine will accept
                    basisSet.max_l()         // max angular momentum of shells this engine will accept
    );
    s_engine.set_params(make_point_charges(atoms));
    const auto& res = s_engine.results();
    auto shell2bf = basisSet.shell2bf();
    double** overlapMatrix = new double*[7];
    for (int m = 0; m < 7; ++m) {
        overlapMatrix[m] = new double[7];
    }

    for (int i = 0; i != basisSet.size(); ++i) {
        for (int k = 0; k != basisSet.size() ; ++k) {
            s_engine.compute(basisSet[i], basisSet[k]);
            auto shellSizeI = basisSet[i].size();
            auto shellSizeK = basisSet[k].size();
            auto bf1 = shell2bf[i];
            auto bf2 = shell2bf[k];
            auto integral_shell = res[0];
            for (auto j = 0; j < shellSizeI; ++j) {
                for (auto l = 0; l < shellSizeK; ++l) {
//                    cout << j+bf1 << " " << l+bf2 << " " << integral_shell[j*shellSizeK+l] << endl;
                    overlapMatrix[j+bf1][l+bf2] = integral_shell[j*shellSizeK+l];
                }
            }
        }
    }


    for (int n = 0; n < 7; ++n) {
        for (int i = 0; i < 7; ++i) {
            cout << overlapMatrix[n][i] << " ";
        }
        cout << endl;
    }
}

int main() {
    libint2::initialize();

//    read the geometry input
    string xyzfilename = "/home/max/ClionProjects/QMax/water.xyz";
    ifstream input_file(xyzfilename);
    vector<Atom> atoms = read_dotxyz(input_file);

//    assign a basis set & print number of BF
    BasisSet basisSet("STO-3G", atoms);
    cout << "Number of basis functions: " << basisSet.nbf() << endl;

//    compute the core Hamiltonian
    SmallMatrix<double> S = computeOneBodyIntegrals(Operator::overlap, basisSet,atoms);
    SmallMatrix<double> NRep = computeOneBodyIntegrals(Operator::nuclear, basisSet,atoms);
    cout << endl << endl;
    cout << S << endl;
    cout << NRep << endl;
//    don't use libint after this!
    libint2::finalize();
    return 0;
}