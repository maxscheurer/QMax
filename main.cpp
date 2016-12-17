#define SRCDATADIR "/usr/local/libint/2.3.0-beta.1/share/libint/2.3.0-beta.1/basis"
#define LINALGWRAP_HAVE_ARMADILLO 1
#include <iostream>
#include <libint2.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/eigensystem.hh>

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

double computeNuclearRepulsionEnergy(vector<Atom> atoms) {
    double nrep = 0.0;
    size_t natoms = atoms.size();
    for (int i = 0; i < natoms; ++i) {
        for (int j = i+1; j < natoms; ++j) {
            auto xij = atoms[i].x - atoms[j].x;
            auto yij = atoms[i].y - atoms[j].y;
            auto zij = atoms[i].z - atoms[j].z;
            auto r2 = xij * xij + yij * yij + zij * zij;
            auto r = sqrt(r2);
            nrep += atoms[i].atomic_number * atoms[j].atomic_number / r;
        }
    }
    return nrep;
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
    int nbf = basisSet.nbf();
    int nshells = basisSet.size();

//    count number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
        nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;

//    compute the nuclear repulsion energy
    double nrep = computeNuclearRepulsionEnergy(atoms);

//    compute the basis set overlap
    SmallMatrix<double> S = computeOneBodyIntegrals(Operator::overlap, basisSet,atoms);

//    compute the core Hamiltonian
    SmallMatrix<double> NAttr = computeOneBodyIntegrals(Operator::nuclear, basisSet,atoms);
    SmallMatrix<double> T = computeOneBodyIntegrals(Operator::kinetic, basisSet,atoms);
//    cout << endl << endl;
//    cout << S << endl;
//    cout << NAttr << endl;
    cout << "Nuclear repulsion energy is: " << nrep << endl;

    SmallMatrix<double> HCore = T + NAttr;
    cout << "Core Hamiltonian: " << endl << HCore << endl;

    bool converged = false;
    int scfiteration = 0;

    SmallMatrix<double> F(nbf,nbf);
    SmallMatrix<double> D(nbf,nbf);
    double escf = 0.0;
    while (!converged) {
        cout << "SCF iteration: " << scfiteration << endl;
        if (scfiteration == 0) {
            F = HCore;
            auto sol = eigensystem_hermitian(F, S);
            auto evals = sol.evalues();
            auto evecs = sol.evectors();
            auto v=evecs.subview(krims::range((unsigned  long) ndocc));
//            cout << evals[0] << endl;
//            cout << evecs << endl;
            D = outer_prod_sum(v,v);
            cout << "Density matrix" << endl << D << endl;
            SmallMatrix<double> EEl = D*(HCore+F);
            cout << EEl << endl;
            double elEnergy = 0.0;
            for (size_t col = 0; col < EEl.n_cols(); ++col) {
                for (size_t row = 0; row < EEl.n_rows(); ++row) {
                    elEnergy += EEl(row, col);
                }
            }
            escf = elEnergy + nrep;
            cout << "Initial electronic energy: " << elEnergy << endl;
            cout << "Initial SCF energy: " << escf << endl;

            SmallMatrix<double> FNew(nbf,nbf);
            for (int i = 0; i < nshells; ++i) {
                for (int j = 0; j < nshells; ++j) {
                    FNew(i,j)=HCore(i,j);
                    for (int k = 0; k < nshells; ++k) {
                        for (int l = 0; l < nshells; ++l) {
//                            FNew(i,j) += D(i,j)*
                        }
                    }
                }
            }

        } else {




            converged = true;
        }
        ++scfiteration;
    }

//    don't use libint after this!
    libint2::finalize();
    return 0;
}