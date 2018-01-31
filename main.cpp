#define SRCDATADIR "/Users/maxscheurer/Software/libint_final/share/libint/2.4.2/basis/"
// #define LINALGWRAP_HAVE_ARMADILLO 1

#include <iostream>
#include <libint2.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <lazyten/SmallMatrix.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/eigensystem.hh>
#include <deque>
using std::deque;

#include "cppe.hpp"

#include "OptionParser.h"

using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using libint2::BasisSet;
using libint2::Atom;
using namespace std;
using libint2::read_dotxyz;
using libint2::make_point_charges;
using namespace lazyten;

#define THRESH 1e-10
#define EN_THRESH 1e-6


Options_t parseCommandline(int argc, char **argv, deque<string> &arguments) {
	OptionParser parser("%prog [options]");
	parser.add_string_option("xyzFile", "-x", "--xyz-filename",
		"The input geometry as xyz file.", "mol.xyz",
		"XYZ_FILENAME");
  parser.add_string_option("basisSet", "-b", "--basis",
  		"Basis set for the calculation.", "sto-3g",
  		"BASIS_SET");
	parser.add_int_option("maxSCFcycles", "-c", "", "Maximum number of SCF cycles"
                        , 50, "MAX_SCF_CYCLES");
	parser.add_bool_option("pe", "-pe", "--pol_embed", "Switch on PE calculation.", false);
	// parser.add_float_option("initialValue", "", "--initial-value",
	// 	"The initial value to use for the calculation.", 1.1);
	// parser.add_count_option("verbosity", "-v", "--verbose", "The verbosity level."
	// 	" Repeating this option increases verbosity.", 0);
	// parser.add_vector_option("stringList", "-s", "--string-item",
	// 	"A string to process.  This option may be specified more than once.",
	// 	"STRING");
	// parser.add_bool_option("fancy", "-f", "--fancy", "Switch on fancy formatting.");
	//cout << parser.get_by_option("-i").info() << endl;
	//cout << parser.get_by_name("help").info() << endl;
	//cout << parser.get_by_name("counter").info() << endl;
	Options_t options = parser.parse_args(argc, argv, arguments);
	//cout << parser.get_by_name("counter").info() << endl;
	// parser.print_info() << endl;
	//parser.print_help();
	//cout << "arguments: ";
	//for (auto el : arguments) {
		//cout << el << ", ";
	//}
	//cout << endl;
	return options;
}

double *computeMullikenCharges(const SmallMatrix<double> D, const SmallMatrix<double> S, vector<Atom> atoms,
                               const BasisSet basisSet) {
    SmallMatrix<double> P = D * S;
    auto shell2bf = basisSet.shell2bf();
    double *mullikenCharges = new double[atoms.size()];
    fill_n(mullikenCharges, atoms.size(), 0);
    vector<long> shell2atom = basisSet.shell2atom(atoms);
    for (vector<long>::iterator i = shell2atom.begin(); i != shell2atom.end(); ++i) {
        auto currentShell = i - shell2atom.begin();
        for (auto bf = 0; bf != basisSet[currentShell].size(); ++bf) {
            auto idx = bf + shell2bf[currentShell];
            mullikenCharges[*i] -= 2 * P(idx, idx);
        }
    }

    for (vector<Atom>::iterator at = atoms.begin(); at != atoms.end(); ++at) {
        auto idx = at - atoms.begin();
        mullikenCharges[idx] += (*at).atomic_number;
    }
    return mullikenCharges;
}

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

    const auto &res = s_engine.results();
    auto shell2bf = basisSet.shell2bf();

    for (int i = 0; i != basisSet.size(); ++i) {
        for (int k = 0; k != basisSet.size(); ++k) {
            s_engine.compute(basisSet[i], basisSet[k]);
            auto shellSizeI = basisSet[i].size();
            auto shellSizeK = basisSet[k].size();
            auto bf1 = shell2bf[i];
            auto bf2 = shell2bf[k];
            const auto *integral_shell = res[0];
            for (auto j = 0, intcount = 0; j < shellSizeI; ++j) {
                for (auto l = 0; l < shellSizeK; ++l, ++intcount) {
                    S(j + bf1, l + bf2) = integral_shell[intcount];
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
        for (int j = i + 1; j < natoms; ++j) {
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

BasisSet sortBasis(const BasisSet basisSet, const std::vector<Atom> atoms) {
	auto shell2atom = basisSet.shell2atom(atoms);
  auto atom2shell = basisSet.atom2shell(atoms);

  // loop over all atoms
	BasisSet newBas = basisSet;
	newBas.clear();
	for (auto a : atom2shell) {
		std::vector<std::vector<Shell>> newShellsForAtom(12);
		int size = 0;
		for (auto shellIdx : a) {
			Shell currentShell = basisSet[shellIdx];
			int max_l = 0;
      for (const auto& cr: currentShell.contr) {
        if (cr.l > max_l) {
          max_l = cr.l;
        }
      }
			if (max_l > size) {
				size = max_l;
			}
			newShellsForAtom[max_l].push_back(currentShell);
		}
		for (size_t i = 0; i <= size; ++i) {
			for (auto sh : newShellsForAtom[i]) {
				newBas.push_back(sh);
				// std::cout << sh << std::endl;
			}
		}
	}
	// basisSet.clear();
	// for (size_t i = 0; i < basisSet.size(); ++i) {
	// 	basisSet[i] = newBas[i];
	// }
	// basisSet = newBas;
	// std::cout << basisSet.size() << std::endl;
	return newBas;
}

int main(int argc, char **argv) {
    deque<string> arguments;
    Options_t options = parseCommandline(argc, argv, arguments);

		bool peCalculation = options["pe"].get_bool();

    libint2::initialize();
    std::cout << "libint initialized" << std::endl;
	  std::cout.precision(10);
//    read the geometry input
    string xyzfilename = options["xyzFile"].get_string();
    ifstream input_file(xyzfilename);
    vector<Atom> atoms = read_dotxyz(input_file);

//    assign a basis set & print number of BF
    BasisSet basisSet(options["basisSet"].get_string(), atoms);

		bool pure = false;
		for(auto s: basisSet) {
    	if(s.contr[0].pure) {
				pure = true;
			}
    }

    // sort the basis set by atom and angular momentum of shells
		basisSet = sortBasis(basisSet, atoms);
		basisSet.set_pure(pure);
		std::cout << "Pure: " << pure << std::endl;

    cout << "Number of basis functions: " << basisSet.nbf() << endl;
    int nbf = basisSet.nbf();
    int nshells = basisSet.size();

		std::cout << "Number of shells: " << nshells << std::endl;

		if (peCalculation) {
			std::string potfile("pehf_cpp.pot");
			call_pe_init(potfile, atoms, basisSet);
		}
		int embeddingEnergy = 0.0;

//    count number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
        nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;

//    compute the nuclear repulsion energy
    double nrep = computeNuclearRepulsionEnergy(atoms);

//    compute the basis set overlap
    SmallMatrix<double> S = computeOneBodyIntegrals(Operator::overlap, basisSet, atoms);

//    compute the core Hamiltonian
    SmallMatrix<double> NAttr = computeOneBodyIntegrals(Operator::nuclear, basisSet, atoms);
    SmallMatrix<double> T = computeOneBodyIntegrals(Operator::kinetic, basisSet, atoms);
//    cout << endl << endl;
//    cout << S << endl;
//    cout << endl << endl;
//    cout << NAttr << endl;
//    cout << endl << endl;
//    cout << T << endl;
//    cout << endl << endl;
    cout << "Nuclear repulsion energy is: " << nrep << endl;

    SmallMatrix<double> HCore = T + NAttr;
//    cout << "Core Hamiltonian: " << endl << HCore << endl;


    bool converged = false;
    int scfiteration = 0;
    SmallMatrix<double> finalD(nbf, nbf);
    SmallMatrix<double> F(nbf, nbf);
    double escf_old = 0.0;
    double escf = 0.0;
    double oldescf = 0.0;
    int maxSCFiterations = options["maxSCFcycles"].get_int();
    SmallMatrix<double> oldD(nbf, nbf);
    while (!converged) {
        Engine c_engine(Operator::coulomb,  // will compute overlap ints
                        basisSet.max_nprim(),    // max # of primitives in shells this engine will accept
                        basisSet.max_l(), 0         // max angular momentum of shells this engine will accept
        );
        const auto &res = c_engine.results();
        auto shell2bf = basisSet.shell2bf();

        cout << "SCF iteration: " << scfiteration << endl;
        if (scfiteration == 0) {
            F = HCore;
        }
        auto sol = eigensystem_hermitian(F, S);
        auto evals = sol.evalues();
        auto evecs = sol.evectors();
        auto v = evecs.subview(krims::range((unsigned long) ndocc));
        SmallMatrix<double> D = outer_prod_sum(v, v);
//        cout << "Density matrix" << endl << D << endl;
        double elEnergy = 0.0;
        for (int i = 0; i < nbf; ++i) {
            for (int j = 0; j < nbf; ++j) {
                elEnergy += D(i, j) * (HCore(i, j) + F(i, j));
            }
        }
				SmallMatrix<double> P = S*oldD*F-F*oldD*S;
				double pError = norm_frobenius_squared(P);
        double rmsd = norm_frobenius(D - oldD);
        oldD = D;
        escf_old = escf;
        escf = elEnergy + nrep;
        cout << scfiteration << " : Electronic energy: " << elEnergy << endl;
        cout << scfiteration << " : SCF energy: " << escf << endl;
				cout << "Pulay error: " << pError << endl;
        if (pError < THRESH && scfiteration && abs(escf-escf_old) < EN_THRESH) {
            converged = true;
            finalD = D;
	    			vector<double> energies;
	    			for (int o = 0; o < ndocc; ++o) {
		    			cout << evals[o] << endl;
							energies.push_back(evals[o]);
						}
            break;
        }
        oldescf = escf;

        SmallMatrix<double> G(nbf, nbf);
        for (int i = 0; i < nshells; ++i) {
            auto bf1 = shell2bf[i];
            for (int j = 0; j < nshells; ++j) {
                auto bf2 = shell2bf[j];
                for (int k = 0; k < nshells; ++k) {
                    auto bf3 = shell2bf[k];
                    for (int l = 0; l < nshells; ++l) {
                        auto bf4 = shell2bf[l];
                        c_engine.compute(basisSet[i], basisSet[j], basisSet[k], basisSet[l]);
                        const auto *integral_shell1 = res[0];

                        if (integral_shell1 == nullptr)
                            continue;

                        for (int s1 = 0, integralIndex = 0; s1 < basisSet[i].size(); ++s1) {
                            const auto i1 = s1 + bf1;
                            for (int s2 = 0; s2 < basisSet[j].size(); ++s2) {
                                const auto i2 = s2 + bf2;
//                                FNew(s1 + bf1, s2 + bf2) = 0;
                                for (int s3 = 0; s3 < basisSet[k].size(); ++s3) {
                                    const auto i3 = s3 + bf3;
                                    for (int s4 = 0; s4 < basisSet[l].size(); ++s4, ++integralIndex) {
                                        const auto i4 = s4 + bf4;
                                        G(i1, i2) += D(i3, i4) * 2.0 * integral_shell1[integralIndex];
                                        G(i1, i4) -= D(i2, i3) * integral_shell1[integralIndex];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
				std::vector<double> dmat;
				for (size_t i = 0; i < nbf; i++) {
					for (size_t j = 0; j < nbf; j++) {
						// cout << finalD(i,j) << endl;
						dmat.push_back(D(i,j));
					}
				}

				if (peCalculation) {
					SmallMatrix<double> fockPE = call_pe_fock(dmat);
	        F = HCore + G + fockPE;
				} else {
					F = HCore + G;
				}

        // cout << F.is_hermitian(1e-13) << " " << F.is_symmetric(1e-13) << endl;
        ++scfiteration;
        if (scfiteration >= maxSCFiterations) {
          std::cout << "Number of SCF iteractions exceeded: " << scfiteration << std::endl;
          break;
        }
    }

		// cout << finalD << endl;
		if (peCalculation) {
			std::vector<double> dmat;
			for (size_t i = 0; i < nbf; i++) {
				for (size_t j = 0; j < nbf; j++) {
					// cout << finalD(i,j) << endl;
					dmat.push_back(finalD(i,j));
				}
			}
			call_pe_energy(dmat);
		}

    // cout << "--- Computing Mulliken Point Charges ---" << endl;
    // double *mulliken = computeMullikenCharges(finalD, S, atoms, basisSet);
    // for (auto i = 0; i < atoms.size(); ++i) {
    //     cout << mulliken[i] << endl;
    // }
		std::cout << "--- SCF Summary ---" << std::endl;
		std::cout << "SCF Energy: " << escf << std::endl;
		// std::cout << "Embedding Energy: " << embeddingEnergy << std::endl;
		std::cout << "Total Energy: " << escf + embeddingEnergy << std::endl;
		std::cout << "-------------------" << std::endl;

//    don't use libint after this!
    libint2::finalize();
    cout << "--- Uh pick up a pancake  ---" << endl;
    return 0;
}
