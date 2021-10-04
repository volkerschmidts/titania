
#ifndef UNIT_STANIPC_H_
#define UNIT_STANIPC_H_

#ifndef MOLECULE_HPP_
class Molecule;
#endif


Molecule *ipcStan();

unsigned int ipc_stan_bonds();
unsigned int ipc_stan_angles();
unsigned int ipc_stan_torsions();
unsigned int ipc_stan_rdcs();

void outputStanMol(Molecule *);

#endif /* UNIT_STANIPC_H_ */
