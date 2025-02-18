/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond_epslog.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondEpslog::BondEpslog(LAMMPS *_lmp) : Bond(_lmp)
{
  born_matrix_enable = 0;
}

/* ---------------------------------------------------------------------- */

BondEpslog::~BondEpslog()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(eps);
  }
}

/* ---------------------------------------------------------------------- */

void BondEpslog::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, rx, ebond, fx;

  ebond = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;


  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];

    rx=abs(delx);
    
    fx=(eps[type]/rx-0.5)/delx;

    
    if (eflag) ebond = eps[type]/rx+0.5*log(rx);
    
    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fx;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= fx;
    }

    
    if (evflag) ev_tally_xyz(i1, i2, nlocal, newton_bond, ebond, fx,0.,0., delx, 1.0, 1.0);
  }
}

/* ---------------------------------------------------------------------- */

void BondEpslog::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(eps, np1, "bond:eps");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondEpslog::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double eps_one = utils::numeric(FLERR, arg[1], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    eps[i] = eps_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondEpslog::equilibrium_distance(int i)
{
  // -eps/d^2+0.5/d == 0 <-> d=2*eps
  return 2*eps[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondEpslog::write_restart(FILE *fp)
{
  fwrite(&eps[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondEpslog::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &eps[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&eps[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondEpslog::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++) fprintf(fp, "%d %g\n", i, eps[i]);
}

/* ---------------------------------------------------------------------- */

double BondEpslog::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
{
  double r = sqrt(rsq);
  fforce=eps[type]/rsq-0.5/r;
  return eps[type]/r+0.5*log(r);;
}

/* ---------------------------------------------------------------------- */

// void BondEpslog::born_matrix(int type, double rsq, int /*i*/, int /*j*/, double &du, double &du2)
// {
//   double r = sqrt(rsq);
//   double dr = r - r0[type];
//   du2 = 0.0;
//   du = 0.0;
//   du2 = 2 * k[type];
//   if (r > 0.0) du = du2 * dr;
// }

/* ----------------------------------------------------------------------
   return ptr to internal members upon request
------------------------------------------------------------------------ */

void *BondEpslog::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "eps") == 0) return (void *) eps;
  return nullptr;
}
