/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "FieldField.H" //added
#include "OFstream.H" //added
// #include <armadillo> //added
#include <cmath> //added
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

double getG (int a, int b, int p, List<List<scalar> > A)
{
  double sumG;
  sumG = 0;
//   Info << "getG runned!"<< nl << endl;
  for (int gi=0; gi<p; gi++)
  {
    sumG += A[a][gi]*A[gi][b];
  }
  return sumG;
}

double getD (int a, int b, int p, List<List<scalar> > A)
{
  double sumD;
//   Info << "getD runned!"<< nl << endl;
  sumD = 0;
  if (a==b)
  {
    for (int di=0; di<p; di++)
    {
      return 1/sqrt(float(getG(a,b,p,A)));
    }    
  }
  else
  {
    return 0;
  }
}

double getR (int a, int b, int p, List<List<scalar> > A)
{
  double sumR;
//   Info << "getR runned!"<< nl << endl;
  sumR = 0;
  for (int t=0; t<p; t++)
  {
    for (int k=0; k<p; k++)
    {
      sumR+=(getD(a,k,p, A)*getG(k,t,p, A))*getD(t,b,p, A);

    }
  }
  return sumR;
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        p.storePrevIter();

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"

        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;

}


// ************************************************************************* //