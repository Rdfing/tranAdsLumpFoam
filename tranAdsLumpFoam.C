/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    tranAdsLumpFoam

Group
    grpBasicSolvers

Description

    Chemical transport and adsorption solver based on LFD and 
    langmuir isotherm.

    Adapted from combustion solver.
    
    Featured with time splitting of Strang scheme, stiff ODE solver, 

    Personally used for fixed bed simulation with fast reaction.

    Author: Randolph

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class myODE
:
public ODESystem
{
    // some parameters
    const scalar alpha_;
    const scalar beta_;
    const scalar qMax_;
    const scalar K_;

    public:

    myODE(const scalar& alpha_1, const scalar& beta_1, const scalar& qMax_1, const scalar& K_1)
    :
    ODESystem(),
    alpha_(alpha_1),
    beta_(beta_1),
    qMax_(qMax_1),
    K_(K_1)
    {}

    label nEqns() const
    {
        return 2;
    }

    void derivatives
    (
        const scalar x,
        const scalarField& y,
        scalarField& dydx
    ) const
    {
        dydx[0] = -alpha_*beta_*((qMax_*K_*y[0])/(1+K_*y[0])-y[1]);
        dydx[1] = beta_*((qMax_*K_*y[0])/(1+K_*y[0])-y[1]);
    }

    void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
        dfdx[0] = 0.0;
        dfdx[1] = 0.0;

        dfdy(0, 0) = -alpha_*beta_*qMax_*K_/Foam::pow(K_*y[0]+1,2);
        dfdy(0, 1) = alpha_*beta_;

        dfdy(1, 0) = beta_*qMax_*K_/Foam::pow(K_*y[0]+1,2);
        dfdy(1, 1) = -beta_;
    }
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Passive scalar transport equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    // constant parameter to simplify the expression
    const dimensionedScalar alpha = ((1-porosity)/porosity)*rhop;
    const dimensionedScalar beta = 15*De/(Rp*Rp);

    // Create the selected ODE system solver
    // ODE system parameter
    const scalar alpha_1 = alpha.value();
    const scalar beta_1 = beta.value(); 

    // langmuir isotherm parameters, pay attention to the unit!
    const scalar qMax_1 = Langmuir_qMax.value();
    const scalar K_1 = Langmuir_b.value();

    // Initial condition
    //const scalar y0 = 2e-20;

    // Create the ODE system
    myODE ode(alpha_1, beta_1, qMax_1, K_1);

    // construct the solver
    word ODESolverName ("rodas23");
    dictionary dict;
    dict.add("solver", ODESolverName);
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);

    // Required to store initial condition 
    scalarField yStart(ode.nEqns());
    // Required to store dydx
    scalarField dyStart(ode.nEqns());

    // some solving controls
    
    //odeSolver->maxSteps() = 20000; //10000 is the default
    //odeSolver->absTol() = 1e-12; //
    odeSolver->relTol() = 1e-4; //1e-4; is the default



    // Strang splitting scheme 
    // You can not use the CK time scheme!!!

    // transport equation time loop
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        const dimensionedScalar deltaT(mesh.time().deltaT()); // extract the time

        // initialize the varibles that pass the data to ODE
        volScalarField q_ODE = q; //current q
        volScalarField Y_ODE = Y; //current Y
        volScalarField Y0 = Y; //current Y

        while (simple.correctNonOrthogonal())
        {

            // Step 1 solve the ODE over deltaT/2
            Info << "Step 1 ODE: Solving for reaction" << endl;
            #include "solveAdsorption.H"

            // Step 2 solve the transport equation over deltaT
            Info << "Step 2 PDE: Solving for transport" << endl;
            fvScalarMatrix YEqn
            (
                fvm::ddt(Y)
              + fvm::div(phi, Y)
              - fvm::laplacian(DT, Y)
             //==
              //  dYdt
              // fvOptions(T)
            );

            YEqn.relax();
            fvOptions.constrain(YEqn);
            YEqn.solve();
            fvOptions.correct(Y);

            // Step 3 solve the ODE equation over deltaT/2
            Info << "Step 3 ODE: Solving for reaction" << endl;
            #include "solveAdsorption.H"
            
        }  

        runTime.write();
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
