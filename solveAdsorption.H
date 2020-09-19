//solve the ODE over deltaT/2

// update the varible
q_ODE = q; //current q
Y_ODE = Y; //current Y
Y0 = Y; //current Y

// loop over all the cells
//Info << "ODE: Solving for reaction " << endl;
forAll(Y_ODE,celli)
{
    // ****************invoke the ODE solver********************
    // Start time
    scalar xStart = 0;

    // Initial condition 
    //yStart[0] = y0; //debug
    //yStart[1] = 0.0; //debug

    yStart[0] = Y[celli]; //current Y
    yStart[1] = q[celli]; //current q

    // define where the time need
    //scalar xEnd = xStart + 10; //debug
    scalar xEnd = xStart + deltaT.value()/2; 
    scalarField y(yStart);

    // Initial ODE step size (does not affect the final solution)
    // but will affect the speed to some extend
    scalar dxEst =  1e-7;

    //- Solve the ODE system from xStart to xEnd, update the state
    //  and return an estimate for the next step in dxTry
    ode.derivatives(xStart, yStart, dyStart);
    odeSolver->solve(xStart, xEnd, y, dxEst);

    //Info<< "Numerical:  y = " << y << ", dxEst = " << dxEst << endl;

    // store the final results of ODE
    Y_ODE[celli] = y[0]; 
    q_ODE[celli] = y[1];
    // *****************end of the ODE solver*******************
}

//updated the q to keep the tab
q_ODE.max(0.0); //clip
q = q_ODE;

// overall reaction rate over the transport step
Y_ODE.max(0.0); //clip
//volScalarField dYdt = (Y_ODE-Y0)/deltaT/2;

// update the equilibirum concentration (postProcess)
forAll(qStar,celli)
{
    qStar[celli] = Langmuir_qMax.value()*Langmuir_b.value()*Y_ODE[celli]/
    (1.0+Langmuir_b.value()*Y_ODE[celli]);
}

// update
Y = Y_ODE;
Y.correctBoundaryConditions();