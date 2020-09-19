Info<< "Reading field Y\n" << endl;


volScalarField Y
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// q is the adsorbed phases concentration [kg adsorbed/kg adsorbent]
// so is dimensionless
volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField qStar // equilibirum q
(
    IOobject
    (
        "qStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// volScalarField Omega // reaction rate
// (
//     IOobject
//     (
//         "Omega",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );



Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading transportProperties\n" << endl;

// number of iteration in the ODE solving
//dimensionedScalar nODEstep_("nODEstep", dimless, transportProperties);

// bulk diffusivity
dimensionedScalar DT("DY", dimViscosity, transportProperties);
// intrapellet diffusivity
dimensionedScalar De("De", dimViscosity, transportProperties);

// adsorbent particle radius
dimensionedScalar Rp("Rp", dimLength, transportProperties);

// adsorbent density
dimensionedScalar rhop("rhop", dimDensity, transportProperties);
// porosity
dimensionedScalar porosity("porosity", dimless, transportProperties);

// Freundlich isotherm parameters
// dimensionedScalar Freundlich_K("Freundlich_K", dimless, transportProperties);
// dimensionedScalar Freundlich_n("Freundlich_n", dimless, transportProperties);

// Langmuir isotherm parameters
// the isotherm plot must be in the SI unit, C (kg/m3) and Qe (g/g)
dimensionedScalar Langmuir_qMax("Langmuir_qMax", dimless, transportProperties);
dimensionedScalar Langmuir_b("Langmuir_b", dimless, transportProperties);

#include "createPhi.H"

#include "createFvOptions.H"
