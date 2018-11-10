//Reading fluid Information
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

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
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


Info<< "Reading field f\n" << endl;
volVectorField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);





#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

//Reading solid Information

Info<< "Reading cloud_of_points from constant\n" << endl;

pointIOField pointcloud_
    (
        IOobject
        (
            "cloud_of_points",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    );

Info<<"Reading currentcloudofpoints from time directory\n"<<endl;

pointIOField current_pointcloud_
    (
        IOobject
        (
            "currentcloudofpoints",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    );


 
 

Info<<"Reading connectivity matrix\n"<<endl;
labelListIOList connectivity_matrix_
   (
        IOobject
        (
            "connectivity",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE                    //NO_WRITE
        )
    );








