/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
     ibFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "labelListIOList.H"
#include "tensor2D.H" 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//TODO check whether the cloud point is in mesh or not



Foam::labelHashSet findNeighbourCells(vector probePoint, fvMesh& mesh)
{
    
     

      // find cell containing this point
      label celli = mesh.findCell(probePoint);

      // container for neighbours set by dumping the cell containing it
      labelHashSet neighbourCellSet(0);
      neighbourCellSet.set(celli);

     
						// number of layers
						int nLayers = 2;
						for (int n = 0; n < nLayers; n++)
						{
								// make a copy of marked cells
								labelHashSet markedNeighbours = neighbourCellSet;

								// loop over all marked cells
								forAllConstIter(labelHashSet, markedNeighbours,iter)
								{
											celli = iter.key();
											 
											// get points of celli
											labelList celliPoints = mesh.cellPoints()[celli];

											forAll(celliPoints,j)
											{
												// get neighbor cells of j th point
												labelList cellJNeighbours = mesh.pointCells()[celliPoints[j]];

												// append these cells in neighbourCellSet
												forAll(cellJNeighbours, k)
												{
													neighbourCellSet.set(cellJNeighbours[k]);
												 
												}
				
											}
								 }  

                         }
 
    return neighbourCellSet;

}
 
Foam::scalar weight(Foam::vector xcloud, Foam::vector xgrid)
{

	scalar hx=0.015625;  //TODO read from mesh
    scalar hy=0.015625;
	scalar phi1;
	scalar phi2;
     

///computing scaled distance between grid point and cloud point 
   scalar r1 = (1/hx)*(xgrid[0]-xcloud[0]);
   scalar r2 = (1/hy)*(xgrid[1]-xcloud[1]);
    

	 r1=fabs(r1);
	 r2=fabs(r2);
     
     //Info<<"r1 ="<<r1<<endl;
     //Info<<"r2 ="<<r2<<endl;
      

///computing phi1

	if (r1<=1.0)
	{
		phi1=(1.0/8.0)*(3.0-2.0*r1+Foam::sqrt(1.0+4.0*r1-4.0*r1*r1));
	} 
	else if(r1>=1.0 && r1<=2.0)
	{
		phi1=(1.0/8.0)*(5.0-2.0*r1-Foam::sqrt(-7.0+12.0*r1-4.0*r1*r1));
	} 
	else
	{
		phi1=0.0;
	} 


///computing phi2

	if (r2<=1.0)
	{
		phi2=(1.0/8.0)*(3.0-2.0*r2+Foam::sqrt(1.0+4.0*r2-4.0*r2*r2));
    } 
	else if(r2>=1.0 && r2<=2.0)
	{
		phi2=(1.0/8.0)*(5.0-2.0*r2-Foam::sqrt(-7.0+12.0*r2-4.0*r2*r2));
	} 
	else
	{
		phi2=0.0;
	}

 


 scalar w=(phi1/hx)*(phi2/hy);

 return w;

}







 
//------------------------


int main(int argc, char *argv[])
{
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"
     


    // * * * * * * * * * * * * * * GEOMETRY* * * * * * * * * * * * * * * * * * * //

      #include "area.H"                     //volume of solid mesh


        scalar hx=0.015625;                          //TODO read from mesh
        scalar hy=0.015625;
         
        scalar delta= 0.001;                     //TODO read time step

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * *INITIALISING VELOCITY FIELD * * * * * * * * * * * * * //
     
        forAll(mesh.C(), cellI)
			
			{
			scalar X = mesh.C()[cellI].component(0);
			//scalar Y = mesh.C()[cellI].component(1);

            scalar lambda=1;   //domain size

			 
		    scalar Uy = Foam::sin((2*(constant::mathematical::pi/lambda))*X);

                   Uy=0.3*Uy;
                    
			scalar Ux =  0.0;
			scalar Uz =  0.0;

			U[cellI] = vector(Ux,Uy,Uz);

			}

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  






    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

         #include "CourantNo.H"


     pointField temporary_pointcloud_(pointcloud_.size());                                                    //initialising temporary points X(n+1/2)
     

     Info<< "\n deform solid for first half time step\n" << endl;
         #include "deform1.H"                                                             //deform solid        X(n+1/2) 


     Info<< "\n Computes Lagrangian nodal force at half time step\n" << endl;      
         #include "Force.H"                                                               //compute solid force F(n+1/2) 



     Info<< "\n Spread the forces to the fluid gird\n" << endl;       
         #include "spread.H"                                                              //spread force        f(n+1/2) 
  


     Info<< "\n Solve for velocity and pressure at first half time step\n" << endl; 
         #include "predictor.H"                                                           //fluid solver        v(n+1/2) and p(n+1/2)   



     Info<< "\n Computing flux\n" << endl;                       
         #include "createPhi.H"                                                           //convective term     v(n+1/2)   



     Info<< "\n Deform the solid to end of the timestep\n" << endl;       
         #include "deform2.H"                                                              //deform solid        X(n+1)



     Info<< "\n Solve for velocity and pressure and complete the timestep \n" << endl;              
         #include "corrector.H"                                                            //fluid solver        v(n+1)   and p(n+1) 
   

        runTime.write();

        runTime.printExecutionTime(Info);

		Info<<"fluid solver end"<<endl; 
    }

    Info<< "Hey Don't worry I'm running\n" << endl;

    return 0;
}


// ************************************************************************* //
