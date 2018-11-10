

forAll(current_pointcloud_,i)
{

 
//COMPUTING VELOCITY AT THE PROBE POINT (Xn)

   vector probePoint=current_pointcloud_[i];               
   vector velocity(0.0,0.0,0.0);                                     //initialising probe point velocity

   //const volVectorField& U =mesh_.lookupObject<volVectorField>("U");   
   labelHashSet cell_neighbors=findNeighbourCells(probePoint,mesh);

   


	forAllConstIter(labelHashSet, cell_neighbors,iter) 
	{
		  label celli = iter.key();
		  vector cellicentre=mesh.cellCentres()[celli];
		  scalar W=weight(probePoint,cellicentre);

	 	  //Info<<"U of neighbour cell"<<celli<<"is"<<U[celli]<<endl;
		  //Info<<"Weight of neighbour cell"<<celli<<"is"<<W<<endl;

		   

		  velocity=velocity+(W*U[celli]*hx*hy);

	}
 
  

    //MOVE THE POINT FOR THE FIRST HALF TIME STEP AND SAVE INTO TEMPORARY POINT CLOUD//

    temporary_pointcloud_[i]=probePoint+ velocity*(delta/2);     
    //Info<<"current point"<<current_pointcloud_[i]<<endl;
    //Info<<"velocity of the cloud point"<<velocity<<endl;
 }


