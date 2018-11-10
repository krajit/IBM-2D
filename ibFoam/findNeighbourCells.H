 






//TODO check whether the cloud point is in mesh or not



Foam::labelHashSet findNeighbourCells(vector probePoint)
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
 
