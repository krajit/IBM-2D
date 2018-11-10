scalarField area_(connectivity_matrix_.size(),0.0) ;	 
{
    

	forAll(area_,i)
	{
		vector  p0=pointcloud_[connectivity_matrix_[i][0]];
	  	vector  p1=pointcloud_[connectivity_matrix_[i][1]];
	  	vector  p2=pointcloud_[connectivity_matrix_[i][2]];
	  	 

		area_[i]= (0.5)*mag( ((p1-p0)^(p2-p0)));
	}
}

