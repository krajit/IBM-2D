/*********************************************2D FORCE**************************************************/

vectorField force(pointcloud_.size(),Zero);           ///Initialising Force at each vertex

forAll(area_,i)                                     ///For loop runs over each element
	{

        //Info<<"element"<<i<<endl;

        //Coordinates of reference triangle:
		vector  s0=pointcloud_[connectivity_matrix_[i][0]];   
		vector  s1=pointcloud_[connectivity_matrix_[i][1]];
		vector  s2=pointcloud_[connectivity_matrix_[i][2]];

        //Coordinates of deformed triangle:
		vector  x0=current_pointcloud_[connectivity_matrix_[i][0]];   
		vector  x1=current_pointcloud_[connectivity_matrix_[i][1]];
		vector  x2=current_pointcloud_[connectivity_matrix_[i][2]];

        //difference between coordinates of reference nodes:
        vector  ds1=s1-s0;
        vector  ds2=s2-s0;
     

        //difference between coordinates of reference nodes:
        vector  dx1=x1-x0;
        vector  dx2=x2-x0;


        //Deformation gradient:

        //First row:
        vector2D  a0( (dx1[0]*ds2[1]-dx2[0]*ds1[1]) , (-dx1[0]*ds2[0]+dx2[0]*ds1[0]) );
        
        //Second row:
        vector2D  a1( (dx1[1]*ds2[1] - dx2[1]*ds1[1]) , (-dx1[1]*ds2[0] + dx2[1]*ds1[0]) );
         
        //Tensor s:

        tensor2D s(ds1[0],ds2[0],ds1[1],ds2[1]);

        //Info<<"area"<<area_[i]<<endl;

        a0=(1.0/(det(s)))*a0;
        a1=(1.0/(det(s)))*a1;
         

        //Info<<"def_grad"<<a0<<endl;
        //Info<<"def_grad"<<a1<<endl;

        //Piola Kirchoff Stress:

        //For simplicity we set deformationgradient=Piola Kirchoff stress

        vector2D P0=a0;
        vector2D P1=a1;
        
        //Derivative of deformation gradient wrt position:
       vector2D dela_delx0( ds2[1]-ds1[1] , -ds2[0]+ds1[0]);
       vector2D dela_delx1( ds2[1] ,  -ds2[0] );
       vector2D dela_delx2(-ds1[1] ,  ds1[0] );


       dela_delx0 = -(1.0/(det(s)))*dela_delx0;
       dela_delx1 =  (1.0/(det(s)))*dela_delx1;
       dela_delx2 =  (1.0/(det(s)))*dela_delx2;
 

       vector f0( P0&dela_delx0, P1&dela_delx0, 0.0);
	   vector f1( P0&dela_delx1, P1&dela_delx1, 0.0);																	    
	   vector f2( P0&dela_delx2, P1&dela_delx2, 0.0);		 
																			 
																				   

       f0 = -area_[i]*f0;
       f1 = -area_[i]*f1;
       f2 = -area_[i]*f2;
       
       //Info<<"f0"<<f0<<endl;
       //Info<<"f1"<<f1<<endl;
       //Info<<"f2"<<f2<<endl;
      

       force[connectivity_matrix_[i][0]]=force[connectivity_matrix_[i][0]] + f0;
       force[connectivity_matrix_[i][1]]=force[connectivity_matrix_[i][1]] + f1;
	   force[connectivity_matrix_[i][2]]=force[connectivity_matrix_[i][2]] + f2;

	}  

Info<<force<<endl;

