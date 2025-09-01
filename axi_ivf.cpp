class INI	//initial volume fraction and level sets
{
	double *X,*Y,*CX,*CY,**Fa,**Phia;
	double Xc,Yc;	//center points of the ellipse
	double a,b;	//major and minor radius
	double H;	//height of liquid pool
	public:
		INI(double *x,double *y,double *cx,double *cy,double **fa,double **phia,double xc,double yc,double A,double B,double h);	//initialization of the input variables
		void inter();	//initial interface is defined here
		void VF();	//determine the initial volume fractions
		void ini_err();	//calculate initialization error in volume fractions
		void LS();	//determine the initial level set function
		void reinit(double **F,double **Phi,double **tag);	//determine the initial level set function from volume fractions
		void int_normal();	//export interface normal vector
};
INI::INI(double *x,double *y,double *cx,double *cy,double **fa,double **phia,double xc,double yc,double A,double B,double h)
{
	X=x; Y=y; CX=cx; CY=cy;
	Fa=fa; Phia=phia;
	Xc=xc; Yc=yc; a=A; b=B; H=h;
}
void INI::inter()
{
	const int nos=200;	//no of co-ordinates points
	double TH,pts[2][nos+1];
	for(int i=0;i<=nos;i++)	//generate ellipse centered at (Rc,Zc)
	{
		TH=2.0*M_PI*i/nos;
		pts[0][i]=Xc+a*cos(TH);	//R co-ordinates
		pts[1][i]=Yc+b*sin(TH);	//Z co-ordinates
	}
	ofstream p_out("ini_inter.dat");
	p_out<<"TITLE = \"INITIAL INTERFACE\""<<endl;
	p_out<<"FILETYPE = FULL"<<endl;
	p_out<<"VARIABLES = \"R\", \"Z\""<<endl;
	p_out<<"ZONE I = "<<nos+1<<", DATAPACKING = BLOCK"<<endl<<endl;
	for(int i=0;i<=nos;i++)	//print R coordinates
		p_out<<" "<<pts[0][i];
	p_out<<endl<<endl;
	for(int i=0;i<=nos;i++)	//print Z coordinates
		p_out<<" "<<pts[1][i];
	p_out.close();
	cout<<"INI: INITIAL INTERFACE FILE OUTPUT SUCCESSFUL"<<endl;
}
void INI::VF()
{
	double dist[4];	//stores the distance between the ellipse center and cell corner points
	double rad[4];	//stores the radius of ellipse at that particular angle
	double theta;	//calculate angle of line
	int p_out=0,p_in=0;	//flag variables
	double mYc;	//center of ellipse WRT shifted origin
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
//---------------------------DETERMINATION OF CELLS THROUGH WHICH INTERFACE PASSES-------------------
			dist[0]=sqrt(pow((Y[j-1]-Yc),2.0)+pow((X[i-1]-Xc),2.0));	//point a
			theta=atan2((Y[j-1]-Yc),(X[i-1]-Xc));
			rad[0]=a*b/sqrt(pow((b*cos(theta)),2.0)+pow((a*sin(theta)),2.0));

			dist[1]=sqrt(pow((Y[j]-Yc),2.0)+pow((X[i-1]-Xc),2.0));	//point b
			theta=atan2((Y[j]-Yc),(X[i-1]-Xc));
			rad[1]=a*b/sqrt(pow((b*cos(theta)),2.0)+pow((a*sin(theta)),2.0));

			dist[2]=sqrt(pow((Y[j]-Yc),2.0)+pow((X[i]-Xc),2.0));	//point c
			theta=atan2((Y[j]-Yc),(X[i]-Xc));
			rad[2]=a*b/sqrt(pow((b*cos(theta)),2.0)+pow((a*sin(theta)),2.0));

			dist[3]=sqrt(pow((Y[j-1]-Yc),2.0)+pow((X[i]-Xc),2.0));	//point d
			theta=atan2((Y[j-1]-Yc),(X[i]-Xc));
			rad[3]=a*b/sqrt(pow((b*cos(theta)),2.0)+pow((a*sin(theta)),2.0));

			for(int temp=0;temp<4;temp++)	//determining the cells containing the interface
			{
				if(abs(dist[temp]/rad[temp]-1.0)<=EPS) continue;	//point coincides with the interface
				if(dist[temp]>rad[temp]) p_out=1;	//point lies in region of light fluid
				else if(dist[temp]<rad[temp]) p_in=1;	//point lies in region of dark fluid
			}
			if((p_out==0)&&(p_in==1)) Fa[j][i]=1.0;	//entire cell lies in dark fluid
			else if((p_out==1)&&(p_in==0)) Fa[j][i]=0.0;	//entire cell lies in light fluid
			else	//interface passes through the cell
			{
//-------------------DETERMINATION OF THE ELLIPSE CENTER WRT THE NEW COORDINATE SYSTEM-------------------
				mYc=Yc-Y[j-1];	//center of ellipse WRT shifted origin
				double A,B,C,Z;	//discriminant at each cell face
				double x0=X[i-1],x1=100.0,x2=100.0,x3=X[i];
				double y0=0.0,y1=100.0,y2=100.0,y3=Y[j]-Y[j-1];	//intersecting coordinates (initialized with a flag value)
				double root1,root2;	//intersection points
				int intg_flag=0,f_flag[2],vol_flag;	//initialized with the default value
				f_flag[0]=f_flag[1]=5;	//cell face flag (default value)

				//f_flag is used to skip cell faces which are not intersected by the interface
				//intg_flag is used to mark the method of integration(for horizontal formulation(1) and for vertical formulation (0))

//-------------------DETERMINING THE ORIENTATION OF THE INTERFACE AND CALCULATING THE INTERSECTION POINTS-------------
				for(int temp=0;temp<4;temp++)	//determining if the curve passes through the corner points
				{
					switch(temp)
					{
						case 0:	//point a
							if(abs((b*b*pow((x0-Xc),2.0)+a*a*mYc*mYc)/(a*a*b*b)-1.0)<=EPS)
							{
								x1=x0; y1=y0;
								f_flag[0]=0; f_flag[1]=1;	//interface does not intersect west and south faces
								intg_flag=1;
							}
							break;
						case 1:	//point b
							if(abs((b*b*pow((x0-Xc),2.0)+a*a*pow((y3-mYc),2.0))/(a*a*b*b)-1.0)<=EPS)
							{
								x1=x0; y1=y3;
								f_flag[0]=0; f_flag[1]=2;	//interface does not intersect west and north faces
							}
							break;
						case 2:	//point c
							if(abs((b*b*pow((x3-Xc),2.0)+a*a*pow((y3-mYc),2.0))/(a*a*b*b)-1.0)<=EPS)
							{
								x2=x3; y2=y3;
								f_flag[0]=2; f_flag[1]=3;	//interface does not intersect north and east faces
							}
							break;
						case 3:	//point d
							if(abs((b*b*pow((x3-Xc),2.0)+a*a*mYc*mYc)/(a*a*b*b)-1.0)<=EPS)
							{
								x2=x3; y2=y0;
								f_flag[0]=1; f_flag[1]=3;	//interface does not intersect south and east faces
								intg_flag=1;
							}
							break;
					}
				}
				for(int temp=0;temp<4;temp++)	//determining the intersecting points between the curve and the cell faces
				{
					if((x1!=100.0)&&(x2!=100.0)) break;	//both the intersection points are obtained
					else if((temp==f_flag[0])||(temp==f_flag[1])) continue;	//skip faces which are not intersected by the interface
					switch(temp)	//setup the quadratic equation
					{
						case 0: //west face
							A=pow(a,2.0);
							B=-2.0*mYc*pow(a,2.0);
							C=pow((a*mYc),2.0)+pow((b*(x0-Xc)),2.0)-pow((a*b),2.0);
							break;
						case 1:	//south face
							A=pow(b,2.0);
							B=-2.0*Xc*pow(b,2.0);
							C=pow((b*Xc),2.0)+pow((a*mYc),2.0)-pow((a*b),2.0);
							break;
						case 2:	//north face
							A=pow(b,2.0);
							B=-2.0*Xc*pow(b,2.0);
							C=pow((b*Xc),2.0)+pow((a*(y3-mYc)),2.0)-pow((a*b),2.0);
							break;
						case 3:	//east face
							A=pow(a,2.0);
							B=-2.0*mYc*pow(a,2.0);
							C=pow((b*(x3-Xc)),2.0)+pow((a*mYc),2.0)-pow((a*b),2.0);
							break;
					}
					Z=pow(B,2.0)-4.0*A*C;	//discriminant calculation
					if(Z<0.0) continue;	//no intersection points on this face
					Z=-0.5*(B+SGN(B)*sqrt(Z));
					root1=Z/A;  root2=C/Z;	//calculate both the roots
					switch(temp)	//determine the required intersection point
					{
						case 0:	//west face
							if((root1>0.0)&&(root1<y3))	//root1 lies in the domain
							{
								x1=x0; y1=root1;	//(x2,y2) cannot lie in west face
							}
							else if((root2>0.0)&&(root2<y3))	//root2 lies in the domain
							{
								x1=x0; y1=root2;	//(x2,y2) cannot lie in west face
							}
							break;
						case 1:	//south face
							if((root1>x0)&&(root1<x3))	//root1 lies in the domain
							{
								if(x1==100.0) { x1=root1; y1=0.0; }
								else if(x2==100.0) { x2=root1; y2=0.0; }
								intg_flag=1;
							}
							else if((root2>x0)&&(root2<x3))	//root2 lies in the domain
							{
								if(x1==100.0) { x1=root2; y1=0.0; }
								else if(x2==100.0) { x2=root2; y2=0.0; }
								intg_flag=1;
							}
							break;
						case 2:	//north face
							if((root1>x0)&&(root1<x3))	//root1 lies in the domain
							{
								if(x1==100.0) { x1=root1; y1=y3; }
								else if(x2==100.0) { x2=root1; y2=y3; }
							}
							else if((root2>x0)&&(root2<x3))	//root2 lies in the domain
							{
								if(x1==100.0) { x1=root2; y1=y3; }
								else if(x2==100.0) { x2=root2; y2=y3; }
							}
							break;
						case 3:	//east face
							if((root1>0.0)&&(root1<y3))	//root1 lies in the domain
							{
								x2=x3; y2=root1;	//(x1,y1) cannot lie in east face
							}
							else if((root2>0.0)&&(root2<y3))	//root2 lies in the domain
							{
								x2=x3; y2=root2;	//(x1,y1) cannot lie in east face
							}
							break;
					}
				}	//end of temp loop
				if(intg_flag==0)	//arranging coordinates for vertical integration. Convention: x1<=x2
				{
					if(x1>=x2)	//interchange coordinates
					{
						A=x1; B=y1;	//variables are reused
						x1=x2; y1=y2;
						x2=A; y2=B;
					}
					y3=y0;
				}
				else if(intg_flag==1)	//arranging coordinates for horizontal integration. Convention: y1<=y2
				{
					if(y1>=y2)	//interchange coordinates
					{
						A=x1; B=y1;	//variables are reused
						x1=x2; y1=y2;
						x2=A; y2=B;
					}
					x3=x0;
				}
//-----------------CALCULATING THE VOLUME FRACTIONS-------------------------------
				if(intg_flag==1)	//volume fraction for horizontal configuration
				{
					vol_flag=(x2>Xc)?1:-1;	//resolve the plus or minus ambiguity
					A=x1*x1*(y1-y0)+x2*x2*(y3-y2)-x0*x0*(y3-y0)+Xc*Xc*(y2-y1)+a*a*(y2-y1)
						-a*a/(3.0*b*b)*(pow((y2-mYc),3.0)-pow((y1-mYc),3.0));
					if((b-abs(y2-mYc))<=EPS) B=b*b*asin(SGN(y2-mYc));
					else B=(y2-mYc)*pow((b*b-pow((y2-mYc),2.0)),0.5)+b*b*asin((y2-mYc)/b);
					if((b-abs(y1-mYc))<=EPS) C=b*b*asin(SGN(y1-mYc));
					else C=(y1-mYc)*pow((b*b-pow((y1-mYc),2.0)),0.5)+b*b*asin((y1-mYc)/b);
					Z=(A+vol_flag*a*Xc/b*(B-C))/(2.0*CX[i]*(X[i]-X[i-1])*(Y[j]-Y[j-1]));
				}
				else if(intg_flag==0)	//volume fraction for vertical configuration
				{
					vol_flag=(y2>mYc)?1:-1;	//resolve the plus or minus ambiguity
					A=0.5*(y1*(x1*x1-x0*x0)+y2*(x3*x3-x2*x2)+mYc*(x2*x2-x1*x1));
					if((a-abs(x1-Xc))<=EPS) B=-0.5*a*a*Xc*asin(SGN(x1-Xc));
					else B=pow((a*a-pow((x1-Xc),2.0)),1.5)/3.0-0.5*Xc*((x1-Xc)*pow((a*a-pow((x1-Xc),2.0)),0.5)+a*a*asin((x1-Xc)/a));
					if((a-abs(x2-Xc))<=EPS) C=-0.5*a*a*Xc*asin(SGN(x2-Xc));
					else C=pow((a*a-pow((x2-Xc),2.0)),1.5)/3.0-0.5*Xc*((x2-Xc)*pow((a*a-pow((x2-Xc),2.0)),0.5)+a*a*asin((x2-Xc)/a));
					Z=(A+vol_flag*b/a*(B-C))/(CX[i]*(X[i]-X[i-1])*(Y[j]-Y[j-1]));
				}
				if(vol_flag==1) Fa[j][i]=Z;
				else if(vol_flag==-1) Fa[j][i]=1.0-Z;
			}	//end of else
			p_out=0; p_in=0;	//reinitialization of the flag variables
		}
	}
	//------------------volume fractions of liquid pool---------------------
	for(int j=1;j<=J;j++)	//check location of pool interface
	{
		if(abs(Y[j]-H)<=EPS) p_in=j;
		else if((H>Y[j-1])&&(H<Y[j])) p_in=j;	//locate interface cell of pool
	}
	theta=(H-Y[p_in-1])/(Y[1]-Y[0]);	//volume fraction of interface cells
	for(int j=1;j<=p_in;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(j<p_in) Fa[j][i]=1.0;
			else if(j==p_in) Fa[j][i]+=theta;
		}
	}
}
void INI::ini_err()
{
	double V=0.0,Vt,err;
	Vt=Xc*M_PI*a*b;	//actual volume of torroidal ellipse
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			V+=Fa[j][i]*CX[i]*(X[i]-X[i-1])*(Y[j]-Y[j-1]);	//volume calculated from F
	err=abs(V-Vt)/Vt;
	cout<<"INI: INITIALIZATION ERROR = "<<err<<endl;
}
void INI::LS()
{
	double Dc,rad,theta;	//distance of a point from ellipse center, distance of a point on the curve from center and angle of the line WRT x axis
	double dx=X[2]-X[1];
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			Dc=sqrt(pow((CX[i]-Xc),2.0)+pow((CY[j]-Yc),2.0));	//calculate distance
			theta=atan2((CY[j]-Yc),(CX[i]-Xc));
			rad=a*b/sqrt(pow((a*sin(theta)),2.0)+pow((b*cos(theta)),2.0));
			if(abs(rad-Dc)<abs(H-CY[j])) Phia[j][i]=rad-Dc;
			else Phia[j][i]=H-CY[j];
		}
	}
}
void INI::reinit(double **F,double **Phi,double **tag)
{
	double h=MAX2((X[2]-X[1]),(Y[2]-Y[1]));
	double temp,a,b,gamma=0.5*h;	//gamma is the distance parameter
	double ray[2];
//--------------------INITIALIZATION SCHEME (TRACING ALGORITHM)--------------------------------------------------
	for(int j=1;j<=J;j++)	//tracing in horizontal direction
	{
		for(int i=1;i<=I;i++)
		{
			tag[j][i]=0;	//reinitialize tag value in the inner domain
		//-----------------VALUE OF RAY---------------------------------------------------
			if((abs(F[j][i])<1e-3)||(F[j][i]<0.0)) ray[1]=-1;	//empty cell
			else if(abs(1.0-F[j][i])<1e-3) ray[1]=1;	//completely filled cell
			else { ray[1]=0; tag[j][i]=1; }	//tag half-filled cell (case A)
		//-------------CELL TAGGING AND ADVANCE RAY---------------------------------------
			if((i!=1)&&(i!=I))	//for inner domain only
			{
				if((ray[0]*ray[1])==-1) { tag[j][i-1]=1; tag[j][i]=1; }	//case B
			}
			ray[0]=ray[1];	//advance the ray for the next cell
		}
	}
	for(int i=1;i<=I;i++)	//tracing in vertical direction
	{
		for(int j=1;j<=J;j++)
		{
		//-----------------VALUE OF RAY---------------------------------------------------
			if((abs(F[j][i])<1e-3)||(F[j][i]<0.0)) ray[1]=-1;	//empty cell
			else if(abs(1.0-F[j][i])<1e-3) ray[1]=1;	//completely filled cell
			else { ray[1]=0; tag[j][i]=1; }	//tag half-filled cell (case A)
		//-------------CELL TAGGING AND ADVANCE RAY---------------------------------------
			if((j!=1)&&(j!=J))	//for inner domain only
			{
				if((ray[0]*ray[1])==-1) { tag[j-1][i]=1; tag[j][i]=1; }	//case B
			}
			ray[0]=ray[1];	//advance the ray for the next cell
		}
	}
	for(int j=1;j<=J;j++)	//initialize level sets
	{
		for(int i=1;i<=I;i++)
		{
			if(tag[j][i]==1)
				Phi[j][i]=(2.0*F[j][i]-1.0)*gamma;	//calculate level set directly from volume fractions
			else
			{
				if((abs(F[j][i])<1e-3)||(F[j][i]<0.0)) Phi[j][i]=-100.0;	//unshaded region
				else Phi[j][i]=100.0;	//shaded region
			}
		}
	}
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(left and right boundaries)
	{
		Phi[j][0]=-100.0;	//ghost cells lie in the unshaded region
		Phi[j][I+1]=-100.0;
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(bottom and top boundaries)
	{
		Phi[0][i]=-100.0;
		Phi[J+1][i]=-100.0;
	}
//---------------SOLUTION OF DISCRETE EQUATIONS(including the ghost cells)-----------------------------
	for(int sweep=1,i_ini,j_ini,di,dj;sweep<=4;sweep++)	//Gauss-Siedel sweeps
	{
		switch(sweep)	//direction of each Gauss-Siedel sweep
		{
			case 1: j_ini=0; i_ini=0;
					dj=1; di=1;
					break;
			case 2: j_ini=0; i_ini=I+1;
					dj=1; di=-1;
					break;
			case 3: j_ini=J+1; i_ini=I+1;
					dj=-1; di=-1;
					break;
			case 4: j_ini=J+1; i_ini=0;
					dj=-1; di=1;
					break;
			default: break;
		}
		for(int j=j_ini;((j>=0)&&(j<=J+1));j+=dj)	//sweep the domain in the required direction (SMART LOOPS!)
		{
			for(int i=i_ini;((i>=0)&&(i<=I+1));i+=di)
			{
				if(tag[j][i]==1) continue;	//interface cells are not updated
				if(i==0) a=Phi[j][i+1];	//left boundary
				else if(i==I+1) a=Phi[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) a=MIN2(Phi[j][i+1],Phi[j][i-1]);
					else a=MAX2(Phi[j][i+1],Phi[j][i-1]);
				}
				if(j==0) b=Phi[j+1][i];	//bottom boundary
				else if(j==J+1) b=Phi[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) b=MIN2(Phi[j+1][i],Phi[j-1][i]);
					else b=MAX2(Phi[j+1][i],Phi[j-1][i]);
				}
				if(SGN(Phi[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-EPS) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MIN2(Phi[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-EPS) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MAX2(Phi[j][i],temp);
				}
			}
		}
	}
}
void INI::int_normal()
{
	double nx,ny;
	ofstream p_out("int_normal.dat");
	p_out<<"TITLE = \"INTERFACE NORMAL VECTOR\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"nx\",\"ny\""<<endl;
	p_out<<"ZONE I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=CELLCENTERED)"<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			nx=0.5*(Phia[j][i+1]-Phia[j][i-1])/(X[1]-X[0]);
			ny=0.5*(Phia[j+1][i]-Phia[j-1][i])/(Y[1]-Y[0]);
			nx/=(sqrt(nx*nx+ny*ny));
			p_out<<" "<<nx;
		}
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			nx=0.5*(Phia[j][i+1]-Phia[j][i-1])/(X[1]-X[0]);
			ny=0.5*(Phia[j+1][i]-Phia[j-1][i])/(Y[1]-Y[0]);
			ny/=(sqrt(nx*nx+ny*ny));
			p_out<<" "<<ny;
		}
		p_out<<endl;
	}
	p_out.close();
	cout<<"INI: INTERFACE NORMAL FILE OUTPUT SUCCESSFULL"<<endl;
}
