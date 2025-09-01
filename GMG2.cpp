//GEOMETRIC MULTIGRID METHOD FOR CELL-CENTERED DISCRETIZATION
//UPTO 5-GRID CYCLE
//PIECEWISE CONSTANT INTERPOLATION IS USED FOR THE PROLONGATION OPERATOR
//BILINEAR INTERPOLATION IS USED FOR THE RESTRICTION OPERATOR
//OPERATORS ARE DESIGNED FOR PURE DIRICHLET BOUNDARY CONDITION
//GALERKIN COARSE GRID APPROXIMATION OF THE COEFFICIENT MATRICES
inline int IND(int I,int i,int j) { return ((j-1)*I+i-1); }	//matrix indices to vector index
class STENCIL16	//16 point restriction stencil (bilinear interpolation)
{
	public:
		double r[16];
		STENCIL16()	//initialization
		{
			for(int i=0;i<16;i++) r[i]=0.0;
		}
		void create(const int I,const int J,int i,int j);	//create restriction stencil
};
void STENCIL16::create(const int I,const int J,int i,int j)
{
	for(int i=0;i<16;i++) r[i]=0.0;	//reinitialization
	int e=1,w=1,n=1,s=1;	//boundary flags set to default value
	if(i==2) w=0;	//boundary modifications in X direction
	else if(i==I) e=0;
	if(j==2) s=0;	//boundary modifications in Y direction
	else if(j==J) n=0;
	r[0]=s*w/16.0; r[1]=(2+w)*s/16.0; r[2]=(2+e)*s/16.0; r[3]=s*e/16.0;
	r[4]=(2+s)*w/16.0; r[5]=(2+w)*(2+s)/16.0; r[6]=(2+e)*(2+s)/16.0; r[7]=(2+s)*e/16.0;
	r[8]=(2+n)*w/16.0; r[9]=(2+w)*(2+n)/16.0; r[10]=(2+e)*(2+n)/16.0; r[11]=(2+n)*e/16.0;
	r[12]=n*w/16.0; r[13]=(2+w)*n/16.0; r[14]=(2+e)*n/16.0; r[15]=n*e/16.0;
}
class STENCIL9	//9 point stencil of coefficient matrix
{
	public:
		double a[9];
		STENCIL9()	//initialization
		{
			for(int i=0;i<9;i++) a[i]=0.0;
		}
		void create(const int I,const int J,int i,int j,int *C,int *R,double *A,double W=1.0);	//create weighted stencil from coefficient matrix
};
void STENCIL9::create(const int I,const int J,int i,int j,int *C,int *R,double *A,double W)
{
	for(int i=0;i<9;i++) a[i]=0.0;	//reinitialization
	if((i<1)||(i>I)||(j<1)||(j>J)) return;	//check index boundedness
	for(int t=R[IND(I,i,j)];t<R[IND(I,i,j)+1];t++)
	{
		if(C[t]==IND(I,i-1,j-1)) a[0]=W*A[t];
		else if(C[t]==IND(I,i,j-1)) a[1]=W*A[t];
		else if(C[t]==IND(I,i+1,j-1)) a[2]=W*A[t];
		else if(C[t]==IND(I,i-1,j)) a[3]=W*A[t];
		else if(C[t]==IND(I,i,j)) a[4]=W*A[t];
		else if(C[t]==IND(I,i+1,j)) a[5]=W*A[t];
		else if(C[t]==IND(I,i-1,j+1)) a[6]=W*A[t];
		else if(C[t]==IND(I,i,j+1)) a[7]=W*A[t];
		else if(C[t]==IND(I,i+1,j+1)) a[8]=W*A[t];
	}
}
class GMG2
{
	protected:
			double *A0; int *C0,*R0;	//CSR arrays for the coefficient matrix (level 0)
			double *X0,*b0,*r0;	//solution, RHS, and residual vectors (level 0)
			double *A1; int *C1,*R1;	//CSR arrays for the coefficient matrix (level 1)
			double *X1,*b1,*r1;	//solution, RHS, and residual vectors (level 1)
			double *A2; int *C2,*R2;	//CSR arrays for the coefficient matrix (level 2)
			double *X2,*b2,*r2;	//solution, RHS, and residual vectors (level 2)
			double *A3; int *C3,*R3;	//CSR arrays for the coefficient matrix (level 3)
			double *X3,*b3,*r3;	//solution, RHS, and residual vectors (level 3)
			double *A4; int *C4,*R4;	//CSR arrays for the coefficient matrix (level 4)
			double *X4,*b4,*r4;	//solution, RHS, and residual vectors (level 4)

			void setup(const int I,const int J,int *C0,int *R0,double *A0,int *C1,int *R1,double *A1);	//setup the multigrid system (Galerkin coarse grid approximation)
			double solve_SSOR(const int NITER,const int N,int *C,int *R,double *A,double *X,double *b,double *r);	//SSOR solver
			void restrict(const int I,const int J,double *r0,double *b1);	//restriction operator
			void prolong(const int I,const int J,double *X1,double *X0);	//prolongation operator
			void reini();	//reinitialize the vectors
	public:
			GMG2(); ~GMG2();	//memory management
};
GMG2::GMG2()
{
	int NNZ,N;	//no of nonzero elements in the sparse matrix and no of elements in the solution vector
	//---------level 0 (finest level)------------------
	NNZ=5*(I-2)*(J-2)+8*(J-2)+8*(I-2)+12;
	N=I*J;
	C0=new int[NNZ];
	A0=new double[NNZ];
	R0=new int[N+1];
	X0=new double[N];
	b0=new double[N];
	r0=new double[N];
	//---------level 1 (standard cell-centered coarsening)------------
	NNZ=9*(I/2-2)*(J/2-2)+12*(J/2-2)+12*(I/2-2)+16;	//GCA 9 point stencil
	N/=4;
	C1=new int[NNZ];
	A1=new double[NNZ];
	R1=new int[N+1];
	X1=new double[N];
	b1=new double[N];
	r1=new double[N];
	//---------level 2 (standard cell-centered coarsening)------------
	NNZ=9*(I/4-2)*(J/4-2)+12*(J/4-2)+12*(I/4-2)+16;	//GCA 9 point stencil
	N/=4;
	C2=new int[NNZ];
	A2=new double[NNZ];
	R2=new int[N+1];
	X2=new double[N];
	b2=new double[N];
	r2=new double[N];
	//---------level 3 (standard cell-centered coarsening)------------
	NNZ=9*(I/8-2)*(J/8-2)+12*(J/8-2)+12*(I/8-2)+16;	//GCA 9 point stencil
	N/=4;
	C3=new int[NNZ];
	A3=new double[NNZ];
	R3=new int[N+1];
	X3=new double[N];
	b3=new double[N];
	r3=new double[N];
	//---------level 4 (standard cell-centered coarsening)------------
	NNZ=9*(I/16-2)*(J/16-2)+12*(J/16-2)+12*(I/16-2)+16;	//GCA 9 point stencil
	N/=4;
	C4=new int[NNZ];
	A4=new double[NNZ];
	R4=new int[N+1];
	X4=new double[N];
	b4=new double[N];
	r4=new double[N];
	cout<<"GMG2: MEMORY ALLOCATED"<<endl;
}
GMG2::~GMG2()
{
	//---------level 0------------------
	delete[] C0;
	delete[] A0;
	delete[] R0;
	delete[] X0;
	delete[] b0;
	delete[] r0;
	//---------level 1------------------
	delete[] C1;
	delete[] A1;
	delete[] R1;
	delete[] X1;
	delete[] b1;
	delete[] r1;
	//---------level 2------------------
	delete[] C2;
	delete[] A2;
	delete[] R2;
	delete[] X2;
	delete[] b2;
	delete[] r2;
	//---------level 3------------------
	delete[] C3;
	delete[] A3;
	delete[] R3;
	delete[] X3;
	delete[] b3;
	delete[] r3;
	//---------level 4------------------
	delete[] C4;
	delete[] A4;
	delete[] R4;
	delete[] X4;
	delete[] b4;
	delete[] r4;
	cout<<"GMG2: MEMORY RELEASED"<<endl;
}
void GMG2::setup(const int I,const int J,int *C0,int *R0,double *A0,int *C1,int *R1,double *A1)
{
	STENCIL9 A[16],Ab;
	STENCIL16 R;
	int cnt=0; R1[cnt]=0;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			R.create(2*I,2*J,2*i,2*j);
			for(int m=-2,t=0;m<=1;m++)	//create 16 fine grid stencils
			{
				for(int n=-2;n<=1;n++)
				{
					A[t].create(2*I,2*J,(2*i+n),(2*j+m),C0,R0,A0,R.r[t]);
					t++;
				}
			}

			Ab.a[0]= A[0].a[0]+A[0].a[1]+A[0].a[3]+A[0].a[4]		//stencil A[0]
					+A[1].a[0]+A[1].a[3]							//stencil A[1]
					+A[4].a[0]+A[4].a[1]							//stencil A[4]
					+A[5].a[0];										//stencil A[5]

			Ab.a[1]= A[0].a[2]+A[0].a[5]							//stencil A[0]
					+A[1].a[1]+A[1].a[2]+A[1].a[4]+A[1].a[5]		//stencil A[1]
					+A[2].a[0]+A[2].a[1]+A[2].a[3]+A[2].a[4]		//stencil A[2]
					+A[3].a[0]+A[3].a[3]							//stencil A[3]
					+A[4].a[2]										//stencil A[4]
					+A[5].a[1]+A[5].a[2]							//stencil A[5]
					+A[6].a[0]+A[6].a[1]							//stencil A[6]
					+A[7].a[0];										//stencil A[7]

			Ab.a[2]= A[2].a[2]+A[2].a[5]							//stencil A[2]
					+A[3].a[1]+A[3].a[2]+A[3].a[4]+A[3].a[5]		//stencil A[3]
					+A[6].a[2]										//stencil A[6]
					+A[7].a[1]+A[7].a[2];							//stencil A[7]

			Ab.a[3]= A[0].a[6]+A[0].a[7]							//stencil A[0]
					+A[1].a[6]										//stencil A[1]
					+A[4].a[3]+A[4].a[4]+A[4].a[6]+A[4].a[7]		//stencil A[4]
					+A[5].a[3]+A[5].a[6]							//stencil A[5]
					+A[8].a[0]+A[8].a[1]+A[8].a[3]+A[8].a[4]		//stencil A[8]
					+A[9].a[0]+A[9].a[3]							//stencil A[9]
					+A[12].a[0]+A[12].a[1]							//stencil A[12]
					+A[13].a[0];									//stencil A[13]

			Ab.a[4]= A[0].a[8]										//stencil A[0]
					+A[1].a[7]+A[1].a[8]							//stencil A[1]
					+A[2].a[6]+A[2].a[7]							//stencil A[2]
					+A[3].a[6]										//stencil A[3]
					+A[4].a[5]+A[4].a[8]							//stencil A[4]
					+A[5].a[4]+A[5].a[5]+A[5].a[7]+A[5].a[8]		//stencil A[5]
					+A[6].a[3]+A[6].a[4]+A[6].a[6]+A[6].a[7]		//stencil A[6]
					+A[7].a[3]+A[7].a[6]							//stencil A[7]
					+A[8].a[2]+A[8].a[5]							//stencil A[8]
					+A[9].a[1]+A[9].a[2]+A[9].a[4]+A[9].a[5]		//stencil A[9]
					+A[10].a[0]+A[10].a[1]+A[10].a[3]+A[10].a[4]	//stencil A[10]
					+A[11].a[0]+A[11].a[3]							//stencil A[11]
					+A[12].a[2]										//stencil A[12]
					+A[13].a[1]+A[13].a[2]							//stencil A[13]
					+A[14].a[0]+A[14].a[1]							//stencil A[14]
					+A[15].a[0];									//stencil A[15]

			Ab.a[5]= A[2].a[8]										//stencil A[2]
					+A[3].a[7]+A[3].a[8]							//stencil A[3]
					+A[6].a[5]+A[6].a[8]							//stencil A[6]
					+A[7].a[4]+A[7].a[5]+A[7].a[7]+A[7].a[8]		//stencil A[7]
					+A[10].a[2]+A[10].a[5]							//stencil A[10]
					+A[11].a[1]+A[11].a[2]+A[11].a[4]+A[11].a[5]	//stencil A[11]
					+A[14].a[2]										//stencil A[14]
					+A[15].a[1]+A[15].a[2];							//stencil A[15]

			Ab.a[6]= A[8].a[6]+A[8].a[7]							//stencil A[8]
					+A[9].a[6]										//stencil A[9]
					+A[12].a[3]+A[12].a[4]+A[12].a[6]+A[12].a[7]	//stencil A[12]
					+A[13].a[3]+A[13].a[6];							//stencil A[13]

			Ab.a[7]= A[8].a[8]										//stencil A[8]
					+A[9].a[7]+A[9].a[8]							//stencil A[9]
					+A[10].a[6]+A[10].a[7]							//stencil A[10]
					+A[11].a[6]										//stencil A[11]
					+A[12].a[5]+A[12].a[8]							//stencil A[12]
					+A[13].a[4]+A[13].a[5]+A[13].a[7]+A[13].a[8]	//stencil A[13]
					+A[14].a[3]+A[14].a[4]+A[14].a[6]+A[14].a[7]	//stencil A[14]
					+A[15].a[3]+A[15].a[6];							//stencil A[15]

			Ab.a[8]= A[10].a[8]										//stencil A[10]
					+A[11].a[7]+A[11].a[8]							//stencil A[11]
					+A[14].a[5]+A[14].a[8]							//stencil A[14]
					+A[15].a[4]+A[15].a[5]+A[15].a[7]+A[15].a[8];	//stencil A[15]

			for(int m=-1,t=0;m<=1;m++)	//update coefficient matrix
			{
				for(int n=-1;n<=1;n++)
				{
					if(abs(Ab.a[t])>SMALL)
					{
						C1[cnt]=(j+m-1)*I+i+n-1;
						A1[cnt]=Ab.a[t];
						cnt++;
					}
					t++;
				}
			}
			R1[(j-1)*I+i]=cnt;
		}
	}
}
double GMG2::solve_SSOR(const int NITER,const int N,int *C,int *R,double *A,double *X,double *b,double *r)
{
	double sum,res,diag;
	int cnt=0;	//iteration counter
	do
	{
		for(int i=0;i<N;i++)	//matrix vector multiplication, forward sweep
		{
			sum=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++)
			{
				if(C[d]==i) diag=A[d];
				else sum+=A[d]*X[C[d]];
			}
			X[i]=(b[i]-sum)/diag;
		}
		for(int i=N-1;i>=0;i--)	//matrix vector multiplication, backward sweep
		{
			sum=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++)
			{
				if(C[d]==i) diag=A[d];
				else sum+=A[d]*X[C[d]];
			}
			X[i]=(b[i]-sum)/diag;
		}
		res=0.0;	//reinitialization
		for(int i=0;i<N;i++)	//residual vector and 2-norm
		{
			sum=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++) sum+=A[d]*X[C[d]];
			r[i]=b[i]-sum;
			res+=pow(r[i],2.0);
		}
		res=sqrt(res);
		cnt++;
	}
	while((res>TOL)&&(cnt<NITER));
	return res;
}
void GMG2::restrict(const int I,const int J,double *r0,double *b1)
{
	for(int j=2;j<J;j++)	//inner domain
		for(int i=2;i<I;i++)
			b1[IND(I,i,j)]=0.0625*(r0[IND(2*I,2*i-2,2*j+1)]+3.0*r0[IND(2*I,2*i-1,2*j+1)]+3.0*r0[IND(2*I,2*i,2*j+1)]+r0[IND(2*I,2*i+1,2*j+1)]+
							3.0*r0[IND(2*I,2*i-2,2*j)]+9.0*r0[IND(2*I,2*i-1,2*j)]+9.0*r0[IND(2*I,2*i,2*j)]+3.0*r0[IND(2*I,2*i+1,2*j)]+
							3.0*r0[IND(2*I,2*i-2,2*j-1)]+9.0*r0[IND(2*I,2*i-1,2*j-1)]+9.0*r0[IND(2*I,2*i,2*j-1)]+3.0*r0[IND(2*I,2*i+1,2*j-1)]+
							r0[IND(2*I,2*i-2,2*j-2)]+3.0*r0[IND(2*I,2*i-1,2*j-2)]+3.0*r0[IND(2*I,2*i,2*j-2)]+r0[IND(2*I,2*i+1,2*j-2)]);
	for(int j=2;j<J;j++)	//left and right boundary
	{
		b1[IND(I,1,j)]=0.0625*(2.0*r0[IND(2*I,1,2*j+1)]+3.0*r0[IND(2*I,2,2*j+1)]+r0[IND(2*I,3,2*j+1)]+
						6.0*r0[IND(2*I,1,2*j)]+9.0*r0[IND(2*I,2,2*j)]+3.0*r0[IND(2*I,3,2*j)]+
						6.0*r0[IND(2*I,1,2*j-1)]+9.0*r0[IND(2*I,2,2*j-1)]+3.0*r0[IND(2*I,3,2*j-1)]+
						2.0*r0[IND(2*I,1,2*j-2)]+3.0*r0[IND(2*I,2,2*j-2)]+r0[IND(2*I,3,2*j-2)]);
		b1[IND(I,I,j)]=0.0625*(r0[IND(2*I,2*I-2,2*j+1)]+3.0*r0[IND(2*I,2*I-1,2*j+1)]+2.0*r0[IND(2*I,2*I,2*j+1)]+
						3.0*r0[IND(2*I,2*I-2,2*j)]+9.0*r0[IND(2*I,2*I-1,2*j)]+6.0*r0[IND(2*I,2*I,2*j)]+
						3.0*r0[IND(2*I,2*I-2,2*j-1)]+9.0*r0[IND(2*I,2*I-1,2*j-1)]+6.0*r0[IND(2*I,2*I,2*j-1)]+
						r0[IND(2*I,2*I-2,2*j-2)]+3.0*r0[IND(2*I,2*I-1,2*j-2)]+2.0*r0[IND(2*I,2*I,2*j-2)]);
	}
	for(int i=2;i<I;i++)	//bottom and top boundary
	{
		b1[IND(I,i,1)]=0.0625*(r0[IND(2*I,2*i-2,3)]+3.0*r0[IND(2*I,2*i-1,3)]+3.0*r0[IND(2*I,2*i,3)]+r0[IND(2*I,2*i+1,3)]+
						3.0*r0[IND(2*I,2*i-2,2)]+9.0*r0[IND(2*I,2*i-1,2)]+9.0*r0[IND(2*I,2*i,2)]+3.0*r0[IND(2*I,2*i+1,2)]+
						2.0*r0[IND(2*I,2*i-2,1)]+6.0*r0[IND(2*I,2*i-1,1)]+6.0*r0[IND(2*I,2*i,1)]+2.0*r0[IND(2*I,2*i+1,1)]);
		b1[IND(I,i,J)]=0.0625*(2.0*r0[IND(2*I,2*i-2,2*J)]+6.0*r0[IND(2*I,2*i-1,2*J)]+6.0*r0[IND(2*I,2*i,2*J)]+2.0*r0[IND(2*I,2*i+1,2*J)]+
						3.0*r0[IND(2*I,2*i-2,2*J-1)]+9.0*r0[IND(2*I,2*i-1,2*J-1)]+9.0*r0[IND(2*I,2*i,2*J-1)]+3.0*r0[IND(2*I,2*i+1,2*J-1)]+
						r0[IND(2*I,2*i-2,2*J-2)]+3.0*r0[IND(2*I,2*i-1,2*J-2)]+3.0*r0[IND(2*I,2*i,2*J-2)]+r0[IND(2*I,2*i+1,2*J-2)]);
	}
	b1[IND(I,1,1)]=0.0625*(2.0*r0[IND(2*I,1,3)]+3.0*r0[IND(2*I,2,3)]+r0[IND(2*I,3,3)]+
					6.0*r0[IND(2*I,1,2)]+9.0*r0[IND(2*I,2,2)]+3.0*r0[IND(2*I,3,2)]+
					4.0*r0[IND(2*I,1,1)]+6.0*r0[IND(2*I,2,1)]+2.0*r0[IND(2*I,3,1)]);	//bottom left corner cell
	b1[IND(I,I,1)]=0.0625*(r0[IND(2*I,2*I-2,3)]+3.0*r0[IND(2*I,2*I-1,3)]+2.0*r0[IND(2*I,2*I,3)]+
					3.0*r0[IND(2*I,2*I-2,2)]+9.0*r0[IND(2*I,2*I-1,2)]+6.0*r0[IND(2*I,2*I,2)]+
					2.0*r0[IND(2*I,2*I-2,1)]+6.0*r0[IND(2*I,2*I-1,1)]+4.0*r0[IND(2*I,2*I,1)]);	//bottom right corner cell
	b1[IND(I,1,J)]=0.0625*(4.0*r0[IND(2*I,1,2*J)]+6.0*r0[IND(2*I,2,2*J)]+2.0*r0[IND(2*I,3,2*J)]+
					6.0*r0[IND(2*I,1,2*J-1)]+9.0*r0[IND(2*I,2,2*J-1)]+3.0*r0[IND(2*I,3,2*J-1)]+
					2.0*r0[IND(2*I,1,2*J-2)]+3.0*r0[IND(2*I,2,2*J-2)]+r0[IND(2*I,3,2*J-2)]);	//top left corner cell
	b1[IND(I,I,J)]=0.0625*(2.0*r0[IND(2*I,2*I-2,2*J)]+6.0*r0[IND(2*I,2*I-1,2*J)]+4.0*r0[IND(2*I,2*I,2*J)]+
					3.0*r0[IND(2*I,2*I-2,2*J-1)]+9.0*r0[IND(2*I,2*I-1,2*J-1)]+6.0*r0[IND(2*I,2*I,2*J-1)]+
					r0[IND(2*I,2*I-2,2*J-2)]+3.0*r0[IND(2*I,2*I-1,2*J-2)]+2.0*r0[IND(2*I,2*I,2*J-2)]);	//top right corner cell
}
void GMG2::prolong(const int I,const int J,double *X1,double *X0)
{
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			X0[IND(2*I,2*i,2*j)]+=X1[IND(I,i,j)];
			X0[IND(2*I,2*i-1,2*j)]+=X1[IND(I,i,j)];
			X0[IND(2*I,2*i,2*j-1)]+=X1[IND(I,i,j)];
			X0[IND(2*I,2*i-1,2*j-1)]+=X1[IND(I,i,j)];
		}
	}
}
void GMG2::reini()
{
	int N=I*J/4;	//level 1
	for(int i=0;i<N;i++)
		X1[i]=b1[i]=r1[i]=0.0;
	N/=4;	//level 2
	for(int i=0;i<N;i++)
		X2[i]=b2[i]=r2[i]=0.0;
	N/=4;	//level 3
	for(int i=0;i<N;i++)
		X3[i]=b3[i]=r3[i]=0.0;
	N/=4;	//level 4
	for(int i=0;i<N;i++)
		X4[i]=b4[i]=r4[i]=0.0;
}
