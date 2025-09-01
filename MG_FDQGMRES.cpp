//MULTIGRID PRECONDITIONED FDQGMRES(K)
class MG_FDQGMRES:public GMG2
{
	double **V,**p;	//truncated Arnoldi basis, and search vector
	double *w,*c,*s,*H;	//rotation vectors and last column of Hessenberg matrix
	double *Z;	//solution of preconditioner
	double var1,var2,gama0,gama1;	//variables required
	double *temp;	//pointer required for shifting of vectors
	int K,N;	//K = dimension of Krylov subspace and N = number of elements in the solution vector
	int MAX;	//max no of iterations
	void SSOR(int CNT0,int M=1);	//SSOR(M) preconditioner
	void MG_SSOR(int CNT0,int M=1);	//multigrid SSOR(M) preconditioner
	protected:
			void solve(int *C,int *R,double *A,double *X,double *b);	//FDQGMRES algorithm
	public:
			MG_FDQGMRES(int k,int n,int max); ~MG_FDQGMRES();
};
MG_FDQGMRES::MG_FDQGMRES(int k,int n,int max)
{
	K=k; N=n; MAX=max;
	V=new double*[K];
	p=new double*[K];
	w=new double[N];
	Z=new double[N];
	c=new double[K];
	s=new double[K];
	H=new double[K+1];
	for(int i=0;i<K;i++)
	{
		V[i]=new double[N];
		p[i]=new double[N];
	}
	cout<<"MG_FDQGMRES: MEMORY ALLOCATED"<<endl;
}
MG_FDQGMRES::~MG_FDQGMRES()
{
	for(int i=0;i<K;i++)
	{
		delete[] V[i];
		delete[] p[i];
	}
	delete[] V;
	delete[] p;
	delete[] w;
	delete[] Z;
	delete[] c;
	delete[] s;
	delete[] H;
	cout<<"MG_FDQGMRES: MEMORY RELEASED"<<endl;
}
void MG_FDQGMRES::solve(int *C,int *R,double *A,double *X,double *b)
{
	int cnt=0;
	gama0=0.0;
	for(int i=0;i<N;i++)	//matrix vector multiplication
	{
		var1=0.0;
		for(int d=R[i];d<R[i+1];d++) var1+=A[d]*X[C[d]];
		V[0][i]=b[i]-var1;	//initial residual vector
		gama0+=V[0][i]*V[0][i];
	}
	gama0=sqrt(gama0);	//initial residual norm
	if(abs(gama0)<=TOL)
	{
		//cout<<"FDQGMRES: BASIS FOR KRYLOV SUBSPACE CANNOT BE OBTAINED!"<<endl;
		return;
	}
	for(int d=0;d<N;d++) V[0][d]/=gama0;	//normalization of initial residual
	for(int M=0,CNT0,CNT1,CNT2;M<MAX;M++)	//convergence loop
	{
		CNT0=(M<=K-1) ? M : K-1;	//special iteration counters
		CNT1=(M+1<=K-1) ? M+1 : K-1;
		CNT2=(M<=K) ? M : K;
		//SSOR(CNT0,8);	//compute solution of preconditioner
		MG_SSOR(CNT0);
		for(int i=0;i<N;i++)	//matrix vector multiplication
		{
			w[i]=0.0;
			for(int d=R[i];d<R[i+1];d++) w[i]+=A[d]*Z[C[d]];
		}
		for(int i=0;i<=CNT0;i++)
		{
			H[i]=0.0;	//reinitialization
			for(int d=0;d<N;d++) H[i]+=w[d]*V[i][d];
			for(int d=0;d<N;d++) w[d]-=H[i]*V[i][d];
		}
		H[CNT0+1]=0.0;	//reinitialization
		for(int d=0;d<N;d++) H[CNT0+1]+=pow(w[d],2.0);
		H[CNT0+1]=sqrt(H[CNT0+1]);
		if(M>=K-1)	//shift V matrix
		{
			temp=V[0];
			for(int i=0;i<K-1;i++) V[i]=V[i+1];
			V[K-1]=temp;
		}
		if(abs(H[CNT0+1])<=TOL) { cout<<"FDQGMRES: LUCKY BREAKDOWN!"<<endl; }
		else
		{
			for(int d=0;d<N;d++) { V[CNT1][d]=w[d]/H[CNT0+1]; w[d]=0.0; }
		}
		if(M<K)	//rotations for full orthogonalization
		{
			for(int i=0;i<M;i++)	//application of previous rotation matrices
			{
				var1=H[i]; var2=H[i+1];
				H[i]=c[i]*var1+s[i]*var2;
				H[i+1]=-s[i]*var1+c[i]*var2;
			}
			if(abs(H[M+1])>abs(H[M])) { var1=H[M]/H[M+1]; s[M]=pow((1.0+var1*var1),-0.5); c[M]=s[M]*var1; }	//rotation vector calculation
			else { var1=H[M+1]/H[M]; c[M]=pow((1.0+var1*var1),-0.5); s[M]=c[M]*var1; }
			H[M]=c[M]*H[M]+s[M]*H[M+1];
		}
		else	//rotations for truncated orthogonalization
		{
			var1=0.0;
			for(int i=0;i<K;i++)
			{
				var2=H[i];
				H[i]=c[i]*var1+s[i]*var2;
				var1=-s[i]*var1+c[i]*var2;
			}
			for(int i=0;i<K-1;i++) { c[i]=c[i+1]; s[i]=s[i+1]; }	//shift rotation vectors
			if(abs(H[K])>abs(var1)) { var2=var1/H[K]; s[K-1]=pow((1.0+var2*var2),-0.5); c[K-1]=s[K-1]*var2; }	//rotation vector calculation
			else { var2=H[K]/var1; c[K-1]=pow((1.0+var2*var2),-0.5); s[K-1]=c[K-1]*var2; }
			H[K]=c[K-1]*var1+s[K-1]*H[K];
		}
		gama1=-s[CNT0]*gama0; gama0*=c[CNT0];	//application of final rotation matrix
		for(int j=0;j<CNT2;j++)	//calculation of new search vector and update of solution vector
			for(int d=0;d<N;d++) w[d]+=H[j]*p[j][d];
		if(M>=K)	//shift p matrix
		{
			temp=p[0];
			for(int i=0;i<K-1;i++) p[i]=p[i+1];
			p[K-1]=temp;
		}
		for(int d=0;d<N;d++)
		{
			p[CNT0][d]=(Z[d]-w[d])/H[CNT2];
			X[d]+=gama0*p[CNT0][d];
			Z[d]=0.0;	//reinitialization
		}
		cnt++;
		if((abs(gama1)<=TOL)||(abs(H[CNT0+1])<=TOL)) break;	//solution converged
		else gama0=gama1;	//update gamma
		//cout<<"MG_FDQGMRES: "<<cnt<<" "<<abs(gama1)<<endl;
	}
	//cout<<"MG_FDQGMRES: "<<cnt<<" "<<abs(gama1)<<endl;
}
void MG_FDQGMRES::SSOR(int CNT0,int M)
{
	double diag,sum;	//diagonal element of the matrix
	for(int m=0;m<M;m++)
	{
		for(int i=0;i<N;i++)	//matrix vector multiplication, forward sweep
		{
			sum=0.0;	//reinitialization
			for(int d=R0[i];d<R0[i+1];d++)
			{
				if(C0[d]==i) diag=A0[d];
				else sum+=A0[d]*Z[C0[d]];
			}
			Z[i]=(V[CNT0][i]-sum)/diag;
		}
		for(int i=N-1;i>=0;i--)	//matrix vector multiplication, backward sweep
		{
			sum=0.0;	//reinitialization
			for(int d=R0[i];d<R0[i+1];d++)
			{
				if(C0[d]==i) diag=A0[d];
				else sum+=A0[d]*Z[C0[d]];
			}
			Z[i]=(V[CNT0][i]-sum)/diag;
		}
		for(int i=0;i<N;i++)	//residual vector
		{
			sum=0.0;	//reinitialization
			for(int d=R0[i];d<R0[i+1];d++) sum+=A0[d]*Z[C0[d]];
			r0[i]=V[CNT0][i]-sum;
		}
	}
}
void MG_FDQGMRES::MG_SSOR(int CNT0,int M)
{
	double res_0,res_1,res_2,res_3,res_4;	//residuals at each level
	for(int cnt=1;cnt<=M;cnt++)	//5-grid V-cycle
	{
		SSOR(CNT0,1);
		restrict(I/2,J/2,r0,b1);
		res_1=solve_SSOR(2,I*J/4,C1,R1,A1,X1,b1,r1);
		restrict(I/4,J/4,r1,b2);
		res_2=solve_SSOR(8,I*J/16,C2,R2,A2,X2,b2,r2);
		restrict(I/8,J/8,r2,b3);
		res_3=solve_SSOR(32,I*J/64,C3,R3,A3,X3,b3,r3);
		restrict(I/16,J/16,r3,b4);
		res_4=solve_SSOR(128,I*J/256,C4,R4,A4,X4,b4,r4);

		prolong(I/16,J/16,X4,X3);
		res_3=solve_SSOR(1,I*J/64,C3,R3,A3,X3,b3,r3);
		prolong(I/8,J/8,X3,X2);
		res_2=solve_SSOR(1,I*J/16,C2,R2,A2,X2,b2,r2);
		prolong(I/4,J/4,X2,X1);
		res_1=solve_SSOR(1,I*J/4,C1,R1,A1,X1,b1,r1);
		prolong(I/2,J/2,X1,Z);
		SSOR(CNT0,1);
		reini();
		//cout<<res_1<<"\t"<<res_2<<"\t"<<res_3<<"\t"<<res_4<<endl;
	}
}
