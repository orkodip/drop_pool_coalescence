//MULTIGRID PRECONDITIONED BiCG-stab
class MG_BICGSTAB:public GMG2
{
	double *r0,*r,*p,*s,*m,*n;	//vectors required
	double *y,*t;	//solution of preconditioner
	double alpha,beta,gama,omega,sum,norm;	//variables required
	int N,MAX;	//number of elements in the vectors, and max no of iterations
	void SSOR(double *xp,double *bp,int M=1);	//SSOR(M) preconditioner
	void MG_SSOR(double *xp,double *bp,int M=1);	//multigrid SSOR(M) preconditioner
	protected:
			void solve(int *C,int *R,double *A,double *X,double *b);	//BiCG-stab algorithm
	public:
			MG_BICGSTAB(int n,int max); ~MG_BICGSTAB();
};
MG_BICGSTAB::MG_BICGSTAB(int n1,int max)
{
	N=n1; MAX=max;
	r0=new double[N];
	r=new double[N];
	p=new double[N];
	s=new double[N];
	m=new double[N];
	n=new double[N];
	y=new double[N];
	t=new double[N];
	cout<<"MG_BICGSTAB: MEMORY ALLOCATED"<<endl;
}
MG_BICGSTAB::~MG_BICGSTAB()
{
	delete[] r0;
	delete[] r;
	delete[] p;
	delete[] s;
	delete[] m;
	delete[] n;
	delete[] y;
	delete[] t;
	cout<<"MG_BICGSTAB: MEMORY RELEASED"<<endl;
}
void MG_BICGSTAB::solve(int *C,int *R,double *A,double *X,double *b)
{
	int M;	//iteration counter
	norm=0.0;	//reinitialization
	for(int i=0;i<N;i++)	//matrix vector multiplication
	{
		sum=0.0;
		for(int d=R[i];d<R[i+1];d++) sum+=A[d]*X[C[d]];
		r[i]=b[i]-sum;	//initial residual vector
		r0[i]=r[i]; p[i]=r[i];
		norm+=r[i]*r[i];
	}
	gama=norm;
	norm=sqrt(norm);	//initial residual norm
	if(abs(norm)<=TOL) return;
	for(M=0;M<MAX;M++)	//convergence loop
	{
		//SSOR(y,p,8);	//compute solution of 1st preconditioner
		MG_SSOR(y,p);
		norm=0.0;	//reinitialization
		for(int i=0;i<N;i++)	//matrix vector multiplication
		{
			m[i]=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++) m[i]+=A[d]*y[C[d]];
			norm+=m[i]*r0[i];
		}
		alpha=gama/norm;
		norm=0.0;	//reinitialization
		for(int d=0;d<N;d++)
		{
			s[d]=r[d]-alpha*m[d]; norm+=s[d]*s[d];
			t[d]=0.0;	//reinitialization
		}
		norm=sqrt(norm);
		if(norm<=TOL)	//convergence check1
		{
			for(int d=0;d<N;d++) X[d]+=alpha*y[d];
			break;
		}
		//SSOR(t,s,8);	//compute solution of 2nd preconditioner
		MG_SSOR(t,s);
		sum=norm=0.0;	//reinitialization
		for(int i=0;i<N;i++)	//matrix vector multiplication
		{
			n[i]=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++) n[i]+=A[d]*t[C[d]];
			sum+=n[i]*s[i]; norm+=n[i]*n[i];
		}
		omega=sum/norm;
		norm=0.0;	//reinitialization
		for(int d=0;d<N;d++)
		{
			X[d]+=alpha*y[d]+omega*t[d];
			r[d]=s[d]-omega*n[d];
			norm+=r[d]*r[d];
		}
		norm=sqrt(norm);
		//cout<<"MG_BICGSTAB: "<<M+1<<" "<<norm<<endl;
		if(norm<=TOL) break;	//convergence check2
		sum=gama;	//sum is used to store gama_M
		gama=0.0;	//reinitialization
		for(int d=0;d<N;d++) gama+=r[d]*r0[d];	//calculation of gama_M+1
		if(abs(gama)<=TOL)	//restart condition
		{
			gama=0.0;	//reinitialization
			for(int d=0;d<N;d++)
			{
				r0[d]=r[d];
				p[d]=r[d];
				gama+=r[d]*r[d];
			}
		}
		else
		{
			beta=alpha*gama/(omega*sum);
			for(int d=0;d<N;d++)
			{
				p[d]=r[d]+beta*(p[d]-omega*m[d]);
				y[d]=0.0;	//reinitialization
			}
		}
	}
	//cout<<"MG_BICGSTAB: "<<M+1<<endl;
}
void MG_BICGSTAB::SSOR(double *xp,double *bp,int M)
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
				else sum+=A0[d]*xp[C0[d]];
			}
			xp[i]=(bp[i]-sum)/diag;
		}
		for(int i=N-1;i>=0;i--)	//matrix vector multiplication, backward sweep
		{
			sum=0.0;	//reinitialization
			for(int d=R0[i];d<R0[i+1];d++)
			{
				if(C0[d]==i) diag=A0[d];
				else sum+=A0[d]*xp[C0[d]];
			}
			xp[i]=(bp[i]-sum)/diag;
		}
		for(int i=0;i<N;i++)	//residual vector
		{
			sum=0.0;	//reinitialization
			for(int d=R0[i];d<R0[i+1];d++) sum+=A0[d]*xp[C0[d]];
			GMG2::r0[i]=bp[i]-sum;
		}
	}
}
void MG_BICGSTAB::MG_SSOR(double *xp,double *bp,int M)
{
	double res_0,res_1,res_2,res_3,res_4;	//residuals at each level
	for(int cnt=1;cnt<=M;cnt++)	//5-grid V-cycle
	{
		SSOR(xp,bp,1);
		restrict(I/2,J/2,GMG2::r0,b1);
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
		prolong(I/2,J/2,X1,xp);
		SSOR(xp,bp,1);
		reini();
		//cout<<res_1<<"\t"<<res_2<<"\t"<<res_3<<"\t"<<res_4<<endl;
	}
}
