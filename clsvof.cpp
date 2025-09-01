/*CONSISTENT TRANSPORT AXISYMMETRIC CLSVOF ALGORITHM USING PLIC SCHEME FOR VOLUME FRACTIONS AND ENO2 SCHEME FOR LEVEL SETS.
BOUNDARY CONDITIONS ARE AS FOLLOWS.
LEFT BOUNDARY - AXISYMMETRIC (FREE SLIP AND NO PENETRATION)
RIGHT BOUNDARY - FREE SLIP AND NO PENETRATION
BOTTOM BOUNDARY - NO SLIP AND NO PENETRATION
TOP BOUNDARY - OUTFLOW CONDITION
*/
class CLSVOF:public MBASE
{
	int **tag;	//tags of interfacial cells
	double **Ft,**Phit,**rhot;	//intermediate volume fractions, level set, and density field
	double **F_trc,**Phi_trc;	//passive tracers of the drop
	double **A_xt,**A_yt;	//intermediate advection terms for R and Z momentum equations
	double mass_act;	//actual mass

	void updt_ghost(double **Phi);	//update the ghost cells of cell centered field
	double ENO2(int flag,int i,int j,double **Phi,double V);	//ENO2 scheme
	double UP1(int flag,int i,int j,double **Phi,double V);	//1st order upwind scheme
	double recon_int(int i,int j,double Fa,double **Phia,VECTOR *N);	//reconstruct the interface
	double vol_frac_flux(int flag,int i,int j,double Fa,double **Phia,double V);	//calculate the volume fraction advection flux
	void reinit(double **Fa,double **Phia);	//LS reinitialization algorithm
	void adv_X();	//solve the advection equation (x direction followed by y direction)
	void adv_Y();	//solve the advection equation (y direction followed by x direction)
	void adv_trc();	//solve the advection equation for tracer fields (y direction followed by x direction)
	void dis_den();	//discard density field after advection
	void LS_bsd_trunc();	//LS based truncation of volume fraction field
	public:
			CLSVOF(); ~CLSVOF();
			void ini(double xc,double yc,double a,double b,double h);
			void den_ini();	//initialize density field at n+1 time step
			void prop_updt();	//update density field for next time step
			void solve(int n);	//CLSVOF advection algorithm
			void mass_err();	//calculate mass error
			void write_advt(int t);	//tecplot file output
			void lsvf_write(int t);	//tecplot file output
			void trc_write(int t);	//tecplot file output for the tracer fields
			void ls_complete(int t);	//ls file output including the ghost cells
};
CLSVOF::CLSVOF()
{
	tag=new int*[J+2];
	Phit=new double*[J+2];
	Ft=new double*[J+2];
	Phi_trc=new double*[J+2];
	F_trc=new double*[J+2];
	rhot=new double*[J+1];
	A_xt=new double*[J+1];
	A_yt=new double*[J+1];
	for(int i=0;i<J+2;i++)
	{
		tag[i]=new int[I+2];
		Phit[i]=new double[I+2];
		Ft[i]=new double[I+2];
		Phi_trc[i]=new double[I+2];
		F_trc[i]=new double[I+2];
		if(i<J+1)
		{
			rhot[i]=new double[I+1];
			A_xt[i]=new double[I+1];
			A_yt[i]=new double[I+1];
		}
	}
	cout<<"CLSVOF: MEMORY ALLOCATED"<<endl;
}
CLSVOF::~CLSVOF()
{
	for(int i=0;i<J+2;i++)
	{
		delete[] tag[i];
		delete[] Phit[i];
		delete[] Ft[i];
		delete[] Phi_trc[i];
		delete[] F_trc[i];
		if(i<J+1)
		{
			delete[] rhot[i];
			delete[] A_xt[i];
			delete[] A_yt[i];
		}
	}
	delete[] tag;
	delete[] Phit;
	delete[] Ft;
	delete[] Phi_trc;
	delete[] F_trc;
	delete[] rhot;
	delete[] A_xt;
	delete[] A_yt;
	cout<<"CLSVOF: MEMORY RELEASED"<<endl;
}
void CLSVOF::ini(double xc,double yc,double a,double b,double h)
{
	INI ms(Xm,Ym,CX,CY,F,Phi,xc,yc,a,b,h);
	ms.VF(); ms.LS();	//initial exact volume fraction and level set field
	INI ms_trc(Xm,Ym,CX,CY,F_trc,Phi_trc,xc,yc,a,b,0.0);
	ms_trc.VF(); ms_trc.LS();	//volume fraction and level set field of the tracer fields
	mass_act=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass_act+=F[j][i]*CX[i]*dx*dy;
	updt_ghost(F);
	for(int j=1;j<=J;j++)	//initialize density, viscosity and advection field
	{
		for(int i=1;i<=I;i++)
		{
			rho_np1[j][i]=rho_n[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			A_x[j][i]=u[j][i];
			A_y[j][i]=v[j][i];
		}
	}
	updt_ghost(rho_n); updt_ghost(rho_np1); updt_ghost(mu);
}
void CLSVOF::den_ini()
{
	for(int j=0;j<=J+1;j++)	//calculate density field at n+1 step (including ghost nodes)
		for(int i=0;i<=I+1;i++)
			rho_np1[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
}
void CLSVOF::dis_den()
{
	for(int j=0;j<=J+1;j++)	//calculate density field at n+1 step (including ghost nodes)
		for(int i=0;i<=I+1;i++)
			rho_np1[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
}
void CLSVOF::prop_updt()
{
	for(int j=0;j<=J+1;j++)	//including ghost nodes
	{
		for(int i=0;i<=I+1;i++)
		{
			rho_n[j][i]=rho_np1[j][i];
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			if((i>=1)&&(i<=I)&&(j>=1)&&(j<=J))	//only inner domain
			{
				A_x[j][i]=u[j][i];
				A_y[j][i]=v[j][i];
			}
		}
	}
}
void CLSVOF::updt_ghost(double **Phia)
{
	for(int j=1;j<=J;j++)	//left and right ghost nodes (Neumann bc)
	{
		Phia[j][0]=Phia[j][1];
		Phia[j][I+1]=Phia[j][I];
	}
	for(int i=1;i<=I;i++)	//bottom and top ghost nodes (Neumann bc)
	{
		Phia[0][i]=Phia[1][i];
		Phia[J+1][i]=Phia[J][i];
	}
}
double CLSVOF::ENO2(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	double plus,minus;	//plus and minus flux
	if(flag==0)	//X direction flux
	{
		plus=Phia[j][i+1]-0.5*MINMOD((Phia[j][i+1]-Phia[j][i]),(Phia[j][i+2]-Phia[j][i+1]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j][i-1]),(Phia[j][i+1]-Phia[j][i]));
	}
	else if(flag==1)	//Y direction flux
	{
		plus=Phia[j+1][i]-0.5*MINMOD((Phia[j+1][i]-Phia[j][i]),(Phia[j+2][i]-Phia[j+1][i]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j-1][i]),(Phia[j+1][i]-Phia[j][i]));
	}
	if(V>0.0) return minus;
	else return plus;
}
double CLSVOF::UP1(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	if(flag==0)	//X direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j][i+1];
	}
	else if(flag==1)	//Y direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j+1][i];
	}
	else { cout<<"CLSVOF: ERROR IN UPWIND SCHEME!"<<endl; return 0; }
}
void CLSVOF::reinit(double **Fa,double **Phia)
{
	double temp,a,b,h=dx;
	VECTOR N;	//interface normal vector
	double gamma=0.5*h;	//gamma is the distance parameter
//--------------------INITIALIZATION SCHEME (GEOMETRIC+ALGEBRAIC)--------------------------------------------------
	for(int j=1;j<=J;j++)	//initialize the LS and tag values and calculate LS for interfacial cells
	{
		for(int i=1;i<=I;i++)
		{
			if((abs(Fa[j][i])<=1e-2)||(abs(Fa[j][i]-1.0)<=1e-2))	//for single phase cells
			{
				tag[j][i]=0;
				Phit[j][i]=(2.0*Fa[j][i]-1.0)*gamma*1e6;
				continue;
			}
			tag[j][i]=1;
			a=recon_int(i,j,Fa[j][i],Phia,&N);
			b=0.5*(abs(N.x)*dx+abs(N.y)*dy);
			Phit[j][i]=a-b;
		}
	}
	for(int j=1;j<=J;j++)	//discard old LS field
		for(int i=1;i<=I;i++)
			Phia[j][i]=Phit[j][i];
	for(int j=1;j<=J;j++)	//checking in horizontal direction
	{
		for(int i=1;i<I;i++)	//excluding right boundary cells
		{
			if((SGN(Phia[j][i]*Phia[j][i+1])!=1)&&((tag[j][i]+tag[j][i+1])==0))	//LS changes sign in untagged cells
			{
				Phia[j][i]=(2.0*Fa[j][i]-1.0)*gamma; tag[j][i]=1;
			}
		}
	}
	for(int i=1;i<=I;i++)	//checking in vertical direction
	{
		for(int j=1;j<J;j++)	//excluding top boundary cells
		{
			if((SGN(Phia[j][i]*Phia[j+1][i])!=1)&&((tag[j][i]+tag[j+1][i])==0))	//LS changes sign in untagged cells
			{
				Phia[j][i]=(2.0*Fa[j][i]-1.0)*gamma; tag[j][i]=1;
			}
		}
	}
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(left and right boundaries)
	{
		Phia[j][0]=Phia[j][1]; tag[j][0]=tag[j][1];
		Phia[j][I+1]=Phia[j][I]; tag[j][I+1]=tag[j][I];
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(bottom and top boundaries)
	{
		Phia[0][i]=Phia[1][i]; tag[0][i]=tag[1][i];
		Phia[J+1][i]=Phia[J][i]; tag[J+1][i]=tag[J][i];
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
				if(i==0) a=Phia[j][i+1];	//left boundary
				else if(i==(I+1)) a=Phia[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phia[j][i])==1.0) a=MIN2(Phia[j][i+1],Phia[j][i-1]);
					else a=MAX2(Phia[j][i+1],Phia[j][i-1]);
				}
				if(j==0) b=Phia[j+1][i];	//bottom boundary
				else if(j==(J+1)) b=Phia[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phia[j][i])==1.0) b=MIN2(Phia[j+1][i],Phia[j-1][i]);
					else b=MAX2(Phia[j+1][i],Phia[j-1][i]);
				}
				if(SGN(Phia[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-EPS) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phia[j][i]=MIN2(Phia[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-EPS) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phia[j][i]=MAX2(Phia[j][i],temp);
				}
			}
		}
	}
}
double CLSVOF::recon_int(int i,int j,double Fa,double **Phia,VECTOR *N)
{
	//---------------------INTERFACE NORMAL CALCULATION----------------
	(*N).x=0.5*(Phia[j][i+1]-Phia[j][i-1])/dx;
	(*N).y=0.5*(Phia[j+1][i]-Phia[j-1][i])/dy;
	if((abs((*N).x)<=EPS)&&(abs((*N).y)<=EPS))
	{
		(*N).x=0.707106781; (*N).y=-0.707106781;	//FAIL-SAFE VALUE
		cout<<"CLSVOF: INTERFACE NORMAL NOT FOUND FROM LS FIELD, FAIL-SAFE VALUE USED!"<<endl;
		cout<<"i = "<<i<<", j = "<<j<<", Phia[j][i] = "<<Phia[j][i]<<endl;
		cout<<"Phia[j][i+1] = "<<Phia[j][i+1]<<", Phia[j][i-1] = "<<Phia[j][i-1]<<endl;
		cout<<"Phia[j+1][i] = "<<Phia[j+1][i]<<", Phia[j-1][i] = "<<Phia[j-1][i]<<endl;
		sleep(3);
	}
	*N=(*N).unit();
	//------------------------------------------------------------------

	if(abs(Fa)<=EPS) return 0.0;	//empty cell
	else if(abs(1.0-Fa)<=EPS) return (abs((*N).x)*dx+abs((*N).y)*dy);	//completely filled cell

	double sx,sy,rx,F_min,F_max,B,C;
	int f;
	sx=abs((*N).x)*dx; sy=abs((*N).y)*dy; rx=abs((*N).x)*CX[i]; f=-SGN((*N).x);
	if(sy<sx)	//case A
	{
		if(abs((*N).y)<=EPS) { F_min=0.0; F_max=1.0; }
		else
		{
			F_min=sy/sx*(0.5+0.08333333*f*(2.0*sy-3.0*sx)/rx);	//Fy
			F_max=1.0+F_min-sy/sx;	//Fx
		}
		if(abs(Fa-F_min)<=EPS) return sy;
		else if(abs(Fa-F_max)<=EPS) return sx;
		else if((Fa>F_min)&&(Fa<F_max))	//solve quadratic equation
		{
			B=f*rx-0.5*(sx+sy);
			C=2.0*(Fa-F_min)*f*rx*sx+pow((f*rx-0.5*(sx-sy)),2.0);
			return (-B+f*sqrt(C));
		}
	}
	else if(sx<sy)	//case D
	{
		if(abs((*N).x)<=EPS) { F_min=0.0; F_max=1.0; }
		else
		{
			F_min=sx/sy*(0.5-0.08333333*f*sx/rx);	//Fx
			F_max=1.0+F_min-sx/sy;	//Fy
		}
		if(abs(Fa-F_min)<=EPS) return sx;
		else if(abs(Fa-F_max)<=EPS) return sy;
		else if((Fa>F_min)&&(Fa<F_max))	return ((Fa-F_min)*sy+sx);	//solve linear equation
	}

	//----------------SOLUTION OF CUBIC EQUATION (NRBS METHOD)--------------------

	double func,dfunc,dx_new,dx_old,xh=MIN2(sx,sy),fh,xl=0.0,fl,sc;	//variables are reused:C=Fc,B=temp

	if(Fa<F_min) { C=Fa; }
	else if(Fa>F_max) { C=1.0-Fa; f=-f; }

	fl=-2.0*C*rx*sx*sy;
	fh=0.33333333*f*pow(xh,3.0)+(rx-0.5*f*sx)*pow(xh,2.0)+fl;

	sc=xh*C/F_min;	//initial guess
	if((rx-0.5*f*sx)<(0.33333333*f*sc)) sc=MIN2(xh,(pow((-fl/(0.33333333*f)),0.33333333)));
	else sc=MIN2(xh,(sqrt(-fl/(rx-0.5*f*sx+0.33333333*f*sc))));

	if(fl>0.0) { B=xl; xl=xh; xh=B; }	//orient the search so that f(xl)<0
	dx_new=dx_old=abs(xh-xl);	//intial bracket sizes

	func=0.33333333*f*pow(sc,3.0)+(rx-0.5*f*sx)*pow(sc,2.0)+fl;	//initial function value
	dfunc=f*pow(sc,2.0)+2.0*(rx-0.5*f*sx)*sc;	//initial derivative value

	for(int t=1;t<=20;t++)
	{
		if((((sc-xh)*dfunc-func)*((sc-xl)*dfunc-func)>=0.0)||(abs(2.0*func)>abs(dx_old*dfunc)))
		{	//use BS if NR is out of bounds or interval is not decreasing fast enough
			dx_old=dx_new;
			dx_new=0.5*(xh-xl);
			sc=xl+dx_new;
			if(abs(sc-xl)<=cc_NRBS) break;	//negligible change in root
		}
		else	//use NR
		{
			dx_old=dx_new;
			dx_new=func/dfunc;
			B=sc;
			sc-=dx_new;
			if(abs(B-sc)<=cc_NRBS) break;	//negligible change in root
		}
		if(abs(dx_new)<=cc_NRBS) break;	//convergence criteria

		func=0.33333333*f*pow(sc,3.0)+(rx-0.5*f*sx)*pow(sc,2.0)+fl;	//evaluate function
		dfunc=f*pow(sc,2.0)+2.0*(rx-0.5*f*sx)*sc;	//evaluate derivative

		if(abs(func)<=cc_NRBS) break;	//convergence criteria
		else if(func<0.0) xl=sc;	//reduce the bracket size
		else if(func>0.0) xh=sc;
		else
		{
			VECTOR m;
			cout<<"CLSVOF: ERROR IN NRBS ITERATIONS, func = "<<func<<endl;
			cout<<"i = "<<i<<", j = "<<j<<", Phia[j][i] = "<<Phia[j][i]<<endl;
			cout<<"Phia[j][i+1] = "<<Phia[j][i+1]<<", Phia[j][i-1] = "<<Phia[j][i-1]<<endl;
			cout<<"Phia[j+1][i] = "<<Phia[j+1][i]<<", Phia[j-1][i] = "<<Phia[j-1][i]<<endl;
			cout<<"Nx = "<<(*N).x<<", Ny = "<<(*N).y<<endl;
			m.x=0.5*(Phia[j][i+1]-Phia[j][i-1])/dx;
			m.y=0.5*(Phia[j+1][i]-Phia[j-1][i])/dy;
			cout<<"mx = "<<m.x<<", my = "<<m.y<<endl;
			sleep(3);
		}
		if(t==20) cout<<"CLSVOF: NOT CONVERGED WITHIN ITERATION LIMIT"<<endl;
	}
	if(Fa<F_min) return sc;
	else if(Fa>F_max) return (sx+sy-sc);
	else { cout<<"CLSVOF: ERROR IN INTERFACE RECONSTRUCTION"<<endl; return 0; }
}
double CLSVOF::vol_frac_flux(int flag,int i,int j,double Fa,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;	//zero velocity
	else if((abs(Fa)<=EPS)||(Fa<0.0)) return 0.0;	//empty cell or over-empty cell
	else if((abs(Fa-1.0)<=EPS)||(Fa>1.0))	//completely filled cell or over-filled cell
	{
		if(flag==0) return ((CX[i]+0.5*SGN(V)*dx)*V);	//advection flux in r direction
		else return V;	//advection flux in z direction
	}
	VECTOR N;
	double s,r_adv,rx,Dx,Dy,Vol_adv,dsx,dsy,dsm,ds_max,ds_min,ds_c,F_c,dv0;
	int flag1=0,f;
	s=recon_int(i,j,Fa,Phia,&N);	//reconstruct interface
	if(flag==0)	//setup for advection in r direction
	{
		r_adv=CX[i]+0.5*SGN(V)*dx;	//r_a
		Dx=SGN(V)*(r_adv-sqrt(r_adv*r_adv-2.0*V*dt*r_adv));	//size of advection cell
		r_adv-=0.5*SGN(V)*Dx;	//r_adv
		if((N.x*V)<0.0)	//complementary part is considered
		{
			Dx=dx-Dx; flag1=1;
			r_adv-=0.5*SGN(V)*dx;	//complementary r_adv
		}
		Dy=dy;
	}
	else	//setup for advection in z direction
	{
		r_adv=CX[i];
		Dy=abs(V)*dt;	//size of advection cell
		if((N.y*V)<0.0) { Dy=dy-Dy; flag1=1; }	//complementary part is considered
		Dx=dx;
	}
	dsx=abs(N.x)*Dx; dsy=abs(N.y)*Dy; dsm=dsx+dsy; rx=r_adv*abs(N.x); f=-SGN(N.x);
	ds_max=MAX2(dsx,dsy); ds_min=MIN2(dsx,dsy);
	if((s>ds_min)&&(s<ds_max))	//calculate considered volume fraction
	{
		if(dsy<dsx) F_c=0.5*((rx-0.5*f*dsx)*(2.0*s-dsy)+f*(s*s-s*dsy+0.33333333*dsy*dsy))/(rx*dsx);
		else if(dsx<dsy) F_c=(s-0.5*dsx*(1.0+0.16666666*f*dsx/(rx+1e-16)))/dsy;
	}
	else
	{
		if(s<ds_min) ds_c=s;
		else if(s>ds_max) { ds_c=MAX2(0,dsm-s); f=-f; }
		F_c=0.5*((rx-0.5*f*dsx)*ds_c*ds_c+0.33333333*f*ds_c*ds_c*ds_c)/(rx*dsx*dsy+1e-16);
		if(s>ds_max) F_c=1.0-F_c;
	}
	dv0=r_adv*Dx*Dy*F_c;	//considered fluid volume
	if(flag1==1) Vol_adv=Fa*CX[i]*dx*dy-dv0;	//adjustment with complementary flag
	else Vol_adv=dv0;
	if(flag==0) return (SGN(V)*Vol_adv/(dt*dy));	//advection flux in r direction
	else return (SGN(V)*Vol_adv/(CX[i]*dx*dt));	//advection flux in z direction
}
void CLSVOF::adv_X()
{
	double F_adv[2],LS_adv[2],rho_adv[2],u_f[2],v_f[2];	//advection fluxes for volume fractions, LS, density, and momentum
	int up;	//donor cell index
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(0,1,j,Phi,u_EW[j][0]);	//upwind scheme for left boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no penetration
				v_f[0]=A_y[j][1];	//free slip
			}
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phi,u_EW[j][I]);	//upwind scheme for right boundary flux
				rho_adv[1]=0.0;
				u_f[1]=0.0;	//no penetration
				v_f[1]=A_y[j][I];	//free slip
			}
			else	//inner domain
			{
				up=(u_EW[j][i]<0.0)?(i+1):i;
				F_adv[1]=vol_frac_flux(0,up,j,F[j][up],Phi,u_EW[j][i]);
				LS_adv[1]=ENO2(0,i,j,Phi,u_EW[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*Xm[i]*u_EW[j][i];
				if((i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use ENO2 scheme
				{
					u_f[1]=ENO2(0,i,j,A_x,u_EW[j][i]);
					v_f[1]=ENO2(0,i,j,A_y,u_EW[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(0,i,j,A_x,u_EW[j][i]);
					v_f[1]=UP1(0,i,j,A_y,u_EW[j][i]);
				}
			}
			Ft[j][i]=(F[j][i]-dt/(CX[i]*dx)*(F_adv[1]-F_adv[0]))/(1.0-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]));
			rhot[j][i]=(rho_n[j][i]-dt/(CX[i]*dx)*(rho_adv[1]-rho_adv[0]))/(1.0-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]));
			Phit[j][i]=(Phi[j][i]-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]*LS_adv[1]-Xm[i-1]*u_EW[j][i-1]*LS_adv[0]))/(1.0-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]));
			A_xt[j][i]=(rho_n[j][i]*A_x[j][i]-dt/(CX[i]*dx)*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]))/(1.0-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]));
			A_yt[j][i]=(rho_n[j][i]*A_y[j][i]-dt/(CX[i]*dx)*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]))/(1.0-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]));
			F_adv[0]=F_adv[1];	//updation for the next cell
			rho_adv[0]=rho_adv[1];
			LS_adv[0]=LS_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_xt[j][i]/=rhot[j][i];	//extract intermediate velocity field
			A_yt[j][i]/=rhot[j][i];
			if(abs(Ft[j][i])<=TRUNC) { Ft[j][i]=0.0; rhot[j][i]=rho_0; }	//empty cell
			else if(abs(Ft[j][i]-1.0)<=TRUNC) { Ft[j][i]=1.0; rhot[j][i]=rho_1; }	//completely filled cell
			if(Ft[j][i]<0.0) { Ft[j][i]=0.0; rhot[j][i]=rho_0; }
			else if(Ft[j][i]>1.0) { Ft[j][i]=1.0; rhot[j][i]=rho_1; }
		}
	}
	for(int i=1;i<=I;i++)
	{
		for(int j=1;j<=J;j++)
		{
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phit,v_NS[0][i]);	//upwind scheme for bottom boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no slip
				v_f[0]=0.0;	//no penetration
			}
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phit,v_NS[J][i]);	//upwind scheme for top boundary flux
				u_f[1]=A_xt[J][i];	//outflow
				v_f[1]=A_yt[J][i];	//outflow
				rho_adv[1]=rho_0*v_NS[J][i];	//outflow
			}
			else
			{
				up=(v_NS[j][i]<0.0)?(j+1):j;
				F_adv[1]=vol_frac_flux(1,i,up,Ft[up][i],Phit,v_NS[j][i]);
				LS_adv[1]=ENO2(1,i,j,Phit,v_NS[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*v_NS[j][i];
				if((i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use ENO2 scheme
				{
					u_f[1]=ENO2(1,i,j,A_xt,v_NS[j][i]);
					v_f[1]=ENO2(1,i,j,A_yt,v_NS[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(1,i,j,A_xt,v_NS[j][i]);
					v_f[1]=UP1(1,i,j,A_yt,v_NS[j][i]);
				}
			}
			F[j][i]=Ft[j][i]*(1.0+dt/dy*(v_NS[j][i]-v_NS[j-1][i]))-dt/dy*(F_adv[1]-F_adv[0]);
			rho_np1[j][i]=rhot[j][i]*(1.0+dt/dy*(v_NS[j][i]-v_NS[j-1][i]))-dt/dy*(rho_adv[1]-rho_adv[0]);
			Phi[j][i]=Phit[j][i]*(1.0+dt/dy*(v_NS[j][i]-v_NS[j-1][i]))-dt/dy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]);
			A_x[j][i]=rhot[j][i]*A_xt[j][i]*(1.0+dt/dy*(v_NS[j][i]-v_NS[j-1][i]))-dt/dy*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]);
			A_y[j][i]=rhot[j][i]*A_yt[j][i]*(1.0+dt/dy*(v_NS[j][i]-v_NS[j-1][i]))-dt/dy*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]);
			F_adv[0]=F_adv[1];	//updation for the next cell
			rho_adv[0]=rho_adv[1];
			LS_adv[0]=LS_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_x[j][i]/=rho_np1[j][i];	//extract velocity field
			A_y[j][i]/=rho_np1[j][i];
			if(abs(F[j][i])<=TRUNC) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }	//empty cell
			else if(abs(F[j][i]-1.0)<=TRUNC) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }	//completely filled cell
			if(F[j][i]<0.0) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }
			else if(F[j][i]>1.0) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }
		}
	}
	updt_ghost(F);	//update ghost cells of F
	updt_ghost(rho_np1);	//update ghost cells of rho
}
void CLSVOF::adv_Y()
{
	double F_adv[2],LS_adv[2],rho_adv[2],u_f[2],v_f[2];	//advection fluxes for volume fractions, LS, density, and momentum
	int up;	//donor cell index
	for(int i=1;i<=I;i++)
	{
		for(int j=1;j<=J;j++)
		{
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phi,v_NS[0][i]);	//upwind scheme for bottom boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no slip
				v_f[0]=0.0;	//no penetration
			}
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phi,v_NS[J][i]);	//upwind scheme for top boundary flux
				u_f[1]=A_x[J][i];	//outflow
				v_f[1]=A_y[J][i];	//outflow
				rho_adv[1]=rho_0*v_NS[J][i];	//outflow
			}
			else
			{
				up=(v_NS[j][i]<0.0)?(j+1):j;
				F_adv[1]=vol_frac_flux(1,i,up,F[up][i],Phi,v_NS[j][i]);
				LS_adv[1]=ENO2(1,i,j,Phi,v_NS[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*v_NS[j][i];
				if((i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use ENO2 scheme
				{
					u_f[1]=ENO2(1,i,j,A_x,v_NS[j][i]);
					v_f[1]=ENO2(1,i,j,A_y,v_NS[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(1,i,j,A_x,v_NS[j][i]);
					v_f[1]=UP1(1,i,j,A_y,v_NS[j][i]);
				}
			}
			Ft[j][i]=(F[j][i]-dt/dy*(F_adv[1]-F_adv[0]))/(1.0-dt/dy*(v_NS[j][i]-v_NS[j-1][i]));
			rhot[j][i]=(rho_n[j][i]-dt/dy*(rho_adv[1]-rho_adv[0]))/(1.0-dt/dy*(v_NS[j][i]-v_NS[j-1][i]));
			Phit[j][i]=(Phi[j][i]-dt/dy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]))/(1.0-dt/dy*(v_NS[j][i]-v_NS[j-1][i]));
			A_xt[j][i]=(rho_n[j][i]*A_x[j][i]-dt/dy*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]))/(1.0-dt/dy*(v_NS[j][i]-v_NS[j-1][i]));
			A_yt[j][i]=(rho_n[j][i]*A_y[j][i]-dt/dy*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]))/(1.0-dt/dy*(v_NS[j][i]-v_NS[j-1][i]));
			F_adv[0]=F_adv[1];	//updation for the next cell
			rho_adv[0]=rho_adv[1];
			LS_adv[0]=LS_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_xt[j][i]/=rhot[j][i];	//extract intermediate velocity field
			A_yt[j][i]/=rhot[j][i];
			if(abs(Ft[j][i])<=TRUNC) { Ft[j][i]=0.0; rhot[j][i]=rho_0; }	//empty cell
			else if(abs(Ft[j][i]-1.0)<=TRUNC) { Ft[j][i]=1.0; rhot[j][i]=rho_1; }	//completely filled cell
			if(Ft[j][i]<0.0) { Ft[j][i]=0.0; rhot[j][i]=rho_0; }
			else if(Ft[j][i]>1.0) { Ft[j][i]=1.0; rhot[j][i]=rho_1; }
		}
	}
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(0,1,j,Phit,u_EW[j][0]);	//upwind scheme for left boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no penetration
				v_f[0]=A_yt[j][1];	//free slip
			}
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phit,u_EW[j][I]);	//upwind scheme for right boundary flux
				rho_adv[1]=0.0;
				u_f[1]=0.0;	//no penetration
				v_f[1]=A_yt[j][I];	//free slip
			}
			else	//inner domain
			{
				up=(u_EW[j][i]<0.0)?(i+1):i;
				F_adv[1]=vol_frac_flux(0,up,j,Ft[j][up],Phit,u_EW[j][i]);
				LS_adv[1]=ENO2(0,i,j,Phit,u_EW[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*Xm[i]*u_EW[j][i];
				if((i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use ENO2 scheme
				{
					u_f[1]=ENO2(0,i,j,A_xt,u_EW[j][i]);
					v_f[1]=ENO2(0,i,j,A_yt,u_EW[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(0,i,j,A_xt,u_EW[j][i]);
					v_f[1]=UP1(0,i,j,A_yt,u_EW[j][i]);
				}
			}
			F[j][i]=Ft[j][i]*(1.0+dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]))-dt/(CX[i]*dx)*(F_adv[1]-F_adv[0]);
			rho_np1[j][i]=rhot[j][i]*(1.0+dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]))-dt/(CX[i]*dx)*(rho_adv[1]-rho_adv[0]);
			Phi[j][i]=Phit[j][i]*(1.0+dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]))-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]*LS_adv[1]-Xm[i-1]*u_EW[j][i-1]*LS_adv[0]);
			A_x[j][i]=rhot[j][i]*A_xt[j][i]*(1.0+dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]))-dt/(CX[i]*dx)*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]);
			A_y[j][i]=rhot[j][i]*A_yt[j][i]*(1.0+dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]))-dt/(CX[i]*dx)*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]);
			F_adv[0]=F_adv[1];	//updation for the next cell
			rho_adv[0]=rho_adv[1];
			LS_adv[0]=LS_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_x[j][i]/=rho_np1[j][i];	//extract velocity field
			A_y[j][i]/=rho_np1[j][i];
			if(abs(F[j][i])<=TRUNC) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }	//empty cell
			else if(abs(F[j][i]-1.0)<=TRUNC) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }	//completely filled cell
			if(F[j][i]<0.0) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }
			else if(F[j][i]>1.0) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }
		}
	}
	updt_ghost(F);	//update ghost cells of F
	updt_ghost(rho_np1);	//update ghost cells of rho
}
void CLSVOF::adv_trc()
{
	double F_adv[2],LS_adv[2];	//advection fluxes for volume fractions and LS fields
	int up;	//donor cell index
	for(int i=1;i<=I;i++)
	{
		for(int j=1;j<=J;j++)
		{
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phi_trc,v_NS[0][i]);	//upwind scheme for bottom boundary flux
			}
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phi_trc,v_NS[J][i]);	//upwind scheme for top boundary flux
			}
			else
			{
				up=(v_NS[j][i]<0.0)?(j+1):j;
				F_adv[1]=vol_frac_flux(1,i,up,F_trc[up][i],Phi_trc,v_NS[j][i]);
				LS_adv[1]=ENO2(1,i,j,Phi_trc,v_NS[j][i]);	//ENO2 scheme for inner domain
			}
			Ft[j][i]=(F_trc[j][i]-dt/dy*(F_adv[1]-F_adv[0]))/(1.0-dt/dy*(v_NS[j][i]-v_NS[j-1][i]));
			Phit[j][i]=(Phi_trc[j][i]-dt/dy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]))/(1.0-dt/dy*(v_NS[j][i]-v_NS[j-1][i]));
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			if(abs(Ft[j][i])<=TRUNC) { Ft[j][i]=0.0; }	//empty cell
			else if(abs(Ft[j][i]-1.0)<=TRUNC) { Ft[j][i]=1.0; }	//completely filled cell
			if(Ft[j][i]<0.0) { Ft[j][i]=0.0; }
			else if(Ft[j][i]>1.0) { Ft[j][i]=1.0; }
		}
	}
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(0,1,j,Phit,u_EW[j][0]);	//upwind scheme for left boundary flux
			}
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phit,u_EW[j][I]);	//upwind scheme for right boundary flux
			}
			else	//inner domain
			{
				up=(u_EW[j][i]<0.0)?(i+1):i;
				F_adv[1]=vol_frac_flux(0,up,j,Ft[j][up],Phit,u_EW[j][i]);
				LS_adv[1]=ENO2(0,i,j,Phit,u_EW[j][i]);	//ENO2 scheme for inner domain
			}
			F_trc[j][i]=Ft[j][i]*(1.0+dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]))-dt/(CX[i]*dx)*(F_adv[1]-F_adv[0]);
			Phi_trc[j][i]=Phit[j][i]*(1.0+dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]))-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]*LS_adv[1]-Xm[i-1]*u_EW[j][i-1]*LS_adv[0]);
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			if(abs(F_trc[j][i])<=TRUNC) { F_trc[j][i]=0.0; }	//empty cell
			else if(abs(F_trc[j][i]-1.0)<=TRUNC) { F_trc[j][i]=1.0; }	//completely filled cell
			if(F_trc[j][i]<0.0) { F_trc[j][i]=0.0; }
			else if(F_trc[j][i]>1.0) { F_trc[j][i]=1.0; }
		}
	}
}
void CLSVOF::LS_bsd_trunc()
{
	double h=dx;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if((F[j][i]>0.0)&&(F[j][i]<1.0))	//interfacial cell
			{
				if(Phi[j][i]>h) F[j][i]=1.0;
				else if(Phi[j][i]<-h) F[j][i]=0.0;
			}
		}
	}
}
void CLSVOF::solve(int n)
{
	if((n%2)==0) adv_Y();	//Strang splitting
	else adv_X();
	reinit(F,Phi);
	adv_trc();
	reinit(F_trc,Phi_trc);
}
void CLSVOF::mass_err()
{
	double mass=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass+=F[j][i]*CX[i]*dx*dy;
	cout<<"CLSVOF: MASS ERROR = "<<(mass-mass_act)/mass_act<<endl;
}
void CLSVOF::write_advt(int t)
{
	string fname="advt_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"INTERMEDIATE ADVECTION FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"adv_xt\",\"adv_yt\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<A_xt[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<A_yt[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: INTERMEDIATE ADVECTION FIELD FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
void CLSVOF::lsvf_write(int t)
{
	string fname="ls_vol_frac_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS AND VOLUME FRACTIONS\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"F\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<F[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: LEVEL SETS AND VOLUME FRACTIONS FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
void CLSVOF::ls_complete(int t)
{
	string fname="ls_comlpete_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS INCLUDING GHOST CELLS\""<<endl;
	p_out<<"FILETYPE = FULL"<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+2<<", J="<<J+2<<", DATAPACKING=BLOCK, SOLUTIONTIME="<<t*dt<<endl;
	double x=0.0-0.5*dx,y=0.0-0.5*dy;
	for(int j=0;j<=J+1;j++)	//X coordinates
	{
		for(int i=0;i<=I+1;i++)
		{
			p_out<<" "<<x;
			x+=dx;
		}
		p_out<<endl;
		x=0.0-0.5*dx;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)	//Y coordinates
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<y;
		p_out<<endl;
		y+=dy;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: COMPLETE LEVEL SET FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
void CLSVOF::trc_write(int t)
{
	string fname="trc_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS AND VOLUME FRACTIONS OF THE TRACERS\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"F_trc\",\"Phi_trc\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<F_trc[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<Phi_trc[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: TRACER FIELD FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
