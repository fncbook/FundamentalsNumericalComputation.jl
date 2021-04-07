using FundamentalsNumericalComputation
using Test,LinearAlgebra,OffsetArrays

@testset "Chapter 1" begin
	@test FNC.horner([-1,3,-3,1],1.6) ≈ 0.6^3
end

@testset "Chapter 2" begin
	A = [ 1 2 3 0; -1 1 2 -1; 3 1 2 4; 1 1 1 1 ]
	L,U = FNC.lufact(A)
	@test norm(L*U - A) < 100eps()
	@test norm(U - triu(U)) < 100eps()
	@test norm(L - tril(L)) < 100eps()
	b = [1,10,0,-1] / 5;
	@test norm(L\b - FNC.forwardsub(L,b)) < 100eps()
	@test norm(U\b - FNC.backsub(U,b)) < 100eps()
end

@testset "Chapter 3" begin
	A = [3 4 5;-1 0 1;4 2 0; 1 1 2; 3 -4 1]
	b = 5:-1:1
	@test FNC.lsnormal(A,b) ≈ A\b
	@test FNC.lsqrfact(A,b) ≈ A\b
	Q,R = qr(A)
	QQ,RR = FNC.qrfact(A)
	@test Q ≈ QQ
	@test R ≈ RR[1:3,:]
end

@testset "Chapter 4" begin

	for c = [2,4,7.5,11]
		f = x -> exp(x) - x - c;
		dfdx = x -> exp(x) - 1;
		x = FNC.newton(f,dfdx,1.0);  r = x[end];
		@test abs(f(r)) < 100eps()
	end

	for c = [2,4,7.5,11]
		f = x -> exp(x) - x - c;
		dfdx = x -> exp(x) - 1;
		x = FNC.secant(f,3,0.5);  r = x[end];
		@test abs(f(r)) < 100eps()
	end

	function nlfun(x)
		f = zeros(3)  
		f[1] = exp(x[2]-x[1]) - 2;
		f[2] = x[1]*x[2] + x[3];
		f[3] = x[2]*x[3] + x[1]^2 - x[2];
		return f
	end
	function nljac(x)
		J = zeros(3,3)
		J[1,:] = [-exp(x[2]-x[1]),exp(x[2]-x[1]), 0]
		J[2,:] = [x[2], x[1], 1]
		J[3,:] = [2*x[1], x[3]-1, x[2]]
		return J
	end

	x = FNC.newtonsys(nlfun,nljac,[0,0,0]);
	@test norm(nlfun(x[end])) < 100eps()
	x = FNC.newtonsys(nlfun,nljac,[1,2,3]);
	@test norm(nlfun(x[end])) < 100eps()

	x = FNC.levenberg(nlfun,[10,-4,-3])
	@test norm(nlfun(x[end])) < 1e-12

end

@testset "Chapter 5" begin
	f = t->cos(5t)
	Q,t = FNC.intadapt(f,-1,3,1e-8)
	@test Q ≈ (sin(15)+sin(5))/5 rtol = 1e-5
	T,_ = FNC.trapezoid(f,-1,3,820)
	@test T ≈ (sin(15)+sin(5))/5 rtol = 1e-4
	
	t = [-2,-0.5,0,1,1.5,3.5,4]/10
	S = FNC.spinterp(t,exp.(t))
	@test S(0.33) ≈ exp(0.33) rtol = 1e-5
	w = FNC.fdweights(t.-0.12,2)
	f = x->cos(3x)
	@test dot(w,f.(t)) ≈ -9cos(0.36) rtol = 1e-3
	y = FNC.hatfun(0.22,t,5)
	@test y ≈ (0.22-t[5])/(t[6]-t[5])
	@test FNC.hatfun(0.6,t,5)==0
	p = FNC.plinterp(t,f.(t)) 
	@test p(0.22) ≈ f(t[5]) + (f(t[6])-f(t[5]))*(0.22-t[5])/(t[6]-t[5])	
end

@testset "Chapter 6" begin
	f = (u,p,t) -> u + p*t^2
	û = exp(1.5) - 2*(-2 + 2*exp(1.5) - 2*1.5 - 1.5^2)
	ivp = ODEProblem(f,1,(0,1.5),-2)
	t,u = FNC.euler(ivp,4000)
	@test û ≈ u[end] rtol = 0.005
	t,u = FNC.am2(ivp,4000)
	@test û ≈ u[end] rtol = 0.005

	g = (u,p,t) -> [t+p-sin(u[2]),u[1]]
	ivp = ODEProblem(g,[-1.,4],(1.,2.),-6)
	sol = solve(ivp,Tsit5())
	t,u = FNC.euler(ivp,4000)
	@test u[end] ≈ sol.u[end] rtol=0.004
	t,u = FNC.ie2(ivp,4000)
	@test u[end] ≈ sol.u[end] rtol=0.0005
	t,u = FNC.rk4(ivp,800)
	@test u[end] ≈ sol.u[end] rtol=0.0005
	t,u = FNC.ab4(ivp,800)
	@test u[end] ≈ sol.u[end] rtol=0.0005
	t,u = FNC.rk23(ivp,1e-4)
	@test u[end] ≈ sol.u[end] rtol=0.0005
	t,u = FNC.am2(ivp,2000)
	@test u[end] ≈ sol.u[end] rtol=0.0005
end

@testset "Chapter 8" begin
	V = randn(4,4)
	D = diagm([-2,0.4,-0.1,0.01])
	A = V*D/V;
	
	γ,x = FNC.poweriter(A,30)
	@test γ[end] ≈ -2 rtol=1e-10
	@test abs( dot(x,V[:,1])/(norm(V[:,1])*norm(x)) ) ≈ 1 rtol=1e-10

	γ,x = FNC.inviter(A,0.37,15)
	@test γ[end] ≈ 0.4 rtol=1e-10
	@test abs( dot(x,V[:,2])/(norm(V[:,2])*norm(x)) ) ≈ 1 rtol=1e-10

	Q,H = FNC.arnoldi(A,ones(4),4)
	@test A*Q[:,1:4] ≈ Q*H

	x,res = FNC.gmres(A,ones(4),3)
	@test norm(ones(4) - A*x) ≈ res[end]
	x,res = FNC.gmres(A,ones(4),4)
	@test A*x ≈ ones(4)
end

@testset "Chapter 9" begin
	f = x -> exp(sin(x)+x^2)
	t = OffsetArray([-cos(k*π/40) for k in 0:40 ],0:40)
	p = FNC.polyinterp(t,f.(t))
	@test p(-0.12345) ≈ f(-0.12345)

	f = x -> exp(sin(π*x))
	n = 30
	t = [ 2k/(2n+1) for k in -n:n ]
	p = FNC.triginterp(t,f.(t))
	@test p(-0.12345) ≈ f(-0.12345)
	t = [ k/n for k in -n:n-1 ]
	p = FNC.triginterp(t,f.(t))
	@test p(-0.12345) ≈ f(-0.12345)

	F = x -> tan(x/2-0.2)
	f = x -> 0.5*sec(x/2-0.2)^2
	@test FNC.ccint(f,40)[1] ≈ F(1)-F(-1)
	@test FNC.glint(f,40)[1] ≈ F(1)-F(-1)

	f = x -> 1/(32+2x^4)
	@test FNC.intde(f,.2,20)[1] ≈ sqrt(2)*π/32 rtol=1e-5

	f = x -> 1/( sin(1+x)^0.5*(1-x)^0.25 )
	@test FNC.intsing(f,.1,1e-7)[1] ≈ 3.16762 rtol=1e-4
end

@testset "Chapter 10" begin
	λ = 0.6
	phi = (r,w,dwdr) -> λ/w^2 - dwdr/r;
	a = eps();  b = 1;
	lval = [];  lder = 0;   # w(a)=?, w'(a)=0
	rval = 1;   rder = [];  # w(b)=1, w'(b)=?

	r,w,dwdx = FNC.shoot(phi,(a,b),lval,lder,rval,rder,0.8)
	@test w[1] ≈ 0.78775 rtol=1e-4

	f = x -> exp(x^2-3x)
	df = x -> (2x-3)*f(x)
	ddf = x -> ((2x-3)^2+2)*f(x)

	t,D,DD = FNC.diffmat2(400,(-0.5,2))
	@test df.(t) ≈ D*f.(t) rtol=1e-3
	@test ddf.(t) ≈ DD*f.(t) rtol=1e-3
	t,D,DD = FNC.diffcheb(80,(-0.5,2))
	@test df.(t) ≈ D*f.(t) rtol=1e-7
	@test ddf.(t) ≈ DD*f.(t) rtol=1e-7

	exact = x -> exp(sin(x));
	p = x -> -cos(x);
	q = sin;
	r = x -> 0; 
	x,u = FNC.bvplin(p,q,r,[0,pi/2],1,exp(1),300);
	@test u ≈ exact.(x) rtol=1e-3

	ϕ = (t,θ,ω) -> -0.05*ω - sin(θ);
	init = collect(LinRange(2.5,-2,101));
	lval,lder = 2.5,[]
	rval,rder = -2,[]

	t,θ = FNC.bvp(ϕ,[0,5],lval,lder,rval,rder,init)
	@test θ[end][7] ≈ 2.421850016880724 rtol=1e-10

	c = x -> x^2;
	q = x -> 4;
	f = x -> sin(π*x);
	x,u = FNC.fem(c,q,f,0,1,100)
	@test u[33] ≈ 0.1641366907307196 rtol=1e-10
end

@testset "Chapter 11" begin
	s = x -> sin(π*(x-0.2))
	c = x -> cos(π*(x-0.2))
	f = x -> 1 + s(x)^2
	df = x -> 2π*s(x)*c(x)
	ddf = x -> 2π^2*(c(x)^2 - s(x)^2)

	t,D,DD = FNC.diffper(400,(0,2))
	@test df.(t) ≈ D*f.(t) rtol=1e-3
	@test ddf.(t) ≈ DD*f.(t) rtol=1e-3
end

@testset "Chapter 13" begin 
	X,Y = FNC.ndgrid(-3:3,0:4)
	@test [X[3,2],Y[3,2]] ≈ [-1,1]

	f = x -> exp(x^2-3x)
	df = x -> (2x-3)*f(x)
	ddf = x -> ((2x-3)^2+2)*f(x)
	X,Y,d = FNC.rectdisc(100,(0.1,0.5),80,(-0.3,-0.2))
	t = X[:,1]
	@test df.(t) ≈ d.Dx*f.(t) rtol=1e-4
	t = Y[1,:]
	@test df.(t) ≈ d.Dy*f.(t) rtol=1e-3

	f = (x,y) -> -sin(3*x.*y-4*y)*(9*y^2+(3*x-4)^2);
	g = (x,y) -> sin(3*x*y-4*y);
	xspan = [0,1];  yspan = [0,2];
	U,X,Y = FNC.poissonfd(f,g,50,xspan,80,yspan);
	@test g.(X,Y) ≈ U rtol=1e-3

	λ = 1.5
	function pde(U,X,Y,d)
		LU = d.Dxx*U + U*d.Dyy';     # apply Laplacian
		F = @. LU - λ/(U+1)^2   # residual
	
		L = kron(d.Dyy,d.Ix) + kron(d.Iy,d.Dxx)  # Laplacian matrix
		u = d.vec(U)
		J = L + spdiagm( @. 2λ/(u+1)^3 ) 
		return F,J
	end      
	g = (x,y) -> 0     # boundary condition
	U,X,Y = FNC.newtonpde(pde,g,100,[0,2.5],80,[0,1]);
	@test U[20,45] ≈ -0.25732335309199816
end