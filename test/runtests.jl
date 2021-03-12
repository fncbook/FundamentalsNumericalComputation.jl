using FundamentalsNumericalComputation
using Test

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