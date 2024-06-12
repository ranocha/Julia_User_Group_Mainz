using Plots, SpecialFunctions
#Solve the 1d Elliptic Problem
#d(K(u) du/dx)/dx=Q
#with dirichlet BC
#homogeneous source (Q)
#Take K = K0 + gam + u^npow as effective K 
function av1d!(uc,u)
    for i in eachindex(uc)
        uc[i] = 0.5*u[i] + 0.5*u[i+1];
    end
end
function analytic(xan,K,Q,uBC)
    Lt  = xan[end]-xan[1];
    C1  = (uBC[2]-uBC[1])/Lt;
    C2  = 0.125*(4.0*K*uBC[1] + 4.0*K*uBC[2] + Lt^2*Q)/K;
    uan = -0.5*Q/K .*xan.^2 .+ C1 .* xan .+ C2;
    return uan
end
function define_problem(L,nx,K,Q,uBC,gam,npow)
    #Impose types
    L    = Float64(L)
    nx   = Int64(nx)
    K    = Float64(K)
    Q    = Float64(Q)
    uBC  = Float64.(uBC)
    gam  = Float64(gam)
    npow = Float64(npow)
    #Create Grid
    dx   = L/(nx-1);
    x    = -L/2:dx:L/2;
    #Initialize Variables
    q    = zeros(Float64,nx-1)
    u    = zeros(Float64,nx)
    dudt = zeros(Float64,nx-2)
    uc   = zeros(Float64,nx-1)
    Keff = K .+ gam .* uc.^npow;
    return x,q,u,dudt,dx,uBC,K,Q,nx,uc,gam,npow, Keff
end
function ptr1!()
#Introduce Problem parameters
L    = 10;
nx   = 100;
K    = 1;
Q    = 1;
gam  = 1;
npow = 2;
uBC  = [1,2];
etol = 1e-9;
x,q,u,dudt,dx,uBC,K,Q,nx,uc, gam, npow, Keff   = define_problem(L,nx,K,Q,uBC,gam,npow)
#Numerical parameters
maxit  = 100000
nout   = 100;
CFL    = 0.3;
dt     = dx^2/K*CFL;
damp   = 1.0 - 5.0/nx
#Initialize BC
u[1]   = uBC[1];
u[end] = uBC[2];
err    = 1e23;
logErr = log10(err)
it     = 0.0;
while it < maxit
    it += 1
    av1d!(uc,u);
    Keff        =  K .+ gam .* uc.^npow;
    dt          =  dx^2/maximum(Keff)*CFL;
    q          .= .- Keff .* diff(u)./dx;
    dudt       .= .- diff(q)/dx   .+ Q + (damp .* dudt);
    u[2:end-1] .= u[2:end-1]      .+ dt .*dudt
    err         = maximum(abs.(dudt))
    logErr      = log10(err)
    if rem(it,nout) == 1
        println("Iteration:  $it  - Err: $logErr")
    end
    if err<etol
        return x,u,err,it
        break
    end
end
end
x,u,err,it= ptr1!()
uan = analytic(x,1,1,[1,2])
plot(x,[u,uan])