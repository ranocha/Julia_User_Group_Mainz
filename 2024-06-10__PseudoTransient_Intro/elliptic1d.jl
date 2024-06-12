using Plots, SpecialFunctions
#Solve the 1d Elliptic Problem
#K*d2u/dx2=Q
#with dirichlet BC
#and homogeneous source (Q)
#Examples of usage are the steady state heat equation
#or the steady state chemical diffusion equation etc
function analytic(xan,K,Q,uBC) #analytic solution
    Lt  = xan[end]-xan[1];
    C1  = (uBC[2]-uBC[1])/Lt;
    C2  = 0.125*(4.0*K*uBC[1] + 4.0*K*uBC[2] + Lt^2*Q)/K;
    uan = -0.5*Q/K .*xan.^2 .+ C1 .* xan .+ C2;
    return uan
end
function define_problem(L,nx,K,Q,uBC)
    #Impose types - this is not necessary but nsures that 
    #the variables are of the correct type
    L    = Float64(L)       #Length of the domain
    nx   = Int64(nx)        #Number of grid points
    K    = Float64(K)       #Diffusion Coefficient
    Q    = Float64(Q)       #Source term
    uBC  = Float64.(uBC)    #Boundary conditions
    #Create Grid
    dx   = L/(nx-1);        #Grid spacing
    x    = -L/2:dx:L/2;     #Grid
    #Initialize Variables
    q    = zeros(Float64,nx-1)  #Fluxes
    u    = zeros(Float64,nx)    #Solution (Numerical)
    dudt = zeros(Float64,nx-2)  #Derivative of the residual
    return x,q,u,dudt,dx,uBC,K,Q,nx
end
function ptr1!()
#Introduce Problem parameters
L    = 10;      #Length of the domain
nx   = 100;     #Number of grid points
K    = 1;       #Diffusion Coefficient
Q    = 1;       #Source term
uBC  = [1,2];   #Boundary conditions
#Define the problem
x,q,u,dudt,dx,uBC,K,Q,nx    = define_problem(L,nx,K,Q,uBC)
#Numerical parameters
maxit  = 100000         #Maximum number of iterations
nout   = 100;           #Output every nout iterations
etol   = 1e-10;         #Tolerance
CFL    = 0.4;           #CFL number
dt     = dx^2/K*CFL;    #Time step
#Initialize BC
u[1]   = uBC[1];        #Left BC
u[end] = uBC[2];        #Right BC
err    = 1e23;          #Error initialization
logErr = log10(err)     #Log of the error
it     = 0;
while it < maxit
    it += 1 #Update iteration
    q          .= .- K.* diff(u)./dx;           #Fluxes
    dudt       .= .- diff(q)/dx   .+ Q;         #time derivative of residual
    u[2:end-1] .= u[2:end-1]      .+ dt .*dudt  #Update solution
    err         = maximum(abs.(dudt))           #Error
    logErr      = log10(err)                    #Log of the error
    if rem(it,nout) == 1
        println("Iteration:  $it  - Err: $logErr") #Output
    end
    if err<etol
        return x,u,err,it   #Exit if error is below tolerance
        break               #Break the loop
    end
end
end
x,u,err,it= ptr1!()             #Call the function
uan = analytic(x,1,1,[1,2])     #Analytic solution
#plot(x,[u,uan])                #Plot the solution