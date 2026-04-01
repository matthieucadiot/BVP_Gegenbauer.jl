using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, JLD2, MAT, MATLAB
include("list_functions.jl")

############ Initialization ##########

N = 100
α = 1
ξ1 = eval_boundary1(N,0) ### evaluation of the boundary at x=1
ξ2 = eval_boundarym1(N,0) ### evaluation of the boundary at x=-1
U = Sequence(Chebyshev(N),zeros(N+1))

############ Newton #######################

U = Sequence(Chebyshev(N),zeros(N+1))
ϵ = 1e-12
kmax = 50
U = Newton_cheb(U,α,N)


########### Existence Proof ##################

conv = mid.(conversion(1,N))*mid.(conversion(0,N))
D2 =  mid.(Deriv(1,N))*mid.(Deriv(0,N)) 

D2 = inv(conv)*D2
V = U 
U = LinearOperator(Chebyshev(N),Chebyshev(N),D2)*V

U = interval.(mid.(U))
α = interval(α)
ν = interval(3.5)
#### aproximate inverse #######################

DF = DF_KS(U,α,N)
A = interval.(inv(mid.(project(DF,Chebyshev(N),Chebyshev(N)))))

####### bound Y ##################

F = F_KS(U,α,N)
Y = norm(weight_nu(ν,N).*(A*F),1) + norm(weight_nu(ν,N+3).*(F[N+1:end]),1)
println("Y = $Y")

########### bound Z1 ####################

Aext = [coefficients(A) interval.(zeros(N+1,N+4));interval.(zeros(N+4,N+1)) interval.(I(N+4))] 
Aext = LinearOperator(Chebyshev(2*N+4), Chebyshev(2*N+4), Aext)

ADF = Aext*DF
Z10_Z11 = opnorm(weight_nu(ν,2N+4).*(LinearOperator(Chebyshev(N+2), Chebyshev(N+2), interval.(I(N+3))) - ADF) ./ weight_nu(ν,N+2)',1)

V0 = interval(2)*α*operator_Linv(N)*U

η0N = interval(0.25)/(interval(N-1)*interval(N)*ν) + interval(0.5)/(interval(N+1)^2-interval(1)) +  interval(0.25)*ν/(interval(N+2)*interval(N+3)) + interval(3)*ν^(-interval(N))/((interval(N+1)^2-interval(4))*(interval(0.5)*interval(N+1) - interval(0.5)))

norm_A = maximum([interval(1) opnorm(weight_nu(ν,N).*A./weight_nu(ν,N)',1)]) 
Z12_Z13 = norm(weight_nu(ν,N+2).*V0,1)*η0N*(interval(1) + norm_A)
#### the above is less sharp than the formulas provided in the paper, but easier to compute

Z1 = maximum([Z10_Z11 Z12_Z13])
println("Z1 = $Z1")


############ bound Z2 ####################
Z2 = interval(0.5)*abs(α)*norm_A
println("Z2 = $Z2")


########## verification of the conditions ###########
if inf((interval(1)-Z1)^2) > sup(interval(2)*Z2*Y)
    r0 = (interval(1)-interval(sup(Z1)) - sqrt((interval(1)-interval(sup(Z1)))^2-interval(2)*Z2*Y))/Z2
    if sup(Z1 + Z2*abs(r0)) < 1
        display("proof successful for r =")
        display(sup(r0))
        return sup(r0)
    else
        display("second condition not verified")
        return naN
    end
else
    display("first condition not verified")
    return naN
end  



################## STABILITY PART ##########################

Aext = [coefficients(A) interval.(zeros(N+1,2));interval.(zeros(2,N+1)) interval.(I(2))] 
Aext = LinearOperator(Chebyshev(N+2), Chebyshev(N+2), Aext)
M = Aext*operator_Linv(N)
M = M[0:N,0:N]


D,P = eigen(mid.(M))

P = interval.(mid.(P))
Pinv = interval.(inv(mid.(P)))

𝒟 = Pinv*M*P 

Linv = operator_Linv(N)
ADFP = ADF[:,0:N]*P
ALinvP = (Aext*Linv)[:,0:N]*P 

nP = opnorm(P,1)
nPinv = opnorm(Pinv,1)

#### we compute an upper bound and a lower bound on a finite number of eigenvalues, that we call max_vp and min_vp
max_vp = interval(-2)
min_vp = interval(100)
    for j=1:N+1
        global V, max_vp, min_vp
        V = interval.(mid.(𝒟[:,j]))
        Cj = [ALinvP[:,j]; interval.(zeros(N+2))]
        Yj = norm(ADFP*V - Cj,1) + Z2*r0*norm(P*V,1)
        ϵj = nPinv*Yj/(interval(1)-(Z1+Z2*r0))
        vj = V[j]
       
        V[j] = interval(0)
        rad = interval(2)*norm(V,1) + ϵj
        max_vp = maximum([vj+rad max_vp])
        min_vp = minimum([vj-rad min_vp])
    end 


############ tail estimate ##############    
norm_A = maximum([interval(1) opnorm(A,1)])

### we compute very rough estimations for the bounds of Proposition 4.5. It turns out that N is big enough for these quantities to be small anyways.
ϵ1 = norm_A*norm(V0,1)*η0N/(interval(N+1)*interval(N+2)) 
ϵ2 = norm_A*norm(V0,1)*opnorm(Linv[1:2,1:2],1)*η0N

ϵN = interval(0.5)*nPinv/(interval(1)-(Z1+Z2*r0))*(Z2*r0*η0N + ϵ1 + ϵ2)
rN = nPinv*η0N 

dist_zero = maximum([rN+ϵN max_vp]) #### we compute the upper bound on the spectrum 
println("value to zero for the tail :  $dist_zero")

C = abs(α)*(norm(U,1)+r0)

if inf(interval(1)/C) > sup(dist_zero)
    display("All the eigenvalues are stricly negative")
    display("Moreover, an upper bound on the eigenvalues is")
    display(sup(interval(1)/min_vp))
else
    display("cannot conclude for the sign of the eigenvalues")
end


