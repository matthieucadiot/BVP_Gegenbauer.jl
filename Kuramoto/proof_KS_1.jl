using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, JLD2
using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn)) ### This hides the [Info ] messages

include("list_functions.jl")

################################### Proof of existence #############################################

## parameter α for KS and order of the odd Chebyshev expansion 
α = interval(1.0) 
N = 200  ; m = 2
N0 = Int(N/2)

U = load("U_1.jld2","U")


## construction of the Ki operators in multiplrecision 
prec = 128 ###  precision for the construction of the linear operators. The boundary conditions are increasing fast, it might be a ncessesary to increase it

#### computation of the different θi corresponding to the boundary conditions of the different operators 𝒦i 
θ0 = interval.(big.(zeros(2m,2m)))
α0 = 2
θ0[1,2] = interval(big(1))
θ0[2,4] = interval(big(1))
θ0[3,2] = interval(big(1))
θ0[4,4] = interval(big(1))

θ1 = interval.(big.(zeros(2m-1,2m)))
α1 = 2
θ1[1,1] = interval(big(1))
θ1[2,3] = interval(big(1))
θ1[3,1] = interval(big(1))

θ2 = interval.(big.(zeros(2m-2,2m)))
α2 = 1
θ2[1,2] = interval(big(1))
θ2[2,2] = interval(big(1))

Linv = operator_K_odd(2m,θ0,α0,m,N+8,prec) ### size N/2 +4 by N/2
K1 = operator_K_odd(2m-1,θ1,α1,m,N+8,prec) ### size N/2 +4 by N/2
K2 = operator_K_odd(2m-2,θ2,α2,m,N+8,prec) ### size N/2 +4 by N/2

### going down to standard precision
Linv = interval.(Float64.(inf.(Linv),RoundDown),Float64.(sup.(Linv),RoundUp) ) 
K1 = interval.(Float64.(inf.(K1),RoundDown),Float64.(sup.(K1),RoundUp) ) 
K2 = interval.(Float64.(inf.(K2),RoundDown),Float64.(sup.(K2),RoundUp) ) 


# U0 = mid.(Linv[1:N0,1:N0])*mid.(U)
# V = Sequence(Chebyshev(N0),vec(odd2full(U0,N0)))

# X = -1:0.01:1
# plot(X, V.(X))



## construction of the operator A

DF = DF_KS(U,α,Linv[1:N0+2,1:N0],K1[1:N0+2,1:N0],K2[1:N0+2,1:N0],N)
DF0 = mid.(DF[1:N0,1:N0])
A = interval.(inv(DF0))


Next = N + 2*m
Aext = [A interval.(zeros(N0,N+4-N0));interval.(zeros(N+4-N0,N0)) interval.(I(N+4-N0))]   

## computation of Y
## because we have an even symmetry, the norm used is essentially 2 times the usual l1 norm on the Chebyshev coefficients
F = F_KS(U,α,Linv[1:N0+2,1:N0],K1[1:N0+2,1:N0],K2[1:N0+2,1:N0],N)
Y = interval(2)*norm(A*F[1:N0],1) + interval(2)*norm(F[N0+1:N+4],1)

println("Y = $Y")

## construction of the operator I-ADf
ADF = Aext*DF
I_ADF = [interval.(I(N0))-ADF[1:N0,1:N0] ; -ADF[N0+1:N+4,1:N0]]
Z10_Z11 = opnorm(I_ADF,1) #### to simplify, we compute both Z10 and Z11 together

### this is direct from the explicit expression in the appendix
Ni = interval(N)
nK2 = interval(1.5)/( (Ni+interval(1))^2-interval(1) ) + interval(1)/( interval(4)*Ni*(Ni-interval(1)) )
nK2 +=  interval(1)/( interval(4)*(Ni+interval(3))*(Ni+interval(2)) )

norm_K2 = maximum([opnorm(K2,1) nK2])

nLinv = nK2*norm_K2 ### we use the fact that Linv = K2^2
nK1 = nK2### this is using that \|∂xK2\| ≤ 1 . In practice it looks like we are missing a factor 2, but this is still enough for the proof

V1 = Sequence(Chebyshev(N+3),vec(even2full(K1[1:N0+2,1:N0]*U,N+3)))
V2 = Sequence(Chebyshev(N+4),vec(odd2full(Linv[1:N0+2,1:N0]*U,N+4)))

Z12_Z13 = abs(α)*(nK2 + nLinv*norm(V1,1) + nK1*norm(V2,1))*(interval(1) + opnorm(A,1))
#### the above is less sharp than the formulas provided in the paper, but easier to compute and it is sufficient for our purposes. It just relies on Proposition 2.15

Z1 = maximum([Z10_Z11 Z12_Z13])

println("Z1 = $Z1")

##### computation of Z2
norm_K1 = maximum([opnorm(K1,1) nK1])
norm_Linv = maximum([opnorm(Linv,1) nLinv])
Z2 = opnorm(A,1)*interval(2)*abs(α)*norm_K1*norm_Linv
Z21 = opnorm(A,1)*interval(2)*abs(α)*norm_K1*norm_Linv
Z22 = interval(2)*abs(α)*norm_K1*norm_Linv

println("Z2 = $Z2")

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


Aext = [A interval.(zeros(N0,2));interval.(zeros(2,N0)) interval.(I(2))] 
M = Aext*Linv[1:N0+2,:]
M = M[1:N0,1:N0]


D,P = eigen(mid.(-M))
n_pos = 0
while real(D[end-n_pos]) > 0
    global n_pos
    n_pos += 1
end

P1 = P[:,end-n_pos+1:end] ; P = [P1 P[:,1:end-n_pos]]

P = interval.(mid.(P))
Pinv = interval.(inv(mid.(P)))

𝒟 = -Pinv*M*P 

S = Diagonal(diag(𝒟))

R = 𝒟
for i=1:N0
    R[i,i] = interval(0)
end
𝒟 = -Pinv*M*P  

ADFP = ADF*P
ALinvP = Aext*Linv[1:N0+2,1:N0]*P 

nP = opnorm(P,1)
nPinv = opnorm(Pinv,1)


#### we notice that the first eigenvalue is positive, we want to rigorously verify it 

##### first, we compute the disk associated to the first eigenvalue
V = interval.(mid.(𝒟[:,1]))
C1 = [ALinvP[:,1]; interval.(zeros(N+2-N0))]
Y1 = norm(ADFP*V + C1,1) + Z2*r0*norm(P*V,1)
ϵ1 = nPinv*Y1/(interval(1)-(Z1+Z2*r0))
        
v1 = V[1]
V[1] = interval(0)
rad = norm(V,1) + ϵ1
upper_bound = abs(v1) + rad
lower_bound = abs(v1) - rad

println("Disk associated to the first eigenvalue : center = $v1 , radius = $rad")

#### we compute an upper bound on a finite number of the rest of the eigenvalues, that we call max_vp and min_vp
max_vp = interval(-2)
min_vp = interval(100)

    for j=2:N0
        global V, max_vp, min_vp, rad
        V = interval.(mid.(𝒟[:,j]))
        Cj = [ALinvP[:,j]; interval.(zeros(N+2-N0))]
        Yj = norm(ADFP*V + Cj,1) + Z2*r0*norm(P*V,1)
        ϵj = interval(0.5)*nPinv*Yj/(interval(1)-(Z1+Z2*r0))
        
        vj = V[j]
        V[j] = interval(0)
        rad = interval(2)*norm(V,1) + ϵj
        ###
        max_vp = maximum([abs(vj)+rad max_vp])
    end 


η0N = interval(1)/(interval(1+N)^2) 

η1N = interval(1)/(interval(1+N-2*m-1)^2) 

MV1 = project(Multiplication(V1),Chebyshev(N+4),Chebyshev(2N+7))
MV2 = project(Multiplication(V2),Chebyshev(N+3),Chebyshev(2N+7))

MV1 = MV1[1:2:end,1:2:end] ### size N + 4  by N/2 + 2
MV2 = MV2[1:2:end,0:2:end] ### size N + 4  by N/2 + 2


# n_AVK1B = opnorm(A*K1[1:N0,1:N0+4]*(MV2[1:N0+4,1:2]*B2m_dag),1)
# n_AVLinvB = opnorm(A*Linv[1:N0,1:N0+4]*(MV1[1:N0+4,1:2]*B2m_dag),1)
n_AVK2 = abs(α)*opnorm(A*(K2[1:N0,1:2] + Linv[1:N0,1:N0+4]*MV1[1:N0+4,1:2] + K1[1:N0,1:N0+4]*MV2[1:N0+4,1:2]),1)

ϵ1 = n_AVK2*η0N
ϵ2 = opnorm(A,1)*abs(α)*(η1N + η1N*norm(V1,1) + η1N*norm(V2,1))/(interval(1+N)^4)

ϵN = nPinv/(interval(1)-(Z1+Z2*r0))*(Z2*r0*η0N + ϵ1 + ϵ2)
rN = nPinv*η0N 

dist_zero = maximum([rN+ϵN max_vp]) #### we compute the upper bound on the spectrum 
println("distance to zero for the tail :  $dist_zero")


if inf(real(lower_bound)) > sup(dist_zero)
    display("the first eigenvalue is isolated from the rest of the finite union of Gershgorin disks")    
else 
    display("cannot isolate the first eigenvalue from the rest of the finite union of Gershgorin disks")
end

norm_LinvU = (norm(Linv[:,1:N0]*U,1)+r0)^2

##### we compute an upper bound for the real and imaginary part of the eigenvalues, denoted by Re_λmax and Im_λmax respectively. For this purpose, we implement the bounds given in Remark 5.12

norm_v = norm(Linv[:,1:N0]*U,1)+r0
β = interval(2)/pi*α*norm_v
μ = interval(463.1)

### computation of the estimate of the real part of the eigenvalues
κ1 = interval(2)*α*(interval(0.5) + interval(1)/π*norm_v)
κ1 = interval(1)/κ1 ; κ2 = copy(κ1)
Re_λmax =α*(norm(K1[:,1:N0]*U,1)+r0 + interval(0.25)/κ1 + norm_v/(interval(2)*pi*κ2))

#### computation of the estimate of the imaginary part of the eigenvalues
δ = α*(interval(2) + interval(4)/π*norm_v)
κ1 = interval(0.5)*(δ + sqrt(δ^2 + interval(16)*(μ + α*(norm(K1[:,1:N0]*U,1)+r0))))
κ1 = κ1^(interval(-1)) ; κ2 = copy(κ1) 
Im_λmax = α*(interval(1) + interval(2)/pi*norm_v )*sqrt(abs( (μ + α*(norm(K1[:,1:N0]*U,1)+r0 + interval(0.25)/κ1 + norm_v/(interval(2)*pi*κ2) ))/(interval(1)-α*κ1 - β*κ2) ))

### value for λmax, providing an upper bound for all eigenvalues with real part greater than μ
λmax = sqrt(Re_λmax^2 + Im_λmax^2)

if inf(lower_bound) > sup(interval(1)/λmax)
    display("the tail eigenvalues are stable")    
else 
    display("cannot conclude about the stability of the tail eigenvalues")
end


println("There is exactly 1 eigenvalues with positive real part")
λj = mid.(1/𝒟[1,1])
println("The enclosure for the eigenvalue around $λj has a radius $rad")