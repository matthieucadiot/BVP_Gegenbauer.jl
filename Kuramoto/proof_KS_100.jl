using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, JLD2, MATLAB
using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn)) ### This hides the [Info ] messages

include("list_functions.jl")

################################### Proof of existence #############################################

## parameter α for KS and order of the odd Chebyshev expansion 
# α = interval(1.0) 
# N = 200  ; m = 2
# N0 = Int(N/2)

# U = load("U_1.jld2","U")


## parameter α for KS and order of the odd Chebyshev expansion 
α =  interval(100.0)
N = 200  ; m = 2
N0 = Int(N/2)

U = load("U_001.jld2","U")



## construction of the Ki operators in multiplrecision 
prec = 256 ###  precision of the computation

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


U0 = mid.(Linv[1:N0,1:N0])*mid.(U)
V = Sequence(Chebyshev(N0),vec(odd2full(U0,N0)))



# SaveMatlabFig(V, "plot_KS100.fig")
    

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


Ni = interval(N)
nK2 = interval(1.5)/( (Ni+interval(1))^2-interval(1) ) + interval(1)/( interval(4)*Ni*(Ni-interval(1)) )
nK2 +=  interval(1)/( interval(4)*(Ni+interval(3))*(Ni+interval(2)) )

norm_K2 = maximum([opnorm(K2,1) nK2])

nLinv = nK2*norm_K2 ### we use the fact that Linv = K2^2
nK1 = nK2### this is using that \|∂xK2\| ≤ 1.

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



M = A*Linv[1:N0,1:N0]

D,P = eigen(mid.(-M))
n_pos = 0
while real(D[end-n_pos]) > 0
    global n_pos
    n_pos += 1
end

#### we reorganize the eigenvalues and eigenvectors so that the ones with positive real part are first
P1 = P[:,end-n_pos+1:end] ; P = [P1 P[:,1:end-n_pos]]
P = interval.(mid.(P))
Pinv = interval.(inv(mid.(P)))

𝒟 = -Pinv*M*P 

R = 𝒟
for i=1:N0
    R[i,i] = interval(0)
end
𝒟 = -Pinv*M*P  

##### we compute an upper bound for the real and imaginary part of the eigenvalues, denoted by Re_λmax and Im_λmax respectively. For this purpose, we implement the bounds given in Remark 5.12

norm_v = norm(Linv[:,1:N0]*U,1)+r0
β = interval(2)/pi*α*norm_v
μ = interval(377)

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


############ Other proof of stability 

ϵ =   abs(α)*(interval(1)/interval((N+2)*(N+3)) + norm(V1,1)/interval((N+2)*(N+3)*(N+4)*(N+5)) + norm(V2,1)/interval((N+2)*(N+3)*(N+4)) + interval(2)*nLinv*norm(V1[N-4:end],1) + interval(2)*nK1*norm(V2[N-4:end],1) ) + λmax/interval((N+2)*(N+3)*(N+4)*(N+5))


MV1 = project(Multiplication(V1),Chebyshev(N+4),Chebyshev(2N+7))
MV2 = project(Multiplication(V2),Chebyshev(N+3),Chebyshev(2N+7))

MV1 = MV1[1:2:end,1:2:end] ### size N + 4  by N/2 + 2
MV2 = MV2[1:2:end,0:2:end] ### size N + 4  by N/2 + 2

K2_ext = [K2[1:N0+2,1:N0];interval.(zeros(N+3 - (N÷2 +1),N÷2))] ### size N + 4 by N/2
DK_ext =  α*K2_ext + α*(MV2*K1[1:N0+2,1:N0]) + α*(MV1*Linv[1:N0+2,1:N0])

err_Pinv = opnorm(interval.(I(N0)) - Pinv*P,1) ## Pinv is not the exact inverse of P, consequently we propagate this error in the estimates. We use a Neumann series to bound the error
norm_P = opnorm(P,1)
norm_Pinv = opnorm(Pinv,1)*interval(1)/(interval(1) - err_Pinv)
norm_DK_P = opnorm(DK_ext[N0+1:end,1:N0]*P,1)
β020 = opnorm(Pinv*A,1)*abs(α)*(nK2 + nLinv*norm(V1,1) + nK1*norm(V2,1))

β01 = norm_DK_P + norm_P*Z22*r0
β02 = β020*β01
β11 = nLinv*opnorm(Pinv*A,1)*β01*interval(1)/(interval(1) - err_Pinv)
β12 = opnorm(Linv[N0+1:N0+4,1:N0]*P,1)
β13 = β020*β12
β14 = β020 + norm_Pinv*Z21*r0
β21 = opnorm(Pinv*A*Linv[1:N0,N0+1:N0+4]*Linv[N0+1:N0+4,1:N0]*P,1)*interval(1)/(interval(1) - err_Pinv)
β22 = nLinv*opnorm(Pinv*A,1)*interval(1)/(interval(1) - err_Pinv)

Z1P = opnorm(interval.(I(N0)) - Pinv*A*DF[1:N0,1:N0]*P,1)*interval(1)/(interval(1) - err_Pinv)
Z2P = norm_P*norm_Pinv*Z2


β0 = β02 + norm_Pinv*Z21*r0*β01 + ϵ/(interval(1)-ϵ)*β01*β14
β1 = β11 + β13 + norm_Pinv*Z21*r0*β12 + ϵ/(interval(1)-ϵ)*(β12*β14+ β01*β22)
β2 = β21 + ϵ/(interval(1)-ϵ)*β12*β22

δ = λmax*opnorm(R,1) + β0 + β1*λmax + β2*λmax^2 + Z1P + Z2P*r0

##### if λi > λmax/(1-R), then the corresponding disk contains a negative eigenvalue

### we loop over the eigenvalues and count how many have a positive real part such that their respective Gershgorin disk is completely contained in the right half plane
min_pos = interval(0)
nb_pos = 0
index_pos = []
   if sup(δ) < 1
        for j=1:N0
         global D, nb_pos, min_pos, δ, λmax, index_pos, λj
         μj = 𝒟[j,j]
            if abs(μj) >= (interval(1)-δ)/λmax
            λj = interval(1)/μj
                if inf(real(λj)-δ*abs(λj)) > 0
                nb_pos += 1 
                push!(index_pos,j)
                min_pos = minimum([min_pos real(λj)-δ*abs(λj)])
                end
            end
        end
    end

    ### for the above obtained positive Gershogorin disks, we ensure that they are disjoint from the rest of the disks. We also ensure that they are contained in the ball of radius λmax
      if sup(δ) < 1
        for j=1:N0
         global D, nb_pos, min_pos, δ, λmax, index_pos, λj
            if j ∉ index_pos
            μj = 𝒟[j,j]
                if abs(μj) >= (interval(1)-δ)/λmax
                    λj = interval(1)/μj
                    if inf(real(λj)+δ*abs(λj)) > 0
                        println("Could not conclude about the sign of the eigenvalue $j")
                        break 
                        brk = brk
                    end

                    if sup(abs(λj)+δ*abs(λj)) > sup(λmax)
                        println("The disk $j crossed the boundary of the ball of radius λmax")
                        break 
                        brk = brk
                    end
                end
            end
        end
    end

    println("There are exactly $nb_pos eigenvalues with positive real part")

    
    for j ∈ index_pos 
        global  λj, δ, Rj, λmax
        λj = 1/𝒟[j,j]
        Rj = δ*abs(λj)
        λmax = abs(λj) + Rj #### we recompute λmax to get a better enclosure of the eigenvalue of interest
        δ = λmax*opnorm(R,1) + β0 + β1*λmax + β2*λmax^2 + Z1P + Z2P*r0
        Rj = sup(δ*abs(λj))
        λj = mid.(λj)
        println("Obtained an enclosure for an eigenvalue around $λj with radius $Rj")
    end

