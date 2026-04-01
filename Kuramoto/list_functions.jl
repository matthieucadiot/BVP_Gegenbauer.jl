
function Deriv(k,N)
    R = interval.(zeros(N+1,N+1))
    if k==0
        R[1,2] = interval(2) 
        for i=1:N-1
            R[i+1,i+2] = interval(i+1)
        end
        return R   
    else
        R[1,2] = interval(4*k)
        for i=1:N-1
            R[i+1,i+2] = interval(2*k)
        end
        return R   
    end    
end 



function Int_k(k,N)
    R = interval.(zeros(N+1,N+1))
    if k==0
        R[1,2] = interval(0.5) 
        for i=1:N-1
            R[i+1,i+2] = interval(1)/interval(i+1)
        end
        return R   
    else
        R[1,2] = interval(1)/interval(4*k)
        for i=1:N-1
            R[i+1,i+2] = interval(1)/interval(2*k)
        end
        return R   
    end    
end 



function conversion(k,N,pres)
    setprecision(pres)
    R = interval.(big.(zeros(N+1,N+1)))
    if k==0
        R[1,1] = interval(big(1)) ; R[1,3] = -interval(big(1))
        for i=1:N-2
            R[i+1,i+1] = interval(big(0.5))
            R[i+1,i+3] = -interval(big(0.5))
        end
        for i=N-1:N 
            R[i+1,i+1] = interval(big(0.5))
        end
        return R   
    else
        R[1,1] = interval(1) ; R[1,3] = -interval(big(2))*interval(big(k))/(interval(big(2))+interval(big(k)))
        for i=1:N-2
            R[i+1,i+1] = interval(big(k))/(interval(big(k))+interval(big(i)))
            R[i+1,i+3] = -interval(big(k))/(interval(big(k))+interval(big(i))+interval(big(2)))
        end
        for i=N-1:N 
            R[i+1,i+1] = interval(big(k))/(interval(big(k))+interval(big(i)))
        end 
        return R   
    end    
end


# function PlotCoeffs1D(U0)
#     dx = 2/201
#     x = collect(-1:dx:1) 

#      m=length(x)
    
#      U = zeros(m)
#  for b₁ = 1:m
#         U[b₁] = U0(x[b₁])
#  end
    
#     mat"plot($x, $U)"
# end

using MATLAB
using MATLAB

function SaveMatlabFig(U0, filename="my_plot.fig")
    dx = 2/201
    x = collect(-1:dx:1) 
    U = U0.(x)
    
    # Get the full absolute path so MATLAB knows exactly where to write
    full_path = joinpath(pwd(), filename)
    
    put_variable(:x_data, x)
    put_variable(:y_data, U)
    put_variable(:file_path, full_path)
    
    println("Generating MATLAB figure file...")

   mat"""
    % 1. Create it invisible so Java doesn't crash Julia
    f = figure('Visible', 'off'); 
    plot(x_data, y_data, 'LineWidth', 1.5);

    
    % 2. CRITICAL: Set the property to 'on' BEFORE saving
    % This doesn't try to render it now, it just writes 'Visible=on' into the file
    set(f, 'Visible', 'on');
    
    % 3. Save the file
    savefig(f, file_path);
    close(f);
    """
    
    if isfile(full_path)
        println("Success! Figure saved to: $full_path")
        println("You can now open this file directly in your MATLAB R2021b app.")
    else
        println("Error: The .fig file was not created.")
    end
end

# Example run:
# SaveMatlabFig(x -> sin(5*x), "SimulationResults.fig")

# Example usage:
# PlotCoeffs1D(x -> sin(pi * x))

# function PlotCoeffs1D(U0)
#     # #U0 is a sequence in 2D
#     # #a,b,c,d are the endpoints of the interval [a,b] × [c,d]

#     # dx = 2/201
#     # x = -1:dx:1
    
#     #
#     # mat"
#     # h = plot($x,$x) ; drawnow;"

#     dx = 0.1
#     x = collect(-1:dx:1) 
    
#     mat"
#     figure;       % FIX 2: Explicitly open a figure container
#     h = plot($x, $x); 
#     shg;          % FIX 3: 'Show Graph' - forces window to front
#     drawnow;      % Process the drawing events
#     "
# end


function translation(N)
    R = interval.(zeros(N+1,N+1))
        for i=0:N-1
            R[i+2,i+1] = interval(1)
        end
        return R 
end



function odd2full(U::Vector{Float64},N)
    V = zeros(N+1)
    for n=1:2:N 
        V[n+1] = U[Int((n+1)/2)]
    end 
    return V 
end 

function odd2full(U::Vector{BigFloat},N)
    V = zeros(N+1)
    for n=1:2:N 
        V[n+1] = U[Int((n+1)/2)]
    end 
    return V 
end 

function odd2full(U::Vector{Interval{Float64}},N)
    V = interval.(zeros(N+1))
    for n=1:2:N 
        V[n+1] = U[Int((n+1)/2)]
    end 
    return V 
end 

function odd2full(U::Vector{Interval{BigFloat}},N)
    V = interval.(big.(zeros(N+1)))
    for n=1:2:N 
        V[n+1] = U[Int((n+1)/2)]
    end 
    return V 
end 


function even2full(U::Vector{Interval{Float64}},N)
    V = interval.(zeros(N+1))
    for n=0:2:N 
        V[n+1] = U[Int((n+2)/2)]
    end 
    return V 
end 


function even2full(U::Vector{Interval{BigFloat}},N)
    V = interval.(big.(zeros(N+1)))
    for n=0:2:N 
        V[n+1] = U[Int((n+2)/2)]
    end 
    return V 
end 

function even2full(U::Vector{Float64},N)
    V = zeros(N+1)
    for n=0:2:N 
        V[n+1] = U[Int((n+2)/2)]
    end 
    return V 
end 

function even2full(U::Vector{BigFloat},N)
    V = zeros(N+1)
    for n=0:2:N 
        V[n+1] = U[Int((n+2)/2)]
    end 
    return V 
end 






function Fcheb(U,α,N)

    V = Sequence(Chebyshev(N),vec(odd2full(U,N)))

    D4 =  mid.(translation(N))^4*mid.(Deriv(3,N))*mid.(Deriv(2,N))*mid.(Deriv(1,N))*mid.(Deriv(0,N)) 
    D2 =  mid.(Deriv(1,N))*mid.(Deriv(0,N))
    D1 =  mid.(Deriv(0,N))

    conv2 = mid.(translation(N))^4*mid.(conversion(3,N,pres))*mid.(conversion(2,N,pres))
    conv3 = mid.(translation(N))^4*mid.(conversion(3,N,pres))*mid.(conversion(2,N,pres))*mid.(conversion(1,N,pres))

    L = α*D4 + conv2*D2

    F = L*coefficients(V) + 0.5*conv3*D1*coefficients(project(V*V,Chebyshev(N))) 
    # F[1] = sum(vec(mid.(ξ1)).*V)  
    F[2] = sum(vec(mid.(ξ1[1:N+1])).*V)  
    # F[3] = sum(vec(mid.(ξ3)).*V)  
    F[4] = sum(vec(mid.(ξ3[1:N+1])).*V)  

    return F[2:2:N]
end 


function DFcheb(U,α,N)
    V = Sequence(Chebyshev(N),vec(odd2full(U,N)))

    D4 =  mid.(translation(N))^4*mid.(Deriv(3,N))*mid.(Deriv(2,N))*mid.(Deriv(1,N))*mid.(Deriv(0,N)) 
    D2 =  mid.(Deriv(1,N))*mid.(Deriv(0,N))
    D1 =  mid.(Deriv(0,N))

    conv2 = mid.(translation(N))^4*mid.(conversion(3,N,pres))*mid.(conversion(2,N,pres))
    conv3 = mid.(translation(N))^4*mid.(conversion(3,N,pres))*mid.(conversion(2,N,pres))*mid.(conversion(1,N,pres))

    L = α*D4 + conv2*D2

    DF = L + conv3*D1*coefficients(project(Multiplication(V),Chebyshev(N),Chebyshev(N))) 
    # DF[1,:] = mid.(ξ1) 
    DF[2,:] = mid.(ξ1[1:N+1]) 
    # DF[3,:] = mid.(ξ3) 
    DF[4,:] = mid.(ξ3[1:N+1]) 

    return DF[2:2:N,2:2:N]
end 



function eval_boundary1(N,k,pres)

    setprecision(pres)
    if k==0
        ξ = interval.(big.(2*ones(N+1)))
        ξ[1] = interval(big(1))
        return ξ
    else

        n = interval.(big.(1:N+1)) 
        ξ = ( n .+ interval(big(1))).*(n .+ interval(big(k)))
        for j=2:2k-1
            ξ = ξ.*(n .+ interval(big(j)))
        end 
        ξ = [interval(big.(zeros(k+1)));interval(big(2^k))*interval(big(factorial(k-1)))*ξ[1:N-k]/interval(big(factorial(2*k-1)))]
        ξ[k+1] = interval(big(2^k))*interval(big(factorial(k)))
        return ξ
    end
end 

function eval_boundarym1(N,k,pres)

    setprecision(pres)
    if k==0
        ξ = interval.(big.(2*ones(N+1)))
        ξ[1] = interval(big(1))
        return ξ.*(interval.(big.( (-1).^(0:N) )))
    else
        n = interval.(big.(1:N+1)) 
        ξ = ( n .+ interval(big(1))).*(n .+ interval(big(k)))
        for j=2:2k-1
            ξ = ξ.*(n .+ interval(big(j)))
        end 
        ξ = ξ.*(interval.(big.( (-1).^(0:N) )))
        ξ = [interval.(big.(zeros(k+1)));interval(big(2^k))*interval(big(factorial(k-1)))*ξ[1:N-k]/interval(big(factorial(2*k-1)))]
        ξ[k+1] = interval(big(2^k))*interval(big(factorial(k)))
        return ξ
    end
end 





function Newton_cheb(U,α,N)

    nf = norm(mid.(Fcheb(U,α,N)),1)
    k=0
    display(nf)
    while (nf>ϵ)&&(k<kmax)
        U = U - mid.(DFcheb(U,α,N))\mid.(Fcheb(U,α,N))
        k = k+1
        nf = norm(mid.(Fcheb(U,α,N)),1)
        display(nf)
    end

    return U, nf

end


function operator_E(k,θ,α,m,N,pres)
    E = interval.(big.(zeros(k,N+1)))

    for j = 1:α 
        for i=0:2*m-1
            E[j,:] = E[j,:] + θ[j,i+1]*eval_boundary1(N,i,pres)
        end
    end 

    for j = α+1:k
        for i=0:2*m-1
            E[j,:] = E[j,:] + θ[j,i+1]*eval_boundarym1(N,i,pres)
        end
    end 

    return E
end 



function operator_E_odd(k,θ,α,m,N,pres)
    E = operator_E(k,θ,α,m,N,pres)
    return E[1:α,2:2:end]
end 

function operator_E_even(k,θ,α,m,N,pres)
    E = operator_E(k,θ,α,m,N,pres)
    return E[1:α,1:2:end]
end 


function operator_B_inv(k,θ,α,m,N,pres)
    E = operator_E(k,θ,α,m,N,pres)
    B = Int_k(0,N)
    for i=1:k-1
        B = Int_k(i,N)*B
    end
    B = translation(N)^k*B

    B11 = inv(E[1:k,1:k])
    B[1:k,:] = -B11*E .* diag(B)'
    
    return B
end


function operator_B_inv_odd(k,θ,α,m,N,pres)
    
    B = Int_k(0,N)
    for i=1:k-1
        B = Int_k(i,N)*B
    end
    
    B = translation(N)^k*B
    b = diag(interval.(big.(mid.(B))))

    if mod(k,2)==0
        E = operator_E_odd(k,θ,α,m,N,pres)
        B = B[2:2:end,2:2:end]
        b = b[2:2:end]
        B11 = inv(E[1:α,1:α])
        B[1:α,:] = -(B11*E) .* b'
        return B
    else
        E = operator_E_even(k,θ,α,m,N,pres)
        B = B[1:2:end,1:2:end]
        b = b[1:2:end]
        B11 = inv(E[1:α,1:α])
        B[1:α,:] = -(B11*E) .* b'
        return B
    end
    
end


function operator_K(k,θ,α,m,N,pres)
    
    K = conversion(0,N+k,pres)
    for i=1:k-1
        K = conversion(i,N+k,pres)*K
    end
    K = interval.(big.(mid.(translation(N+k)^k)))*K
    K = operator_B_inv(k,θ,α,m,N+k,pres)*K[:,1:N+1]
    
    return K
end


function operator_K_odd(k,θ,α,m,N,pres)
    
    K = conversion(0,N+k,pres)
    for i=1:k-1
        K = conversion(i,N+k,pres)*K
    end
    K = interval.(big.(mid.(translation(N+k)^k)))*K

    if mod(k,2)==0
        K = operator_B_inv_odd(k,θ,α,m,N+k,pres)*K[2:2:end,2:2:N+1]
        return K
    else
        K = operator_B_inv_odd(k,θ,α,m,N+k,pres)*K[1:2:end,2:2:N+1]
        return K
    end
    
end





function Newton_KS(U,ν,Linv,K1,K2,N)

    nf = norm(mid.(F_KS(U,ν,Linv,K1,K2,N)),1)
    k=0
    display(nf)
    while (nf>ϵ)&&(k<kmax)
        U = U - mid.(DF_KS(U,ν,Linv,K1,K2,N))\mid.(F_KS(U,ν,Linv,K1,K2,N))
        k = k+1
        nf = norm(mid.(F_KS(U,ν,Linv,K1,K2,N)),1)
        display(nf)
    end

    return U, nf

end




function F_KS(U,ν::Float64,Linv,K1,K2,N)

    U_DU = Sequence(Chebyshev(N+3),vec(even2full(K1*U,N+3)))*Sequence(Chebyshev(N+4),vec(odd2full(Linv*U,N+4)))
    return U  + ν*(K2*U)[1:N÷2] + ν*U_DU[1:2:end][1:N÷2] 
end




function DF_KS(U,ν::Float64,Linv,K1,K2,N)

    V1 = Sequence(Chebyshev(N+3),vec(even2full(K1*U,N+3)))
    V2 = Sequence(Chebyshev(N+4),vec(odd2full(Linv*U,N+4)))
    MV1 = project(Multiplication(V1),Chebyshev(N+4),Chebyshev(2N+7))
    MV2 = project(Multiplication(V2),Chebyshev(N+3),Chebyshev(2N+7))

    MV1 = MV1[1:2:end,1:2:end] ### size N + 4  by N/2 + 2
    MV2 = MV2[1:2:end,0:2:end] ### size N + 4  by N/2 + 2

    ##### need to padd K2 with a bunch of zeros at the bottom to make the sizes work
    K2 = [K2;zeros(N+3 - (N÷2 +1),N÷2)] ### size N + 4 by N/2
    DF =  ν*K2 + ν*MV2*K1 + ν*MV1*Linv
    return DF[1:N÷2,1:N÷2] + I(N÷2)  #### add the identity
end





function F_KS(U,ν::Interval{Float64},Linv,K1,K2,N)

    U_DU = Sequence(Chebyshev(N+3),vec(even2full(K1*U,N+3)))*Sequence(Chebyshev(N+4),vec(odd2full(Linv*U,N+4)))
    N0 = length(U_DU[1:2:end])
    N1 = length(U)

    return [U ; interval.(zeros(N0 - N1))]  + ν*[K2*U; interval.(zeros(N0 - N1 - 2))] + ν*U_DU[1:2:end]
end


function DF_KS(U,ν::Interval{Float64},Linv,K1,K2,N)

    V1 = Sequence(Chebyshev(N+3),vec(even2full(K1*U,N+3)))
    V2 = Sequence(Chebyshev(N+4),vec(odd2full(Linv*U,N+4)))
    MV1 = project(Multiplication(V1),Chebyshev(N+4),Chebyshev(2N+7))
    MV2 = project(Multiplication(V2),Chebyshev(N+3),Chebyshev(2N+7))

    MV1 = MV1[1:2:end,1:2:end] ### size N + 4  by N/2 + 2
    MV2 = MV2[1:2:end,0:2:end] ### size N + 4  by N/2 + 2
    
    ##### need to padd K2 with a bunch of zeros at the bottom to make the sizes work
    K2 = [K2;interval.(zeros(N+3 - (N÷2 +1),N÷2))] ### size N + 4 by N/2
    DF =  ν*K2 + ν*(MV2*K1) + ν*(MV1*Linv)
    DF[1:N÷2,:] = DF[1:N÷2,:] + interval.(I(N÷2))  #### add the identity
    return DF
end



function tail_estimates(N,pres)

    θ0 = interval.(zeros(6,6))
    α0 = 3
    θ0[1,2] = interval(1)
    θ0[2,4] = interval(1)
    θ0[3,6] = interval(1)
    θ0[4,2] = interval(1)
    θ0[5,4] = interval(1)
    θ0[6,6] = interval(1)

    θ1 = interval.(zeros(7,8))
    α1 = 4
    θ1[1,1] = interval(1)
    θ1[2,3] = interval(1)
    θ1[3,5] = interval(1)
    θ1[4,7] = interval(1)
    θ1[5,1] = interval(1)
    θ1[6,3] = interval(1)
    θ1[7,5] = interval(1)

    K0 = operator_K_odd(6,θ0,α0,3,N,pres)
    K1 = operator_K_odd(7,θ1,α1,4,N,pres)

    return K0, K1
end 