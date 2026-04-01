
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



function conversion(k,N)
    R = interval.(zeros(N+1,N+1))
    if k==0
        R[1,1] = interval(1) ; R[1,3] = -interval(1)
        for i=1:N-2
            R[i+1,i+1] = interval(0.5)
            R[i+1,i+3] = -interval(0.5)
        end
        for i=N-1:N 
            R[i+1,i+1] = interval(0.5)
        end
        return R   
    else
        R[1,1] = interval(1) ; R[1,3] = -interval(2)*interval(k)/(interval(2)+interval(k))
        for i=1:N-2
            R[i+1,i+1] = interval(k)/(interval(k)+interval(i))
            R[i+1,i+3] = -interval(k)/(interval(k)+interval(i)+interval(2))
        end
        for i=N-1:N 
            R[i+1,i+1] = interval(k)/(interval(k)+interval(i))
        end 
        return R   
    end    
end


function PlotCoeffs1D(U0)
    dx = 2/201
    x = collect(-1:dx:1) 

     m=length(x)
    
     U = zeros(m)
 for b₁ = 1:m
        U[b₁] = U0(x[b₁])
 end
    
    # We use 'visible', 'off' to tell MATLAB not to try opening a window
    # Then we save the plot as 'test_plot.png' in your current folder
    mat" 
    f = plot($x, $U);
    title('Test Plot');
    saveas(f, 'plot_KS.fig');
    close(f);
    "
end



function translation(N)
    R = interval.(zeros(N+1,N+1))
        for i=0:N-1
            R[i+2,i+1] = interval(1)
        end
        return R 
end



function Fcheb(U,α,N)

    L =  mid.(translation(N))^2*mid.(Deriv(1,N))*mid.(Deriv(0,N))

    conv2 = mid.(translation(N))^2*mid.(conversion(1,N))*mid.(conversion(0,N))

    F = L*coefficients(U) + α*conv2*coefficients(project(U*U,Chebyshev(N)) + 1) 
    # F[1] = sum(vec(mid.(ξ1)).*V)  
    F[1] = sum(vec(mid.(ξ1[1:N+1])).*U)  
    # F[3] = sum(vec(mid.(ξ3)).*V)  
    F[2] = sum(vec(mid.(ξ2[1:N+1])).*U)  

    return Sequence(Chebyshev(N),F)
end 



function weight_nu(ν,N)
    w = interval.(zeros(N+1))
    w[1] = interval(1)
    for i=1:N
        w[i+1] = interval(2)*ν^interval(i)
    end
    return w
end



function DFcheb(U,α,N)

    L =  mid.(translation(N))^2*mid.(Deriv(1,N))*mid.(Deriv(0,N))
    conv2 = mid.(translation(N))^2*mid.(conversion(1,N))*mid.(conversion(0,N))
    V = coefficients(project(Multiplication(U),Chebyshev(N),Chebyshev(N))) 

    DF = L + 2*α*conv2*V
    # DF[1,:] = mid.(ξ1) 
    DF[1,:] = mid.(ξ1[1:N+1]) 
    # DF[3,:] = mid.(ξ3) 
    DF[2,:] = mid.(ξ2[1:N+1]) 

    return LinearOperator(Chebyshev(N), Chebyshev(N), DF)
end 



function eval_boundary1(N,k)

    if k==0
        ξ = interval.(2*ones(N+1))
        ξ[1] = interval(1)
        return ξ
    else

        n = interval.(1:N+1) 
        ξ = ( n .+ interval(1)).*(n .+ interval(k))
        for j=2:2k-1
            ξ = ξ.*(n .+ interval(j))
        end 
        ξ = [interval(zeros(k+1));interval(2^k)*interval(factorial(k-1))*ξ[1:N-k]/interval(factorial(2*k-1))]
        ξ[k+1] = interval(2^k)*interval(factorial(k))
        return ξ
    end
end 

function eval_boundarym1(N,k)

    if k==0
        ξ = interval.(2*ones(N+1))
        ξ[1] = interval(1)
        return ξ.*(interval.( (-1).^(0:N) ))
    else
        n = interval.(1:N+1) 
        ξ = ( n .+ interval(1)).*(n .+ interval(k))
        for j=2:2k-1
            ξ = ξ.*(n .+ interval(j))
        end 
        ξ = ξ.*(interval.( (-1).^(0:N) ))
        ξ = [interval(zeros(k+1));interval(2^k)*interval(factorial(k-1))*ξ[1:N-k]/interval(factorial(2*k-1))]
        ξ[k+1] = interval(2^k)*interval(factorial(k))
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

    return U

end



function operator_Linv(N)
    #### the argument is of size N, the output is of size N+2 by N+2
    ξ1 = eval_boundary1(N+2,0) ### evaluation of the boundary at x=1
    ξ2 = eval_boundarym1(N+2,0) ### evaluation of the boundary at x=-1

    Linv = Int_k(1,N+2)*Int_k(0,N+2)
    Linv = translation(N+2)^2*Linv

    B11 = inv([ξ1[1:2]';ξ2[1:2]'])
    Linv[1:2,:] = -B11*[ξ1';ξ2'] .* diag(Linv)'
    Linv = Linv*translation(N+2)^2*conversion(1,N+2)*conversion(0,N+2) 

    return LinearOperator(Chebyshev(N), Chebyshev(N+2), Linv[1:N+3,1:N+1])
end


function F_KS(U,α,N)
    Linv = operator_Linv(N)
    e0 = Sequence(Chebyshev(N),interval.(zeros(N+1)))
    e0[0] = interval(1)
    return  U + α*(Linv*U)*(Linv*U) + α*e0
end




function DF_KS(U,α,N)
    Linv = operator_Linv(N)
    V = project(Multiplication(Linv*U),Chebyshev(N+2),Chebyshev(2*N+4))

    return  interval(2)*α*V*Linv + LinearOperator(Chebyshev(N+2), Chebyshev(N+2), interval.(I(N+3)))
end



function PlotCoeffs1D(U0)
    dx = 2/201
    x = collect(-1:dx:1) 

     m=length(x)
    
     U = zeros(m)
 for b₁ = 1:m
        U[b₁] = U0(x[b₁])
 end
    
    # We use 'visible', 'off' to tell MATLAB not to try opening a window
    # Then we save the plot as 'test_plot.png' in your current folder
    mat" 
    f = plot($x, $U);
    title('Test Plot');
    saveas(f, 'plot_toy_model.fig');
    close(f);
    "
end