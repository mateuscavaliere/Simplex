

function SimplexFase1(A::Array{Float64},b::Array{Float64},c::Array{Float64},maxIter::Int32,epslon = 10.0^-9)
    
    A = hcat(A,-1*ones(1,size(A)[1])')
    c = vcat(zeros(length(c)),-1)

    m = size(A)[1] + 1
    n = size(A)[2] - m

    vN = vcat([i for i in 1:n],size(A)[2])
    vB  = [i for i in (n+1):(n+m-1)]

    iter    = 1

    B = A[:,vB]
    N = A[:,vN]
    xB = B \ b

    if iter == 1    # devemos inserir a variavel artificial (w) na base
        j          = length(vN)    # indice da variavel artificial (w) inserida no problema
        k          = indmin(xB)     # indice da variavel basica com o menor valor
        varR_Base  = vB[k]   # coluna associada a variavel que deve ser retirada da base
        varR_nBase = vN[j] # coluna associada a variavel que deve ser retirada da nao base
        vB[k]      = varR_nBase  # nova Base
        vN[j]      = varR_Base  # nova nBase
    end

    for i = 2:maxIter

        println("Iteracao # ",i)

        B = A[:,vB]
        N = A[:,vN]
        xB = B \ b
        dB = -B \ N
        cB = c[vB]
        cN = c[vN]
        cR = (cN' + cB'*dB)'
        j  = indmax(cR)

        x      = zeros(n+m)
        x[vB]  = xB

        println("x     = ", x)
        println("Base  = ", vB)
        println("nBase = ", vN)
        println("")

        if all(cR .<= epslon)
            vB
            vN = vN[vN .!= size(A)[2]]
            status = 1
            # return (sort(vB), sort(vN), status)
            return (vB, vN, status)
        end

        if all(dB[:,j] .>= epslon)
            status = -2
            return (status);
        end

        r  = xB./(dB[:,j])
        k  = indmax(r[r .< 0])

        auxB  = vB[k]
        vB[k] = vN[j]
        vN[j] = auxB
    end
end

function Simplex(A::Array{Float64},b::Array{Float64},c::Array{Float64},maxIter::Int32,epslon = 10.0^-9)

    m = size(A)[1]
    n = size(A)[2] - m

    vB = [i for i in (n+1):(n+m)]
    vN = [i for i in 1:n]

    iter    = 1

    B = A[:,vB]
    N = A[:,vN]

    xB = B \ b

    if any(xB .< 0)
        vB,vN,flag = SimplexFase1(A,b,c,maxIter)
    end

    
    println("")
    println("##############################")
    println("------- Simplex Fase 2 -------")
    println("##############################")
    println("")

    for i = 1:maxIter

        println("Iteracao # ",i)

        B = A[:,vB]
        N = A[:,vN]
        xB = B \ b
        dB = -B \ N
        cB = c[vB]
        cN = c[vN]
        cR = (cN' + cB'*dB)'
        j  = indmax(cR)

        x      = zeros(n+m)
        x[vB]  = xB

        println("x     = ", x)
        println("Base  = ", vB)
        println("nBase = ", vN)
        println("")
        
        if all(cR .<= epslon)
            z      = c'*x;
            status = 1;
            println("Solucao obtida:")
            println("x      = ", x)
            println("z      = ", z)
            println("status = ", status)
            return (x,z,status);
        end
        
        if all(dB[:,j] .>= epslon)
            d      = zeros(n+m)
            d[vN]  = 1
            d[vB]  = dB[:,j]
            z      = Inf
            status = -1
            println("Problema ilimitado!")
            println("d      = ", d)
            println("z      = ", z)
            println("status = ", status)
            return (d,z,status);
        end
        
        r  = xB./(dB[:,j])
        k  = indmax(r[r .< 0])
        
        auxB  = vB[k]
        vB[k] = vN[j]
        vN[j] = auxB
    end
end

#--- Problema 1

A = float([2 1 1 0 ; 1 2 0 1]);
b = float([4;4]);
c = float([4;3;0;0]);
maxIter = Int32(1000);

x,z,status = Simplex(A,b,c,maxIter);

#---- Problema 2

A = float([0.5 -1 1 0; -4 1 0 1])
b = float([0.5 ; 1])
c = float([1 ; 1; 0; 0])
maxIter = Int32(1000);

x,z,status = Simplex(A,b,c,maxIter);


#---- Problema 3

A = float([2 1 1 0 0; 1 2 0 1 0; -1 -1 0 0 1])
b = float([4, 4, -1])
c = float([4 3 0 0 0])
maxIter = Int32(1000);

x,z,status = Simplex(A,b,c,maxIter);
