#HAAR STATE
using Random
using LinearAlgebra
using Statistics
using Dates
using CSV, DataFrames
using DelimitedFiles
using DataFrames, GLM

function f(x)
        if x == 1
            return [0,1]
        end
        if x == 0 
            return [1,0]
        end
end

# Turns the binary list to the base ten number number
function tob10(vecs)
    join(vecs)
    parse(Int, join(vecs); base = 2)
end

# Splits the Vector
function splitleft(state, dimstate)
    statenew = zeros(dimstate)
    for i in 1:dimstate
        statenew[i] = Int64(state[i])
    end
    return round.(Int,statenew)
end

function splitright(state,dimstate)
    statenew = zeros(length(state)-dimstate)
    for i in dimstate+1:length(state)
        statenew[i-dimstate] = Int64(state[i])
    end
    return round.(Int,statenew)
end

px = [0 1;1 0]
py = [0 -1im; 1im 0]
pz = [1 0; 0 -1]

function ranm()
    alp = rand(Float64)*2*pi
    psi = rand(Float64)*2*pi
    chi = rand(Float64)*2*pi
    eps = rand(Float64)
    phi = asin(sqrt(eps))
    return exp(1im*alp)*[exp(1im*psi)*cos(phi) exp(1im*chi)*sin(phi); -1*exp(-1im*chi)*sin(phi) exp(-1im*psi)*cos(phi)]
end

function rangate()
    phig = rand(Float64)*2*pi
    J = rand(Float64)*2*pi
    up = ranm()
    um = ranm()
    vp = ranm()
    vm = ranm()
    return exp(1im*phig)*kron(up,um)*exp(-1im*(pi/4)*(kron(px,px)+kron(py,py))-1im*J*kron(pz,pz))*kron(vp,vm)
end

function binaryvec(n)
    ket0 = string(n, base=2, pad = nspin)
    for k in range(1,length(ket0))      
        kettemp[k] = parse(Int64,ket0[k:k]);
    end
    return kettemp
end

function reyni(state,nth,cut)
    # Generate the basis kets, and the Matrix A for SVD
    n = nth
    diml = cut
    dimr = nspin-cut
    A = zeros(ComplexF64,2^dimr,2^diml)
    srn = 0 
    for j in range(0,dim-1)
        l = tob10(splitright(ket[j+1],cut))
        m = tob10(splitleft(ket[j+1],cut))
        A[l+1,m+1] = state[j+1] 
    end
    # Compute the Von-Neumann Entropy and the nth Reyni Entropy
    evls = svdvals(A)
        for i in 1:length(evls)
            srn += ((norm(evls[i]))^2n)
        end
        srni = (1/(1-n))*log2(srn)
        return(srni)
end


function gate(statein,i)
    L = nspin
    U = rangate()
    psi_temp = zeros(ComplexF64, 4)
    if i!=nspin
        for nL=0:2^(L-1-i)-1, nR=0:2^(i-1)-1
            ind1 = nR + nL*2^(i+1) + 1
            ind2 = nR + nL*2^(i+1) + 2^(i-1) + 1
            ind3 = nR + nL*2^(i+1) + 2^(i) + 1
            ind4 = nR + nL*2^(i+1) + 2^(i) + 2^(i-1) + 1
            psi_temp[1] = statein[ind1]
            psi_temp[2] = statein[ind2]
            psi_temp[3] = statein[ind3]
            psi_temp[4] = statein[ind4]
            psi_temp = U*psi_temp
            statein[ind1] = psi_temp[1]
            statein[ind2] = psi_temp[2]
            statein[ind3] = psi_temp[3]
            statein[ind4] = psi_temp[4]
        end
    end
    return(statein)
end

function even_projmsm(statein)
# Projective Measurements
    L = nspin
    for i=1:L
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR+2^(i)*nL + 2^(i-1)+1
                pup += statein[n]*conj(statein[n])
            end
            p2 = rand()
        # Measure up 
            if p2<Real(pup)
                FET[t:end] .+= -log(pup)
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL +1
                    statein[n] = 0 
                end
        # Measure down 
            else 
                FET[t:end] .+= -log(1-pup)
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL + 2^(i-1)+1
                    statein[n] = 0
                end
            end
        end
        statein = statein/norm(statein)
    end
    return(statein)
end


function odd_projmsm(statein)
# Projective Measurements
    L = nspin
    for i=2:L-1
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR+2^(i)*nL + 2^(i-1)+1
                pup += statein[n]*conj(statein[n])
            end
            p2 = rand()
        # Measure up 
            if p2<Real(pup)
                FET[t:end] .+= -log(pup)
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL +1
                    statein[n] = 0 
                end
        # Measure down 
            else 
                FET[t:end] .+= -log(1-pup)
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL + 2^(i-1)+1
                    statein[n] = 0
                end
            end
        end
        statein = statein/norm(statein)
    end
    return(statein)
end

function runy(nspins, iteration, probability, filename) 
    global state
    #Generate Initial State
    prob = probability
    global prob
    nspin = nspins
    dimen = 2^nspin
    global nspin
    global dimen
    ket = [zeros(Int64,nspin)];
    global ket
    kettemp = zeros(Int64,nspin)
    global kettemp
    # Initialization conditions for Gates
    timestep = Int64(2*nspin + 14*nspin)
    FET = zeros(2*timestep)
    global FET
    global timestep
    iterations = iteration
    vals1 = zeros(timestep)
    errors = zeros(iteration)
    x = zeros(2*timestep)
    cut = Int64(nspin/2)
    fout = zeros(2*timestep)
    t = 0
    global t
    state = zeros(ComplexF64,dimen)
    global state
    for j in 1:dimen-1
        temp = Vector(binaryvec(j))
        push!(ket,temp)
    end
    for k=1:iterations
        t = 0
        FET .= 0 
        # Generate a new Haar random state for each trajectory
        state = zeros(Complex{Float64}, dimen)
        for i=1:dimen
            state[i] = rand(Normal(0,1))[1] + rand(Normal(0,1))[1]im
        end
        state = state/norm(state)
        for i=1:timestep
            for j=1:2:nspin
               state = gate(state,Int64(j))
            end
            t = t + 1
            x[t] = t
            state = even_projmsm(state)
            for j = 2:2:nspin-1
                state = gate(state,Int64(j))
            end
            t = t + 1
            x[t] = t
            state = odd_projmsm(state)
        end
        fout .+= (FET)/(iterations)
    end
    writedlm(string(filename,"_FET"),fout,", ")
    return(fout)
end
runy(8,1,0.14,"test")
