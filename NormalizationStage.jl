# This computes the Normalization stage of the model

using LsqFit
using DifferentialEquations
using Sundials
using MAT
using JLD
using Plots
using IterableTables, DataFrames
pyplot()


# Spatiotemporal Input for channel 1
function SgnCh1(I, pad, t)
    # Inside this function you can manually select the stimulus size and
    # the desired contrast value for the input signal. We picked this method
    # because it is the fastest.
    s = 390 # Stimulus size
    middle = pad + 1*1800
    base = 60.5 # background luminance
    # Select the contrast stimulus here is list. You can always use this formula
    # Delta = (MC * 2 * base)/(1 - MC) + base ; where MC is the Michelsen Contrast and base is the background luminance
    # 61.81 --> 1.07%
    # 62.1 ---> 1.3%
    # 62.19 ---> 1.38%
    # 62.23 --> 1.4%
    # 62.34 --> 1.5%
    # 62.47 --> 1.6%
    # 62.72 --> 1.8%
    # 62.8435 --> 1.9%
    # 63.5 --> 2.43%
    # 63.6026 ---> 2.5%
    # 63.73 ---> 2.6%
    # 63.98602880658436 --> 2.8%
    # 64.88 ---> 3.5%
    # 65.56 --> 4%
    # 67.54232804232804 --> 5.5%
    # 69.02173913043478 --> 8%
    # 70.59866901330986 --> 9.09%
    # 73.94444444444444 --> 10%
    # 75.65393271070103 --> 11.13%
    # 88.94991108476586 --> 17.65%
    # 104.62820512820514 --> 24%
    # 163.57407407407408 --> 46%
    # 1450.0000000000007 --> 92%
    contrast = 62.19 # change this depending on the contrast value
    I[1:(middle - s)] .= base
    if t >= 0.0 && t <= 11.8
        I[(middle - s):(middle + s)] .= contrast
    else
        I[(middle - s):(middle + s)] .= base
    end
    I[(middle + s):(end)] .= base
    return I
end

# Spatiotemporal Input for channel 2
function SgnCh2(I, pad, t )

    s = 279 # 392
    middle = pad + 1*1800 #- 3*50
    startL = middle - 2*s -48#+ 46 # marks the start of the left signal
    endL = startL + 2*s  # left signal
    startR = middle + 0*s +48#- 46# marks the start of the right signal
    endR = startR + 2*s

    base = 60.5
    # Select the contrast stimulus here is list. You can always use this formula
    # Delta = (MC * 2 * base)/(1 - MC) + base ; where MC is the Michelsen Contrast and base is the background luminance
    # 61.81 --> 1.07%
    # 63.98602880658436 --> 2.8%
    # 67.54232804232804 --> 5.5%
    # 69.02173913043478 --> 8%
    # 70.59866901330986 --> 9.09%
    # 73.94444444444444 --> 10%
    # 75.65393271070103 --> 11.13%
    # 88.94991108476586 --> 17.65%
    # 104.62820512820514 --> 24%
    # 163.57407407407408 --> 46%
    # 1450.0000000000007 --> 92%

    contrast = 163.5740740    # change this depending on the contrast value
    I[1:startL] .= base
    I[endL:startR] .= base
    if t >= 0.5 && t <= 5.0
        I[startL:endL] .= contrast
        I[startR:endR] .= contrast
    else
        I[startL:endL] .= base
        I[startR:endR] .= base

    end
    I[endR:(end)] .= base # 7198 - 1799 = 5399
    return I

end

# Input for the motion discrimination task
function SgnCh3(I, pad, t )

    s = 350 # 392
    middle = pad + 1*1800 - 3*50
    startL = middle - 2*s + 46 # marks the start of the left signal
    endL = startL + 2*s  # left signal
    startR = middle + s -46# marks the start of the right signal
    endR = startR + 2*s

    base = 60.5
    # 65.54232804232804 --> 5.5%
    # 69.02173913043478 --> 8%
    # 73.94444444444444 --> 10%
    # 73.65393271070103 --> 11.13%
    # 80.94991108476586 --> 15.65%
    # 94.62820512820514 --> 22%
    #161.57407407407408 --> 46%
    # 1450.0000000000007 --> 92%
    contrast = 1450.0000000000007 # change this depending on the contrast value
    I[1:startL] .= base #1799
    if t >= 0.5 && t <= 2.0
        I[startL:endR] .= contrast
    else
        I[startL:endR] .= base

    end
    I[endR:(end)] .= base # 7198 - 1799 = 5399
    return I

end

# Main program starts here
function main()

    # Kernels
    n = 3599 # number of Kernel cells
    half_kernel = Int((n+1)/2)
    C_width = 2 # half of center's width
    S_width = 2 # half of surround's width
    Lc = half_kernel - C_width
    Rc = half_kernel + C_width
    Ls = Lc - S_width
    Rs = Rc + S_width
    Kc = zeros(n)
    Kc[Lc:Rc] .= 2.0
    Ks = zeros(n)
    Ks[Lc:Rc] .= 0.5
    Ks[Ls:Lc] .= 0.5
    Ks[Rc:Rs] .= 0.5
    Ks[1:Ls] .= 0.009 # inhibition for the rest of the values
    Ks[Rs:n] .= 0.009
    DoK = Kc - Ks

    fullsize = 1*3600
    # Initial Conditions of the system
    L0 = 60.5
    I = L0*ones(fullsize) # note: the signal starts at position i = 7

    # Parameters of the equation
    A, B, D = 0.5, 1.0, 0.0
    # initial conditions for neuron layer x
    I0 = L0*sum(Kc) # first value after zero padding
    Ij = L0*sum(Ks)


    X0 = B*I0/(A + I0 + Ij)
    println(X0)
    MC = 0.11 # Michelsen Contrast
    deltaI = (MC * 2 * L0)/(1 - MC) + L0
    println(deltaI)
    # Defining the function
    function f(du, u, p, t)
        A, B, RFc, RFs, Inp, numcell, filter, pad = p

        println(t)

        I = SgnCh1(Inp, pad, t) # call the input value at time t



        ###### Middle of the x cell array ######
        for i = 1:numcell
            #println(i)
            esum, isum = 0.0, 0.0
            esum = sum(I[i:i+filter-1].*RFc)
            isum = sum(I[i:i+filter-1].*RFs)
            du[i] = -A*u[i] + (B - u[i])*esum - u[i]*isum
        end

    end # end of function f

    # ODE
    x0 = Array{Float64}(undef, (fullsize))
    # fill array with initial values
    for i =1:fullsize
        x0[i] = X0
    end


    # Input Stimuli S
    Input_1 = zeros(1*3600) #arr_46
    input_r = size(Input_1)[1]

    # Zero padding
    filter = size(Kc)[1]  # Size of the Kernel
    pad = (filter - 1) ÷ 2 # Integer division. Number of zeros to be added at both edges.
    input_padded = zeros(input_r+(2*pad))

    Input_1 = input_padded
    input_r = size(Input_1)[1]
    result = zeros(input_r-2*pad) # Size of the output signal

    result_r = size(result)[1]
    display(input_r)
    p = [A, B, Kc, Ks, Input_1, result_r, filter, pad]
    tspan = (0.0,16.0)
    prob = ODEProblem(f,x0,tspan,p) # solver
    sol = solve(prob, alg_hints=[:stiff]) # use :stiff for stiff problems

    return sol # main return
end

out = @time main() # run the main function
times = [t for t in out.t] # total number of time steps
fullsize = 1*3600
xi = Array{Float64}(undef, (fullsize, length(times)))

for i = 1:fullsize
    xi[i, :] = [u[i] for u in out.u] # for each cell i store all temporal stamps (:)
end

ch = "ch1"
if ch == "ch1"
    save("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/contrast13_390.jld", "xi_ss",  xi)
    save("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/contrast13_time_390.jld", "xi_ss", times)
    println(size(xi))
    println(xi[1800, floor(Int8, length(times)/2)])
    println(xi[980, floor(Int8, length(times)/2)])
    display(plot(times, xi[1800, :]))
    display(plot(xi[:, floor(Int8, length(times)/2)]))
else
    save("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/contrast46_MStask_ch2.jld", "xi_ss",  xi)
    save("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/contrast46_time_MStask_ch2.jld", "xi_ss", times)
    println(size(xi))
    println(xi[1600, floor(Int8, length(times)/2)])
    println(xi[980, floor(Int8, length(times)/2)])
    display(plot(times, xi[1600, :]))
    display(plot(xi[:, floor(Int8, length(times)/2)]))
end

# save the file into .mat format
#filem = matopen("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/AdaptiveCS_outputStage/contrast92_MDtask_etapa1.mat", "w")
#write(filem, "Xi_MDtask_92_etapa1", xi)
#close(filem)

#file_time = matopen("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/AdaptiveCS_outputStage/contrast92_time_MDtask_etapa1.mat", "w")
#write(file_time, "t_MDtask_92_etapa1", times)
#close(file_time)
