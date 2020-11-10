# This file contains the code to compute the Adaptive center-surround stage of the model.
# Vectorized version of the Adaptive_CS_output stage with Interpolated Spatiotemporal input signal
# I added the interpolation function Spline!
using LsqFit
using DifferentialEquations
using Dierckx
using Sundials
using MAT
using JLD
using Plots
using IterableTables, DataFrames
pyplot()


# creating Gaussian function
function Gaussian(G, Amp, var, halfsize, fullsize)
    for i = 1:fullsize
        G[i] = Amp*exp(-(i - halfsize)*(i - halfsize)/var^2)
    end
    return G
end

# Center non-linearity
function SigmoidC_unit(x)
    b, center, amp = [75.251, 0.261, 7.] #23
    NonL_C = amp/(1 + exp(-b*(x - center  ))) #pass all vector x through the sigmoid
    return NonL_C
end

# Surround non-linearity
function SigmoidS_unit(x)
    b, center, amp = [80.7925, 0.2685, 6.5] # In case of emergency use 0.235
    NonL_S = amp/(1 + exp(-b*(x - center  ))) #pass all vector x through the sigmoid
    return NonL_S
end

# Sigmoid Non-linearity for the center
function SigmoidC(x)
    b, center, amp = [75.251, 0.261, 7.] #23
    NonL_C = amp./(1 .+ exp.(-b*(x .- center  ))) #pass all vector x through the sigmoid
    NonL_C[NonL_C .< 0.3060265654022103] .= 0.23 # change the values of neurons with x < 0.22
    return NonL_C
end

# Sigmoid Non-linearity for the surround
function SigmoidS(x)
    b, center, amp = [80.7925, 0.2685, 6.5] # In case of emergency use 0.235
    NonL_S = amp./(1 .+ exp.(-b*(x .- center  ))) #pass all vector x through the sigmoid
    NonL_S[NonL_S .< 0.3493847156007452] .= 0.253 # change the values of neurons with x < 0.23
    return NonL_S
end

# Creates the input in space...
function Inp_norm(I, S, contrast, base, middle)

    s = 392 # 392
    I[1:(middle - s)] .= base
    I[(middle - s):(middle + s)] .= contrast # 5.00 deg
    I[(middle + s):end] .= base
return I

end

# Spatiotemporal Input for channel 1
function SgnCh1(spl, t )

    s = collect(1:7198) #7198
    tt = ones(7198).*t
    I = spl(s, tt)
    return I
end

# Spatiotemporal Input for channel 2
function SgnCh2(spl, t )

    s = collect(1:7198) #7198
    tt = ones(7198).*t
    I = spl(s, tt)
    return I

end


# Main program starts here
function main()

    # Kernels
    # Determining the Gaussian RFs
    fullsize = 3599
    halfsize = Int((fullsize+1)/2)
    Gc = Array{Float64}(undef, fullsize)
    Gs = Array{Float64}(undef, fullsize)
    C = 67 #35
    S = 300 #94
    RFc = Gaussian(Gc, 2.0/1, C, halfsize, fullsize)
    RFs = Gaussian(Gs, 1.0/1, S, halfsize, fullsize)
    Dog = RFc - RFs
    display(plot(Dog, color=:black, linewidth=3, grid=false, legend=false, border=:none))
    #savefig("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Plots/DoG")

    fullsize = 3600
    # Initial Conditions of the system
    L0 = 60.5
    I = L0*ones(fullsize) # note: the signal starts at position i = 7
    sumGc = sum(RFc)
    sumGs = sum(RFs)
    # Parameters of the equation
    B, A, By, Dy = 1., 1., 100.0, 50.0
    # initial conditions for neuron layer x
    I0 = L0*10 # first value after zero padding
    Ij = L0*35.837#
    X0 = B*I0/(0.5 + I0 + Ij)
    println(X0)
    #println(sumGs)

    # initial conditions for neuron layer y
    # for the initial condition use the value assigned to the output of the nonlinearity for the baseline X0 = 0.22
    Y0 = (By * sumGc * 0.23- Dy * sumGs * 0.253)/(A + sumGc * 0.23 + sumGs * 0.253)

    println(Y0)

    # Defining the function
    function f(du, u, p, t)
        #x, y = u
        A, B, D, RFc, RFs, Inp, numcell, filter = p

        println(t)

        I = SgnCh1(Inp, t) # obtain the input signal at time t

        ###### Middle of the x cell array ######
        for i = 1:numcell
            #println(i)
            esum, isum = 0.0, 0.0
            esum = sum(SigmoidC(I[i:i+filter-1]).*RFc)
            isum = sum(SigmoidS(I[i:i+filter-1]).*RFs)

            du[i] = -A*u[i] + (B - u[i])*esum - (D + u[i])*isum
        end

    end # end of function f

    # ODE
    y0 = Array{Float64}(undef, (fullsize))
    # fill array with initial values
    for i =1:fullsize
        y0[i] = Y0
    end


    # Input Stimuli S
    # load the input from the Normalization Stage
    ch = "ch1"
    if ch == "ch1"
        In = load("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/contrast13_390.jld")
        In_norm = In["xi_ss"]
        t = load("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/contrast13_time_390.jld")
        time_1 = t["xi_ss"]
    else
        In = load("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/contrast46_MStask_ch2.jld")
        In_norm = In["xi_ss"]
        t = load("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/contrast46_time_MStask_ch2.jld")
        time_1 = t["xi_ss"]
    end

    input_r = size(In_norm)[1]
    input_c = size(In_norm)[2]

    # Zero padding
    filter = size(RFc)[1]  # Size of the Kernel
    pad = (filter - 1) ÷ 2 # Integer division. Number of zeros to be added at both edges.
    input_padded = zeros(input_r+(2*pad), input_c)
    # fill the new padded input
    for i in 1:input_r
        for j in 1:input_c
            input_padded[i+pad, j] = In_norm[i, j]
        end
    end
    Input_1 = input_padded
    # Interpolation function
    s = collect(1:size(Input_1)[1]) # spatial array
    spl = Spline2D(s, time_1, Input_1) # creates and interpolation object

    display(plot(time_1, Input_1[3600, :]))
    input_r = size(Input_1)[1]
    result = zeros(input_r-filter+1) # Size of the output signal
    result_r = size(result)[1]
    display(input_r)
    p = [A, By, Dy, RFc, RFs, spl, result_r, filter]
    tspan = (0.0, 16.0)
    prob = ODEProblem(f,y0,tspan,p)
    sol = solve(prob, alg_hints=[:stiff])

    return sol # main return
end


out = @time main() # run the main function
times = [t for t in out.t] # total number of time steps
fullsize = 3600
xi = Array{Float64}(undef, (fullsize, length(times)))

for i = 1:fullsize
    xi[i, :] = [u[i] for u in out.u] # for each cell i store all temporal stamps (:)
end



save("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics13_390.jld", "xi_ss", xi)
save("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics13_390_times_.jld", "xi_ss", times)
# save the file into .mat format
ch = "ch1"
if ch == "ch1"
    filem = matopen("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics13_390.mat", "w")
    write(filem, "MDynamics13_390", xi)
    close(filem)

    file_time = matopen("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics13_390_times_.mat", "w")
    write(file_time, "MDynamics13_390_times_", times)
    close(file_time)
    display(plot(times, xi[1800, :]))
    display(plot(xi[:, floor(Int8, length(times)/2)]))
else
    filem = matopen("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Modeling paper/Data/MStask/Trial2/MStask46_ch2.mat", "w")
    write(filem, "MStask46_ch2", xi)
    close(filem)

    file_time = matopen("/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Modeling paper/Data/MStask/Trial2/MStask46_times_ch2.mat", "w")
    write(file_time, "MStask46_times_ch2", times)
    close(file_time)
    display(plot(times, xi[1800, :]))
    display(plot(xi[:, floor(Int8, length(times)/2)]))
end
