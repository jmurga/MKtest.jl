using Comonicon

"
	julia --threads N abcmk_cli.jl rates 10000 661 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000 --gH 200,1000 --gL 1,10 --gam_flanking -1000,-500 --gam_dfe -1000,-100 --shape 0.184 --alpha 0.1,0.9 --iterations 1000000 --output rates.jld2

Function to solve analytical fixation rates and the expected SFS. The function will create N random models from prior values. Use the arguments to defined the input range for each parameter.

If rho and/or theta are not set, default values will be used (0.001).

To parallelize the estimations please be sure you start up Julia using --threads/-t option and set the number of cores.

The function returns a HDF5 file containing models solved and rates. The rates will be used to compute summary statistics required at ABC.

Please check the documentation to get more info about models parameters or detailed arguments description https://jmurga.github.io/MKtest.jl/dev/cli/ to check model
"
@lazyload using MKtest @cast function rates(N::Int64,
		samples::Int64;
		alpha::String="0.1,0.9",
		gam_neg::String="-1000,-200",
		gam_flanking::String="-1000,-500",
		gL::String="5,10",
		gH::String="400,1000",
		dac::String="1,2,4,5,10,20,50,100,200,400,500,661,925,1000",
		shape::Float64=0.184,
		rho::Float64=0.001,
		theta::Float64=0.001,
		iterations::Int64=100000,
		cutoff::String="0.0,1.0",
		output::String="rates.jld2")
    
	alpha, cutoff = map(x-> parse.(Int64,x),split.([alpha, cutoff],","))
    gam_neg, gam_flanking, gL, gH, dac = map(x-> parse.(Int64,x),split.([gam_neg,gam_flanking,gL,gH,dac],","))


    adap = MKtest.parameters(N=N,n=samples,dac=dac,cutoff=cutoff)

    println(adap)
    # MKtest.rates(adap,gam_neg=gam_neg,gam_flanking=gam_flanking,gL=gL,gH=gH,alpha=alpha,iterations=iterations,output=output)
end

@main