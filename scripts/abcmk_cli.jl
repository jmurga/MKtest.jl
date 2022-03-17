using Fire, Distributed, MKtest

"
	julia abcmk_cli.jl rates --samples 661 --gam_neg -1000,-200 --gL 1,10 --gH 400,1000 --rho 0.001 --theta 0.001 --iterations 1000000 --output rates.jld2 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000 --nprocs 7

Function to solve fixation and polymorphic rates analitically. The function will create N random models from prior values. Use the arguments to defined the input range for each parameter.

If rho and/or theta are set to nothing, the function will input random values given the range 0.0005:0.0005:0.01. Otherwise you can fix the values.

If gL is set to nothing, the function will not account the role of the weakly selected alleles in the estimation.

If scheduler is set to threads the software will be run using multi-threading, not distributed computing. Please be sure you start up Julia using the same number of threads as the argument nprocs using the option julia -t nprocs

The function returns a HDF5 file containing models solved and rates. The rates will be used to compute summary statistics required at ABC.

Please check the documentation to get more info about models parameters or detailed arguments description https://jmurga.github.io/MKtest.jl/dev/cli/ to check model
"
@main function rates(;ne::Int64=1000,
		samples::Int64=661,
		alpha::String="0.1,0.9",
		gam_neg::String="-1000,-200", 
		gL::String="5,10", 
		gH::String="400,1000",
		dac::String="1,2,4,5,10,20,50,100,200,400,500,661,925,1000",
		shape::Float64=0.184,
		rho::String="nothing",
		theta::String="nothing",
		iterations::Int64=100000,
		output::String="/home/jmurga/rates.jld2",
		scheduler::String="local",
		nprocs::Int64=1
	)

	α          = parse.(Float64,split(alpha,","));
	tmp_neg    = parse.(Int,split(gam_neg,","));
	tmp_strong = parse.(Int,split(gH,","));
	dac        = parse.(Int,split(dac,","));

	if (gL == "nothing")
		tmp_weak = nothing
	else
		tmp_weak = parse.(Int,split(gL,","))
	end

	if (rho == "nothing")
		rho = nothing
	else
		rho = parse(Float64,rho)
	end

	if (theta == "nothing")
		theta = nothing
	else
		theta = parse(Float64,theta)
	end

	threads=false;
	if (scheduler == "slurm")
		@eval using ClusterManagers
		@eval addprocs_slurm($nprocs)
	elseif (scheduler == "htcondor")
		@eval using ClusterManagers
		@eval addprocs_htc($nprocs)
	elseif (scheduler == "local")
		@eval addprocs($nprocs)
	else
		threads=true;
	end
	
	@eval @everywhere using MKtest
	@eval adap = MKtest.parameters(N=$ne,n=$samples,dac=$dac,al=$shape);

	@eval MKtest.rates(
		param      = $adap, 
		α          = $α,
		gH         = $tmp_strong[1]:$tmp_strong[2],
		gL         = $tmp_weak[1]:$tmp_weak[2],
		gam_neg    = $tmp_neg[1]:$tmp_neg[2],
		iterations = $iterations,
		rho        = $rho,
		theta      = $theta,
		output     = $output,
		threads    = $threads
	);

	if scheduler != "threads"
		for i in workers()
			rmprocs(i)
		end
	end
end

"
	julia abcmk_cli.jl parse_data --analysis_folder analysis/ --gene_list analysis/dnaVipsList.txt

Function to parse polymorphic and divergence data from Uricchio et. al (2019) and Murga-Moreno et al (2019). Please input a path to create a new analysis folder. You can filter the dataset using a file containing a list of Ensembl IDs. 

The function returns files containing raw polymorphic and divergence data, parsed SFS and parsed divegence required to estimate summary statistics.	

Please check the documentation to get more info https://jmurga.github.io/MKtest.jl/dev/cli/
"
@main function data(;
	analysis_folder::String="<folder>",
	dataset::String="tgp",
	gene_list::String="false",
	bins::String="false")
	
	@eval using DataFrames, CSV

	mkpath(analysis_folder)

	dataset = lowercase(dataset)
	data    = analysis_folder * "/" * dataset * ".txt"
	
	download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/"* dataset * ".txt",data)

	# Check if bins or genelist are defined
	@eval if $gene_list != "false"
		@eval g_list = String.(Array(CSV.read($gene_list,DataFrame,header=false)))[:,2:end]
	else
		g_list = nothing
	end

	@eval if $bins != "false"
		@eval bins_size = parse(Int,$bins)
	else
		bins_size = nothing
	end

	# Parsing TGP data
	if dataset == "tgp"
		@eval α,sfs, divergence = MKtest.parse_sfs(sample_size=661,data=$data,gene_list=$g_list,bins=$bins_size)
	elseif occursin("zi",dataset)
		@eval α,sfs, divergence = MKtest.parseSfs(sample_size=154,data=$data,gene_list=$g_list,bins=$bins_size,isolines=true)
	elseif occursin("ral",dataset)
		@eval α,sfs, divergence = MKtest.parseSfs(sample_size=160,data=$data,gene_list=$g_list,bins=$bins_size,isolines=true)
	end

	# Writting data to folder
	@eval s_name = $analysis_folder .* "/sfs".*string.(1:length(sfs)).*".tsv"
	@eval d_name = $analysis_folder .* "/div".*string.(1:length(sfs)).*".tsv"

	@eval CSV.write.($s_name,DataFrame.($sfs,:auto),delim='\t',header=false)
	@eval CSV.write.($d_name,DataFrame.($divergence,:auto),delim='\t',header=false)

end

"
	julia abcmk_cli.jl summaries --analysis_folder analysis/ --rates analysis/rates.jld2 --samples 661 --dac 2,4,5,10,20,50,200,661,925 --summstat_size 1000000 --nprocs 7

Estimate summary statistics from analytical rates. You must provide a path containing the parsed SFS and divergence file.

The function returns files containing bootstrapped datasets (alphas.txt) and summary statistics (summstat.txt)

Check the documentation to get more info https://jmurga.github.io/MKtest.jl/dev/cli"
@main function summaries(;analysis_folder::String="<folder>",h5_file::String="rates.jld2",ne::Int64=1000, samples::Int64=661,dac::String="2,4,5,10,20,50,200,661,925",summstat_size::Int64=100000,bootstrap::Int64=0)
	
	@eval  using JLD2, DataFrames, CSV
	
	sfs     = CSV.read.(filter(x -> occursin("sfs",x), readdir(analysis_folder,join=true)),DataFrame,header=false) .|> Array;
	divergence     = CSV.read.(filter(x -> occursin("div",x), readdir(analysis_folder,join=true)),DataFrame,header=false) .|> Array;

	@eval adap      = MKtest.parameters(N=$ne,n=$samples,dac =parse.(Int,split($dac,",")))

	@eval if $bootstrap != 0
		@eval summstat  = MKtest.summary_statistics(
            param           = adap,
            h5_file         = $h5_file,
            sfs             = $sfs,
            divergence      = $divergence,
            analysis_folder = $analysis_folder,
            summstat_size   = $summstat_size,
            bootstrap       = $bootstrap
		);
	else
		@eval summstat  = MKtest.summary_statistics(
            param           = $adap,
            h5_file         = $h5_file,
            sfs             = $sfs,
            divergence      = $divergence,
            analysis_folder = $analysis_folder,
            summstat_size   = $summstat_size
		);
	end
end

"ABCreg inference.

The function returns posterior distributions from ABC inference. Each posterior file contains information about alpha_w, alpha_s, alpha, gam_neg and shape parameter. The number of posterior distributions will depend on the number of bootstrap replicas.

Check the documentation to get more info https://jmurga.github.io/MKtest.jl/dev/cli
"
@main function inference(;analysis_folder::String="<folder>",P::Int64=5,S::Int64=9,tol::Float64=0.025,ABCreg::String="/home/jmurga/ABCreg/src/reg")
	
	@eval MKtest.ABCreg(
        analysis_folder = $analysis_folder,
        P               = $P,
        S               = $S,
        tol             = $tol,
        abcreg          = $ABCreg
	);

end

#="Plot Maximum a posterior distribution"
@main function plotMap(;analysis_folder::String="<folder>")
	try
		@eval using RCall, GZip, DataFrames, CSV
		
		@eval MKtest.sourcePlotMapR(script=$analysis_folder)
		@eval MKtest.plotMap(analysis_folder=$analysis_folder)
		@eval RCall.endEmbeddedR()
	catch
		println("Please install R, ggplot2, data.table and locfit in your system before execute this function")
	end
end=#

