# project + reduce
function project_moments(sfs::Union{DataFrame,Matrix{Float64}},n::Int64)

    moments = PythonCall.pyimport("moments")

    @assert size(sfs,2) > 2 "Check SFS"
    
    _nn = 2*n
    _n = _nn - 1

    tmp = zeros(_n,3)
    tmp[:,1] .= (1:_n)
    for i=2:3
        tmp[:,i] = pyconvert(Vector{Float64},moments.Spectrum(view(sfs,:,i)).project([_nn]))[2:end-1]
    end

    tmp[:,2:3] .= floor.(tmp[:,2:3])
    return tmp
end


true_alpha = DataFrame[]
true_omega = DataFrame[]
grapes     = DataFrame[]
abc_mk     = DataFrame[]

function read_simulations(param,folder)

    # mkpath(output);
    tmp = split(folder,"/")[end]
    
    df_div = CSV.read.(filter(x->occursin("div",x),readdir(folder,join=true)),DataFrame,header=true);
    
    ln,ls = [7.5e7,2.5e7]

    df_div = vcat(df_div...)
    alpha = DataFrame([df_div[:,3]./df_div[:,1],df_div[:,4]./df_div[:,1],sum.(eachrow(df_div[:,3:4]))./df_div[:,1],fill(analysis,size(df_div,1)),fill("slim",size(df_div,1))],[:α_weak,:α_strong,:α,:analysis,:method])

    dₙ = @. (df_div[:,1]/ln)
    dₛ = @. (df_div[:,2]/ls)
    D₋ = @. df_div[:,1] - df_div[:,3] - df_div[:,4]

    ω  = @. dₙ / dₛ
    ωₙ = @. (D₋/ln) / dₛ
    ωₐ_weak   = @. (df_div[:,3]/ln) / dₛ
    ωₐ_strong =  @. (df_div[:,4]/ln) / dₛ
    ωₐ   = @. ((df_div[:,3]+df_div[:,4])/ln) / dₛ

    omega_a = DataFrame([ωₐ_weak,ωₐ_strong,ωₐ,fill(analysis,size(df_div,1)),fill("slim",size(df_div,1))],[:ωₐ_weak,:ωₐ_strong,:ωₐ,:analysis,:method])

    divergence = Int64.(Array(df_div)[:,1:2])
    divergence = hcat(divergence,repeat([ln ls],size(divergence,1)))

    divergence  = [permutedims(divergence[i,:]) for i in 1:size(divergence,1)]

    sfs = CSV.read.(filter(x->occursin("sfs",x),readdir(folder,join=true)),DataFrame,header=true);
    sfs = map(x -> Array(x)[:,1:3],sfs)
    map(x -> x[:,1] .= 1:(param.nn-1),sfs)
    return(sfs,divergence)
end

analysis = filter(x-> isdir(x),readdir(labstorage * "raw_data/no_demog",join=true))

adap = MKtest.parameters(N=500,n=500,shape=0.184,dac=[1,2,4,5,10,20,50,100,200,400,500,700])
adap_20 = MKtest.parameters(N=500,n=20)



for folder in analysis
    @show folder
    
    sfs,divergence = read_simulations(adap,folder)
    sfs_projected = map(x-> project_moments(x,20),sfs)


    for i in [0.5,0.25,0.1,0.01,0.001,0.0001]
        sfs_tmp = map(x->sum(x[:,2:3].*i),sfs)
        sfs_projected_tmp = map(x->sum(x[:,2:3].*i),sfs_projected)
        divergence_tmp = map(x->sum(x.*i),divergence)

        MKtest.imputedMK(adap,sfs,divergence,cutoff=0.3)[1]
        MKtest.imputedMK(adap,sfs_projected,divergence,cutoff=0.3)
        MKtest.fwwMK(adap,sfs,divergence,cutoff=0.3)
        MKtest.fwwMK(adap_20,sfs_projected,divergence,cutoff=0.3)

        MKtest.aMK(adap,sfs,divergence,na_rm=true)
        MKtest.aMK(adap_20,sfs_projected,divergence,na_rm=true)
    end

end
