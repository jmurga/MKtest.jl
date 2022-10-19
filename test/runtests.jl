using MKtest
using Test

# Running MKtest approximation with default parameters
@testset "MKtest.jl" begin

	adap  = MKtest.parameters(N=1000,n=661,gam_dfe=-457,gam_flanking=-1000,al_low=0.2,B_bins=[0.2,0.999])
	binom = MKtest.binom_op(adap)
	MKtest.set_θ!(adap)
	MKtest.set_ppos!(adap)

	@test adap.ppos_h == 0.00020733640485526953
	@test adap.θ_flanking == 3.0015010007506417e-6

end