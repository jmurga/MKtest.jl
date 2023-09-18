using MKtest
using Test

# Running MKtest approximation with default parameters
@testset "MKtest.jl" begin

	adap  = MKtest.parameters(N=1000,n=661,gam_dfe=-457,gam_flanking=-1000,al_low=0.2,B=0.999)
	binom = MKtest.binom_op(adap)
	MKtest.set_θ!(adap)
	MKtest.set_ppos!(adap)

	@test adap.ppos_h == 0.00020733640485526956
	@test adap.θ_flanking == 3.001501000750676e-6

end
