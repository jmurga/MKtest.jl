using MKtest
using Test

# Running MKtest approximation with default parameters
@testset "MKtest.jl" begin

	adap  = MKtest.parameters(N=1000,n=661,al_low=0.2,B_bins=[0.2,0.999])
	binom = MKtest.binom_op!(adap)
	MKtest.set_θ!(adap)
	MKtest.set_ppos!(adap)

	@test adap.ppos_h == 0.00020805325103594085
	@test adap.θ_noncoding == 1.6433217979109007e-6

end