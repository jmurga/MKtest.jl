module MKtest

	using Parameters, SparseArrays, Distributed, CSV, JLD2, DataFrames, ProgressMeter, Quadmath, GZip, StatsBase, Unzip, ThreadsX, ParallelUtilities

	# Analytical solutions
	import Roots: find_zero
	import NLsolve: nlsolve
	import SpecialFunctions: polygamma, zeta
	import PoissonRandom: pois_rand
	import Distributions: Binomial, pdf

	# Parse data
	import GZip: open
	import Parsers: parse
	import OrderedCollections: OrderedDict
	import FastaIO: readfasta
	import Random: randstring
	import Suppressor: @suppress_out

	# MK-approaches
	import LsqFit: curve_fit, confidence_interval
	import HypothesisTests: pvalue,FisherExactTest

	include("parameters.jl")
	include("fixations.jl")
	include("polymorphism.jl")
	include("summary_statistics.jl")
	include("rates.jl")
	include("infer_tools.jl")
	include("read_fasta.jl")
	include("methods.jl")

end
