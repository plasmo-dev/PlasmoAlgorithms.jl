module PlasmoSchwarz

using Printf
using DataStructures
using LinearAlgebra
using Metis
import Base.@kwdef

using MathOptInterface, Plasmo, JuMP

export SchwarzAlgorithm, run_algorithm!, calculate_primal_feasibility, calculate_dual_feasibility

include("utils.jl")

include("algorithm.jl")

include("interface.jl")

end
