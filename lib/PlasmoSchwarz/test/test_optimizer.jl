include((@__DIR__) * "/../examples/optimal_control.jl")

@test Plasmo.termination_status(optimizer) == MOI.LOCALLY_SOLVED
