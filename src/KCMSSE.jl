module KCMSSE
using LinearAlgebra, FFTW, Statistics

include("lattice.jl")
include("Leg.jl")
include("Estimator.jl")
include("update.jl")

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end # module KCMSSE
