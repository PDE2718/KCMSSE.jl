using Revise
using KCMSSE, Statistics

function onesimu(L::Tuple{Int,Int}, β::Float64, μ::Float64;
    t_bin::Int=1000, n_bin::Int=32,
    t_th::Int=30000, PBC=false
)
    N::Int = prod(L)
    Λ0::Int = ceil(Int, β * N)
    ψ00::Matrix{Bool} = zeros(Bool, L)
    ψ00[1] = true

    X = Estimator(ψ00, Λ0; PBC=PBC)
    X.β = β
    X.μ = μ
    X.ξ = 0.0
    M = [Obs(X) for t ∈ 1:n_bin]

    for t ∈ 1:t_th
        bisweep!(X)
        increase_estimator!(X, ψ00, 1.5)
    end
    for b ∈ M
        for t ∈ 1:t_bin
            bisweep!(X)
            onestep_measure!(X, b; track_trace=false, take_snapshot=iszero(t % 10))
        end
    end
    res = [dumptoNamedTuple(m) for m ∈ M]
    return res
end

res = onesimu((5, 5), 32.0, 0.0)

mean(x -> mean(x.E), res) # ≈ -0.394742
mean(x -> mean(x.ρ), res) # ≈  0.374215


















