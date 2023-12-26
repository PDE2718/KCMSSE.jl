mutable struct Estimator
    const PBC::Bool # if this is finite temp
    n::Int
    β::f64
    ξ::f64
    μ::f64
    const H::HString
    const legs0::Matrix{Leg}
    const legsβ::Matrix{Leg}
    const ψ0::Matrix{Bool}
    const ψt::Matrix{Bool}
    const ψm::Matrix{Bool}
    ######################## observable buffers
    const ψk::Matrix{ComplexF64}
    const Sk::Matrix{f64}
    ######################## other constants
end
function Estimator(ψini::Matrix{Bool}, Λ0::Int; PBC::Bool=true)
    @assert Λ0 |> iseven
    H = HString(Λ0)
    ψ0 = deepcopy(ψini)
    ψt = PBC ? ψ0 : deepcopy(ψini)
    ψm = deepcopy(ψini)
    X = Estimator(PBC,
        0, 0.0, 0.0, 0.0,
        H,
        [Leg(H, 0, 0) for _ ∈ ψ0],
        [Leg(H, 0, 0) for _ ∈ ψ0],
        ψ0, ψt, ψm,
        zeros(ComplexF64, size(ψ0)),
        zeros(Float64, size(ψ0)),
    )
    return X
end
function inspect_wl(X::Estimator)
    for p ∈ string_orders(X.H)
        Leg(X.H, p, 5) |> println
    end
end
function inspect_kinks(X::Estimator)
    for p ∈ string_orders(X.H)
        l = Leg(X.H, p, 5)
        v = get_flag(l)
        if v == σx_
            l |> println
        end
    end
end

function increment!(X::Estimator, ψ0::Matrix{Bool}, r::Real=1.5)::Bool
    @assert r > 1
    H::HString = X.H
    n::Int = X.n
    Λ0::Int = string_length(H)
    Λ1::Int = ceil(Int, n * r)
    Λp::Int = ceil(Int, n * abs2(r))
    Λ1 += isodd(Λ1)
    Λp += isodd(Λp)
    if Λ1 > Λ0
        resize!(H.buf, Λp)
        ini_string!(H)
        X.ψ0 .= ψ0
        X.ψt .= ψ0
        X.ψm .= ψ0
        X.n = 0
        return true
    else
        return false
    end
end