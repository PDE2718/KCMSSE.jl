μshift(μ::f64)::f64 = 1.0 + abs(μ)
Iμ_weight(ψ::Bool, μ::f64)::f64 = μshift(μ) + μ * ψ
σx_weight(cnt::Integer, ξ::f64)::f64 = cnt == 1 ? (1.0 - ξ) : ξ

function diagonal_update!(X::Estimator)
    X.n = diagonal_update!(X.PBC, X.H, X.n, X.β, X.ξ, X.μ, X.legs0, X.legsβ, X.ψ0, X.ψt, X.ψm)
end
function diagonal_update!(PBC::Bool, H::HString,
    n::Int, β::f64, ξ::f64, μ::f64,
    legs0::Matrix{Leg}, legsβ::Matrix{Leg},
    ψ0::Matrix{Bool}, ψt::Matrix{Bool}, ψm::Matrix{Bool},
    )::Int
    
    L::Tuple{Int,Int} = size(ψ0)
    N::Int = length(ψ0)
    lattice = eachindex(ψ0)
    Λ::Int = string_length(H)
    M::Int = Λ ÷ 2
    ψt .= ψ0

    fill!(legs0, nullleg(H))
    fill!(legsβ, nullleg(H))

    l::Leg = nullleg(H)
    ic::Int = 0
    w::f64 = 0.0
    v::Vert = I0_
    @inbounds for p ∈ string_orders(H)
        l = center_leg(H, p)
        v = get_flag(l)
        if v == I0_
            ic = rand(lattice)
            w = Iμ_weight(ψt[ic], μ)
            if metro((N * w * β) / (Λ - n))
                set_flag!(l, Iμ_)
                set_j!(l, ic)
                n += 1
            end

        elseif v == Iμ_
            ic = get_j(l)
            w = Iμ_weight(ψt[ic], μ)
            if metro((Λ - n + 1) / (N * w * β))
                set_flag!(l, I0_)
                # set_j!(l, 0)
                n -= 1
            end

        elseif v == σx_
            ic = get_j(l)
            ψt[ic] ⊻= true
        end

        if p == M
            ψm .= ψt
        end

        if get_flag(l) ≠ I0_
            ic = get_j(l)
            for (r, j) ∈ enumerate(udlrx(ic, L))
                l = l[r]
                set_j!(l, j)
                set_ψ!(l, ψt[j])
                if legs0[j] |> isnull
                    legs0[j] = l
                    legsβ[j] = l
                else
                    set_prev!(l, legsβ[j])
                    set_next!(legsβ[j], l)
                    legsβ[j] = l
                end
            end
        end
    end
    for (l0, lβ) ∈ zip(legs0, legsβ)
        if ~isnull(l0) && ~isnull(lβ)
            set_prev!(l0, PBC ? lβ : nullleg(l))
            set_next!(lβ, PBC ? l0 : nullleg(l))
        end
    end
    return n
end

function wf2wi_tail(l::Leg, ξ::f64, μ::f64)::f64
    ψ::Bool = l |> get_ψ
    c::i32 = l |> get_count
    v::Vert = l |> get_flag
    if v == Iμ_
        return σx_weight(c, ξ) / Iμ_weight(ψ, μ)
    elseif v == σx_
        return Iμ_weight(~ψ, μ) / σx_weight(c, ξ)
    end
end

function wf2wi_head(l::Leg, ξ::f64, μ::f64)::f64
    ψ::Bool = l |> get_ψ
    c::i32 = l |> get_count
    v::Vert = l |> get_flag
    if v == Iμ_
        return σx_weight(c, ξ) / Iμ_weight(ψ, μ)
    elseif v == σx_
        return Iμ_weight(ψ, μ) / σx_weight(c, ξ)
    end
end

function wf2wi_body(l::Leg, ξ::f64, μ::f64)::f64
    v::Vert = l |> get_flag
    if v == Iμ_
        return 1.0
    elseif v == σx_
        if iszero(ξ)
            return 0.0
        end
        ψ::Bool = l |> get_ψ
        ci::i32 = l |> get_count
        cf::i32 = ci + (~ψ) - ψ
        return σx_weight(cf, ξ) / σx_weight(ci, ξ)
    end
end

function wf2wi_cyclic(l::Leg, ξ::f64, μ::f64)::f64
    v::Vert = l |> get_flag
    if v == Iμ_
        ψ::Bool = l |> get_ψ
        if ~iszero(σx_weight(l |> get_count, ξ))
            return Iμ_weight(~ψ, μ) / Iμ_weight(ψ, μ)
        else
            return 0.0
        end
    elseif v == σx_
        return 1.0
    end
end

include("update_ahead.jl")
include("update_aback.jl")

function update_ψ0t!(X::Estimator)
    for (i, l) ∈ enumerate(X.legs0)
        if ~isnull(l)
            X.ψ0[i] = get_ψ_before(l)
        end
    end
    for (i, l) ∈ enumerate(X.legsβ)
        if ~isnull(l)
            X.ψt[i] = get_ψ(l)
        end
    end
    return nothing
end

function offdiag_update!(H::HString, ξ::Float64, μ::Float64, PBC::Bool, update_forward::Bool)
    Λ = string_length(H)
    p_orders = string_orders(H)
    if update_forward
        for i ∈ p_orders
            update_ahead!(center_leg(H, i), ξ, μ, PBC)
        end
    else
        for i ∈ reverse(p_orders)
            update_aback!(center_leg(H, i), ξ, μ, PBC)
        end
    end
    return nothing
end

function offdiag_update!(X::Estimator)
    offdiag_update!(X.H, X.ξ, X.μ, X.PBC, rand(Bool))
    update_ψ0t!(X)
    return nothing
end

function sweep!(X::Estimator)
    diagonal_update!(X)
    offdiag_update!(X)
    return nothing
end

function bisweep!(X::Estimator)
    ξ::f64, μ::f64, H::HString, PBC::Bool = X.ξ, X.μ, X.H, X.PBC
    Λ::Int = string_length(H)
    p_orders = string_orders(H)
    
    diagonal_update!(X)
    for i ∈ p_orders
        update_ahead!(center_leg(H, i), ξ, μ, PBC)
    end
    update_ψ0t!(X)

    diagonal_update!(X)
    for i ∈ reverse(p_orders)
        update_aback!(center_leg(H, i), ξ, μ, PBC)
    end
    update_ψ0t!(X)

    return nothing
end
