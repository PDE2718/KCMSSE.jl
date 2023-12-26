include("accumulator.jl")
Base.@kwdef struct Obs
    t::Int64
    L::Tuple{Int64,Int64}
    β::Float64
    ξ::Float64
    μ::Float64
    E::Accum{Float64} = Accum(0.0)
    ρ::Accum{Float64} = Accum(0.0)
    nwl1::Accum{Int} = Accum(0)
    nwl2::Accum{Int} = Accum(0)
    Sz1::Accum{Int} = Accum(0)
    Sz2::Accum{Int} = Accum(0)
    Sk::Accum{Matrix{Float64}} = Accum(zeros(Float64, L), 0)
    ψ̄::Accum{Matrix{Float64}} = Accum(zeros(Float64, L), 0)
    ρtrace::Vector{Float64} = Float64[]
    ψsnapshots::Vector{BitMatrix} = BitMatrix[]
end
Obs(X::Estimator, t::Int64=0) = Obs(
    t=t, L=size(X.ψ0), β=X.β, ξ=X.ξ, μ=X.μ
)

function onestep_measure!(X::Estimator, O::Obs;
    track_trace::Bool=false,
    take_snapshot::Bool=false
)
    ψm, ψk, Sk, β, ξ, μ = X.ψm, X.ψk, X.Sk, X.β, X.ξ, X.μ

    N = length(ψm)
    # record the world line number and energy
    nwl1::Int64 = X.n
    push!(O.nwl1, nwl1)
    E = -nwl1 / β / N + μshift(μ)
    push!(O.E, E)
    # for specific heat
    nwl2::Int64 = abs2(nwl1)
    push!(O.nwl2, nwl2)

    # magnetization
    Sz1 = sum(ψm)
    push!(O.Sz1, Sz1)
    push!(O.Sz2, Sz1 |> abs2)

    # density
    ρ = Sz1 / N
    push!(O.ρ, ρ)

    # real space mean and the structure factor
    ψk .= ψm
    push!(O.ψ̄, ψm)
    fft!(ψk)
    map!(abs2, Sk, ψk)
    push!(O.Sk, Sk)

    # record optional data
    if track_trace
        push!(O.ρtrace, ρ)
    end
    if take_snapshot
        push!(O.ψsnapshots, X.ψm .== true)
    end

    return nothing
end

obs_ninj(Sk::Matrix{Float64})::Matrix{Float64} = real.(bfft(Sk)) ./ (length(Sk)^2)

function obs_cicj(Sk::Matrix{Float64})::Matrix{Float64}
    S̃ = deepcopy(Sk)
    S̃[1] = 0
    return real.(bfft(S̃)) ./ (length(S̃)^2)
end
# Auxiliary functions for visualization
function pbcfill(A)
    B = fftshift(A)
    if iseven(size(B, 1))
        B = vcat(B, B[1:1, 1:end])
    end
    if iseven(size(B, 2))
        B = hcat(B, B[1:end, 1:1])
    end
    return B
end
function pbcgrid(A)
    @assert size(A) |> allequal
    n = size(A, 1)
    @assert isodd(n)
    return -(n ÷ 2):(n÷2)
end