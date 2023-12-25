mutable struct Accum{T}
    X::T
    num::Int
end
function Accum(x)
    return Accum(zero(x), 0)
end
import Base: push!, empty!
function push!(A::Accum{T}, x::T) where {T<:Number}
    A.X += x
    A.num += 1
    return nothing
end
function push!(A::Accum{T1}, x::T2) where {T1<:AbstractArray,T2<:AbstractArray}
    A.X .+= x
    A.num += 1
    return nothing
end
function empty!(A::Accum{T}) where {T}
    A.X = 0
    A.num = 0
    return nothing
end
import Statistics: mean
function mean(A::Accum{T}) where {T}
    return A.X ./ A.num
end
