using Base: unsafe_getindex, unsafe_setindex!
import Base: getindex, setindex!, length
import Base: >, <, ≥, ≤, ==, ≠
const i32 = Int32
const f64 = Float64

struct HString
    buf::Vector{i32}
end
Base.getindex(H::HString, i::Integer) = unsafe_getindex(H.buf, i)
Base.setindex!(H::HString, v::Integer, i::Integer) = unsafe_setindex!(H.buf, v, i)
Base.getindex(H::HString, p::Integer, q::Integer) = unsafe_getindex(H.buf, ((p - 1) << 4) + q)
Base.setindex!(H::HString, v::Integer, p::Integer, q::Integer) = unsafe_setindex!(H.buf, v, ((p - 1) << 4) + q)
Base.length(H::HString) = length(H.buf) >> 4
string_length(H::HString) = length(H)
string_orders(H::HString) = Base.OneTo(length(H))
function HString(N::Integer)
    H = HString(zeros(i32, N << 4))
    for p ∈ string_orders(H)
        H[p, 16] = p << 4
    end
    return H
end

struct Leg
    H::HString
    p::i32
    r::i32
end

id_self(p::Integer, r::Integer)::i32 = ((p - 1) << 4) + r
id_op(p::Integer)::i32 = ((p - 1) << 4) + 16
const prev_shift_::i32 = 5
const next_shift_::i32 = 10

mod16(x::Integer)::i32 = mod(x, 16)
parse_p(p_r::Integer)::i32 = p_r >> 4
parse_r(p_r::Integer)::i32 = mod16(p_r)
parse_p_r(p_r::Integer)::Tuple{i32,i32} = (parse_p(p_r), parse_r(p_r))
zip_p_r(p::Integer, r::Integer)::i32 = (p << 4) + r

Leg(l::Leg, p::Integer, r::Integer) = Leg(l.H, p, r)
Leg(l::Leg, p_r::Integer) = Leg(l.H, parse_p(p_r), parse_r(p_r))
Base.getindex(l::Leg, rp::Integer) = Leg(l.H, l.p, rp)
get_self_p_r(l::Leg)::i32 = zip_p_r(l.p, l.r)

function get_next(l::Leg)::Leg
    next_p_r = l.H[id_self(l.p, l.r)+next_shift_]
    return Leg(l, next_p_r)
end
function set_next!(l::Leg, p_r::i32)::i32
    l.H[id_self(l.p, l.r)+next_shift_] = p_r
end
set_next!(l::Leg, p::i32, r::i32)::i32 = set_next!(l, zip_p_r(p, r))
set_next!(l::Leg, p_r::Tuple{i32,i32})::i32 = set_next!(l, zip_p_r(p_r[1], p_r[2]))
set_next!(l::Leg, lp::Leg) = set_next!(l, lp |> get_self_p_r)

function get_prev(l::Leg)::Leg
    prev_p_r = l.H[id_self(l.p, l.r)+prev_shift_]
    return Leg(l, prev_p_r)
end
function set_prev!(l::Leg, p_r::i32)::i32
    l.H[id_self(l.p, l.r)+prev_shift_] = p_r
end
set_prev!(l::Leg, p::i32, r::i32)::i32 = set_prev!(l, zip_p_r(p, r))
set_prev!(l::Leg, p_r::Tuple{i32,i32})::i32 = set_prev!(l, zip_p_r(p_r[1], p_r[2]))
set_prev!(l::Leg, lp::Leg) = set_prev!(l, lp |> get_self_p_r)

checkout_p(l::Leg)::i32 = l.H[id_op(l.p)] >> 4
checkout(l::Leg)::Bool = (l.H[id_op(l.p)] >> 4) == l.p

@enum Vert::Int8 I0_ = 0 Iμ_ = 1 σx_ = 2
get_flag(l::Leg)::Vert = l.H[id_op(l.p)] |> mod16 |> Vert
function set_flag!(l::Leg, v::Vert)::Vert
    pid = id_op(l.p)
    l.H[pid] = ((l.H[pid] >> 4) << 4) + i32(v)
    return v
end
function flip_flag!(l::Leg)::Vert
    v = get_flag(l)
    set_flag!(l, v == Iμ_ ? σx_ : Iμ_)
end
function get_self(l::Leg)::i32
    return l.H[id_self(l.p, l.r)]
end
function set_self!(l::Leg, j_ψ::Integer)::i32
    return l.H[id_self(l.p, l.r)] = j_ψ
end

parse_j(j_ψ::Integer)::i32 = j_ψ >> 4
parse_ψ(j_ψ::Integer)::Bool = isodd(j_ψ)
parse_j_ψ(j_ψ::Integer)::Tuple{i32,Bool} = (j_ψ >> 4, isodd(j_ψ))
zip_j_ψ(j::Integer, ψ::Bool)::i32 = (j << 4) + ψ

get_j(l::Leg)::i32 = l.H[id_self(l.p, l.r)] |> parse_j
get_ψ(l::Leg)::Bool = l.H[id_self(l.p, l.r)] |> parse_ψ
get_j_ψ(l::Leg)::Tuple{i32,Bool} = l.H[id_self(l.p, l.r)] |> parse_j_ψ

function get_ψ_before(l::Leg)
    ψ::Bool = get_ψ(l)
    if (l |> iscenter) && ((l |> get_flag) == σx_)
        return ~ψ
    else
        return ψ
    end
end

function set_j!(l::Leg, j::Integer)::i32
    id = id_self(l.p, l.r)
    ψ = l.H[id] |> parse_ψ
    l.H[id] = zip_j_ψ(j, ψ)
    return j
end
function set_ψ!(l::Leg, ψ::Bool)::Bool
    id = id_self(l.p, l.r)
    j = l.H[id] |> parse_j
    l.H[id] = zip_j_ψ(j, ψ)
    return ψ
end
function set_j_ψ!(l::Leg, j::Integer, ψ::Bool)::Tuple{i32,Bool}
    l.H[id_self(l.p, l.r)] = zip_j_ψ(j, ψ)
    return (j, ψ)
end
function flip_ψ!(l::Leg)::Bool
    id = id_self(l.p, l.r)
    j_ψ = l.H[id]
    if j_ψ |> parse_ψ
        l.H[id] = j_ψ - 1
        return false
    else
        l.H[id] = j_ψ + 1
        return true
    end
end
function get_count(l::Leg)::i32
    id::i32 = id_self(l.p, 0)
    cnt::i32 = +(
        parse_ψ(l.H[id+1]),
        parse_ψ(l.H[id+2]),
        parse_ψ(l.H[id+3]),
        parse_ψ(l.H[id+4]),
    )
    return cnt
end
function get_count_flip(l::Leg)::Tuple{i32,i32}
    id::i32 = id_self(l.p, 0)
    nbconf = (
        parse_ψ(l.H[id+1]),
        parse_ψ(l.H[id+2]),
        parse_ψ(l.H[id+3]),
        parse_ψ(l.H[id+4]),
    )
    cnt::i32 = sum(nbconf)
    if l.r == 5
        return (cnt, cnt)
    else
        ψt::Bool = nbconf[l.r]
        ψf::Bool = ~ψt
        return (cnt, cnt + ψf - ψt)
    end
end
function get_ψ_udlrx(l::Leg)
    id::i32 = id_self(l.p, 0)
    return (
        parse_ψ(l.H[id+1]),
        parse_ψ(l.H[id+2]),
        parse_ψ(l.H[id+3]),
        parse_ψ(l.H[id+4]),
        parse_ψ(l.H[id+5]),
    )
end

isnull(l::Leg) = iszero(l.r)
nullleg(l::Leg) = Leg(l.H, 0, 0)
nullleg(H::HString) = Leg(H, 0, 0)
iscenter(l::Leg)::Bool = l.r == 5
==(l1::Leg, l2::Leg)::Bool = (l1.p == l2.p)
≠(l1::Leg, l2::Leg)::Bool = (l1.p ≠ l2.p)
# ==(l1::Leg, l2::Leg)::Bool = l1.p == l2.p
<(l1::Leg, l2::Leg)::Bool = l1.p < l2.p
>(l1::Leg, l2::Leg)::Bool = l1.p > l2.p
≤(l1::Leg, l2::Leg)::Bool = l1.p ≤ l2.p
≥(l1::Leg, l2::Leg)::Bool = l1.p ≥ l2.p

function center_leg(H, p)
    return Leg(H, p, 5)
end

ψ2char(ψ::Bool) = ψ ? '◼' : '□'
ψ2charflip(ψ::Bool) = ψ ? '⬓' : '⬒'
function Base.show(io::IO, ll::Leg)
    if ll |> isnull
        print(io, "Leg{undef} at p=$(ll.p)")
    elseif ll |> get_flag == I0_
        print(io, "⊠⊠⊠⊠ ⊠ Leg{I} at p=$(ll.p)")
    elseif ll |> get_flag == Iμ_
        udlrx_ = get_ψ_udlrx(ll)
        printstyled(io, ψ2char(udlrx_[1]), color=ll.r == 1 ? :red : :cyan)
        printstyled(io, ψ2char(udlrx_[2]), color=ll.r == 2 ? :red : :cyan)
        printstyled(io, ψ2char(udlrx_[3]), color=ll.r == 3 ? :red : :cyan)
        printstyled(io, ψ2char(udlrx_[4]), color=ll.r == 4 ? :red : :cyan)
        printstyled(io, ' ', ψ2char(udlrx_[5]), color=ll.r == 5 ? :red : :cyan)
        print(io, " Leg{μ} at p=$(ll.p), c=$(get_j(ll[5]))")
    elseif ll |> get_flag == σx_
        udlrx_ = get_ψ_udlrx(ll)
        printstyled(io, ψ2char(udlrx_[1]), color=ll.r == 1 ? :red : :cyan)
        printstyled(io, ψ2char(udlrx_[2]), color=ll.r == 1 ? :red : :cyan)
        printstyled(io, ψ2char(udlrx_[3]), color=ll.r == 1 ? :red : :cyan)
        printstyled(io, ψ2char(udlrx_[4]), color=ll.r == 1 ? :red : :cyan)
        printstyled(io, ' ', ψ2charflip(udlrx_[5]), color=ll.r == 1 ? :red : :cyan)
        print(io, " Leg{σx} at p=$(ll.p), c=$(get_j(ll[5]))")
    end
end
