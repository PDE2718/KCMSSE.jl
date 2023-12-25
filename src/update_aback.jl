function update_aback!(l0::Leg, ξ::f64, μ::f64, PBC::Bool)::Bool
    # if l0 |> get_flag == I0_ || ~iscenter(l0)
    #     return false
    # end
    tail::Leg = head::Leg = l0
    wr::f64 = 1.0
    across::Bool = cyclic::Bool = false
    while true
        if ~PBC && (tail |> get_prev |> isnull)
            across = true
            wr *= wf2wi_head(head, ξ, μ)
            break
        end
        tail = tail |> get_prev

        if PBC && (tail == head)
            cyclic = true
            wr *= wf2wi_cyclic(tail, ξ, μ)
            break
        elseif tail |> iscenter
            wr *= wf2wi_tail(tail, ξ, μ)
            wr *= wf2wi_head(head, ξ, μ)
            break
        else # normal continue
            wr *= wf2wi_body(tail, ξ, μ)
            if iszero(wr)
                return false
            end
        end
    end
    @assert (across && cyclic) == false
    if wr > 0 && metro(wr)
        if across
            flip_flag!(head)
            flip_ψ!(head) # one more flip compare to ahead!!!!
            while true
                flip_ψ!(head)
                head = head |> get_prev
                if head |> isnull
                    return true
                end
            end
        else
            flip_flag!(tail)
            while true
                flip_ψ!(tail)
                tail = tail |> get_next
                if tail == head
                    flip_flag!(tail)
                    return true
                end
            end
        end
    else
        return false
    end
end