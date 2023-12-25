function update_ahead!(l0::Leg, ξ::f64, μ::f64, PBC::Bool)::Bool
    # if l0 |> get_flag == I0_ || ~iscenter(l0)
    #     return false
    # end
    tail::Leg = head::Leg = l0
    wr::f64 = 1.0
    across::Bool = cyclic::Bool = false
    while true
        if ~PBC && (head |> get_next |> isnull)
            across = true
            wr *= wf2wi_tail(tail, ξ, μ)
            break
        end
        head = head |> get_next

        if PBC && (head == tail)
            cyclic = true
            wr *= wf2wi_cyclic(head, ξ, μ)
            break
        elseif head |> iscenter
            wr *= wf2wi_tail(tail, ξ, μ)
            wr *= wf2wi_head(head, ξ, μ)
            break
        else # normal continue
            wr *= wf2wi_body(head, ξ, μ)
            if iszero(wr)
                return false
            end
        end
    end
    @assert (across && cyclic) == false
    if wr > 0 && metro(wr)
        if across
            flip_flag!(tail)
            # one less flip compare to aback!!!!
            while true
                flip_ψ!(tail)
                tail = tail |> get_next
                if tail |> isnull
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

