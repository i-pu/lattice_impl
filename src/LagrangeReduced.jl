#=
アルゴリズム42 Lagrange簡約
=#
function ragrange_reduced(b1, b2)
    if norm(b1) > norm(b2)
        b1, b2 = b2, b1
    end

    while true
        q = -round(b1 ⋅ b2 / norm(b1) ^ 2)
        b2, b1 = b1, b2 + q * b1
        if !(norm(b1) < norm(b2))
            break
        end
    end
    b2, b1
end