function Weight(rs::RootSystem, coordinate_vector, basis::Symbol)
    if basis == :functionals
        return Weight(rs, missing, missing, coordinate_vector)
    elseif basis == :fundamental_weights
        return Weight(rs, missing, coordinate_vector, missing)
    elseif basis == :roots
        return Weight(rs, coordinate_vector, missing, missing)
    end

    throw("Invalid basis")
end

function as_root(w::Weight)
    if ismissing(w.as_root)
        if !ismissing(w.as_functional)
            w.as_root = get_functional_to_root_matrix(w.rs) * w.as_weight   
        elseif !ismissing(w.as_fundamental_weight)
            w.as_root = get_cartan_matrix_inv(w.rs) * w.as_fundamental_weight
        end
        # unreachable
    end

    return w.as_root
end

function as_fundamental_weight(w::Weight)
    if ismissing(w.as_fundamental_weight)
        if !ismissing(w.as_functional)
            w.as_fundamental_weight = get_functional_to_fundamental_weight_matrix(w.rs) * w.as_functional
        elseif !ismissing(w.as_root)
            w.as_fundamental_weight = get_cartan_matrix(w.rs) * w.as_root
        end
        # unreachable
    end

    return w.as_fundamental_weight
end

function as_functional(w::Weight)
    if ismissing(w.as_functional)
        if !ismissing(w.as_fundamental_weight)
            w.as_functional = get_fundamental_weight_to_functional_matrix(w.rs) * w.as_fundamental_weight
        elseif !ismissing(w.as_root)
            w.as_functional = get_root_to_functional_matrix(w.rs) * w.as_root
        end
        # unreachable
    end

    return w.as_functional
end

function convert(w::Weight, type::Symbol)
    if type == :functionals
        return as_functional(w)
    elseif type == :fundamental_weights
        return as_fundamental_weight(w)
    elseif type == :roots
        return as_root(w)
    end
end

# returns <w, a>, only valid if a is a root! 
function lf(w::Weight, a::Weight)
    w_fw = as_fundamental_weight(w)
    a_r = as_root(a)

    
    #println("LF: ", as_fundamental_weight(w), " * ", as_root(a), " = ", dot(w_fw, a_r))

    return Int(dot(w_fw, a_r))
end

function Base.:+(w1::Weight, w2::Weight)
    @assert(w1.rs == w2.rs)
    # make sure that not all sums will be missing
    w1_fw = as_fundamental_weight(w1)
    w2_fw = as_fundamental_weight(w2)
    return Weight(w1.rs, w1.as_root .+ w2.as_root, w1_fw .+ w2_fw, w1.as_functional .+ w2.as_functional )
end

function Base.:-(w1::Weight, w2::Weight)
    @assert(w1.rs == w2.rs)
    # make sure that not all sums will be missing
    w1_fw = as_fundamental_weight(w1)
    w2_fw = as_fundamental_weight(w2)
    return Weight(w1.rs, w1.as_root .- w2.as_root, w1_fw .- w2_fw, w1.as_functional .- w2.as_functional )
end

function Base.:*(m::T, w::Weight) where {T <: Integer}
    return Weight(w.rs, m .* w.as_root, m .* w.as_fundamental_weight, m .* w.as_functional )
end