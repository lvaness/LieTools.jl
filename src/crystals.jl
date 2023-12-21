function Base.show(io::IO, cb::CrystalBase)
    print(io, collect(cb.B))
end

function vector_rep_crystal_type_a(n)
    B = collect(1:n+1)
    e = function(i, j) # e_i * v_j
        return j == i+1 ? B[i] : 0
    end
    f = function(i, j)
        return j == i ? B[i+1] : 0
    end
    eps = function(i, j)
        return j == i+1 ? 1 : 0
    end
    phi = function(i, j)
        return j == i ? 1 : 0
    end

    return CrystalBase(B, n+1, n, e, f, eps, phi, missing) # TODO replace missing?
end

function vector_rep_crystal_type_c(n)
    B = collect(-n:1:n)
    e = function(i, j)
        if i == n
            if -j == n
                return n
            else
                return 0
            end
        else
            if i + 1 == j || i == -j
                return j-1
            else
                return 0
            end
        end
    end
    f = function(i, j)
        if i == n
            if j == n
                return -n
            else
                return 0
            end
        else
            if i == j || i + 1 == -j
                return j+1
            else
                return 0
            end
        end
    end
    eps = function(i, j)
        return (i == -j == n || i+1 == j || i == -j) ? 1 : 0
    end
    phi = function(i, j)
        return (i == j == n || i == j || i + 1 == -j) ? 1 : 0
    end

    return CrystalBase(B, 2*n, n, e, f, eps, phi, missing)
end

function i_sgn(v, crystal, fn_index)
    # - ≡ 0, + ≡ 1
    sign_vec = Vector{Int64}()
    indicator_vec = Vector{Int64}() # tells to which tensor index each sign belongs
    for (i,v_i) in enumerate(v)
        eps_val = crystal.eps(fn_index, v_i)
        phi_val = crystal.phi(fn_index, v_i)
        append!(sign_vec, repeat([0], eps_val))
        append!(sign_vec, repeat([1], phi_val))

        append!(indicator_vec, repeat([i], eps_val + phi_val))
    end

    for j in length(sign_vec)-1:-1:1
        if sign_vec[j] == 1 && length(sign_vec) >= j+1 && sign_vec[j+1] == 0
            deleteat!(sign_vec, [j, j+1])
            deleteat!(indicator_vec, [j, j+1])
        end
    end

    return sign_vec, indicator_vec
end

function tensor_product_helper(v, crystal, fn_index, is_f::Bool)
    res = copy(v)
    sgn, indicator = i_sgn(v, crystal, fn_index)
    if is_f
        first_plus_index = findfirst(x -> x == 1, sgn)
        if first_plus_index === nothing
            return repeat([0], length(v))
        end
        tensor_index = indicator[first_plus_index]
        res[tensor_index] = crystal.f(fn_index, res[tensor_index])
    else # is_e
        last_minus_index = findlast(x -> x == 0, sgn)
        if last_minus_index === nothing
            return repeat([0], length(v))
        end
        tensor_index = indicator[last_minus_index]
        res[tensor_index] = crystal.e(fn_index, res[tensor_index])
    end
    return res
end

# TODO iterative solution
function tensor_product_helper_phieps(v, crystal, fn_index, is_phi::Bool)
    if length(v) == 1
        # recursion anchor
        if is_phi
            return crystal.phi(fn_index, first(v))
        else
            return crystal.eps(fn_index, first(v))
        end
    else
        #recursion step
        if is_phi
            # cut off the last index and recurse
            v1 = v[begin:end-1]
            v2 = v[end]
            return crystal.phi(fn_index, v2) + max(0, tensor_product_helper_phieps(v1, crystal, fn_index, is_phi) - crystal.eps(fn_index, v2))
        else
            # cut off the first index and recurse
            v1 = v[begin]
            v2 = v[begin+1:end]
            return crystal.eps(fn_index, v1) + max(0, tensor_product_helper_phieps(v2, crystal, fn_index, is_phi) - crystal.phi(fn_index, v1))
        end
    end
end

#= TODO: wollen wir das hier vorher checken bevor wir den Crystal manuell aufbauen?
function crystal_operands(hwm::HighestWeightModule)
    if hwm.rs.type == :A
        shape = as_functional(hwm.weight)
        shape = deleteat!(shape, shape .== 0) # Oscar doesn't want trailing zeroes for the partition
        vrc = get_vector_representation_crystal(hwm.rs)
        return map(far_eastern_reading, semistandard_tableaux(shape, vrc.num_elem)) # todo eigentlich num_elem der vektorrep.
    #elseif hwm.rs.type == :C

    end
end
=#

function _get_maximal_vector(w::Weight)
    shape = as_functional(w)
    shape = deleteat!(shape, shape .== 0)
    col_shape = _transpose_shape(shape)
    mv = Int64[]
    for col_len in reverse(col_shape)
        append!(mv, collect(1:col_len)) # c.f. far-eastern reading
    end

    return mv
end

function get_crystal(mod)
    if ismissing(mod.crystal)
        vrc = get_vector_representation_crystal(mod.rs)

        num_ops = vrc.num_ops
        e = function(i,v)
           return tensor_product_helper(v, vrc, i, false) 
        end
        f = function(i,v)
            return tensor_product_helper(v, vrc, i, true) 
        end
        eps = function(i,v)
            return tensor_product_helper_phieps(v, vrc, i, false)
        end
        phi = function(i,v)
            return tensor_product_helper_phieps(v, vrc, i, true)
        end

        # calculate B from the crystal graph
        crystal_graph = build_crystal_graph(_get_maximal_vector(mod.weight), f, num_ops, get_weyl_group_elem(mod))
    
        mod.crystal = CrystalBase(crystal_graph.nodes, length(crystal_graph.nodes), num_ops, e, f, eps, phi, crystal_graph)
    end

    return mod.crystal
end

function get_crystal_graph(hwm::HighestWeightModule)
    c = get_crystal(hwm)

    return c.graph
end

function get_crystal_graph(dmz::DemazureModule)
    c = get_crystal(dmz)

    return c.graph
end




function get_maximal_vector(cb::CrystalBase)
    for b in cb.B
        if all(i -> cb.eps(i, b) == 0, 1:cb.num_ops)
            return b
        end
    end
    @assert false # there must always be a maximal vector!
end


function weight(R::RootSystem, vec)
    # TODO other types
    # Immer auf die ersten epsilons beschränken?
    if R.type == :A
        return Weight(R, [count(x -> x == i, vec) for i in 1:R.rank] .- count(x -> x == R.rank+1, vec), :functionals)
    end
    # return Weight(R, [count(x -> x == i, vec) for i in 1:R.rank] - [count(x -> x == -i, vec) for i in 1:R.rank], :functionals)
end

### Littlewood-Richardson
function maximal_vectors_of_tensor_product(cb1, cb2)
    @assert cb1.num_ops == cb2.num_ops
    b1 = get_maximal_vector(cb1)
    maximal_vectors = []
    phi_b1 = map(i -> cb1.phi(i, b1), 1:cb1.num_ops)
    for b2 in cb2.B
        eps_b2 = map(i -> cb2.eps(i, b2), 1:cb2.num_ops)
        if all(i -> phi_b1[i] >= eps_b2[i], 1:cb1.num_ops)
            push!(maximal_vectors, vcat(b1, b2))
        end
    end

    # cor 4.4.4 S. 104
    # S. 10 simple reflections
    return maximal_vectors
end