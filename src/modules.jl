using DataStructures

function HighestWeightModule(rs::RootSystem, w::Weight)
    return HighestWeightModule(rs, w, missing, missing, missing, missing, missing, missing)
end

function get_weight(hwm::HighestWeightModule)
    return hwm.weight
end

function DemazureModule(rs::RootSystem, w::Weight, g::WeylGroupElement)
    gens = get_weyl_group(rs).simple_refls
    nf_elem = isempty(g) ? one(get_weyl_group(rs).cox_grp) : short_lex(prod(map(i -> gens[i], g)))
    return DemazureModule(rs, w, nf_elem, missing, missing, missing, missing)
end

function get_weight(dmz::DemazureModule)
    return dmz.weight
end

#
# Characters
#

function _weyl_character_formula(rs::RootSystem, w::Weight)
    weyl_elems = get_weyl_elements(rs)
    parities = map(w -> length(w) % 2 == 0 ? 1 : -1, weyl_elems)

    rho = Weight(rs, repeat([1], rs.rank), :fundamental_weights)
    w_plus_rho = w + rho

    Xs = gens(get_character_ring(rs))

    # TODO build context?
    numerator = sum([p * as_character_ring_elem(apply_weyl_group_element(r, w_plus_rho)) for (p, r) in zip(parities, weyl_elems)])
    denominator = sum([p * as_character_ring_elem(apply_weyl_group_element(r, rho)) for (p, r) in zip(parities, weyl_elems)])

    return numerator/denominator
end

function get_character(hwm::HighestWeightModule, method=:weyl_char; cache=true)
    character = hwm.character

    if ismissing(character) || !cache
        if method == :weyl_char
            character = _weyl_character_formula(hwm.rs, hwm.weight)
        elseif method == :littelmann
            character = _littelmann_demazure_formula(hwm.rs, get_crystal(hwm))
        end
    end

    if cache
        hwm.character = character
    end        

    return character
end

#
# demazure characters
#

function _demazure_operator(root, polynomial)
    ans = zero(parent(polynomial))
    for (coeff, weight_fw) in zip(coefficients(polynomial), exponent_vectors(polynomial))
        w = Weight(root.rs, weight_fw, :fundamental_weights)
        chain_len = lf(w, root)

        if chain_len > -1
            for i in 0:chain_len
                ans += coeff * as_character_ring_elem(w - i * root)
            end
        elseif chain_len < -1
            for i in 0:chain_len
                ans -= coeff * as_character_ring_elem(w + (1 - i) * root)
            end
        end
    end
    return ans
end

function _demazure_character_formula(w, weyl_group_element)
    res = as_character_ring_elem(w)
    for refl_index in reverse(weyl_group_element)
        res = _demazure_operator(Weight(w.rs, [i == refl_index ? 1 : 0 for i in 1:w.rs.rank], :roots), res)
    end
    return res
end

# Use Littelmann's refined Demazure character formula
# cf. https://doi.org/10.1215/S0012-7094-93-07131-1
function _littelmann_demazure_formula(rs::RootSystem, c::CrystalBase)

    return sum(as_character_ring_elem(weight(rs, b)) for b in c.B)
end

function get_character(dmz::DemazureModule, method=:littelmann; use_cached=false)
    if ismissing(dmz.character) || !use_cached
        if method == :littelmann
            dmz.character = _littelmann_demazure_formula(dmz.rs, get_crystal(dmz))
        elseif method == :demazure
            dmz.character = _demazure_character_formula(dmz.weight, dmz.weyl_group_element)
        end
    end

    return dmz.character
end

#
# dimensions
#

function _weyl_dimension_formula(rs::RootSystem, w::Weight)
    root_system_gap = get_gap_root_system(rs)
    n = rs.rank
    # these are in the basis of fundamental weights!
    pos_roots_gap = GAP.Globals.PositiveRoots(root_system_gap)
    pos_roots_fw = GAP.gap_to_julia(Vector{Vector{Int64}}, pos_roots_gap)
    pos_roots = map(w -> Weight(rs, w, :fundamental_weights), pos_roots_fw)

    rho = Weight(rs, repeat([1], rs.rank), :fundamental_weights)

    lambda_plus_rho = w + rho

    numerator = prod([lf(lambda_plus_rho, a) for a in pos_roots])
    denominator = prod([lf(rho, a) for a in pos_roots])

    return Int64(numerator//denominator)
end

function get_dimension(hwm::HighestWeightModule, method=:weyl_dim)
    if ismissing(hwm.dimension)
        if method == :weyl_dim
            return _weyl_dimension_formula(hwm.rs, hwm.weight)
        end
        # TODO: weyl_char, crystal
    end
end
    
function tensor_product(mod1, mod2, method=:crystal)
    @assert mod1.rs == mod2.rs
    if method == :crystal
        return TensorProductWrapper(mod1.rs, maximal_vectors_of_tensor_product(get_crystal(mod1), get_crystal(mod2)), sum(as_functional(mod1.weight)))
    end
end

function show_weights(t::TensorProductWrapper, display_type=:fundamental_weights)
    c = counter(map(v ->convert(weight(t.rs, v), display_type), t.maximal_vectors))
    return c
end

function tail_of_maximal_vectors(t::TensorProductWrapper)
    return map(x -> x[t.first_component_length+1:end], t.maximal_vectors)
end

function degenerate_module(hwm::HighestWeightModule)
    #@assert hwm.rs.type == :A
    rank = hwm.rs.rank
    new_rank = 2 * rank - 1
    new_rs = RootSystem(hwm.rs.type, new_rank)
    
    weight_fw = as_fundamental_weight(hwm.weight)
    new_weight_fw = [ (i % 2 == 0) ? 0 : weight_fw[(i ÷ 2) + 1] for i in 1:new_rank]
    new_weight = Weight(new_rs, new_weight_fw, :fundamental_weights)

    new_weyl_group_elem = reduce(vcat, [collect(0:k-1) .+ k for k in rank:-1:1])
    return DemazureModule(new_rs, new_weight, new_weyl_group_elem)
end

function get_weyl_group_elem(dmz::DemazureModule)
    return dmz.weyl_group_element
end

function get_weyl_group_elem(hwm::HighestWeightModule)
    return short_lex(get_longest_weyl_elem(hwm.rs))
end