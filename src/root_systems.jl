function RootSystem(type::Symbol, rank::Int64)
    valid = false
    if type in [:A, :B, :C, :D]
        valid = rank > 0 # TODO
    elseif type == :E
        valid = rank in [6,7,8]
    elseif type == :F
        valid = rank == 4
    elseif type == :G
        valid = rank == 2
    end

    if !valid
        throw("Invalid Root System type")
    end

    return RootSystem(type, rank, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
end

function get_longest_weyl_elem(rs::RootSystem)
    # TODO different per type?
    weyl_grp = get_weyl_group(rs)
    return longest_element(weyl_grp.cox_grp)
end

# Iterate the group by acting on identity by right multiplication
function get_weyl_elements(rs::RootSystem)
    g = get_weyl_group(rs)
    if ismissing(g.elements)
        orbit = [one(g.cox_grp)]
        j = 1
        l = 1
        while j <= l
            current_elem = orbit[j]
            for refl in g.simple_refls
                prod = current_elem * refl
                if prod ∉ orbit
                    push!(orbit, prod)
                    l += 1
                end
            end
            j += 1
        end
        g.elements = orbit
    end
    return g.elements
end

function get_vector_representation_crystal(rs::RootSystem)
    if ismissing(rs.vector_rep_crystal)
        if rs.type == :A
            rs.vector_rep_crystal = vector_rep_crystal_type_a(rs.rank)
        elseif rs.type == :C
            rs.vector_rep_crystal = vector_rep_crystal_type_c(rs.rank)
        end
    end

    return rs.vector_rep_crystal
end

function get_gap_lalg(rs::RootSystem)
    if ismissing(rs.gap_lalg)
        rs.gap_lalg = GAP.Globals.SimpleLieAlgebra(GAP.GapObj(string(rs.type)), rs.rank, GAP.Globals.Rationals)
    end

    return rs.gap_lalg
end

function get_gap_root_system(rs::RootSystem)
    if ismissing(rs.gap_rs)
        lalg = get_gap_lalg(rs)
        rs.gap_rs = GAP.Globals.RootSystem(lalg)
    end

    return rs.gap_rs
end

function get_cartan_matrix(rs::RootSystem)
    if ismissing(rs.cartan_matrix)
        gap_rs = get_gap_root_system(rs)
        rs.cartan_matrix = Rational.(transpose(GAP.gap_to_julia(Matrix{Int64}, GAP.Globals.CartanMatrix(gap_rs))))
    end

    return rs.cartan_matrix
end

function get_cartan_matrix_inv(rs::RootSystem)
    if ismissing(rs.cartan_matrix_inv)
        cartan_matrix = get_cartan_matrix(rs)
        rs.cartan_matrix_inv = inv(cartan_matrix)
    end

    return rs.cartan_matrix_inv
end

function get_functional_to_root_matrix(rs::RootSystem)
    if ismissing(rs.functional_to_root_matrix)
        n = rs.rank
        if rs.type == :A
            rs.functional_to_root_matrix = [ j > i ? -i//(n+1) : (n+1-i)//(n+1) for i=1:n, j=1:n ]
        end
    end

    return rs.functional_to_root_matrix
end

function get_root_to_functional_matrix(rs::RootSystem)
    if ismissing(rs.root_to_functional_matrix)
        n = rs.rank
        if rs.type == :A
            rs.root_to_functional_matrix = [ 
                if i == j == n 
                    2
                elseif j == n || i == j 
                    1
                elseif i == j+1 
                    -1
                else 
                    0
                end
            for i=1:n, j=1:n ]
        end
    end

    return rs.root_to_functional_matrix
end

function get_functional_to_fundamental_weight_matrix(rs::RootSystem)
    if ismissing(rs.functional_to_fundamental_weight_matrix)
        n = rs.rank
        if rs.type == :A
            rs.functional_to_fundamental_weight_matrix = [
                if i == j
                    1
                elseif i+1 == j
                    -1
                else
                    0
                end
            for i=1:n, j=1:n]
        end
    end

    return rs.functional_to_fundamental_weight_matrix
end

function get_fundamental_weight_to_functional_matrix(rs::RootSystem)
    if ismissing(rs.fundamental_weight_to_functional_matrix)
        n = rs.rank
        if rs.type == :A
            rs.fundamental_weight_to_functional_matrix = [ i <= j ? 1 : 0 for i=1:n, j=1:n]
        end
    end

    return rs.fundamental_weight_to_functional_matrix
end

function get_positive_roots(rs::RootSystem)
    # TODO
end

function get_weyl_group(rs::RootSystem)
    if ismissing(rs.weyl_group)
        rs.weyl_group = WeylGroup(coxeter_group_min(Int.(get_cartan_matrix(rs)))..., missing)
    end

    return rs.weyl_group
end

function apply_simple_reflection(root_index, w::Weight)
    ith_root = Weight(w.rs, [k == root_index ? 1 : 0 for k = 1:w.rs.rank], :roots)
    return w - lf(w, ith_root) * ith_root
end

function apply_weyl_group_element(weyl_group_elem, w::Weight)
    v = short_lex(weyl_group_elem)
    res = w
    for i in reverse(v)
        res = apply_simple_reflection(i, res)
    end
    return res
end

function get_character_ring(rs::RootSystem)
    if ismissing(rs.character_ring)
        indices = [map(x -> x + 0x2080 - 0x30, string(i)) for i in 1:rs.rank] # 0x30 = 0, 0x2080 = subscript 0

        R, _ = LaurentPolynomialRing(AbstractAlgebra.ZZ, "ω" .* indices)

        rs.character_ring = R
    end

    return rs.character_ring
end

function as_character_ring_elem(w::Weight)
    R = get_character_ring(w.rs)
    Xs = gens(R)

    return prod(Xs .^ Int.(as_fundamental_weight(w)))
end
