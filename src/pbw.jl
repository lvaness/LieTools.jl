
"""
Credit: Thomas Breuer
"""
function _sc_table_info(rs)
    type = string(rs.type)
    n = rs.rank

    A = GAP.Globals.SimpleLieAlgebra(GAP.GapObj(type), n, GAP.Globals.Rationals)::GAP.GapObj
    B = GAP.Globals.CanonicalBasis(A)::GAP.GapObj
    T = GAP.Globals.StructureConstantsTable(B)::GAP.GapObj
    T = GAP.Globals.ELMS_LIST(T, GAP.GapObj(1:length(B)))::GAP.GapObj
    sct = Matrix{Tuple{Vector{Int}, Vector{Int}}}(T)
    roots = Vector(map(l -> Vector{Int}(GAP.Globals.List(l, x -> GAP.Globals.Position(B, x))),
                       GAP.Globals.ChevalleyBasis(A)))
    return roots, sct
end

function get_num_pos_roots(rs)
    v = gens(get_uea(rs))
    return (length(v) - rs.rank) ÷ 2
end

function get_uea_vars(rs)
    v = gens(get_uea(rs))
    n_simple = rs.rank
    n_pos = get_num_pos_roots(rs)
    f = v[1:n_pos]
    h = v[(1:n_simple) .+ n_pos]
    e = v[(1:n_pos) .+ (n_simple + n_pos)]

    return f,h,e
end

function set_uea_ordering_matrix(rs, ordering)
    rs.uea_ordering = ordering

    if !ismissing(rs.uea)
        _regenerate_uea(rs)
    end
end

function get_uea(rs)
    if ismissing(rs.uea)
        _regenerate_uea(rs)
    end

    return rs.uea
end

function _regenerate_uea(rs)
    _, sct = _sc_table_info(rs)

    root_system = get_gap_root_system(rs)
    simple_roots = GAP.gap_to_julia(Vector{Vector{Int64}}, GAP.Globals.SimpleSystem(root_system))
    pos_roots = (GAP.gap_to_julia(Vector{Vector{Int64}}, GAP.Globals.PositiveRoots(root_system)))
    n_simple = length(simple_roots)
    n_pos = length(pos_roots)

    # Translate between indexing
    # GAP: pos - neg - h
    # desired: neg - h - pos
    function desired_to_gap(i)
        n_roots = n_pos
        n_h = n_simple
        if i ≤ n_roots # we want negative
            #return i + n_roots
            return 2*n_roots - i + 1 # ordering_adjusted
        elseif i <= n_roots + n_h # we want h
            return i + n_roots
        else # we want pos
            return i - n_roots - n_h
        end
    end

    function gap_to_desired(i)
        n_roots = n_pos
        n_h = n_simple
        if i ≤ n_roots # we want positive
            return i + n_roots + n_h
        elseif i <= 2*n_roots # we want neg
            #return i - n_roots
            return n_roots - (i - n_roots) + 1 
        else # we want h
            return i - n_roots
        end
    end

    # Set up the base symmetric algebra over the Lie algebra
    n⁻ = "f_".*string.(1:n_pos)
    h = "h_".*string.(1:n_simple)
    n⁺ = "e_".*string.(1:n_pos)
    SA, v = polynomial_ring(QQ, [n⁻; h; n⁺]; ordering = :lex)

    # Get the structure constants obtained from GAP and transform them to be usable by Oscar
    C = [v[i] * v[j] for i = 1:length(v), j = 1:length(v)]
    D = [ 
            sct[desired_to_gap(x),desired_to_gap(y)][1] == [] ? 0 : 
                -sum([sct[desired_to_gap(x),desired_to_gap(y)][2][k] * v[gap_to_desired(sct[desired_to_gap(x),desired_to_gap(y)][1][k])]
                for k = 1:length(sct[desired_to_gap(x),desired_to_gap(y)][1])]) 
                for x = 1:length(v), y = 1:length(v)
            
        ] 

    # Set up the universal enveloping algebra
    UA, v = pbw_algebra(SA, C+D, ismissing(rs.uea_ordering) ? :grlex : matrix_ordering(v, rs.uea_ordering))

    rs.uea = UA
end

function get_ideal(hwm)
    if ismissing(hwm.ideal)
        # Set up the ideal that we need to "mod out" for the standard cyclic module
        maximal_weight = as_fundamental_weight(hwm.weight)
        f, h, e = get_uea_vars(hwm.rs)
        gens_1 = e;
        gens_2 = [ h[i] - maximal_weight[i] for i=1:hwm.rs.rank ]
        gens_3 = [ f[-i+1 + length(e)]^(maximal_weight[i]+1) for i in 1:hwm.rs.rank ] # TODO This is jank, we need to make sure the f[i] correspond to h[i] above to get the proper module, we need to fix e to get the expected operation, reverting the order here is a quick fix

        hwm.ideal = left_ideal(get_uea(hwm.rs), [ gens_1 ; gens_2; gens_3 ])
    end

    return hwm.ideal
end

function get_polytope(hwm)
    # TODO: implementation detail
    exps = [Singular.leading_exponent_vector(p)[1:get_num_pos_roots(hwm.rs)] for p in Singular.gens(Singular.kbase(Singular.std(get_ideal(hwm).sdata)))]

    # Construct the associated polytope. This should be the FFLV polytope for types A and C
    return convex_hull(mapreduce(permutedims, vcat, exps))
end

function get_gb(hwm)
    # TODO: implementation detail
    return Singular.gens(Singular.std(LieTools.get_ideal(hwm).sdata, complete_reduction=true))
end