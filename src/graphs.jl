
function build_crystal_graph(starting_node, f, num_ops, dmz_edge_restriction #= longest weyl perm =#)
    image_of_current_refl = Vector{Vector{Int64}}() # after applying the highest necessary exponent of the current reflection, this is what we will apply the next one to
    not_eval_fully_yet_idxs = Set{Int64}() # fill this with the image of each fn application, once this is empty, we can move on to the next one

    g = Graph{Directed}(0)
    incidence_tuples = Set{Tuple{Int64, Int64, Int64}}()
    
    push!(image_of_current_refl, starting_node) # the maximal vector is where we start
    add_vertex!(g)

    for simple_refl_idx in reverse(dmz_edge_restriction) # since we start from the maximal vector, invert the permutation
        not_eval_fully_yet_idxs = Set(1:length(image_of_current_refl)) # gather everything we've reached before any use this to apply to

        while !isempty(not_eval_fully_yet_idxs)
            current_vector_idx = pop!(not_eval_fully_yet_idxs)
            img = f(simple_refl_idx, image_of_current_refl[current_vector_idx])
            if !iszero(img)
                img_idx = findfirst(x -> x == img, image_of_current_refl)
                if isnothing(img_idx)
                    push!(image_of_current_refl, img)
                    add_vertex!(g)
                    img_idx = length(image_of_current_refl)
                end
                push!(not_eval_fully_yet_idxs, img_idx)
                push!(incidence_tuples, (current_vector_idx, img_idx, simple_refl_idx))
                add_edge!(g, current_vector_idx, img_idx)
            end
        end
    end

    # now we have all nodes, but we might be missing some edges
    for (node_idx, node) in enumerate(image_of_current_refl)
        for op_idx in 1:num_ops
            img = f(op_idx, node)
            if iszero(img) continue end
            img_idx = findfirst(x -> x == img, image_of_current_refl)
            if isnothing(img_idx) continue end

            add_edge!(g, node_idx, img_idx)
            push!(incidence_tuples, (node_idx, img_idx, op_idx))
        end
    end


    incidence_tuples = collect(incidence_tuples)

    incidence_i = getindex.(incidence_tuples, 1)
    incidence_j = getindex.(incidence_tuples, 2)
    incidence_v = getindex.(incidence_tuples, 3)

    incidence_mat = sparse(incidence_i, incidence_j, incidence_v)

    return CrystalGraph(g, image_of_current_refl, incidence_mat)
end

function show_crystal_graph(mod, as_tableaux=false)
    cg = get_crystal_graph(mod)
    crystal = get_crystal(mod)
    edge_colors = Colors.distinguishable_colors(crystal.num_ops, Colors.colorant"red")
    grg = GR.SimpleDiGraph(map(x-> GR.Edge((Oscar.src(x), Oscar.dst(x))), Oscar.edges(cg.graph)))
    vert_labels = nothing
    if as_tableaux
        shape = as_functional(get_weight(mod)) # the weight of the maximal vector is the shape of the tableaux
        vert_labels = _print_nice_tableau.(map(b -> tableau_from_far_eastern(shape, b), cg.nodes))
    else
        vert_labels = string.(cg.nodes)
    end
    vert_colors = [:black for _ in 1:nv(cg.graph)]
    vert_colors[1] = :red
    f, ax, p = graphplot(grg, nlabels=vert_labels, elabels=collect(map(x->string(cg.incidence[GR.src(x), GR.dst(x)]),GR.edges(grg))), edge_color=collect(map(x->edge_colors[cg.incidence[GR.src(x), GR.dst(x)]], GR.edges(grg))), node_color=vert_colors, nlabels_attr=(;font = "Mono"))
    Makie.hidespines!(ax)
    Makie.hidedecorations!(ax)
    return f
end