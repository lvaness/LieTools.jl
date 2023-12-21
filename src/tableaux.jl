
"""
    middle_eastern_reading(tableau)

Compute the middle eastern reading of a (Young) tableau.

This reading returns a vector iterating the values of the tableau,
moving algong the rows right-to-left from top to bottom.

# Examples
```julia-repl
julia> tab=Tableau([[1,2,3],[4,5],[6]])
[[1, 2, 3], [4, 5], [6]]

julia> middle_eastern_reading(tab)
6-element Vector{Int64}:
 3
 2
 1
 5
 4
 6
```
"""
function middle_eastern_reading(y::Tableau)
    v = Int64[]
    for (row,row_length) in enumerate(shape(y)) 
        for col in row_length:-1:1
            push!(v, y[row, col])
        end
    end
    return v
end

function tableau_from_middle_eastern(shape, v)
    new_rows = Vector{Vector{Int64}}()
    # now flip row-wise
    current_row_start = 1
    for row_length in shape
        push!(new_rows, reverse(v[current_row_start:current_row_start + row_length - 1]))
        current_row_start += row_length
    end

    return Tableau(new_rows)
end

function _transpose_shape(row_shape)
    return [findlast(j -> j >= i, row_shape) for i in 1:maximum(collect(row_shape))]
end

function far_eastern_reading(y::Tableau)
    v = Int64[]
    column_shape = _transpose_shape(shape(y))
    for (col, col_length) in reverse(collect(enumerate(column_shape)))
        for row in 1:col_length
            push!(v, y[row, col])
        end
    end
    return v
end

function tableau_from_far_eastern(shape, v)
    column_shape = _transpose_shape(shape)
    mat = zeros(Int64, column_shape[1], shape[1])
    i = 1
    for (col, col_length) in reverse(collect(enumerate(column_shape)))
        for row in 1:col_length
            mat[row, col] = v[i]
            i += 1
        end
    end
    vecs = [mat[row, 1:shape[row]] for row in 1:column_shape[1]]

    return Tableau(vecs)
end


tikz_colors= ["red", "green", "blue", "cyan", "magenta", "yellow", "black", "gray", "darkgray", "lightgray", "brown", "lime", "olive", "orange", "pink", "purple", "teal", "violet"]

# requires tikz-cd, ytableau
#=
function get_latex_string(t::Tableau)
    s = raw"\begin{ytableau}"
    s *= "\n"
    for row in t
        row_len = length(row)
        for i in 1:row_len
            s *= string(row[i], " ")
            if i == row_len
                s *= string(raw"\\ ", "\n")
            else
                s *= "& "
            end
        end
    end
    s *= raw"\end{ytableau}"

    return s
end
=#

function _print_nice_tableau(y::Tableau)
    str = ""
    for line in y
        str *= "["
        for (i,number) in enumerate(line)
            str *= "$number"
            if i != lastindex(line)
                str *= ","
            end
        end
        str *= "]\n"
    end
    return str
end