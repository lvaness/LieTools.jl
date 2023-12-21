using LieTools
using LinearAlgebra
using Oscar

# We want to compute Gröbner bases for highest weight modules over Lie algebras of type A and C
R1 = LieTools.RootSystem(:C, 3)

# The monomial ordering is given by https://arxiv.org/abs/1010.2321 Defn. 2.6
# This ordering can be realized as a matrix ordering (see https://docs.oscar-system.org/stable/CommutativeAlgebra/GroebnerBases/orderings/#Matrix-Orderings)
# Each column references a particular positive root, in descending order based on height (the exact order is currently an implementation detail)
const c3_fflv_top_left =
[1  1   1   1   1   1   1   1   1;
 0  0   0   0   0   0   -1  0   0;
 0  0   -1  0   -1  0   0   -1  0;
 -1 -1  0   -1   0  -1  0   0  -1;
 0  0   0   0   0   0   1   0   0;
 0  0   1   0   0   0   0   0   0;
 0  0   0   0   1   0   0   0   0;
 0  0   0   0   0   0   0   1   0;
 1  0   0   0   0   0   0   0   0;
 0  1   0   0   0   0   0   0   0;
 0  0   0   1   0   0   0   0   0;
 0  0   0   0   0   1   0   0   0;
 0  0   0   0   0   0   0   0   1]
# All other root vectors are sorted arbitrarily (taking a grlex ordering)
const c3_fflv_top_right = [transpose(repeat([1], 12)); zeros(Int64, 12, 12)]
c3_fflv_ordering_matrix = [c3_fflv_top_left c3_fflv_top_right; zeros(Int64, 12,9) I(12)]

LieTools.set_uea_ordering_matrix(R1, c3_fflv_ordering_matrix)

M1 = LieTools.HighestWeightModule(R1, LieTools.Weight(R1, [1,1,0], :fundamental_weights))

# Show the polytope associated to the basis of the module. This should be the FFLV polytope
println("C3 - Weight ω_1 + ω_2\nPolytope:")
display(facets(LieTools.get_polytope(M1)))

# Show the Gröbner basis associated to the ideal
println("Gröbner basis:")
display(LieTools.get_gb(M1))

# Type A
R2 = LieTools.RootSystem(:A, 3)

const a3_fflv_top_left =
[1  1   1   1   1   1;
 0  0   0   1   0   0;
 0  1   0   0   0   0;
 0  0   0   0   1   0;
 1  0   0   0   0   0;
 0  0   1   0   0   0;
 0  0   0   0   0   1];
const a3_fflv_top_right = [transpose(repeat([1], 9)); zeros(Int64, 6, 9)];
a3_fflv_ordering_matrix = [a3_fflv_top_left a3_fflv_top_right; zeros(Int64, 9,6) I(9)];

LieTools.set_uea_ordering_matrix(R2, a3_fflv_ordering_matrix)
M2 = LieTools.HighestWeightModule(R2, LieTools.Weight(R2, [1,1,0], :fundamental_weights))

println("\nA3 - Weight ω_1 + ω_2\nPolytope:")
display(facets(LieTools.get_polytope(M2)))
println("Gröbner basis:")
display(LieTools.get_gb(M2))