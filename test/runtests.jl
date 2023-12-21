using LieTools
using SparseArrays
using Test
using Oscar

@testset "LieTools.jl" begin
    # dominant weights
    # TODO

    @testset "Generalized Cartan Matrices" begin
        generalized_cartan_matrix = [2 -1; -1 2]
        @test LieTools.is_generalized_cartan(generalized_cartan_matrix)
        @test LieTools.is_generalized_cartan(sparse(generalized_cartan_matrix))

        not_generalized_cartan_matrix_non_square = [2 -1 0; -1 2 0]
        @test !LieTools.is_generalized_cartan(not_generalized_cartan_matrix_non_square)
        @test !LieTools.is_generalized_cartan(sparse(not_generalized_cartan_matrix_non_square))

        not_generalized_cartan_matrix_non_integral = [2 -0.5; -0.5 2]
        @test !LieTools.is_generalized_cartan(not_generalized_cartan_matrix_non_integral)
        @test !LieTools.is_generalized_cartan(sparse(not_generalized_cartan_matrix_non_integral))

        not_generalized_cartan_matrix_diagonal_wrong = [1 0; 0 2]
        @test !LieTools.is_generalized_cartan(not_generalized_cartan_matrix_diagonal_wrong)
        @test !LieTools.is_generalized_cartan(sparse(not_generalized_cartan_matrix_diagonal_wrong))

        not_generalized_cartan_matrix_positive_nondiagonal = [2 1; 1 2]
        @test !LieTools.is_generalized_cartan(not_generalized_cartan_matrix_positive_nondiagonal)
        @test !LieTools.is_generalized_cartan(sparse(not_generalized_cartan_matrix_positive_nondiagonal))

        not_generalized_cartan_matrix_asymmetric_zeroes = [2 0; -1 2]
        @test !LieTools.is_generalized_cartan(not_generalized_cartan_matrix_asymmetric_zeroes)
        @test !LieTools.is_generalized_cartan(sparse(not_generalized_cartan_matrix_asymmetric_zeroes))
    end

    @testset "Tableau Reading" begin
        young_tableau = Oscar.Tableau([[1,2,3], [4,5], [6]])
        @test LieTools.tableau_from_middle_eastern(Oscar.shape(young_tableau), LieTools.middle_eastern_reading(young_tableau)) == young_tableau

        reading = [3,2,1,5,4,6]
        @test LieTools.middle_eastern_reading(LieTools.tableau_from_middle_eastern([3,2,1], reading)) == reading
    end
end
