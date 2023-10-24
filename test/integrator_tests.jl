using Test

function run_tests()
  @testset "integrators" begin
  end
  
  @testset "parser" begin
    @test parse_data("data\\2_body.csv") == [[(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0, 0], [(0.0, 0.0, 0.0), (0.0, -1.0, 0.0), 1.0, 0]]
  end
end