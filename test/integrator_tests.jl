using Test

function run_tests()
  @testset "integrators" begin
    @test 1==1
    @test 2+2 ≈ 4.000000000000001
  end
  
  @testset "parser" begin
    @test parse_data("data\\2_body.csv") == [[(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0], [(0.0, 0.0, 0.0), (0.0, -1.0, 0.0), 1.0]]
  end
end