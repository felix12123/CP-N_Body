using Test

function run_tests()
  function gravity(body, bodys)
    G = 1
    F = (0, 0, 0)
    for b in bodys
      F = F .+ (G * b[3] * body[3] / sqrt(sum((b[1] .- body[1]) .^ 2))^3) .* (b[1] .- body[2])
    end 
    return F
  end
  bs = [[(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0], [(0.0, 0.0, 0.0), (0.0, -1.0, 0.0), 1.0]]
  b = [(-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1]

  @testset "integ_euler" begin
    
    @test isa(integ_euler(b, x -> gravity(x, bs), 0.01), Vector)
  end
  
  @testset "parser" begin
    @test parse_data("data\\2_body.csv") == [[(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0], [(0.0, 0.0, 0.0), (0.0, -1.0, 0.0), 1.0]]
  end
end


# if string(@__MODULE__) === "Main"
#   run_tests()
# end