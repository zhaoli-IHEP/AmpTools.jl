using BenchmarkTools, Combinatorics, LinearAlgebra, SymEngine, FeynUtils, Dates, Test  

@info "basicTest starts @ $(now())"

@testset "drop_coeff" begin
  @vars shat, x1, x2
  @test drop_coeff( 2*shat*x1 ) == shat*x1
  @test drop_coeff( 2//3*shat*x1 ) == shat*x1
  @test drop_coeff( 2*im*shat*x1 ) == shat*x1
  @test drop_coeff( 2//3*im*shat*x1 ) == shat*x1
  @test drop_coeff( im*shat*x1 ) == shat*x1
end # @testset

@testset "drop_coeff_keep_im" begin
  @vars shat, x1, x2
  @test drop_coeff_keep_im( 2*shat*x1 ) == shat*x1
  @test drop_coeff_keep_im( 2//3*shat*x1 ) == shat*x1
  @test drop_coeff_keep_im( 2*im*shat*x1 ) == im*shat*x1
  @test drop_coeff_keep_im( 2//3*im*shat*x1 ) == im*shat*x1
  @test drop_coeff_keep_im( im*shat*x1 ) == im*shat*x1
end # @testset


@testset "get_Groebner_basis" begin

  @vars x, y
  f = x^3 − 2*x*y 
  g = x^2*y − 2*y^2 + x

  #gb_bench = Basic[(-1/2)*x + y^2, x*y, x^2]
  gb_bench = Basic[y^3, x - 2*y^2]
  gb = get_Groebner_basis_v2( [f, g], [x,y] )
  @test (isempty∘setdiff)( gb, gb_bench )

  #gb_bench = Basic[x^2, x*y, (-1/2)*x + y^2]
  gb_bench = Basic[y^3, x - 2*y^2]
  gb = get_Groebner_basis( [f, g], [x,y] )
  @test (isempty∘setdiff)( gb, gb_bench )

end # @testset


@testset "get_degree" begin
  @vars x, y
  poly = 1+2*x+y*3*x^3
  @test FeynUtils.get_degree( poly, x ) == 3
end # @testset


@testset "make_SP" begin
  @funs SP
  @vars k1, k2, k3, k4, q1, q2
  @test FeynUtils.make_SP(k1,k2) == SP(k1,k2)
  @test FeynUtils.make_SP(k1,-k2) == -SP(k1,k2)
  @test FeynUtils.make_SP(-k1,-k2) == SP(k1,k2)
  @test FeynUtils.make_SP(k1-k3,k2+k4) == SP(k1, k2) + SP(k1, k4) - SP(k2, k3) - SP(k3, k4)

  @test FeynUtils.make_SP( k1, k2 ) == SP(k1, k2)
  @test FeynUtils.make_SP( k1, k2+2*k3 ) == SP(k1, k2) + 2*SP(k1, k3)
  test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
  @test FeynUtils.make_SP( test_mom, test_mom ) == 4*SP(k1, k1) + 12*SP(k1, k2) + 4*SP(k1, k3) + 20*SP(k1, k4) + 24*SP(k1, q1) + 28*SP(k1, q2) + 9*SP(k2, k2) + 6*SP(k2, k3) + 30*SP(k2, k4) + 36*SP(k2, q1) + 42*SP(k2, q2) + SP(k3, k3) + 10*SP(k3, k4) + 12*SP(k3, q1) + 14*SP(k3, q2) + 25*SP(k4, k4) + 60*SP(k4, q1) + 70*SP(k4, q2) + 36*SP(q1, q1) + 84*SP(q1, q2) + 49*SP(q2, q2)
end # @testset





@testset "get_det_naive" begin
  for dim ∈ 1:8
    if dim == 1
      tmp_mat = Matrix{Basic}(undef, 1, 1)
      tmp_mat[1, 1] = Basic("x_1_1")
    else
      tmp_mat = reduce(hcat, [Basic("x_$(ii)_$(jj)") for ii ∈ 1:dim] for jj ∈ 1:dim)
    end

    diff_result = FeynUtils.get_det_naive(tmp_mat)
    for perm ∈ permutations(1:dim)
      σ = (iseven ∘ parity)(perm) ? 1 : -1
      diff_result -= σ * reduce(*, tmp_mat[ii, perm[ii]] for ii ∈ 1:dim)
    end
    @test (iszero ∘ expand)(diff_result)
  end

  for dim ∈ 1:10
    tmp_mat = rand(Basic.(0:100), dim, dim)
    @test det(tmp_mat) == FeynUtils.get_det_naive(tmp_mat)
  end
end


@testset "get_det" begin
  for dim ∈ 1:8
    if dim == 1
      tmp_mat = Matrix{Basic}(undef, 1, 1)
      tmp_mat[1, 1] = Basic("x_1_1")
    else
      tmp_mat = reduce(hcat, [Basic("x_$(ii)_$(jj)") for ii ∈ 1:dim] for jj ∈ 1:dim)
    end

    diff_result = FeynUtils.get_det(tmp_mat)
    for perm ∈ permutations(1:dim)
      σ = (iseven ∘ parity)(perm) ? 1 : -1
      diff_result -= σ * reduce(*, tmp_mat[ii, perm[ii]] for ii ∈ 1:dim)
    end
    @test (iszero ∘ expand)(diff_result)
  end

  for dim ∈ 1:10
    tmp_mat = rand(Basic.(0:100), dim, dim)
    @test det(tmp_mat) == FeynUtils.get_det(tmp_mat)
  end
end


@testset "`get_det_naive` vs. `get_det`" begin
  for dim ∈ 1:10
    tmp_mat = rand(Basic.(-100:100), dim, dim)
    println("For $dim dimensional numeric matrix.")
    println("Original `get_det_naive`:")
    @btime FeynUtils.get_det_naive($tmp_mat)
    println("Bisection `get_det`:")
    @btime FeynUtils.get_det($tmp_mat)
    println()
  end

  for dim ∈ 1:8
    tmp_mat = if dim == 1
      mat = Matrix{Basic}(undef, 1, 1)
      mat[1, 1] = Basic("x_1_1")
      mat
    else
      reduce(hcat, [Basic("x_$(ii)_$(jj)") for ii ∈ 1:dim] for jj ∈ 1:dim)
    end
    println("For $dim dimensional symbolic matrix.")
    println("Original `get_det_naive`:")
    @btime FeynUtils.get_det_naive($tmp_mat)
    println("Bisection `get_det`:")
    @btime FeynUtils.get_det($tmp_mat)
    println()
  end
  @test true
end


@testset "get_adj" begin
  for dim ∈ 1:10
    tmp_mat = rand(Basic.(0:100), dim, dim)
    det_tmp_mat = det(tmp_mat)
    while iszero(det_tmp_mat)
      tmp_mat = rand(Basic.(0:100), dim, dim)
      det_tmp_mat = det(tmp_mat)
    end

    inv_tmp_mat = inv(tmp_mat)
    @test get_adj(tmp_mat) == det_tmp_mat * inv_tmp_mat
  end
end



@testset "gen_sorted_str" begin
  @vars k1, k2, k3, k4, q1, q2
  test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
  expr = expand(test_mom^3+test_mom)
  @test FeynUtils.gen_sorted_str(expr) == ":Add( :Mul( 108,:Pow( q1, 2 ),k3 ),:Mul( 108,k2,k3,q1 ),:Mul( 12,:Pow( k1, 2 ),k3 ),:Mul( 125,:Pow( k4, 3 ) ),:Mul( 126,k2,k3,q2 ),:Mul( 1260,k4,q1,q2 ),:Mul( 135,:Pow( k2, 2 ),k4 ),:Mul( 147,:Pow( q2, 2 ),k3 ),:Mul( 15,:Pow( k3, 2 ),k4 ),:Mul( 150,:Pow( k4, 2 ),k1 ),:Mul( 162,:Pow( k2, 2 ),q1 ),:Mul( 18,:Pow( k3, 2 ),q1 ),:Mul( 180,k1,k2,k4 ),:Mul( 180,k3,k4,q1 ),:Mul( 189,:Pow( k2, 2 ),q2 ),:Mul( 2,k1 ),:Mul( 21,:Pow( k3, 2 ),q2 ),:Mul( 210,k3,k4,q2 ),:Mul( 216,:Pow( q1, 2 ),k1 ),:Mul( 216,:Pow( q1, 3 ) ),:Mul( 216,k1,k2,q1 ),:Mul( 225,:Pow( k4, 2 ),k2 ),:Mul( 252,k1,k2,q2 ),:Mul( 252,k3,q1,q2 ),:Mul( 27,:Pow( k2, 2 ),k3 ),:Mul( 27,:Pow( k2, 3 ) ),:Mul( 294,:Pow( q2, 2 ),k1 ),:Mul( 3,k2 ),:Mul( 324,:Pow( q1, 2 ),k2 ),:Mul( 343,:Pow( q2, 3 ) ),:Mul( 36,:Pow( k1, 2 ),k2 ),:Mul( 36,k1,k2,k3 ),:Mul( 360,k1,k4,q1 ),:Mul( 420,k1,k4,q2 ),:Mul( 441,:Pow( q2, 2 ),k2 ),:Mul( 450,:Pow( k4, 2 ),q1 ),:Mul( 5,k4 ),:Mul( 504,k1,q1,q2 ),:Mul( 525,:Pow( k4, 2 ),q2 ),:Mul( 54,:Pow( k2, 2 ),k1 ),:Mul( 540,:Pow( q1, 2 ),k4 ),:Mul( 540,k2,k4,q1 ),:Mul( 6,:Pow( k3, 2 ),k1 ),:Mul( 6,q1 ),:Mul( 60,:Pow( k1, 2 ),k4 ),:Mul( 60,k1,k3,k4 ),:Mul( 630,k2,k4,q2 ),:Mul( 7,q2 ),:Mul( 72,:Pow( k1, 2 ),q1 ),:Mul( 72,k1,k3,q1 ),:Mul( 735,:Pow( q2, 2 ),k4 ),:Mul( 75,:Pow( k4, 2 ),k3 ),:Mul( 756,:Pow( q1, 2 ),q2 ),:Mul( 756,k2,q1,q2 ),:Mul( 8,:Pow( k1, 3 ) ),:Mul( 84,:Pow( k1, 2 ),q2 ),:Mul( 84,k1,k3,q2 ),:Mul( 882,:Pow( q2, 2 ),q1 ),:Mul( 9,:Pow( k3, 2 ),k2 ),:Mul( 90,k2,k3,k4 ),:Pow( k3, 3 ),k3 )" 

end # @testset



@testset "gen_SPcombo" begin
  @vars q1,q2,q3,k1,k2,k3
  @funs SP
  @test gen_SPcombo( [2,2,0], [k1,k2] ) == sort( [ 
          make_SP(q2, k1)^2*make_SP(q1, k2)*make_SP(q1, k1), 
          make_SP(q2, k2)^2*make_SP(q1, k1)^2, 
          make_SP(q2, k2)^2*make_SP(q1, k2)*make_SP(q1, k1), 
          make_SP(q2, k2)^2*make_SP(q1, k2)^2, 
          make_SP(q2, k1)*make_SP(q2, k2)*make_SP(q1, k2)^2, 
          make_SP(q2, q2)*make_SP(q1, k2)*make_SP(q1, k1), 
          make_SP(q2, k1)*make_SP(q2, k2)*make_SP(q1, k1)^2, 
          make_SP(q2, q2)*make_SP(q1, k1)^2, 
          make_SP(q2, k1)^2*make_SP(q1, k2)^2, 
          make_SP(q2, k1)*make_SP(q1, k2)*make_SP(q1, q2), 
          make_SP(q2, k2)^2*make_SP(q1, q1), 
          make_SP(q2, q2)*make_SP(q1, q1), 
          make_SP(q2, k1)^2*make_SP(q1, q1), 
          make_SP(q2, k1)*make_SP(q2, k2)*make_SP(q1, k2)*make_SP(q1, k1), 
          make_SP(q2, q2)*make_SP(q1, k2)^2, 
          make_SP(q1, q2)^2, 
          make_SP(q2, k2)*make_SP(q1, k1)*make_SP(q1, q2), 
          make_SP(q2, k1)*make_SP(q1, k1)*make_SP(q1, q2), 
          make_SP(q2, k1)^2*make_SP(q1, k1)^2, 
          make_SP(q2, k1)*make_SP(q2, k2)*make_SP(q1, q1), 
          make_SP(q2, k2)*make_SP(q1, k2)*make_SP(q1, q2) 
        ], by=gen_sorted_str )

end # @testset


@info "basicTest ends @ $(now())"

