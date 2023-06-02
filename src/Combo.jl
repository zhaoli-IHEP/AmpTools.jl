


##########################
function split_coeff( 
    expr::Basic 
)::Tuple{Basic,Basic}
##########################

  expr_class = SymEngine.get_symengine_class( expr )
  if expr_class != :Mul
    if expr_class ∈ [:Complex,:Integer,:Rational] 
      return expr, one(Basic)
    else
      return one(Basic), expr
    end # if
  else # == :Mul
    factor_list = get_args(expr)
    coeff_expr = one(Basic)
    trimmed_expr = one(Basic)
    for one_factor in factor_list 
      if SymEngine.get_symengine_class(one_factor) ∈ [:Complex,:Integer,:Rational]
        coeff_expr *= one_factor
      else
        trimmed_expr *= one_factor
      end # if
    end # for one_factor
    return coeff_expr, trimmed_expr
  end # if

end # function split_coeff


##########################
function split_coeff( 
    expr_list::Vector{Basic} 
)::Tuple{Vector{Basic},Vector{Basic}}
##########################

  coeff_list = zeros( Basic, length(expr_list) )
  remnant_list = zeros( Basic, length(expr_list) )

  for index in 1:length(expr_list)
    coeff_list[index], remnant_list[index] = split_coeff(expr_list[index])
  end # for index

  return coeff_list, remnant_list

end # function split_coeff

#########################################
function drop_coeff( expr::Basic )::Basic
#########################################

  return is_class( :Mul, expr ) ?
         (prod∘filter)( x_-> SymEngine.get_symengine_class(x_) ∉ [:Complex,:Integer,:Rational], get_args(expr) ) : expr

end # function drop_coeff



#########################################
function drop_coeff_keep_im( expr::Basic )::Basic
#########################################

  if is_class( :Mul, expr ) 
    remove_integer_list = filter( x_-> SymEngine.get_symengine_class(x_) ∉ [:Integer,:Rational], get_args(expr) )
    keep_im_list = map( x_-> SymEngine.get_symengine_class(x_) == :Complex ? im : x_, remove_integer_list )
    return prod(keep_im_list)
  else
    return expr
  end # if

  return expr

end # function drop_coeff_keep_im





#########################################
function generate_SPcombo(
    rank_str::String,
    independent_mom_list::Array{Basic}
)::Array{Basic}
#########################################

  # e.g. generate_SP([2,2],[k1,k2,K3])
  # only for two-loop level

  @vars q1,q2
  @funs SP

  ###### r1 is the rank of q1, r2 is the rank of q2
  r1 = count( "q1", rank_str )
  r2 = count( "q2", rank_str )
  num_mom = length(independent_mom_list)

  # term_q1_k = SP(q1,k1)+SP(q1,k2)+...
  term_q1_k = 0
  for i in 1:num_mom
    term_q1_k += make_SP(q1,independent_mom_list[i])
  end # for

  # term_q2_k = SP(q2,k1)+SP(q2,k2)+...
  term_q2_k = 0
  for i in 1:num_mom
    term_q2_k += make_SP(q2,independent_mom_list[i])
  end # for


  #q1_q2_exp is the exponent of SP(q1,q2)
  #q1_q1_exp is the exponent of SP(q1,q1)
  #q2_q2_exp is the exponent of SP(q2,q2)
  #q1_remain = total q1 rank - q1 rank in SP(q1,q2)
  #q2_remain = total q2 rank - q2 rank in SP(q1,q2)
  #q1_k_exp is the exponent of (SP(q1,k1)+SP(q1,k2)+...)
  #q2_k_exp is the exponent of (SP(q2,k1)+SP(q2,k2)+...)
  #total_term is the summation of all possible terms
  #expand the total_term and drop the coefficients

  total_term = 0
  for q1_q2_exp in 0:min(r1,r2)
    q1_remain = r1-q1_q2_exp
    q2_remain = r2-q1_q2_exp
    max_q1_q1_exp = floor(Int,q1_remain/2)
    max_q2_q2_exp = floor(Int,q2_remain/2)
    for q1_q1_exp in 0:max_q1_q1_exp
      q1_k_exp = q1_remain - 2*q1_q1_exp
      for q2_q2_exp in 0:max_q2_q2_exp
        q2_k_exp = q2_remain - 2*q2_q2_exp
        total_term += term_q1_k^q1_k_exp*term_q2_k^q2_k_exp*
                      make_SP(q1,q1)^q1_q1_exp*make_SP(q2,q2)^q2_q2_exp*make_SP(q1,q2)^q1_q2_exp
      end # for q2_q2_exp
    end # for q1_q1_exp
  end # for q1_q2_exp

  SP_list = get_add_vector_expand(total_term)
  SP_list = sort( map(drop_coeff,SP_list), by=gen_sorted_str )

  return SP_list

end # function generate_SPcombo





#########################################
function gen_SPcombo_v2(
    rank_list::Vector{Int64},
    ind_mom_list::Array{Basic}
)::Array{Basic}
#########################################

  # e.g. gen_SPcombo_v2( [1,2,0],[k1,k2,K3])
  # only for three-loop level

  @vars q1,q2,q3
  @funs SP

  @assert length(rank_list) == 3
  ori_n_q1 = rank_list[1] # count( "q1", rank_str )
  ori_n_q2 = rank_list[2] # count( "q2", rank_str )
  ori_n_q3 = rank_list[3] # count( "q3", rank_str )

  # term_q1ki = SP(q1,k1)+SP(q1,k2)+...
  term_q1ki = sum( mom -> make_SP(q1,mom), ind_mom_list )
  # term_q2ki = SP(q2,k1)+SP(q2,k2)+...
  term_q2ki = sum( mom -> make_SP(q2,mom), ind_mom_list )
  # term_q3ki = SP(q3,k1)+SP(q3,k2)+...
  term_q3ki = sum( mom -> make_SP(q3,mom), ind_mom_list )

  #q1q2_xpt is the exponent of SP(q1,q2)
  #q2q3_xpt is the exponent of SP(q2,q3)
  #q1q3_xpt is the exponent of SP(q1,q3)
  #q1q1_xpt is the exponent of SP(q1,q1)
  #q2q2_xpt is the exponent of SP(q2,q2)
  #q3q3_xpt is the exponent of SP(q3,q3)
  #q1ki_xpt is the exponent of SP(q1,k1)+SP(q1,k2)+...
  #q2ki_xpt is the exponent of SP(q2,k1)+SP(q2,k2)+...
  #q3ki_xpt is the exponent of SP(q3,k1)+SP(q3,k2)+...
  #total_term is the summation of all possible terms
  #expand the total_term and drop the coefficients

  total_term = zero(Basic)
  for q1q2_xpt in 0:min(ori_n_q1,ori_n_q2)
    n_q1 = ori_n_q1
    n_q2 = ori_n_q2
    n_q3 = ori_n_q3

    n_q1 -= q1q2_xpt
    n_q2 -= q1q2_xpt
    for q1q3_xpt in 0:min(n_q1,n_q3)
      n_q1 -= q1q3_xpt
      n_q3 -= q1q3_xpt
      for q2q3_xpt in 0:min(n_q2,n_q3)
        n_q2 -= q2q3_xpt
        n_q3 -= q2q3_xpt

        max_q1q1_xpt = floor(Int,n_q1/2)
        max_q2q2_xpt = floor(Int,n_q2/2)
        max_q3q3_xpt = floor(Int,n_q3/2)

        for q1q1_xpt in 0:max_q1q1_xpt, q2q2_xpt in 0:max_q2q2_xpt, q3q3_xpt in 0:max_q3q3_xpt
          q1ki_xpt = n_q1-2*q1q1_xpt
          q2ki_xpt = n_q2-2*q2q2_xpt
          q3ki_xpt = n_q3-2*q3q3_xpt

          total_term += term_q1ki^q1ki_xpt * term_q2ki^q2ki_xpt * term_q3ki^q3ki_xpt *
                  make_SP(q1,q1)^q1q1_xpt * make_SP(q2,q2)^q2q2_xpt * make_SP(q3,q3)^q3q3_xpt *
                  make_SP(q1,q2)^q1q2_xpt * make_SP(q1,q3)^q1q3_xpt * make_SP(q2,q3)^q2q3_xpt 
            
        end # for q1q1_xpt, q2q2_xpt, q3q3_xpt
        
      end # for q2q3_xpt
    end # for q1q3_xpt
  end # for q1q2_xpt

  SP_list = get_add_vector_expand(total_term)
  SP_list = sort( map(drop_coeff,SP_list), by=gen_sorted_str )

  return SP_list

end # function gen_SPcombo_v2





#########################################
# Created by Quan-feng WU
# Mar. 28, 2023
function gen_SPcombo_v3(
  rank_list::Vector{Int64},
  ind_mom_list::Array{Basic}
)::Array{Basic}
#########################################

  q_list = reduce(vcat,
    [ Basic("q$ii") for _ ∈ 1:rank_list[ii] ]
      for ii ∈ eachindex(rank_list)
  ) # end reduce

  SP_list = Basic[]

  for num_single_q ∈ 0:length(q_list)
    if isodd(length(q_list) - num_single_q)
      continue
    end # if

    for single_q_list ∈ multiset_combinations( q_list, num_single_q )

      pair_q_list = copy(q_list)
      for single_q ∈ single_q_list
        q_index = findfirst( isequal(single_q), pair_q_list )
        @assert !isnothing(q_index)
        deleteat!( pair_q_list, q_index )
      end # for single_q

      single_q_SP_list = Basic[]
      for ind_mom_comb ∈ with_replacement_combinations(ind_mom_list, num_single_q)
        for ind_mom_ord ∈ multiset_permutations(ind_mom_comb, num_single_q)
          push!( single_q_SP_list, prod( make_SP.( ind_mom_ord, single_q_list ) ) )
        end # for ind_mom_ord
      end # for ind_mom_comb

      pair_q_SP_list = Basic[]
      if isempty(pair_q_list)
        push!( pair_q_SP_list, one(Basic) )
      else
        pair_q_SP_list = Basic[]
        for pair_q_partition ∈ partitions( pair_q_list, 2 )
          if (length∘first)(pair_q_partition) * 2 != length(pair_q_list)
            continue
          end # if
          left_pair_q, right_pair_q = pair_q_partition

          for right_pair_q_ord ∈ multiset_permutations( right_pair_q, length(right_pair_q) )
            push!( pair_q_SP_list, prod( make_SP.( left_pair_q, right_pair_q_ord ) ) )
          end # for right_pair_q_ord
        end # for pair_q_partition
      end # if

      unique!(single_q_SP_list)
      unique!(pair_q_SP_list)
      for single_q_SP ∈ single_q_SP_list, pair_q_SP ∈ pair_q_SP_list
        push!( SP_list, single_q_SP * pair_q_SP )
      end # for single_q_SP
    end # for single_q_list
  end # for num_single_q

  unique!( SP_list )
  sort!( SP_list; by=gen_sorted_str )

  return SP_list
end # function gen_SPcombo_v3

