


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
# Created by Quan-feng WU
# Mar. 28, 2023
function gen_SPcombo(
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

end # function gen_SPcombo

