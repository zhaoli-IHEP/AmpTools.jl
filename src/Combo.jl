
#########################################
function drop_coeff( expr::Basic )::Basic
#########################################

  return SymEngine.get_symengine_class( expr ) == :Mul ?
         (prod∘filter)( x_-> SymEngine.get_symengine_class(x_) != :Integer, get_args(expr) ) : expr

end # function drop_coeff

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
    term_q1_k = SP(q1,independent_mom_list[i])+term_q1_k
  end # for

  # term_q2_k = SP(q2,k1)+SP(q2,k2)+...
  term_q2_k = 0
  for i in 1:num_mom
    term_q2_k = SP(q2,independent_mom_list[i])+term_q2_k
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
                      SP(q1,q1)^q1_q1_exp*SP(q2,q2)^q2_q2_exp*SP(q1,q2)^q1_q2_exp
      end # for q2_q2_exp
    end # for q1_q1_exp
  end # for q1_q2_exp

  SP_list = (get_add_vector∘expand)(total_term)
  SP_list = map(drop_coeff,SP_list)

  return SP_list

end # function generate_SPcombo

