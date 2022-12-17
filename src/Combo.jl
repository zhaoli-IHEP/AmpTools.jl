
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
    term_q1_k += SP(q1,independent_mom_list[i])
  end # for

  # term_q2_k = SP(q2,k1)+SP(q2,k2)+...
  term_q2_k = 0
  for i in 1:num_mom
    term_q2_k += SP(q2,independent_mom_list[i])
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
  SP_list = sort( map(drop_coeff,SP_list), by=gen_sorted_str )

  return SP_list

end # function generate_SPcombo





#########################################
function gen_SPcombo_v2(
    rank_str::String,
    ind_mom_list::Array{Basic}
)::Array{Basic}
#########################################

  # e.g. generate_SP("q1q1q2q2",[k1,k2,K3])
  # only for three-loop level

  @vars q1,q2,q3
  @funs SP

  n_q1 = count( "q1", rank_str )
  n_q2 = count( "q2", rank_str )
  n_q3 = count( "q3", rank_str )

  # term_q1ki = SP(q1,k1)+SP(q1,k2)+...
  term_q1ki = (sum∘map)( mom -> SP(q1,mom), ind_mom_list )
  # term_q2ki = SP(q2,k1)+SP(q2,k2)+...
  term_q2ki = (sum∘map)( mom -> SP(q2,mom), ind_mom_list )
  # term_q3ki = SP(q3,k1)+SP(q3,k2)+...
  term_q3ki = (sum∘map)( mom -> SP(q3,mom), ind_mom_list )

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

  total_term = 0
  for q1q2_xpt in 0:min(n_q1,n_q2)
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
                  SP(q1,q1)^q1q1_xpt * SP(q2,q2)^q2q2_xpt * SP(q3,q3)^q3q3_xpt *
                  SP(q1,q2)^q1q2_xpt * SP(q1,q3)^q1q3_xpt * SP(q2,q3)^q2q3_xpt 
            
        end # for q1q1_xpt, q2q2_xpt, q3q3_xpt
        
      end # for q2q3_xpt
    end # for q1q3_xpt
  end # for q1q2_xpt

  SP_list = (get_add_vector∘expand)(total_term)
  SP_list = sort( map(drop_coeff,SP_list), by=gen_sorted_str )

  return SP_list

end # function gen_SPcombo_v2

