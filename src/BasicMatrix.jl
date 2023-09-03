


############################################
# Created by Quan-feng Wu, Feb 21, 2023
function get_det_naive( 
    MM::Matrix{Basic}, 
    shape_mat::Matrix{Int64} = Matrix{Int64}(undef,0,0) 
)::Basic
############################################

  nr, nc = size(MM)
  @assert nr == nc

  if nr == 1
    return MM[1,1]
  end # if

  if nr == 2 
    return MM[1,1]*MM[2,2] - MM[1,2]*MM[2,1]
  end # if

  if isempty(shape_mat)
    shape_mat = zeros( Int64, nr, nc )
    for rr in 1:nr, cc in 1:nc
      shape_mat[rr,cc] = iszero(MM[rr,cc]) ? zero(Int64) : one(Int64)
    end # for rr, cc
  end # if
  
  # larger than 2x2 uses the Laplace expansion
  detMM = zero(Basic)

  row_count_list = map( rr -> sum(shape_mat[rr,1:nc]), 1:nr )
  min_sum, chosen_rr = findmin(row_count_list)
  if iszero(min_sum)
    return zero(Basic)
  end # if

  for cc in 1:nc
    element = MM[chosen_rr,cc]
    if iszero(element)
      continue
    end # if
    cofactor = (-1)^(chosen_rr+cc)
    sub_mat = MM[ setdiff(1:nr,chosen_rr), setdiff(1:nc,cc) ]
    sub_shape_mat = shape_mat[ setdiff(1:nr,chosen_rr), setdiff(1:nc,cc) ]
    sub_det = get_det_naive(sub_mat,sub_shape_mat)
    detMM += element * cofactor * sub_det
  end # for cc

  return expand(detMM)

end # function get_det_naive


############################################
# Created by Quan-feng Wu, Feb 21, 2023
function get_det( 
    MM::Matrix{Basic} 
)::Basic
############################################

  nr, nc = size(MM)
  @assert nr == nc

  if nr == 1
    return MM[1, 1]
  end

  if nr == 2
    return MM[1, 1] * MM[2, 2] - MM[1, 2] * MM[2, 1]
  end

  det = 0
  sub_dim = (Int∘floor)(nr / 2)
  for selected_cols in combinations(1:nc, sub_dim)
    cofactor = (-1)^(sum(1:sub_dim) + sum(selected_cols))
    fac1 = get_det( MM[1:sub_dim, selected_cols] )
    if iszero(fac1)
      continue
    end # if
    fac2 = get_det( MM[(sub_dim+1):end, setdiff(1:nc, selected_cols)] ) 
    det += fac1 * cofactor * fac2
  end # for selected_cols

  return  det

end # function get_det


#############################
function lu_decomposition(
    A::Matrix{Basic} 
)::Tuple{Matrix{Basic},Matrix{Basic}} 
#############################

  n = size(A, 1)
  L = one( Matrix{Basic}(undef, n, n) )
  U = deepcopy(A)
    
  for i in 1:n
    for j in 1:i-1
      sum = 0.0
      for k in 1:j-1
        sum += L[i, k] * U[k, j]
      end # for k
      U[i, j] -= sum
    end # for j
        
    for j in i+1:n
        sum = 0.0
        for k in 1:i-1
          sum += L[i, k] * U[k, j]
        end # for k
        L[i, j] = (A[i, j] - sum) / U[i, i]
    end # for j
  end # for i
    
  return L, U

end # function lu_decomposition



#################################
function get_adj_naive( 
    MM::Matrix{Basic} 
)::Matrix{Basic}
#################################

  nr, nc = size( MM )
  @assert nr == nc

  adjMM = Matrix{Basic}(undef,nr,nc)
  if nr == 1

    adjMM[1,1] = Basic(1)

  elseif nr == 2

    adjMM[1,1] = MM[2,2]
    adjMM[1,2] = -MM[1,2]
    adjMM[2,1] = -MM[2,1]
    adjMM[2,2] = MM[1,1]

  elseif nr == 3

    adjMM[1,1] = -MM[2,3]*MM[3,2] + MM[2,2]*MM[3,3] 
    adjMM[1,2] = MM[1,3]*MM[3,2] - MM[1,2]*MM[3,3]
    adjMM[1,3] = -MM[1,3]*MM[2,2] + MM[1,2]*MM[2,3] 
    adjMM[2,1] = MM[2,3]*MM[3,1] - MM[2,1]*MM[3,3]
    adjMM[2,2] = -MM[1,3]*MM[3,1] + MM[1,1]*MM[3,3]
    adjMM[2,3] = MM[1,3]*MM[2,1] - MM[1,1]*MM[2,3]
    adjMM[3,1] = -MM[2,2]*MM[3,1] + MM[2,1]*MM[3,2]
    adjMM[3,2] = MM[1,2]*MM[3,1] - MM[1,1]*MM[3,2]
    adjMM[3,3] = -MM[1,2]*MM[2,1] + MM[1,1]*MM[2,2]

  else
    error( "Larger nr is not expected!" )
  end # if

  return adjMM

end # function get_adj_naive







####################################################
function get_adj( MM::Matrix{Basic} )::Matrix{Basic}
####################################################

  nr, nc = size( MM )
  @assert nr == nc

  if nr == 1
    return ones(Basic, 1, 1)
  end # if

  if nr == 2
    adjMM = Matrix{Basic}(undef,nr,nc)
    adjMM[1,1] = MM[2,2]
    adjMM[1,2] = -MM[1,2]
    adjMM[2,1] = -MM[2,1]
    adjMM[2,2] = MM[1,1]
    return adjMM
  end # if

  # for larger matrix, notice there should a transpose.
  adjMM = Matrix{Basic}(undef,nr,nc)
  for rr in 1:nr, cc in 1:nc
    cofactor = (-1)^(rr+cc)
    sub_mat = MM[ setdiff(1:nr,cc), setdiff(1:nc,rr) ]
    adjMM[rr,cc] = cofactor * get_det(sub_mat)
  end # for rr, cc
  
  return adjMM

end # function get_adj





################################
function get_dot( 
    m1::Matrix{Basic}, 
    m2::Matrix{Basic} 
)::Basic
################################

  nr1, nc1 = size(m1)
  nr2, nc2 = size(m2)

  @assert nr1 == nr2 && nc1 == nc2
  nr = nr1
  nc = nc1

  if nr == 1
    vec1 = collect( view( m1, 1, : ) )
    vec2 = collect( view( m2, 1, : ) )
  elseif nc == 1
    vec1 = collect( view( m1, :, 1 ) )
    vec2 = collect( view( m2, :, 1 ) )
  else
    error( "nr, nc is not expected: $(nr), $(nc)" )
  end # if
 
  return sum( vec1 .* vec2 )

end # function get_dot



########################
function rref(
    A::Matrix{Basic}
)::Matrix{Basic}
########################

  m, n = size(A)
  lead = 1

  for r in 1:m
    if n < lead
      return A
    end

    i = r
    while (iszero∘expand)( A[i, lead] )
      i += 1
      if m < i
        i = r
        lead += 1
        if n < lead
          return A
        end
      end
    end

    if i != r
      A[[i, r], :] = A[[r, i], :]
    end

    pivot = A[r, lead]
    A[r, :] = A[r, :] / pivot


    for i in 1:m
      if i != r
        factor = A[i, lead]
        A[i, :] = A[i, :] - factor * A[r, :] 
      ##A[i, :] = map( x->begin
      ##  numer, denom = get_numer_denom_fermat(x)
      ##  numer/denom
      ##end, A[i, :] - factor * A[r, :] )
      end
    end

    A = thread_rationalize_fermat(A)

    lead += 1
  end

  return A

end # function rref





############################
function calc_null_space(
    mat::Matrix{Basic}
)::Matrix{Basic}
############################

  mat = Matrix(mat)  # Convert input matrix to a numeric matrix
  nr, nc = size(mat)
  
  if iszero(mat)
    return one( Matrix{eltype(mat)}(undef, nc, nc) )
  end
  
  rref_mat = rref(mat)

  pivot_cols = Vector{Int64}()
  for row in 1:nr
    if all(iszero∘expand, rref_mat[row, :])
      continue
    end # if
    for col in 1:nc
      if !(iszero∘expand)(rref_mat[row, col])
        push!(pivot_cols, col)
        break
      end # if
    end # for col
  end # for row
    
  null_cols = setdiff(collect(1:nc), pivot_cols)
  
  if length(null_cols) == 0
    return Matrix{eltype(mat)}(undef, nc, 0)
  end # if

  ns_mat = zeros(eltype(mat), nc, length(null_cols))
  for (i, null_col) in enumerate(null_cols)
    ns_mat[null_col,i] = one(Basic)
    for rref_rr in 1:length(pivot_cols)
      pivot_col = pivot_cols[rref_rr]
      ns_mat[pivot_col,i] = -rref_mat[rref_rr,null_col]
    end # for pivot_col
  end # for

  @assert all( iszero∘expand, mat*ns_mat )

  return ns_mat

end # function calc_null_space



#################################
function get_matrix_shape_str(
    mat::Matrix{Basic}
)::String
#################################

  nr, nc = size(mat)
  result_str = string()
  for rr in 1:nr
    for cc in 1:nc
    ele = mat[rr,cc]
      if iszero(ele)
        result_str *= " "
      else
        result_str *= "⋅"
      end # if
    end # for cc
    result_str *= "\n"
  end # for rr

  return result_str

end # function get_matrix_shape_str


