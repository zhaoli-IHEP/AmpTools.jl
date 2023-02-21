


############################################
# Created by Quan-feng Wu, Feb 21, 2023
function get_det_old( MM::Matrix{Basic} )::Basic
############################################

  nr, nc = size(MM)
  @assert nr == nc

  if nr == 1
    return MM[1,1]
  end # if

  if nr == 2 
    return MM[1,1]*MM[2,2] - MM[1,2]*MM[2,1]
  end # if
  
  # larger than 2x2 uses the Laplace expansion
  detMM = zero(Basic)

  for cc in 1:nc
    element = MM[1,cc]
    cofactor = (-1)^(1+cc)
    sub_mat = MM[ 2:nr, setdiff(1:nc,cc) ]
    detMM += element * cofactor * get_det_old(sub_mat)
  end # for cc

  return detMM

end # function get_det_old


############################################
# Created by Quan-feng Wu, Feb 21, 2023
function get_det( MM::Matrix{Basic} )::Basic
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
    det += get_det( MM[1:sub_dim, selected_cols] ) * cofactor *
      get_det( MM[sub_dim+1:end, setdiff(1:nc, selected_cols)] )
  end # for selected_cols

  return  det

end # function get_det


####################################################
function get_adj_old( MM::Matrix{Basic} )::Matrix{Basic}
####################################################

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

end # function get_adj_old







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





###############################################################
function get_dot( m1::Matrix{Basic}, m2::Matrix{Basic} )::Basic
###############################################################

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
    error( "nr, nc is not expected: "*string(nr)*", "*string(nc) )
  end # if
 
  return sum( vec1.*vec2 )

end # function get_dot












