



#######################################
function get_rref_fermat( 
    mat::Matrix{Basic}, 
    name_str::String = "rref"
)::Matrix{Basic}
#######################################

  if all(iszero,mat)
    return mat
  end # if

  symbol_str_list = map( string, free_symbols(mat) )
  symbol_str_list = filter( x->x!="im", symbol_str_list )
  capital_symbol_list = (to_Basic∘filter)( x->x!=lowercase(x), symbol_str_list )
  lowered_symbol_list = (to_Basic∘map)( lowercase∘string, capital_symbol_list )

  symbol_decl = join( [ "&(J=$(x));" for x in map(lowercase,symbol_str_list) ], "\n" )

  lower_dict = Dict{Basic,Basic}( capital_symbol_list .=> lowered_symbol_list )
  recover_dict = Dict{Basic,Basic}( lowered_symbol_list .=> capital_symbol_list )

  nr, nc = size(mat)

  init_str = string()
  for rr in 1:nr, cc in 1:nc
    ele = mat[rr,cc]
    if !iszero(ele) && !isempty(lower_dict)
      ele = subs( ele, lower_dict )
    end # if
    init_str *= "mat[$rr,$cc]:=$(ele);\n"
  end # for rr, cc

  # The reshape in julia is col-wise
  out_str = string()
  for cc in 1:nc, rr in 1:nr 
    out_str *= "!!(&o, mat[$rr,$cc]);\n"
    if rr != nr || cc != nc 
      out_str *= "!!(&o, ',');\n"
    end # if
  end # for rr, cc

  fer_script_str = """
  &(N=0);
  &(t=0);
  &(_t=0);
  &(_d=100000);
  &(M='');
  &(_s=0);
  
  &(S='$(name_str).out');
  
  &(J=im);
  $(symbol_decl)
  &(P=im^2+1,1);

  Array mat[$nr,$(nc+1)] Sparse;
  
  $(init_str) 

  Redrowech([mat]);
  
  &(U=1);
  
  $(out_str)
  
  &(U=0);
  
  &q;
  """

  file = open( "$(name_str).fer", "w" )
  write( file, fer_script_str )
  close(file)

  run( pipeline( `fermat`, stdin="$(name_str).fer", stdout="$(name_str).log" ) )

  out_file = open( "$(name_str).out", "r" )
  result_str = read( out_file, String ) 
  close( out_file )

  str_list = split( replace( result_str, "\n"=>"" ), "," )
  @assert length(str_list) == nr*nc
  ele_list = map( to_Basic, str_list )
  if !isempty(recover_dict)
    ele_list = map( x -> subs( x, recover_dict ), ele_list )
  end # if
  rref_mat = reshape( ele_list, nr, nc )

  rm( "$(name_str).fer" )
  rm( "$(name_str).out" )
  rm( "$(name_str).log" )

  return rref_mat

end # function get_rref_fermat



#############################################
function rational_function_simplify( 
    rat::Basic, 
    surfix::String = "" 
)::Basic
#############################################

  symbol_str_list = map( string, free_symbols(rat) )
  symbol_str_list = filter( x->x!="im", symbol_str_list )
  capital_symbol_list = (to_Basic∘filter)( x->x!=lowercase(x), symbol_str_list )
  lowered_symbol_list = (to_Basic∘map)( lowercase∘string, capital_symbol_list )

  symbol_decl = join( [ "&(J=$(x));" for x in map(lowercase,symbol_str_list) ], "\n" )

  lower_dict = Dict{Basic,Basic}( capital_symbol_list .=> lowered_symbol_list )
  recover_dict = Dict{Basic,Basic}( lowered_symbol_list .=> capital_symbol_list )
  if !isempty(lower_dict)
    rat = subs( rat, lower_dict )
  end # if

  fer_script_str = """
  &(N=0);
  &(t=0);
  &(_t=0);
  &(_d=100000);
  &(M='');
  &(_s=0);
  
  &(S='rational$(surfix).out');
  
  &(J=im);
  $(symbol_decl)
  &(P=im^2+1,1);
  
  expr:=$(rat);
  
  
  &(U=1);
  
  !!(&o, expr);
  
  &(U=0);
  
  &q;
  """

  file = open( "rational$(surfix).fer", "w" )
  write( file, fer_script_str )
  close(file)

  run( pipeline( `fermat`, stdin="rational$(surfix).fer", stdout="rational$(surfix).log" ) )

  out_file = open( "rational$(surfix).out", "r" )
  result_expr = (to_Basic∘read)( out_file, String ) 
  close( out_file )

  if !isempty(recover_dict)
    result_expr = subs( result_expr, recover_dict )
  end # if

  rm( "rational$(surfix).fer" )
  rm( "rational$(surfix).out" )
  rm( "rational$(surfix).log" )

  return result_expr

end # function rational_function_simplify



#############################################
function rational_function_simplify( 
    mat::Matrix{Basic}, 
    surfix::String = "";
    verbose = true
)::Matrix{Basic}
#############################################

  cart = CartesianIndices(mat)
  result_mat = zero(mat)
  Threads.@threads for index in 1:length(cart)
    if verbose == true
      printstyled( ".", color=:light_yellow )
    end # if
    rr = cart[index][1]
    cc = cart[index][2]
    if !iszero(mat[rr,cc])
      result_mat[rr,cc] = rational_function_simplify( mat[rr,cc], "r$(rr)c$(cc)" )
    end # if
    if verbose == true
      printstyled( ".", color=:light_green )
    end # if
  end # for index

  if verbose == true
    println()
  end # if

  return result_mat

end # function rational_function_simplify









#############################################
function numer_denom_fermat( 
    rat::Basic, 
    surfix::String = "" 
)::Tuple{Basic,Basic}
#############################################

  symbol_str_list = map( string, free_symbols(rat) )
  symbol_str_list = filter( x->x!="im", symbol_str_list )
  capital_symbol_list = (to_Basic∘filter)( x->x!=lowercase(x), symbol_str_list )
  lowered_symbol_list = (to_Basic∘map)( lowercase∘string, capital_symbol_list )

  symbol_decl = join( [ "&(J=$(x));" for x in map(lowercase,symbol_str_list) ], "\n" )

  lower_dict = Dict{Basic,Basic}( capital_symbol_list .=> lowered_symbol_list )
  recover_dict = Dict{Basic,Basic}( lowered_symbol_list .=> capital_symbol_list )
  if !isempty(lower_dict)
    rat = subs( rat, lower_dict )
  end # if

  fer_script_str = """
  &(N=0);
  &(t=0);
  &(_t=0);
  &(_d=100000);
  &(M='');
  &(_s=0);
  
  &(S='numerdenom$(surfix).out');
  
  &(J=im);
  $(symbol_decl)
  &(P=im^2+1,1);
  
  expr:=$(rat);
  
  &(U=1);
  
  !!(&o, Numer(expr));
  !!(&o, ',');
  !!(&o, Denom(expr));
  
  &(U=0);
  
  &q;
  """

  file = open( "numerdenom$(surfix).fer", "w" )
  write( file, fer_script_str )
  close(file)

  run( pipeline( `fermat`, stdin="numerdenom$(surfix).fer", stdout="numerdenom$(surfix).log" ) )

  out_file = open( "numerdenom$(surfix).out", "r" )
  result_str = read( out_file, String ) 
  close( out_file )

  str_list = split( result_str, "," )
  @assert length(str_list) == 2
  numer = (to_Basic∘string)(str_list[1]) 
  denom = (to_Basic∘string)(str_list[2]) 

  if !isempty(recover_dict)
    numer = subs( numer, recover_dict )
    denom = subs( denom, recover_dict )
  end # if

  rm( "numerdenom$(surfix).fer" )
  rm( "numerdenom$(surfix).out" )
  rm( "numerdenom$(surfix).log" )

  return numer, denom

end # function numer_denom_fermat 










