



#######################################
function get_rref_fermat( 
    mat::Matrix{Basic}, 
    name_str::String = "rref"
)::Matrix{Basic}
#######################################

  if all(iszero,mat)
    return mat
  end # if

  symbol_list = filter( s_->s_!="im", string.(free_symbols(mat)) )
  symbol_decl = join( [ "&(J=$(str));" for str in symbol_list ], "\n" )

  nr, nc = size(mat)

  init_str = string()
  for rr in 1:nr, cc in 1:nc
    init_str *= "mat[$rr,$cc]:=$(mat[rr,cc]);\n"
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

  symbol_list = filter( x->x!="im", string.( free_symbols(rat) ) )
  symbol_decl = join( [ "&(J=$(str));" for str in symbol_list ], "\n" )

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
  result_expr = Basic( read( out_file, String ) )
  close( out_file )

  result_expr = subs( result_expr, Basic("im")=>im )

  rm( "rational$(surfix).fer" )
  rm( "rational$(surfix).out" )
  rm( "rational$(surfix).log" )

  return result_expr

end # function rational_function_simplify
