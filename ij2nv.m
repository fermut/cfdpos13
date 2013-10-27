function ng = ij2nv (i,j)
  global NI NJ
  if(NI < NJ+1)
    ng = i + (j-1)*NI;
  else
    ng = j + (i-1)*(NJ+1);
  endif
#if (ng > NI*(NJ+1))
#  disp("OUT OF RAGE i,j=",i,j); fflush(stdout);
#endif
endfunction

