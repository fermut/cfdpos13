function ng = ij2n (i,j)
  global NI NJ
  if(NI < NJ)
    ng = i + (j-1)*NI;
  else
    ng = j + (i-1)*NJ;
  endif
endfunction

