function ng = ij2nu (i,j)
  global NI NJ
  if(NI+1 < NJ)
    ng = i + (j-1)*(NI+1);
  else
    ng = j + (i-1)*NJ;
  endif
#if (ng > (NI+1)*NJ)
#  disp("OUT OF RAGE i,j=",i,j); fflush(stdout);
#endif
endfunction

