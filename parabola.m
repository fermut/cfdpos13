function f = parabola(N)
  x = linspace(0,1,N);
  f = zeros(N,1);
  for i=1:N
    f(i) = 4*x(i)*(1-x(i));
  endfor
endfunction
