#Problem definition
global NI=3*12*2+1
global NJ=3* 2*2+2
Lx = 12.0
Ly = 2.0
dx = Lx/(NI-1)
dy = Ly/(NJ-2)
h = (dx+dy)/2

#Material properties
nu = 1/100;
rho= 1;

#Temporal integration
dt = 0.001;
T0 = 0;
T1 = 100;
NT = (T1-T0)/dt;

#Arrays
dimP = NI*NJ;
dimU = (NI+1)*NJ;
dimV = NI*(NJ+1);

Ap = spalloc(dimP,dimP,5*dimP);
Dx = spalloc(dimU,dimU,5*dimU);
Dy = spalloc(dimV,dimV,5*dimV);

Apu = spalloc(dimP,dimU,4*dimP);
Apv = spalloc(dimP,dimV,4*dimP);

AxuE = spalloc(dimU,dimU,2*dimU);
AxuW = spalloc(dimU,dimU,2*dimU);
AxuN = spalloc(dimU,dimU,2*dimU);
AxuS = spalloc(dimU,dimU,2*dimU);
AxvN = spalloc(dimU,dimV,2*dimU);
AxvS = spalloc(dimU,dimV,2*dimU);

AyvE = spalloc(dimV,dimV,2*dimV);
AyvW = spalloc(dimV,dimV,2*dimV);
AyvN = spalloc(dimV,dimV,2*dimV);
AyvS = spalloc(dimV,dimV,2*dimV);
AyuE = spalloc(dimV,dimU,2*dimV);
AyuW = spalloc(dimV,dimU,2*dimV);

gradU = spalloc(dimP,dimU,2*dimP);
gradV = spalloc(dimP,dimV,2*dimP);
gradPU = spalloc(dimU,dimP,2*dimU);
gradPV = spalloc(dimV,dimP,2*dimV);

#Mask definition
INNER=0;
WALL=1;
INFLOW=2;
OUTFLOW=3;
mask=INNER*ones(dimP,1);
disp("creating mask ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     ###DOMAIN BOUNDARY
     if (i == 1)
       mask(ij2n(i,j))=INFLOW;
     endif
     if (i == NI)
       mask(ij2n(i,j))=OUTFLOW;
     endif
     if ((j == 1)||(j == NJ))
       mask(ij2n(i,j))=WALL;
     endif
     ###INNER BOUNDARY
     if ((i <= (NI-1)/3+1)&&(j <= (NJ-2)/2+1))
       mask(ij2n(i,j))=WALL;
     endif
   endfor
endfor

#Mask transcription to U-,V-cells
masku=INNER*ones(dimU);
maskv=INNER*ones(dimV);
disp("transcriving mask to u,v cells ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     if (mask(ij2n(i,j)) != INNER)
       masku(ij2nu(i,j))  =mask(ij2n(i,j)); 
       masku(ij2nu(i+1,j))=mask(ij2n(i,j)); 
       maskv(ij2nv(i,j))  =mask(ij2n(i,j)); 
       maskv(ij2nv(i,j+1))=mask(ij2n(i,j)); 
     endif
     if (mask(ij2n(i,j)) == OUTFLOW)
       masku(ij2nu(i,j)) = INNER;
     endif
   endfor
endfor

#Boundary Conditions
pbc = zeros(dimP,1);
ubc = zeros(dimU,1);
vbc = zeros(dimV,1);
N=0;
for j=1:NJ
  if (masku(ij2nu(2,j)) == INFLOW)
    N = N + 1;
  endif
endfor
UI=parabola(N);
N=1;
for j=1:NJ
  if (masku(ij2nu(2,j)) == INFLOW)
    ubc(ij2nu(2,j)) = UI(N);
    N = N + 1;
  endif
endfor

#P-cells
disp("P-cells: pressure laplacian and velocity gradients ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     ######################################
     #### PRESSURE LAPLACIAN
     if (mask(ij2n(i,j)) == INNER)
       #-- east face
       if (masku(ij2nu(i+1,j)) == INNER)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i+1,j)) =  1/h/h;
       else
         Apu(ij2n(i,j),ij2nu(i+1,j)) = 1;
       endif
       #-- west face
       if (masku(ij2nu(i,j)) == INNER)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i-1,j)) =  1/h/h;
       else
         Apu(ij2n(i,j),ij2nu(i,j)) = 1;
       endif
       #-- north face
       if (maskv(ij2nv(i,j+1)) == INNER)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j+1)) =  1/h/h;
       else
         Apv(ij2n(i,j),ij2nv(i,j+1)) = 1;
       endif
       #-- south face
       if (maskv(ij2nv(i,j)) == INNER)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j-1)) =  1/h/h;
       else
         Apv(ij2n(i,j),ij2nv(i,j)) = 1;
       endif
     else
       #not inner: left unchanged
       Ap(ij2n(i,j),ij2n(i,j)) = 1;
     endif
     ######################################
     #### VELOCITY GRADIENTS
     if (mask(ij2n(i,j)) == INNER)
       gradU(ij2n(i,j),ij2nu(i,j)) = -1/h;
       gradU(ij2n(i,j),ij2nu(i+1,j)) = 1/h;
       gradV(ij2n(i,j),ij2nv(i,j)) = -1/h;
       gradV(ij2n(i,j),ij2nv(i,j+1)) = 1/h;
     endif
   endfor
endfor

#U-cells
disp("U-cells: face velocities and viscous laplacian ..."); fflush(stdout);
for i=1:NI+1
   for j=1:NJ
     ######################################
     #### FACE VELOCITIES (NORMAL)
     ####   ux{E,W,N,S} = Axu{E,W,N,S}*u
     if (masku(ij2nu(i,j)) == INNER)
       #-- east face
       if (masku(ij2nu(i+1,j)) == OUTFLOW)
         # zero gradient
         AxuE(ij2nu(i,j),ij2nu(i,j)) = 1;  
       else
         AxuE(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuE(ij2nu(i,j),ij2nu(i+1,j)) = 1/2;
       endif
       #-- west face
       AxuW(ij2nu(i,j),ij2nu(i,j)) = 1/2;
       AxuW(ij2nu(i,j),ij2nu(i-1,j)) = 1/2;
       #--north face
       if (masku(ij2nu(i,j+1)) != WALL)
         AxuN(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuN(ij2nu(i,j),ij2nu(i,j+1)) = 1/2;
       endif
       #--south face
       if (masku(ij2nu(i,j-1)) != WALL)
         AxuS(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuS(ij2nu(i,j),ij2nu(i,j-1)) = 1/2;
       endif
     ######################################
     #### FACE VELOCITIES (TANGENCIAL)
     ####   vx{N,S}     = Axv{N,S}*v
       #--north face
       if (masku(ij2nu(i,j+1)) == INNER)
         AxvN(ij2nu(i,j),ij2nv(i,j+1)) = 1/2;
         AxvN(ij2nu(i,j),ij2nv(i-1,j+1)) = 1/2;
       endif
       #--south face
       if (masku(ij2nu(i,j-1)) == INNER)
         AxvS(ij2nu(i,j),ij2nv(i-1,j)) = 1/2;
         AxvS(ij2nu(i,j),ij2nv(i,j)) = 1/2;
       endif
     endif
     ######################################
     #### VISCOUS LAPLACIAN
     if (masku(ij2nu(i,j)) == INNER)
       #-- east face
       if (masku(ij2nu(i+1,j)) == INNER)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i+1,j)) =  1/h/h;
       elseif (masku(ij2nu(i+1,j)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
       elseif (masku(ij2nu(i+1,j)) == OUTFLOW)
         #zero gradient
         #do nothing
       endif
       #-- west face
       if (masku(ij2nu(i-1,j)) == INNER)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i-1,j)) =  1/h/h;
       elseif (masku(ij2nu(i-1,j)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
       elseif (masku(ij2nu(i-1,j)) == INFLOW)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i-1,j)) =  1/h/h;
       endif
       #-- north face
       if (masku(ij2nu(i,j+1)) == INNER)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i,j+1)) =  1/h/h;
       elseif (masku(ij2nu(i,j+1)) == WALL)
         #ghost point
         Dx(ij2nu(i,j),ij2nu(i,j))  += -2/h/h;
       endif
       #-- south face
       if (masku(ij2nu(i,j-1)) == INNER)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i,j-1)) =  1/h/h;
       elseif (masku(ij2nu(i,j-1)) == WALL)
         #ghost point
         Dx(ij2nu(i,j),ij2nu(i,j))  += -2/h/h;
       endif
     else
       #not inner: left unchanged
       Dx(ij2nu(i,j),ij2nu(i,j)) = 1;
     endif
     ######################################
     #### PRESSURE GRADIENT
     if (masku(ij2nu(i,j)) == INNER)
       gradPU(ij2nu(i,j),ij2n(i,j)) = 1/h;
       gradPU(ij2nu(i,j),ij2n(i-1,j)) = -1/h;
     endif
   endfor
endfor

#V-cells
disp("V-cells: face velocities and viscous laplacian ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ+1
     ######################################
     #### FACE VELOCITIES (NORMAL)
     ####   vy{N,S,E,W} = Ayv{N,S,E,W}*v
     if (maskv(ij2nv(i,j)) == INNER)
       #--north face
       AyvN(ij2nv(i,j),ij2nv(i,j)) = 1/2;
       AyvN(ij2nv(i,j),ij2nv(i,j+1)) = 1/2;
       #--south face
       AyvS(ij2nv(i,j),ij2nv(i,j)) = 1/2;
       AyvS(ij2nv(i,j),ij2nv(i,j-1)) = 1/2;
       #--east face
       if (maskv(ij2nv(i+1,j)) != WALL)
         AyvE(ij2nv(i,j),ij2nv(i,j)) = 1/2;
         AyvE(ij2nv(i,j),ij2nv(i+1,j)) = 1/2;
       endif
       #--west face
       if ((maskv(ij2nv(i-1,j)) != WALL) &&
           (maskv(ij2nv(i-1,j)) != INFLOW))
         AyvW(ij2nv(i,j),ij2nv(i,j)) = 1/2;
         AyvW(ij2nv(i,j),ij2nv(i-1,j)) = 1/2;
       endif
     ######################################
     #### FACE VELOCITIES (TANGENCIAL)
     ####   uy{E,W}     = Ayu{E,W}*u
       #--east face
       if (maskv(ij2nv(i+1,j)) != WALL)
         AyuE(ij2nv(i,j),ij2nu(i+1,j)) = 1/2;
         AyuE(ij2nv(i,j),ij2nu(i+1,j-1)) = 1/2;
       endif
       #--west face
       if (maskv(ij2nv(i-1,j)) != WALL)
         AyuW(ij2nv(i,j),ij2nu(i,j-1)) = 1/2;
         AyuW(ij2nv(i,j),ij2nu(i,j)) = 1/2;
       endif
     endif
     ######################################
     #### VISCOUS LAPLACIAN
     if (maskv(ij2nv(i,j)) == INNER)
       #-- east face
       if (maskv(ij2nv(i+1,j)) == INNER)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i+1,j)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         #ghost point
         Dy(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
         #do nothing
       elseif (maskv(ij2nv(i+1,j)) == OUTFLOW)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
       #-- west face
       if (maskv(ij2nv(i-1,j)) == INNER)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i-1,j)) =  1/h/h;
       elseif (maskv(ij2nv(i-1,j)) == WALL)
         #ghost point
         Dy(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
         #do nothing
       elseif (maskv(ij2nv(i-1,j)) == INFLOW)
         #ghost point
         Dy(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
         #do nothing
       endif
       #-- north face
       if (maskv(ij2nv(i,j+1)) == INNER)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i,j+1)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
       #-- south face
       if (maskv(ij2nv(i,j-1)) == INNER)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i,j-1)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
     else
       #not inner: left unchanged
       Dy(ij2nv(i,j),ij2nv(i,j)) = 1;
     endif
     ######################################
     #### PRESSURE GRADIENT
     if (maskv(ij2nv(i,j)) == INNER)
       gradPV(ij2nv(i,j),ij2n(i,j)) = 1/h;
       gradPV(ij2nv(i,j),ij2n(i,j-1)) = -1/h;
     endif
   endfor
endfor

### END INITIALIZATION
####################################################################

####################################################################
### BEGIN TIME INTEGRATION

p = zeros(dimP,1);
u = zeros(dimU,1);
v = zeros(dimV,1);

u = ubc;

i = 0;
t = T0;
while ((t < T1)&&(i < NT))

  #advection operator
  uE = AxuE*u; uW = AxuW*u;
  uN = AxuN*u; uS = AxuS*u;
  vN = AxvN*v; vS = AxvS*v;
  Ax = (uE.*uE - uW.*uW + uN.*vN - uS.*vS)/h/h;

  vN = AyvN*v; vS = AyvS*v;
  vE = AyvE*v; vW = AyvW*v;
  uE = AyuE*u; uW = AyuW*u;
  Ay = (vN.*vN - vS.*vS + vE.*uE - vW.*uW)/h/h;

  #prediction step
  ustar = u + 1*dt*(-Ax + nu*Dx*u);
  vstar = v + 1*dt*(-Ay + nu*Dy*v);

  #pressure poisson equation
  rhs = (rho/dt)*(gradU*ustar + gradV*vstar) - Apu*ubc - Apv*vbc;
  pnew = Ap\rhs;

  #correction step
  unew = ustar - (dt/rho)*gradPU*pnew;
  vnew = vstar - (dt/rho)*gradPV*pnew;

  #finish timestep
  u = unew;
  v = vnew;
  p = pnew;
  t = t + dt;
  i = i + 1;
  fprintf("t= %12.5E  |u|= %12.5E |v|= %12.5E |p|= %12.5E\n",
          t,norm(u),norm(v),norm(p));
  fflush(stdout);
endwhile

### END INITIALIZATION
####################################################################

####################################################################
### BEGIN POST-PROCESSING

XU=linspace(0,(NI+1)*dx,NI+1);
YU=linspace(dy/2,NJ*dy,NJ);
XV=linspace(dx/2,NI*dx,NI);
YV=linspace(0,(NJ+1)*dy,NJ+1);
XP=linspace(dx/2,NI*dx,NI);
YP=linspace(dy/2,NJ*dy,NJ);

U=reshape(u,NJ,NI+1);
V=reshape(v,NJ+1,NI);
P=reshape(p,NJ,NI);
M=reshape(mask,NJ,NI);

#surf(XP,YP,P);
contourf(XP,YP,P);
colorbar();
title(sprintf("Pressure at T=%g",t));
figure();
#surf(XU,YU,U);
contourf(XU,YU,U);
colorbar();
title(sprintf("U-Velocity at T=%g",t));
figure();
#surf(XV,YV,V);
contourf(XV,YV,V);
colorbar();
title(sprintf("V-Velocity at T=%g",t));

### END POST-PROCESSING
####################################################################
