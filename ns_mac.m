#Problem definition
global NI=10*12*2+1
global NJ=10* 2*2+2
Lx = 12.0
Ly = 2.0
dx = Lx/(NI-1)
dy = Ly/(NJ-2)
h = (dx+dy)/2
umax = 1

#Material properties
nu = 1/100;
rho= 1;

#Temporal integration
dtd = h*h/(nu*4)
dtc = 2*nu/umax
dt = 0.5*min(dtd,dtc)
T0 = 0;
T1 = 1000;
NT = (T1-T0)/dt;
dtol = 1.0E-6;

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

avgU = spalloc(dimP,dimU,2*dimP);
avgV = spalloc(dimP,dimV,2*dimP);

#Mask definition
INTERIOR=0;
WALL=1;
INFLOW=2;
OUTFLOW=3;
mask=INTERIOR*ones(dimP,1);
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
     ###INTERIOR BOUNDARY
     if (1)
     if ((i <= (NI-1)/3+1)&&(j <= (NJ-2)/2+1))
       mask(ij2n(i,j))=WALL;
     endif
     endif
   endfor
endfor

#Mask transcription to U-,V-cells
masku=INTERIOR*ones(dimU);
maskv=INTERIOR*ones(dimV);
disp("transcriving mask to u,v cells ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     if (mask(ij2n(i,j)) != INTERIOR)
       masku(ij2nu(i,j))  =mask(ij2n(i,j)); 
       masku(ij2nu(i+1,j))=mask(ij2n(i,j)); 
       maskv(ij2nv(i,j))  =mask(ij2n(i,j)); 
       maskv(ij2nv(i,j+1))=mask(ij2n(i,j)); 
     endif
     if (mask(ij2n(i,j)) == OUTFLOW)
       masku(ij2nu(i,j)) = INTERIOR;
     endif
   endfor
endfor

#Boundary Conditions
disp("generating boundary conditions ..."); fflush(stdout);
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
    ubc(ij2nu(2,j)) = umax*UI(N);
    N = N + 1;
  endif
endfor

######################################
#### PRESSURE LAPLACIAN
disp("pressure laplacian ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     if (mask(ij2n(i,j)) == INTERIOR)
       #-- east face
       if (masku(ij2nu(i+1,j)) == INTERIOR)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i+1,j)) =  1/h/h;
       else
         Apu(ij2n(i,j),ij2nu(i+1,j)) = 1;
       endif
       #-- west face
       if (masku(ij2nu(i,j)) == INTERIOR)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i-1,j)) =  1/h/h;
       else
         Apu(ij2n(i,j),ij2nu(i,j)) = -1;
       endif
       #-- north face
       if (maskv(ij2nv(i,j+1)) == INTERIOR)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j+1)) =  1/h/h;
       else
         Apv(ij2n(i,j),ij2nv(i,j+1)) = 1;
       endif
       #-- south face
       if (maskv(ij2nv(i,j)) == INTERIOR)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j-1)) =  1/h/h;
       else
         Apv(ij2n(i,j),ij2nv(i,j)) = -1;
       endif
     else
       #not inner: left unchanged
       Ap(ij2n(i,j),ij2n(i,j)) = 1;
     endif
   endfor
endfor

######################################
#### FACE VELOCITIES (U-Cells)
disp("face velocities (u-cells) ..."); fflush(stdout);
for i=1:NI+1
   for j=1:NJ
     ######################################
     #### NORMAL VELOCITIES
     ####   ux{E,W,N,S} = Axu{E,W,N,S}*u
     if (masku(ij2nu(i,j)) == INTERIOR)
       #-- east face
       if (masku(ij2nu(i+1,j)) == INTERIOR)
         AxuE(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuE(ij2nu(i,j),ij2nu(i+1,j)) = 1/2;
       elseif (masku(ij2nu(i+1,j)) == WALL)
         AxuE(ij2nu(i,j),ij2nu(i,j)) = 1/2;
       elseif (masku(ij2nu(i+1,j)) == OUTFLOW)
         # zero gradient
         AxuE(ij2nu(i,j),ij2nu(i,j)) = 1;  
       endif
       #-- west face
       if (masku(ij2nu(i-1,j)) == INTERIOR)
         AxuW(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuW(ij2nu(i,j),ij2nu(i-1,j)) = 1/2;
       elseif (masku(ij2nu(i-1,j)) == WALL)
         AxuE(ij2nu(i,j),ij2nu(i,j)) = 1/2;
       elseif (masku(ij2nu(i-1,j)) == INFLOW)
         AxuW(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuW(ij2nu(i,j),ij2nu(i-1,j)) = 1/2;
       endif
       #--north face
       if (masku(ij2nu(i,j+1)) == INTERIOR)
         AxuN(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuN(ij2nu(i,j),ij2nu(i,j+1)) = 1/2;
       elseif (masku(ij2nu(i,j+1)) == WALL)
         #do nothing
       endif
       #--south face
       if (masku(ij2nu(i,j-1)) == INTERIOR)
         AxuS(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuS(ij2nu(i,j),ij2nu(i,j-1)) = 1/2;
       elseif (masku(ij2nu(i,j-1)) == WALL)
         #do nothing
       endif
     ######################################
     #### TANGENCIAL VELOCITIES
     ####   vx{N,S}     = Axv{N,S}*v
       #--north face
       if (masku(ij2nu(i,j+1)) == INTERIOR)
         AxvN(ij2nu(i,j),ij2nv(i,j+1)) = 1/2;
         AxvN(ij2nu(i,j),ij2nv(i-1,j+1)) = 1/2;
       elseif (masku(ij2nu(i,j+1)) == WALL)
         #do nothing
       endif
       #--south face
       if (masku(ij2nu(i,j-1)) == INTERIOR)
         AxvS(ij2nu(i,j),ij2nv(i-1,j)) = 1/2;
         AxvS(ij2nu(i,j),ij2nv(i,j)) = 1/2;
       elseif (masku(ij2nu(i,j-1)) == WALL)
         #do nothing
       endif
     endif
   endfor
endfor

######################################
#### FACE VELOCITIES (V-Cells)
disp("face velocities (v-cells) ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ+1
     ######################################
     #### NORMAL VELOCITIES
     ####   vy{E,W,N,S} = Ayv{E,W,N,S}*v
     if (maskv(ij2nv(i,j)) == INTERIOR)
       #--east face
       if (maskv(ij2nv(i+1,j)) == INTERIOR)
         AyvE(ij2nv(i,j),ij2nv(i,j)) = 1/2;
         AyvE(ij2nv(i,j),ij2nv(i+1,j)) = 1/2;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         #do nothing
       elseif (maskv(ij2nv(i+1,j)) == OUTFLOW)
         AyvE(ij2nv(i,j),ij2nv(i,j)) = 1/2;
       endif
       #--west face
       if (maskv(ij2nv(i-1,j)) == INTERIOR)
         AyvW(ij2nv(i,j),ij2nv(i,j)) = 1/2;
         AyvW(ij2nv(i,j),ij2nv(i-1,j)) = 1/2;
       elseif (maskv(ij2nv(i-1,j)) == WALL)
         #do nothing
       elseif (maskv(ij2nv(i-1,j)) == INFLOW)
         #do nothing
       endif
       #--north face
       if (maskv(ij2nv(i,j+1)) == INTERIOR)
         AyvN(ij2nv(i,j),ij2nv(i,j)) = 1/2;
         AyvN(ij2nv(i,j),ij2nv(i,j+1)) = 1/2;
       elseif (maskv(ij2nv(i,j+1)) == WALL)
         AyvN(ij2nv(i,j),ij2nv(i,j)) = 1/2;
       endif
       #--south face
       if (maskv(ij2nv(i,j-1)) == INTERIOR)
         AyvS(ij2nv(i,j),ij2nv(i,j)) = 1/2;
         AyvS(ij2nv(i,j),ij2nv(i,j-1)) = 1/2;
       elseif (maskv(ij2nv(i,j-1)) == WALL)
       endif
     ######################################
     #### TANGENCIAL VELOCITIES
     ####   uy{E,W}     = Ayu{E,W}*u
       #--east face
       if (maskv(ij2nv(i+1,j)) == INTERIOR)
         AyuE(ij2nv(i,j),ij2nu(i+1,j)) = 1/2;
         AyuE(ij2nv(i,j),ij2nu(i+1,j-1)) = 1/2;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         #do nothing
       elseif (maskv(ij2nv(i+1,j)) == OUTFLOW)
         AyuE(ij2nv(i,j),ij2nu(i+1,j)) = 1/2;
         AyuE(ij2nv(i,j),ij2nu(i+1,j-1)) = 1/2;
       endif
       #--west face
       if (maskv(ij2nv(i-1,j)) == INTERIOR)
         AyuW(ij2nv(i,j),ij2nu(i,j)) = 1/2;
         AyuW(ij2nv(i,j),ij2nu(i,j-1)) = 1/2;
       elseif (maskv(ij2nv(i-1,j)) == WALL)
         #do nothing
       elseif (maskv(ij2nv(i-1,j)) == INFLOW)
         AyuW(ij2nv(i,j),ij2nu(i,j)) = 1/2;
         AyuW(ij2nv(i,j),ij2nu(i,j-1)) = 1/2;
       endif
     endif
   endfor
endfor

######################################
#### VISCOUS LAPLACIAN (U-Cells)
disp("viscous laplacian (u-cells) ..."); fflush(stdout);
for i=1:NI+1
   for j=1:NJ
     if (masku(ij2nu(i,j)) == INTERIOR)
       #-- east face
       if (masku(ij2nu(i+1,j)) == INTERIOR)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i+1,j)) =  1/h/h;
       elseif (masku(ij2nu(i+1,j)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
       elseif (masku(ij2nu(i+1,j)) == OUTFLOW)
         #zero gradient
         #do nothing
       endif
       #-- west face
       if (masku(ij2nu(i-1,j)) == INTERIOR)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i-1,j)) =  1/h/h;
       elseif (masku(ij2nu(i-1,j)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
       elseif (masku(ij2nu(i-1,j)) == INFLOW)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i-1,j)) =  1/h/h;
       endif
       #-- north face
       if (masku(ij2nu(i,j+1)) == INTERIOR)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i,j+1)) =  1/h/h;
       elseif (masku(ij2nu(i,j+1)) == WALL)
         #ghost point
         Dx(ij2nu(i,j),ij2nu(i,j))  += -2/h/h;
       endif
       #-- south face
       if (masku(ij2nu(i,j-1)) == INTERIOR)
         Dx(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dx(ij2nu(i,j),ij2nu(i,j-1)) =  1/h/h;
       elseif (masku(ij2nu(i,j-1)) == WALL)
         #ghost point
         Dx(ij2nu(i,j),ij2nu(i,j))  += -2/h/h;
       endif
     else
       #not inner: left unchanged
       Dx(ij2nu(i,j),ij2nu(i,j)) = 0;
     endif
   endfor
endfor

######################################
#### VISCOUS LAPLACIAN (V-Cells)
disp("viscous laplacian (v-cells) ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ+1
     if (maskv(ij2nv(i,j)) == INTERIOR)
       #-- east face
       if (maskv(ij2nv(i+1,j)) == INTERIOR)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i+1,j)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         #ghost point
         Dy(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
       elseif (maskv(ij2nv(i+1,j)) == OUTFLOW)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
       #-- west face
       if (maskv(ij2nv(i-1,j)) == INTERIOR)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i-1,j)) =  1/h/h;
       elseif (maskv(ij2nv(i-1,j)) == WALL)
         #ghost point
         Dy(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
       elseif (maskv(ij2nv(i-1,j)) == INFLOW)
         #ghost point
         Dy(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
       endif
       #-- north face
       if (maskv(ij2nv(i,j+1)) == INTERIOR)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i,j+1)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
       #-- south face
       if (maskv(ij2nv(i,j-1)) == INTERIOR)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dy(ij2nv(i,j),ij2nv(i,j-1)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
     else
       #not inner: left unchanged
       Dy(ij2nv(i,j),ij2nv(i,j)) = 0;
     endif
   endfor
endfor

######################################
#### VELOCITY GRADIENTS (P-Cells)
disp("velocity gradients ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     if (mask(ij2n(i,j)) == INTERIOR)
       gradU(ij2n(i,j),ij2nu(i+1,j)) = 1/h;
       gradU(ij2n(i,j),ij2nu(i,j)) = -1/h;
       gradV(ij2n(i,j),ij2nv(i,j+1)) = 1/h;
       gradV(ij2n(i,j),ij2nv(i,j)) = -1/h;
     endif
   endfor
endfor

######################################
#### PRESSURE GRADIENT (U-Cells)
disp("pressure gradient (u-cells) ..."); fflush(stdout);
for i=1:NI+1
   for j=1:NJ
     if (masku(ij2nu(i,j)) == INTERIOR)
       gradPU(ij2nu(i,j),ij2n(i,j)) = 1/h;
       gradPU(ij2nu(i,j),ij2n(i-1,j)) = -1/h;
     endif
   endfor
endfor

######################################
#### PRESSURE GRADIENT (V-Cells)
disp("pressure gradient (v-cells) ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ+1
     if (maskv(ij2nv(i,j)) == INTERIOR)
       gradPV(ij2nv(i,j),ij2n(i,j)) = 1/h;
       gradPV(ij2nv(i,j),ij2n(i,j-1)) = -1/h;
     endif
   endfor
endfor

######################################
#### VELOCITY AVERAGES (for visualization)
disp("velocity averages ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     if (mask(ij2n(i,j)) == INTERIOR)
       avgU(ij2n(i,j),ij2nu(i,j)) = 1/2;
       avgU(ij2n(i,j),ij2nu(i+1,j)) = 1/2;
       avgV(ij2n(i,j),ij2nv(i,j)) = 1/2;
       avgV(ij2n(i,j),ij2nv(i,j+1)) = 1/2;
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
  Ax = (uE.*uE - uW.*uW + uN.*vN - uS.*vS)/h;

  vN = AyvN*v; vS = AyvS*v;
  vE = AyvE*v; vW = AyvW*v;
  uE = AyuE*u; uW = AyuW*u;
  Ay = (vN.*vN - vS.*vS + vE.*uE - vW.*uW)/h;

  #prediction step
  ustar = u + dt*(-Ax + nu*Dx*u);
  vstar = v + dt*(-Ay + nu*Dy*v);

  #pressure poisson equation
  rhs = (rho/dt)*(gradU*ustar + gradV*vstar) + Apu*ubc + Apv*vbc;
  pnew = Ap\rhs;

  #projection step
  unew = ustar - (dt/rho)*gradPU*pnew;
  vnew = vstar - (dt/rho)*gradPV*pnew;

  #deltas
  dp=norm(pnew-p);
  du=norm(unew-u);
  dv=norm(vnew-v);

  #finish timestep
  u = unew;
  v = vnew;
  p = pnew;
  t = t + dt;
  i = i + 1;
  if (mod(i,100) == 0)
  fprintf("t= %12.5E  |u|= %12.5E |v|= %12.5E |p|= %12.5E\n",
          t,norm(u),norm(v),norm(p));
  fflush(stdout);
  endif

  #convergence ?
  if ((dp < dtol)&&(du < dtol)&&(dv < dtol))
    break;
  endif
endwhile

### END INITIALIZATION
####################################################################

####################################################################
### BEGIN POST-PROCESSING

XP=linspace(dx/2,(NI-0.5)*dx,NI);
YP=linspace(dy/2,(NJ-0.5)*dy,NJ);

XU=linspace(0,NI*dx,NI+1);
YU=linspace(dy/2,(NJ-0.5)*dy,NJ);

XV=linspace(dx/2,(NI-0.5)*dx,NI);
YV=linspace(0,NJ*dy,NJ+1);

P=reshape(p,NJ,NI);
U=reshape(u,NJ,NI+1);
V=reshape(v,NJ+1,NI);

#surf(XP,YP,P);
contourf(XP(2:NI),YP(2:NJ-1),P(2:NJ-1,2:NI));
#axis([XP(2) XP(NI) YP(2) YP(NJ-1)])
axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
colorbar();
title(sprintf("Pressure at T=%g",t));
figure();

#surf(XU,YU,U);
contourf(XU(2:NI),YU(2:NJ-1),U(2:NJ-1,2:NI));
#axis([XU(2) XU(NI) YU(2) YU(NJ-1)])
hold
czero=contour(XU(2:NI),YU(2:NJ-1),U(2:NJ-1,2:NI),[0,0]);
axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
colorbar();
xzero=czero(end-1)-dx;
title(sprintf("U-Velocity at T=%g, Vortex Attachment Point at X=%g",t,xzero));
figure();

#surf(XV,YV,V);
contourf(XV(2:NI),YV(2:NJ),V(2:NJ,2:NI));
#axis([XV(2) XV(NI) YV(2) YV(NJ)])
axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
colorbar();
title(sprintf("V-Velocity at T=%g",t));
figure();

up = avgU*u;
vp = avgV*v;
UP=reshape(up,NJ,NI);
VP=reshape(vp,NJ,NI);

h=quiver(XP(2:NI),YP(2:NJ-1),UP(2:NJ-1,2:NI),VP(2:NJ-1,2:NI));
axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
set (h, "maxheadsize", 0.33);

### END POST-PROCESSING
####################################################################
