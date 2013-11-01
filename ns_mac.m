#Problem definition
global NI NJ
Lx = 12.0
Ly = 2.0
h = 0.25
NI = floor(Lx/h)
NJ = floor(Ly/h)
dx = Lx/NI
dy = Ly/NJ
h = (dx+dy)/2
NI = NI + 1       #left layer
NJ = NJ + 2       #top and bottom layer
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

#Postprocessing
NPRIN = 10;
NSAVE = 500;

#Arrays
dimP = NI*NJ
dimU = (NI+1)*NJ
dimV = NI*(NJ+1)

Ap = spalloc(dimP,dimP,5*dimP);
Dxu = spalloc(dimU,dimU,5*dimU);
Dyv = spalloc(dimV,dimV,5*dimV);

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
       mask(ij2n(i,j)) = INFLOW;
     endif
     if (i == NI)
       mask(ij2n(i,j)) = OUTFLOW;
     endif
     if ((j == 1)||(j == NJ))
       mask(ij2n(i,j)) = WALL;
     endif
     ###INTERIOR WALLS
     xp = (i-1.5)*dx;
     yp = (j-1.5)*dy;
     if ((xp < 4)&&(yp < 1))
       mask(ij2n(i,j)) = WALL;
     endif
   endfor
endfor

#Mask transcription to U-,V-cells
masku=INTERIOR*ones(dimU);
maskv=INTERIOR*ones(dimV);
disp("transcribing mask to u,v cells ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     if (mask(ij2n(i,j)) != INTERIOR)
       masku(ij2nu(i,j))   = mask(ij2n(i,j)); 
       masku(ij2nu(i+1,j)) = mask(ij2n(i,j)); 
       maskv(ij2nv(i,j))   = mask(ij2n(i,j)); 
       maskv(ij2nv(i,j+1)) = mask(ij2n(i,j)); 
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
       endif
       #-- west face
       if (masku(ij2nu(i,j)) == INTERIOR)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i-1,j)) =  1/h/h;
       endif
       #-- north face
       if (maskv(ij2nv(i,j+1)) == INTERIOR)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j+1)) =  1/h/h;
       endif
       #-- south face
       if (maskv(ij2nv(i,j)) == INTERIOR)
         Ap(ij2n(i,j),ij2n(i,j))  += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j-1)) =  1/h/h;
       endif
     else
       #not interior: left unchanged
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
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dxu(ij2nu(i,j),ij2nu(i+1,j)) =  1/h/h;
       elseif (masku(ij2nu(i+1,j)) == WALL)
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
       elseif (masku(ij2nu(i+1,j)) == OUTFLOW)
         #zero gradient
         #do nothing
       endif
       #-- west face
       if (masku(ij2nu(i-1,j)) == INTERIOR)
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dxu(ij2nu(i,j),ij2nu(i-1,j)) =  1/h/h;
       elseif (masku(ij2nu(i-1,j)) == WALL)
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
       elseif (masku(ij2nu(i-1,j)) == INFLOW)
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dxu(ij2nu(i,j),ij2nu(i-1,j)) =  1/h/h;
       endif
       #-- north face
       if (masku(ij2nu(i,j+1)) == INTERIOR)
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dxu(ij2nu(i,j),ij2nu(i,j+1)) =  1/h/h;
       elseif (masku(ij2nu(i,j+1)) == WALL)
         #ghost point
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -2/h/h;
       endif
       #-- south face
       if (masku(ij2nu(i,j-1)) == INTERIOR)
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -1/h/h;
         Dxu(ij2nu(i,j),ij2nu(i,j-1)) =  1/h/h;
       elseif (masku(ij2nu(i,j-1)) == WALL)
         #ghost point
         Dxu(ij2nu(i,j),ij2nu(i,j))  += -2/h/h;
       endif
     else
       #not interior: left unchanged
       Dxu(ij2nu(i,j),ij2nu(i,j)) = 1;
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
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dyv(ij2nv(i,j),ij2nv(i+1,j)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         #ghost point
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
       elseif (maskv(ij2nv(i+1,j)) == OUTFLOW)
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
       #-- west face
       if (maskv(ij2nv(i-1,j)) == INTERIOR)
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dyv(ij2nv(i,j),ij2nv(i-1,j)) =  1/h/h;
       elseif (maskv(ij2nv(i-1,j)) == WALL)
         #ghost point
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
       elseif (maskv(ij2nv(i-1,j)) == INFLOW)
         #ghost point
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -2/h/h;
       endif
       #-- north face
       if (maskv(ij2nv(i,j+1)) == INTERIOR)
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dyv(ij2nv(i,j),ij2nv(i,j+1)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
       #-- south face
       if (maskv(ij2nv(i,j-1)) == INTERIOR)
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
         Dyv(ij2nv(i,j),ij2nv(i,j-1)) =  1/h/h;
       elseif (maskv(ij2nv(i+1,j)) == WALL)
         Dyv(ij2nv(i,j),ij2nv(i,j))  += -1/h/h;
       endif
     else
       #not interior: left unchanged
       Dyv(ij2nv(i,j),ij2nv(i,j)) = 1;
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
### POST-PROCESSING INITIALIZATION
XP=linspace(dx/2-dx,(NI-0.5)*dx-dx,NI);
YP=linspace(dy/2-dy,(NJ-0.5)*dy-dy,NJ);

XU=linspace(0-dx,NI*dx-dx,NI+1);
YU=linspace(dy/2-dy,(NJ-0.5)*dy-dy,NJ);

XV=linspace(dx/2-dx,(NI-0.5)*dx-dx,NI);
YV=linspace(0-dy,NJ*dy-dy,NJ+1);

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

  #viscous operator
  Dx = Dxu*u-ubc;
  Dy = Dyv*v-vbc;

  #prediction step
  ustar = u + dt*(-Ax + nu*Dx);
  vstar = v + dt*(-Ay + nu*Dy);

  #pressure poisson equation
  rhs = (rho/dt)*(gradU*ustar + gradV*vstar);
  pnew = Ap\rhs;

  #correction step
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

  #output ?
  if (mod(i,NPRIN) == 0)
    fprintf("t= %12.5E  |u|= %12.5E |v|= %12.5E |p|= %12.5E\n",
            t,norm(u),norm(v),norm(p));
    fflush(stdout);
  endif

  #save image ?
  if (mod(i,NSAVE) == 0)
    fprintf("saving solution at T= %g\n", t);
    fflush(stdout);

    P=reshape(p,NJ,NI);
    U=reshape(u,NJ,NI+1);
    V=reshape(v,NJ+1,NI);

    figure(1,'visible','off');
    clf();
    contourf(XP(2:NI),YP(2:NJ-1),P(2:NJ-1,2:NI));
    axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
    colorbar();
    title(sprintf("Pressure at T=%g",t));
    print(sprintf("press%08d.png",i),"-dpng");
    close();

    figure(2,'visible','off');
    clf();
    contourf(XU(2:NI),YU(2:NJ-1),U(2:NJ-1,2:NI));
    hold
    czero=contour(XU(2:NI),YU(2:NJ-1),U(2:NJ-1,2:NI),[0,0]);
    axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
    colorbar();
    xzero=czero(end-1)-dx;
    title(sprintf("U-Velocity at T=%g, Vortex Attachment Point at X=%g",
                  t,xzero));
    print(sprintf("velou%08d.png",i),"-dpng");
    close();

    figure(3,'visible','off');
    clf();
    contourf(XV(2:NI),YV(2:NJ),V(2:NJ,2:NI));
    axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
    colorbar();
    title(sprintf("V-Velocity at T=%g",t));
    print(sprintf("velov%08d.png",i),"-dpng");
    close();
  endif

  #convergence ?
  if ((dp < dtol)&&(du < dtol)&&(dv < dtol))
    break;
  endif
endwhile

####################################################################
### POST-PROCESSING

figure();
contourf(XP(2:NI),YP(2:NJ-1),P(2:NJ-1,2:NI));
axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
colorbar();
title(sprintf("Pressure at T=%g",t));

figure();
contourf(XU(2:NI),YU(2:NJ-1),U(2:NJ-1,2:NI));
hold
czero=contour(XU(2:NI),YU(2:NJ-1),U(2:NJ-1,2:NI),[0,0]);
axis([XU(1) XU(NI+1) YV(1) YV(NJ+1)])
colorbar();
xzero=czero(end-1)-dx;
title(sprintf("U-Velocity at T=%g, Vortex Attachment Point at X=%g",
              t,xzero));

figure();
contourf(XV(2:NI),YV(2:NJ),V(2:NJ,2:NI));
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
