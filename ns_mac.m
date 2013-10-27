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
dt = 0.00001;
T0 = 0;
T1 = 1;
NT = 1;

#Arrays
dimP = NI*NJ;
dimU = (NI+1)*NJ;
dimV = NI*(NJ+1);

Ap = spalloc(dimP,dimP,5*dimP);
Dx = spalloc(dimU,dimU,5*dimU);
Dy = spalloc(dimV,dimV,5*dimV);

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
     if (i == 1)
       mask(ij2n(i,j))=INFLOW;
     endif
     if (i == NI)
       mask(ij2n(i,j))=OUTFLOW;
     endif
     if ((j == 1)||(j == NJ))
       mask(ij2n(i,j))=WALL;
     endif
     if ((i <= NI/3)&&(j <= NJ/2))
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


#P-cells
disp("assembling pressure laplacian and velocity gradients ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ
     ######################################
     #### PRESSURE LAPLACIAN
     if (mask(ij2n(i,j)) == INNER)
       # east face
       if ((mask(ij2n(i+1,j)) == INNER) ||
           (mask(ij2n(i+1,j)) == OUTFLOW))
         Ap(ij2n(i,j),ij2n(i,j)) += -1/h/h;
         Ap(ij2n(i,j),ij2n(i+1,j))=  1/h/h;
       endif
       #west face
       if (mask(ij2n(i-1,j)) == INNER)
         Ap(ij2n(i,j),ij2n(i,j)) += -1/h/h;
         Ap(ij2n(i,j),ij2n(i-1,j))=  1/h/h;
       endif
       #north face
       if (mask(ij2n(i,j+1)) == INNER)
         Ap(ij2n(i,j),ij2n(i,j)) += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j+1))=  1/h/h;
       endif
       #south face
       if (mask(ij2n(i,j-1)) == INNER)
         Ap(ij2n(i,j),ij2n(i,j)) += -1/h/h;
         Ap(ij2n(i,j),ij2n(i,j-1))=  1/h/h;
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
disp("assembling u-face veloc and viscous laplacian ..."); fflush(stdout);
for i=1:NI+1
   for j=1:NJ
     ######################################
     #### FACE VELOCITIES
     ####   ux{E,W,N,S} = Axu{E,W,N,S}*u
     ####   vx{N,S}     = Axv{N,S}*v
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
       if (masku(ij2nu(i,j+1)) == INNER)
         AxuN(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuN(ij2nu(i,j),ij2nu(i,j+1)) = 1/2;
         AxvN(ij2nu(i,j),ij2nv(i,j+1)) = 1/2;
         AxvN(ij2nu(i,j),ij2nv(i-1,j+1)) = 1/2;
       endif
       #--south face
       if (masku(ij2nu(i,j-1)) == INNER)
         AxuS(ij2nu(i,j),ij2nu(i,j)) = 1/2;
         AxuS(ij2nu(i,j),ij2nu(i,j-1)) = 1/2;
         AxvS(ij2nu(i,j),ij2nv(i-1,j)) = 1/2;
         AxvS(ij2nu(i,j),ij2nv(i,j)) = 1/2;
       endif
     endif
     ######################################
     #### VISCOUS LAPLACIAN
     if (masku(ij2nu(i,j)) == INNER)
       Dx(ij2nu(i,j),ij2nu(i,j)) = -4/h/h;
       #-- east face
       if (masku(ij2nu(i+1,j)) == OUTFLOW)
         #zero gradient
         Dx(ij2nu(i,j),ij2nu(i,j)) +=  1/h/h;
       elseif (masku(ij2nu(i+1,j)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j)) += -1/h/h;
       else
         Dx(ij2nu(i,j),ij2nu(i+1,j)) +=  1/h/h;
       endif
       #-- west face
       if (masku(ij2nu(i-1,j)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j)) += -1/h/h;
       else
         Dx(ij2nu(i,j),ij2nu(i-1,j)) +=  1/h/h;
       endif
       #-- north face
       if (masku(ij2nu(i,j+1)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j)) += -1/h/h;
       else
         Dx(ij2nu(i,j),ij2nu(i,j+1)) +=  1/h/h;
       endif
       #-- south face
       if (masku(ij2nu(i,j-1)) == WALL)
         Dx(ij2nu(i,j),ij2nu(i,j)) += -1/h/h;
       else
         Dx(ij2nu(i,j),ij2nu(i,j-1)) +=  1/h/h;
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
disp("assembling v-face veloc and viscous laplacian ..."); fflush(stdout);
for i=1:NI
   for j=1:NJ+1
     ######################################
     #### ADVECTION AVERAGES 
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
         AyuE(ij2nv(i,j),ij2nu(i+1,j)) = 1/2;
         AyuE(ij2nv(i,j),ij2nu(i+1,j-1)) = 1/2;
       endif
       #--west face
       if (maskv(ij2nv(i+1,j)) != WALL)
         AyvW(ij2nv(i,j),ij2nv(i,j)) = 1/2;
         AyvW(ij2nv(i,j),ij2nv(i-1,j)) = 1/2;
         AyuW(ij2nv(i,j),ij2nu(i,j-1)) = 1/2;
         AyuW(ij2nv(i,j),ij2nu(i,j)) = 1/2;
       endif
     endif
     ######################################
     #### VISCOUS LAPLACIAN
     if (maskv(ij2nv(i,j)) == INNER)
       Dy(ij2nv(i,j),ij2nv(i,j)) = -4/h/h;
       #-- east face
       if (maskv(ij2nv(i+1,j)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j)) += -1/h/h;
       else
         Dy(ij2nv(i,j),ij2nv(i+1,j)) +=  1/h/h;
       endif
       #-- west face
       if (maskv(ij2nv(i-1,j)) == INFLOW)
         Dy(ij2nv(i,j),ij2nv(i,j)) += -1/h/h;
       elseif (maskv(ij2nv(i-1,j)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j)) += -1/h/h;
       else
         Dy(ij2nv(i,j),ij2nv(i-1,j)) +=  1/h/h;
       endif
       #-- north face
       if (maskv(ij2nv(i,j+1)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j)) += -1/h/h;
       else
         Dy(ij2nv(i,j),ij2nv(i,j+1)) +=  1/h/h;
       endif
       #-- south face
       if (maskv(ij2nv(i,j-1)) == WALL)
         Dy(ij2nv(i,j),ij2nv(i,j)) += -1/h/h;
       else
         Dy(ij2nv(i,j),ij2nv(i,j-1)) +=  1/h/h;
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

#Inflow
N = (NJ-1)-(NJ/2+1)+1
UI = parabola(N);
for k=1:N
  j=k+(NJ/2);
  u(ij2nu(2,j)) = UI(k);
endfor

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
  ustar = u + 0*dt*(-Ax + nu*Dx*u);
  vstar = v + 0*dt*(-Ay + nu*Dy*v);

  #pressure poisson equation
  rhs = (rho/dt) * (gradU*ustar + gradV*vstar)
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

surf(XP,YP,P);
figure();
surf(XU,YU,U);
figure();
surf(XV,YV,V);

### END POST-PROCESSING
####################################################################
