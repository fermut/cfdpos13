# Finite Volume Solver for the convection-reaction-diffusion equation on a
# two-dimensional rectangle of dimensions Lx and Ly.
# The integral (conservation) equation being solve is as follows:
#
# /                      /                                /
# |  dU                  |                                | 
# | (-- + r u) dOmega  + | (-k grad(U) + U V).nG dGamma = | f dOmega
# |  dt                  |                                |
# / Omega                / Gamma                          / Omega
#
# where U(x,y,t) is the conserved scalar (unknown), r is the reaction coeff.,
# k(x,y) is the diffusion coeff., V is the velocity, nG is the outgoing normal,
# and f(x,y) the source term. For this problem, V = (vx, 0) is assumed.
# A theta-scheme for time discretization leads to the following system:
#
#  ((1/dt+theta*r)*Am + theta*Af).W(n+1) = 
#                     ((1/dt-(1-theta)*r)*Am - (1-theta)*Af)W(n) + Bf
#
# U(x,y) is prescribed to zero (Dirichlet boundary conditions) at the four
# walls of the cross section: ((0,y),(Lx,y),(x,0),(x,Ly))
#
# Author: Fernando Mut && Gustavo Buscaglia
# Last Modification: 09/09/2013

#BEGIN PROBLEM DEFINITION
#---------------------------------------------------------------------------
#Convention used
# NJ   x  x  x  ... x 
# .    .
# .    .
# .    . 
# 3    x  x  x  ... x  
# 2    x  x  x  ... x  
# 1    x  x  x  ... x  
#
# j/i  1  2  3  ... NI

# Global index if(NI < NJ) n = ij2n(i,j) = i + (j-1)*NI
# Global index if(NI > NJ) n = ij2n(i,j) = j + (i-1)*NJ

global NI = 100;   #X Direction
global NJ = 50;    #Y Direction
Lx = 2.0;
Ly = 1.0;
reac = 0.;         #Reaction
diff = 0.000025;   #Diffusion
gamma = 1.;        #0<=g<=1, 1: full upwinding, 0: centered
vx = 1.;
#uniform grid spacing
dx = (Lx/(NI))*ones(NI,1);
dy = (Ly/(NJ))*ones(NJ,1);
#if dx or dy are variable, program the required values here
#dx( 1: 40) = 0.015;
#dx(41: 60) = 0.040;
#dx(61:100) = 0.015;
#dy( 1: 20) = 0.015;
#dy(21: 30) = 0.040;
#dy(31: 50) = 0.015;
#uniform diffusion
k = diff*ones(NI,NJ);
#if difussion is variable, program the required values here
#k(NI/2:NI,:)=3*diff;
#uniform source
fs = 1*ones(NI*NJ,1);
#if difussion is variable, program the required values here
for i=1:NI
   for j=1:NJ
       if (j > NJ/2)
          fs(ij2n(i,j),1) = 3.;
       endif
       if (i > NI/2)
          fs(ij2n(i,j),1) = 0.;
       endif
    endfor
endfor
#temporal integration
theta = 0.5;
dt = .05;
T0 = 0.;
T1 = 4.;
#plotting
NP = 1;                 #timesteps increment for plots
NT = (T1-T0)/dt+1;      #number of stored timesteps
#working space
nunk = NI*NJ;
fprintf(stdout, "  * PROBLEM DEFINITION DONE -> Unknowns = %d\n", nunk);
fflush(stdout);
#matrices and arrays
Af = zeros(nunk,nunk);   # flux matrix
Am = zeros(nunk,nunk);   # mass matrix
bm = zeros(nunk,1);      # mass RHS
#arrays for diagnostic purposes
U0 = zeros(nunk,1);      # initial condition
UN = zeros(nunk,NT/NP);  # stored solutions
TN = zeros(NT/NP);       # times for stored solutions
diag = zeros(10,NT);     # diagnostics
#---------------------------------------------------------------------------
#END PROBLEM DEFINITION

#BEGIN PRE-PROCESSING STEP
#---------------------------------------------------------------------------
#-- Assembly: loop over control volumes
for i=1:NI
   for j=1:NJ
#-- Mass Matrix and Mass RHS: Am and Bm
      Am(ij2n(i,j),ij2n(i,j))=dx(i)*dy(j);
      bm(ij2n(i,j))=dx(i)*dy(j);
#-- Flux Matrix: Af
# east face
      if (i<NI)
      kface=(k(i,j)*dx(i+1)+k(i+1,j)*dx(i))/(dx(i)+dx(i+1));
      aux1 = kface*2*dy(j)/(dx(i)+dx(i+1)) + gamma*vx*dy(j);
      aux2 = kface*2*dy(j)/(dx(i)+dx(i+1)) - (1-gamma)*vx*dy(j);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux1;
      Af(ij2n(i,j),ij2n(i+1,j))=Af(ij2n(i,j),ij2n(i+1,j))-aux2;
      else
      kface=k(i,j);
      aux = 2*kface*dy(j)/dx(i) + max(2*(gamma-0.5),0.)*vx*dy(j);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      endif
# west face
      if (i>1)
      kface=(k(i,j)*dx(i-1)+k(i-1,j)*dx(i))/(dx(i)+dx(i-1));
      aux1 = kface*2*dy(j)/(dx(i)+dx(i-1)) - (1-gamma)*vx*dy(j);
      aux2 = kface*2*dy(j)/(dx(i)+dx(i-1)) + gamma*vx*dy(j);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux1;
      Af(ij2n(i,j),ij2n(i-1,j))=Af(ij2n(i,j),ij2n(i-1,j))-aux2;
      else
      kface=k(i,j);
      aux = 2*kface*dy(j)/dx(i);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      endif
# north face
      if (j<NJ)
      kface=(k(i,j)*dy(j+1)+k(i,j+1)*dy(j))/(dy(j)+dy(j+1));
      aux = kface*2*dx(i)/(dy(j)+dy(j+1));
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      Af(ij2n(i,j),ij2n(i,j+1))=Af(ij2n(i,j),ij2n(i,j+1))-aux;
      else
      kface=k(i,j);
      aux = 2*kface*dx(i)/dy(j);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      endif
# south face
      if (j>1)
      kface=(k(i,j)*dy(j-1)+k(i,j-1)*dy(j))/(dy(j)+dy(j-1));
      aux = kface*2*dx(i)/(dy(j)+dy(j-1));
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      Af(ij2n(i,j),ij2n(i,j-1))=Af(ij2n(i,j),ij2n(i,j-1))-aux;
      else
      kface=k(i,j);
      aux = 2*kface*dx(i)/dy(j);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      endif
   endfor
endfor

#-- Timestepping Matrices: Ai, Ae
Ai = ((1/dt)+theta*reac)*Am + theta*Af;
Ae = ((1/dt)-(1-theta)*reac)*Am - (1-theta)*Af;

#-- Switch to sparse storage ?
flag_sparse = 1;
fprintf(stdout, "  * CREATING FINAL MATRICES -> flag_sparse = %d\n",
        flag_sparse);
fflush(stdout);

if(flag_sparse == 1) 
  Ai = sparse(Ai);
  Ae = sparse(Ae);
endif
#---------------------------------------------------------------------------
#END PRE-PROCESSING STEP

#BEGIN SOLUTION STEP
#---------------------------------------------------------------------------
#-- Temporal integration loop
nt = 1;
np = 0;
tn = T0;
un = U0;
dnorm = 0;
fprintf(stdout, "     SOLVING ...\n");
exit_flag = 0;
while 1
   #-- Diagnostics

   #store 'some' solutions for visualization
   if rem(nt-1,NP) == 0
      np=np+1;
      TN(np)=tn;
      UN(:,np)=un;
   endif

   #time(1)
   diag(1,nt)=tn;

   #concentration(2)
   diag(2,nt)=bm'*un;

   #screen output
   unorm=norm(un);
   fprintf(stdout, "t= %12.5E |u|= %12.5E |du|= %12.5E\n", tn, unorm, dnorm);
   fflush(stdout);

   #should exit ?
   if exit_flag > 0
     break;
   endif 

   #-- Advance in time
   tm=tn+dt;

   #assemble RHS
   b=Ae*un + bm.*fs;

   #solve linear system
   um = Ai\b;

   #delta
   dm = (um.-un)/dt;
   dnorm=norm(dm);

   #end of integration ?
   if tn >= T1
     exit_flag = 2;
   endif

   #end timestep
   un = um;
   tn = tm;
   nt = nt + 1;
endwhile
#---------------------------------------------------------------------------
#END SOLUTION STEP

#BEGIN POST-PROCESSING STEP
#---------------------------------------------------------------------------
#-- Integrals
plot(diag(1,:),diag(2,:))
title("Concentration")
figure()

#-- Meshgrid: points are placed in the center of each cell
X=zeros(NI,1);
X(1)=dx(1)/2;
for i=2:NI
  DX=0.5*(dx(i-1)+dx(i));
  X(i)=X(i-1)+DX;
endfor
Y=zeros(NJ,1);
Y(1)=dy(1)/2;
for j=2:NJ
  DY=0.5*(dy(j-1)+dy(j));
  Y(j)=Y(j-1)+DY;
endfor

#-- Limits
xmin=min(X);
xmax=max(X);
ymin=min(Y);
ymax=max(Y);
umin=min(UN(:,np));
umax=max(UN(:,np));
for i=1:np
  umin = min(umin,min(UN(:,i)));
  umax = max(umax,max(UN(:,i)));
endfor

#-- Plot Solution 
for i=1:np
  U = reshape(UN(:,i),NJ,NI);
  surf(X,Y,U);
  axis([xmin xmax ymin ymax umin umax]);
  title(sprintf("T=%g",TN(i)))
  #figure()
  #contourf(X,Y,U);
  pause();
endfor
#---------------------------------------------------------------------------
#END POST-PROCESSING STEP
