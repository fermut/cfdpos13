# Finite Volume Solver for a fully-developed flow inside a rectangular
# pipe of dimensions Lx and Ly (UNSTEADY VERSION).
# The integral (conservation) equation being solve is as follows:
#
# /                  /                        /
# |     dW           |                        |
# | rho -- dOmega  - | mu grad(W).nG dGamma = | g dOmega
# |     dt           |                        |
# / Omega            / Gamma                  / Omega
#
# where W(x,y,t) is the velocity along z (logitudinal direction), 
# rho is the density, nG is the outgoing normal, mu(x,y) is the viscosity
# and g is the resulting (unknown) pressure gradient.
# The equation for the flow rate reads as follows:
#
# /
# |
# | W dOmega = Q
# |
# / Omega
#
# where Q(t) is the given flow rate.
#
# A theta-scheme for time discretization leads to the following system:
#
#  (1/dt*Am + theta*Af + Aq).W(n+1) = (1/dt*Am - (1-theta)*Af)W(n) = Bq(n+1)
#
#
# W(x,y) is prescribed to zero (Dirichlet boundary conditions) at the four
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

global NI = 100;  #X Direction
global NJ = 50;   #X Direction
Lx = 2.0;
Ly = 1.0;
rho = 1.;
visc = 1.;
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
#uniform viscosity
mu = visc*ones(NI,NJ);
#if viscosity is variable, program the required values here
#mu(NI/2:NI,:)=3*visc;
#temporal integration
theta = 1.0;
dt = 0.01;
T0 = 0.;
T1 = 1.;
#plotting
IP = 1;                 #timesteps increment for plots
NP = (T1-T0)/dt+1;      #number of stored timesteps
#working space
nunk = NI*NJ+1;
fprintf(stdout, "  * PROBLEM DEFINITION DONE -> Unknowns = %d\n", nunk);
fflush(stdout);
#matrices and arrays
Aq = zeros(nunk,nunk);   # flow matrix
Af = zeros(nunk,nunk);   # flux matrix
Am = zeros(nunk,nunk);   # mass matrix
bm = zeros(nunk,1);      # mass RHS
bq = zeros(nunk,1);      # flow RHS
#arrays for diagnostic purposes
W0 = zeros(nunk,1);      # initial condition
WN = zeros(nunk,NP/IP);  # stored solutions
TN = zeros(NP/IP);       # times for stored solutions
diag = zeros(10,NP);     # diagnostics
#---------------------------------------------------------------------------
#END PROBLEM DEFINITION

#BEGIN PRE-PROCESSING STEP
#---------------------------------------------------------------------------
#-- Assembly: loop over control volumes
for i=1:NI
   for j=1:NJ
#-- Mass Matrix and RHS: Am and Bm
      Am(ij2n(i,j),ij2n(i,j))=dx(i)*dy(j);
      bm(ij2n(i,j))=dx(i)*dy(j);
#-- Flow Matrix: Aq
      Aq(ij2n(i,j),NI*NJ+1)=dx(i)*dy(j);
#-- Flux Matrix: Af
# east face
      if (i<NI)
      muface=(mu(i,j)*dx(i+1)+mu(i+1,j)*dx(i))/(dx(i)+dx(i+1));
      aux = muface*2*dy(j)/(dx(i)+dx(i+1));
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      Af(ij2n(i,j),ij2n(i+1,j))=Af(ij2n(i,j),ij2n(i+1,j))-aux;
      else
      muface=mu(i,j);
      aux = muface*dy(j)/dx(i);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+2*aux;
      endif
# west face
      if (i>1)
      muface=(mu(i,j)*dx(i-1)+mu(i-1,j)*dx(i))/(dx(i)+dx(i-1));
      aux = muface*2*dy(j)/(dx(i)+dx(i-1));
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      Af(ij2n(i,j),ij2n(i-1,j))=Af(ij2n(i,j),ij2n(i-1,j))-aux;
      else
      muface=mu(i,j);
      aux = muface*dy(j)/dx(i);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+2*aux;
      endif
# north face
      if (j<NJ)
      muface=(mu(i,j)*dy(j+1)+mu(i,j+1)*dy(j))/(dy(j)+dy(j+1));
      aux = muface*2*dx(i)/(dy(j)+dy(j+1));
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      Af(ij2n(i,j),ij2n(i,j+1))=Af(ij2n(i,j),ij2n(i,j+1))-aux;
      else
      muface=mu(i,j);
      aux = muface*dx(i)/dy(j);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+2*aux;
      endif
# south face
      if (j>1)
      muface=(mu(i,j)*dy(j-1)+mu(i,j-1)*dy(j))/(dy(j)+dy(j-1));
      aux = muface*2*dx(i)/(dy(j)+dy(j-1));
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+aux;
      Af(ij2n(i,j),ij2n(i,j-1))=Af(ij2n(i,j),ij2n(i,j-1))-aux;
      else
      muface=mu(i,j);
      aux = muface*dx(i)/dy(j);
      Af(ij2n(i,j),ij2n(i,j))=Af(ij2n(i,j),ij2n(i,j))+2*aux;
      endif
   endfor
endfor

#-- Timestepping Matrices: Ai, Ae
Ai = (rho/dt)*Am + theta*(Af - Aq) + Aq';
Ae = (rho/dt)*Am - (1-theta)*(Af - Aq);

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
wn = qflow(tn)*bm/(Lx*Ly);
printf("Qini= %g %g\n",qflow(tn),bm'*wn);
anorm = 0;
fprintf(stdout, "     SOLVING ...\n");
exit_flag = 0;
while 1
   #-- Diagnostics

   #store 'some' solutions for visualization
   if rem(nt-1,IP) == 0
      np=np+1;
      TN(np)=tn;
      WN(:,np)=wn;
   endif

   #time(1)
   diag(1,nt)=tn;

   #kinetic energy(2)
   diag(2,nt)=0.5*rho*(bm.*wn)'*wn;

   #flow rate(3)
   diag(3,nt)=bm'*wn;

   #dp/dz(4)
   diag(4,nt)=wn(NI*NJ+1);

   #screen output
   wnorm=norm(wn);
   fprintf(stdout, "t= %12.5E |w|= %12.5E |a|= %12.5E\n", tn, wnorm, anorm);
   fflush(stdout);

   #should exit ?
   if exit_flag > 0
     break;
   endif 

   #-- Advance in time
   tm=tn+dt;

   #assemble RHS
   bq(NI*NJ+1)=qflow(tm);
   b=Ae*wn + bq;

   #solve linear system
   wm = Ai\b;

   #acceleration
   am = (wm.-wn)/dt;
   anorm=norm(am);

   #end of integration ?
   if tn >= T1
     exit_flag = 1;
   endif

   #end timestep
   wn = wm;
   tn = tm;
   nt = nt + 1;
endwhile
#---------------------------------------------------------------------------
#END SOLUTION STEP

#BEGIN POST-PROCESSING STEP
#---------------------------------------------------------------------------
#-- Integrals
plot(diag(1,:),diag(2,:))
title("Kinetic Energy")
figure()
plot(diag(1,:),diag(3,:))
title("Flow Rate")
figure()
plot(diag(1,:),diag(4,:))
title("Pressure Gradient")
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
wmin=min(WN(1:NI*NJ,np));
wmax=max(WN(1:NI*NJ,np));
for i=1,np
  wmin = min(wmin,min(WN(1:NI*NJ,i)));
  wmax = max(wmax,max(WN(1:NI*NJ,i)));
endfor

#-- Plot Solution 
for i=1:np
  U = reshape(WN(1:NI*NJ,i),NJ,NI);
  surf(X,Y,U);
  axis([xmin xmax ymin ymax wmin wmax]);
  title(sprintf("T=%g",TN(i)))
  #figure()
  #contourf(X,Y,U);
  pause();
endfor
#---------------------------------------------------------------------------
#END POST-PROCESSING STEP
