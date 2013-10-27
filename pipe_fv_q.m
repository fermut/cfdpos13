# Finite Volume Solver for a fully-developed flow inside a rectangular
# pipe of dimensions Lx and Ly (STEADY VERSION).
# The integral (conservation) equation being solve is as follows:
#
#   /                          /
#   |                          | dp
# - | mu grad(W).nG dGamma = - | -- dOmega     =>     Af.W = bp
#   |                          | dz
#   / Gamma                    / Omega
#
# where W(x,y) is the velocity along z (logitudinal direction), 
# nG is the outgoing normal, mu(x,y) is the viscosity and dp/dz is the
# applied pressure gradient. The LHS (including the minus sign, so Af is
# defined positive) is discretized into the matrix Af (flux matrix).
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
global NJ = 50;   #Y Direction
Lx = 2.0;
Ly = 1.0;
visc = 1.;
qf = 1.5;
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
#working space
nunk = NI*NJ+1;
fprintf(stdout, "  * PROBLEM DEFINITION DONE -> Unknowns = %d\n", nunk);
fflush(stdout);
#matrices and arrays
Af = zeros(nunk,nunk);   # flux matrix
Aq = zeros(nunk,nunk);   # flow matrix
bq = zeros(nunk,1);      # RHS
#---------------------------------------------------------------------------
#END PROBLEM DEFINITION

#BEGIN PRE-PROCESSING STEP
#---------------------------------------------------------------------------
#-- Assembly: loop over control volumes
for i=1:NI
   for j=1:NJ
#-- Flow Matrix: Aw
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

#-- Final System: A,b
A = Af+Aq+Aq';
bq(NI*NJ+1) = qf;

#-- Switch to sparse storage
flag_sparse = 1;
fprintf(stdout, "  * CREATING FINAL MATRICES -> flag_sparse = %d\n",
        flag_sparse);
fflush(stdout);

if (flag_sparse == 1) 
  A = sparse(A);
endif
#---------------------------------------------------------------------------
#END PRE-PROCESSING STEP

#BEGIN SOLUTION STEP
#---------------------------------------------------------------------------
fprintf(stdout, "     SOLVING ...");
fflush(stdout);
t0 = time;

#-- Solve linear system
W = A\bq;
tsol = time - t0;

fprintf(stdout, "\n      -> Solution time = %f\n", tsol);
fflush(stdout);
#---------------------------------------------------------------------------
#END SOLUTION STEP

#BEGIN POST-PROCESSING STEP
#---------------------------------------------------------------------------
#-- Splitting
G = W(NI*NJ+1);
W = W(1:NI*NJ);
printf("dpdt= %g\n",G);
#-- Integrals
bm=Aq(1:NI*NJ,NI*NJ+1);
printf("FlowRate= %g\n",bm'*W);
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
#-- Plot solution 
U = reshape(W, NJ, NI);
surf(X,Y,U);
#figure()
#contourf(X,Y,U);
#---------------------------------------------------------------------------
#END POST-PROCESSING STEP
