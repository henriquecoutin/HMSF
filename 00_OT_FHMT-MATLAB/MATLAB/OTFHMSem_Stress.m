%% TOPOLOGY OPTIMIZATION WITH HYBRID-MIXED STRESS ELEMENT %% 
%  BY LUIS ARMANDO, AUGUST 2017 % 
function TOHMSF (nelx,nely,volfrac,penal,rmin,ft)
% INITIALIZE
x(1:nely,1:nelx) = volfrac;
change = 1.; c= 0.; loop = 0; dv = ones(nely,nelx);
[H,Hs]   = check(nelx,nely,rmin);  
while change > 0.01  
xold = x; cold = c; loop = loop + 1;
% HMSF-ANALYSIS   
[STRESS] = FHMT(nelx,nely,x,penal);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [~,~,~,~,~,~,~,~,~,FEL] = EL;
  c = 0.; dc = zeros(nely,nelx);
  for ely = 1:nely
    for elx = 1:nelx
        n1 = elx + (nelx+1)*(ely-1);
        n2 = elx + (nelx+1)*ely;
        S = STRESS([n1*3-2; n1*3-1; n1*3; n1*3+1; n1*3+2; n1*3+3;
            n2*3+1; n2*3+2; n2*3+3; n2*3-2; n2*3-1; n2*3],1);
      c = c + x(ely,elx)^(-penal)*S'*FEL*S;
      dc(ely,elx) = penal*x(ely,elx)^(-penal-1)*S'*FEL*S;
    end
  end
% FILTERING 
if ft == 1 % ON SENSITIVITIES
    dc(:) = H*(x(:).*dc(:))./Hs;
elseif ft == 2 % ON DENSITY 
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
end
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
[x] = OC(nelx,nely,x,volfrac,dc,dv,ft,H,Hs);
change = max(max(abs(x-xold))); change = min(change,abs(cold-c));
% PRINT RESULTS
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES  
  colormap(gray); imagesc(-flip(x)); axis equal; axis tight; axis off;pause(1e-6);
end 
end
% ELEMENT MATRICES
function [iF,jF,sF,iAe,jAe,sAe,iAg,jAg,sAg,FEL] = EL
F = 1e-8*[6250 -1875 0; -1875 6250 0; 0 0 16250];
FEL = [F F F F; F F F F; F F F F; F F F F]; [iF,jF,sF] = find(FEL);
AEL = [-1/8 0 -1/8 0 -1/8 0 -1/8 0; 0 -1/8 0 -1/8 0 -1/8 0 -1/8;
    -1/8 -1/8 -1/8 -1/8 -1/8 -1/8 -1/8 -1/8; 1/8 0 1/8 0 1/8 0 1/8 0;
    0 -1/8 0 -1/8 0 -1/8 0 -1/8; -1/8 1/8 -1/8 1/8 -1/8 1/8 -1/8 1/8;
    -1/8 0 -1/8 0 -1/8 0 -1/8 0; 0 1/8 0 1/8 0 1/8 0 1/8;
    1/8 -1/8 1/8 -1/8 1/8 -1/8 1/8 -1/8; 1/8 0 1/8 0 1/8 0 1/8 0;
    0 1/8 0 1/8 0 1/8 0 1/8; 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8]; [iAe,jAe,sAe] = find(AEL);
AGEL = [-1/3 0 0 0 -1/6 0 0 0; 0 -1/3 0 -1/6 0 0 0 0;
    -1/3 -1/3 -1/6 0 0 -1/6 0 0; 0 0 1/3 0 0 0 1/6 0;
    0 -1/6 0 -1/3 0 0 0 0; -1/6 0 -1/3 1/3 0 0 0 1/6;
    -1/6 0 0 0 -1/3 0 0 0; 0 0 0 0 0 1/3 0 1/6;
    0 -1/6 0 0 1/3 -1/3 1/6 0; 0 0 1/6 0 0 0 1/3 0;
    0 0 0 0 0 1/6 0 1/3; 0 0 0 1/6 1/6 0 1/3 1/3]; [iAg,jAg,sAg] = find(AGEL);
end
% HMSF ANALYSIS 
function [STRESS] = FHMT(nelx,nely,x,penal)
[iF,jF,sF,iAe,jAe,sAe,iAg,jAg,sAg,~] = EL;
ntrip = 0; ndof = (nelx+1)*(nely+1)*7; nodes = (nelx+1)*(nely+1);
I = zeros((nelx)*(nely)*(size(iF,1)+size(iAg,1)+size(iAg,1)),1);
J = zeros((nelx)*(nely)*(size(iF,1)+size(iAg,1)+size(iAg,1)),1);
S = zeros((nelx)*(nely)*(size(iF,1)+size(iAg,1)+size(iAg,1)),1);
for ely = 1:nely
    for elx = 1:nelx
        n1 = elx + (nelx+1)*(ely-1);
        n2 = elx + (nelx+1)*ely;
        fdof = [n1*3-2 n1*3-1 n1*3 n1*3+1 n1*3+2 n1*3+3 ...
            n2*3-2 n2*3-1 n2*3 n2*3+1 n2*3+2 n2*3+3];
        dof = [n1*2-1 n1*2 n1*2+1 n1*2+2 n2*2-1 n2*2 n2*2+1 n2*2+2];  
        aedof = dof + nodes*3; agdof = dof + nodes*5;
        for fi = 1:size(iF,1)
            ntrip = ntrip + 1;
            I(ntrip) = fdof(iF(fi)); J(ntrip) = fdof(jF(fi)); S(ntrip) = x(ely,elx)^(-penal)*sF(fi);
        end
        for ei = 1:size(iAe,1)
            ntrip = ntrip + 1;
            I(ntrip) = fdof(iAe(ei)); J(ntrip) = aedof(jAe(ei)); S(ntrip) = sAe(ei);
            ntrip = ntrip + 1;
            I(ntrip) = aedof(jAe(ei)); J(ntrip) = fdof(iAe(ei)); S(ntrip) = sAe(ei);
        end
        for gi = 1:size(iAg,1)
            ntrip = ntrip + 1;
            I(ntrip) = fdof(iAg(gi)); J(ntrip) = agdof(jAg(gi)); S(ntrip) = -sAg(gi);
            ntrip = ntrip + 1;
            I(ntrip) = agdof(jAg(gi)); J(ntrip) = fdof(iAg(gi)); S(ntrip) = -sAg(gi);
        end
    end
end
K = sparse(I,J,S,ndof,ndof);  
% APPLIED LOAD
FI =  5*nodes+(nelx+1)*(nely/2+1)*2; FJ = 1; FS = -10; F = sparse(FI,FJ,-FS,ndof,1);
% BOUNDARY CONDITION
nfixed = 1:nelx+1:(nelx+1)*nely+1;
doffixed = sort([2*nfixed 2*nfixed-1]) + 5*nodes;
% MATRIX MODIFICATION FOR SOLVING SYSTEM
iB = zeros(3*nodes,1); jB = zeros(3*nodes,1); sB = zeros(3*nodes,1); 
for bi = 1:3*nodes
    iB(bi) = bi; jB(bi) = bi; sB(bi) = K(bi,bi)^ (-1/2);
end
B = sparse(iB,jB,sB,3*nodes,3*nodes);
K(:,1:3*nodes) = K(:,1:3*nodes)*B;
K(1:3*nodes,:) = B*K(1:3*nodes,:);
K = K + 10^-8*speye(7*nodes);

K(doffixed,:) = 0; K(doffixed,doffixed) = eye(max(size(doffixed)));
F(doffixed) = 0;
% MULTIFRONTAL METHOD FOR SOLVING LINEAR SYSTEM
struct = ma57_factor(K);
U = ma57_solve(K, full(F), struct); 
STRESS = B*U(1:3*nodes);
end
% H/Hs MATRIX DEFINITION FOR FILTER APPLICATION        
function [H,Hs]=check(nelx,nely,rmin)
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nely
    for j1 = 1:nelx
        e1 = (j1-1)*nely+i1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nely)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nelx)
                e2 = (j2-1)*nely+i2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
end
% OPTIMALITY CRITERIA METHOD 
function [xnew]=OC(nelx,nely,x,volfrac,dc,dv,ft,H,Hs)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(dc./dv/lmid)))));
  if ft == 2
      xnew(:) =  (H*xnew(:))./Hs;
  end
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
end