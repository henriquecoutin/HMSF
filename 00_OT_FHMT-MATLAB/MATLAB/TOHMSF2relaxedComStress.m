%% TOPOLOGY OPTIMIZATION WITH HYBRID-MIXED STRESS ELEMENT %% 
%  BY LUIS ARMANDO, AUGUST 2017 % 
function TOHMSF2relaxed ()
% INITIALIZE
nelx = 80; nely = 80; volfrac = 0.3; 
tx = 10; ty = 10; dx = tx/nelx; dy = ty/nely;
penal = 3.0; rmin = 0.02*nelx; ft = 2;
x(1:nely,1:nelx) = volfrac;
change = 1.; loop = 0;vol = ones(nely,nelx);
energy1 = 1000000000;
[H,Hs] = check(nelx,nely,rmin);  
[I,J,S,F,local,matrix_stress] = PREK(nelx,nely,dx,dy);
[IJS,Msize] = dkdp(nelx,nely); %%%%%%%%%%%%%%%%%%% ATENCAO
low = ones(nelx*nely,1); upp = low; % MMA
%%%%% L-Beam %%%% 
index = 0;
for j = 0.4*nely:nely
    for i = 0.4*nelx:nelx
        index = index + 1;
        passive(index) = (j-1)*nelx+i;
        x(j,i) = 0.001;
    end
end
all = 1:1:nelx*nely; free = setdiff(all,passive);
xold1 = reshape(x',nelx*nely,1); xold2 = reshape(x',nelx*nely,1);
xmma = xold1;
 while change > 0.01 && loop < 101  
xold = x; loop = loop + 1;
% HMSF-ANALYSIS   
[STRESS,K,U,struct,B] = FHMT(nelx,nely,x,penal,I,J,S,F,local);
% GLOBAL STRESS
[stress_el,von_mises] = EL_STRESS(STRESS,matrix_stress,nelx,nely,x);
[adjoint,sum_stress,second,third,mises] = ADJOINT(x,struct,K,von_mises,stress_el,matrix_stress,nelx,nely,B,passive);
[IJS,Msize] = dkdp(nelx,nely);
[dgdp] = SENSITIVITIES(adjoint,U,nelx,nely,x,penal,Msize,IJS);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [~,~,~,~,~,~,~,~,~,FEL] = EL(dx,dy);
  energy = 0; dc = zeros(nely,nelx);
  for ely = 1:nely
    for elx = 1:nelx
        n1 = elx + (nelx+1)*(ely-1);
        n2 = elx + (nelx+1)*ely;
        St = STRESS([n1*3-2; n1*3-1; n1*3; n1*3+1; n1*3+2; n1*3+3;
            n2*3+1; n2*3+2; n2*3+3; n2*3-2; n2*3-1; n2*3],1);
      energy = energy + x(ely,elx)^(-penal)*St'*FEL*St;
      dc(ely,elx) = -penal*x(ely,elx)^(-penal-1)*St'*FEL*St;
    end
  end
% FILTERING 
dv = ones(nely,nelx);
if ft == 1 % ON SENSITIVITIES
    dc(:) = H*(x(:).*dc(:))./Hs;
elseif ft == 2 % ON DENSITY 
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
end
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
%[x] = OC(nelx,nely,x,volfrac,dc,dv,ft,H,Hs);
% DESIGN UPDATE BY MMA 
 m  = 2; % GENERAL CONSTRAINTS 
 n = size(free,2);% NUMBER OF VARIABLES 
 xval = reshape(x',nelx*nely,1);  % COLUMN VECTOR CURRENT VALUES OF X
 xmin = 0.001*ones(nelx*nely,1); xmax = ones(nelx*nely,1);
 fval(1,1) = energy/30-1; %VALUE OF OBJECTIVE FUNCTION
 fval(2,1) = sum_stress-1; 
 dfdx = dc'; dfdx = dfdx(:)'/30;
 %dvp = flip(dv); dvp = dvp';
 dfdx(2,:) = dgdp'+second+third; %SENSITIVITIES OF OBJECTIVE FUNCITON
 dfdx(2,mises<0)=0;
 f0val = sum(xval); % /(nelx*nely*volfrac)-1; %CONSTRAINT VALUE AT ITER
 df0dx = reshape(dv',nelx*nely,1); %'/(nelx*nely*volfrac); %
 a0 = 1; a = zeros(m,1); c = 1000*ones(m,1); d = zeros(m,1);
 [xmma(free),~,~,~,~,~,~,~,~,low(free),upp(free)] = ...
 mmasub(m,n,loop,xval(free),xmin(free),xmax(free),xold1(free),xold2(free), ...
 f0val,df0dx(free),fval,dfdx(:,free),low(free),upp(free),a0,a,c,d);
 x = reshape(xmma,nelx,nely)';
 x(:) = H*(x(:)./Hs);
 xold2 = xold1; xold1 = reshape(x',nelx*nely,1);

change = max(max(abs(x-xold)));
change = min(change,abs(energy1-energy));
energy1 = energy;
% PRINT RESULTS
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
       ' stress.: ' sprintf('%6.3f',fval(2,1)) ...
       ' compl.: ' sprintf('%6.3f',fval(1,1)) ...
        ' ch.: ' sprintf('%6.3f',change)])
% PLOT DENSITIES  
colormap(gray); imagesc(-flip(x)); axis equal; axis tight; axis off;pause(1e-6);

 end

%%%%%%%%%%
% STRESS %
%%%%%%%%%%
[~,von_mises] = EL_STRESS(STRESS,matrix_stress,nelx,nely,x);
%von_mises = reshape(von_mises,nelx,nely)';
k = 0;
for i = 1:nely
    for j = 1:nelx
        k = k + 1;
        X(k) = j - 0.5;
        Y(k) = i - 0.5;
    end
end
fileID = fopen('stress.txt','w');
fprintf(fileID,'%6.2f %6.2f %12.8f\n',[X; Y; von_mises']);
fclose(fileID);
end

% ELEMENT MATRICES
function [iF,jF,sF,iAe,jAe,sAe,iAg,jAg,sAg,FEL] = EL(dx,dy)
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
function [STRESS,K,U,struct,B] = FHMT(nelx,nely,x,penal,I,J,S,F,local)

ndof = (nelx+1)*(nely+1)*7; nodes = (nelx+1)*(nely+1);
xtemp = 0.;
xtemp = reshape(x',1,nelx*nely); xtemp = repmat(xtemp,80,1);
xtemp = reshape(xtemp,80*nelx*nely,1);
Stemp = S;
Stemp(local) = Stemp(local).*xtemp.^(-penal);

K = sparse(I,J,Stemp,ndof,ndof);  
% BOUNDARY CONDITION
% nfixed = 1:nelx+1:(nelx+1)*nely+1; % left
nfixed = (nelx+1)*nely:(nelx+1)*(nely+1);
doffixed = sort([2*nfixed 2*nfixed-1]) + 5*nodes; % LEFT
% doffixed = sort([2*nfixed-1 2*(nelx+1)]) + 5*nodes; % MBB
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
U(1:3*nodes) = B*U(1:3*nodes);
%U = zeros(7*nodes,1);
%U = babuska(K,F,U);
STRESS = U(1:3*nodes);
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
% BABUSKA METHOD
function [x] = babuska(K,F,U) 
        j = 0;
        x = U;
        r = F; tol = 10^-3; err = 1.0;
        den = transpose(r)*r;
        while err > tol
            a = K\r;
            x = x + a;
            r = r - K*a; num = transpose(r)*r;
            j = j + 1;
            err = sqrt(num/den);
        end
end
% PRE-ASSEMBLY OF COEFFICIENT MATRIX
function [I,J,S,F,local,matrix_stress] = PREK(nelx,nely,dx,dy)
[iF,jF,sF,iAe,jAe,sAe,iAg,jAg,sAg,~] = EL(dx,dy);
ntrip = 0; ndof = (nelx+1)*(nely+1)*7; nodes = (nelx+1)*(nely+1);
I = zeros((nelx)*(nely)*(size(iF,1)+size(iAg,1)+size(iAg,1)),1);
J = zeros((nelx)*(nely)*(size(iF,1)+size(iAg,1)+size(iAg,1)),1);
S = zeros((nelx)*(nely)*(size(iF,1)+size(iAg,1)+size(iAg,1)),1);
local = zeros(80*nelx*nely,1); trip = 0; 
count = 0; It = zeros(12*nelx*nely,1); Jt = zeros(12*nelx*nely,1);
Sst = zeros(12*nelx*nely,1); cout = 0;
for ely = 1:nely
    for elx = 1:nelx
        n1 = elx + (nelx+1)*(ely-1); n2 = elx + (nelx+1)*ely;
        fdof = [n1*3-2 n1*3-1 n1*3 n1*3+1 n1*3+2 n1*3+3 ...
            n2*3-2 n2*3-1 n2*3 n2*3+1 n2*3+2 n2*3+3];
        dof = [n1*2-1 n1*2 n1*2+1 n1*2+2 n2*2-1 n2*2 n2*2+1 n2*2+2];
        aedof = dof + nodes*3; agdof = dof + nodes*5;
        
        It(count+1:count+4) = cout+1; Jt(count+1:count+4) = [n1*3-2; n1*3+1; n2*3-2; n2*3+1];
        It(count+5:count+8) = cout+2; Jt(count+5:count+8) = [n1*3-1; n1*3+2; n2*3-1; n2*3+2];
        It(count+9:count+12) = cout+3; Jt(count+9:count+12) = [n1*3; n1*3+3; n2*3; n2*3+3];
        Sst(count+1:count+12) = 1;
        count = count + 12; cout = cout + 3;
        
        for fi = 1:size(iF,1)
            ntrip = ntrip + 1; trip = trip + 1;
            I(ntrip) = fdof(iF(fi)); J(ntrip) = fdof(jF(fi)); S(ntrip) = sF(fi);
            local(trip) = ntrip;
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
% APPLIED LOAD % MIDDLE
% FI =  5*nodes+(nelx+1)*(nely/2+1)*2; FJ = 1; FS = -10; F = sparse(FI,FJ,-FS,ndof,1);
% MBB
% FI =  5*nodes+(nelx+1)*(nely)*2+2; FJ = 1; FS = -10; F = sparse(FI,FJ,-FS,ndof,1);
% L-BEAM
FI =  5*nodes+(nelx+1)*(0.4*nely/2+1)*2; FJ = 1; FS = -10; F = sparse(FI,FJ,-FS,ndof,1);

% STRESS MATRIX CALCULATOR
matrix_stress = sparse(It,Jt,Sst,3*nelx*nely,nodes*3);
end
% ELEMENT STRESSES AND VON MISES STRESS 
function [stress_el,von_mises] = EL_STRESS(STRESS,matrix_stress,nelx,nely,x)
stress_el = 1/4*matrix_stress*STRESS; 
stress_el = reshape(stress_el,3,nelx*nely);
stress_el = stress_el';
pho = reshape(x',nelx*nely,1);
% for i = 1:nelx*nely
%     stress_el(i,:) = stress_el(i,:)*pho(i)^0.5; 
% end
T = [1 -1/2 0; -1/2 1 0; 0 0 3];
von_mises = dot(stress_el*T,stress_el,2);
von_mises = sqrt(von_mises);
end
% STRESS CONSTRAINT SENSITIVITY 
function [dtalndpe,sum_stress] = CONSTRAINT_SENS(stress_sens,matrix_stress,von_mises,stress_el,nelx,nely)
der = 1/4*matrix_stress*stress_sens; 
der = reshape(der,3,nelx*nely);
der = der';
der(:,1) = dxdpe; der(:,2) = dydpe; der(:,3) = dtaldpe;

p = 6; % P-NORM FACTOR 
sum_stress = (sum(von_mises.^p))^(1/p-1);

temp = von_mises'.^(p-2)*1/2*((2*repmat(stress_el(:,1),[1,nelx*nely])-repmat(stress_el(:,2),[1,nelx*nely])).*dxdpe +...
   (2*repmat(stress_el(:,2),[1,nelx*nely])-repmat(stress_el(:,1),[1,nelx*nely])).*dydpe +...
   6*repmat(stress_el(:,3),[1,nelx*nely]).*dtaldpe);

dtalndpe = sum_stress*temp;
end
% ADJOINT EQUATION 
function [adjoint,sum_stress,second,third,mises] = ADJOINT(x,struct,K,von_mises,stress_el,matrix_stress,nelx,nely,B,passive)
% DERIVATE P-NORM GLOBAL STRESS CONSTRAINT (IN SIGMA)
%von_mises(abs(von_mises)<1e-3) = 0;
%stress_el(abs(stress_el)<1e-3) = 0;
N = nelx*nely;
q = 2.5; % Stress Relax Factor
sl = 5.5; % Stress Limit
epsi = 0.01;
p = 100; % P-NORM FACTOR 
pho = reshape(x',nelx*nely,1);
von_mises(passive) = 0.;
vm = von_mises;
sum_stress = (sum(von_mises+epsi-epsi./pho).^p)^(1/p-1);
von_mises = dot(von_mises,pho.^(-q),2)/sl; % Writing Gvm/(Gl*x^q)
for ij = 1:nelx*nely
    mises(ij,1) = von_mises(ij) + epsi - epsi/pho(ij);
end
mises(mises<0)=0.;
sum_stress = (sum(mises).^p)^(1/p-1);
von_mises(mises<0)=0.;
rows = 1:3:nelx*nely*3-2; columns = 1:3:(nelx+1)*(nely+1)*3-2;
mstress = matrix_stress(rows,columns); % nelx*nely by (nelx+1)*(nely+1) 
sel_stress = mstress'.*repmat(mises',[(nelx+1)*(nely+1) 1]);
sx = mstress'.*repmat(stress_el(:,1)',[(nelx+1)*(nely+1) 1]);
sy = mstress'.*repmat(stress_el(:,2)',[(nelx+1)*(nely+1) 1]);
txy = mstress'.*repmat(stress_el(:,3)',[(nelx+1)*(nely+1) 1]);
dgdx = zeros(3*(nelx+1)*(nely+1),1);
stg = zeros(size(sel_stress));
for iii = 1:nelx*nely
    stg(:,iii) = sel_stress(:,iii).^(p-1)*(von_mises(iii)*sl*pho(iii))^(-1);
end
for i = 1:(nelx+1)*(nely+1)
dgdx(3*i-2) = 1/4*1/2*sum_stress*stg(i,:)*(2*sx(i,:)'-sy(i,:)'); 
dgdx(3*i-1) = 1/4*1/2*sum_stress*stg(i,:)*(2*sy(i,:)'-sx(i,:)');
dgdx(3*i) = 1/4*1/2*sum_stress*stg(i,:)*6*txy(i,:)';
end
rhs = zeros((nelx+1)*(nely+1)*7,1);
rhs(1:(nelx+1)*(nely+1)*3,1) = dgdx;
adjoint = ma57_solve(K, full(rhs), struct);
adjoint(1:(nelx+1)*(nely+1)*3,1) = B*adjoint(1:(nelx+1)*(nely+1)*3,1);
for ii = 1:nelx*nely
    second(ii) = sum_stress*von_mises(ii)^(p-1)*(-q)*vm(ii)*x(ii)^(3-q-1)*(1/sl);
    third(ii) = sum_stress*mises(ii)^(p-1)*epsi/(pho(ii)^(-2));
end
sum_stress = (sum(mises.^p))^(1/p);
end
% STRESS SENSITIVITIES 
function [dgdp] = SENSITIVITIES(adjoint,U,nelx,nely,x,penal,Msize,dkdp)
dgdp = zeros(nelx*nely,1);
n = nelx*nely;
x = reshape(x',nelx*nely,1);
% U(abs(U)<1e-3) = 0; U = sparse(U);
    for i=1:n
        I = 0; J = 0; S = 0; dk = 0;
        l = (i-1)*Msize+1; u = i*Msize;
        I = dkdp(l:u,1); J = dkdp(l:u,2); S = -penal*x(i)^(-1-penal)*dkdp(l:u,3);
        dk = sparse(I,J,S,7*(nelx+1)*(nely+1),7*(nelx+1)*(nely+1));
        dgdp(i,1) = -adjoint'*dk*U;
    end
end
% MATRICES dkdp 
function [dkdp,Msize] = dkdp(nelx,nely)
[iF,jF,sF,~,~,~,~,~,~,~] = EL(1,1);
ntrip = 0; i = 0; Msize = size(iF,1);
I = zeros((size(iF,1)*nelx*nely),1);
J = zeros((size(iF,1)*nelx*nely),1);
S = zeros((size(iF,1)*nelx*nely),1); 
for ely = 1:nely
    for elx = 1:nelx
        i = i + 1; 
        n1 = elx + (nelx+1)*(ely-1); n2 = elx + (nelx+1)*ely;
        fdof = [n1*3-2 n1*3-1 n1*3 n1*3+1 n1*3+2 n1*3+3 ...
            n2*3-2 n2*3-1 n2*3 n2*3+1 n2*3+2 n2*3+3];        
        for fi = 1:size(iF,1)
            ntrip = ntrip + 1;
            I(ntrip) = fdof(iF(fi)); J(ntrip) = fdof(jF(fi)); 
            S(ntrip) = sF(fi);
        end
    end
end
dkdp = [I J S];
end
% ACTIVE SET OF CONSTRAINTS (q DEFINED)
function [active,cons_val] = ACTIVE_CONSTRAINT(x,von_mises,loop)
n = size(von_mises,1);
i = 0; sigmaL = 6.0; q = 3.0; k = 0;
von_mises(abs(von_mises)<1e-2)= 0.0;
lhs = von_mises./(sigmaL*x.^q);
% thre = linspace(0.75,0.95,15);
% if loop < 16
%     threshold = thre(loop);
% else
     threshold = 0.95;
% end
for j = 1:n
    i = i+1;
    if lhs(j) >= threshold
        k = k + 1;
        active(k,1) = i;
    end
end
cons_val = lhs(active)-1;
end
% ADJOINT FOR LOCAL STRESS CONSTRAINT 
function [lambda] = ADJOINT_LOCAL(struct,K,B,von_mises,active,stress_el,matrix_stress)
n = size(active,1); node = size(K,1)/7;
von_mises = von_mises.^(-1/2);
T = [2 -1 0; -1 2 0; 0 0 6];
for i = 1:n
    rhs = zeros(node*7,1);
    j = active(i);
    dsduk = 1/8*von_mises(j)*stress_el(j,:)*T;
    rhs(1:node*3,1) = (matrix_stress(3*j-2,:)*dsduk(1)+matrix_stress(3*j-1,:)*dsduk(2)+...
        matrix_stress(3*j,:)*dsduk(3))';
    lambda(:,i) = ma57_solve(K, full(rhs), struct);
    lambda(1:3*node,i) = B*lambda(1:3*node,i); 
end
end
% SENSITIVITIES OF LOCAL STRESS CONSTRAINT (q DEFINED)
function [deqdx] = LOCAL_SENS(x,U,lambdas,penal,dkdp,Msize,nelx,nely,active,von_mises)
n = nelx*nely; nact = size(lambdas,2);
node = size(U,1); q = 3.0; E = eye(n);
U(abs(U)<1e-3) = 0; U = sparse(U);
    for i=1:n
        l = (i-1)*Msize+1; u = i*Msize;
        I = dkdp(l:u,1); J = dkdp(l:u,2); S = -penal*x(i)^(-1-penal)*dkdp(l:u,3);
        dk = sparse(I,J,S,7*(nelx+1)*(nely+1),7*(nelx+1)*(nely+1));
        %for k=1:nact
            dtaldx(:,i) = lambdas(:,:)'*dk*U;
        %end
    end
    for k=1:nact
        deqdx(k,:) = x(active(k))^(-q)*dtaldx(k,:)-...
            q*E(active(k),:)*x(active(k))^(-q-1)*von_mises(active(k));
    end
    deqdx = deqdx/6;
end
