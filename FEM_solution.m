
%*********************************************************************
% Computation of the homogenized solution
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Sol,x1] = FEM_solution(A,T,ele1,ele2,...
    A_EfectivePerm,field,p_prev,~)
%%
%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************
global Macro_geo N_real b_nolin
global L_scheme bp_nolin Lref 

% Assemble matrix B
B    = sparse(Macro_geo.(field).noedges, Macro_geo.(field).noedges);

Sol      = struct();
Sol.Pres = zeros(Macro_geo.(field).nElement,1);
Sol.Vel  = zeros(Macro_geo.(field).noedges,1);

%% REFERENCE MESH FOR PERMEABILITY
[meshx,meshy] = meshgrid(linspace(0,N_real(1),N_real(1)+1),...
    linspace(0,N_real(2),N_real(2)+1));
meshx = meshx/Lref; meshy = meshy/Lref;

%% Matrix B
for j = 1:Macro_geo.(field).nElement
    coord = Macro_geo.(field).coordinate(Macro_geo.(field).element(j,:),:)';
    I = diag(Macro_geo.(field).nodes2edge(Macro_geo.(field).element(j,[2 3 1]),Macro_geo.(field).element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==Macro_geo.(field).edge2element(I,4)))= -1;
    
    aux11 = interp2(meshx,meshy,A_EfectivePerm.(field)(:,:,1)',...
        Macro_geo.(field).bari(j,1),Macro_geo.(field).bari(j,2),'nearest');
    aux22 = interp2(meshx,meshy,A_EfectivePerm.(field)(:,:,2)',...
        Macro_geo.(field).bari(j,1),Macro_geo.(field).bari(j,2),'nearest');
    aux12 = interp2(meshx,meshy,A_EfectivePerm.(field)(:,:,3)',...
        Macro_geo.(field).bari(j,1),Macro_geo.(field).bari(j,2),'nearest');
    aux21 = interp2(meshx,meshy,A_EfectivePerm.(field)(:,:,4)',...
        Macro_geo.(field).bari(j,1),Macro_geo.(field).bari(j,2),'nearest');
    
    %     perm_eval2 = [1/aux11 0; 0 1/aux22];
    perm_eval = inv([aux11 aux12; aux21 aux22]);
    
    B(I,I)= B(I,I)+ diag(signum)*...
        stimaB(coord,perm_eval)*diag(signum);
end
% Complete solution
A(1:Macro_geo.(field).noedges,1:Macro_geo.(field).noedges) = B;

%% BOUNDARY CONDITIONS
k     = boundary(Macro_geo.(field).coordinate(:,1),Macro_geo.(field).coordinate(:,2));
Gamma = [k(1:end-1) k(2:end)];

tmp =zeros(Macro_geo.(field).noedges+Macro_geo.(field).nElement,1);
tmp(diag(Macro_geo.(field).nodes2edge(Gamma(:,1),Gamma(:,2))))=...
    ones(size(diag(Macro_geo.(field).nodes2edge(Gamma(:,1),Gamma(:,2))),1),1);

FreeEdge = find(~tmp); % + nElemnt numeration
freePOS  = (setdiff(FreeEdge,Macro_geo.(field).noedges+[ele1;ele2]))';

it = 1;
residual = zeros(0);
residual(1) = inf;
bprev = b_nolin(p_prev);
pprev_it = p_prev;
bprev_it = b_nolin(p_prev);
change_nonLinear = 0;

ind_max = zeros(size(pprev_it,1),1);
while residual(it)>1E-8
    if residual(it) < 1E-1 && change_nonLinear==0
        change_nonLinear = 1;
    end
    if change_nonLinear == 1
        LL = bp_nolin(pprev_it);
    else
        LL = L_scheme;
    end
    
    b_it = [zeros(Macro_geo.(field).noedges,1); ...
        T*(bprev - bprev_it + LL.*pprev_it)];
    
    
    %% SOURCE TERM - Impossing Pressure
    x_it=zeros(Macro_geo.(field).noedges+Macro_geo.(field).nElement,1);
    x_it(Macro_geo.(field).noedges+ele1,1) = 1;
    x_it(Macro_geo.(field).noedges+ele2,1) = 0;
    
    A(end-Macro_geo.(field).nElement+1:end,...
        end-Macro_geo.(field).nElement+1:end) = LL.*T;
    
    b_it = b_it-A*x_it;
    
    % LINEAR SOLUTION (DIRECT SOLVER)
    x_it(freePOS,1) = A(freePOS,freePOS)\b_it(freePOS);
    
    it = it+1;
    
    differ = (x_it(end-Macro_geo.(field).nElement+1:end)- pprev_it);
    residual(it) = norm(differ);
    %         if sign(residual(it) - residual(it-1))==1
    %             disp('up')
    %         end
    
    pprev_it = x_it(end-Macro_geo.(field).nElement+1:end);
    ind_max = max([ind_max,differ./abs(pprev_it)],[],2);
    bprev_it = b_nolin(pprev_it);
    %     disp(field)
    %     disp(it)
    %     disp(residual(it))
end

%% POST PROCESSING
% X1 PRESURE =(1*NELEMENT)
x1 = x_it(end-Macro_geo.(field).nElement+1:end);
% X2 VELOCITY =(2*nedges)
Sol.flux = x_it(1:end-Macro_geo.(field).nElement);
Sol.residual = residual;
Sol.ind_mx  = ind_max;

Sol.PresCont  = zeros(size(Macro_geo.(field).coordinate,1),1);
contador_Cont = zeros(size(Macro_geo.(field).coordinate,1),1);

for j = 1:Macro_geo.(field).nElement
    pos = 3*(j-1)+1:3*j;
    % Sol.Pres= 3*Nelement
    
    Sol.Pres(pos,1)= x1(j);
    % Sol.Pres= NCOORDINATES
    
    Sol.PresCont(Macro_geo.(field).element(j,:),1)= ...
        Sol.PresCont(Macro_geo.(field).element(j,:),1)+x1(j);
    
    contador_Cont(Macro_geo.(field).element(j,:),1) = contador_Cont(Macro_geo.(field).element(j,:),1) +1;
end

Sol.PresCont = Sol.PresCont./contador_Cont;

% Sol.flux = x2;
[Sol.Vel,Sol.VelCont,Sol.diverg]=...
    fluxEB(Macro_geo.(field).element,Macro_geo.(field).coordinate,...
    -Sol.flux,Macro_geo.(field).nodes2edge,Macro_geo.(field).edge2element);

end
