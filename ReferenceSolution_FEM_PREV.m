%*********************************************************************
% Preprocess matrices for computing the reference solution
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Macro_geoREF,Macro_SolREF,A,ele1,ele2,T] = ReferenceSolution_FEM_PREV(K_perm)
%%
%*********************************************************************
%*                                                                   *
%*                 TRIANGULAR MESH FOR THE MACRO SOLUTION           *
%*                                                                   *
%*********************************************************************
global N_real Time nondimC Lref ISO

[xx,yy] = meshgrid(0:N_real(1),0:N_real(2));
xx = xx/Lref; yy=yy/Lref;

Macro_geoREF         = struct();
Macro_geoREF.coordinate = [xx(:),yy(:)];
Macro_geoREF.element    = delaunay(xx,yy);

% nElement -> number of element at each mesh
Macro_geoREF.nElement  = size(Macro_geoREF.element,1);
% nnodes -> number of nodes at each mesh
Macro_geoREF.nnodes = size(Macro_geoREF.coordinate,1);

% MACRO features
[~,Macro_geoREF.nodes2edge,Macro_geoREF.noedges,Macro_geoREF.edge2element,Macro_geoREF.interioredge] = ...
    edge(Macro_geoREF.element,Macro_geoREF.coordinate);

%%
%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************

% Assemble matrices B and C
B = sparse(Macro_geoREF.noedges,Macro_geoREF.noedges);
C = sparse(Macro_geoREF.noedges,Macro_geoREF.nElement);
T = sparse(Macro_geoREF.nElement,Macro_geoREF.nElement);

Macro_geoREF.size = 0;
Macro_SolREF      = struct();

field     = sprintf('time%i',0);
Macro_SolREF.(field).Pres = zeros(Macro_geoREF.nElement,1);
Macro_SolREF.(field).Vel  = zeros(Macro_geoREF.noedges,1);
Macro_SolREF.(field).PresCont  = zeros(Macro_geoREF.nnodes,1);

% To impose pressure
infl = 8;
P  = [ 0 0;infl 0;0 infl
    N_real(1) N_real(2); N_real(1)-infl N_real(2); N_real(1) N_real(2)-infl];
P = P/Lref;
sT = [1 2  3; 4 5 6];
TR = triangulation(sT,P);
ele1 = []; ele2 = [];

for j = 1:Macro_geoREF.nElement
    
    coord = Macro_geoREF.coordinate(Macro_geoREF.element(j,:),:)';
    I = diag(Macro_geoREF.nodes2edge(Macro_geoREF.element(j,[2 3 1]),Macro_geoREF.element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==Macro_geoREF.edge2element(I,4)))= -1;
    
    bari = sum(coord,2)/3;
    area(j) = det([1,1,1; coord])/2;
    
    n = coord(:,[3,1,2])-coord(:,[2,3,1]);
    C(I,j) = diag(signum)*[norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]';
    
    T(j,j) = area(j);
    
    stbary= cartesianToBarycentric(TR,1,bari');
    if stbary(1)>0 && stbary(2)>0
        ele1 = [ele1;j];
    end
    stbary= cartesianToBarycentric(TR,2,bari');
    if stbary(1)>0 && stbary(2)>0
        ele2 = [ele2;j];
    end
    
    Macro_geoREF.size = max(Macro_geoREF.size,max([norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]));
    
%     if ISO ==1
%         aux1 = interp2(xx,yy,K_perm,bari(1),bari(2),'nearest');
%         %         aux2 = interp2(CoarseX,CoarseY,A_Efective(:,:,2)',bari(1),bari(2));
%         B(I,I)= B(I,I)+ diag(signum)*...
%             stimaB(coord,[1/aux1 0;0 1/aux1])*diag(signum);
%     else
        aux11 = interp2(xx,yy,K_perm(:,:,1),bari(1),bari(2),'nearest');
        aux22 = interp2(xx,yy,K_perm(:,:,2),bari(1),bari(2),'nearest');
        aux12 = interp2(xx,yy,K_perm(:,:,3),bari(1),bari(2),'nearest');
        aux21 = interp2(xx,yy,K_perm(:,:,4),bari(1),bari(2),'nearest');
        b_aux = inv([aux11 aux12; aux21 aux22]);
        B(I,I)= B(I,I)+ diag(signum)*...
            stimaB(coord,b_aux)*diag(signum);
%     end
end
A = [B ,         C,      ;
    -(Time.dt/nondimC)*C', T];

end

function B=stimaB(coord,A)

%N=coord(:)*ones(1,3)-[coord;coord;coord];

N=coord(:)*ones(1,3)-repmat(coord,3,1);

D=diag([norm(N([5,6],2)) norm(N([1,2],3)) norm(N([1,2],2))]);

M=spdiags([ones(6,1),ones(6,1),2*ones(6,1),ones(6,1),ones(6,1)],...
          [-4,-2,0,2,4],6,6);

% N_aux = (repmat(diag(A),3,3).*N)'; % ONLY DIAGONAL

coord_aux = A*coord;
N_aux=(coord_aux(:)*ones(1,3)-repmat(coord_aux,3,1))';

B = D*N_aux*M*N*D/(24*det([1,1,1;coord])); 
end