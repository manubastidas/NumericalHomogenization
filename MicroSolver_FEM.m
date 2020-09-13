%*********************************************************************
% Micro scale solver (MFEM) 
%
% This code is based on:
% Bahriawati, C., & Carstensen, C. (2005). 
% Three MATLAB implementations of the lowest-order Raviart-Thomas 
% MFEM with a posteriori error control. 
% Computational Methods in Applied Mathematics, 5(4), 333-361.
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Micro_sol1, Micro_sol2] = MicroSolver_FEM(Micro_geo,xx,yy,K_micro11,K_micro12,K_micro21,K_micro22,Kx,Ky)

element      = Micro_geo.element;
coordinate   = Micro_geo.coordinate;
noedges      = Micro_geo.noedges;
nodes2edge   = Micro_geo.nodes2edge;
edge2element = Micro_geo.edge2element;

% Assemble matrices B and C
B = sparse(noedges, noedges);
C = sparse(noedges,Micro_geo.nElement);

% Volume force
b1 = sparse(noedges+size(element,1),1);
b2 = sparse(noedges+size(element,1),1);

for j = 1:Micro_geo.nElement
    coord = coordinate(element(j,:),:)';
    I = diag(nodes2edge(element(j,[2 3 1]),element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==edge2element(I,4)))= -1;
    
    bari = sum(coord,2)./3;
    
    aux11 = interp2(xx,yy,K_micro11',bari(1),bari(2),'nearest');
    aux12 = interp2(xx,yy,K_micro12',bari(1),bari(2),'nearest');
    aux21 = interp2(xx,yy,K_micro21',bari(1),bari(2),'nearest');
    aux22 = interp2(xx,yy,K_micro22',bari(1),bari(2),'nearest');
    
    K_aux = inv([aux11 aux12; aux21 aux22]);
    B(I,I)= B(I,I)+ diag(signum)*...
        stimaB(coord,K_aux)*diag(signum);
    
    n = coord(:,[3,1,2])-coord(:,[2,3,1]);
    C(I,j) = diag(signum)*[norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]';
    
    b1(noedges+j)= det([1,1,1; coordinate(element(j,:),:)']) * ...
        interp2(xx,yy,Kx',bari(1),bari(2))/2;
    
    b2(noedges+j)= det([1,1,1; coordinate(element(j,:),:)']) * ...
        interp2(xx,yy,Ky',bari(1),bari(2))/2;
end

% Global stiffness matrix
M = [B ,         C,      ;
    C', sparse(size(C,2),size(C,2))];

[x1,x2] = Micro_Boundary(M,b1,b2,Micro_geo);


Micro_sol1 = x1(1:end-Micro_geo.nElement);
Micro_sol2 = x2(1:end-Micro_geo.nElement);
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

