
%*********************************************************************
% Preprocess matrices for computing the homogenized solution
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium


function [A,elegidos,T] = FEM_solutionPreprocess(field)
%%
%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************
global Macro_geo N_real Time nondimC Lref

[~,Macro_geo.(field).nodes2edge,Macro_geo.(field).noedges,...
    Macro_geo.(field).edge2element,~] = ...
    edge(Macro_geo.(field).element,Macro_geo.(field).coordinate);

% Assemble matrices B and C
% B = sparse(Macro_geo.(field).noedges, Macro_geo.(field).noedges);
C = sparse(Macro_geo.(field).noedges,Macro_geo.(field).nElement);
T = sparse(Macro_geo.(field).nElement,Macro_geo.(field).nElement);

Macro_geo.(field).size = [0 inf];

%% SOURCE TERM
% Auxiliar code
int = 8;
P = [ 0 0;int 0;0 int
    N_real(1) N_real(2); N_real(1)-int N_real(2); N_real(1) N_real(2)-int];
P = P/Lref;
TR = triangulation([1 2  3; 4 5 6],P);
ele1 = []; ele2=[];

%% MASS MATRIX ENSAMBLE

for j = 1:Macro_geo.(field).nElement
    
    coord = Macro_geo.(field).coordinate(Macro_geo.(field).element(j,:),:)';
    I = diag(Macro_geo.(field).nodes2edge(Macro_geo.(field).element(j,[2 3 1]),Macro_geo.(field).element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==Macro_geo.(field).edge2element(I,4)))= -1;
    
    bari = Macro_geo.(field).bari(j,:);
%      sum(coord,2)/3; = bari;
    
    %     perm_eval = K_permLevel(A_EfectivePerm,bari,Nreal);
    %
    %     B(I,I)= B(I,I)+ diag(signum)*...
    %         stimaB(coord,perm_eval)*diag(signum);
    
    n      = coord(:,[3,1,2])-coord(:,[2,3,1]);
    C(I,j) = diag(signum)*[norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]';
    
    T(j,j) = det([1,1,1; coord])/2;
    
    %% Source 
    stbary= cartesianToBarycentric(TR,1,bari);
    if stbary(1)>0 && stbary(2)>0
        ele1 = [ele1;j];
    end
    stbary= cartesianToBarycentric(TR,2,bari);
    if stbary(1)>0 && stbary(2)>0
        ele2 = [ele2;j];
    end
    Macro_geo.(field).area(j) = T(j,j);
    Macro_geo.(field).size(1) = max(Macro_geo.(field).size(1),max([norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]));
    Macro_geo.(field).size(2) = min(Macro_geo.(field).size(2),min([norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]));
end
% Global stiffness matrix A
A        = [zeros(Macro_geo.(field).noedges) , C ;
            -(Time.dt/nondimC)*C', T];
elegidos = {ele1,ele2};

end
