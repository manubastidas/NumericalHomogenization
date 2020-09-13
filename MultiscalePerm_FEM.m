%*********************************************************************
% Calculation of the Effective permeability. 
% This function calls the Micro-scale solvers and the integral calcul
% and assamble the effective parameter.
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [A_Efective,A_Efective_anis] = MultiscalePerm_FEM(K_perm,NCoarse)

%% Create Micro mesh
% Micro features
global N_real Lref Harm

Nmicro = N_real./NCoarse;
[xx_micro,yy_micro] = meshgrid(0:Nmicro(1),0:Nmicro(2));
xx_micro = xx_micro/Lref; yy_micro = yy_micro/Lref;

Micro_geo            = struct();
Micro_geo.coordinate = [xx_micro(:),yy_micro(:)];
Micro_geo.element    = delaunay(xx_micro,yy_micro);

% nElement -> number of element at each mesh
Micro_geo.nElement  = size(Micro_geo.element,1);
% nnodes -> number of nodes at each mesh
Micro_geo.nnodes    = size(Micro_geo.coordinate,1);

[~,Micro_geo.nodes2edge,Micro_geo.noedges,Micro_geo.edge2element,Micro_geo.interioredge] = ...
    edge(Micro_geo.element,Micro_geo.coordinate);

Micro_geo = PreProcess_Identification(Micro_geo);
Micro_geo = Position_indicators(Micro_geo,0);

%%
%*********************************************************************
%*                                                                   *
%*                    EFFECTIVE DIFUSSION TENSOR - COARSE MESH          *
%*                                                                   *
%*********************************************************************

A_Efective      = zeros(NCoarse);
A_Efective_anis = zeros(NCoarse);

numberCoarse = numel(A_Efective);
for e = 1:numberCoarse
    
    fprintf('\n Cell %i/%i',e,numberCoarse);
    
    % Choose diffusion tensor
    %     i  = ceil(e/NCoarse(2));
    %     j  = e-(i-1)*NCoarse(2);
    [j,i] = ind2sub(NCoarse,e);
    
    Nmicro_aux = floor(Nmicro);
    
    K_micro11 = K_perm(((i-1)*Nmicro_aux(1)+1):i*Nmicro_aux(1)+1,...
        ((j-1)*Nmicro_aux(2)+1):j*Nmicro_aux(2)+1,1);
    K_micro22 = K_perm(((i-1)*Nmicro_aux(1)+1):i*Nmicro_aux(1)+1,...
        ((j-1)*Nmicro_aux(2)+1):j*Nmicro_aux(2)+1,2);
    K_micro12 = K_perm(((i-1)*Nmicro_aux(1)+1):i*Nmicro_aux(1)+1,...
        ((j-1)*Nmicro_aux(2)+1):j*Nmicro_aux(2)+1,3);
    K_micro21 = K_perm(((i-1)*Nmicro_aux(1)+1):i*Nmicro_aux(1)+1,...
        ((j-1)*Nmicro_aux(2)+1):j*Nmicro_aux(2)+1,4);
    
    if Harm == 1
        %     m = harmmean(X,'all')
        A_Efective(j,i,1) = harmmean(K_micro11,'all');
        A_Efective(j,i,2) = harmmean(K_micro22,'all');
        A_Efective_anis(j,i,1) = harmmean(K_micro12,'all');
        A_Efective_anis(j,i,2) = harmmean(K_micro21,'all');
    else
        % NUMERICAL GRADIENT :(
        [Kx11,~] = gradient(K_micro11);
        [~,Ky21] = gradient(K_micro21);
        [Kx12,~] = gradient(K_micro12);
        [~,Ky22] = gradient(K_micro22);
        
        fx = Kx11+Ky21;
        fy = Kx12+Ky22;
        [Vel1, Vel2] = MicroSolver_FEM(Micro_geo,xx_micro,yy_micro,K_micro11,K_micro12,K_micro21,K_micro22,fx,fy);
        
        [K_Efective] = EfectivePermTensor(Micro_geo,xx_micro,yy_micro,K_micro11,K_micro12,K_micro21,K_micro22,Vel1,Vel2);
        
        A_Efective(j,i,1) = abs(K_Efective(1,1));
        A_Efective(j,i,2) = abs(K_Efective(2,2));
        A_Efective_anis(j,i,1) = abs(K_Efective(1,2));
        A_Efective_anis(j,i,2) = abs(K_Efective(2,1));
    end
end
A_Efective = [A_Efective(:,:,:) A_Efective(:,end,:)];
A_Efective = [A_Efective(:,:,:);A_Efective(end,:,:)];

A_Efective_anis = [A_Efective_anis(:,:,:) A_Efective_anis(:,end,:)];
A_Efective_anis = [A_Efective_anis(:,:,:);A_Efective_anis(end,:,:)];