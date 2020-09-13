
%*********************************************************************
% Offline computation of the different permeability levels. 
%*********************************************************************
%
%***------------------------------------
%***Inputs: Heterogeneos permeability field and partition
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [A_Efective,numl] = permeability_levels(K_perm,NCoarse)
% RESTRICTIONS:  Nmicro = [speX speY]./NCoarse -> Even Integer
%                Nmicro(1) ~ Nmicro(2) -> Square domains

%%
global N_real Lref 
Nmicro  = N_real./NCoarse;

A_Efective = struct();
% A_Efective_anis = struct();
rr=0;

while min(Nmicro)-1>=1
    field = sprintf('ref%i',rr);
    
    [x_coarse, y_coarse]        = meshgrid(0:Nmicro(1):N_real(1),0:Nmicro(2):N_real(2));
    x_coarse = x_coarse/Lref; y_coarse= y_coarse/Lref;
    
    A_Efective.(field).gridmesh = {x_coarse,y_coarse};
    
    [A_Efective.(field).tensor,A_Efective.(field).anisot] = MultiscalePerm_FEM(K_perm,NCoarse);
    
    NCoarse = NCoarse*2;
    Nmicro = N_real./NCoarse;
    rr = rr+1;
end

%% Ultimo nivel Si la division es par
if isequal(Nmicro,[1 1])
    numl = rr;
    % Last one is the original field
    field = sprintf('ref%i',numl);
    [x_coarse, y_coarse]        = meshgrid(0:N_real(1),0:N_real(2));
     x_coarse = x_coarse/Lref; y_coarse= y_coarse/Lref;
     
    A_Efective.(field).gridmesh = {x_coarse,y_coarse};
    [A_Efective.(field).tensor(:,:,1)] =  K_perm(:,:,1)';
    [A_Efective.(field).tensor(:,:,2)] =  K_perm(:,:,2)';
    [A_Efective.(field).anisot(:,:,1)] =  K_perm(:,:,3)';
    [A_Efective.(field).anisot(:,:,2)] =  K_perm(:,:,4)';
else 
    numl = rr-1;
end



