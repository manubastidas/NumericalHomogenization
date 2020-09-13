
%*********************************************************************
% Computation of the reference solution
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Macro_SolREF] = ReferenceSolution_FEM(Macro_geoREF,Macro_SolREF,A,ele1,ele2,T)

global b_nolin L_scheme Time bp_nolin

%%
%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************
%
field     = sprintf('time%i',0);

k     = boundary(Macro_geoREF.coordinate(:,1),Macro_geoREF.coordinate(:,2));
Gamma = [k(1:end-1) k(2:end)];

tmp =zeros(Macro_geoREF.noedges+Macro_geoREF.nElement,1);
tmp(diag(Macro_geoREF.nodes2edge(Gamma(:,1),Gamma(:,2))))=...
    ones(size(diag(Macro_geoREF.nodes2edge(Gamma(:,1),Gamma(:,2))),1),1);
% temp = 1 at Neumann edges

FreeEdge=find(~tmp); % + nElemnt numeration
freePOS =(setdiff(FreeEdge,Macro_geoREF.noedges+[ele1;ele2]))';

p_prev = Macro_SolREF.(field).Pres;
Macro_SolREF.norm = zeros(Time.tnSteps,1);
for tt = 1:Time.tnSteps
    fprintf('\n Time %i/%i',tt,Time.tnSteps);
    
    %     b = sparse(Macro_geoREF.noedges+Macro_geoREF.nElement,1);
    %     temp = Time.time_vec(tt);
    it = 1;
    residual = zeros(0);
    residual(1) = inf;
    bprev = b_nolin(p_prev);
    pprev_it = p_prev;
    bprev_it = b_nolin(p_prev);
    change_nonLinear = 0;
    
    while residual(it)>1E-8
        if residual(it) < 1E-1 && change_nonLinear==0
            change_nonLinear = 1;
        end
        if change_nonLinear == 1
            LL = bp_nolin(pprev_it);
        else
            LL = L_scheme;
        end
        
        b_it = [zeros(Macro_geoREF.noedges,1); ...
            T*(bprev - bprev_it + LL.*pprev_it)];
        
        x_it = zeros(Macro_geoREF.noedges+Macro_geoREF.nElement,1);
        x_it(Macro_geoREF.noedges+ele1,1) = 1;
        x_it(Macro_geoREF.noedges+ele2,1) = 0;
        
        A(end-Macro_geoREF.nElement+1:end,...
            end-Macro_geoREF.nElement+1:end) = LL.*T;
        
        b_it = b_it-A*x_it;
        x_it(freePOS) = A(freePOS,freePOS)\b_it(freePOS);
        
        it = it+1;
        residual(it) = norm(x_it(end-Macro_geoREF.nElement+1:end)- pprev_it);
        pprev_it= x_it(end-Macro_geoREF.nElement+1:end);
        bprev_it = b_nolin(pprev_it);
        
    end
    
    x1 = x_it(end-Macro_geoREF.nElement+1:end);
    x2 = x_it(1:end-Macro_geoREF.nElement);
    
    field     = sprintf('time%i',tt);
    
    contador = zeros(size(Macro_geoREF.coordinate,1),1);
    
    Macro_SolREF.(field).PresCont = zeros(size(Macro_geoREF.coordinate,1),1);
    for j = 1:Macro_geoREF.nElement
        
        pos = 3*(j-1)+1:3*j;
        Macro_SolREF.(field).Pres(pos,1)= x1(j);
        Macro_SolREF.(field).PresCont(Macro_geoREF.element(j,:),1) =...
            Macro_SolREF.(field).PresCont(Macro_geoREF.element(j,:),1) + Macro_SolREF.(field).Pres(pos,1);
        contador(Macro_geoREF.element(j,:),1)= contador(Macro_geoREF.element(j,:),1) + ones(3,1);
    end
    Macro_SolREF.(field).PresCont = Macro_SolREF.(field).PresCont./contador;
    
    Macro_SolREF.(field).flux = x2;
    [Macro_SolREF.(field).Vel,Macro_SolREF.(field).VelCont,Macro_SolREF.(field).diverg]=fluxEB(Macro_geoREF.element,Macro_geoREF.coordinate,...
        -x2,Macro_geoREF.nodes2edge,Macro_geoREF.edge2element);
    
    p_prev = x1;
    
    [~, Macro_SolREF.norm(tt), ~] =...
        errorL2_Macro(0.*Macro_SolREF.(field).Pres,...
        Macro_SolREF.(field).Pres,Macro_geoREF,Macro_geoREF);
    Macro_SolREF.(field).residual = residual;
end
