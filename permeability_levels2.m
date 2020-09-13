%*********************************************************************
% Auxiliar code to copy the values of the coarse scale permeabilities 
% to a fine scale resolution
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium


function [Klevels2] = permeability_levels2(Klevels,numl)

for rr=0:numl-1
    field = sprintf('ref%i',rr);
    
    A1 = Klevels.(field).tensor(1:end-1,1:end-1,1);
    A2 = Klevels.(field).tensor(1:end-1,1:end-1,2);
    A3 = Klevels.(field).anisot(1:end-1,1:end-1,1);
    A4 = Klevels.(field).anisot(1:end-1,1:end-1,2);
    
    A_EfectivePerm1 = kron(A1,ones(2^(numl-rr)));
    A_EfectivePerm2 = kron(A2,ones(2^(numl-rr)));
    A_EfectivePerm3 = kron(A3,ones(2^(numl-rr)));
    A_EfectivePerm4 = kron(A4,ones(2^(numl-rr)));

    A_EfectivePerm1 = [A_EfectivePerm1(:,:) A_EfectivePerm1(:,end)];
    A_EfectivePerm1 = [A_EfectivePerm1(:,:);A_EfectivePerm1(end,:)];

    A_EfectivePerm2 = [A_EfectivePerm2(:,:) A_EfectivePerm2(:,end)];
    A_EfectivePerm2 = [A_EfectivePerm2(:,:);A_EfectivePerm2(end,:)];
    
    A_EfectivePerm3 = [A_EfectivePerm3(:,:) A_EfectivePerm3(:,end)];
    A_EfectivePerm3 = [A_EfectivePerm3(:,:);A_EfectivePerm3(end,:)];
   
    A_EfectivePerm4 = [A_EfectivePerm4(:,:) A_EfectivePerm4(:,end)];
    A_EfectivePerm4 = [A_EfectivePerm4(:,:);A_EfectivePerm4(end,:)];
    
    Klevels2.(field)(:,:,1) = A_EfectivePerm1;
    Klevels2.(field)(:,:,2) = A_EfectivePerm2;
    Klevels2.(field)(:,:,3) = A_EfectivePerm3;
    Klevels2.(field)(:,:,4) = A_EfectivePerm4;
% malla refinada (artificial)
% Klevels2.(field) = cat(3,A_EfectivePerm1(1:end-(2*(numl-rr)-1),1:end-(2*(numl-rr)-1),1),...
%     A_EfectivePerm2(1:end-(2*(numl-rr)-1),1:end-(2*(numl-rr)-1),1));

end

if isempty(rr)
    field= 'ref0';
else
field = sprintf('ref%i',rr+1);
end
Klevels2.(field)(:,:,1) = Klevels.(field).tensor(:,:,1);
Klevels2.(field)(:,:,2) = Klevels.(field).tensor(:,:,2);
Klevels2.(field)(:,:,3) = Klevels.(field).anisot(:,:,1);
Klevels2.(field)(:,:,4) = Klevels.(field).anisot(:,:,2);