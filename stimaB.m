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