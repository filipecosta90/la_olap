
A_s=sparse(3,2); 

A_s(1,1)=5; 
A_s(1,2)=2; 
A_s(2,1)=1; 
A_s(2,2)=-1; 
A_s(3,2)=3; 

A = full (A_s);

B_s=sparse(2,2); 

B_s(1,1)=7; 
B_s(1,2)=2; 
B_s(2,1)=3; 
B_s(2,2)=-1; 

B = full (B_s);
C = kron ( A, B);



D = sparse(6,6);
E = sparse(6,2);

D_E = sparse(6,2);

F = sparse(6,6);
G = sparse(6,2);

F_G = sparse(6,2);

D(1,1)=5; D(2,2)=5; D(3,3)=1; D(4,4)=1;D(5,5)=0;D(6,6)=0;
E(1,1)=7; E(2,1)=3; E(3,1)=7; E(4,1)=3; E(5,1)=7; E(6,1)=3;
%E(1,2)=2; E(2,2)=-1; E(3,2)=2; E(4,2)=-1; E(5,2)=2; E(6,2)=-1;

F(1,1)=2; F(2,2)=2; F(3,3)=-1; F(4,4)=-1;F(5,5)=3;F(6,6)=3;
%G(1,3)=7; G(2,3)=3; G(3,3)=7; G(4,3)=3; G(5,3)=7; G(6,3)=3;
G(1,2)=2; G(2,2)=-1; G(3,2)=2; G(4,2)=-1; G(5,2)=2; G(6,2)=-1;

D_E = D*E;
F_G = F*G;
D_E_plus_F_G = full( D_E + F_G );
 
matrix2latexmatrix(A, 'krao.tex');
matrix2latexmatrix(B, 'krao.tex');
matrix2latexmatrix(C, 'krao.tex');


matrix2latexmatrix(full(D), 'krao.tex');
matrix2latexmatrix(full(E), 'krao.tex');
matrix2latexmatrix(full(D_E), 'krao.tex');
matrix2latexmatrix(full(F), 'krao.tex');
matrix2latexmatrix(full(G), 'krao.tex');
matrix2latexmatrix(full(F_G), 'krao.tex');

matrix2latexmatrix( D_E_plus_F_G , 'krao.tex');






