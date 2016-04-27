
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

matrix2latexmatrix(A, 'a.tex');
matrix2latexmatrix(B, 'b.tex');
matrix2latexmatrix(C, 'c.tex');
