////////////////////////////////////////////////// 
// This code is part of the pcc_q MAGMA library //
//                                              //
// copyright (c) 2017 Jan Tuitman               //
//////////////////////////////////////////////////


ord_0:=function(f)

  // Compute ord_0(f), where f is a rational function.

  return Valuation(Numerator(f))-Valuation(Denominator(f));
end function;


ord_0_mat:=function(A)

  // Compute ord_0(A), where A is a matrix of rational functions.

  v:=ord_0(A[1,1]);
  for i:=1 to NumberOfRows(A) do
    for j:=1 to NumberOfColumns(A) do
      if ord_0(A[i,j]) lt v then
        v:=ord_0(A[i,j]);
      end if;
    end for;
  end for;
  return v;
end function;


ord_inf:=function(f)

  // Compute ord_inf(f), where f is a rational function.

  return -Degree(Numerator(f))+Degree(Denominator(f));
end function;


ord_inf_mat:=function(A);

  // Compute ord_inf(A), where A is a matrix of rational functions.

  v:=ord_inf(A[1,1]);
  for i:=1 to NumberOfRows(A) do
    for j:=1 to NumberOfColumns(A) do
      if ord_inf(A[i,j]) lt v then
        v:=ord_inf(A[i,j]);
      end if;
    end for;
  end for;
  return v;
end function;


ord_r:=function(f,r)

  // Compute ord_r(f), where f is a rational function.
 
  if f eq 0 then
    return Infinity();
  end if;
  fac:=Factorization(r);
  vlist:=[];
  for i:=1 to #fac do
    v:=0;
    while (Numerator(f) mod fac[i][1] eq 0) do
      f:=f/fac[i][1];
      v:=v+1;
    end while;
    while (Denominator(f) mod fac[i][1] eq 0) do
      f:=f*fac[i][1];
      v:=v-1;
    end while;
    Append(~vlist,v);
  end for; 
  min:=Minimum(vlist);
  return min;
end function;


ord_r_mat:=function(A,r) 

  // Compute ord_r(A), where f is a matrix of rational functions.

  v:=ord_r(A[1,1],r);
  for i:=1 to NumberOfRows(A) do
    for j:=1 to NumberOfColumns(A) do
      if ord_r(A[i,j],r) lt v then
        v:=ord_r(A[i,j],r);
      end if;
    end for;
  end for;
  return v;
end function;


Zax_to_Kx:=function(f,Kx);

  // Convert from Zax to Kx.
 
  Zax:=Parent(f); Za:=BaseRing(Zax);
  K:=BaseRing(Kx);
  phi:=hom<Za->K|K.1>; phi:=hom<Zax->Kx|phi,Kx.1>;
  fK:=phi(f);
  
  return fK;
end function;


Zaxy_to_Kxy:=function(f,Kxy);

  // Convert from Zaxy to Kxy.

  C:=Coefficients(f);
  D:=[];
  for i:=1 to #C do
    D[i]:=Zax_to_Kx(C[i],BaseRing(Kxy));
  end for;
  return Kxy!D;
end function;


push_to_Kt:=function(f,K);

  // Push element to rational function field over K
  
  if f eq 0 then
    return f;
  end if;
  
  KT:=PolynomialRing(K);
  
  C:=Coefficients(Numerator(f));
  C:=[Eltseq(coef):coef in C];
  C:=[K!elt:elt in C];
  
  D:=Coefficients(Denominator(f)); 
  D:=[Eltseq(coef) : coef in D];
  D:=[K!elt:elt in D];
  
  return Polynomial(C)/Polynomial(D);
end function;


push_to_Kt_mat:=function(A,K);

  // Push matrix to the rational function field over K
  
  Kt:=RationalFunctionField(K);
  B:=ZeroMatrix(Kt,NumberOfRows(A),NumberOfColumns(A));
  for i:=1 to NumberOfRows(A) do
    for j:=1 to NumberOfColumns(A) do
      B[i,j]:=push_to_Kt(A[i,j],K);
    end for;
  end for;
  return B;
end function;


mat_W0:=function(Q,Kxy)

  // Compute the matrix W0 using MaximalOrderFinite

  Kx<x>:=RationalFunctionField(CoefficientRing(CoefficientRing(Kxy)));
  fun_field:=ext<Kx|Zaxy_to_Kxy(Q,Kxy)>;
  b0:=Basis(MaximalOrderFinite(fun_field));
  d:=Degree(Q);
  mat:=ZeroMatrix(Kx,d,d);
  for i:=1 to d do
    for j:=1 to d do
      mat[i,j]:=Eltseq(fun_field!b0[i])[j];
    end for;
  end for;
  return mat;

end function;


mat_Winf:=function(Q,Kxy);

  // Compute the matrix Winf using MaximalOrderFinite

  Kx<x>:=RationalFunctionField(CoefficientRing(CoefficientRing(Kxy)));
  Kxy<y>:=PolynomialRing(Kx);
  Qnew:=Kxy!0;
  C:=Coefficients(Zaxy_to_Kxy(Q,Kxy));
  for i:=1 to #C do
    Qnew:=Qnew+Evaluate(C[i],1/x)*y^(i-1);
  end for;
  Q:=Qnew;
  fun_field:=ext<Kx|Q>;
  binf:=Basis(MaximalOrderFinite(fun_field));
  d:=Degree(Q);
  mat:=ZeroMatrix(Kx,d,d);
  for i:=1 to d do
    for j:=1 to d do
      mat[i,j]:=Eltseq(fun_field!binf[i])[j];
    end for;
  end for;
  return Evaluate(mat,1/x);

end function;


mat_Winf_alternative:=function(Q,W0,Kxy) 

  // Compute the matrix Winf using W0, MaximalOrderInfinite
  // and some linear algebra.

  K:=CoefficientRing(CoefficientRing(Kxy));
  Kx<x>:=RationalFunctionField(K);
  Kxy<y>:=PolynomialRing(Kx);
  d:=Degree(Q);

  Qnew:=Kxy!0;
  C:=Coefficients(Zaxy_to_Kxy(Q,Kxy));
  for i:=1 to #C do
    Qnew:=Qnew+Evaluate(C[i],x)*y^(i-1);
  end for;
  Q:=Qnew;
  
  fun_field:=ext<Kx|Q>;
  binf:=Basis(MaximalOrderInfinite(fun_field));
  Winf:=ZeroMatrix(Kx,d,d);
  for i:=1 to d do
    for j:=1 to d do
      Winf[i,j]:=Eltseq(fun_field!binf[i])[j];;
    end for;
  end for;
  
  W:=W0*Winf^(-1);
  W:=Evaluate(W,1/x);
  Kt<t>:=PolynomialRing(K);
  
  denom:=Kt!1;
  for i:=1 to d do
    for j:=1 to d do
      denom:=LCM(denom,Denominator(W[i,j]));
    end for;
  end for;

  A:=ZeroMatrix(Kt,d,d);
  for i:=1 to d do
    for j:=1 to d do
      A[i,j]:=Kt!(W[i,j]*denom);
    end for;
  end for;  

  S,P1,P2:=SmithForm(A);
  
  S:=ChangeRing(S,Kx);
  S:=S/denom;
  Sinv:=S^(-1);
  for i:=1 to d do
    Sinv[i,i]:=x^(ord_0(Sinv[i,i]));
  end for;

  Winf:=Evaluate(Sinv*P1,1/x)*W0;

  return Winf;

end function;


ddx:=function(f);

  // Differentiate polynomial f(x)(y) with respect to x.

  C:=Coefficients(f);
  for i:=1 to #C do
    C[i]:=Derivative(C[i]);
  end for;
  return Parent(f)!C;
end function;


ddx_mat:=function(A)

  // Differentiate matrix of rational functions.

  for i:=1 to NumberOfRows(A) do
    for j:=1 to NumberOfColumns(A) do
      A[i,j]:=Derivative(A[i,j]);
    end for;
  end for;
  return A;
end function;


ddx_vec:=function(v)

  // Differentiate vector of rational functions.

  for i:=1 to #Eltseq(v) do
    v[i]:=Derivative(v[i]);
  end for;
  return v;
end function;


reduce_mod_Q:=function(f,Q);

  // Eliminate powers of y >= d_x.

  while Degree(f) gt Degree(Q)-1 do
    f:=f-LeadingCoefficient(f)*(Parent(f).1)^(Degree(f)-Degree(Q))*Q;
  end while;
  return f;
end function;


con_mat:=function(Q,Delta,s,K);

  // Compute the connection matrix G.

  d:=Degree(Q);
  Kx<x>:=RationalFunctionField(K);
  Kxy<y>:=PolynomialRing(Kx);
  Delta:=Zax_to_Kx(Delta,Kx);
  Q:=Zaxy_to_Kxy(Q,Kxy);
  s:=Zaxy_to_Kxy(s,Kxy);
  list:=[];
  list[1]:=Kxy!0;
  for i:=2 to d do
    list[i]:=-(i-1)*y^(i-2)*(s/Delta)*ddx(Q);
  end for;
  for i:=1 to #list do
    list[i]:=reduce_mod_Q(list[i],Q);
  end for;
  G:=ZeroMatrix(Kx,d,d);
  for i:=1 to d do
    for j:=1 to d do
      G[i,j]:=Coefficient(list[i],j-1); // G acts on the right on row vectors
    end for;
  end for;
  return(G);
end function;


basis_kernel_mod_pN:=function(A,p,n,m,N)

  // Compute a basis for the kernel of the matrix A modulo p^N, where A is a
  // matrix over the numberfield K defined by m, which is unramified at p.

  // Get rid of p-adic denominator

  val:=Minimum([0] cat [Valuation(j,p) : j in &cat([Eltseq(i) : i in Eltseq(A)])]);
  A:=p^(-val)*A; 
  N:=N-val;

  // Push matrix into a p-adic ring

  Zp:=pAdicQuotientRing(p,N); 
  Zq:=UnramifiedExtension(Zp,m);
  row:=NumberOfRows(A);
  col:=NumberOfColumns(A);
  matpN:=ZeroMatrix(Zq,row,col);
  for i:=1 to row do
    for j:=1 to col do
      matpN[i,j]:=Zq!Eltseq(A[i,j]);
    end for;
  end for;

  // Compute a basis for the kernel

  S,P1,P2:=SmithForm(matpN);
  b:=[];
  for i:=Rank(S)+1 to row do
    Append(~b,P1[i]);
  end for;
  if #b gt 0 then
    b:=RowSequence(HermiteForm(Matrix(b)));
  end if;

  // Push the basis back into the number field K

  K:=BaseRing(A);
  bnew:=[];
  for i:=1 to #b do
    v:=[];
    for j:=1 to #(b[i]) do
      v[j]:=K!Eltseq(b[i][j]);
    end for;
    Append(~bnew,v);
  end for;
  b:=bnew;

  return b;
end function;


basis_kernel_mod_pN_Ki:=function(A,p,n,m,N)

  // Compute a basis for the kernel of the matrix A modulo p^N, where A is a
  // matrix over the numberfield Ki over the Numberfield K defined by m, and K 
  // is unramified at p. Does the linear algebra over K, which works better.

  Ki:=BaseRing(A); 
  K:=BaseRing(Ki);
  
  // Find the matrix over K

  row:=NumberOfRows(A);
  col:=NumberOfColumns(A);
  degKi:=Degree(Ki);
  B:=ZeroMatrix(K,degKi*row,degKi*col);
  
  for i:=1 to row do
    for j:=1 to degKi do
      v:=[Ki|];
      for k:=1 to row do
        v[k]:=Ki!0;
      end for;
      v[i]:=Ki.1^(j-1);
      w:=Vector(v)*A;
      for k:=1 to col do
        for l:=1 to degKi do
          B[(i-1)*degKi+j,(k-1)*degKi+l]:=Eltseq(w[k])[l];
        end for;
      end for;
    end for;
  end for;

  // Compute a basis over K modulo p^N

  bK:=basis_kernel_mod_pN(B,p,n,m,N);

  // Pick subset which is a basis over Ki modulo p^N

  bKi:=[];
  for i:=1 to (#bK div degKi) do
    v:=bK[(i-1)*degKi+1]; // guaranteed to be a basis over Ki?
    w:=[];
    for j:=1 to row do
      w[j]:=Ki!0;
    end for;
    for j:=1 to row do
      for k:=1 to degKi do
        w[j]:=w[j]+(v[(j-1)*degKi+k])*Ki.1^(k-1); 
      end for;
    end for;
    Append(~bKi,w);
  end for;

  return bKi;
end function;


invert_matrix_mod_pN:=function(A,p,n,m,N)

  // Compute the inverse of A modulo p^N, where A is a matrix over 
  // the numberfield K defined by m, which is unramified at p.

  // Get rid of p-adic denominator
  
  val:=Minimum([0] cat [Valuation(j,p) : j in &cat([Eltseq(i) : i in Eltseq(A)])]);
  A:=p^(-val)*A; 
  N:=N-2*val;

  // Push matrix into a p-adic ring

  Zp:=pAdicQuotientRing(p,N); 
  Zq:=UnramifiedExtension(Zp,m);
  row:=NumberOfRows(A);
  matpN:=ZeroMatrix(Zq,row,row);
  for i:=1 to row do
    for j:=1 to row do
      matpN[i,j]:=Zq!Eltseq(A[i,j]);
    end for;
  end for;

  // Invert

  S,P,Q:=SmithForm(matpN);
  maxvalS:=0;
  for i:=1 to row do
    if Valuation(S[i,i]) gt maxvalS then
      maxvalS:=Valuation(S[i,i]);
    end if;
  end for;

  Sinv:=ZeroMatrix(Zq,row,row);
  for i:=1 to row do
    Sinv[i,i]:=(Zq!p)^(maxvalS-Valuation(S[i,i])); // is invariant factor always power of p?
  end for;

  B:=Q*Sinv*P;

  // Push the result back into the number field K

  K:=BaseRing(A); 
  Bnew:=ZeroMatrix(K,row,row);
  for i:=1 to row do
    for j:=1 to row do
      Bnew[i,j]:=K!Eltseq(B[i,j]);
    end for;
  end for;
  B:=Bnew;
  B:=p^(-val-maxvalS)*B;

  return B;
end function;


invert_matrix_mod_pN_Ki:=function(A,p,n,m,N)

  // Compute the inverse of A modulo p^N, where A is a matrix over 
  // the numberfield Ki over the numberfield K defined by m, where
  // K is unramified at p. Does the linear algebra over K, which
  // works better.

  Ki:=BaseRing(A); 
  K:=BaseRing(Ki);

  // Find the matrix over K
  
  row:=NumberOfRows(A);
  degKi:=Degree(Ki);
  B:=ZeroMatrix(K,degKi*row,degKi*row); 
  for i:=1 to row do
    for j:=1 to degKi do
      v:=[Ki|];
      for k:=1 to row do
        v[k]:=Ki!0;
      end for;
      v[i]:=Ki.1^(j-1);
      w:=Vector(v)*A;
      for k:=1 to row do
        for l:=1 to degKi do
          B[(i-1)*degKi+j,(k-1)*degKi+l]:=Eltseq(w[k])[l];
        end for;
      end for;
    end for;
  end for;

  // Invert

  C:=invert_matrix_mod_pN(B,p,n,m,N);

  // Push the result back into Ki

  Cnew:=ZeroMatrix(Ki,row,row);
  for i:=1 to row do
    for j:=1 to row do
      for k:=1 to degKi do
        Cnew[i,j]:=Cnew[i,j]+C[(i-1)*degKi+1,(j-1)*degKi+k]*Ki.1^(k-1);
      end for;
    end for;
  end for;

  C:=Cnew;

  return C;
end function; 


fin_ram_ind:=function(r,G0,Kx)

  // Compute the finite ramification indices

  d:=NumberOfRows(G0);
  K:=BaseRing(Kx);
  rK:=Zax_to_Kx(r,Kx);
  fac:=Factorization(rK);
  M0:=G0*rK;
  e0:=1;
  e0list:=[**];
  resG0list:=[**];
  for i:=1 to #fac do
    e0listi:=[];
    ri:=fac[i][1];
    if Degree(ri) eq 1 then
      Ki:=K;
      s:=-Evaluate(fac[i][1],0);
    else
      Ki<s>:=ext<K|fac[i][1]>; 
    end if;
    resG0:=Evaluate(M0,s)/Evaluate(Derivative(rK),s);
    Append(~resG0list,resG0);
    for j:=1 to d do                                          
      if Determinant(resG0-(1/j)*IdentityMatrix(Ki,d)) eq 0 then 
        Append(~e0listi,j);   
        e0:=Maximum(e0,j);                          
      end if;                                                 
    end for;
    Append(~e0list,e0listi);                                                 
  end for;

  return e0,e0list,resG0list;
end function;


inf_ram_ind:=function(Ginf,Kx)

  // Compute the infinite ramification indices

  d:=NumberOfRows(Ginf);
  K:=BaseRing(Kx);

  resGinf:=-Evaluate((1/Kx.1)*Evaluate(Ginf,1/Kx.1),0);
  einflist:=[];
  for j:=1 to d do
    if Determinant(resGinf-(1/j)*IdentityMatrix(K,d)) eq 0 then
      Append(~einflist,j);
    end if;
  end for;
  if #einflist eq 0 then
    einf:=1;
  else
    einf:=Maximum(einflist);
  end if;  

  return einf,einflist,resGinf;
end function;


jordan_0:=function(p,n,m,r,e0list,resG0list,Kx)

  // Precompute diagonalisations of finite residue matrices

  d:=NumberOfRows(resG0list[1]);
  J0:=[**];
  T0:=[**];
  T0inv:=[**];

  rK:=Zax_to_Kx(r,Kx);
  fac:=Factorisation(rK);
  
  for i:=1 to #fac do
    Ki:=Parent(resG0list[i][1,1]); // what if resG0list[i] is defined over K....
    J0i:=ZeroMatrix(Ki,d,d); 
    //if exactcoho then
      b:=Basis(Kernel(resG0list[i]));
    //else
     // if Degree(fac[i][1]) eq 1 then
     //   b:=basis_kernel_mod_pN(resG0list[i],p,n,m,3*N);
     // else
     //   b:=basis_kernel_mod_pN_Ki(resG0list[i],p,n,m,3*N); 
     // end if;
    //end if;
    cnt:=#b+1;
    for j:=1 to #e0list[i] do 
      for k:=1 to (e0list[i][j]-1) do
      //  if exactcoho then
          b1:=Basis(Kernel(resG0list[i]-k/e0list[i][j]*IdentityMatrix(Ki,d)));  
      //  else
       //   if Degree(fac[i][1]) eq 1 then
       //     b1:=basis_kernel_mod_pN(resG0list[i]-k/e0list[i][j]*IdentityMatrix(Ki,d),p,n,m,3*N);
       //   else
       //     b1:=basis_kernel_mod_pN_Ki(resG0list[i]-k/e0list[i][j]*IdentityMatrix(Ki,d),p,n,m,3*N);
       //   end if;
       // end if;
        for l:=1 to #b1 do
          J0i[cnt,cnt]:=k/e0list[i][j];
          cnt:=cnt+1;
        end for;
        b:=b cat b1;
      end for;
    end for;
    Append(~J0,J0i);
    Append(~T0,Matrix(b));
    //if exactcoho then
      Append(~T0inv,T0[i]^(-1));
   // else 
   //   if Degree(fac[i][1]) eq 1 then
   //     Append(~T0inv,invert_matrix_mod_pN(T0[i],p,n,m,3*N));
    //  else
    //    Append(~T0inv,invert_matrix_mod_pN_Ki(T0[i],p,n,m,3*N));
    //  end if;
    //end if;
  end for;

  return J0,T0,T0inv; 
end function;


jordan_inf:=function(p,n,m,einflist,resGinf)

  // Precompute diagonalisation of infinite residue matrix

  d:=NumberOfRows(resGinf);
  K:=Parent(resGinf[1,1]);
  Jinf:=ZeroMatrix(K,d,d);
 
  //if exactcoho then 
  b:=Basis(Kernel(resGinf));
  //else
 //   b:=basis_kernel_mod_pN(resGinf,p,n,m,3*N); 
 // end if;
  cnt:=#b+1;
  for i:=1 to #einflist do
    for j:=1 to (einflist[i]-1) do
      //if exactcoho then
        b1:=Basis(Kernel(resGinf-j/einflist[i]*IdentityMatrix(K,d)));
      //else
      //  b1:=basis_kernel_mod_pN(resGinf-j/einflist[i]*IdentityMatrix(K,d),p,n,m,3*N);
      //end if;
      for k:=1 to #b1 do
        Jinf[cnt,cnt]:=j/einflist[i];
        cnt:=cnt+1;
      end for;
      b:=b cat b1;
    end for;
  end for;
  Tinf:=Matrix(b);
  //if exactcoho then
    Tinfinv:=Tinf^(-1);
 // else
  //  Tinfinv:=invert_matrix_mod_pN(Tinf,p,n,m,3*N);
  //end if;
  return Jinf,Tinf,Tinfinv; 
end function;


res_0:=function(w,QK,rK,J0,T0inv)

  // Compute res_0(\sum w_i b^0_i dx/r).

  d:=Degree(QK);
  fac:=Factorization(rK);
  reslist:=[];
  for i:=1 to #fac do
    Ki:=Parent(T0inv[i][1,1]);
    if Degree(fac[i][1]) eq 1 then
      s:=Ki!(-Coefficient(fac[i][1],0));
    else
      s:=Ki.1;
    end if;
    v:=Vector(Evaluate(w,s));
    v:=v*T0inv[i];
    for j:=1 to d do
      if J0[i][j,j] eq 0 then
        if Degree(fac[i][1]) eq 1 then   
          reslist:=reslist cat [v[j]]; 
        else                             
          reslist:=reslist cat Eltseq(v[j]);
        end if;
      end if; 
    end for; 
  end for;  

  return Vector(reslist);

end function;


val_Kttinv_d:=function(v)

  // Compute the valuation of an element of Kttinvd.  

  val:=Valuation(v[1]);
  for i:=2 to #Eltseq(v) do
    if Valuation(v[i]) lt val then
      val:=Valuation(v[i]);
    end if;
  end for;
  return val;
end function;


res_inf:=function(w,QK,rK,W0,Winf,Ginf,Jinf,Tinfinv,Kxy)

  // Compute res_inf(\sum w_i b^0_i dx/r).

  Kx:=BaseRing(Kxy);
  K:=BaseRing(Kx); 
  d:=Degree(QK);
  Kd:=RSpace(K,d);
  degr:=Degree(rK);
  Kttinv<t>:=LaurentSeriesRing(K);
  Kttinvd:=RSpace(Kttinv,d);

  W:=Winf*W0^(-1);
  Winv:=W^(-1);
  w:=Kttinvd!Evaluate(w,t^(-1));
  w:=w*Evaluate(Winv,t^(-1));

  resGinf:=-Evaluate((1/Parent(Winf[1,1]).1)*Evaluate(Ginf,1/Parent(Winf[1,1]).1),0);

  // reduce to a cohomologous 1-form that is logarithmic at all points lying over x=inf:

  while val_Kttinv_d(w) lt -degr+1 do
    m:=-val_Kttinv_d(w)-degr+1;
    mat:=resGinf-m*IdentityMatrix(K,d);
    rhs:=Kd!0;
    for i:=1 to d do
      rhs[i]:=rhs[i]+Coefficient(-w[i],-m-degr+1)/LeadingCoefficient(rK);
    end for;
    vbar:=rhs*mat^(-1);
    w:=w-ChangeRing(vbar,Kttinv)*t^(-m)*Evaluate(rK*Ginf,t^(-1))-Evaluate(rK,1/t)*m*t^(1-m)*ChangeRing(vbar,Kttinv);  
  end while;

  // now sum w_i b^{inf}_i dx/r is logarithmic at all points lying over x=inf

  w:=w*t^(degr-1);
  v:=Kd!0;
  for i:=1 to d do
    v[i]:=Coefficient(w[i],0);
  end for;

  // project v onto the eigenspace of res_Ginf of eigenvalue 0 

  v:=v*Tinfinv;

  res:=[];
  for i:=1 to d do
    if Jinf[i,i] eq 0 then
      Append(~res,v[i]);
    end if;
  end for;
  
  return Vector(res);
end function;

polys_to_vec:=function(polys,degx,K);

  // Converts a sequence of polynomials to a vector  

  dim:=#polys*(degx+1);
  v:=[];
  cnt:=1;
  for i:=1 to #polys do
    for j:=0 to degx do
      v[cnt]:=Coefficient(polys[i],j);
      cnt:=cnt+1;
    end for;
  end for;

  V:=VectorSpace(K,dim);

  return V!v;
end function;



basis_coho:=function(Q,p,r,W0,Winf,G0,Ginf,J0,Jinf,T0inv,Tinfinv,useU,basis0,basis1,basis2,K,Kx,Kxy)
  
  // Compute a basis for H^1(X).

  //Qx<x>:=PolynomialRing(RationalField());
 // Qxy<y>:=PolynomialRing(Qx);
  x:=Kx.1;
  y:=Kxy.1;
  d:=Degree(Q);
  Kxd:=RSpace(Kx,d);
  degr:=Degree(r);
  pIdeal:=ideal<RingOfIntegers(K)|p>;
 
  rK:=Zax_to_Kx(r,Kx);
  QK:=Zaxy_to_Kxy(Q,Kxy);
 
  W:=Winf*W0^(-1);

  Winv:=W^(-1);
  ord0W:=ord_0_mat(W);
  ordinfW:=ord_inf_mat(W);
  ord0Winv:=ord_0_mat(Winv);
  ordinfWinv:=ord_inf_mat(Winv);

  // Compute a basis for E0

  deg_bound_E0:=degr-ord0W-ordinfW-2; 
  basisE0:=[];
  for i:=0 to d-1  do 
    for j:=0 to deg_bound_E0 do
      basisE0:=Append(basisE0,[i,j]);
    end for;
  end for;
  dimE0:=#basisE0;
  E0:=VectorSpace(K,dimE0);

  // Compute a matrix with kernel (E0 intersect Einf).

  matE0nEinf:=ZeroMatrix(K,dimE0,d*(-ordinfW-ordinfWinv));
  for i:=1 to dimE0 do
    temp:=RowSequence(x^(basisE0[i][2])*Winv)[basisE0[i][1]+1];
    for j:=0 to d-1 do
      for k:=0 to (-ordinfW-ordinfWinv-1) do
        matE0nEinf[i,j*(-ordinfW-ordinfWinv)+k+1]:=Coefficient(Numerator((Parent(W[1,1]).1)^(-ord0Winv)*temp[j+1]),k-ord0W-ord0Winv+degr-1);
      end for;
    end for;
  end for;  

  E0nEinf:=Kernel(matE0nEinf);

  // Compute a matrix with kernel the elements of E0 logarithmic at infinity.

  matlogforms:=ZeroMatrix(K,dimE0,d*(-ord0W-ordinfW-ordinfWinv-1));
  for i:=1 to dimE0 do
    temp:=RowSequence(x^(basisE0[i][2])*Winv)[basisE0[i][1]+1];
    for j:=0 to d-1 do
      for k:=0 to (-ord0W-ordinfW-ordinfWinv-2) do
        matlogforms[i,j*(-ord0W-ordinfW-ordinfWinv-1)+k+1]:=Coefficient(Numerator((Parent(W[1,1]).1)^(-ord0Winv)*temp[j+1]),k-ord0Winv+degr);
      end for;
    end for;
  end for;

  logforms:=E0nEinf meet Kernel(matlogforms);

  // Compute the finite residues.
  
  w:=Kxd!0;
  w[1]:=1;
  res0dim:=Dimension(Parent(res_0(w,QK,rK,J0,T0inv)));
  matres0:=ZeroMatrix(K,dimE0,res0dim); 
  for i:=1 to dimE0 do
    w:=Kxd!0;
    w[basisE0[i][1]+1]:=x^(basisE0[i][2]);
    coefs:=res_0(w,QK,rK,J0,T0inv); 
    for j:=1 to res0dim do
        matres0[i,j]:=coefs[j];
    end for;
  end for;

  // Compute the infinite residues. 

  w:=Kxd!0;
  w[1]:=1;
  resinfdim:=Dimension(Parent(res_inf(w,QK,rK,W0,Winf,Ginf,Jinf,Tinfinv,Kxy)));
  matresinf:=ZeroMatrix(K,dimE0,resinfdim); 
  for i:=1 to dimE0 do
    w:=Kxd!0;
    w[basisE0[i][1]+1]:=x^(basisE0[i][2]);
    coefs:=res_inf(w,QK,rK,W0,Winf,Ginf,Jinf,Tinfinv,Kxy); 
    for j:=1 to resinfdim do
        matresinf[i,j]:=coefs[j];
    end for;
  end for;

  forms2ndkind:=Kernel(matres0) meet Kernel(matresinf);
  cocyc:=E0nEinf meet forms2ndkind;
  forms1stkind:=logforms meet forms2ndkind;
  
  // Compute a matrix with kernel (B0 intersect Binf)

  deg_bound_B0:=-ord0W-ordinfW-1;
  basisB0:=[];
  for i:=0 to d-1  do 
    for j:=0 to deg_bound_B0 do
      basisB0:=Append(basisB0,[i,j]);
    end for;
  end for;
  dimB0:=#basisB0;
  B0:=VectorSpace(K,dimB0);

  matB0nBinf:=ZeroMatrix(K,dimB0,d*(-ordinfW-ordinfWinv));
  for i:=1 to dimB0 do
    power_x:=basisB0[i][2];
    power_y:=basisB0[i][1];
    temp:=RowSequence(x^(power_x)*Winv)[power_y+1];
    for j:=0 to d-1 do
      for k:=0 to (-ordinfW-ordinfWinv-1) do
        matB0nBinf[i,j*(-ordinfW-ordinfWinv)+k+1]:=Coefficient(Numerator((Parent(W[1,1]).1)^(-ord0Winv)*temp[j+1]),k-ord0W-ord0Winv);
      end for;
    end for;
  end for;

  // Compute d(B0 intersect Binf).

  B0nBinf:=Kernel(matB0nBinf);
  basisB0nBinf:=Basis(B0nBinf);
  dimB0nBinf:=#basisB0nBinf;  

  list:=[];
  for i:=1 to dimB0nBinf do
    vec:=basisB0nBinf[i];
    vecKxd:=Kxd!0;
    for j:=1 to dimB0 do
      vecKxd[basisB0[j][1]+1]:=vecKxd[basisB0[j][1]+1]+vec[j]*x^(basisB0[j][2]);
    end for;
    vecKxd:=vecKxd*ChangeRing(rK*G0,Kx)+rK*ddx_vec(vecKxd);
    coefs:=[];
    for j:=1 to dimE0 do
      power_x:=basisE0[j][2];
      power_y:=basisE0[j][1];
      coefs[j]:=Coefficient(vecKxd[power_y+1],power_x);  
    end for;
    list:=Append(list,E0!coefs);
  end for;
  matd:=Matrix(list);

  // Compute bases

  cobound:=sub<cocyc|list>;

  if basis0 eq [] then
    b0:=Basis(forms1stkind);
  else
    b0:=[];
    for i:=1 to #basis0 do
      b0[i]:=polys_to_vec(basis0[i],deg_bound_E0,K);
    end for;
  end if;  

  b5:=[];
  for i:=2 to dimB0nBinf do
    b5:=Append(b5,E0!list[i]);
  end for;
  
  dualspace:=Complement(cocyc,forms1stkind+cobound); // Take the dual w.r.t. cup product? Right now just any complement of forms of the 1st kind in H^1(X).
  if basis1 eq [] then
    b1:=Basis(dualspace);
  else
    b1:=[];
    for i:=1 to #basis1 do
      b1[i]:=polys_to_vec(basis1[i],deg_bound_E0,K);
    end for;
  end if;  

  basisH1X:=b0 cat b1;
  dimH1X:=#basisH1X;

  finiteregularlogarithmic:=logforms meet Kernel(matres0); // 1-forms that generate H^1(Y) over H^1(X), where Y=X-x^{-1}(\infty)
  H1YmodH1X:=Complement(finiteregularlogarithmic,forms1stkind);

  if basis2 eq [] then
    b2:=Basis(H1YmodH1X);
  else
    b2:=[];
    for i:=1 to #basis2 do
      b2[i]:=polys_to_vec(basis2[i],deg_bound_E0,K);
    end for;
  end if;

  b3:=Basis(Complement(E0nEinf,cocyc+H1YmodH1X));

  b4:=Basis(Complement(E0,E0nEinf));  
  
  b:=b0 cat b1 cat b2 cat b3 cat b4 cat b5;

  dimH1U:=#b0+#b1+#b2+#b3;

  if useU then
    dim:=dimH1U;
  else
    dim:=dimH1X;
  end if;

  for i:=1 to dim do
    valdenom:=0;
    for j:=1 to dimE0 do 
      valdenom:=Minimum(valdenom,Valuation(b[i][j],pIdeal));
    end for;
    b[i]:=p^(-valdenom)*b[i];
  end for; 

  matb:=Matrix(b);
  quo_map:=matb^(-1);

  integrals:=[Kxd|];
  for i:=2 to dimB0nBinf do
    vec:=Kxd!0;
    for j:=1 to dimB0 do
      vec[basisB0[j][1]+1]:=vec[basisB0[j][1]+1]+(basisB0nBinf[i][j])*x^(basisB0[j][2]);
    end for;
    integrals:=Append(integrals,LeadingCoefficient(r)*vec); // factor lc(r) here, since working with dx/z basis instead of dx/r
  end for;

  Zaxd:=RSpace(Zax,d);
  basis:=[Zaxd|];
  
  for i:=1 to dim do
    vec:=Zaxd!0;
    for j:=1 to dimE0 do
      vec[basisE0[j][1]+1]:=vec[basisE0[j][1]+1]+(Za!Eltseq(b[i][j]))*(Zax.1)^(basisE0[j][2]);
    end for;
    basis:=Append(basis,vec);
  end for;

  return basis,integrals,quo_map;
    
end function;

