WeakValuation:=function(f)
  coeffs:=Coefficients(f);
  m:=Minimum([i : i in [0..Degree(f)] | not IsWeaklyZero(coeffs[i+1])]);
  return m;
end function;

// TODO: go through Tuitman's proof
max_prec:=function(Q,p,N,g,W0,Winf,e0,einf)

  // Compute the p-adic precision required for provable correctness

  d:=Degree(Q);
  W:=Winf*W0^(-1);
  
  Nmax:=N+Floor(log(p,-p*(ord_0_mat(W)+1)*einf));
  while (Nmax-Floor(log(p,p*(Nmax-1)*e0))-Floor(log(p,-(ord_inf_mat(W^(-1))+1)*einf)) lt N) do 
    Nmax:=Nmax+1;
  end while;

  Nmax:=Maximum(Nmax,2); 

  return Nmax; // precision to be used in the computations
end function;


frobmatrix:=function(Q,p,n,m,N,Nmax,g,r,W0,Winf,G0,Ginf,frobmatb0r,red_list_fin,red_list_inf,basis,integrals,quo_map,verbose,Kx);

  // Compute the matrix of F_p on H^1(X) mod p^N with respect to 'basis'.
  K:=BaseRing(Kx);
  F:=ZeroMatrix(K,#basis,#basis);  
  f0list:=[];
  finflist:=[];
  fendlist:=[];

  for i:=1 to #basis do
    dif:=frobenius(basis[i],Q,p,n,m,Nmax,r,frobmatb0r);
    dif:=convert_to_Kxzzinvd(dif,Q,Kx);

    coefs,f0,finf,fend:=reduce_with_fs(dif,Q,p,N,Nmax,r,W0,Winf,G0,Ginf,red_list_fin,red_list_inf,basis,integrals,quo_map,Kx);

    for j:=1 to #basis do
      F[i,j]:=coefs[j];
    end for;
    
    f0list:=Append(f0list,f0);
    finflist:=Append(finflist,finf);
    fendlist:=Append(fendlist,fend);

    if verbose then
      printf ".";
    end if;

  end for;
 
  return F,f0list,finflist,fendlist;

end function;


coleman_data:=function(Q,p,m,N:useU:=false,basis0:=[],basis1:=[],basis2:=[],verbose:=false,W0:=0,Winf:=0)


  // Takes a polynomial Q in three variables a,x,y over the rationals which is monic in y, where a represents a generator of K over Q.
  // Returns the Coleman data of (the projective nonsingular model of) the curve defined
  // by Q at p to p-adic precision N.
  
  exactcoho:=false;

  if not IsPrime(p) then
    error "p is not prime";
  end if;

  if not IsIrreducible(Q) then
    error "Curve is not irreducible";
  end if;
  
  //Kxy:=Parent(Q);
 // Kx:=Parent(Coefficients(Q)[1]);
 // K:=BaseRing(Kx);
 // OK:=RingOfIntegers(K);
 // n:=InertiaDegree(ideal<OK|p>);
  //m:=DefiningPolynomial(FiniteField(q));
  
  K:=NumberField(m);
  if Type(K) eq Type(Rationals()) then
	K:=QNF();
  end if;
  OK := Integers(K);
  pfrs := Decomposition(OK, p);
  n, pos := Min([Degree(pfr[1]) : pfr in pfrs | pfr[2] eq 1]); // pick an unramified prime ideal pfr of minimal degree n
  pfr := pfrs[pos][1];
  Kp, ip := Completion(K, pfr);
  //OKp := RingOfIntegers(Kp);
  //Fq, mod_pfr := ResidueClassField(OKp);
  Kx<x> := PolynomialRing(K);
  Kxy<y> := PolynomialRing(Kx);
  Kxyz<z> := LaurentSeriesRing(Kxy);

  ResetMaximumMemoryUsage();
  t0:=Cputime();

  if verbose then
    print "Computing Delta,r,s,W^0,W^inf,G:";
  end if;

  d:=Degree(Q);
  g:=genus(Q,p,n,m);
  
  r,Delta,s:=auxpolys(Q,p,n,K);

  if W0 eq 0 then
    W0:=mat_W0(Q,Kxy);
  end if;
  if Winf eq 0 then
    Winf:=mat_Winf(Q,Kxy);
  end if;
  W0inv:=W0^(-1); 
  W:=Winf*W0inv;
  Winfinv:=Winf^(-1);

  assert Type(LeadingCoefficient(Delta)) cmpeq RngIntElt;
  // TODO: use the "smooth" function corretly. We need a function to push polynomials from Zax to Kpx. I think Zax_to_Kpx will not work for us, since our prime is not inert and hence the minimal polynomial of K/Q is not one of Kp/Qp (the latter might have lower degree).
  if (FiniteField(p)!LeadingCoefficient(Delta) eq 0) or (Degree(r) lt 1) or (not smooth(r,Kp)) or (not (is_integral(W0,p,n) and is_integral(W0inv,p,n) and is_integral(Winf,p,n) and is_integral(Winfinv,p,n))) then
    error "bad prime";
  end if;

  G:=con_mat(Q,Delta,s,K);
  G0:=W0*Evaluate(G,Parent(W0[1,1]).1)*W0^(-1)+ddx_mat(W0)*W0^(-1);
  Ginf:=Winf*Evaluate(G,Parent(Winf[1,1]).1)*Winf^(-1)+ddx_mat(Winf)*Winf^(-1);

  e0,e0list,resG0list := fin_ram_ind(r,G0,Kx);
  einf,einflist,resGinf := inf_ram_ind(Ginf,Kx);
  Jinf,Tinf,Tinfinv:=jordan_inf(p,n,m,einflist,resGinf);
  J0,T0,T0inv:=jordan_0(p,n,m,r,e0list,resG0list,Kx);
  //e0,einf:=ram(J0,Jinf);
 
  delta:=Floor(log(p,-(ord_0_mat(W)+1)*einf))+Floor(log(p,(Floor((2*g-2)/d)+1)*einf));

  if verbose then
    print "Time (s) :    ", Cputime(t0);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing basis H^1(X):"; 
  end if;

  basis,integrals,quo_map:=basis_coho(Q,p,r,W0,Winf,G0,Ginf,J0,Jinf,T0inv,Tinfinv,useU,basis0,basis1,basis2,K,Kx,Kxy);
  
 
  if verbose then 
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing Frobenius lift:";
  end if;

  Nmax:=max_prec(Q,p,N,g,W0,Winf,e0,einf);
  SetPrecision(Kp, Nmax);

  frobmatb0r:=froblift(Q,p,n,m,Nmax-1,r,Delta,s,W0);

  if verbose then 
    print "Time (s) :    ", Cputime(t);
    print "Memory (MiB) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  if verbose then
    print "Computing reduction matrices:";
  end if;

  red_list_fin,red_list_inf:=red_lists(Q,p,n,m,Nmax,r,W0,Winf,G0,Ginf,e0,einf,J0,Jinf,T0,Tinf,T0inv,Tinfinv,Kx);

  if verbose then
    print "Time (s) :    ", Cputime(t);
    print "Memory (MiB) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing Frobenius matrix:";
  end if;

  F,f0list,finflist,fendlist:=frobmatrix(Q,p,n,m,N,Nmax,g,r,W0,Winf,G0,Ginf,frobmatb0r,red_list_fin,red_list_inf,basis,integrals,quo_map,verbose,Kx);

  if verbose then
    print "";
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  rK:=Zax_to_Kx(r,Kx);
  rK:=rK/LeadingCoefficient(rK);

  // formatting the output into a record:
 
  format:=recformat<Q,p,m,N,g,W0,Winf,r,Delta,s,G0,Ginf,e0,einf,delta,basis,quo_map,integrals,F,f0list,finflist,fendlist,Nmax,red_list_fin,red_list_inf,minpolys,K,n,Kx,Kxy,Kp,ip,rK>;
  out:=rec<format|>;
  out`Q:=Q; out`p:=p; out`m:=m; out`N:=N; out`g:=g; out`W0:=W0; out`Winf:=Winf; out`r:=r; out`Delta:=Delta; out`s:=s; out`G0:=G0; out`Ginf:=Ginf; 
  out`e0:=e0; out`einf:=einf; out`delta:=delta; out`basis:=basis; out`quo_map:=quo_map; out`integrals:=integrals; out`F:=F; out`f0list:=f0list; 
  out`finflist:=finflist; out`fendlist:=fendlist; out`Nmax:=Nmax; out`red_list_fin:=red_list_fin; out`red_list_inf:=red_list_inf;
  out`K:=K; out`n:=n; out`Kx:=Kx; out`Kxy:=Kxy; out`Kp:=Kp; out`ip:=ip; out`rK:=rK;

  return out;
end function;


set_point:=function(x0,y0,data)

  // Constructs a point from affine coordinates x0,y0. 

  Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0;
  K:=data`K; Kp:=data`Kp; n:=data`n; Kx:=data`Kx; r:=data`r;
  d:=Degree(Q);
  rK:=Zax_to_Kx(r,Kx);
  Kpx:=PolynomialRing(Kp);
  x0:=Kp!x0; y0:=Kp!y0;

  if Valuation(x0) lt 0 then
    error "x0 has negative valuation";
  end if;
  
  if (not(W0 eq IdentityMatrix(BaseRing(W0),d))) and (Valuation(Evaluate(rK,x0)) gt 0) then
    error "W0 is not the identity and r(x0) is zero mod p";
  end if;
  
  format:=recformat<x,b,inf,xt,bt,index>;
  P:=rec<format|>;
  P`inf:=false;
  P`x:=x0;

  y0powers:=[];
  y0powers[1]:=Kp!1;
  for i:=2 to d do
    y0powers[i]:=(y0)^(i-1);
  end for;
  y0powers:=Vector(y0powers);
  W0Kp:=matrix_push_to_Kp(W0,d,Kpx);
  W0x0:=Transpose(Evaluate(W0Kp,x0));

  P`b:=Eltseq(y0powers*W0x0); // the values of the b_i^0 at P

  return P;
end function;


set_bad_point:=function(x,b,inf,data)

  Q:=data`Q; p:=data`p; N:=data`N; 
  K:=data`K; Kp:=data`Kp; d:=Degree(Q);

  format:=recformat<x,b,inf,xt,bt,index>;
  P:=rec<format|>;
  P`inf:=inf;
  P`x:=Kp!x;
  P`b:=[Kp!b[i]:i in [1..d]];

  return P; 

end function;


is_bad:=function(P,data)

  // check whether the point P is bad

  x0:=P`x; r:=data`r; 
  Kp:=Parent(x0);
  Kpx:=PolynomialRing(Kp);
  rKp:=Zax_to_Kpx(r,Kpx);

  if P`inf then // infinite point
    return true;
  end if;

  if Valuation(Evaluate(rKp,x0)) gt 0 then // finite bad point
    return true;
  end if;

  return false;
  
end function;


is_very_bad:=function(P,data)

  // check whether the point P is very bad

  x0:=P`x; r:=data`r; N:=data`N; 
  Kp:=Parent(x0);
  Kpx:=PolynomialRing(Kp);
  rKp:=Zax_to_Kpx(r,Kpx);

  if P`inf then // infinite point
    if Valuation(x0) ge N then // infinite very bad point
      return true;
    end if;
  else // finite point
    if Valuation(Evaluate(rKp,x0)) ge N then // finite very bad point
      return true;
    end if;
  end if;

  return false;

end function;


lie_in_same_disk:=function(P1,P2,data)

  // check whether two points P1,P2 lie in the same residue disk

  x1:=P1`x; b1:=P1`b; x2:=P2`x; b2:=P2`b; Q:=data`Q;
  d:=Degree(Q);
  
  if P1`inf ne P2`inf then // one point infinite, other not
    return false;
  else
    for i:=1 to d do
        if Valuation(b1[i]-b2[i]) lt 1 then
          return false;
        end if;
      end for;
      return true;
  end if;
 
end function;


minpoly:=function(f1,f2)

  // computes the minimum polynomial of f2 over K(f1), where
  // f1,f2 are elements of a 1 dimensional function field over K

  FF:=Parent(f1);
  K:=BaseRing(FF);

  d:=5;  

  done:=false;
  while not done do

    S:=[];
    for i:=0 to d do
      for j:=0 to d do
        S:=Append(S,f1^j*f2^i);
      end for;
    end for;

    denom:=1;
    for i:=1 to #S do
      E:=Eltseq(S[i]);
      for j:=1 to #E do
        denom:=LCM(denom,Denominator(E[j]));
      end for;
    end for;
  
    maxdeg:=0;
    for i:=1 to #S do
      E:=Eltseq(S[i]);
      for j:=1 to #E do
        deg:=Degree(Numerator(denom*E[j]));
        if  deg gt maxdeg then
          maxdeg:=deg;
        end if;
      end for;
    end for;

    T:=[];
    for i:=1 to #S do
      E:=Eltseq(S[i]);
      v:=[];
      for j:=1 to #E do
        for k:=0 to maxdeg do
          v[(j-1)*(maxdeg+1)+k+1]:=Coefficient(Numerator(denom*E[j]),k);
        end for;  
      end for;
      T:=Append(T,v);
    end for;

    b:=Basis(NullSpace(Matrix(T)));  

    if #b gt 0 then
      done:=true;
      R:=b[1];
    else
      d:=d+3;
    end if;
  
  end while;

  Kx:=PolynomialRing(K);
  Kxy:=PolynomialRing(Kx);
  poly:=Kxy!0;
  for i:=0 to d do
    for j:=0 to d do
      poly:=poly+R[i*(d+1)+j+1]*Kx.1^j*Kxy.1^i;
    end for;
  end for;

  fac:=Factorisation(poly);

  for i:=1 to #fac do
    factor:=fac[i][1];
    test:=FF!0;
    for j:=0 to Degree(factor) do
      test:=test+Evaluate(Coefficient(factor,j),f1)*f2^j;
    end for;
    if test eq 0 then
      poly:=factor;
    end if;
  end for;

  return poly;
 
end function;


update_minpolys:=function(data,inf,index);

  // TODO comment

  Q:=data`Q; W0:=data`W0; Winf:=data`Winf; 
  d:=Degree(Q); Kxy:=data`Kxy; K:=data`K;

  if not assigned data`minpolys then
    data`minpolys:=[ZeroMatrix(Kxy,d+2,d+2),ZeroMatrix(Kxy,d+2,d+2)];
  end if;
  minpolys:=data`minpolys;

  Kt:=RationalFunctionField(K); Kty:=PolynomialRing(Kt);
  
  QK:=Zaxy_to_Kxy(Q,Kxy);
  f:=Kty!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(QK,i)) do
      f:=f+Coefficient(Coefficient(QK,i),j)*Kty.1^i*Kt.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  if inf then 
    W:=Winf;
  else
    W:=W0;
  end if;

  bfun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+W[i,j]*FF.1^(j-1);
    end for;
    bfun[i]:=bi;
  end for;

  if inf then // b=b^{infty}
    if index eq 0 then
       for i:=1 to d do
         if minpolys[2][1,i+1] eq 0 then
           minpolys[2][1,i+1]:=minpoly(FF!(1/Kt.1),bfun[i]);
         end if;
       end for;
    else
      if minpolys[2][index+1,1] eq 0 then
        minpolys[2][index+1,1]:=minpoly(bfun[index],FF!(1/Kt.1));
      end if;
      for i:=1 to d do
        if minpolys[2][index+1,i+1] eq 0 then
          minpolys[2][index+1,i+1]:=minpoly(bfun[index],bfun[i]);
        end if;
      end for;
    end if;
  else // b=b^0
    if index eq 0 then
      for i:=1 to d do
        if minpolys[1][1,i+1] eq 0 then
          minpolys[1][1,i+1]:=minpoly(FF!Kt.1,bfun[i]);
        end if;
      end for;
    else
      if minpolys[1][index+1,1] eq 0 then
        minpolys[1][index+1,1]:=minpoly(bfun[index],FF!Kt.1);
      end if;
      for i:=1 to d do
        if minpolys[1][index+1,i+1] eq 0 then
          minpolys[1][index+1,i+1]:=minpoly(bfun[index],bfun[i]);
        end if;
      end for;
    end if;
  end if;

  data`minpolys:=minpolys;

  return data;

end function;


frobenius_pt:=function(P,data);

  // Computes the image of P under Frobenius

  x0:=P`x; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; Kp:=data`Kp;
  K:=data`K; Kxy:=data`Kxy; ip:=data`ip;
  n:=data`n;
  q:=p^n;
  d:=Degree(Q);  Kpy:=PolynomialRing(Kp);
    
  W0Kp:=matrix_push_to_Kp(W0,d,Kpy);
  
  x0q:=x0^q;
  b:=P`b;

  Kt:=RationalFunctionField(K); Kty:=PolynomialRing(Kt);
  QK:=Zaxy_to_Kxy(Q,Kxy);
  f:=Kty!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(QK,i)) do
      f:=f+Coefficient(Coefficient(QK,i),j)*Kty.1^i*Kt.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  if not is_bad(P,data) then // finite good point
    
    W0invx0:=Transpose(Evaluate(W0Kp^(-1),x0));
    
    ypowers:=Vector(b)*W0invx0;
    y0:=ypowers[2];
  
    C:=[Kx_to_Kpt(c,ip,Kpy) : c in Coefficients(QK)];
    D:=[];
    for i:=1 to #C do
      D[i]:=Evaluate(C[i],x0q);
    end for;
    fy:=Kpy!D;

    y0q:=HenselLift(fy,y0^q); // Hensel lifting
  
    y0qpowers:=[];
    y0qpowers[1]:=Kp!1;
    for i:=2 to d do
      y0qpowers[i]:=y0q^(i-1);
    end for;
    y0qpowers:=Vector(y0qpowers);

    W0x0:=Transpose(Evaluate(W0Kp,x0));
  
    b:=Eltseq(y0qpowers*W0x0);

  elif P`inf then // infinite point
    error "INFINITE POINTS NOT IMPLEMENTED!";
    /*
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+Winf[i,j]*FF.1^(j-1);
      end for;

      if assigned data`minpolys and data`minpolys[2][1,i+1] ne 0 then
        poly:=data`minpolys[2][1,i+1];
      else
        poly:=minpoly(FF!(1/Kt.1),bi);
      end if;

      C:=[Kx_to_Kpt(c,ip,Kpy) : i in Coefficients(poly)];
      D:=[];
      for i:=1 to #C do
        D[i]:=Evaluate(C[i],x0q); 
      end for;
      fy:=Kpy!D;

      fac:=Factorisation(fy); // Roots has some problems that Factorisation does not
      zeros:=[];
      for j:=1 to #fac do
        if Degree(fac[j][1]) eq 1 then
          zeros:=Append(zeros,-Coefficient(fac[j][1],0)/Coefficient(fac[j][1],1));
        end if;
      end for;
      
      done:=false;
      j:=1;
      while not done and j le #zeros do
        if Valuation(zeros[j]-b[i]^p) gt p then // was gt 0 before 
          done:=true;
          b[i]:=zeros[j];
        end if;
        j:=j+1;
      end while;
      if not done then
        error "Frobenius does not converge at P";
      end if;
    end for;
   */
  else // finite bad point
   error "BAD POINTS NOT SUPPORTED!";
   /*
   for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+W0[i,j]*FF.1^(j-1);
      end for;

      if assigned data`minpolys and data`minpolys[1][1,i+1] ne 0 then
        poly:=data`minpolys[1][1,i+1];
      else
        poly:=minpoly(FF!Qt.1,bi);
      end if;

      C:=Coefficients(poly);
      D:=[];
      for i:=1 to #C do
        D[i]:=Evaluate(C[i],x0p); 
      end for;
      fy:=Ky!D;

      fac:=Factorisation(fy); // Roots has some problems that Factorisation does not
      zeros:=[];
      for j:=1 to #fac do
        if Degree(fac[j][1]) eq 1 then
          zeros:=Append(zeros,-Coefficient(fac[j][1],0)/Coefficient(fac[j][1],1));
        end if;
      end for;

      done:=false;
      j:=1;
      while not done and j le #zeros do
        if Valuation(zeros[j]-b[i]^p) gt p then
          done:=true;
          b[i]:=zeros[j];
        end if;
        j:=j+1;
      end while;
      if not done then
        error "Frobenius does not converge at P";
      end if;
    end for;
    */
  end if;
  
    P`x:=x0q;
    P`b:=b;
    delete P`xt;
    delete P`bt;
    delete P`index;

  return P;
end function;


teichmueller_pt:=function(P,data)

  // Compute the Teichmueller point in the residue disk at a good point P

  x0:=P`x; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q); K:=Parent(x0); Ky:=PolynomialRing(K);

  if is_bad(P,data) then
    error "Point is bad";
  end if;

  x0new:=K!TeichmuellerLift(FiniteField(p)!x0,pAdicQuotientRing(p,N)); 
  b:=P`b; 
  W0invx0:=Transpose(Evaluate(W0^(-1),x0));
  ypowers:=Vector(b)*W0invx0;
  y0:=ypowers[2];
  
  C:=Coefficients(Q);
  D:=[];
  for i:=1 to #C do
    D[i]:=Evaluate(C[i],x0new);
  end for;
  fy:=Ky!D;

  y0new:=HenselLift(fy,y0); // Hensel lifting
  y0newpowers:=[];
  y0newpowers[1]:=K!1;
  for i:=2 to d do
    y0newpowers[i]:=y0newpowers[i-1]*y0new;
  end for;
  y0newpowers:=Vector(y0newpowers);

  W0x0:=Transpose(Evaluate(W0,x0));
  b:=Eltseq(y0newpowers*W0x0);

  P`x:=x0new;
  P`b:=b;
  delete P`xt;
  delete P`bt;
  delete P`index;

  return P;

end function;


local_data:=function(P,data)

  // For a point P, returns the ramification index of the map x on the residue disk at P

  Q:=data`Q; p:=data`p; W0:=data`W0; Winf:=data`Winf; x0:=P`x; b:=P`b; d:=Degree(Q); n:=data`n;
  Kp:=data`Kp;

  if not is_bad(P,data) then
    eP:=1;
    index:=0;
    return eP,index;
  else     
  Op:=RingOfIntegers(Kp); Fp:=ResidueClassField(Op); Fpx:=RationalFunctionField(Fp); Fpxy:=PolynomialRing(Fpx);
    f:=Fpxy!0;
    C:=Coefficients(Q);
    for i:=1 to #C do
      D:=Coefficients(C[i]);
      for j:=1 to #D do
        E:=Coefficients(D[j]);
        for k:=1 to #E do
          f:=f+(Fp!E[k])*Fp.1^(k-1)*Fpxy.1^(i-1)*Fpx.1^(j-1);
        end for;
      end for;
    end for;  
    FFp:=FunctionField(f); // function field of curve mod p
    
    if P`inf then
      places:=InfinitePlaces(FFp); // infinite places of function field of curve mod p
      W:=Winf;
    else
      Px0:=Zeros(Fpx.1-Fp!x0)[1]; 
      places:=Decomposition(FFp,Px0); // places of function field of curve mod p lying over x0 mod p
      W:=W0;
    end if;

    bmodp:=[]; // elements of b mod p, where b is either b^0 or b^inf
    for i:=1 to d do
      f:=FFp!0;
      for j:=1 to d do
        f:=f+(Fpx!W[i,j])*FFp.1^(j-1);
      end for;
      bmodp[i]:=f;
    end for;

    done:=false;

    for i:=1 to #places do
      same:=true;
      for j:=1 to d do
        if Evaluate(bmodp[j],places[i]) ne Fp!b[j] then
          same:=false;
        end if;
      end for;    
      if same then
        place:=places[i];
        done:=true;
      end if;
    end for;

    if not done then
      error "Point does not lie on curve";
    end if;

    eP:=RamificationIndex(place);

    if eP eq 1 then
      index:=0;
    else
      done:=false;
      i:=1;
      while not done do
        ord:=Valuation(bmodp[i]-Evaluate(bmodp[i],place),place);
        if ord eq 1 then
          index:=i;
          done:=true;
        end if;
        i:=i+1;
      end while;
    end if;

    return eP,index,place,bmodp;
  end if;

end function;


hensel_lift:=function(fy,root,data);

  // Finds a root of the polynomial fy over Kp[[t]]
  // by Hensel lifting from an approximate root.
  //
  // Assumes that the starting criterion for Hensel's 
  // lemma is satisfied

  K:=data`K; ip:=data`ip; 
  Kpty:=Parent(fy);
  Kpt:=BaseRing(Kpty);
  Kp:=BaseRing(Kpt);
  tprec:=Precision(Kpt); // t-adic precision
  Kt:=PowerSeriesRing(K,tprec);
  Kty:=PolynomialRing(Kt);
  p:=Prime(Kp);
  pprec:=Precision(Kp);  // p-adic precision
  Op:=RingOfIntegers(Kp);
  Opt:=PowerSeriesRing(Op,tprec);  

  //fy:=Kty!fy;
  derfy:=Derivative(fy);

  if not Valuation(LeadingCoefficient(Evaluate(derfy,root))) eq 0 then
    error "In Hensel lift of power series, derivative has leading term divisible by p";
  end if;

  v1:=WeakValuation(Evaluate(fy,root));
  v2:=WeakValuation(Evaluate(derfy,root));

  if not v1 gt 2*v2 then
    error "Condition Hensel's Lemma not satisfied";
  end if;

  prec_seq:=[];
  k:=tprec;
  
  while k gt v1 do
    prec_seq:=Append(prec_seq,k);
    k:=Ceiling(k/2+v2);
  end while;
  prec_seq:=Reverse(prec_seq);

  for j:=1 to #prec_seq do
    root:=ChangePrecision(root,prec_seq[j]);
    root:=root-(Evaluate(fy,root))/(Evaluate(derfy,root));
    root:=Opt!root;
  end for;

  return root;

end function;


mod_p_prec:=function(fy);

  // Finds the t-adic precision necessary to separate the roots
  // of the polynomial fy over Kp[[t]] modulo p and start Hensel lift.
  //
  // Temporarily uses intrinsic Factorisation instead of 
  // intrinsic Roots because of multiple problems with Roots.
  
  Kpty:=Parent(fy);
  Kpt:=BaseRing(Kpty);
  tprec:=Precision(Kpt);
  Kp:=BaseRing(Kpt);
  Op:=RingOfIntegers(Kp);
  Fp,red:=ResidueClassField(Op);
  Fpt:=PowerSeriesRing(Fp,tprec);
  Fpty:=PolynomialRing(Fpt);

  f:=Fpty!0;
    C:=Coefficients(fy);
    for i:=1 to #C do
      D:=Coefficients(C[i]);
      for j:=1 to #D do
          f:=f+red(D[j])*Fpty.1^(i-1)*Fpt.1^(j-1);
      end for;
    end for;  

  fymodp:=f;
  derfymodp:=Derivative(fymodp);

  zeros:=[];
  fac:=Factorisation(fymodp); // can be slow...
  for i:=1 to #fac do
    if fac[i][2] gt 1 then
      error "t-adic precision not high enough";
    end if;
    factor:=fac[i][1];
    if Degree(factor) eq 1 and LeadingCoefficient(factor) eq 1 then
      zeros:=Append(zeros,-Coefficient(factor,0));
    end if;
  end for;

  modpprec:=1;
  for i:=1 to #zeros do
    done:=false;
    prec:=1;
    while not done do
      v1:=Valuation(Evaluate(fymodp,ChangePrecision(zeros[i],prec)));
      v2:=Valuation(Evaluate(derfymodp,ChangePrecision(zeros[i],prec)));
      if Minimum(prec,v1) gt 2*v2 then
        done:=true;
      end if;
      prec:=prec+1;
    end while;
    modpprec:=Maximum(modpprec,prec);
  end for;

  for i:=1 to #zeros do
    for j:=i+1 to #zeros do
      modpprec:=Maximum(modpprec,Valuation(zeros[i]-zeros[j]));
    end for;
  end for;
 
  return modpprec;
 
end function;


approx_root:=function(fy,y0,modpprec,expamodp,data)

  // Computes an approximation to t-adic precision modpprec of 
  // a root of the polynomial fy over Kp[[t]] which is congruent to:
  //
  // y0 modulo t
  // expamodp modulo p 
  //
  // This approximation is then used as root in hensel_lift.

  K:=data`K;
  Kpty:=Parent(fy);
  Kpt:=BaseRing(Kpty);
  tprec:=Precision(Kpt); // t-adic precision
  Kp:=BaseRing(Kpt);
  Op:=RingOfIntegers(Kp);
  Fp,red:=ResidueClassField(Op);
  p:=Characteristic(Fp);
  pprec:=Precision(Kp);  // p-adic precision
  Opt:=PowerSeriesRing(Op,tprec);
  Opz:=PolynomialRing(Op);

  Kt:=PowerSeriesRing(K,tprec);
  Kty:=PolynomialRing(Kt);
  Kpz:=PolynomialRing(Kp);
  Kpzt:=PowerSeriesRing(Kpz,tprec);
  
  roots:=[[*Kpt!y0,1*]];
  i:=1;
  while i le #roots do
    root:=roots[i][1];
    Nroot:=roots[i][2];
    if Nroot lt modpprec then
      roots:=Remove(roots,i);
      newroot:=root+Kpty.1*Kpt.1^Nroot;
      C:=Coefficients(fy);
      fynewroot:=Kpty!0;
      for j:=1 to #C do
        fynewroot:=fynewroot+(C[j])*newroot^(j-1);
      end for;
      fznewroot:=Kpzt!0;
      for j:=0 to Degree(fynewroot) do
        for k:=0 to tprec-1 do
          fznewroot:=fznewroot+Coefficient(Coefficient(fynewroot,j),k)*(Kpz.1)^j*(Kpzt.1)^k;
        end for;
      end for;
      fac:=Factorisation(Opz!Coefficient(fznewroot,WeakValuation(fznewroot)));
      for j:=1 to #fac do
        if (Degree(fac[j][1]) eq 1) and (Coefficient(fac[j][1],1) eq 1) then
          sol:=-Coefficient(fac[j][1],0); 
          if red(sol) eq Coefficient(expamodp,Nroot) then
            roots:=Insert(roots,i,[*Evaluate(newroot,sol),Nroot+1*]);
          end if;
          roots;
        end if;
      end for;
    else
      i:=i+1;
    end if;
  end while;
  roots;
  if #roots ne 1 then
    error "something is wrong, number of approximate roots is different from 1";
  end if;

  root:=roots[1][1];
  root:=Opt!root;

  v1:=Valuation(Evaluate(fy,root));
  v2:=Valuation(Evaluate(Derivative(fy),root));

  if v1 le 2*v2 then
    error "something is wrong, approximate root not good enough for Hensel lift";
  end if;

  return root;

end function;


mod_p_expansion:=function(f,place,tmodp,modpprec);

  // Finds the power series expansion of f in the function field
  // modulo p at place with respect to local parameter tmodp to
  // absolute precision modpprec.

  FFp:=Parent(f);
  Fpx:=BaseRing(FFp);
  Fp:=BaseRing(Fpx);
  Fpt:=PowerSeriesRing(Fp,modpprec);

  expamodp:=Fpt!0;
  for i:=0 to modpprec-1 do
    val:=Evaluate(f,place);
    expamodp:=expamodp+val*Fpt.1^i;
    f:=(f-val)/tmodp;
  end for;
  
  return expamodp;
  
end function;


local_coord:=function(P,prec,data);

  // Computes powerseries expansions of x and
  // the b^0_i or b^infty_i (depending on whether
  // P is infinite or not) in terms of the local
  // coordinate computed by local_data.

  if assigned P`xt and Precision(Parent(P`xt)) ge prec then
    xt:=P`xt;
    bt:=P`bt;
    index:=P`index;
    return xt,bt,index;
  end if;

  if is_bad(P,data) and not is_very_bad(P,data) then
    error "Cannot compute local parameter at a bad point which is not very bad";
  end if;

  x0:=P`x; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; d:=Degree(Q); b:=P`b;
  K:=data`K; Kxy:=data`Kxy; Kp:=data`Kp; ip:=data`ip;
  //Kp:=Parent(x0);
  Kpt<t>:=PowerSeriesRing(Kp,prec); Kpty:=PolynomialRing(Kpt);
  KptF:=FieldOfFractions(Kpt);
  Kt:=RationalFunctionField(K); Kty:=PolynomialRing(Kt);
  Kpx:=PolynomialRing(Kp);
  Fq:=ResidueClassField(RingOfIntegers(Kp));

  f:=Kty!0;
  QK:=Zaxy_to_Kxy(Q,Kxy);
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(QK,i)) do
      f:=f+Coefficient(Coefficient(QK,i),j)*Kty.1^i*Kt.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve
  
  W0Kp:=matrix_push_to_Kp(W0,d,Kpx);
  if not is_bad(P,data) then // finite good point

    xt:=t+x0;

    W0invx0:=Transpose(Evaluate(W0Kp^(-1),x0));
    ypowers:=Vector(b)*W0invx0;
    y0:=ypowers[2];

    C:=[Kx_to_Kpt(c,ip,Kpt) : c in Coefficients(QK)];
    D:=[];
    for i:=1 to #C do
      D[i]:=Evaluate(C[i],xt); 
    end for;
    fy:=Kpty!D;
    derfy:=Derivative(fy);

    yt:=hensel_lift(fy,Kpt!y0,data);

    ypowerst:=[];
    ypowerst[1]:=KptF!1;
    ypowerst[2]:=yt;
    for i:=3 to d do
      ypowerst[i]:=ypowerst[i-1]*yt;
    end for;
    
    bt:=Eltseq(Vector(ypowerst)*ChangeRing(Transpose(Evaluate(W0Kp,xt)),KptF));
    
    btnew:=[];
    for i:=1 to d do
      btnew[i]:=Kpt!bt[i];
    end for;
    bt:=btnew;

    index:=0;

  elif P`inf then // infinite point

    eP,index,place,bmodp:=local_data(P,data);
    FFp:=Parent(bmodp[1]);
    Fpx:=BaseRing(FFp);

    bfun:=[];
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+Winf[i,j]*FF.1^(j-1);
      end for;
      bfun[i]:=bi;
    end for;
    
    if eP eq 1 then // P is an infinite point that is not ramified
      
      xt:=t+x0;
      bt:=[];

      for i:=1 to d do

        if assigned data`minpolys and data`minpolys[2][1,i+1] ne 0 then
          poly:=data`minpolys[2][1,i+1]; 
        else 
          poly:=minpoly(FF!(1/Kt.1),bfun[i]);
        end if;

        C:=[Kx_to_Kpt(c,ip,Kpt) : c in Coefficients(poly)];
        D:=[];
        for j:=1 to #C do
          D[j]:=Evaluate(C[j],xt); 
        end for;
        fy:=Kpty!D;
        derfy:=Derivative(fy);

        modpprec:=mod_p_prec(fy);

        if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
          approxroot:=P`bt[i];
        else
          tmodp:=1/Fpx.1-Fq!x0;
          expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
          approxroot:=approx_root(fy,b[i],modpprec,expamodp);
        end if;

        bti:=hensel_lift(fy,approxroot,data);
        bt[i]:=bti;

      end for;

    else // P is an infinite point that is ramified

      if assigned data`minpolys and data`minpolys[2][index+1,1] ne 0 then
        poly:=data`minpolys[2][index+1,1];
      else
        poly:=minpoly(bfun[index],FF!1/(Kt.1));
      end if;

      C:=[Kx_to_Kpt(c,ip,Kpt) : c in Coefficients(poly)];
      D:=[];
      for j:=1 to #C do
        D[j]:=Evaluate(C[j],t+b[index]); 
      end for;
      fy:=Kpty!D;
      derfy:=Derivative(fy);

      modpprec:=mod_p_prec(fy);

      if assigned P`xt and Precision(Parent(P`xt)) ge modpprec then
        approxroot:=P`xt;
      else
        tmodp:=bmodp[index]-Fq!b[index];
        expamodp:=mod_p_expansion(FFp!1/Fpx.1,place,tmodp,modpprec);
        approxroot:=approx_root(fy,x0,modpprec,expamodp);
      end if;

      xt:=hensel_lift(fy,approxroot,data);

      bt:=[];
      for i:=1 to d do 
      
        if i eq index then
          bt[i]:=t+b[i];
        else
          
          if assigned data`minpolys and data`minpolys[2][index+1,i+1] ne 0 then
            poly:=data`minpolys[2][index+1,i+1];
          else
            poly:=minpoly(bfun[index],bfun[i]);
          end if;

          C:=[Kx_to_Kpt(c,ip,Kpt) : c in Coefficients(poly)];
          D:=[];
          for j:=1 to #C do
            D[j]:=Evaluate(C[j],t+b[index]); 
          end for;

          fy:=Kpty!D;
          derfy:=Derivative(fy);

          modpprec:=mod_p_prec(fy);

          if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
            approxroot:=P`bt[i];
          else
            tmodp:=bmodp[index]-Fq!b[index];
            expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
            approxroot:=approx_root(fy,b[i],modpprec,expamodp);
          end if;

          bti:=hensel_lift(fy,approxroot,data);
          bt[i]:=bti;

        end if;
 
      end for;

    end if;

  else // finite bad point

    eP,index,place,bmodp:=local_data(P,data);
    FFp:=Parent(bmodp[1]);
    Fpx:=BaseRing(FFp);

    bfun:=[];
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+W0[i,j]*FF.1^(j-1);
      end for;
      bfun[i]:=bi;
    end for;

    if eP eq 1 then // P is a finite point that is not ramified

      xt:=t+x0;
      bt:=[];
      for i:=1 to d do
        
        if assigned data`minpolys and data`minpolys[1][1,i+1] ne 0 then
          poly:=data`minpolys[1][1,i+1];
        else
          poly:=minpoly(FF!Kt.1,bfun[i]);
        end if;

        C:=[Kx_to_Kpt(c,ip,Kpt) : c in Coefficients(poly)];
        D:=[];
        for j:=1 to #C do
          D[j]:=Evaluate(C[j],xt); 
        end for;
        fy:=Kpty!D;
        derfy:=Derivative(fy);

        modpprec:=mod_p_prec(fy);

        if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
          approxroot:=P`bt[i];
        else
          tmodp:=Fpx.1-Fq!x0;
          expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
          approxroot:=approx_root(fy,b[i],modpprec,expamodp);
        end if;

        bti:=hensel_lift(fy,approxroot,data);
        bt[i]:=bti;

      end for;

    else // P is a finite point that ramifies

      if assigned data`minpolys and data`minpolys[1][index+1,1] ne 0 then
        poly:=data`minpolys[1][index+1,1];
      else
        poly:=minpoly(bfun[index],FF!Kt.1);
      end if;

      C:=[Kx_to_Kpt(c,ip,Kpt) : c in Coefficients(poly)];
      D:=[];
      for j:=1 to #C do
        D[j]:=Evaluate(C[j],t+b[index]); 
      end for;
      fy:=Kpty!D;
      derfy:=Derivative(fy);

      modpprec:=mod_p_prec(fy);

      if assigned P`xt and Precision(Parent(P`xt)) ge modpprec then
        approxroot:=P`xt;
      else
        tmodp:=bmodp[index]-Fq!b[index];
        expamodp:=mod_p_expansion(FFp!Fpx.1,place,tmodp,modpprec);
        approxroot:=approx_root(fy,x0,modpprec,expamodp,data);
      end if;

      xt:=hensel_lift(fy,approxroot,data);  

      bt:=[];
      for i:=1 to d do 
      
        if i eq index then
          bt[i]:=t+b[i];
        else
          
          if assigned data`minpolys and data`minpolys[1][index+1,i+1] ne 0 then
            poly:=data`minpolys[1][index+1,i+1];
          else
            poly:=minpoly(bfun[index],bfun[i]);
          end if;

          C:=[Kx_to_Kpt(c,ip,Kpt) : c in Coefficients(poly)];
          D:=[];
          for j:=1 to #C do
            D[j]:=Evaluate(C[j],t+b[index]);
          end for;

          fy:=Kpty!D;

          derfy:=Derivative(fy);

          modpprec:=mod_p_prec(fy);

          if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
            approxroot:=P`bt[i];
          else
            tmodp:=bmodp[index]-Fq!b[index];
            expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
            approxroot:=approx_root(fy,b[i],modpprec,expamodp);
          end if;

          bti:=hensel_lift(fy,approxroot,data);
          bt[i]:=bti;

        end if;
 
      end for;

    end if;

  end if;

  return xt,bt,index;

end function;


tiny_integral_prec:=function(prec,e,maxpoleorder,maxdegree,mindegree,val,data);

  // Determines the p-adic precision to which tiny_integrals_on_basis is correct.

  N:=data`N; p:=data`p;

  // Precision loss from terms of positive order we do consider:

  m1:=N*e-val;
  for i:=1 to maxdegree do
    m1:=Minimum(m1,N*e+i-e*Floor(Log(p,i+1)));
  end for;  

  // Precision loss from terms we omit:

  m2:=mindegree+2-e*Floor(Log(p,mindegree+2));
  for i:=mindegree+2 to Ceiling(e/Log(p)) do
    m2:=Minimum(m2,i+1-e*Floor(Log(p,i+1)));
  end for;

  // Precision loss from terms of negative order

  m3:=N*e-val;
  if maxpoleorder ge 2 then
    m3:=N*e-val-maxpoleorder*val-e*Floor(Log(p,maxpoleorder-1));
  end if;

  m:=Minimum([m1,m2,m3]);

  return m/e;

end function;


find_bad_point_in_disk:=function(P,data);

  // Find the very bad point in the residue disk of a bad point P.

  x0:=P`x; b:=P`b; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; r:=data`r;
  d:=Degree(Q); K:=Parent(x0); Ky:=PolynomialRing(K); 

  if not is_bad(P,data) then
    error "Residue disk does not contain a bad point";
  end if;

  if P`inf then
    x0:=K!0;
  else
    rQp:=Ky!Coefficients(r);
    x0:=HenselLift(rQp,x0);
  end if;

  Qt:=RationalFunctionField(RationalField()); Qty:=PolynomialRing(Qt);

  f:=Qty!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+Coefficient(Coefficient(Q,i),j)*Qty.1^i*Qt.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  eP,index:=local_data(P,data);

  if P`inf then
    W:=Winf;
  else
    W:=W0;
  end if;

  bfun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+W[i,j]*FF.1^(j-1);
    end for;
    bfun:=Append(bfun,bi);
  end for;

  if index eq 0 then
    if P`inf then
      xfun:=FF!(1/Qt.1);
    else
      xfun:=FF!(Qt.1);
    end if;

    for i:=1 to d do
      poly:=minpoly(xfun,bfun[i]);
      C:=Coefficients(poly);
      D:=[];
      for i:=1 to #C do
        D[i]:=Evaluate(C[i],x0); 
      end for;
      fy:=Ky!D;
      fac:=Factorisation(fy);
      done:=false;
      j:=1;
      while not done and j le #fac do
        if Degree(fac[j][1]) eq 1 and Valuation(-Coefficient(fac[j][1],0)-b[i]) gt 0 then
          done:=true;
          b[i]:=-Coefficient(fac[j][1],0);
        end if;
        j:=j+1;
      end while;
    end for;
  else
   bindex:=bfun[index];
   if P`inf then
      xfun:=FF!(1/Qt.1);
    else
      xfun:=FF!(Qt.1);
    end if;
    poly:=minpoly(xfun,bindex);
    C:=Coefficients(poly);
    D:=[];
    for i:=1 to #C do
      D[i]:=Evaluate(C[i],x0); 
    end for;
    fy:=Ky!D;
    fac:=Factorisation(fy);
    done:=false;
    j:=1;
    while not done and j le #fac do
      if Degree(fac[j][1]) eq 1 and Valuation(-Coefficient(fac[j][1],0)-b[index]) gt 0 then
        done:=true;
        b[index]:=-Coefficient(fac[j][1],0);
      end if;
      j:=j+1;
    end while;
    for i:=1 to d do
      if i ne index then
        poly:=minpoly(bindex,bfun[i]);
        C:=Coefficients(poly);
        D:=[];
        for i:=1 to #C do
          D[i]:=Evaluate(C[i],b[index]); 
        end for;
        fy:=Ky!D;
        fac:=Factorisation(fy); // Roots has some problems that Factorisation does not
        done:=false;
        j:=1;
        while not done and j le #fac do
          if Degree(fac[j][1]) eq 1 and Valuation(-Coefficient(fac[j][1],0)-b[i]) gt 0 then
            done:=true;
            b[i]:=-Coefficient(fac[j][1],0);
          end if;
          j:=j+1;
        end while;
      end if;
    end for; 
  end if;

  Pbad:=set_bad_point(x0,b,P`inf,data);

  return Pbad;

end function;


tadicprec:=function(data,e);

  // Compute the necessary t-adic precision to compute tiny integrals

  p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf;
  W:=Winf*W0^(-1);

  prec:=1;
  while Floor(prec/e)+1-Floor(Log(p,prec+1)) lt N do
    prec:=prec+1;
  end while;
  prec:=Maximum([prec,100]); // 100 arbitrary, avoid problems with small precisions 

  return prec;

end function;


tiny_integrals_on_basis:=function(P1,P2,data:prec:=0,P:=0);

  // Compute tiny integrals of basis elements from P1 to P2.
  // If P1 is not defined over Qp (but a totally ramified 
  // extension) then a point P defined over Qp in the same
  // residue disk as P1 has to be specified.

  x1:=P1`x; x2:=P2`x; b1:=P1`b; b2:=P2`b; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; Kx:=data`Kx; r:=data`r; basis:=data`basis; N:=data`N;
  d:=Degree(Q); W:=Winf*W0^(-1); Kp:=data`Kp;  rK:=data`rK; K:=data`K;

  if not lie_in_same_disk(P1,P2,data) then
    error "the points do not lie in the same residue disk";
  end if;

  //if ((x1 eq x2) and (b1 eq b2)) then 
  //  return RSpace(Kp,#basis)!0, N*Degree(Kp);
  //end if;

  if (Valuation(x1-x2)/Valuation(Parent(x1-x2)!p) ge N) and (Minimum([Valuation(b1[i]-b2[i])/Valuation(Parent(b1[i]-b2[i])!p):i in [1..d]]) ge N) then
    return RSpace(Kp,#basis)!0, N*Degree(Kp);
  end if; 


  if AbsoluteRamificationDegree(Kp) gt 1 then // P1 needs to be defined over Qq
    tinyPtoP2,NtinyPtoP2:=$$(P,P2,data);
    tinyPtoP1,NtinyPtoP1:=$$(P,P1,data);
    return tinyPtoP2-tinyPtoP1,Minimum(NtinyPtoP2,NtinyPtoP1);
  end if;

  if not Type(P) eq Rec then
    P:=P1;
  end if;

  if is_bad(P,data) and not is_very_bad(P,data) then // on a bad disk P1 needs to be very bad
    P:=find_bad_point_in_disk(P,data);
    tinyPtoP2,NtinyPtoP2:=$$(P,P2,data);
    tinyPtoP1,NtinyPtoP1:=$$(P,P1,data);
    return tinyPtoP2-tinyPtoP1,Minimum(NtinyPtoP2,NtinyPtoP1);
  end if;

  e:=AbsoluteRamificationDegree(Parent(x2));

  if prec eq 0 then // no t-adic precision specified
    prec:=tadicprec(data,e);
  end if;

  Kpt:=LaurentSeriesRing(Kp,prec);
  Op:=RingOfIntegers(Kp);
  Opt:=LaurentSeriesRing(Op,prec);

  xt,bt,index:=local_coord(P1,prec,data);

  Kt<t>:=LaurentSeriesRing(K,prec);
  /*xt:=Qt!xt;
  btnew:=[Qt|];
  for i:=1 to d do
    btnew[i]:=Qt!bt[i];
  end for;
  bt:=Vector(btnew);*/
  
  // WE ONLY CONSIDER FINITE GOOD POINTS

  if P1`inf then
    error "We only consider finite good points";
/*    xt:=1/xt;
    xt:=Qt!Kpt!xt; 
    Winv:=W0*Winf^(-1);          
    bt:=bt*Transpose(Evaluate(Winv,xt));
    for i:=1 to d do
      bt[i]:=Qt!(Kpt!bt[i]);
    end for; */
  end if;
  rKpt:=Kx_to_Kpt(rK,data`ip,Kpt);
  if P1`inf or not is_bad(P1,data) then 
    denom:=(1/Evaluate(rKpt,xt));
  else
    error "We only consider finite good points";
    /*
    Qp:=pAdicField(p,N);
    Qpx:=PolynomialRing(Qp);
    rQp:=Qpx!r;
    zero:=HenselLift(rQp,x1);
    sQp:=rQp div (Qpx.1-zero);
    denom:=Qt!Kt!((Qt!OKt!(xt-Coefficient(xt,0)))^(-1)*(Qt!Kt!(1/Evaluate(sQp,xt))));*/
  end if;

  derxt:=Derivative(xt); 
  diffs:=[];
  basisKpt:=[Vector([Kx_to_Kpt(Zax_to_Kx(basis[i][j],Kx),data`ip,Kpt): j in [1..d]]) : i in [1..#basis]];
  for i:=1 to #basis do
    basisxt:=Evaluate(basisKpt[i],xt);
    /*for j:=1 to d do
      basisxt[1][j]:=Qt!Kt!basisxt[1][j];
    end for;*/
    diffs[i]:=Kpt!(InnerProduct(Vector(basisxt*derxt*denom),Vector(bt)));
    //diffs[i]:=Qt!Kt!diffs[i];
    if Coefficient(diffs[i],-1) ne 0 then
      diffs[i]:=diffs[i]-Coefficient(diffs[i],-1)*Kpt.1^(-1); // temporary, deal with logs later, important for double integrals
    end if;
  end for;
  A:=[WeakValuation(diffs[i]): i in [1..#basis] | not IsWeaklyZero(diffs[i])];
  if IsEmpty(A) then
    maxpoleorder:=-N;
  else
    maxpoleorder:=-(Minimum([WeakValuation(diffs[i]): i in [1..#basis] | not IsWeaklyZero(diffs[i])]));
  end if;
  B:=[Degree(diffs[i]): i in [1..#basis] | not IsWeaklyZero(diffs[i])];
  if IsEmpty(B) then
    maxdegree:=0;
    mindegree:=1000;
  else
    maxdegree:=Maximum([Degree(diffs[i]): i in [1..#basis] | not IsWeaklyZero(diffs[i])]);
    mindegree:=Minimum([Degree(diffs[i]): i in [1..#basis] | not IsWeaklyZero(diffs[i])]);
  end if;

  indefints:=[];
  for i:=1 to #basis do
    indefints := Append(indefints, Integral(diffs[i]));
  end for;

  tinyP1toP2:=[];
  for i:=1 to #basis do
    if index eq 0 then // x-x(P1) is the local coordinate
      tinyP1toP2[i]:=Evaluate(indefints[i],x2-x1);
      val:=Valuation(x2-x1);
    else // b[index]-b[index](P1) is the local coordinate
      tinyP1toP2[i]:=Evaluate(indefints[i],b2[index]-b1[index]);
      val:=Valuation(b2[index]-b1[index]);
    end if;
  end for;

  NtinyP1toP2:=tiny_integral_prec(prec,e,maxpoleorder,maxdegree,mindegree,val,data);

  return Vector(tinyP1toP2),NtinyP1toP2;

end function;


tiny_integrals_on_basis_to_z:=function(P,data:prec:=0);

  // Compute tiny integrals of basis elements from P to an
  // arbitrary point in its residue disk as a power series
  // in the local parameter there. The series expansions xt
  // and bt of the coordinates on the curve in terms of this 
  // local parameter are also returned.

  x0:=P`x; b:=P`b; Q:=data`Q; p:=data`p; N:=data`N; basis:=data`basis; r:=data`r; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q); lc_r:=LeadingCoefficient(r); W:=Winf*W0^(-1); K:=Parent(x0);

  if is_bad(P,data) and not is_very_bad(P,data) then // on a bad disk P needs to be very bad
    P1:=find_bad_point_in_disk(P,data);  
  else
    P1:=P;
  end if;
  x1:=P1`x;

  IPP1,NIPP1:=tiny_integrals_on_basis(P,P1,data:prec:=prec);

  if prec eq 0 then // no t-adic precision specified
    prec:=tadicprec(data,1);
  end if;

  Kt<t>:=LaurentSeriesRing(K,prec);
  OK:=RingOfIntegers(K);
  OKt:=LaurentSeriesRing(OK,prec);

  xt,bt,index:=local_coord(P1,prec,data);

  xtold:=xt;
  btold:=bt;

  Qt<t>:=LaurentSeriesRing(RationalField(),prec);
  xt:=Qt!xt;
  btnew:=[Qt|];
  for i:=1 to d do
    btnew[i]:=Qt!bt[i];
  end for;
  bt:=Vector(btnew);

  if P1`inf then
    xt:=1/xt;
    xt:=Qt!Kt!xt; 
    Winv:=W0*Winf^(-1);          
    bt:=bt*Transpose(Evaluate(Winv,xt));
    for i:=1 to d do
      bt[i]:=Qt!(Kt!bt[i]);
    end for; 
  end if;

  if P1`inf or not is_bad(P1,data) then 
    denom:=Qt!Kt!(1/Evaluate(r,xt));
  else
    Qp:=pAdicField(p,N);
    Qpx:=PolynomialRing(Qp);
    rQp:=Qpx!r;
    zero:=HenselLift(rQp,x1);
    sQp:=rQp div (Qpx.1-zero);
    denom:=Qt!Kt!((Qt!OKt!(xt-Coefficient(xt,0)))^(-1)*(Qt!Kt!(1/Evaluate(sQp,xt))));
  end if;

  derxt:=Qt!Kt!Derivative(xt); 
  diffs:=[];
  for i:=1 to #basis do
    basisxt:=Evaluate(basis[i],xt);
    for j:=1 to d do
      basisxt[1][j]:=Qt!Kt!basisxt[1][j];
    end for;
    diffs[i]:=InnerProduct(Vector(basisxt*derxt*lc_r*denom),bt);
    diffs[i]:=Qt!Kt!diffs[i];
    if Coefficient(diffs[i],-1) ne 0 then
      diffs[i]:=diffs[i]-Coefficient(diffs[i],-1)*t^(-1); // temporary, TODO deal with logs later, important for double integrals
    end if;
  end for;

  indefints:=[];
  for i:=1 to #basis do
    indefints := Append(indefints, Integral(diffs[i]));
  end for;

  xt:=xtold;
  bt:=Vector(btold);

  return Vector(indefints)+IPP1,xt,bt,NIPP1;

end function;


pow:=function(x,k);

  if k eq 0 then
    return Parent(x)!1;
  else
    return x^k;
  end if;

end function;


evalf0:=function(f0,P,data);

  // Evaluate vector of functions f0 at P.
 
  x0:=P`x; b:=P`b; Q:=data`Q; rK:=data`rK; W0:=data`W0; Winf:=data`Winf; N:=data`N; Nmax:=data`Nmax; p:=data`p;
  d:=Degree(Q); Kp:=data`Kp;
  Kpx:=PolynomialRing(Kp);
  Kpt:=LaurentSeriesRing(Kp);
  rKp:=Kx_to_Kpt(rK,data`ip,Kpt);
  valf0:=0;

  if P`inf then 
    Winv:=W0*Winf^(-1); 
    WinvKp:=matrix_push_to_Kp(Winv,d,Kpx);
    b:=Vector(b)*Transpose(Evaluate(Evaluate(Winv,1/Kpt.1),x0)); // values of the b_i^0 at P
    
    z0:=Evaluate(rKp,1/x0);
    invz0:=1/z0;
    invz0pow:=[Kp!1];
    for i:=1 to p*(Nmax-1) do
      invz0pow[i+1]:=invz0pow[i]*invz0;
    end for;
    
    invx0:=1/x0;
    invx0pow:=[Kp!1];
    for i:=1 to Degree(rK)-1 do
      invx0pow[i+1]:=invx0pow[i]*invx0;
    end for;

    f0P:=Kp!0;
    for i:=1 to d do
      f0i:=f0[i];
      C:=Coefficients(f0i);
      val:=Valuation(f0i);
      for j:=1 to #C do
        D:=Coefficients(C[j]);
        for k:=1 to #D do
          f0P:=f0P+(Kp!D[k])*invx0pow[k]*invz0pow[2-j-val]*b[i];
          valf0:=Minimum(valf0,Valuation(K!D[k]));
        end for;
      end for;
    end for;
    Nf0P:=N*Degree(Kp)+(ord_inf_mat(Winv)+1)*Valuation(x0)+valf0;

  else
    
    z0:=Evaluate(rKp,x0);  
    invz0:=1/z0;
    invz0pow:=[Kp!1];
    for i:=1 to p*(Nmax-1) do
      invz0pow[i+1]:=invz0pow[i]*invz0;
    end for;

    x0pow:=[Kp!1];
    for i:=1 to Degree(rKp)-1 do
      x0pow[i+1]:=x0pow[i]*x0;
    end for;  
 
    f0P:=Kp!0;
    for i:=1 to d do
      f0i:=f0[i];
      C:=Coefficients(f0i);
      val:=Valuation(f0i);
      for j:=1 to #C do
        D:=Coefficients(C[j]);
        for k:=1 to #D do
          f0P:=f0P+(Kp!D[k])*x0pow[k]*invz0pow[2-j-val]*b[i];
          valf0:=Minimum(valf0,Valuation(Kp!D[k]));
        end for;
      end for;
    end for;
    Nf0P:=N*Degree(Kp)-p*(Nmax-1)*Valuation(z0)+valf0; // TODO this is error of terms we did consider, take error of terms we ignored into account as well
  end if;

  return f0P,Nf0P/Degree(Kp);

end function;


evalfinf:=function(finf,P,data);

  // Evaluate vector of functions finf at P.

  x0:=P`x; b:=P`b; Q:=data`Q; W0:=data`W0; Winf:=data`Winf; N:=data`N; p:=data`p;
  d:=Degree(Q); K:=Parent(x0); 

  W:=Winf*W0^(-1); 

  valfinf:=0;

  if P`inf then
    finfP:=K!0;
    for i:=1 to d do
      finfi:=finf[i];
      C:=Coefficients(finfi);
      val:=Valuation(finfi);
      for j:=1 to #C do
        finfP:=finfP+(K!C[j])*pow(1/x0,val+j-1)*b[i];
        valfinf:=Minimum(valfinf,Valuation(K!C[j]));
      end for;
    end for;
    NfinfP:=N*Degree(K)+p*(ord_0_mat(W)+1)*Valuation(x0)+valfinf;
  else 
    finf:=finf*ChangeRing(W,BaseRing(finf));
    finfP:=K!0;
    for i:=1 to d do
      finfi:=finf[i];
      C:=Coefficients(finfi);
      val:=Valuation(finfi);
      for j:=1 to #C do
        finfP:=finfP+(K!C[j])*pow(x0,val+j-1)*b[i];
        valfinf:=Minimum(valfinf,Valuation(K!C[j]));
      end for;
    end for;
    NfinfP:=N*Degree(K)+valfinf;
  end if;

  return finfP, NfinfP/Degree(K);

end function;


evalfend:=function(fend,P,data);

  // Evaluate vector of functions fend at P.

  x0:=P`x; b:=P`b; Q:=data`Q; W0:=data`W0; Winf:=data`Winf; N:=data`N;
  d:=Degree(Q);
  K:=Parent(x0);

  valfend:=0;

  if P`inf then
    Winv:=W0*Winf^(-1);
    Qt:=BaseRing(Winv);
    b:=Vector(b)*Transpose(Evaluate(Evaluate(Winv,1/Qt.1),x0)); // values of the b_i^0 at P
    fendP:=K!0;
    for i:=1 to d do
      fendi:=fend[i];
      C:=Coefficients(fendi);
      for j:=1 to #C do
        fendP:=fendP+(K!C[j])*pow(1/x0,j-1)*b[i];
        valfend:=Minimum(valfend,Valuation(K!C[j]));
      end for;
    end for;
    NfendP:=N*Degree(K)+(ord_0_mat(Winf)+1)*Valuation(x0)+valfend;
  else
    fendP:=K!0;
    for i:=1 to d do
      fendi:=fend[i];
      C:=Coefficients(fendi);
      for j:=1 to #C do
        fendP:=fendP+(K!C[j])*pow(x0,j-1)*b[i];
        valfend:=Minimum(valfend,Valuation(K!C[j]));
      end for;
    end for;
    NfendP:=N*Degree(K)+valfend;
  end if;

  return fendP, NfendP/Degree(K);

end function;


round_to_Qp:=function(L)

  // Rounds a vector over a totally ramified extension of Qp to one over Qp.

  K:=CoefficientRing(L);
  deg:=Degree(K);
  e:=Precision(K);

  l:=[];
  for i:=1 to #Eltseq(L) do
    l[i]:=Eltseq(L[i])[1];  
    e:=Minimum(e,Valuation(L[i]-l[i]));
  end for;

  return Vector(l),e/deg;

end function;


coleman_integrals_on_basis:=function(P1,P2,data:e:=1)

  // Integrals of basis elements from P1 to P2. 

  F:=data`F; Q:=data`Q; basis:=data`basis; x1:=P1`x; f0list:=data`f0list; finflist:=data`finflist; fendlist:=data`fendlist; p:=data`p; N:=data`N; delta:=data`delta;
  d:=Degree(Q); K:=Parent(x1); 

  // First make sure that if P1 or P2 is bad, then it is very bad

  
  if is_bad(P1,data) and not is_very_bad(P1,data) then
    S1:=find_bad_point_in_disk(P1,data);
    _,index:=local_data(S1,data);
    data:=update_minpolys(data,S1`inf,index);
    xt,bt,index:=local_coord(S1,tadicprec(data,e),data);
    S1`xt:=xt;
    S1`bt:=bt;
    S1`index:=index;
    IS1P1,NIS1P1:=tiny_integrals_on_basis(S1,P1,data:prec:=tadicprec(data,e));
    IS1P2,NIS1P2:=$$(S1,P2,data:e:=e);
    IP1P2:=IS1P2-IS1P1;
    NIP1P2:=Ceiling(Minimum([NIS1P1,NIS1P2]));
    return IP1P2,NIP1P2;
  end if;

  if is_bad(P2,data) and not is_very_bad(P2,data) then
    S2:=find_bad_point_in_disk(P2,data);
    _,index:=local_data(S2,data);
    data:=update_minpolys(data,S2`inf,index);
    xt,bt,index:=local_coord(S2,tadicprec(data,e),data);
    S2`xt:=xt;
    S2`bt:=bt;
    S2`index:=index;
    IP1S2,NIP1S2:=$$(P1,S2,data:e:=e);
    IP2S2,NIP2S2:=tiny_integrals_on_basis(P2,S2,data:prec:=tadicprec(data,e));
    IP1P2:=IP1S2-IP2S2;
    NIP1P2:=Ceiling(Minimum([NIP1S2,NIP2S2]));
    return IP1P2,NIP1P2;
  end if;

  // If P1,P2 is bad (hence very bad), use a near boundary point.

  _,index:=local_data(P1,data);
  data:=update_minpolys(data,P1`inf,index);
  _,index:=local_data(P2,data);
  data:=update_minpolys(data,P2`inf,index);

  if is_bad(P1,data) then
    xt,bt,index:=local_coord(P1,tadicprec(data,e),data);
    P1`xt:=xt;       
    P1`bt:=bt;       
    P1`index:=index; 
    Qp:=Parent(P1`x);
    Qpa<a>:=PolynomialRing(Qp);
    K<a>:=TotallyRamifiedExtension(Qp,a^e-p);
    format:=recformat<x,b,inf,xt,bt,index>;
    S1:=rec<format|>;                                                    
    S1`inf:=P1`inf;
    S1`x:=Evaluate(xt,a);
    S1`b:=[Evaluate(bt[i],a):i in [1..d]];
  else
    xt,bt,index:=local_coord(P1,tadicprec(data,1),data);
    P1`xt:=xt;       
    P1`bt:=bt;       
    P1`index:=index; 
    S1:=P1;
  end if;

  if is_bad(P2,data) then
    xt,bt,index:=local_coord(P2,tadicprec(data,e),data);
    P2`xt:=xt;       
    P2`bt:=bt;       
    P2`index:=index; 
    if not is_bad(P1,data) then
      Qp:=Parent(P2`x);
      Qpa<a>:=PolynomialRing(Qp);
      K<a>:=TotallyRamifiedExtension(Qp,a^e-p);
    end if;
    format:=recformat<x,b,inf,xt,bt,index>;
    S2:=rec<format|>;                                                    
    S2`inf:=P2`inf;
    S2`x:=Evaluate(xt,a);
    S2`b:=[Evaluate(bt[i],a):i in [1..d]];
  else
    xt,bt,index:=local_coord(P2,tadicprec(data,1),data);
    P2`xt:=xt;       
    P2`bt:=bt;       
    P2`index:=index; 
    S2:=P2;
  end if;

  // Split up the integral and compute the tiny ones.

  tinyP1toS1,NP1toS1:=tiny_integrals_on_basis(P1,S1,data);
  tinyP2toS2,NP2toS2:=tiny_integrals_on_basis(P2,S2,data);

  FS1:=frobenius_pt(S1,data);
  FS2:=frobenius_pt(S2,data);

  tinyS1toFS1,NS1toFS1:=tiny_integrals_on_basis(S1,FS1,data:P:=P1); 
  tinyS2toFS2,NFS2toS2:=tiny_integrals_on_basis(S2,FS2,data:P:=P2); 

  NIP1P2:=Minimum([NP1toS1,NP2toS2,NS1toFS1,NFS2toS2]);  

  // Evaluate all functions.

  I:=[];
  for i:=1 to #basis do
    f0iS1,Nf0iS1:=evalf0(f0list[i],S1,data);
    f0iS2,Nf0iS2:=evalf0(f0list[i],S2,data);
    finfiS1,NfinfiS1:=evalfinf(finflist[i],S1,data);
    finfiS2,NfinfiS2:=evalfinf(finflist[i],S2,data);
    fendiS1,NfendiS1:=evalfend(fendlist[i],S1,data);
    fendiS2,NfendiS2:=evalfend(fendlist[i],S2,data);
    NIP1P2:=Minimum([NIP1P2,Nf0iS1,Nf0iS2,NfinfiS1,NfinfiS2,NfendiS1,NfendiS2]);
    I[i]:=(K!f0iS1)-(K!f0iS2)+(K!finfiS1)-(K!finfiS2)+(K!fendiS1)-(K!fendiS2)-(K!tinyS1toFS1[i])+(K!tinyS2toFS2[i]);
  end for; 

  valIP1P2:=Minimum([Valuation(I[i])/Valuation(K!p):i in [1..#basis]]);

  mat:=(F-IdentityMatrix(RationalField(),#basis));
  valdet:=Valuation(Determinant(mat),p);
  mat:=mat^-1;
  Nmat:=N-valdet-delta;
  valmat:=Minimum([Valuation(e,p):e in Eltseq(mat)]);

  NIP1P2:=Minimum([NIP1P2+valmat,Nmat+valIP1P2]);                            
  
  IS1S2:=Vector(I)*Transpose(ChangeRing(mat,K));    // Solve the linear system.
  IP1P2:=IS1S2+ChangeRing(tinyP1toS1,K)-ChangeRing(tinyP2toS2,K);
  IP1P2,Nround:=round_to_Qp(IP1P2);

  assert Nround ge NIP1P2;                          // Check that rounding error is within error bound.
  
  NIP1P2:=Ceiling(NIP1P2);

  for i:=1 to #basis do
    IP1P2[i]:=IP1P2[i]+O(Parent(IP1P2[i])!p^(NIP1P2));
  end for;

  return IP1P2,NIP1P2;
end function;




coleman_integrals_on_basis_between_good_points:=function(P1,P2,data:e:=1)

  // Integrals of basis elements from P1 to P2. We ASSUME that P1 and P2 are good. Else, behavior is undefined. 

  F:=data`F; Q:=data`Q; basis:=data`basis; x1:=P1`x; f0list:=data`f0list; finflist:=data`finflist; fendlist:=data`fendlist; p:=data`p; N:=data`N; delta:=data`delta;
  d:=Degree(Q); K:=Parent(x1); 


//  WE IGNORE THE FOLLOWING PART WHICH DEALS WITH BAD POINTS
  /*
  // First make sure that if P1 or P2 is bad, then it is very bad

  if is_bad(P1,data) and not is_very_bad(P1,data) then
    S1:=find_bad_point_in_disk(P1,data);
    _,index:=local_data(S1,data);
    data:=update_minpolys(data,S1`inf,index);
    xt,bt,index:=local_coord(S1,tadicprec(data,e),data);
    S1`xt:=xt;
    S1`bt:=bt;
    S1`index:=index;
    IS1P1,NIS1P1:=tiny_integrals_on_basis(S1,P1,data:prec:=tadicprec(data,e));
    IS1P2,NIS1P2:=$$(S1,P2,data:e:=e);
    IP1P2:=IS1P2-IS1P1;
    NIP1P2:=Ceiling(Minimum([NIS1P1,NIS1P2]));
    return IP1P2,NIP1P2;
  end if;

  if is_bad(P2,data) and not is_very_bad(P2,data) then
    S2:=find_bad_point_in_disk(P2,data);
    _,index:=local_data(S2,data);
    data:=update_minpolys(data,S2`inf,index);
    xt,bt,index:=local_coord(S2,tadicprec(data,e),data);
    S2`xt:=xt;
    S2`bt:=bt;
    S2`index:=index;
    IP1S2,NIP1S2:=$$(P1,S2,data:e:=e);
    IP2S2,NIP2S2:=tiny_integrals_on_basis(P2,S2,data:prec:=tadicprec(data,e));
    IP1P2:=IP1S2-IP2S2;
    NIP1P2:=Ceiling(Minimum([NIP1S2,NIP2S2]));
    return IP1P2,NIP1P2;
  end if;
*/
  // If P1,P2 is bad (hence very bad), use a near boundary point.

  _,index:=local_data(P1,data);
  data:=update_minpolys(data,P1`inf,index);
  _,index:=local_data(P2,data);
  data:=update_minpolys(data,P2`inf,index);
/*
  if is_bad(P1,data) then
    xt,bt,index:=local_coord(P1,tadicprec(data,e),data);
    P1`xt:=xt;       
    P1`bt:=bt;       
    P1`index:=index; 
    Qp:=Parent(P1`x);
    Qpa<a>:=PolynomialRing(Qp);
    K<a>:=TotallyRamifiedExtension(Qp,a^e-p);
    format:=recformat<x,b,inf,xt,bt,index>;
    S1:=rec<format|>;                                                    
    S1`inf:=P1`inf;
    S1`x:=Evaluate(xt,a);
    S1`b:=[Evaluate(bt[i],a):i in [1..d]];
  else
*/
    xt,bt,index:=local_coord(P1,tadicprec(data,1),data);
    P1`xt:=xt;       
    P1`bt:=bt;       
    P1`index:=index; 
    S1:=P1;
//  end if;

/*
  if is_bad(P2,data) then
    xt,bt,index:=local_coord(P2,tadicprec(data,e),data);
    P2`xt:=xt;       
    P2`bt:=bt;       
    P2`index:=index; 
    if not is_bad(P1,data) then
      Qp:=Parent(P2`x);
      Qpa<a>:=PolynomialRing(Qp);
      K<a>:=TotallyRamifiedExtension(Qp,a^e-p);
    end if;
    format:=recformat<x,b,inf,xt,bt,index>;
    S2:=rec<format|>;                                                    
    S2`inf:=P2`inf;
    S2`x:=Evaluate(xt,a);
    S2`b:=[Evaluate(bt[i],a):i in [1..d]];
  else
*/
    xt,bt,index:=local_coord(P2,tadicprec(data,1),data);
    P2`xt:=xt;       
    P2`bt:=bt;       
    P2`index:=index; 
    S2:=P2;
 // end if;

  // Split up the integral and compute the tiny ones.

  tinyP1toS1,NP1toS1:=tiny_integrals_on_basis(P1,S1,data);
  tinyP2toS2,NP2toS2:=tiny_integrals_on_basis(P2,S2,data);

  FS1:=frobenius_pt(S1,data);
  FS2:=frobenius_pt(S2,data);

  tinyS1toFS1,NS1toFS1:=tiny_integrals_on_basis(S1,FS1,data:P:=P1); 
  tinyS2toFS2,NFS2toS2:=tiny_integrals_on_basis(S2,FS2,data:P:=P2); 

  NIP1P2:=Minimum([NP1toS1,NP2toS2,NS1toFS1,NFS2toS2]);  

  // Evaluate all functions.

  I:=[];
  for i:=1 to #basis do
    f0iS1,Nf0iS1:=evalf0(f0list[i],S1,data);
    f0iS2,Nf0iS2:=evalf0(f0list[i],S2,data);
    finfiS1,NfinfiS1:=evalfinf(finflist[i],S1,data);
    finfiS2,NfinfiS2:=evalfinf(finflist[i],S2,data);
    fendiS1,NfendiS1:=evalfend(fendlist[i],S1,data);
    fendiS2,NfendiS2:=evalfend(fendlist[i],S2,data);
    NIP1P2:=Minimum([NIP1P2,Nf0iS1,Nf0iS2,NfinfiS1,NfinfiS2,NfendiS1,NfendiS2]);
    I[i]:=(K!f0iS1)-(K!f0iS2)+(K!finfiS1)-(K!finfiS2)+(K!fendiS1)-(K!fendiS2)-(K!tinyS1toFS1[i])+(K!tinyS2toFS2[i]);
  end for; 

  valIP1P2:=Minimum([Valuation(I[i])/Valuation(K!p):i in [1..#basis]]);

  mat:=(F-IdentityMatrix(RationalField(),#basis));
  valdet:=Valuation(Determinant(mat),p);
  mat:=mat^-1;
  Nmat:=N-valdet-delta;
  valmat:=Minimum([Valuation(e,p):e in Eltseq(mat)]);

  NIP1P2:=Minimum([NIP1P2+valmat,Nmat+valIP1P2]);                            
  
  IS1S2:=Vector(I)*Transpose(ChangeRing(mat,K));    // Solve the linear system.
  IP1P2:=IS1S2+ChangeRing(tinyP1toS1,K)-ChangeRing(tinyP2toS2,K);
  IP1P2,Nround:=round_to_Qp(IP1P2);

  assert Nround ge NIP1P2;                          // Check that rounding error is within error bound.
  
  NIP1P2:=Ceiling(NIP1P2);

  for i:=1 to #basis do
    IP1P2[i]:=IP1P2[i]+O(Parent(IP1P2[i])!p^(NIP1P2));
  end for;

  return IP1P2,NIP1P2;
end function;





coleman_integral:=function(P1,P2,dif,data:e:=1,IP1P2:=0,NIP1P2:=0);

  // Integral of 1-form dif from P1 to P2.

  Q:=data`Q; p:=data`p; N:=data`N; Nmax:=data`Nmax; r:=data`r; W0:=data`W0; Winf:=data`Winf;
  G0:=data`G0; Ginf:=data`Ginf; red_list_fin:=data`red_list_fin; red_list_inf:=data`red_list_inf;
  basis:=data`basis; integrals:= data`integrals; quo_map:=data`quo_map;

  coefs,f0,finf,fend:=reduce_with_fs(dif,Q,p,N,Nmax,r,W0,Winf,G0,Ginf,red_list_fin,red_list_inf,basis,integrals,quo_map); // TODO: handle precision here

  if NIP1P2 eq 0 then 
    IP1P2,NIP1P2:=coleman_integrals_on_basis(P1,P2,data:e:=e);
  end if;
  
  f0P1,Nf0P1:=evalf0(f0,P1,data);
  f0P2,Nf0P2:=evalf0(f0,P2,data);
  finfP1,NfinfP1:=evalfinf(finf,P1,data);
  finfP2,NfinfP2:=evalfinf(finf,P2,data);
  fendP1,NfendP1:=evalfend(fend,P1,data);
  fendP2,NfendP2:=evalfend(fend,P2,data);

  IdifP1P2:=f0P2-f0P1+finfP2-finfP1+fendP2-fendP1;
  NIdifP1P2:=Minimum([Nf0P1,Nf0P2,NfinfP1,NfinfP2,NfendP1,NfendP2]);
  
  for i:=1 to #coefs do
    IdifP1P2:=IdifP1P2+coefs[i]*IP1P2[i];
    NIdifP1P2:=Minimum(NIdifP1P2,NIP1P2+Valuation(coefs[i],p));
  end for;

  NIdifP1P2:=Ceiling(NIdifP1P2);
  IdifP1P2:=IdifP1P2+O(Parent(IdifP1P2)!p^NIdifP1P2);

  return IdifP1P2, NIdifP1P2;

end function;
