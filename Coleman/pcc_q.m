////////////////////////////////////////////////// 
// This code is part of the pcc_q MAGMA library //
//                                              //
// copyright (c) 2017 Jan Tuitman               //
//////////////////////////////////////////////////


print "";
print "*****************************";
print "* pcc_q Magma library v2.17 *";
print "*                           *";
print "* by Jan Tuitman            *";
print "*****************************";
print "";


Z:=IntegerRing();
Za<a>:=PolynomialRing(Z);
Zax<x>:=PolynomialRing(Za);
Zaxy<y>:=PolynomialRing(Zax);


load "auxpolys_q.m";
load "coho_q.m";
load "froblift_q.m";
load "reductions_q.m";
load "zeta_q.m";
print "";


num_zeta:=function(Q,q: N:=0, m:=0, verbose:=false, exactcoho:=false, W0:=0, Winf:=0);

  if IsPrimePower(q) then
    p:=Factorization(q)[1][1];
    n:=Factorization(q)[1][2];
  else
    error "q is not a power of a prime";
  end if;

  if m eq 0 then
    m:=Za!DefiningPolynomial(FiniteField(q));
  end if;

  N0:=N;

  ResetMaximumMemoryUsage();
  t0:=Cputime();

  K:=NumberField(m);
  Kx:=PolynomialRing(K);
  Kxy:=PolynomialRing(Kx);
  Kxyz:=LaurentSeriesRing(Kxy);

  if verbose then
    print "Computing Delta,r,s,g,W^0,W^inf,G:";
  end if;

  r,Delta,s:=auxpolys(Q,p,n,K);
  g:=genus(Q,p,n,m);

  if W0 eq 0 then
    W0:=mat_W0(Q,Kxy);
  else
    W0:=push_to_Kt_mat(W0,K);
  end if;
  if Winf eq 0 then
    Winf:=mat_Winf(Q,Kxy);
  else
    Winf:=push_to_Kt_mat(Winf,K);
  end if;

  W0inv:=W0^(-1); 
  Winfinv:=Winf^(-1); 

  if (Degree(r) lt 1) or (not smooth(r,p,n,m)) or (not (is_integral(W0,p,n) and is_integral(W0inv,p,n) and is_integral(Winf,p,n) and is_integral(Winfinv,p,n))) then
    error "Bad model for curve";
  end if;

  G:=con_mat(Q,Delta,s,K);
  G0:=W0*Evaluate(G,Parent(W0[1,1]).1)*W0^(-1)+ddx_mat(W0)*W0^(-1);
  Ginf:=Winf*Evaluate(G,Parent(Winf[1,1]).1)*Winf^(-1)+ddx_mat(Winf)*Winf^(-1);
  e0,e0list,resG0list:=fin_ram_ind(r,G0,Kx);
  einf,einflist,resGinf:=inf_ram_ind(Ginf,Kx);

  if verbose then
    print "Time (s) :    ", Cputime(t0);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing basis H^1(X):";
  end if;

  if N eq 0 then
    prov_prec:=true;
    N:=provable_prec(Q,p,n,g,W0,Winf,e0,einf);
  else
    prov_prec:=false;
  end if;

  J0,T0,T0inv:=jordan_0(p,n,m,N,r,e0list,resG0list,Kx,exactcoho);
  Jinf,Tinf,Tinfinv:=jordan_inf(p,n,m,N,einflist,resGinf,exactcoho);

  basis,quo_map:=basis_coho(Q,p,n,m,N,r,W0,Winf,G0,Ginf,J0,Jinf,T0inv,Tinfinv,Kxy,exactcoho);
  v:=Minimum([0] cat [Valuation(j,p) : j in &cat([Eltseq(i) : i in Eltseq(quo_map)])]);
  if (v lt 0 and prov_prec) then
    N:=provable_prec(Q,p,n,g,W0,Winf,e0,einf:val:=v);
  end if;

  if verbose then
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing Frobenius lift:";
  end if; 

  frobmatb0r:=froblift(Q,p,n,m,N-1,r,Delta,s,W0);

  if verbose then
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing reduction matrices:";
  end if;

  redlistfin,redlistinf:=red_lists(Q,p,n,m,N,r,W0,Winf,G0,Ginf,e0,einf,J0,Jinf,T0,Tinf,T0inv,Tinfinv,Kx);

  if verbose then
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing Frobenius matrix:";
  end if;

  Fp:=frobmatrix(Q,p,n,m,N,r,W0,Winf,G0,Ginf,frobmatb0r,redlistfin,redlistinf,basis,quo_map,Kx,verbose);
  Fq:=galois_norm(Fp,p,n,m,N);

  if verbose then
    print "";
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  chi:=revcharpoly(Fq,p,n);

  if verbose then
    print "Total time (s):", Cputime(t0), "\n";
    print "The p-adic precision N:", N, "\n";
  end if;

  if not test_num_zeta(chi,q) then
    print "Q =", Q, "p =", p, "n =", n;
    print "Output seems to be wrong. Please report this example to jan.tuitman@kuleuven.be.";
    return 0; 
  end if;

  ZT<T>:=PolynomialRing(IntegerRing());
  chi:=ZT!chi;

  return chi;

end function;
