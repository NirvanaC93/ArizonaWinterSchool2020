load "coleman.m";
print"";
print "m:=a^2-10;"; m:=a^2-10;print "K=NumberField(m); i.e. K=Q(sqrt(10))";
print"";


print "Q:=y^2-x^3+10*x-9;"; Q:=y^2-x^3+10*x-9;print""; 

print "p:=17;"; p:=17;print"";

print "N:=10;"; N:=10;print"";

print "data:=coleman_data(Q,p,m,N);"; data:=coleman_data(Q,p,m,N);print"";
print "data`r;"; data`r;print"";

print "data`basis;"; data`basis;print"";

print "data`F;"; data`F;print"";
print "P1:=set_point(0,3,data);"; P1:=set_point(0,3,data);
K:=data`Kp;print"";
print "P1bis:=set_point(K.1,3,data);"; P1bis:=set_point(K.1,3,data);print"";
print "P2:=set_point(1,0,data);"; P2:=set_point(1,0,data);print"";
//print "P3:=set_bad_point(1,[1,0],false,data);"; P3:=set_bad_point(1,[1,0],false,data);print"";
//print "P4:=set_bad_point(0,[1,0],true,data);"; P4:=set_bad_point(0,[1,0],true,data);print"";
print "is_bad(P1,data);"; is_bad(P1,data);print"";
print "is_bad(P2,data);"; is_bad(P2,data);print"";
print "is_very_bad(P2,data);"; is_very_bad(P2,data);print"";
print "lie_in_same_disk(P1,P2,data);"; lie_in_same_disk(P1,P2,data);print"";
print "lie_in_same_disk(P1,P1bis,data);"; lie_in_same_disk(P1,P1bis,data);print"";
//lie_in_same_disk(P1,P2,data);
tiny_integrals_on_basis(P1,P1bis,data);