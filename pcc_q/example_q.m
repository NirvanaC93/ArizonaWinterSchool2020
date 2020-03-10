////////////////////////////////////////////////// 
// This code is part of the pcc_q MAGMA library //
//                                              //
// copyright (c) 2017 Jan Tuitman               //
//////////////////////////////////////////////////


//////////////
// EXAMPLES //
//////////////


/////////////////////////////
// 1. nondegenerate curves //
/////////////////////////////


/////////////////////////////
// 1.1 some plane quartics //
/////////////////////////////


// q:=3^10;	
// q:=5^10; 	
// q:=7^10;	
// q:=11^10;    
// q:=3^20;	
// q:=5^20;	

// Q:=y^4+y^3*((3*a^9 + 5*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 4*a^2 + 8*a + 2)*x+(6*a^9 + a^8 + 3*a^7 + 6*a^6 + 2*a^5 + 5*a^4 + 3*a^3 + 5*a^2 + 3*a + 5))+y^2*((4*a^9 + 3*a^8 + 6*a^7 + 2*a^6 + 6*a^5 + 2*a^3 + 2*a^2)*x^2+(a^9 + a^8 + 2*a^7 + 4*a^6 + 2*a^5 + a^4 + 2*a^3 + 4*a^2 + 6*a + 2)*x+(5*a^8 + 2*a^7 + 2*a^6 + 3*a^2 + a))+y*((7*a^9 + 2*a^8 + 7*a^7 + 10*a^6 + 3*a^5 + 5*a^4 + a^3 + 10*a^2 + 4*a + 6)*x^3+(4*a^9 + 6*a^8 + 4*a^7 + 3*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 6*a^2 + 1)*x^2+(a^9 + 4*a^8 + 3*a^7 + 6*a^6 + 6*a^5 + 6*a^4 + 3*a^3 + a^2 + 9*a + 10)*x+(7*a^9 + 2*a^8 + 7*a^7 + 10*a^6 + 3*a^5 + 5*a^4 + a^3 + 10*a^2 + 4*a + 6))+(a^9 + 6*a^8 + a^7 + a^6 + 6*a^5 + a^3 + 5*a^2 + 2*a + 6)*x^4+(6*a^9 + a^8 + 3*a^7 + 6*a^6 + 2*a^5 + 5*a^4 + 3*a^3 + 5*a^2 + 3*a + 5)*x^3+(2*a^9 + a^8 + 6*a^7 + 3*a^6 + 5*a^5 + 6*a^4 + 3*a^3 + a)*x^2+(5*a^9 + 3*a^8 + a^7 + 2*a^6 + 4*a^5 + a^4 + 3*a^3 + 5*a^2 + 2)*x+(a^9 + a^8 + 2*a^7 + 4*a^6 + 2*a^5 + a^4 + 2*a^3 + 4*a^2 + 6*a + 2);


/////////////////////////////
// 1.2 some C_{a,b} curves //
/////////////////////////////


/////////////
// C_{3,4} //
/////////////


// q:=5^10; 	
// q:=7^10;	
// q:=11^10;		
// q:=5^20; 	
// q:=7^20;     
// q:=11^20;     
// q:=5^40;	

// Q:=y^3 + ((2*a^9 + a^8 + a^7 + 5*a^6 + 4*a^5 + 9*a^3 + 4*a + 1)*x + 4*a^9 + 7*a^8 + 7*a^7 + 2*a^6 + 5*a^5 + 5*a^4 + 9*a^3 + 8*a^2 + 9*a + 5)*y^2 + ((a^9 + 4*a^8 + 3*a^7 + 6*a^6 + 6*a^5 + 6*a^4 + 3*a^3 + a^2 + 9*a + 10)*x^2 + (7*a^9 + 2*a^8 + 7*a^7 + 10*a^6 + 3*a^5 + 5*a^4 + a^3 + 10*a^2 + 4*a + 6)*x + (2*a^9 + 8*a^8 + 6*a^7 + 4*a^6 + a^5 + 2*a^3 + 7*a^2 + 4*a + 3))*y +(4*a^9 + 6*a^8 + 4*a^7 + 3*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 6*a^2 + 1)*x^4 +(5*a^9 + a^8 + 5*a^7 + 3*a^6 + 10*a^5 + 4*a^4 + 3*a^3 + 3*a^2 + 3*a + 6)*x^3 + (3*a^9 + 5*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 4*a^2 + 8*a + 2)*x^2 + (3*a^9 + 6*a^8 + 6*a^7 + 6*a^6 + 5*a^5 + 8*a^4 + 2*a^3 + 4*a^2 + 3*a)*x + 3*a^9 + 3*a^8 + a^7 + 10*a^6 + 3*a^4 + 5*a^3 + 3*a^2 + a + 1;


/////////////
// C_{3,5} //
///////////// 


// q:=5^10;	
// q:=7^10;	
// q:=11^10;		
// q:=13^10;   	
// q:=5^20;	
// q:=7^20;     
// q:=11^20;    

// Q:= y^3 + ((2*a^9 + a^8 + a^7 + 5*a^6 + 4*a^5 + 9*a^3 + 4*a + 1)*x + 4*a^9 + 7*a^8 + 7*a^7 + 2*a^6 + 5*a^5 + 5*a^4 + 9*a^3 + 8*a^2 + 9*a + 5)*y^2 + ((7*a^9 + 2*a^8 + 7*a^7 + 10*a^6 + 3*a^5 + 5*a^4 + a^3 + 10*a^2 + 4*a + 6)*x^3 + (a^9 + 4*a^8 + 3*a^7 + 6*a^6 + 6*a^5 + 6*a^4 + 3*a^3 + a^2 + 9*a + 10)*x^2 + (7*a^9 + 2*a^8 + 7*a^7 + 10*a^6 + 3*a^5 + 5*a^4 + a^3 + 10*a^2 + 4*a + 6)*x + (2*a^9 + 8*a^8 + 6*a^7 + 4*a^6 + a^5 + 2*a^3 + 7*a^2 + 4*a + 3))*y + (2*a^9 + a^8 + 6*a^7 + 3*a^6 + 5*a^5 + 6*a^4 + 3*a^3 + a)*x^5 + (4*a^9 + 6*a^8 + 4*a^7 + 3*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 6*a^2 + 1)*x^4 +(5*a^9 + a^8 + 5*a^7 + 3*a^6 + 10*a^5 + 4*a^4 + 3*a^3 + 3*a^2 + 3*a + 6)*x^3 + (3*a^9 + 5*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 4*a^2 + 8*a + 2)*x^2 + (3*a^9 + 6*a^8 + 6*a^7 + 6*a^6 + 5*a^5 + 8*a^4 + 2*a^3 + 4*a^2 + 3*a)*x + 3*a^9 + 3*a^8 + a^7 + 10*a^6 + 3*a^4 + 5*a^3 + 3*a^2 + a + 1;


///////////////////////////////
// 1.3 d_x=3, d_y=3, genus 4 //
///////////////////////////////   


// q:=5^10;     
// q:=7^10;     
// q:=11^10;    
// q:=5^20;     
// q:=7^20;     
// q:=11^20;    

// Q:= y^3 + ((4*a^9 + 3*a^8 + 6*a^7 + 2*a^6 + 6*a^5 + 2*a^3 + 2*a^2)*x^3 + (a^9 + 6*a^8 + a^7 + a^6 + 6*a^5 + a^3 + 5*a^2 + 2*a + 6)*x^2 + (5*a^8 + 4*a^7 + 5*a^6 + 5*a^5 + a^4 + 4*a^3 + 6*a^2 + 6*a + 1)*x + (3*a^9 + 3*a^8 + 6*a^7 + a^6 + 4*a^5 + 3*a^4 + a^3 + 5*a^2 + 4))*y^2 + ((4*a^9 + a^8 + 2*a^7 + 2*a^6 + 2*a^4 + 3*a^3 + 5*a^2 + 5*a + 2)*x^3 + (2*a^9 + 2*a^8 + 3*a^7 + 6*a^6 + 2*a^5 + 6*a^4 + 6*a^3 + 6*a^2 + 6*a)*x^2 + (6*a^9 + a^8 + 3*a^7 + 6*a^6 + 2*a^5 + 5*a^4 + 3*a^3 + 5*a^2 + 3*a + 5)*x + (5*a^9 + a^8 + a^7 + 3*a^6 + 4*a^5 + 6*a^3 + a^2 + 6*a + 2))*y + (5*a^9 + a^8 + 5*a^7 + 3*a^6 + 3*a^5 + 4*a^4 + 3*a^3 + 3*a^2 + 3*a + 6)*x^3 + (3*a^9 + 5*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 4*a^2 + 4*a + 2)*x^2 + (3*a^9 + 6*a^8 + 6*a^7 + 6*a^6 + 5*a^5 + 3*a^4 + 2*a^3 + 4*a^2 + 3*a)*x + 3*a^9 + 3*a^8 + a^7 + 6*a^6 + 3*a^4 + 5*a^3 + 3*a^2 + a + 1;


///////////////////////////////
// 1.4 d_x=3, d_y=4, genus 6 //
///////////////////////////////

// q:=5^10;  
// q:=7^10;  
// q:=11^10;

// Q:= y^3 + ((6*a^9 + 4*a^8 + 6*a^7 + 3*a^6 + 2*a^5 + 2*a^4 + 2*a^3 + 6*a^2 + a + 1)*x^4 + (6*a^9 + 5*a^8 + 3*a^7 + 2*a^6 + 5*a^5 + 6*a^3 + 4*a^2 + 3)*x^3 + (2*a^9 + 2*a^7 + 2*a^6 + 3*a^5 + 6*a^4 + 6*a^3 + 6*a^2 + 3*a + 4)*x^2 + (3*a^9 + a^8 + 2*a^7 + 5*a^6 + a^4 + 3*a^3 + 6*a^2 + 6*a + 3)*x + (6*a^9 + 6*a^8 + 3*a^6 + 3*a^3 + 2*a^2 + 3*a + 6))*y^2 + ((2*a^8 + 3*a^7 + 4*a^6 + a^5 + 5*a^4 + 3*a^3 + 5*a^2 + 2*a + 4)*x^4 + (2*a^9 + a^8 + 6*a^7 + 3*a^6 + 5*a^5 + 6*a^4 + 3*a^3 + a)*x^3 + (6*a^9 + a^6 + a^5 + 4*a^4 + 5*a^3 + 3*a^2 + a + 3)*x^2 + (5*a^9 + a^8 + 5*a^7 + 5*a^6 + 4*a^5 + 3*a^4 + 4*a^2 + a + 6)*x + (4*a^9 + a^8 + a^7 + 6*a^6 + 6*a^5 + 3*a^4 + 4*a^2 + a + 5))*y + (4*a^9 + a^8 + a^7 + 4*a^6 + 2*a^5 + a^4 + 3*a^3 + 6*a^2 + 3*a + 2)*x^4 + (5*a^9 + 4*a^8 + 6*a^7 + 3*a^6 + 6*a^5 + 6*a^3 + 5*a^2 + 4*a + 1)*x^3 + (2*a^9 + 2*a^8 + 6*a^7 + 2*a^6 + 3*a^5 + 5*a^4 + 2*a^2 + 3*a + 5)*x^2 + (4*a^9 + 6*a^8 + 4*a^7 + 3*a^6 + 3*a^5 + 2*a^4 + 6*a^3 + 6*a^2 + 1)*x + a^9 + 5*a^8 + 5*a^6 + 6*a^4 + 2*a^3 + 6*a^2 + 5*a + 5;


///////////////////////////////
// 1.5 d_x=3, d_y=5, genus 8 //
///////////////////////////////


// q:=5^10; 	
// q:=7^10; 	
// q:=11^10;	

// Q:=y^3 + ((a^9 + 5*a^7 + 3*a^5 + 6*a^4 + 4*a^3 + 2*a^2 + 5*a + 1)*x^5 + (5*a^9 + 5*a^8 + 2*a^7 + 2*a^6 + 3*a^5 + a^4 + 6*a^3 + 4*a + 4)*x^4 + (2*a^9 + 6*a^7 + 6*a^6 + 2*a^5 + 6*a^4 + 5*a^3 + 6)*x^3 + (3*a^9 + 2*a^8 + 3*a^7 + 3*a^6 + a^5 + 4*a^4 + 5*a^3 + 4*a^2 + 3*a + 3)*x^2 + (5*a^9 + 3*a^8 + a^7 + 2*a^6 + 4*a^5 + a^4 + 3*a^3 + 5*a^2 + 2)*x + (4*a^8 + 2*a^7 + 4*a^6 + a^4 + 4*a^3 + a^2 + 2*a + 4))*y^2 + ((2*a^9 + 3*a^8 + 3*a^7 + 6*a^6 + 6*a^5 + 6*a^4 + 4*a^3 + 5*a^2 + 6*a)*x^5 + (5*a^9 + 3*a^8 + 2*a^6 + 2*a^5 + 4*a^4 + 2*a^3 + 4*a^2 + 3*a + 6)*x^4 + (3*a^9 + 3*a^8 + 6*a^7 + 5*a^6 + 3*a^5 + 3*a^4 + 5*a^3 + 4*a^2 + 4*a + 1)*x^3 + (2*a^9 + 2*a^8 + 5*a^7 + 5*a^6 + 5*a^5 + 6*a^4 + a^3 + a^2 + 2*a + 2)*x^2 + (3*a^8 + 3*a^6 + 3*a^5 + 5*a^3 + 4*a^2 + 4*a + 2)*x + (4*a^9 + 2*a^8 + 5*a^7 + 5*a^6 + 2*a^5 + 5*a^4 + 6*a^3 + 2*a + 4))*y + (a^9 + a^8 + 2*a^7 + 4*a^6 + 2*a^5 + a^4 + 2*a^3 + 4*a^2 + 6*a + 2)*x^5 + (4*a^9 + 5*a^7 + a^6 + a^5 + 3*a^4 + 2*a^3 + 6*a + 6)*x^4 + (4*a^9 + 4*a^8 + 4*a^7 + a^6 + a^5 + 5*a^4 + 2*a^3 + a^2 + 2*a)*x^3 + (5*a^9 + 5*a^7 + 6*a^6 + 3*a^5 + 6*a^4 + 4*a^3 + 3*a^2 + 6*a)*x^2 + (5*a^8 + 2*a^7 + 2*a^6 + 3*a^2 + a)*x + a^9 + 6*a^8 + 6*a^7 + 2*a^6 + 6*a^5 + 4*a^4 + 3*a^3 + 5*a + 2;


///////////////////////////////
// 1.6 d_x=4, d_y=4, genus 9 //
///////////////////////////////


// q:=3^5;  
// q:=5^5   
// q:=7^5;  
// q:=3^10;   
// q:=5^10; 
// q:=7^10;  

// Q:=y^4 + ((5*a^4 + 8*a^2 + 4*a + 9)*x^4 + (6*a^4 + 10*a^2 + 3*a + 2)*x^3 + (8*a^4 + 10*a^3 + 8*a^2 + 6*a + 7)*x^2 + (5*a^4 + 7*a^3 + 10*a^2 + 5*a + 8)*x + (4*a^4 + 5*a^3 + 10*a + 4))*y^3 + ((10*a^4 + 3*a^3 + 5*a^2 + 3*a)*x^4 + (7*a^4 + 9*a^3 + 9*a + 1)*x^3 + (10*a^3 + 7*a^2 + 4*a + 5)*x^2 + (5*a^4 + 5*a^3 + 7*a^2 + 2*a + 1)*x + (6*a^4 + 7*a^3 + a^2 + a + 3))*y^2 + ((7*a^4 + 7*a^3 + a^2 + 7*a + 4)*x^4 + (6*a^4 + 10*a^3 + 5*a^2 + 5*a + 10)*x^3 + (4*a^4 + 4*a^3 + 5*a^2 + 4)*x^2 + (2*a^4 + 10*a^3 + 7*a^2 + 5*a + 8)*x + (6*a^4 + 5*a + 6))*y + (7*a^4 + 7*a^3 + 10*a^2 + 6*a + 4)*x^4 + (4*a^4 + 6*a^3 + 9*a^2 + 8*a + 7)*x^3 + (5*a^4 + 8*a^3 + 6*a^2 + 10*a + 2)*x^2 + (4*a^4 + 8*a^3 + 5*a^2 + 5*a + 1)*x + 9*a^4 + 10*a^3 + 2*a + 2;


////////////////////////////////
// 1.7 d_x=4, d_y=5, genus 12 //
////////////////////////////////

// q:=3^5; 
// q:=5^5; 
// q:=7^5; 

// Q:=y^4 + ((2*a^4 + 10*a^3 + 4*a^2 + 6)*x^5 + (4*a^4 + 5*a^3 + 4*a^2 + 5*a + 8)*x^4 + (10*a^4 + 8*a^3 + 9*a^2 + 5*a + 4)*x^3 + (9*a^4 + 2*a^3 + 3*a^2 + 3)*x^2 + (6*a^4 + 7*a^3 + 6*a^2 + 10*a + 9)*x + (a^4 + 7*a^3 + 9*a^2 + 9*a + 6))*y^3 + ((7*a^4 + 3*a^3 + 9)*x^5 + (6*a^4 + 3*a^2 + 8*a + 1)*x^4 + (a^4 + 5*a^3 + 6*a^2 + 2)*x^3 + (9*a^3 + 10*a^2 + 8*a + 5)*x^2 + (a^4 + 5*a^3 + 2*a^2 + 10*a + 5)*x + (4*a^4 +7*a^3 + 5*a^2 + a + 5))*y^2 + ((8*a^4 + 3*a^3 + 3*a^2 + a + 3)*x^5 + (8*a^4 + 2*a^3 + a^2 + 7*a)*x^4 + (3*a^4 + 10*a^3 + 4*a^2 + 8*a + 7)*x^3 + (10*a^4 + 2*a^3 + 10*a^2 + 7*a + 6)*x^2 + (8*a^4 + 9*a^3 + 4*a + 9)*x + (3*a^4 + 8*a^3 + 3*a^2 + 8*a + 7))*y + (6*a^4 + 6*a^3 + 2*a^2 + 5)*x^5 + (5*a^4 + 9*a^3 + 3*a^2 + 3)*x^4 + (6*a^4 + 9*a^3 + a^2 + 10*a + 6)*x^3 + (2*a^4 + 7*a^3 + a^2 + a + 3)*x^2 + (8*a^4 + a^3 + 9*a^2 + 6*a + 2)*x + 6*a^4 + 4*a^3 + 5*a^2 + 4*a + 4;


//////////////////////////////////////////////////////////////////////////////////////////////////
// more general (i.e. degenerate) examples for genus 4,5 can be found in the goodmodels package //
//////////////////////////////////////////////////////////////////////////////////////////////////