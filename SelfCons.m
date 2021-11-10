b=7; 
h=1;
Dev=1;
Temp=0.01;
Lambda=0.1;

Sparam=1/Lambda-0.5*log(Temp);
Op=1;

hh=[];
OOp=[];

for dh=0:100
    
    h=dh*0.1
    
   Dev=1;


while (Dev>0.0001)

OpNew=0;

for n=0:ceil(100/Temp)    

OmegaM=(2*n+1)*pi*Temp; 

Gorkov= OPfun(OmegaM,Op,h,b) -Op/OmegaM ; 

OpNew=OpNew + pi*Lambda*Temp*Gorkov ;

end

OpNew=OpNew+Lambda*Sparam*Op;     

Dev=abs(Op-OpNew)  ;

Op=OpNew;

end 

hh=[hh h];
OOp=[OOp Op]

end

file=['Op'] ;  
save (file,'OOp', 'hh');
