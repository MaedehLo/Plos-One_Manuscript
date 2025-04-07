function S = Index_Orientation(A,alpha)

 
f1=0;
f2=0;
for i=1:size(A,1)
  
    f1=f1+((cosd(A(i,1)-alpha))^2);   
    f2=f2+1; 
end
 S=2*(f1/f2)-1;  % S=0 (totally isotropic)    S=1(perfectly aligned)
 
end