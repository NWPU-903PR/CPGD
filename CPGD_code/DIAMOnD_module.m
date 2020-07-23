function [s,Density]=DIAMOnD_module(A,b1)
%forming the targets module by using DIAMOnD_module
%   A:adjacency matrix;
%   b:targets;
%The neighborhood nodes of Targets 1 for drug A


A(A~=0)=1;
degree=sum(A,2);
N=length(find(degree~=0));

Density=[];All_ss=[];
d0=sum(sum(A));

for y=1:10

    All_ss{y,1}=b1;
   
    
   cand=sum(sum(triu(A(b1,b1))));
   
   
   c1=A(:,b1);
s1=[find(sum(c1,2)~=0);b1];
ss=unique(s1);

 Density(y,1)=cand/(length(setdiff(find(sum(c1,2)~=0),b1)));    

p_value=[];
for i=1:length(ss)
    
    c1=A(:,ss(i));
    t=[find(sum(c1,2)~=0);i];
    tt=unique(t);
    p_value(i,1)=abs(1-sum(hygepdf(0:length(intersect(tt,ss))-1,N,length(ss),length(tt))));

    
end

[a,b]=find(p_value<0.001);


b1=[b1;ss(a)];
b1=unique(b1);

end

[w,p]=max(Density);
s=All_ss{p,1};



end



