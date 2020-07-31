function [scores]=DIAMOnD_combinational_drugs_scores(gene_list,Targets1,Targets2,Breast_genes_set,CC);  

%*******************************************************************
%   function:calculate the combinational scores of two drugs
%   Input:
%         Targets1 and Targets2:Targets for drug 1 and drug 2
%         Breast_genes_set:Breaset cancer genes
%         gene_list:All genes set
%         CC:The adjacency matrix

%   Output:
%          scores:predited scores
%*******************************************************************
A=CC;
k=find(sum(A,2)~=0);


[c,d]=ismember(Breast_genes_set,gene_list);
d(d==0)=[];

%********Targets of Drug 1 overlapping with breast cancer genes******


[a1,b1]=ismember(Targets1,gene_list);
b1(b1==0)=[];
%The neighborhood nodes of Targets 1 for drug A
%c1=A(:,b1);
%s1=[find(sum(c1,2)~=0);b1];

%Module of drug A
[s1,density]=DIAMOnD_module(A,b1);


t1=intersect(s1,d);
AA=sum(sum(A(t1,t1)))/2;

Nmax=100;A0=[];
for i=1:Nmax

t0= intersect(k(randi([1 length(k)],length(s1),1)),d);
A0(i,1)=sum(sum(A(t0,t0)))/2;
    
    
end
z_score=(AA -mean(A0))/(std(A0));
scores1=(1-1+normcdf(z_score));

%********Targets of Drug 2 overlapping with breast cancer genes******

[a2,b2]=ismember(Targets2,gene_list);
b2(b2==0)=[];
%The neighborhood nodes of Targets 2 for drug B

%c2=A(:,b2);
%s2=[find(sum(c2,2)~=0);b2];
[s2,density]=DIAMOnD_module(A,b2);


t2=intersect(s2,d);
AA=sum(sum(A(t2,t2)))/2;

Nmax=100;A0=[];
for i=1:Nmax

t0= intersect(k(randi([1 length(k)],length(s2),1)),d);
A0(i,1)=sum(sum(A(t0,t0)))/2;
    
    
end
z_score=(AA -mean(A0))/(std(A0));
scores2=(1-1+normcdf(z_score));


%******scores of module 1 and module 2******************
t3=intersect(s2,s1);
AA=sum(sum(A(t3,t3)))/2;
t4=unique([s1;s2]);
NN=sum(sum(A(t4,t4)))/2;

scores3=(NN-AA)/NN;

scores=scores1*scores2*scores3;


end


