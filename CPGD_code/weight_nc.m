function [ x,index ] = weight_nc( z,NNN,scores,lamda )
%non-lineaqr controllability of undirected networks
%   input:
%         z:network structure
%         scores:the scores of each node
%         lamda:the parameter of our PDC
%   Output:
%         lc_index:the number of driver nodes

%***********************solve the problem**************************************
%**********************MATLAB2014*****************************

N1=NNN;
[N2,~]=size(z);
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N2,N1);
for i=1:N2
    
        A_adjacent(i,z(i,1))=1;
        A_adjacent(i,z(i,2))=1;
       
end
A=-A_adjacent;
b=-ones(size(A,1),1);
f = lamda*ones(size(A,2),1)-scores;%
intcon = [1:size(A,2)];
lb=zeros(size(A,2),1);
ub=ones(size(A,2),1);
[x,fval,exitflag,output]  = intlinprog(f,intcon,A,b,[],[],lb,ub);

%**********************************************************
index=length(find(x(:,1)~=0));

end

