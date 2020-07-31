function [ new_x,index ] = Nar_Opti_weight_nc( dz,NN,scores,lamda  )
%non-lineaqr controllability of undirected networks
%   input:
%         z:network structure
%         scores:the scores of each node
%         lamda:the parameter of our PDC
%   Output:
%         lc_index:the number of driver nodes

%***********************solve the problem**************************************
%**********************MATLAB2014******************************

z=dz;
N=NN;

%lamda=0;
%N=max(max(z));
%scores=ones(N,1);


N1=N;
[N2,~]=size(z);
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N2,N1);
for i=1:N2
    
        A_adjacent(i,z(i,1))=1;
        A_adjacent(i,z(i,2))=1;
       
end
A=-A_adjacent;
b=-ones(size(A,1),1);
f = ones(size(A,2),1)-lamda*scores;%

model.A = sparse(A);
model.obj = f';
model.rhs = b;
model.sense=[];
for i=1:size(A,1)
model.sense = [model.sense,'<'];;
end

model.modelsense = 'min';


lb=zeros(size(A,2),1);
ub=ones(size(A,2),1);

 model.lb = lb;
 model.ub = ub;

n = size(model.A, 2);
model.vtype = repmat('C', n, 1);

intcon = [1:size(A,2)];
model.vtype(intcon) = 'I';
%model.varnames = names;

%gurobi_write(model, 'mip1.lp');

params.outputflag = 0;

result = gurobi(model, params);

%disp(result);
x=result.x;
object=result.objval;
index=sum(x);

new_x=x;

end

