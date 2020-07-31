function [ CPGD_rank_drugs ] = condorcet( All_parameter_num,All_drugs)
%function:
%         rank the candidate genes with the condorcet method
%Input:
%      new_result_driver_gene_module:the filtered information of individual sample
%      node0:gene information
%Output:
%      SCS_rank_genes:the rank information of genes


%list all the candidate driver mutation
sample_driver_rank=[];
node0=All_drugs;

All_A=zeros(length(node0));

for k=1:5

sample_driver_rank=All_parameter_num{k,1};

    A0=zeros(length(node0));
 
for i=1:length(sample_driver_rank)
    for j=i+1:length(sample_driver_rank)
        
                 
                 x=sample_driver_rank(i,1);
           
                 y=sample_driver_rank(j,1);
            
           
                 value_x=sample_driver_rank(i,2);
           
                 value_y=sample_driver_rank(j,2);


            
             if value_x>value_y
                 A0(x,y)=1;
             end
             
             if value_y>value_x
                 A0(y,x)=1;
             end
             
             
         end
         
        
    end

   

All_A=All_A+A0;


end

ind=sum(All_A,2);
[value,PSCS_rank_genes_ind]=sort(ind,'descend');
value(value~=0)=1;
CPGD_rank_drugs=node0(PSCS_rank_genes_ind,:);




end