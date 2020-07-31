function [Result_predict_drug,All_Sample_drug,Personalized_side_effect] = CPGD( expression_tumor_fileName,expression_normal_fileName,lamda )
%CPGD outputs the patient-specific drug profiles
%   Input:
%         expression_tumor_fileName:the tumor expression data in cancer
%         expression_normal_fileName:the normal expression data in cancer 
%         lamda:the balance parameter between the scores objective and minimum objective,the bigger the value is, the more preferable we choose the minimum objective.
%               the fefault value is 0.01.


%  Output:
%         Output1:
%         Result_drug_combinations:individual combinational drugs
%         (DCDB);the first column is the sample name with paired data
%         (normal and tumor) and the second column is the ranked
%         combinational drug name in descend with defined scores.The
%         defined scores are the number of targeted personalized driver
%         genes.

%         Output2:
%         Result_drug_targets: individual drug-target genes；the 1-4 column is efficacious synergistic drug pairs including Drug A, targets of drug A, Drug B, targets of drug B, and the 5 colunm is the related combinational drugs (name in DCDB), the 6 colunm is the corresponding synergistic scores

%         Output3:
%         Result_patients_side_effect：the first column is the sample name and the second colunm denote the aggrevating effect and the third colunm denotes the improving effect;

%************************part1:LOAD BRCA sample data and network data************************
%********************obtain the paired expression data******************
%read the original text of tumor and normal

%The temple data in our study


%*************************tumor****************************

[tumor,name_tumor,~]=xlsread(expression_tumor_fileName);

gene_list=name_tumor(2:end,1);




Sample_name_tumor=name_tumor(1,2:end);
tumor_data=tumor(:,2:end);

%*************************normal****************************

[normal,name_normal,~]=xlsread(expression_normal_fileName);
Sample_name_normal=name_normal(1,2:end);
normal_data=normal(:,2:end);


%**********GIN interaction Network************************
load('GIN_network_information.mat', 'edge0')
PPI=edge0;

[x1,y1]=ismember(PPI(:,1),gene_list);
[x2,y2]=ismember(PPI(:,2),gene_list);

P=unique(PPI);
[xx,yy]=ismember(P,gene_list);

y=y1.*y2;
z=[y1 y2];
z(find(y==0),:)=[];

N1=length(gene_list);
[N2,~]=size(z);


%calculate the adjacency matrix of PPI 

Net_adjacent=zeros(N1,N1);
for i=1:N2
    i
         Net_adjacent(z(i,2),z(i,1))=1;  %undirected gene-gene interaction network
         Net_adjacent(z(i,1),z(i,2))=1;
       
end




%************Construct Drug-Target interaction network**************

[data]=importdata('Combination_drug_genes.xlsx');
network=data.Sheet1;
Drug_name=unique(network(:,1));

[~,z2]=ismember(network(:,2),gene_list);

[z1]=cellfun('isempty',network(:,1));
z=intersect(find(z1==0),find(z2~=0));
index=unique(z);
Final_net=network(index,:);

Final_drug_name=unique(Final_net(:,1));%drug name

[~,z1]=ismember(Final_net(:,1),Final_drug_name);
[~,z2]=ismember(Final_net(:,2),gene_list);
z=[z1 z2];

Matrix_Drug_target=zeros(length(gene_list),length(Final_drug_name));%drug-target interaction name

for i=1:length(z)
    
    Matrix_Drug_target(z2(i,1),z1(i,1))=1;
end


%************part2:Construct personalized network and the main control procedure****************

new_T=tumor_data;
new_N=normal_data;
new_gene=gene_list;
cx=[1:length(gene_list)]';

%******************Output1:Result_predict_drug,the prediction of personalized drug combinations*******
%read the size of data

[row,colunm]=size(new_T);
Ref=new_N;  
k0=sum(Net_adjacent);%
[b0,a0]=find(k0~=0);
subnetwork_nodes0=a0';


new_T0=new_T(subnetwork_nodes0,:);
new_N0=new_N(subnetwork_nodes0,:);
Ref0=Ref(subnetwork_nodes0,:);

NC_nodes=[];
NC_genes=[];

Final_Sample_Tumor=Sample_name_tumor;

for i=1:size(Final_Sample_Tumor,2)
%for i=1:2
%for i=1:size(new_T,2)
    
    i
        

    %for the i-th sample
    %construct the tumor SSN**********
    
    
    sample_tumor=new_T0(:,i);
    
    
    [R0,P]=SSN(sample_tumor,Ref0);
    
    P(isnan(P))=0;
    P(abs(P)>=0.05)=0;
    P(P~=0)=1;
   
    
    
    R=abs(R0);
    R(isnan(R))=0;
    
    %sample_network{i,1}=C;
    D0=abs(P.*R);
    
    
    %sample_network{i,1}=C;

CC=abs(D0.*Net_adjacent(subnetwork_nodes0,subnetwork_nodes0));

CC(isnan(CC))=0;CC(CC==inf)=0;

[z1,z2]=find(triu(CC)~=0);
p=[];
for ii=1:length(z1)
p(ii,1)=CC(z1(ii,1),z2(ii,1));
end

All_CC{i,1}=[subnetwork_nodes0(z1) subnetwork_nodes0(z2) p];

%the nonlinear control


k=sum(CC);%
scores=k';

subnetwork_adjacency=CC;
[z1,z2]=find(triu(subnetwork_adjacency)~=0);
z=[z1 z2];

%the nonlinear control

NN=length(subnetwork_adjacency);
[ nc_x,nc_index ] = Nar_Opti_weight_nc( z,NN,scores,lamda  );

%[ nc_x,nc_index ] = weight_nc( z,NNN,scores,lamda  );


NC_nodes{i,1}=subnetwork_nodes0(find(nc_x~=0));
NC_genes{i,1}=gene_list(subnetwork_nodes0(find(nc_x~=0)),1);






end



%******part3:Predict the combinational drugs******

Result_predict_drug=[];

Final_Sample_name_normal=Sample_name_tumor';

for i=1:size(NC_genes,1)
%for i=1:1    

control_node=NC_nodes{i,1};
drug_scores_matrix=Matrix_Drug_target(control_node,:);
%Drug_scores(i,:)=sum(drug_scores_matrix);.*repmat(All_sample_scores(control_node,i),1,size(Matrix_Drug_target,2))
k=sum(drug_scores_matrix);
[value,ind]=sort(k,'descend');
%ind(value==0)=[];value(value==0)=[];
predict_drug=Final_drug_name(ind);
Result_predict_drug{i,1}=Final_Sample_name_normal{i,1};
drug=[];
for j=1:length(value)
     
    drug{j,1}=Final_drug_name(ind(1,j));
   
    drug{j,2}=value(1,j);
end

Result_predict_drug{i,2}=gene_list(control_node);
Result_predict_drug{i,3}=drug;
  

end


%*****Output2:All_Sample_drug,Obtain efficacious drug pairs and corresponding drug targets*******


for i=1:size(NC_nodes,1)
i
control_node=NC_nodes{i,1};
drug_scores_matrix=Matrix_Drug_target(control_node,:);
%Drug_scores(i,:)=sum(drug_scores_matrix);.*repmat(All_sample_scores(control_node,i),1,size(Matrix_Drug_target,2))
k=sum(drug_scores_matrix);
[u,v]=sort(k,'descend');
dk_structure=drug_scores_matrix;

predict_drug= Final_drug_name(v(1,1:10),1);
Personalized_drug{i,1}=predict_drug;


end



genes=gene_list(control_node);
[~,data,~]=xlsread('Drug_design_data.xlsx');
drug_id=data(:,5);


All_drug_design_information=[];
for i=1:size(NC_nodes,1)
% for i=1:10

i

z=All_CC{i,1};

CC=zeros(length(gene_list));

for iii=1:size(z,1)
    CC(z(iii,1),z(iii,2))=z(iii,3);
    CC(z(iii,2),z(iii,1))=z(iii,3);
end


control_node=NC_nodes{i,1};
cand=control_node;
%nodes=intersect(cand,cx);
nodes=cand;
Nodes_CC=CC(nodes,nodes);
cand_personalized_drug=Personalized_drug{i,1};

A=Nodes_CC;
drug_design_information=[];
for j=1:length(cand_personalized_drug)
j
cand=cand_personalized_drug{j,1};
[a,b]=ismember(drug_id,cand);
drug_data=data(find(b~=0),1:4);
cand_drug=unique([drug_data(:,2);drug_data(:,4)]);

drug_data1=data(:,1:4);
[~,x1]=ismember(drug_data1(:,1),genes);
[~,x2]=ismember(drug_data1(:,2),cand_drug);
[~,x3]=ismember(drug_data1(:,3),genes);
[~,x4]=ismember(drug_data1(:,4),cand_drug);
x=[x1 x2;x3 x4];
y=x(:,1).*x(:,2);
x(find(y==0),:)=[];
B=zeros(length(genes),length(cand_drug));
for k=1:size(x,1)
B(x(k,1),x(k,2))=1;
end

drug_design_information{j,1}=cand_drug;
drug_design_information{j,2}=genes;
drug_design_information{j,3}=A;
drug_design_information{j,4}=B;
end
All_drug_design_information{i,1}=drug_design_information;



end




[data]=importdata('COVID_cancer_genes.xlsx');
Breast_genes_set=unique(data.textdata.Sheet1(2:end,1));


All_Sample_drug=[];
%for i=1:size(All_drug_design_information,1)
for i=1:size(NC_nodes,1)
  
i
z=All_CC{i,1};

CC=zeros(length(gene_list));

for iii=1:size(z,1)
    CC(z(iii,1),z(iii,2))=z(iii,3);
    CC(z(iii,2),z(iii,1))=z(iii,3);
end



control_node=NC_nodes{i,1};
drug_ID=Personalized_drug{i,1};

cand=All_drug_design_information{i,1};
Sample_drug=[];
Breast_genes_set=gene_list(cx);

Drug_combinations_information=[];
kk=1;

for j=1:size(cand,1)
    
j
Drug_matrix=cand{j,4};drug=cand{j,1};
All_possible_combination=combntns([1:size(Drug_matrix,2)],2);
Net_drug_similarity=[];Targets=[];Attribute=[];

for k=1:size(All_possible_combination,1)
    
k

cand1=Drug_matrix(:,All_possible_combination(k,1));
cand2=Drug_matrix(:,All_possible_combination(k,2));
Targets1=cand{j,2}(find(cand1~=0));
Targets2=cand{j,2}(find(cand2~=0));


%[scores]=combinational_drugs_scores(gene_list,Targets1,Targets2,Breast_genes_set,CC);
[scores]=DIAMOnD_combinational_drugs_scores(gene_list,Targets1,Targets2,Breast_genes_set,CC);

Drug_combinations_information{kk,1}=drug(All_possible_combination(k,1));
Drug_combinations_information{kk,2}=Targets1;
Drug_combinations_information{kk,3}=drug(All_possible_combination(k,2));
Drug_combinations_information{kk,4}=Targets2;
Drug_combinations_information{kk,5}=drug_ID{j,1};
Drug_combinations_information{kk,6}=scores;
kk=kk+1;

end

end
All_Sample_drug{i,1}=Drug_combinations_information;


end

%**************Output3: Personalzied Side effect*************


load('side_effect_deal_data.mat')

for i=1:size(NC_nodes,1)

    i
    
    %breast_NC_genes=intersect(NC_genes{i,1},CGC_genes);
    breast_NC_genes=NC_genes{i,1};
    %breast_NC_genes=intersect(NC_genes{i,1},cancer_genes);
    [aggrevating_side_effect,inprove_side_effect] = function_sample_side_effect(Filter_BEGene,breast_NC_genes,new_Filter_DT_network,new_Filter_DT_network_pharmacological,drug_names_all,target_names_all,target_names_pharmacological);
    sample_inprove_side_effect(i,1)=inprove_side_effect;    
    sample_aggrevating_side_effect(i,1)=aggrevating_side_effect; 
    
   
end

Personalized_side_effect=[sample_aggrevating_side_effect sample_inprove_side_effect];








end


