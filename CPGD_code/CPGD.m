function [Result_predict_drug,All_Sample_drug,Personalized_side_effect] = CPGD( expression_tumor_fileName,expression_normal_fileName,lamda )
%PDC outputs the patient-specific drug profiles
%   Input:
%         expression_tumor_fileName:the tumor expression data in cancer
%         expression_normal_fileName:the normal expression data in cancer 
%         lamda:the balance parameter between the scores objective and minimum objective,the bigger the value is, the more preferable we choose the minimum objective.
%               the fefault value is 100.


%  Output:
%         Output1:
%         Result_drug_combinations:individual combinational drugs
%         (DCDB);the first column is the sample name with paired data
%         (normal and tumor) and the second column is the ranked
%         combinational drug name in descend with defined scores.The
%         defined scores are the number of targeted personalized driver
%         genes.clc

%         Output2:
%         Result_drug_targets: individual drug-target genes；the first column is the sample name with paired data (normal and tumor) and the second column is the predicted personalized drug targets.
%         Result_efficacious_Drug_pairs：the 1-4 column is efficacious synergistic drug pairs
%         and the 5 colunm is the related combinational drugs (name in DCDB),the 6 colunm is the corresponding synergistic scores

%         Output3:
%         Result_patients_side_effect：the first column is the sample name and the second colunm denote the aggrevating effect and the third colunm denotes the improving effect;

%************************part1:LOAD BRCA sample data and network data************************
%********************obtain the paired expression data******************
%read the original text of tumor and normal

%The temple data in our study


%*************************tumor****************************

[tumor,~,name_tumor]=importdata(expression_tumor_fileName);
Sample_gene_list=tumor.textdata(2:end,1);
Sample_name_tumor=tumor.textdata(1,2:end);
tumor_data=tumor.data;

%*************************normal****************************

[normal,~,name_normal]=importdata(expression_normal_fileName);
Sample_name_normal=normal.textdata(1,2:end);
normal_data=normal.data;

%************************paired data**************************
%*****************obtain the tumor sample name**************
filter_name_tumor=[];
for i=1:length(Sample_name_tumor)
    
    
    S=Sample_name_tumor{1,i};
    %[a,b]=find(S=='-');
    %filter_name_tumor{i,1}=S(1:b(1,end)-1);
    filter_name_tumor{i,1}=S;
end

%*****************obtain the normal sample name***************
filter_name_normal=[];
for i=1:length(Sample_name_normal)
    
    
    S=Sample_name_normal{1,i};
    %[a,b]=find(S=='-');
    %filter_name_normal{i,1}=S(1:b(1,end)-1);
    filter_name_normal{i,1}=S;
    
end
%************* matching for the tumor and normal data************

[ind,address]=ismember(filter_name_normal,filter_name_tumor);

k=1;Sample_Tumor=[];Sample_Normal=[];
for i=1:length(address)
    
    i
    if address(i,1)~=0  
       Final_Sample_name_normal{k,1}=filter_name_normal{i,1};
        Sample_Normal(:,k)=normal_data(:,i);
        Sample_Tumor(:,k)=tumor_data(:,address(i,1));
        k=k+1;
        
    end
    
end



%The temple data in our study

%*************************tumor****************************

%unpack zip files
unzip('BRCA_normal.zip') 
unzip('BRCA_tumor.zip') 

expression_tumor_fileName = 'BRCA_tumor.txt';
expression_normal_fileName = 'BRCA_normal.txt';

[tumor,~,name_tumor]=importdata(expression_tumor_fileName);
temple_gene_list=tumor.textdata(2:end,1);
Sample_name_tumor=tumor.textdata(1,2:end);
tumor_data=tumor.data;

%*************************normal****************************

[normal,~,name_normal]=importdata(expression_normal_fileName);
Sample_name_normal=normal.textdata(1,2:end);
normal_data=normal.data;

%************************paired data**************************
%*****************obtain the tumor sample name**************
filter_name_tumor=[];
for i=1:length(Sample_name_tumor)
    
    
    S=Sample_name_tumor{1,i};
    %[a,b]=find(S=='-');
    %filter_name_tumor{i,1}=S(1:b(1,end)-1);
    filter_name_tumor{i,1}=S;
end

%*****************obtain the normal sample name***************
filter_name_normal=[];
for i=1:length(Sample_name_normal)
    
    
    S=Sample_name_normal{1,i};
    %[a,b]=find(S=='-');
    %filter_name_normal{i,1}=S(1:b(1,end)-1);
    filter_name_normal{i,1}=S;
    
end
%************* matching for the tumor and normal data************

[ind,address]=ismember(filter_name_normal,filter_name_tumor);

k=1;temple_Tumor=[];temple_Normal=[];
for i=1:length(address)
    
    i
    if address(i,1)~=0  
       Final_Sample_name_normal{k,1}=filter_name_normal{i,1};
        temple_Normal(:,k)=normal_data(:,i);
        temple_Tumor(:,k)=tumor_data(:,address(i,1));
        k=k+1;
        
    end
    
end

%Final data in the following procedure

gene_list=intersect(Sample_gene_list,temple_gene_list);

[Sample_ind,Sample_address]=ismember(gene_list,Sample_gene_list);
Final_Sample_Tumor=Sample_Tumor(Sample_address,:);
Final_Sample_Normal=Sample_Normal(Sample_address,:);


[temple_ind,temple_address]=ismember(gene_list,temple_gene_list);
Final_temple_Tumor=temple_Tumor(temple_address,:);
Final_temple_Normal=temple_Normal(temple_address,:);

Tumor0=[Final_Sample_Tumor Final_temple_Tumor];
Normal0=[Final_Sample_Normal Final_temple_Normal];

%keep the order of rows
[bTumor,mTumor]=unique(Tumor0.','rows');
[mTumor,mTumor]=sort(mTumor);
Tumor=bTumor(mTumor,:).';

[bNormal,mNormal]=unique(Normal0.','rows');
[mNormal,mNormal]=sort(mNormal);
Normal=bNormal(mNormal,:).';

%**********Synthetic lethality Network************************


Network=[];
D=[];
Result=[];
[weight,PPI,~]=xlsread('network_lethal.xlsx');
value=weight(:,3);
[x1,y1]=ismember(PPI(:,1),gene_list);
[x2,y2]=ismember(PPI(:,2),gene_list);

P=unique(PPI);
[xx,yy]=ismember(P,gene_list);

y=y1.*y2;
z=[y1 y2];
z(find(y==0),:)=[];
value(find(y==0),:)=[];
N1=length(gene_list);
[N2,~]=size(z);


%calculate the adjacency matrix of PPI 

Net_adjacent=zeros(N1,N1);
for i=1:N2
    i
         Net_adjacent(z(i,2),z(i,1))=1;  %undirected gene-gene interaction network
         Net_adjacent(z(i,1),z(i,2))=1;
       
end




%***********construct co-mutation network for BRCA mutation data set***********

[weight,PPI,~]=xlsread('BRCA_lethal_net_co_mutation.xlsx');
value=weight(:,end);
[x1,y1]=ismember(PPI(:,1),gene_list);
[x2,y2]=ismember(PPI(:,2),gene_list);

P=unique(PPI);
[xx,yy]=ismember(P,gene_list);

y=y1.*y2;
z=[y1 y2];
z(find(y==0),:)=[];
value(find(y==0),:)=[];
N1=length(gene_list);
[N2,~]=size(z);


%calculate the adjacency matrix of PPI 

New_A_network=zeros(N1,N1);
for i=1:N2
    i
         New_A_network(z(i,2),z(i,1))=value(i,1);  %undirected gene-gene interaction network
         New_A_network(z(i,1),z(i,2))=value(i,1);
       
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
for i=1:size(Final_Sample_Tumor,2)
%for i=1:2
%for i=1:size(new_T,2)
    
    i
        

    %for the i-th sample
    %construct the tumor SSN**********
    
    
    sample_tumor=new_T0(:,i);
    
    
    [R0,P]=SSN(sample_tumor,Ref0);
    
    P(isnan(P))=0;
    P(abs(P)>=0.05)=2;P(abs(P)<0.05)=1;
    P(P==2)=0;
    %construct the normal SSN
    clear  sample_tumor 
    sample_normal=new_N0(:,i);
    [R1,P1]=SSN(sample_normal,Ref0);
    clear  sample_normal 
    
    P1(isnan(P1))=0;
    P1(abs(P1)>=0.05)=2;P1(abs(P1)<0.05)=1;
    P1(P1==2)=0;
    C=P-P1; 
    
    R=abs(log2(abs(R0./R1)));
    R(isnan(R))=0;
    
    %sample_network{i,1}=C;
    D0=abs(C.*R);
    
    
    %sample_network{i,1}=C;

CC=abs(mapminmax(D0,0,1).*Net_adjacent(subnetwork_nodes0,subnetwork_nodes0).*mapminmax(New_A_network(subnetwork_nodes0,subnetwork_nodes0),0,1));
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

NNN=length(subnetwork_adjacency);
%[ nc_x,nc_index ] = Opti_weight_nc( z,NNN,scores,lamda  );
[ nc_x,nc_index ] = weight_nc( z,NNN,scores,lamda  );


NC_nodes{i,1}=subnetwork_nodes0(find(nc_x~=0));
NC_genes{i,1}=gene_list(subnetwork_nodes0(find(nc_x~=0)),1);






end



%******part3:Predict the combinational drugs******

Result_predict_drug=[];
[data]=importdata('validated_drug_breast_cancer.xlsx');
dc_standard=unique(data.textdata.Sheet1(2:end,1));



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
    [c,d]=ismember(Final_drug_name(ind(1,j)),dc_standard);
    drug{j,2}=c;
    drug{j,3}=value(1,j);
end
Result_predict_drug{i,2}=gene_list(control_node);
Result_predict_drug{i,3}=drug;
  

end


%*****Output2:All_Sample_drug,Obtain efficacious drug pairs and corresponding drug targets*******


data=importdata('GO_data.csv');
data(find(data(:,1)>length(gene_list)),:)=[];
data(find(data(:,2)>length(gene_list)),:)=[];
GO_sim_data0=zeros(length(gene_list));
for i=1:size(data,1)
    i
    GO_sim_data0(data(i,1),data(i,2))=data(i,3);
     GO_sim_data0(data(i,2),data(i,1))=data(i,3);
    
end
v=ones(1,length(gene_list));
D=diag(v);
GO_sim_data=GO_sim_data0+D;


%cmap matrix
[~,data,~]=xlsread('Drug_design_data.xlsx');
all_drugs=unique([data(2:2895,2);data(2:2895,4)]);


unzip('cmap_disimilarity.zip') 
u1=importdata('cmap_disimilarity.csv');
cmap_data=u1.data;
name_u1=upper(u1.textdata(2:end,1));

[c1,d1]=ismember(all_drugs,name_u1);
cmap_dissimilarity_data0=zeros(length(d1),length(cmap_data));

for i=1:length(d1)
    
    if d1(i,1)~=0
        i
        cmap_dissimilarity_data0(i,:)=cmap_data(d1(i,1),:);     
        
    end
    
end

cmap_dissimilarity_data=zeros(length(d1));
for i=1:length(d1)
    
    if d1(i,1)~=0     
        cmap_dissimilarity_data(:,i)=cmap_dissimilarity_data0(:,d1(i,1));        
    end
    
 end
cmap_dissimilarity_data(isnan(cmap_dissimilarity_data))=0;


%structure matrix
u2=importdata('structural_similarity.csv');
name_z1=u2.textdata(2:end,2);
name_z2=u2.textdata(2:end,3);

structure_data=u2.data;
[~,z1]=ismember(name_z1,all_drugs);
[~,z2]=ismember(name_z2,all_drugs);


structural_similarity_matrix=[];
for i=1:length(z1)
    
    structural_similarity_matrix(z1(i,1),z2(i,1))=structure_data(i,1);
    structural_similarity_matrix(z2(i,1),z1(i,1))=structure_data(i,1);
    
end
structural_similarity_matrix(isnan(structural_similarity_matrix))=0;




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
    CC(z(iii,1),z(iii,2))=p(iii,1);
    CC(z(iii,2),z(iii,1))=p(iii,1);
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




[data]=importdata('breast_cancer_genes.xlsx');
Breast_genes_set=unique(data.textdata.Sheet1(2:end,1));


All_Sample_drug=[];
%for i=1:size(All_drug_design_information,1)
for i=1:size(NC_nodes,1)
  
i
z=All_CC{i,1};

CC=zeros(length(gene_list));

for iii=1:size(z,1)
    CC(z(iii,1),z(iii,2))=p(iii,1);
    CC(z(iii,2),z(iii,1))=p(iii,1);
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
[scores]=DIAMOnD_combinational_drugs_scores(GO_sim_data,gene_list,Targets1,Targets2,Breast_genes_set,CC);

C1=drug(All_possible_combination(k,1));
C2=drug(All_possible_combination(k,2));
[~,uu1]=ismember(C1,all_drugs);
[~,uu2]=ismember(C2,all_drugs);

ss=structural_similarity_matrix(uu1,uu2);
sg=cmap_dissimilarity_data(uu1,uu2);
v=ss+1-sg;


Drug_combinations_information{kk,1}=drug(All_possible_combination(k,1));
Drug_combinations_information{kk,2}=Targets1;
Drug_combinations_information{kk,3}=drug(All_possible_combination(k,2));
Drug_combinations_information{kk,4}=Targets2;
Drug_combinations_information{kk,5}=drug_ID{j,1};
Drug_combinations_information{kk,6}=scores+v;
kk=kk+1;

end

end
All_Sample_drug{i,1}=Drug_combinations_information;


end

%**************Output3: Personalzied Side effect*************


[~,CGC,~]=xlsread('CGC.xlsx');
CGC_genes=CGC;
NCG_genes=importdata('NCG_name.txt');
cancer_genes=[CGC_genes;NCG_genes];


load('side_effect_deal_data.mat')

for i=1:size(NC_nodes,1)

    i
    
    %breast_NC_genes=intersect(NC_genes{i,1},CGC_genes);
    %breast_NC_genes=NC_genes{i,1};
    breast_NC_genes=intersect(NC_genes{i,1},cancer_genes);
    [aggrevating_side_effect,inprove_side_effect] = function_sample_side_effect(Filter_BEGene,breast_NC_genes,new_Filter_DT_network,new_Filter_DT_network_pharmacological,drug_names_all,target_names_all,target_names_pharmacological);
    sample_inprove_side_effect(i,1)=inprove_side_effect;    
    sample_aggrevating_side_effect(i,1)=aggrevating_side_effect; 
    
   
end

Personalized_side_effect=[sample_aggrevating_side_effect sample_inprove_side_effect];



end


