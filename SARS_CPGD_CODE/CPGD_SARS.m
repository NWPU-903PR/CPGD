function [CPGD_rank_drugs] = CPGD_SARS(expression_tumor_fileName,expression_normal_fileName,lamda)
%CPGD outputs the drug combinations of SARS-COV2 
%   Input:
%         expression_tumor_fileName:the tumor expression data in cancer
%         expression_normal_fileName:the normal expression data in cancer 
%         lamda:the balance parameter between the scores objective and minimum objective,the bigger the value is, the more preferable we choose the minimum objective.
%               the fefault value is 0.01.


%  Output:
%         Output:CPGD_rank_drugs
%         CPGD_rank_drugs:The ranking pairwise drug combinationas

All_parameter_Result_predict_drug=[];
All_parameter_All_Sample_drug=[];
All_parameter_Personalized_side_effect=[];
All_parameter_Pattern_drugs=[];
All_parameter_Frequency_drugs=[];


for kkk=1:5

kkk
lamda=(10^kkk)*0.001; %the default value is 0.01 in BRCA cancer data set.

tic

[Result_predict_drug,All_Sample_drug,Personalized_side_effect] = CPGD( expression_tumor_fileName,expression_normal_fileName,lamda );

toc

All_parameter_Result_predict_drug{kkk,1}=Result_predict_drug;
All_parameter_All_Sample_drug{kkk,1}=All_Sample_drug;
All_parameter_Personalized_side_effect{kkk,1}=Personalized_side_effect;





[~,data,~]=xlsread('/home/disk1/guoweifeng/old_server_files/BRCA_data_deal/Drug_design_data.xlsx');
All_drugs=unique(data(2:2895,[2,4]));

All_Sample_drug4=All_Sample_drug;
Patern_All_Sample_drug=[];All_patern_drug=[];
All_zz=[];All_value=[];

for i=1:size(All_Sample_drug4,1)
    
    
    i
 
    cand=All_Sample_drug4{i,1};
    value=cell2mat(cand(:,6));
   value(isnan(value))=0;
    
     new_cand=[];
     for k=1:size(cand,1)
        
        new_cand{k,1}=char(cand{k,1}{1,1});
        new_cand{k,2}=char(cand{k,3}{1,1});
     end
    
    [~,z1]=ismember(new_cand(:,1),All_drugs);
    [~,z2]=ismember(new_cand(:,2),All_drugs);
    zz=[z1 z2];
     
    Sample_drugs_id{i,1}=zz;
    Sample_value{i,1}=value;
    
    All_zz=[All_zz;zz];
    All_value=[All_value;value];
    
    Patern_All_Sample_drug{i,1}=cand;
    
    
    
end   
  


u=unique(All_zz,'rows');

%calculate the mean scores of each two drug combinations

N1=max(max(u));
A=zeros(N1);

for i=1:size(Sample_drugs_id,1)
    i
    cand1=Sample_drugs_id{i,1};
    cand2=Sample_value{i,1};
    A0=zeros(N1);
    for j=1:size(cand1,1)
        A0(cand1(j,1),cand1(j,2))=cand2(j,1);
         A0(cand1(j,2),cand1(j,1))=cand2(j,1);
    end
    A=A+A0;
end
A1=A/size(Sample_drugs_id,1);


Frequency_drugs4=[];
for i=1:size(u,1)
    
    i

    a=A1(u(i,1),u(i,2));
    Frequency_drugs4(i,1)=a;
    
end


Pattern_drugs=All_drugs(u);


All_parameter_Pattern_drugs{kkk,1}=Pattern_drugs;
All_parameter_Frequency_drugs{kkk,1}=Frequency_drugs4;





end



all=[];
for i=1:5

cand1=All_parameter_Frequency_drugs{i,1};
cand2=All_parameter_Pattern_drugs{i,1};
Final_All_parameter_Pattern_drugs{i,1}=cand2(find(cand1~=0),:);
Final_All_parameter_Frequency_drugs{i,1}=cand1(find(cand1~=0),:);

all=[all;unique(cand2)];



end


All_x=[];
for i=1:5

cand1=Final_All_parameter_Pattern_drugs{i,1};



[~,x1]=ismember(cand1(:,1),all);
[~,x2]=ismember(cand1(:,2),all);

x=[x1 x2];
All_parameter_id{i,1}=x;

All_x=[All_x;x];

end

All_x=unique(All_x,'rows');
All_drugs=all(All_x);



for i=1:5

cand=All_parameter_id{i,1};
cand0=Final_All_parameter_Frequency_drugs{i,1};


cand_id=[];
for j=1:size(cand,1)

     

[u1,x1]=ismember(All_x(:,1),cand(j,1));
[u2,x2]=ismember(All_x(:,2),cand(j,2));

xx=x1.*x2;

cand_id(j,1)=find(xx~=0);




end

All_parameter_num{i,1}=[cand_id cand0];




end



[ CPGD_rank_drugs ] = condorcet( All_parameter_num,All_drugs);


end

