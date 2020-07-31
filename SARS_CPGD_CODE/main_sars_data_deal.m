
expression_tumor_fileName = 'SARS_TUMOR_F.xlsx';
expression_normal_fileName = 'SARS_NORMAL_F.xlsx';


[tumor,name_tumor,~]=xlsread(expression_tumor_fileName);

genes=name_tumor(2:end,1);




Sample_name_tumor=name_tumor(1,2:end);
tumor_data=tumor(:,2:end);

%*************************normal****************************

[normal,name_normal,~]=xlsread(expression_normal_fileName);
Sample_name_normal=name_normal(1,2:end);
normal_data=normal(:,2:end);


data=[tumor_data normal_data];
samples=[Sample_name_tumor Sample_name_normal]';
load('data_input.mat', 'normals')

normals(1:length(Sample_name_tumor))=0;
normals(length(Sample_name_tumor)+1:length(Sample_name_tumor)+length(Sample_name_normal))=1;

