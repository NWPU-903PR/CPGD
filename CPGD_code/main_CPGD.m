clc
clear
%   $Id: main_CPGD.m Created at 2020-06-22 16:25:22 $
%   by Weifeng Guo, Northwestern Polytechtical University, China
%   Copyright (c) 2014-2019 by Key Laboratory of Information Fusion Technology of Ministry of Education in Northwestern Polytechnical University,
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%**************Part 1:Input the BRCA samples information ****

%***install gurubi software**************
%cd '/opt/gurobi810/linux64/matlab'
%gurobi_setup
%cd '/home/disk1/guoweifeng/CPGD_code'


expression_tumor_fileName = 'Example_tumor.txt';
expression_normal_fileName = 'Example_normal.txt';

%%**************Part 2:PDC outputs the predicted combinational drugs****

%the balance parameter between the scores objective and minimum objective,the bigger the value is, the more preferable we choose the minimum objective.
lamda=0.01; %the default value is 0.01 in BRCA cancer data set.

tic

[Result_predict_drug,All_Sample_drug,Personalized_side_effect] = CPGD( expression_tumor_fileName,expression_normal_fileName,lamda );

toc

%%**************Part 3:save the result****

save('CPGD_BRCA_results.mat','Result_predict_drug','All_Sample_drug','Personalized_side_effect')
