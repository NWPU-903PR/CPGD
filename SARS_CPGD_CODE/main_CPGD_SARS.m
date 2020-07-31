clc
clear
%   $Id: main_CPGD.m Created at 2020-06-22 16:25:22 $
%   by Weifeng Guo, Northwestern Polytechtical University, China
%   Copyright (c) 2014-2019 by Key Laboratory of Information Fusion Technology of Ministry of Education in Northwestern Polytechnical University,
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%**************Part 1:Input the samples information ****

%***install gurubi software**************
%This step should be revised according to the users' enviroments
cd '/opt/gurobi810/linux64/matlab'
gurobi_setup
cd '/home/disk1/guoweifeng/new_server_files/SARS_CODE'


expression_tumor_fileName = 'SARS_TUMOR_F.xlsx';
expression_normal_fileName = 'SARS_NORMAL_F.xlsx';

%%**************Part 2:PDC outputs the predicted combinational drugs****

[CPGD_rank_drugs] = CPGD_SARS(expression_tumor_fileName,expression_normal_fileName,lamda);

%%**************Part 3:save the result****

save('CPGD_results.mat')
