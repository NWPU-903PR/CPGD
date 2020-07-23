# Instruction of CPGD package
This package includes Matlab scripts and several datasets for demo of CPGD approach:
(a)	main_CPGD.m is a Matlab function for the routine of experimental analysis. CPGD aims to screen the control role of personalized drug combinations (for breast cancer patients) including the prioritization of personalized anti-cancer drug combinations, identification of synergistic paiwise drug combinations and evaluation of drug side effect on personalized drug targets.

(b)  main_CPGD.m is the main script to call CPGD by supplying following parameters:
    (1)	expression_tumor_fileName: the directory locating of the gene expression data as the input data.
    (2)	expression_normal_fileName: the directory locating of the copy number variations data as the input data.

(c) Algorithm_CPGD directory includes Matlab scripts for each step of CPGD analysis, and called in main_CPGD.m

(d) The input datasets include:
(1) tumor.txt: the tumor expression data in cancer.
(2) normal.txt: the normal expression data in cancer.
Note: Our CPGD outputs the information of samples with paired data in the both two files.

(e) The analysis results are saved in directory pointed by fileName: The variable “Result_drug_combinations”, “Result_efficacious_Drug_pairs”  and “Result_patients_side_effect” are the output of our CPGD, indicting the predicted individual combinational drugs and the target genes. 
(1) For “Result_drug_combinations”, the first column is the sample name with paired data (normal and tumor) and the second column is the ranked combinational drug name in descend with defined scores. The scores are the number of targeted personalized driver genes.
(2) For “Result_efficacious_Drug_pairs”, the 1-4 column is efficacious synergistic drug pairs including Drug A, targets of drug A, Drug B, targets of drug B, and the 5 colunm is the related combinational drugs (name in DCDB), the 6 colunm is the corresponding synergistic scores

(3) For “Result_patients_side_effect”, the row denote the patients and the first colunm denote the aggrevating effect and the second colunm denotes the improving effect;


(f) As a demo, users can directly run main_CPGD.m in Matlab. We choose a simple data with 5 breast invasive carcinoma (BRCA) patients as a test case in our demo. This package has been tested in different computer environments as: Window 7 or above; Matlab 2014 or above.

(g) When users analyzed yourself new data, please:
   (1) Prepare input datasets as introduced in (d).
   (2) Clear the previous results.
   (3) Set parameters in main_CPGD.m as introduced in (b).
   (4) Run main_CPGD.m.
   (5) Suggest that the users add all fille in our folders to your folder.

%   $Id: main_CPGD.m Created at 2020-06-22 16:25:22 $ 

%   $Copyright (c) 2014-2020 by Key Laboratory of Information Fusion Technology of Ministry of Education in Northwestern Polytechnical University,and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science$; 
%   $If any problem,pleasse contact shaonianweifeng@126.com for help. $
