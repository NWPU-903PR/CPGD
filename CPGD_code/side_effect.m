function [B_effect_value]=side_effect(drug1,drug2,sample_DT,drug_names_all,Filter_BEGene,Filter_genes_pharmacological)
%function:obtain the side effect value for drug1 and drug2
%Input:
%      Drug:drug1 and drug2
%      Network:DT_network,drug_names_all,target_names_all,target_names_pharmacological
%Output:
%      B_effect_value
Filter_DT_network=sample_DT;
[x1,y1]=ismember(drug1,drug_names_all);


[x2,y2]=ismember(drug2,drug_names_all);
v1=Filter_DT_network(:,y1);
v2=Filter_DT_network(:,y2);

v=v1.*v2;
ss=length(setdiff(Filter_BEGene(find(v1~=0)),Filter_genes_pharmacological));

ind=length(setdiff(Filter_BEGene(find(v~=0)),Filter_genes_pharmacological));

%ind=length(find(v~=0))-length(intersect(Filter_BEGene(find(v~=0)),Filter_genes_pharmacological));

B_effect_value=0;

vv=v1+v2;
e_vv=vv.*v;

if ind~=0
    
    s=ind;
    s_coherent=length(setdiff(Filter_BEGene(find(e_vv~=0)),Filter_genes_pharmacological));
    s_incoherent=s-s_coherent;
    B_effect_value=(s_incoherent-s_coherent)/ss;
    
end





end

