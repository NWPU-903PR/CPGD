function [aggrevating_side_effect,inprove_side_effect] = function_sample_side_effect(Filter_BEGene,breast_NC_genes,new_Filter_DT_network,new_Filter_DT_network_pharmacological,drug_names_pharmacological,target_names_all,target_names_pharmacological)
    
    drug_names_all=drug_names_pharmacological;
    [x1,y1]=ismember(Filter_BEGene,breast_NC_genes);
    sample_DT=new_Filter_DT_network(find(x1~=0),:);
    
    k0=sum(sample_DT);
    [~,effect_drug]=find(k0~=0);
    
    B=[];
    
    
    for j=1:size(effect_drug,2)
        
        
        for k=1:size(effect_drug,2)
         if k~=j   
        drug1=drug_names_all(1,effect_drug(1,j));
        drug2=drug_names_all(1,effect_drug(1,k));
        v=[effect_drug(1,j);effect_drug(1,k)];
        s_target_names_pharmacological=target_names_pharmacological(1,find(sum(new_Filter_DT_network_pharmacological(:,v),2)==2));
        %s_target_names_pharmacological=[];
        [B_effect_value]=side_effect(drug1,drug2,sample_DT,drug_names_all,target_names_all,s_target_names_pharmacological);
        B(j,k)=B_effect_value;
        
         end
        end
        
    end
    [I1,J1]=find(B<0);    
    aggrevating_side_effect=length(I1);
    
    [I2,J2]=find(B>0);
    inprove_side_effect=length(I2);

end

