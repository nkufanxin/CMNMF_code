function [ranked_pathways_for_specific_gene] = rankgene(Th,G,G0,target)

%      G=G(:,SelectedPathway_G0_NoLessThan_T_norm);
%      G0 = G0(:,SelectedPathway_G0_NoLessThan_T_norm);

      %index为G每一行降序的索引值
      [~,index]=sort(G,2,'descend');
      [~,G0_col] = size(G0);
      
      %preditcted_pathways_idxes为预测出的基因所属的pathway
      %predicted_pathways_idxes = index(target,1:Th);
      %predicted = zeros(1,G0_col);
      %predicted(1,predicted_pathways_idxes) = 1;
      
%       %C为预测结果和真实结果索引的交集
%       [C,ia,ib] = intersect(predicted_pathways_idxes,find(G0(target,:)==1));
        
%       
%       %ranked_pathways_for_specific_gene为输出的交集向量
        ranked_pathways_for_specific_gene = zeros(1,G0_col);
%       ranked_pathways_for_specific_gene(1,C) = 1; 
        for i = 1:G0_col
            if(G0(target,index(target,i))==1);
            ranked_pathways_for_specific_gene(i) = 1;
            end
        end

        
       
       
end

