function [ranked_pathways_for_specific_gene] = rankgene(Th,G,G0,target)

%      G=G(:,SelectedPathway_G0_NoLessThan_T_norm);
%      G0 = G0(:,SelectedPathway_G0_NoLessThan_T_norm);

      %indexΪGÿһ�н��������ֵ
      [~,index]=sort(G,2,'descend');
      [~,G0_col] = size(G0);
      
      %preditcted_pathways_idxesΪԤ����Ļ���������pathway
      %predicted_pathways_idxes = index(target,1:Th);
      %predicted = zeros(1,G0_col);
      %predicted(1,predicted_pathways_idxes) = 1;
      
%       %CΪԤ��������ʵ��������Ľ���
%       [C,ia,ib] = intersect(predicted_pathways_idxes,find(G0(target,:)==1));
        
%       
%       %ranked_pathways_for_specific_geneΪ����Ľ�������
        ranked_pathways_for_specific_gene = zeros(1,G0_col);
%       ranked_pathways_for_specific_gene(1,C) = 1; 
        for i = 1:G0_col
            if(G0(target,index(target,i))==1);
            ranked_pathways_for_specific_gene(i) = 1;
            end
        end

        
       
       
end

