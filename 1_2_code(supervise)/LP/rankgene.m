function [ranked_pathways_for_specific_gene] = rankgene(G,G0,target)


       ranked_pathways_for_specific_gene = ((G0(target,:)-G(target,:))==0);
       
       
end

