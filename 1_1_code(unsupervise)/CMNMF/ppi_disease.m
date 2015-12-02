load('mmu_pathway_mgi_id_true.mat');
load('ppi_network_4588.mat');
load('selected_pathway.mat');
% for i=1:length(selected_pathway)
% type2=mmu_pathway_mgi_id_true(258,:);
% [type2_gene,type2_ppi] = ppi(type2);
% end
load('cluster_result.mat');
cluster_gene = unique(b);
[~,~,ib]=intersect(unique(b),entrez_id);
a = nnz(ppi_network(ib,ib));
% for i=1:229
%     c=intersect(b(:,i),entrez_id);
%     intersect_ppi(i) = length(c);
% end
% [order,index]=sort(intersect_ppi,'descend');
% select_cluster=cell(229);
% for i=1:229
%     [gene,~,ib]=intersect(b(:,index(i)),entrez_id);
%     select_cluster{i}=ppi_network(ib,ib);
%     x(i)=nnz(select_cluster{i});
% end