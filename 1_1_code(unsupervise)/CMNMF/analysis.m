load('symbol_1274.mat');
load('mgi_id.mat');
a=importdata('cluster-genes.xlsx');
b = a.data.Sheet2;

selected = [6,11,29,39,146,167,173];
result=cell(length(selected));
for i=1:length(selected)
    k = b(:,selected(i));
    [~,~,ib]=intersect(k,mgi_id);
    result{i} = C{1,1}(ib);
end

