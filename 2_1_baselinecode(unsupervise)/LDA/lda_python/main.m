path(path,'../../3_usefuldata');
load('g_p_netwrok_2015_3_4&5.mat');

%max_pathway_gene_numΪɸѡ��G0����ֵ��pathway�ĳ������ֵ��
%TΪɸѡ��G0_NolessThan_T����ֵ��ÿ��pathway������mgi_id�ཻ��pathway����С���ȣ�
max_pathway_gene_num=200;
T=2;

build_G0(max_pathway_gene_num,T);
load('G0_data.mat');
load('mgi_id');
load('go_mgi_network.mat');
load('mmu_pathway_data.mat');
load('gt_network229.mat');


% %G0Ϊgene_mgi_id��pathway�Ķ�Ӧ��ϵ
% A=mmu_pathway_mgi_id;
% [x,~]=size(A);
% 
% G0=zeros(length(mgi_id),x);
% 
% for i=1:x
%     [~,~,ib] = intersect(A(i,:),mgi_id);
%     G0(ib,i) = 1;
% end     

% jaccards=zeros(a*b*c*d,t_circle);
% nonEmptyClusterNum = zeros(a*b*c*d,t_circle);

% %pathway_gene_mgi_idΪԤ�����gene��pathway��ϵ
% RD = zeros(a*b*c*d,t_circle);
% F = zeros(a*b*c*d,t_circle);
% Z_filter = cell(a*b*c*d,t_circle);
% pathway_gene_mgi_id = cell(a*b*c*d,t_circle);
% Precision=zeros(a*b*c*d,t_circle);
% Recall=zeros(a*b*c*d,t_circle);
% jaccard=zeros(a*b*c*d,t_circle);

alpha_ratio = 1;
%TΪzscore��ֵ
T=3;
Z_filter = cell(10);
pathway_gene_mgi_id = cell(10);

% for t = 1:10
% a=['gt_network' num2str(t)];
% b=load('gt_network_real.mat',a);
% c=struct2cell(b);
% [Z_filter{t},pathway_gene_mgi_id{t}]= predicted_pathway(c{1},T,mgi_id);
% [RD(t),F(t),Precision(t),Recall(t),jaccard(t)]=rand_index(Z_filter{t},G0_NoLessThan_T,alpha_ratio);
% end

% for t = 1:10
% a=['gt_network' num2str(t)];
% b=load('gt_network_100.mat',a);
% c=struct2cell(b);
% [Z_filter{t},pathway_gene_mgi_id{t}]= predicted_pathway(c{1},T,mgi_id);
% [RD(t),F(t),Precision(t),Recall(t),jaccard(t)]=rand_index(Z_filter{t},G0_NoLessThan_T,alpha_ratio);
% end

for t = 1:1
[Z_filter{t},pathway_gene_mgi_id{t}]= predicted_pathway(gt_network1,T,mgi_id);
[RD(t),F(t),Precision(t),Recall(t),jaccard(t)]=rand_index(Z_filter{t},G0_NoLessThan_T,alpha_ratio);
end




