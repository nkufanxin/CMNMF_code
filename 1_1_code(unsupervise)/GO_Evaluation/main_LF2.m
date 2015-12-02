path(path,'../../3_usefuldata');
load('g_p_netwrok_2015_3_4&5.mat');

%max_pathway_gene_numΪɸѡ��G0����ֵ��pathway�ĳ������ֵ��
%TΪɸѡ��G0_NolessThan_T����ֵ��ÿ��pathway������mgi_id�ཻ��pathway����С���ȣ�
max_pathway_gene_num=200;
T=5;

build_G0(max_pathway_gene_num,T);
load('G0_data.mat');
load('mgi_id');



%MaxIterΪһ��NMF������InnerMaxIterΪ�ڲ�����
MaxIter=30;
InnerMaxIter=3;

%�������K
K=150;
%t_circleΪ��ʼ������
t_circle=20;

%V1Ϊlevel4 g-p����
%V2Ϊlevel5 g-p����
V1 = g_p_network4;
V2 = g_p_network5;

[m,n1]=size(V1);
[~,n2]=size(V2);

%����
alphas=[0,1];
gamas=[0,0.9];
%lamtas1 ̫��᲻����
lamtas1=[1];
lamtas2=[1];

load('go_mgi_network.mat');
load('mgi_truth_ground.mat');
load('hierarchy_g_p_network.mat','mgi_id');
load('mmu_pathway_data.mat');

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

a = length(alphas);
b = length(gamas);
c = length(lamtas1);
d = length(lamtas2);

% jaccards=zeros(a*b*c*d,t_circle);
% nonEmptyClusterNum = zeros(a*b*c*d,t_circle);

%pathway_gene_mgi_idΪԤ�����gene��pathway��ϵ
RD = zeros(a*b*c*d,t_circle);
F = zeros(a*b*c*d,t_circle);
Z_filter = cell(a*b*c*d,t_circle);
pathway_gene_mgi_id = cell(a*b*c*d,t_circle);
rat = zeros(a*b*c*d,t_circle+4);
%TΪzscore��ֵ
T=3;
for al=1:a
    alpha1=alphas(al);
    for ga=1:b
        gama1=gamas(ga); 
        for j=1:c
            lamta1=lamtas1(j);
            for k=1:d
                lamta2=lamtas2(k);
                for t=1:t_circle
                fn = ['../../5_6_result_2015/CMNMF/CMNMF_LF/CMNMF_simple_alpha' num2str(alpha1)  '_gama&' num2str(gama1) '_lamta1&' num2str(lamta1)  '_lamta2&' num2str(lamta2) '_t' num2str(t) '.mat'];
                if(~exist(fn,'file'))
                    continue;
                end
                load(fn);
                %locationΪĳ�������Ӧ��λ��
                location=(al-1)*b*c*d+(ga-1)*c*d+(j-1)*d+k;
                [Z_filter{location,t},pathway_gene_mgi_id{location,t}]= predicted_pathway(W_out,C,T,mgi_id);
                [RD(location,t),F(location,t)]=rand_index(Z_filter{location,t},G0_NoLessThan_T);
                rat(location,1) = alpha1;
                rat(location,2) = gama1;
                rat(location,3) = lamta1;
                rat(location,4) = lamta2;
                rat(location,t+4) = ratio(pathway_gene_mgi_id{location,t});
                end
            end
        end
    end
end

%save('ratio.mat','F','ratio');
%�洢���
% datetime=fix(clock);
%       s='';
%       for i=1:6
%       s=[s num2str(datetime(i))];
%       end
%       fn2=[ '../../5_6_result_2015/CMNMF/CMNMF_LF/predicted_pathway_LF_' s '.mat'];
%       save(fn2,'pathway_gene_mgi_id','RD','F');

