path(path,'../../3_usefuldata');
load('g_p_netwrok_2015_3_4&5.mat');

%MaxIterΪһ��NMF������InnerMaxIterΪ�ڲ�����
MaxIter=30;
InnerMaxIter=3;

%�������K
K=300;
%t_circleΪ��ʼ������
t_circle=10;

%V1Ϊlevel4 g-p����
%V2Ϊlevel5 g-p����
V1 = g_p_network4;
V2 = g_p_network5;

[m,n1]=size(V1);
[~,n2]=size(V2);

%����
alphas=[0,0.001,0.01,0.1,1];
gamas=[0]; 
lamtas1=[1];
lamtas2=[10];

% for alpha=alphas
%     for gama = gamas
%         for lamta1=lamtas1
%             for lamta2=lamtas2
%                 for t=1:t_circle
%                     W= rand(m,K);
%                     H1 = rand(K,n1);
%                     H2 = rand(K,n2);
%                     
%                     [L,W_out,H1_out,H2_out] =  CMNMF_L1(K,MaxIter,InnerMaxIter,V1,V2,W,H1,H2,M,alpha,gama,lamta1,lamta2);
%                      
%                          directory='../../5_6_result_2015/CMNMF/CMNMF_L1(';
%                          if(~exist(directory,'dir'))
%                              mkdir(directory);
%                          end
%                        
%                          fn = ['../../5_6_result_2015/CMNMF/CMNMF_L1(/CMNMF_simple_alpha' num2str(alpha)  '_gama&' num2str(gama) '_lamta1&' num2str(lamta1)  '_lamta2&' num2str(lamta2) '_t' num2str(t) '.mat'];
%                          disp([datestr(now) ':  '  fn ]);
%                          save(fn,'L','W_out','H1_out','H2_out');
%                    
%                 end
%             end
%         end
%     end
% end

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
pathway_gene_mgi_id = cell(a*b*c*d,t_circle);

alpha=alphas(1);
gama=gamas(1);

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
                fn = ['../../5_6_result_2015/CMNMF/CMNMF_L1(/CMNMF_simple_alpha' num2str(alpha1)  '_gama&' num2str(gama1) '_lamta1&' num2str(lamta1)  '_lamta2&' num2str(lamta2) '_t' num2str(t) '.mat'];
                if(~exist(fn,'file'))
                    continue;
                end
                load(fn);
                %locationΪĳ�������Ӧ��λ��
                location=(al-1)*b*c*d+(ga-1)*c*d+(j-1)*d+k;
                [pathway_gene_mgi_id{location,t}]= predicted_pathway(W_out,C,T,mgi_id);
                end
            end
        end
    end
end

%�洢���
datetime=fix(clock);
      s='';
      for i=1:6
      s=[s num2str(datetime(i))];
      end
      fn2=[ '../../5_6_result_2015/CMNMF/CMNMF_L1(/predicted_pathway_L1_' s '.mat'];
      save(fn2,'pathway_gene_mgi_id');

