path(path,'../../3_usefuldata');
load('g_p_netwrok_2015_3_4&5.mat');

%max_pathway_gene_numΪɸѡ��G0����ֵ��pathway�ĳ������ֵ��
%TΪɸѡ��G0_NolessThan_T����ֵ��ÿ��pathway������mgi_id�ཻ��pathway����С���ȣ�
max_pathway_gene_num=200;
T=5;

build_G0(max_pathway_gene_num,T);
load('G0_data.mat');
load('mgi_id');

%MaxIterΪ��������
%U=Inner_MaxIterΪ�ڲ���������
MaxIter=30;
Inner_MaxIter=2;

[~,K]=size(G0_NoLessThan_T_norm);

%t_circleΪ��ʼ������
t_circle=2;

%V1Ϊlevel4 G-P����
%V2Ϊlevel5 G-P����
V1 = g_p_network4;
V2 = g_p_network5;



W0 = G0_NoLessThan_T_norm;


[m,n1]=size(V1);
[~,n2]=size(V2);

%alphas,beta,gamas,lamtas1,lamtas2Ϊ����
alphas=[1];
gamas=[0.01];
beta=1;
%beta=1-alphas-gamas;
lamtas1=[1];
lamtas2=[1];

L = cell(length(alphas)*length(gamas),t_circle);

[~,col]=size(G0_NoLessThan_T_norm);

for alpha=alphas
    for gama = gamas
          for lamta1=lamtas1
            for lamta2=lamtas2
                for t=1:t_circle
                    for i=1:length(SelectedGene)
                    W= rand(m,K);
                    H1 = rand(K,n1);
                    H2 = rand(K,n2);
                    
                    W0=G0_NoLessThan_T_norm;
                    [~,col]=size(W0);
                   
                    %W0Ϊȥ��һ��SelectedGene��ĵ�Gene-Pathway����
                    W0(mgi_id==SelectedGene(i),SelectedPathway_G0_NoLessThan_T_norm)=0;                   
                    [L{i,t_circle},W_out,H1_out,H2_out] =  SCMNMF_many_L3( MaxIter,Inner_MaxIter,V1,V2,W,W0,H1,H2,M,alpha,beta,gama,lamta1,lamta2);
                    directory='../../5_8_result_2015/SCMNMF/SCMNMF_many_L3';
                        if(~exist(directory,'dir'))
                            mkdir(directory);
                        end
                       
                        fn = ['../../5_8_result_2015/SCMNMF/SCMNMF_many_L3/SCMNMF_simple_alpha' num2str(alpha) '_beta&' num2str(beta) '_gama&' num2str(gama) '_lamta1&' num2str(lamta1)  '_lamta2&' num2str(lamta2) '_t' num2str(t) '_i' num2str(i) '.mat'];
                        disp([datestr(now) ':  '  fn ]);
                        save(fn,'L','W_out','H1_out','H2_out');
                    end
                end
            end
        end
    end
end

 a = length(alphas);
 b = length(gamas);
 c = length(lamtas1);
 d = length(lamtas2);
 
 %Th Ԥ��ʱ��ȡǰTh��pathwayΪ����
 Th=25;
 %roc_n_setΪ����rocʱ�Ĵ���ֵn
 roc_n_set = 168;
 %roc_geneΪ��Ž����������ÿһ��Ԫ�ش���ȥ����i������ѵ�����auc���ֵ
 %roc_gene=zeros(length(roc_n_set),length(SelectedGene));
 roc_gene=zeros(length(SelectedGene),t_circle);
 %ÿһ��roc_n_set��Ӧһ��roc cell��ÿ��cellΪһ�����󣬴洢��ͬ�����£���ͬ��ʼ�����aucֵ
  
  alpha=alphas(1);
  gama=gamas(1);

 for al=1:a
     alpha1=alphas(al);
     for ga=1:b
         gama1=gamas(ga); 
        for j=1:1
            lamta1=lamtas1(j);
            for k=1:1
                lamta2=lamtas2(k);
                for t=1:t_circle                  
                    for i=1:length(SelectedGene)
                           
                            fn = ['../../5_8_result_2015/SCMNMF/SCMNMF_many_L3/SCMNMF_simple_alpha' num2str(alpha1) '_beta&' num2str(beta) '_gama&' num2str(gama1) '_lamta1&' num2str(lamta1)  '_lamta2&' num2str(lamta2) '_t' num2str(t) '_i' num2str(i) '.mat'];
                            if(~exist(fn,'file'))
                                continue;
                            end

                             load(fn);
                             %location Ϊÿһ����������ŵ�λ��
                             location=(al-1)*b+ga;

                             %targetΪ��Ҫȥ���Ļ�������Ӧ��mgi_idλ��
                             target=find(mgi_id==SelectedGene(i));

                             ranked_pathways_for_specific_gene=rankgene(Th,W_out,G0_NoLessThan_T,target);

                             roc_gene(i,t)=ROC(ranked_pathways_for_specific_gene,roc_n_set);  

                    end
                    roc(location,t)=mean(roc_gene(:,t));  
                end
            end
        end
    end
end

datetime=fix(clock);
      s='';
      for i=1:6
      s=[s num2str(datetime(i))];
      end
      fn2=[ '../../5_8_result_2015/SCMNMF/SCMNMF_many_L1/ROC_L3_' s '.mat'];
      save(fn2,'roc');

