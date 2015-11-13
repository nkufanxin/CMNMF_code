
path(path,'../../3_usefuldata');
load('g_p_netwrok_2015_3_4&5.mat');
load('similarity.mat');

%max_pathway_gene_numΪɸѡ��G0����ֵ��pathway�ĳ������ֵ��
%TΪɸѡ��G0_NolessThan_T����ֵ��ÿ��pathway������mgi_id�ཻ��pathway����С���ȣ�
max_pathway_gene_num=200;
T=5;
iter = 300;
build_G0(max_pathway_gene_num,T);
load('G0_data.mat');
load('mgi_id');

%t_circleΪ��ʼ������
t_circle=5;
          
 %roc_n_setΪ����rocʱ�Ĵ���ֵn
 roc_n_set = 168;
 %roc_geneΪ��Ž����������ÿһ��Ԫ�ش���ȥ����i������ѵ�����auc���ֵ
 %roc_gene=zeros(length(roc_n_set),length(SelectedGene));
 roc_gene=zeros(length(SelectedGene),t_circle);
 %ÿһ��roc_n_set��Ӧһ��roc cell��ÿ��cellΪһ�����󣬴洢��ͬ�����£���ͬ��ʼ�����aucֵ
 roc = zeros(1,t_circle);
 for t=1:t_circle                  
      for i=1:length(SelectedGene)
            %targetΪ��Ҫȥ���Ļ�������Ӧ��mgi_idλ��
            target=find(mgi_id==SelectedGene(i));
            
            [W_out,rmse,f_u] = lpa(simi,G0_NoLessThan_T_norm,target,iter);

            ranked_pathways_for_specific_gene=rankgene(W_out,G0_NoLessThan_T,target);


            roc_gene(i,t)=ROC(ranked_pathways_for_specific_gene,roc_n_set);  

      end
       roc(t)=mean(roc_gene(:,t));  
 end

datetime=fix(clock);
      s='';
      for i=1:6
      s=[s num2str(datetime(i))];
      end
      fn2=[ '../../5_8_result_2015/LP/ROC_' s '.mat'];
      save(fn2,'roc');

