[heart_scale_label, heart_scale_inst] = libsvmread('gpp.txt');

max_pathway_gene_num=200;
T=2;
randpick_gene_num=10;

build_G0(max_pathway_gene_num,T,randpick_gene_num);
load('G0_data.mat');
load('mgi_id');

gene_label = load('gene_map.txt');
gene_label = gene_label+1;
heart_scale_label = heart_scale_label+1;
t=1;
testing_label_vector = zeros(length(heart_scale_label),1);
%W0为去掉一个SelectedGene后的的Gene-Pathway矩阵
for i = 1:length(SelectedGene)
    target = find(mgi_id==SelectedGene(i));
    target_label = find(gene_label==target);                   
    heart_scale_label(target_label) = 0;
    model = libsvmtrain(heart_scale_label, heart_scale_inst,'4 2 -c 5');
    [predict,~,~] = libsvmpredict(heart_scale_label, heart_scale_inst, model);
%     directory='../../5_8_result_2015/SVM';
%         if(~exist(directory,'dir'))
%             mkdir(directory);
%         end                  
%     fn = ['../../5_8_result_2015/SVM_t' num2str(t) '_i' num2str(i) '.mat'];
%     disp([datestr(now) ':  '  fn ]);
%     save(fn,'L','W_out','H1_out','H2_out');
end

%[predicted_label, accuracy, decision_values] = svmpredict(testing_label_vector,testing_instance_matrix,model);