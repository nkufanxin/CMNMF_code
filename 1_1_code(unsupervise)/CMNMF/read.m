clc;
clear;
fid = fopen( 'symbols_1274.txt', 'r' );
C = textscan(fid, '%s'); % һ������������
%data = cell2mat(C); % ת������ͨ����plot(data);hold on;