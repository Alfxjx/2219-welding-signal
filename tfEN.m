% ʱƵ�ؼ���
% ����HHT�� Ȼ��Զ�άƽ����зֿ飬
% ���ŶԷֿ���������й�һ��
% ������Ϣ�׵�ʱƵ��

% ��άƽ�� hs
cell = mat2cell(hs,[10,10,10,10,10,10,10,10,10,11],[4096,4096,4096,4096,4096,4096,4096,4096,4096,4097]);
energy_of_each = cellfun(@(x) sum(full(x(:))), cell);
% ���ŶԷֿ���������й�һ��
energy_of_all = sum(energy_of_each(:));
en_q = energy_of_each/energy_of_all;
% ʱƵ��
res_en = 0;
size_en = size(energy_of_each);

for i=1:size_en(1)
    for j=1:size_en(2)
        single_en = -1*en_q(i,j)* log(en_q(i,j));
        res_en = res_en+single_en;
    end
end
