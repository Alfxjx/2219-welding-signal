% disp('����ƽ����ʱƵ��')
% ��������
% importData
% disp(cond);
dadi = importdata(cond,' ',0);
% ����ȡ��λ�ú�ȡ���ĳ���
% start1 = input('input start dadi: ');
start1 = pos;
step1 = 40960;
% start2 = input('input start gm: ');
% step2 = input('input step gm: ');
% ��ʾһ�´�������С��
dadi_volt = dadi(start1:(start1+step1),1);
dadi_curr = dadi(start1:(start1+step1),2);
volt = dadi_volt;
current = dadi_curr;

% emd�ֽ�
imf_curr = eemd(current,0,1);
% HHT�任
[hs,f,t] = hht(imf_curr);

% ����ʱƵ��
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

disp(['start: ',num2str(start1),' end: ', num2str(start1+step1),', ʱƵ��Ϊ: ',num2str(res_en)]);


