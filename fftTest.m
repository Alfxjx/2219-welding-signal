% ��������
data = importdata('1.txt',',',0);
% ������ѹ

% ȡһС�β���һ�� ��Ҫ��Ӧ
current = data(23000:36000,2);
volt = data(23000:36000,1);

% Ԥ�������ݣ�����ѹ/�����ź�ת���ɷǸ��ź�

minVIndex = find(volt==min(volt));

out = volt(minVIndex);

for i=1:len(volt)
    volt(i) = volt(i)+out;
end
for j=1: