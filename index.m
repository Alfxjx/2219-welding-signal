% ����ļ��������������� matlab �Ĺ�������ɣ��ڱ��ļ����´�Ű����ŷָ������ݡ�
% ���������ǵ�һ�е�ѹ���ڶ��е���
% ���ȵ������� �ֱ𱣴�Ϊ dadi_ �� gm_ ��ͷ�ı�����
importData

% �����ݽ���һ���򵥵Ĺ۲�
showbrief

% ����ȡ��λ�ú�ȡ���ĳ���
start1 = input('input start dadi: ');
step1 = input('input step dadi: ');
% start2 = input('input start gm: ');
% step2 = input('input step gm: ');
% ��ʾһ�´�������С��
dadi_volt = dadi(start1:(start1+step1),1);
dadi_curr = dadi(start1:(start1+step1),2);

% gm_volt = gm(start2:(start2+step2),1);
% gm_curr = gm(start2:(start2+step2),2);

showall

% �����ݽ��� emd ����
% �����õ�IMF�о���Ƶ��
volt = dadi_volt;
current = dadi_curr;
emdTest
fftTest
% �����Ӧ��IMF�� HHT �任
hhtTest

% �л�һ��ע��

% volt = gm_volt;
% current = gm_curr;
% emdTest
% % �����Ӧ��IMF�� HHT �任
% hhtTest

% ����ʱƵ��
tfEN


