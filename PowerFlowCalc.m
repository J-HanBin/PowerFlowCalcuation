% PowerFlowCalc
% PowerFlow Calcuation Program by MATLAB Ver. R2018a

% Y Bus Calculator
% 임의의 n개의 bus로 이루어진 전력 시스템에서 선로 정보의 데이터를 통해 어드미턴스 매트릭스를 만드는 코드
% 입력 받는 데이터의 양식은 다음과 같은 형태로 만들어져 있어야 한다
% 
% ex)
% % 3G 9 Bus Power Flow Data                       <- line data format
% % 2004.03.09 Seon-Ju Ahn
% % From	To	      R	          X	      Charging
%    4		1		0		  0.0576		0
%    7		2		0		  0.0625		0
%    9		3		0		  0.0586		0
%    7		8		0.0085    0.072		    0.149
%    9		8		0.0119	  0.1008	    0.209
%    7		5		0.032	  0.161         0.306
%    9		6		0.039	  0.17          0.358
%    5		4		0.01	  0.085         0.176
%    6		4		0.017	  0.092         0.158
 
load Line.dat

% 데이터 유효성 검사 (선로 정보에 대한 중복 입력, From, To가 잘못 입력된 경우 방지)
for i=1:size(Line,1)
    if Line(i,1) == Line(i,2)   % From 숫자와 To 숫자가 같게 입력된 경우
       sprintf('WRONG LINE DATA : FROM NUM CANNOT BE THE SAME WITH TO NUM')
       return
    end
    
    for j=(i+1):size(Line,1)
        if (Line(i,1)==Line(j,1) && Line(i,2)==Line(j,2)) || (Line(i,1)==Line(j,2) && Line(i,2)==Line(j,1))     % 선로 정보가 중복 입력된 경우
            sprintf('WRONG LINE DATA : DUPLICATED LINE DATA CANNOT EXIST')
            return
        end
    end
end

% Admittance Matrix (Yb)
FromTo = Line(:,1:2);
Yb = zeros(max(FromTo(:)));     % definition of the size of Admittance Matrix

for i=1:size(Line,1)
    m1 = Line(i,1);     m2 = Line(i,2); 
    R = Line(i,3);  X = Line(i,4);  B = Line(i,5);
    zl = complex(R,X);   % Series Impedance
    yc = complex(0,B);     % Shunt Admittance
    
    Yb(m1, m1) = Yb(m1, m1) + 1/zl + yc/2;    % Self-Admittance
    Yb(m2, m2) = Yb(m2, m2) + 1/zl + yc/2;    % Self-Admittance
    Yb(m1, m2) = -1/zl;
    Yb(m2, m1) = -1/zl;
end

G = real(Yb);   B = imag(Yb);     % Conductance G & Susceptance S

% Bus Data Input
% 임의의 n개의 bus로 이루어진 전력 시스템에서 각 Bus의 데이터를 입력받는 코드
% 입력 받는 데이터의 양식은 다음과 같은 형태로 만들어져 있어야 한다
% 
% ex)
% % 3G 9 Bus Power Flow Data             <- Bus data format
% % 2004.03.09 Seon-Ju Ahn
% % Bus Type
% %					1: Slack
% %					2: Generator
% %					3: Load
% %					Sbase : 100MVA
% % Num	  Type    PG		  QG		PL	  QL	  VM	   VAngle
%   1	   1       0 	  	   0	 	0	   0       1.04	    0
%   2	   2	  1.63         0		0	   0	   1        0
%   3	   3	   0           0		0	   0	   1        0
%   4	   3	   0		   0		0	   0	   1	 	0
%   5	   3	   0           0        1.25   0.5	   1		0
%   6	   3	   0		   0		0.9    0.3	   1		0
%   7	   3 	   0	       0		0	   0	   1		0
%   8	   3	   0		   0		1      0.35    1		0
%   9	   3	   0		   0		0	   0	   1		0

load Bus.dat

% 데이터 유효성 검사 (버스 정보에 대한 중복 입력, 정보 누락 방지)
Bus = sortrows(Bus);     % 버스 정보 1열 기준으로 오름차순 정렬
for i=1:size(Bus,1)
   if Bus(i,1) ~= i
       sprintf('WRONG BUS DATA EXIST')
       return
   end
end

v = Bus(:,7);   % Voltage Magnitude
theta = Bus(:,8);   % Phase Angle
P = zeros(size(Bus,1),1);   % Active Power
Q = zeros(size(Bus,1),1);   % Reactive Power

% 주어진 유효전력 및 무효전력 값의 버스 넘버를 저장하기 위한 배열 정의
% 모르는 위상각과 전압의 크기의 버스 넘버에 대한 정보도 알 수 있음
% ex) pa에 저장된 버스 -> 유효전력 값이 주어짐 -> 해당 버스의 위상각이 미지수
% ex) qa에 저장된 버스 -> 무효전력 값이 주어짐 -> 해당 버스의 전압의 크기가 미지수
pn = 0;     qn = 0;     % pa와 qa의 원소 개수를 세기 위한 변수

for i=1:size(Bus,1)
    if Bus(i,2) == 2
        pn = pn + 1;
    else if Bus(i,2) == 3
         pn = pn + 1;   qn = qn + 1;   
        end
    end
end

pa = zeros(pn, 1);    qa = zeros(qn, 1);   
pn = 1;     qn = 1;

for i=1:size(Bus,1)
    if Bus(i,2)==1  % Slack Bus & Reference Bus -> 전력 공급 -> P>0, Q>0
        P(i) = Bus(i,3);   Q(i) = Bus(i,4);     % P,Q 주어지지 않음 - 없어도 상관 X
    else if Bus(i,2)==2     % PV Bus -> 전력 공급 -> P>0, Q>0
             P(i) = Bus(i,3);   Q(i) = Bus(i,4);    % Q 주어지지 않음 - 없어도 상관 X
             pa(pn) = i;
             pn = pn + 1;
    else if Bus(i,2)==3     % PQ Bus -> 전력 소비 -> P<0, Q<0
            P(i) = -Bus(i,5);    Q(i) = -Bus(i,6);
            pa(pn) = i;     qa(qn) = i;
            pn = pn + 1;    qn = qn + 1;
        else
            sprintf('WRONG BUS TYPE : TYPE NUM SHOULD BE 1 TO 3')
            return
        end
        end
    end
end

% Initial Condition (미지수 x 초기화)
xp = zeros(length(pa),1);   xq = ones(length(qa),1);
x = [xp; xq];

% for initialization of iteration criterion
del_x = inf;
n = 0;

while abs(del_x) > 0.0001    % iteration criterion
n = n+1;

% Power Flow Equations
fp = zeros(length(pa),1);   fq = zeros(length(qa),1);
for i=1:size(Bus,1)  
   for j=1:length(pa)
       k = pa(j);
       
       if i==1
          fp(j) = fp(j) - P(k); 
       end
       
       fp(j) = fp(j) + v(k)*v(i)*(G(k,i)*cos(theta(k)-theta(i))+B(k,i)*sin(theta(k)-theta(i)));
   end
   
   for j=1:length(qa)
      k = qa(j);
      
      if i==1
          fq(j) = fq(j) - Q(k); 
      end
       
      fq(j) = fq(j) + v(k)*v(i)*(G(k,i)*sin(theta(k)-theta(i))-B(k,i)*cos(theta(k)-theta(i)));
   end
end

f = [fp; fq];

% Jacobian Matrix (J)
J11 = zeros(length(pa));    % dP/d(theta)
for i=1:length(pa)
    k = pa(i);
    
   for j=1:length(pa)
      if k==pa(j)
          J11(i,j) = -Q(k)-B(k,k)*(v(k))^2;
      else  % i~=j
          J11(i,j) = v(k)*v(pa(j))*(G(k,pa(j))*sin(theta(k)-theta(pa(j)))-B(k,pa(j))*cos(theta(k)-theta(pa(j))));
      end
   end
end

J22 = zeros(length(qa));    % dQ/dv
for i=1:length(qa)
   k = qa(i);
   
   for j=1:length(qa)
       if k==qa(j)
          J22(i,j) = Q(k)/v(k)-B(k,k)*v(k);
       else     % i~=j
           J22(i,j) = v(k)*(G(k,qa(j))*sin(theta(k)-theta(qa(j)))-B(k,qa(j))*cos(theta(k)-theta(qa(j))));
       end
   end
end

J12 = zeros(size(J11,1),size(J22,2));   % dP/dv
for i=1:length(pa)
    k = pa(i);
    
   for j=1:length(qa)
       if k==qa(j)  % pa(i) == qa(j)
           J12(i,j) = P(k)/v(k)+G(k,k)*v(k);
       else     % pa(i) ~= qa(j)
           J12(i,j) = v(k)*(G(k,qa(j))*cos(theta(k)-theta(qa(j)))+B(k,qa(j))*sin(theta(k)-theta(qa(j))));
       end
   end
end

J21 = zeros(size(J22,1),size(J11,2));   % dQ/d(theta)
for i=1:length(qa)
    k = qa(i);
    
   for j=1:length(pa)
       if k==pa(j)  % qa(i) == pa(j)
           J21(i,j) = P(k)-G(k,k)*(v(k))^2;
       else     % qa(i) ~= pa(j)
           J21(i,j) = v(k)*v(pa(j))*(-G(k,pa(j))*cos(theta(k)-theta(pa(j)))-B(k,pa(j))*sin(theta(k)-theta(pa(j))));
       end
   end
end

J = [J11 J12; J21 J22];

% Newton Raphson Method
del_x = J\f;
x = x - del_x;

for i=1:length(x)
   if i<=length(pa)
       theta(pa(i)) = x(i);
   else     % i>length(pa)
       v(qa(i-length(pa))) = x(i);
   end
end

end

% Calculate P & Q of slack bus & PV bus
for i=1:size(Bus, 1)
    if Bus(i, 2) == 1   % slack bus
        for j=1:size(Bus, 1)
            P(i) = P(i) + v(i)*v(j)*(G(i,j)*cos(theta(i)-theta(j))+B(i,j)*sin(theta(i)-theta(j)));
            Q(i) = Q(i) + v(i)*v(j)*(G(i,j)*sin(theta(i)-theta(j))-B(i,j)*cos(theta(i)-theta(j)));
        end
    else if Bus(i, 2) == 2  % PV bus
            for j=1:size(Bus, 1)
                Q(i) = Q(i) + v(i)*v(j)*(G(i,j)*sin(theta(i)-theta(j))-B(i,j)*cos(theta(i)-theta(j)));
            end
        end
    end
end

theta_deg = zeros(length(theta), 1);
for i=1:length(theta)
   theta_deg(i)=theta(i)/pi*180; % 위상각 단위를 라디안에서 도로 변환 
end

% PU Method
Sbase = 30;  % 단위 GVA
Vbase = 345;    % 단위 kV
Ibase = Sbase*1000/Vbase;   % 단위 kA

vreal = Vbase*v;
Preal = Sbase*P;    Qreal = Sbase*Q;

% Calculate Power Flow Current
% Ireal(i,j) -> i bus에서 j bus로 흐르는 조류량
Ireal = zeros(size(Bus, 1));
for i=1:size(Bus, 1)
    for j=1:size(Bus, 1)
        Ireal(i, j) = abs(Yb(i, j)) * (vreal(i) - vreal(j));
    end
end

I = Ireal/Ibase;  % p.u.값으로 계산된 조류량

% Results
% 각 버스별 전압(크기, 위상각)과 유효전력 무효전력 값에 대한 계산 결과 값 & 그래프 출력
% 각 선로별 조류 양에 대한 계산 결과 값 출력
fprintf('전압의 크기(Vm) [p.u.]\n')
fprintf('------------------------------------------\n')
for i=1:length(v)
   fprintf('Bus %d : %.4f\n', i, v(i)) 
end
fprintf('\n')

fprintf('전압의 위상각(theta) [rad]\n')
fprintf('------------------------------------------\n')
for i=1:length(theta)
   fprintf('Bus %d : %.4f\n', i, theta(i)) 
end
fprintf('\n')

fprintf('전압의 크기(Vm) [kV]\n')
fprintf('------------------------------------------\n')
for i=1:length(vreal)
   fprintf('Bus %d : %.2f\n', i, vreal(i)) 
end
fprintf('\n')

fprintf('유효전력(P) [GW]\n')
fprintf('------------------------------------------\n')
for i=1:length(Preal)
    if Preal(i) > 0     % 전력 공급
        fprintf('Bus %d : %.1f (유효전력 %.1f [GW] 공급)\n', i, Preal(i), abs(Preal(i))) 
    else    % Preal(i) <= 0     % 전력 소비
        fprintf('Bus %d : %.1f (유효전력 %.1f [GW] 소비)\n', i, Preal(i), abs(Preal(i)))
    end
end
fprintf('\n')

fprintf('무효전력(Q) [GVar]\n')
fprintf('------------------------------------------\n')
for i=1:length(Qreal)
    if Qreal(i) > 0     % 전력 공급
        fprintf('Bus %d : %.1f (무효전력 %.1f [GVar] 공급)\n', i, Qreal(i), abs(Qreal(i)))
    else    % Qreal(i) <= 0     % 전력 소비
        fprintf('Bus %d : %.1f (무효전력 %.1f [GVar] 소비)\n', i, Qreal(i), abs(Qreal(i)))
    end
end
fprintf('\n')

fprintf('각 모선에서 조류량(I(i, j) From i To j) [p.u.]\n')
fprintf('------------------------------------------\n')
for i=1:size(I, 1)
    for j=1:size(I, 2)
        fprintf('%4.2f\t', I(i,j))
    end
    fprintf('\n')
end
fprintf('\n')

fprintf('각 모선에서 조류량(I(i, j) From i To j) [kA]\n')
fprintf('------------------------------------------\n')
for i=1:size(I, 1)
    for j=1:size(I, 2)
        fprintf('%.2f\t', Ireal(i,j))
    end
    fprintf('\n')
end
fprintf('\n')

% Goal of Voltage Control (전압조정목표) [p.u.]
% 기준전압 345kV는 전압조정목표를 336kV ~ 360kV로 한다 by 전력계통 신뢰도 및 전력품질 유지기준
% 336kV (0.9739 p.u.) / 360kV (1.0435 p.u.) 
yline = ones(size(Bus, 1), 1);
yline1 = (1 - (Vbase-336)/Vbase) * yline;
yline2 = (1 + (360-Vbase)/Vbase) * yline;

figure(1);  % First Graph
yyaxis left
plot(Bus(:, 1), v, '-o');
hold on;
plot(Bus(:, 1), yline1, '--');
plot(Bus(:, 1), yline2, '--');
title('Votage Magnitue(Vm) & Phase Angle(Theta)');
ylabel('Voltage[p.u.]');
xlabel('경남(1)  강원(2)  서울(3)  전북(4)  경북(5)  전남(6)  충북(7)  경기(8)  충남(9)');
grid on;

yyaxis right
plot(Bus(:, 1), theta_deg, '-*');
ylabel('Angle[degree]')

figure(2);  % Second Graph
PQ = [Preal Qreal];
bar(Bus(:, 1), PQ);
title('Active Power(P) & Reactive Power(Q)')
ylabel('P/Q[GW/GVar](+:gen,-:consum)');
xlabel('경남(1)  강원(2)  서울(3)  전북(4)  경북(5)  전남(6)  충북(7)  경기(8)  충남(9)');
legend('Active Power', 'Reactive Power');
grid on;
