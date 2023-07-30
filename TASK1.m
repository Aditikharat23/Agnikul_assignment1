%% TASK 1

% Inbuilt function used : 1. factorial(n)

%% 
clear all
% sum = 0;
% tracksum = zeros(1,50);
% for n = 0: 2
%     temp  = (fact(4*n)/(fact(n))^4)*((1103+(26390*n))/(396)^(4*n));
%     sum = sum+temp;
%     tracksum(n+1) = sum;
% end
% temp2 = sqrt(8)*sum/9801;
% 
% pival = 1/temp2;

%%Unable to get accuracy of more than 30 digits
%Inbuilt matlab function used digits to get the precision required.
digits(31416)
pi_val = vpa(pi);
%Convert to character
c = char(pi_val);
freq = zeros(1,10);
for i= 3:length(c)
    temp =str2num(c(i));
    freq(1,temp+1) = freq(1,temp+1) +1;
end
h=histogram('Categories',{'0','1','2','3','4','5','6','7','8','9'},'BinCounts',[freq(1,:)])
title('Frequency of numbers in 31415 decimal places of pi ')
