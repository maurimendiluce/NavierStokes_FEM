%Estimacion de orden

h = [1/16,1/32,1/64,1/128];
error = [0.0952183,0.0673277,0.0476136,0.0336701];

plot(log(h),log(error),'--o',log(h),log(h),'k-.')
xlabel("log(h)")
ylabel("log(error)")
title("Estimación de orden")
legend("Orden P2-P1","O(h)")
%polyfit(log(h),log(error),1)