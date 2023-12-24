close all 
clear 
c_f_range=linspace(0,15);
H=3.49;
K=0.005;
q_Henry=H.*c_f_range;
q_Langmuir=H.*c_f_range./(1+H.*c_f_range);
v=1.2;
q_Freundlich=H.*c_f_range.^(1/v);
% figure(1)
% % plot(c_f_range,q_Henry,'LineWidth',2.0)
% % hold on 
% plot(c_f_range,q_Langmuir)
% % hold on 
% % plot(c_f_range,q_Freundlich,'LineWidth',2.0)
% legend('Henry','Langmuir','Freundlich')
% hold off
% xlabel('c(g/L) concentration in fluid phase')
% ylabel('n(g/L) uptake in solid phase')
% ax = gca; 
% ax.FontSize = 15;