%Quick Feigenbaum calculation, ps=0 for period-4 trajectory
%Taking data from LyapStats and LyapStatsZoom
%Feigenbaum constant for Logistic map universality class: 4.6692...
%Feigenbaum constant for Hénon map universality class: 8.7210...
%\rvwg2022

a1 = cursor_info.Position(1);	%3.6276
a2 = cursor_info.Position(1);	%3.6666
a3 = cursor_info.Position(1);   %3.6725
a4 = cursor_info.Position(1);   %3.6737
a5 = cursor_info.Position(1);   %3.6739

delta1 = (a2-a1)/(a3-a2);   %6.6102
delta2 = (a3-a2)/(a4-a3);   %5.0427
delta3 = (a4-a3)/(a5-a4);   %4.6800