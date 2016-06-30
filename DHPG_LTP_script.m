close all

pulse_interval=20;% in seconds
baseline_duration=1200; % in seconds
stimulus_duration=300;
stimulus='DHPG 5 min';
%%
j=1;
for i = 1:length(Region1IcSlopemVms)
    if ~isnan(Region1IcSlopemVms(i))
        Region1IcSlopemVms_final(j,1)=abs(Region1IcSlopemVms(i));
        j=j+1;
    end
end

Times_final=linspace(-baseline_duration,length(Region1IcSlopemVms_final)*pulse_interval-baseline_duration,length(Region1IcSlopemVms_final));
Times_final_m=minutes(seconds(Times_final)); %to alter and put to work in


%% Plotting when all data collected

% Times_final=linspace(0,length(Region1IcSlopemVms_final)*pulse_interval,length(Region1IcSlopemVms_final));
% Times_final=seconds(Times_final);
% Times_final_m=duration(Times_final,'Format','m'); %to alter and put to work in 

slopes_tma = tsmovavg(Region1IcSlopemVms_final,'t',20,1);
baseline = mean(Region1IcSlopemVms_final(1:baseline_duration/pulse_interval)); %baseline is given by a simple average of the base line period
slope_normalized = Region1IcSlopemVms_final./baseline;
slope_normalized_tma = tsmovavg(slope_normalized,'s',10,1);

%in order to display graphs in minutes
baseline_duration_m=minutes(seconds(baseline_duration));
stimulus_duration_m=minutes(seconds(stimulus_duration));

%% Averaged for 6 points (unlike moving average, points averaged do not repeat)

Averaged_slope=NaN(1,length(slope_normalized));

j=1;
while j<=length(slope_normalized)
    Averaged_slope(j)=mean(slope_normalized(j:j+5));
    j=j+6;
end



%% Plotting
f3=figure
hold on
title({'DHPG induced LTD, 1st try';' '})
scatter(round(Times_final_m,1),Averaged_slope,8)
%plot(Times_final_m,slope_normalized_tma)
line([0 stimulus_duration_m],[1.5 1.5],'LineWidth',4,'Color','r')
end_of_stimulus=baseline_duration_m+stimulus_duration_m;

line([Times_final_m(1) Times_final_m(end)],[1 1],'LineWidth',0.8,'Color','k','LineStyle','--')
line([0 0],[0 1.5],'LineWidth',0.8,'Color','k','LineStyle','-.')
line([stimulus_duration_m stimulus_duration_m],[0 1.5],'LineWidth',0.8,'Color','k','LineStyle','-.')

text(stimulus_duration_m/2,1.51,stimulus,'HorizontalAlignment','center','VerticalAlignment','bottom');
fig3=gca;
fig3.YTick=(0:0.2:1.5);
axis([Times_final_m(1) Times_final_m(end) 0 1.5 ])
xlabel({' ';'Duration (min)'})
ylabel({'fEPSP slope normalized to baseline (mV/s)';' '})

%legend('recorded values','moving average')
hold off

% poly_coeff=[];
% 
% for j=1:length(o)
%      o_temp1=o{j};  
%      o_temp=o_temp1(7:end,:);
%      [m,n]=size(o_temp);
%      
%         for i=1:n
%             [mm,~]=size(o_temp(:,i));
%             o_temp_cell=cell(mm,1);
%             o_temp_cell(:,1)=o_temp(:,i);
% 
%             mat_o_sub=cell2mat(o_temp_cell);
%             xx=(1:length(mat_o_sub));
%             poly_coeff(i)=polyfit(xx',mat_o_sub,1);
%             
%         end
% end