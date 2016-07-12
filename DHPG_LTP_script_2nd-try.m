close all

pulse_interval=20;% in seconds
baseline_duration=1200; % in seconds
stimulus_duration=300;
stimulus='DHPG';
%%
j=1;
for i = 1:length(Region1IcSlopemVms)
    if ~isnan(Region1IcSlopemVms(i))
        Region1IcSlopemVms_final(j,1)=abs(Region1IcSlopemVms(i));
        j=j+1;
    end
end

%% Truncation step
Region1IcSlopemVms_final=Region1IcSlopemVms_final(193:end);

%%
%Times_final=linspace(-baseline_duration,length(Region1IcSlopemVms_final)*pulse_interval-baseline_duration,length(Region1IcSlopemVms_final));

Times_final=linspace(0,length(Region1IcSlopemVms_final)*pulse_interval-baseline_duration,length(Region1IcSlopemVms_final));

Times_final=linspace(-baseline_duration,length(Region1IcSlopemVms_final)*pulse_interval-baseline_duration,length(Region1IcSlopemVms_final));

Times_final_m=minutes(seconds(Times_final)); %to alter and put to work in


%% Plotting when all data collected

% Times_final=linspace(0,length(Region1IcSlopemVms_final)*pulse_interval,length(Region1IcSlopemVms_final));
% Times_final=seconds(Times_final);
% Times_final_m=duration(Times_final,'Format','m'); %to alter and put to work in 

%slopes_tma = tsmovavg(Region1IcSlopemVms_final,'t',30,1);
%baseline = mean(Region1IcSlopemVms_final(1:baseline_duration/pulse_interval)); %baseline is given by a simple average of the base line period
baseline_end=38;
baseline = mean(Region1IcSlopemVms_final(1:baseline_end)); %checked by looking at the graph and then the indexes on the array
%baseline = 1;

%% correction of outliers

slope_normalized = Region1IcSlopemVms_final./baseline;
slope_normalized(13:18)=1;
slope_normalized(49:55)=1;

baseline_of_baseline = mean(slope_normalized(1:baseline_end));

slope_normalized_tma = tsmovavg(slope_normalized,'t',20,1);

%in order to display graphs in minutes
baseline_duration_m=minutes(seconds(baseline_duration));
stimulus_duration_m=minutes(seconds(stimulus_duration));

%% Averaged for 6 points (unlike moving average, points averaged do not repeat)

Averaged_slope=NaN(1,length(slope_normalized));

j=1;
while j+5<=length(slope_normalized) 
    Averaged_slope(j)=mean(slope_normalized(j:j+5));
    j=j+6;
end

LTD_percentage=nanmean(slope_normalized_tma((end-54):end));

%% Truncate slopes

% left_trunc=193;
% Averaged_slope_truncated=Averaged_slope(left_trunc:end);
% Averaged_slope_truncated(13)=mean(Averaged_slope(1:77));
% Averaged_slope_truncated(49)=mean(Averaged_slope(1:77));
% 
% 
% Times_final=linspace(0,length(Averaged_slope_truncated)*pulse_interval,length(Averaged_slope_truncated));
% Times_final_m_trunc=minutes(seconds(Times_final)); %to alter and put to work in

%% Plotting
f3=figure
hold on
title({'DHPG induced LTD, 2nd try';' '})

scatter(round(Times_final_m,1),Averaged_slope,8)

%scatter(round(Times_final_m,1),Averaged_slope,8)


plot(Times_final_m,slope_normalized_tma)

%line([0 stimulus_duration_m],[1.5 1.5],'LineWidth',4,'Color','r')

end_of_stimulus=baseline_duration_m+stimulus_duration_m;

line([Times_final_m(1) Times_final_m(end)],[1 1],'LineWidth',0.8,'Color','k','LineStyle','--')
%line([0 0],[0 1.5],'LineWidth',0.8,'Color','k','LineStyle','-.')
%line([stimulus_duration_m stimulus_duration_m],[0 1.5],'LineWidth',0.8,'Color','k','LineStyle','-.')

%text(stimulus_duration_m/2,1.51,stimulus,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(4,1.52,stimulus,'HorizontalAlignment','center','VerticalAlignment','bottom');
line([0 8],[1.5 1.5],'LineWidth',4,'Color','r')

%text((33-28)/2+28,1.52,stimulus,'HorizontalAlignment','center','VerticalAlignment','bottom');
%line([28 33],[1.5 1.5],'LineWidth',4,'Color','r')

fig3=gca;
fig3.YTick=(0:0.2:4);
axis([Times_final_m(1) Times_final_m(end) 0 2 ])
xlabel({' ';'Duration (min)'})
ylabel({'fEPSP slope normalized to baseline (mV/s)';' '})

%legend('recorded values','moving average')
hold off