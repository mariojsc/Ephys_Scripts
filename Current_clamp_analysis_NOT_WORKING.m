%% Current clamp analysis

    
%% current clamp files organization

types={'CC.txt'}; %defines columns of o (add here for more :))
variable=0;

for i=1:length(CC)
    aC(i)=(height(CC(i).values));
end

c1=char(types);    
CC_T{i,j}=ephyst_extract(CC,c1,variable,aC,0); %turns negative values to absolute



%% Statistics
r=1;
%for total average frequency, each cell
dur=3*60*1000;
cache=o{1,1};

for i=1:length(a)
    mean_frequency(1:4,i)=cache(1:4,i);
    mean_frequency{5,i}=1000*(a(i)-1)/dur;
end

statistics(r).parameter='mean_frequency_each_cell';
statistics(r).values=mean_frequency;
r=r+1;

statistics(r).parameter='mean_frequency_overall';
mean_frequencies=cell2mat(mean_frequency(5,:));
statistics(r).values=mean(mean_frequencies);
r=r+1;

statistics(r).parameter='SEM_frequency_overall';
SEM=std(mean_frequencies)/sqrt(length(mean_frequencies));
statistics(r).values=SEM;
r=r+1;

%average and median amplitude
cache=o{2,1};
for i=1:length(a)
    av_amplitude(1:4,i)=cache(1:4,i);
    array=cell2mat(cache(6:a(i)+5,i));
    av_amplitude{5,i}=-1*sum(array)/a(i);
    av_amplitude{6,i}=-1*median(array);
end

statistics(4).parameter='av_amplitude_each_cell';
statistics(4).values=av_amplitude;

%for Vm and Ic of the first spike during CC

[~,m]=size(oT{1,1});
cacheV=oV{1,1};
cacheIc=oT{1,1};
AP_threshold=cell(size(cacheIc));

for i=1:m
    AP_threshold(1:4,i)=cacheIc(1:4,i);
    AP_threshold{5,i}=cacheV{5+cacheIc{6,i},i};
    AP_threshold{6,i}=I_steps(cacheIc{6,i});
    
end

statistics(2).parameter='AP_threshold';
statistics(2).values='AP_threshold';


%to get absolute values


%statistics will come here

%% exports files and gives names to the csv
 [n,m]=size(o);
 for j=1:m
    for i=1:n
        T=cell2table(o{i,j});
        writetable(T,strcat(char(types(j)),'--',char(variables(i)),'.csv'));
    end
 end

%% APs Analysis

%to get: number of AP per current inj, IV curve, AP threshold

types={'CC1_curve.txt'}; %defines columns of o (add here for more :))
variables={'CC1: IN 0 Mean'};%,'Interevent Interval (ms)'}; %defines rows of o

%for the special case of trace
for j=1:length(types)
    c1=char(types(j));    
    for i=1:length(variables)
        oV{i,j}=ephyst_extract(recA,c1,variables(i),aA,1); %turn into cell to export later
    end
end

variables={'CC1: _Cmd 0 at Cursor 1'}; %now it is necessary for variable str to be inside cells!!! look at the first test inside the function
I_steps=ephyst_extract(recA,c1,variables(1),aA,1); %revise!!!
I_curve=cell2table(I_steps);
I_steps=cell2mat(table2array(I_curve(6:45,1)));

types={'CC1_spiket.txt'}; %defines columns of o (add here for more :))
variables={'Trace'};  %,'Interevent Interval (ms)'}; %defines rows of o

%for the special case of trace
for j=1:length(types)
    c1=char(types(j));    
    for i=1:length(variables)
        oT{i,j}=ephyst_extract(recB,c1,variables(i),aB,1); %turn into cell to export later
    end
end

target=cell(size(oT{1})); % this compromises code!!! because index for cell of interest just appears inside the next loop
target_z=zeros(size(target));

for i=1:numel(target)
target{i}=target_z(i);
end

for jj=1:length(oT)
    origin_cache=oT{jj};
    
    [~,m]=size(origin_cache);
    
    for i=1:m %AP threshold
    
    target(1:5,i)=origin_cache(1:5,i);
    target{6,i}=I_steps(origin_cache{6,i},1);

    x=6;    
    
    while x<length(origin_cache(5:end,i))

            inc=1;
            jjj=x;
            
            while origin_cache{jjj,i}==origin_cache{jjj+1,i}
            inc=inc+1;
            jjj=jjj+1;
            end
            
            if origin_cache{x,i}>=1
                target{origin_cache{x,i}+6,i}=inc;
            end
               x=x+inc;
        end
    end
end

%%
%read Membrane Constant files and organize be same as above order

% filename = '/Users/mariocarvalho/Documents/MATLAB/NCBL/ephys/MSUS2#_ephys_Recs/~Membrane_constants.csv';
% delimiter = ',';
% formatSpec = '%s%s%s%s%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% [m,~]=size(dataArray{1,1});
% [~,n]=size(dataArray);
% Membrane_constants=cell(m,n);
% 
% for i=1:n
%     for j=1:m
%         cache=dataArray{1,i};
%         Membrane_constants{j,i}=cache{j,1};
%     end
% end
% 
% organized_mc=cell(size(Membrane_constants)');
% 
% [~,m_c]=size(o{1,1});
% cache_o=o{1,1};
% id=o(1:4,:);
% j=1;

% while i<=m
%     ido=id{1,j};
%     
%     if strcmp(
%     
%     end
% end
% 

%oT and oV are of different length because not all traces of IV have
%spike?!


% ot=zeros(n,m);

% for i=1:m
%     j=1;
%     jt=1;
%     while j<=n && jt<=n
%         if (o{j,i})>0
%         ot(jt,i)=ot(jt,i)+1;
%         j=j+1;
%         if j<=n
%             while o{j,i}==o{j-1,i} 
%             ot(jt,i)=ot(jt,i)+1;
%             j=j+1;
%             end
%         end
%         jt=jt+1;
%         end
%         j=j+1;
%     end
% end

% statistics will come here

%exports files and gives names to the csv
