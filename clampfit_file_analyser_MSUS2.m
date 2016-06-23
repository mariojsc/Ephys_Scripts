%cd the Rec folder with animals
%the code is A# for animal, s# for slice, c# for cell, 'holding current'.
%txt for file

%there are Nan values on top of values in value cells of o and
%such. This is the result of ephyst_extract starting applying str2double
%from the very top of the files, so variable name gets Nan, and for instant
%frequency, also the first N/A

cd



%% import parameters for trimming

fileID=fopen('~analysis_parameters.txt');
P=textscan(fileID,'%s %d');
fclose(fileID);
parameter_strct=struct('parameter',P{1,1});
values=P{1,2};

for i=1:length(parameter_strct)
    parameter_strct(i).value=values(i);
    parameter_strct(i).exclusion_vector=[];
end


%% import blinding_correspondence

clear P

fileID=fopen('~blinding.txt');
P=textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);



%% import data sets

r=0;
r1=0;
r2=0;
r3=0;
r4=0;
fld1=dir('A*');

for i=1:length(fld1);
%go to animal
oldfolder=cd(fld1(i).name);

fld11=dir('s*'); 
    
    for j=1:length(fld11);
    %go to slice
    oldfolderr=cd(fld11(j).name);    
    
    fld111=dir('c*');
        
        for x=1:length(fld111);
        %go to cell
        oldfolderrr=cd(fld111(x).name);
        
        % Get all files in the current cell
        
        filesc= dir('*.atf');
        if length(filesc)>=1
            for id = 1:length(filesc)
                    % Get the file name (minus the extension)
                    [~, f] = fileparts(filesc(id).name);
                    % rename
                    movefile(filesc(id).name, sprintf('%s.txt',f));
            end
        end
        
        filesc= dir('-55_bl_results_E.txt'); %change name of this type of file
        if length(filesc)>=1
            for id = 1:length(filesc)
                    % Get the file name (minus the extension)
                    [~, f] = fileparts(filesc(id).name);
                    % rename
                    movefile(filesc(id).name,'-55_E.txt');
            end
        end
        
        files = dir('-55_E.txt');
        if (length(files))~=0
            for id = 1:length(files);
                        % creates tables with values from the converted files
                        
                        B=readtable(files(id).name,'HeaderLines',2,'Delimiter','tab','ReadVariableNames',false);
                        %stores in struct
                        rec(id+r).type=files(id).name;
                        rec(id+r).values=B;
                        rec(id+r).animal=fld1(i).name;
                        rec(id+r).slice=fld11(j).name;
                        rec(id+r).cell=fld111(x).name; 
            end
            r=r+id;
        end
        
        files = dir('*.png');
        if (length(files))~=0
            for id = 1:length(files);
                
                        % stores in struct
                        recD(id+r4).type=files(id).name;
                        recD(id+r4).values=files(id).date;
                        recD(id+r4).animal=fld1(i).name;
                        recD(id+r4).slice=fld11(j).name;
                        recD(id+r4).cell=fld111(x).name;
            end
            r4=r4+id;
        end
        
        files = dir('CC1_curve.txt');
        if (length(files))~=0
            for id = 1:length(files);
                        % creates tables with values from the converted files
                        
                        B=readtable('CC1_curve.txt','HeaderLines',2,'Delimiter','tab','ReadVariableNames',false);
                        %stores in struct
                        recA(id+r1).type=files(id).name;
                        recA(id+r1).values=B;
                        recA(id+r1).animal=fld1(i).name;
                        recA(id+r1).slice=fld11(j).name;
                        recA(id+r1).cell=fld111(x).name; 
            end
            r1=r1+id;
        end
        
        files = dir('CC.txt');
        if (length(files))~=0
            for id = 1:length(files);
                        % creates tables with values from the converted files
                        
                        B=readtable('CC.txt','HeaderLines',0,'Delimiter','tab','ReadVariableNames',false);
                        %stores in struct
                        CC(id+r2).type=files(id).name;
                        CC(id+r2).values=B;
                        CC(id+r2).animal=fld1(i).name;
                        CC(id+r2).slice=fld11(j).name;
                        CC(id+r2).cell=fld111(x).name;
                        
            end
            r2=r2+id;
        end
        
        files = dir('CC1_spiket.txt');
            if (length(files))~=0
                for id = 1:length(files);
                            % creates tables with values from the converted files

                            B=readtable(files(id).name,'HeaderLines',2,'Delimiter','tab','ReadVariableNames',false);
                            %stores in struct
                            recB(id+r3).type=files(id).name;
                            recB(id+r3).values=B;
                            recB(id+r3).animal=fld1(i).name;
                            recB(id+r3).slice=fld11(j).name;
                            recB(id+r3).cell=fld111(x).name; 
                end
                r3=r3+id;
            end
        cd(oldfolderrr)
        end
    cd(oldfolderr) 
    end
cd(oldfolder)
end

%take rec structure and generate a blinding.txt with fields for each
%parameter. This text file has to be filed in order to separate groups and
%perform statistics

clear cache
cache=dir;

v=0;
for i=1:length(cache)
    ccache=cache(i).name;
        if strcmp(ccache,'~blinding.txt')
        v=v+1;
        end
end

if v==0

    
fID=fopen('~blinding.txt','wt');
fprintf(fID,'type animal slice cell pharmacology experimental_group\n');

for i=1:length(rec)
    type=rec(i).type;
    animal=rec(i).animal;
    slice=rec(i).slice;
    cell=rec(i).cell;
    fprintf(fID,'%s %s %s %s\n',type,animal,slice,cell);
end
    fclose(fID);
    
end

%% Make file size tables
    %take the size of the table from biggest size peak tables
    for i=1:length(rec)
        [n,~]=size(rec(i).values);
        a(i)=n;
    end 
   
    for i=1:length(recA)
        [n,~]=size(recA(i).values);
        aA(i)=n;
    end     
    
    for i=1:length(recB)
        [n,~]=size(recB(i).values);
        aB(i)=n;
    end     
    
%% Current clamp Analysis

types={'CC.txt'}; %defines columns of o (add here for more :))
variable=0;

for i=1:length(CC)
    aC(i)=(height(CC(i).values));
end

c1=char(types);    
CC_T{i,j}=ephyst_extract(CC,c1,variable,aC,0); %turns negative values to absolute



%% mini events files parsing

types={'-55_E.txt'}; %defines columns of o (add here for more :))
variables={'Inst. Freq. (Hz)','Peak Amp (pA)'}; %defines rows of o

for j=1:length(types)
    
    c1=char(types(j));    
    for i=1:length(variables)
        o{i,j}=ephyst_extract(rec,c1,variables(i),a,1); %turns negative values to absolute
    end

end

%% Trimming: Variation over time

%this part finds an index of time-dependent variability of data (in this
%case through least square fitting) and enables trimming by a coefficient
%cutt-off

C=size(o{1,1}); %the use of this vector admits that all cells inside o are the same size
o_mat=cell(length(variables),C(1,2));

cache=o{j};
[m,n]=size(cache);
master_metadata(1:5,2:n+1)=cache(1:5,:);

for j=1:length(o) %this length is given by the number of element in the variable array
     o_temp1=o{j};
     o_temp=o_temp1(5:end,:); %trimming out the headers
     [m,n]=size(o_temp);
     
        for i=1:n
            [mm,~]=size(o_temp(:,i));
            o_temp_cell=cell(mm,1); %make a temporary cell for storing each column of o
            o_temp_cell(:,1)=o_temp(:,i); %store each c of o...
            mat_o_sub=cell2mat(o_temp_cell); %convert to mat, in this case normal matrix
            
            o_mat{j,i}=mat_o_sub; %store the matrix out, for further use
            
            
            xx=(1:length(mat_o_sub));
            c=polyfit(xx',mat_o_sub,1);
            poly_coeff(j,i)=c(1,1); %this array can be used to get the indexes of the "cells" that display abnormal variability
        end
end

clear cache

[m,n]=size(o_mat);
o_mat_x=o_mat;

parameter=parameter_strct(4).value;
coeff_trim=poly_coeff>parameter;

%% for testing the variation coefficients visually
% Cn=C(1,2);
% CM=jet(Cn);
% for j=1:m
%    f=figure;
%    hold on
%    
%     for i=1:Cn
%         cc=numel(o_mat{j,i});
%         cc=(1:cc);
%         o_mat_x{j,i}=cc';
%         scatter(cc',o_mat{j,i},30,CM(i,1:3));
%         l=lsline;
%         set(l(1),'color',CM(i,1:3),'LineWidth',2);
%     end
%     
%     xlabel('event number')
%     ylabel(variables{j})
%     
%     hold off
% end

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

% %exports files and gives names to the csv
% [n,m]=size(o);
% for j=1:m
%    for i=1:n
%        T=cell2table(o{i,j});
%        writetable(T,strcat(char(types(j)),'--',char(variables(i)),'.csv'));
%    end
% end


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
