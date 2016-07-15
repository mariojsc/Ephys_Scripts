%cd the Rec folder with animals
%the code is A# for animal, s# for slice, c# for cell, 'holding current'.
%txt for file

%there are Nan values on top of values in value cells of o and
%such. This is the result of ephyst_extract starting applying str2double
%from the very top of the files, so variable name gets Nan, and for instant
%frequency, also the first N/A

clear all
cd

event_type='excitatory';

switch event_type
    case 'excitatory'
        file_type='-55_E.txt';
    otherwise
        file_type='-55_I.txt';
end
        


begin_event_number=1 %this parameter can be further used, together with parameter for number of analysed events per recording, to define the interval of events to analyse per recording

%% INPUT parameters for statistics

fileID=fopen('~statistics_settings.txt');
S=textscan(fileID,'%s %s');
fclose(fileID);
Statistics_settings_strct=struct('parameter',S{1,1});
values=S{1,2};

for i=1:length(Statistics_settings_strct)
    Statistics_settings_strct(i).argument=values(i);
end

%special case of comparison factors
A=Statistics_settings_strct(end).argument;
Comparison_factors=strsplit(A{1},',');

%special case of number of events
A=Statistics_settings_strct(1).argument;
Statistics_settings_strct(1).argument=str2num(A{1});

clear S A
%% INPUT parameters for trimming

fileID=fopen('~analysis_parameters.txt');
P=textscan(fileID,'%s %f');
fclose(fileID);
Parameter_strct=struct('parameter',P{1,1});
values=P{1,2};

for i=1:length(Parameter_strct)
    Parameter_strct(i).value=values(i);
    Parameter_strct(i).exclusion_vector=[];
end


%% INPUT blinding_correspondence

clear P
Header=1;
fileID=fopen('~blinding_complete.txt');
P=textscan(fileID,'%s %s %s %s %s %s','HeaderLines',Header);
fclose(fileID);

[m,~]=size(P{1,1});

P1=P{1,1};
P2=P{1,2};
P3=P{1,3};
P4=P{1,4};
P5=P{1,5};
P6=P{1,6};

for i=1:m
   Blinding_strct(i).type=P1{i,1};
   Blinding_strct(i).animal=P2{i,1};
   Blinding_strct(i).slice=P3{i,1};
   Blinding_strct(i).cell=P4{i,1};
   Blinding_strct(i).pharmacology=P5{i,1};
   Blinding_strct(i).phenotype=strtrim(P6{i,1});
end

clear  P P1 P2 P3 P4 P5
%% INPUT recording metadata (membrane constants)

Header=1;
Delimiter=',';
fileID=fopen('~M_constant.txt');
C=textscan(fileID,'%s %s %s %s %s','HeaderLines',Header,'Delimiter',Delimiter);
fclose(fileID);

[m,~]=size(C{1,1});
[~,n]=size(C);
Rec_meta=cell(m,n);

C1=C{1,1};
C2=C{1,2};
C3=C{1,3};
C4=C{1,4};
C5=C{1,5};

for i=1:m
   Rec_meta{i,1}=C1{i,1};
   Rec_meta{i,2}=C2{i,1};
   Rec_meta{i,3}=C3{i,1};
   Rec_meta{i,4}=C4{i,1};
   Rec_meta{i,5}=C5{i,1};
end

j=1;
for i=1:m
    
   if strncmpi(Rec_meta{i,1},'A',1)
    animal=Rec_meta{i,1};
   end
   
   if strncmpi(Rec_meta{i,1},'s',1)
    slices=Rec_meta{i,1};
   end
   
   if strncmpi(Rec_meta{i,1},'c',1)
    cells=Rec_meta{i,1};
    c_counter=0;
   end
   
   if isnan(str2double(Rec_meta{i,1}))
       
   else
       c_counter=c_counter+1;
       
       Membrane_C_strct(j).animal=animal;
       Membrane_C_strct(j).slices=slices;
       Membrane_C_strct(j).cells=cells;
       
       Membrane_C_strct(j).Cm=str2double(Rec_meta{i,1});
       Membrane_C_strct(j).Rm=str2double(Rec_meta{i,2});
       Membrane_C_strct(j).Ra=str2double(Rec_meta{i,3});
       Membrane_C_strct(j).HA=str2double(Rec_meta{i,4});
       Membrane_C_strct(j).Comp=str2double(Rec_meta{i,5});
       Membrane_C_strct(j).sample=c_counter;
       
       j=j+1;
   end
   
end

clear Rec_meta C C1 C2 C3 C4 C5 animal slices cell 
%% INPUT recording data sets

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
        %% 
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
        
        files = dir(file_type);
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
% 
% clear cache
% cache=dir;
% 
% v=0;
% for i=1:length(cache)
%     ccache=cache(i).name;
%         if strcmp(ccache,'~blinding.txt')
%         v=v+1;
%         end
% end
% 
% if v==0
% 
%     
% fID=fopen('~blinding.txt','wt');
% fprintf(fID,'type animal slice cell pharmacology experimental_group\n');
% 
% for i=1:length(rec)
%     type=rec(i).type;
%     animal=rec(i).animal;
%     slice=rec(i).slice;
%     cell=rec(i).cell;
%     fprintf(fID,'%s %s %s %s\n',type,animal,slice,cell);
% end
%     fclose(fID);
%     
% end

for i=1:length(rec)
    rec(i).pharmacology=Blinding_strct(i).pharmacology;
    rec(i).phenotype=Blinding_strct(i).phenotype;
end

clear r1 r2 r3 r4 fld1 fld11 fld111
%% Metadata Headers

for i=1:length(rec)
   Metadata_headers{1,i}=rec(i).phenotype;
   Metadata_headers{2,i}=rec(i).animal ;
   Metadata_headers{3,i}=rec(i).slice;
   Metadata_headers{4,i}=rec(i).cell;
   Metadata_headers{5,i}=rec(i).pharmacology ;
   Metadata_headers{6,i}=rec(i).type;
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

%% mini events files organization

types={file_type}; %defines columns of o (add here for more)
variables={'Inst. Freq. (Hz)','Peak Amp (pA)','Event Start Time (ms)'}; %defines rows of o

for j=1:length(types)
    
    c1=char(types(j));    
    for i=1:length(variables)
        o{i,j}=ephyst_extract(rec,c1,variables(i),a,1); %turns negative values to absolute
    end

end

%% Trimming by Membrane constants 

% variation of Ra of more than X%
% left to be done: implement one step removal of cells using max of all 
% exclusion vectors

clear cache

time_of_membrane_constant_measurement=1;
Ra_max=Parameter_strct(1).value;
Ra_variation_max=Parameter_strct(2).value;
Rm_min=Parameter_strct(3).value;
Comp_crit=Parameter_strct(11).value;

j=1;
cc=1;
ii=1;

for i=1:length(Membrane_C_strct) %creates Ra variation, Ra tables of each Ra sample per each cell
    
   
    
    if Membrane_C_strct(i).sample==cc
        variation_Ra{j,cc}=Membrane_C_strct(i).Ra;
        Rm{j,cc}=Membrane_C_strct(i).Rm;
        Comp{j,cc}=Membrane_C_strct(i).Comp;
        cc=cc+1;
    end
    
    if (i+1)<=length(Membrane_C_strct) && Membrane_C_strct(i+1).sample~=cc
        
            %if cc==time_of_membrane_constant_measurement
                Membrane_C_org(ii)=Membrane_C_strct(i);
                ii=ii+1;
            %end

        cc=1;
        j=j+1;
    end
    
end

[~,n]=size(variation_Ra);
[~,m]=size(rec);

for i=1:m %creates exclusion vectors for Ra,Rm,Comp,...
    
        variation_rate_cache=cell2mat(variation_Ra(i,1:n));
        
        for j=1:length(variation_rate_cache)-1
            variation_Ra_full{i,j}=abs((variation_rate_cache(j)-variation_rate_cache(j+1))/variation_rate_cache(j));
        end
        
        variation_rate_m=cell2mat(variation_Ra_full(i,:));
        Rm_mat=cell2mat(Rm(i,1:end));
        Rm_mean{i,1}=mean(Rm_mat);
        Cm_mat=cell2mat(Comp(i,1:end));
        
        if variation_rate_cache(end)>=Ra_max+1
            exclusion_vector_1(i)=1;
        else
            exclusion_vector_1(i)=0;
        end
        
        if variation_rate_m(end)>Ra_variation_max
            exclusion_vector_2(i)=1;
        else
            exclusion_vector_2(i)=0;
        end
        
        if variation_Ra_full{i,1}<Rm_min
            exclusion_vector_3(i)=1;
        else
            exclusion_vector_3(i)=0;
        end
        
       
        if Comp_crit==1 && Cm_mat(end)==1
            exclusion_vector_11(i)=1;
        else
            exclusion_vector_11(i)=0;
        end
end

Parameter_strct(1).exclusion_vector=exclusion_vector_1;
Parameter_strct(2).exclusion_vector=exclusion_vector_2;
Parameter_strct(3).exclusion_vector=exclusion_vector_3;
Parameter_strct(11).exclusion_vector=exclusion_vector_11;

clear variation_rate variation_c



%% Trimming by Dataset properties, before any statistics

%this part finds an index of time-dependent variability of data sets and
%enables trimming by a threshold over this index


%this part converts the cells to mat; the obtained matrices will continue
%to be used

C=size(o{1,1}); %the use of this vector admits that all cells inside have the same number or columns
o_mat=cell(length(variables),C(1,2));
cache=o{1,1};

[~,n]=size(o{1,1});
master_metadata(1:5,2:n+1)=cache(1:5,:);

for j=1:length(o) 
     o_temp1=o{j};
     o_temp=o_temp1(5:end,:); %trimming out the headers
     [~,n]=size(o_temp);
     
        for i=1:n
            mm=length(o_temp(:,i));
            o_temp_cell=cell(mm,1); %make a temporary cell for storing each column of o
            o_temp_cell(:,1)=o_temp(:,i); %store each c of o...
            mat_o_sub=cell2mat(o_temp_cell); %convert to mat, in this case normal matrix
            
            sign_test=sign(mat_o_sub);
            
            if mean(sign_test(3:end,1))<0
                mat_o_sub=mat_o_sub.*-1;
            end
                
            Datasets_all{j,i}=mat_o_sub; %store the matrix out, for further use
           
        end
end

clear cache

min_event_num=Parameter_strct(5).value; %creates exclusion vector for minimum number of events
[~,n]=size(Datasets_all);
cache=Datasets_all;
j=1;
for i=1:n
    
    if numel(cache{1,i})>min_event_num
        trim_cache(:,j)=cache(:,i);
        %trim_metadata(:,j)=master_metadata(:,i+1); %useless
        exclusion_vector_5(i)=0;
        
        j=j+1;
    else
        exclusion_vector_5(i)=1;
    end
    
end

Parameter_strct(5).exclusion_vector=exclusion_vector_5;

j=1;
for i=1:length(exclusion_vector_5)
   
    if exclusion_vector_5(i)==0
        trim1_rec_meta(j)=Membrane_C_strct(i);
        j=j+1;
    end
    
end

%% Statistics and Trimming

clear cache

end_event_number=Parameter_strct(7).value;
Sample_rate=Parameter_strct(10).value;
Min_event_frequency=Parameter_strct(6).value;
%Min_event_amplitude=Parameter_strct().value;
Max_event_frequency=Parameter_strct(8).value;
Max_event_amplitude=Parameter_strct(9).value;

Datasets_all_trimmed_1=Datasets_all;
[m,n]=size(Datasets_all_trimmed_1);

for j=1:m

   for i=1:n
    cache=Datasets_all_trimmed_1{j,i};

    if length(cache)>=end_event_number+1;
    cache=cache(begin_event_number+1:end_event_number+1);
    end

    Datasets_all_trimmed_1{j,i}=cache;
   end

end

%clear cache cache_T
for i=1:n
    
    cache_T=Datasets_all_trimmed_1{3,i};
    cache_A=Datasets_all_trimmed_1{2,i};
    mean_amp(i)=mean(cache_A(2:end));
    deltaT_events(i)=cache_T(end)-cache_T(2);
    frequency_events_cells(i)=end_event_number/(deltaT_events(i)*10^-4);
    
end

exclusion_vector_6=frequency_events_cells<Min_event_frequency;
exclusion_vector_8=frequency_events_cells>Max_event_frequency;
exclusion_vector_9=mean_amp>Max_event_amplitude;

Parameter_strct(6).exclusion_vector=exclusion_vector_6;
Parameter_strct(8).exclusion_vector=exclusion_vector_8;
Parameter_strct(9).exclusion_vector=exclusion_vector_9;

%% for testing variation index visually
% clear cache
% 
% [m,n]=size(trim_cache);
% CM=jet(n);
% 
% %axis_var=[0 70;5 45];
% axis_var=[-1 2;-2 2]; %for zscore
% 
% for j=1:m
%    f=figure;
%    hold on
%    
%     for i=1:n
%         cache=trim_cache{j,i};
%         cc=numel(cache);
%         cc=(3:cc);
%         o_mat_x{j,i}=cc';
%         
%         zeed=zscore(cache(3:end,1));
%         moving_av=tsmovavg(zeed,'t',40,1);
%         
%         %moving_av=tsmovavg(cache,'t',40,1);
%         
%         scatter(cc',zeed,1,CM(i,1:3));
%         
%         %scatter(cc',cache,1,CM(i,1:3));
%         
%         plot(moving_av,'Color',CM(i,1:3),'LineWidth',0.8);
%         
%         %set(l(1),'color',CM(i,1:3),'LineWidth',2);
%     end
%     
%     axis([0 600 axis_var(j,1) axis_var(j,2)])
%     
%     xlabel('event number')
%     ylabel(variables{j})
%     line([0,600],[0,0],'Color','k')
%     fig=gcf;
%     fig.Position=([1 1 1440 275]);
%     
%     hold off
% end

%% Final Trimming


Results_cells(1,1:length(frequency_events_cells))=num2cell(frequency_events_cells);
Results_cells(2,1:length(frequency_events_cells))=num2cell(mean_amp);
exc_sz=Parameter_strct(1).exclusion_vector;

for i=1:length(Parameter_strct)
    if ~isempty(Parameter_strct(i).exclusion_vector)
        P(i,1:length(exc_sz))=Parameter_strct(i).exclusion_vector;
    end
end

j=1;
for i=1:length(frequency_events_cells)
    
   P_cache=sum(P(1:end,i));
   
   if P_cache==0
    Results_cells_trim_1(:,j)=Results_cells(:,i);
    Metadata_trim(:,j)=Metadata_headers(:,i);
    Datasets_trimmed_2(:,j)=Datasets_all_trimmed_1(:,i);
    Membrane_C_org_trim(j)=Membrane_C_org(i);
    j=j+1;
   end
   
end

[m,n]=size(Results_cells_trim_1);

for i=1:n

    for j=1:m
        Results_cells_trim_1{m+j,i}=Datasets_trimmed_2{j,i};
    end
    
    Results_cells_trim_1{m+j+1,i}=Membrane_C_org_trim(i).Cm;
    Results_cells_trim_1{m+j+2,i}=Membrane_C_org_trim(i).Rm;
    Results_cells_trim_1{m+j+3,i}=Membrane_C_org_trim(i).Ra;
    Results_cells_trim_1{m+j+4,i}=Membrane_C_org_trim(i).HA;
    Results_cells_trim_1{m+j+5,i}=Membrane_C_org_trim(i).Comp;
    
end


%% Final sorting

organizer_index=1;

for i=1:length(Comparison_factors)
    
    C_f=Comparison_factors{1,i};
    switch C_f
        case 'phenotype'
            C_f_index=1;
        case 'animal'
            C_f_index=2;
        case 'slice'
            C_f_index=3;
        case 'cells'
            C_f_index=4;
        case 'pharmacology'
            C_f_index=5;
        case 'type'
            C_f_index=6;
    end
    
    C_f_possible{i,1}=unique(Metadata_trim(C_f_index,:));
    Filler=length(Metadata_trim(:,1));
    Array_possible=C_f_possible{i,1};
    
    for jjj=1:length(Results_cells_trim_1(:,1))
            
        Final_results_cache=cell(length(Results_cells_trim_1)+1,length(C_f_possible{i,1})+length(Metadata_trim(:,1)));
        for j=1:length(Array_possible);
            Final_results_cache{1,length(Metadata_trim(:,1))+j}=Array_possible{1,j};
        end

        inner_j=2;

              for jj=1:length(Array_possible)

                for j=1:length(Results_cells_trim_1(1,:))

                   if strcmp(Metadata_trim{C_f_index,j},Array_possible{jj})

                        Final_results_cache{inner_j,1}=Metadata_trim{1,j};
                        Final_results_cache{inner_j,2}=Metadata_trim{2,j};
                        Final_results_cache{inner_j,3}=Metadata_trim{3,j};
                        Final_results_cache{inner_j,4}=Metadata_trim{4,j};
                        Final_results_cache{inner_j,5}=Metadata_trim{5,j};
                        Final_results_cache{inner_j,6}=Metadata_trim{6,j};

                        Final_results_cache{inner_j,Filler+jj}=Results_cells_trim_1{jjj,j};

                        inner_j=inner_j+1;
                   end

                end

              end

                Final_results{jjj,organizer_index}=Final_results_cache;
                inner_j=1;
    end
      
organizer_index=organizer_index+1;
        
end



%% OUTPUT

[m,n]=size(Datasets_all_trimmed_1); %exclusion report
M=length(Parameter_strct);
[M2,N2]=size(Metadata_headers);

for i=1:M
exclusion_report{i+M2,1}=Parameter_strct(i).parameter;
end

exclusion_report(2:M2+1,2:N2+1)=Metadata_headers;

clear cache
for i=1:M
    cache=Parameter_strct(i).exclusion_vector;
    if ~isempty(cache)
    cache=num2cell(cache)';
    exclusion_report(i+M2,2:end)=cache;
    end
end

E=cell2table(exclusion_report);
writetable(E,'~Exclusion_report.csv');

Variables={'mean_frequency','mean_amplitude','every_mini_frequency','every_mini_amplitude','Cm','Rm','Ra','HA','Comp'};

[n,m]=size(Final_results);
 for j=1:m
    for i=1:n
        T=cell2table(Final_results{i,j});
        writetable(T,char(strcat(Variables(i),'--',Comparison_factors(j),'.csv')));
    end
 end
 