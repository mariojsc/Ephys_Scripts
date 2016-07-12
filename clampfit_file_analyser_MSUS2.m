%cd the Rec folder with animals
%the code is A# for animal, s# for slice, c# for cell, 'holding current'.
%txt for file

%there are Nan values on top of values in value cells of o and
%such. This is the result of ephyst_extract starting applying str2double
%from the very top of the files, so variable name gets Nan, and for instant
%frequency, also the first N/A

cd

%% INPUT parameters for trimming

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
P=textscan(fileID,'%s %s %s %s %s','HeaderLines',Header);
fclose(fileID);

[m,~]=size(P{1,1});
[~,n]=size(P);

P1=P{1,1};
P2=P{1,2};
P3=P{1,3};
P4=P{1,4};
P5=P{1,5};

for i=1:m
   Blinding_strct(i).type=P1{i,1};
   Blinding_strct(i).animal=P2{i,1};
   Blinding_strct(i).slice=P3{i,1};
   Blinding_strct(i).cell=P4{i,1};
   Blinding_strct(i).treatment=P5{i,1};
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

clear r1 r2 r3 r4 fld1 fld11 fld111
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

types={'-55_E.txt'}; %defines columns of o (add here for more)
variables={'Inst. Freq. (Hz)','Peak Amp (pA)'}; %defines rows of o

for j=1:length(types)
    
    c1=char(types(j));    
    for i=1:length(variables)
        o{i,j}=ephyst_extract(rec,c1,variables(i),a,1); %turns negative values to absolute
    end

end

%% Trimming

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
            
%           xx=(1:length(mat_o_sub));
%           c=polyfit(xx',mat_o_sub,1);
%           poly_coeff(j,i)=c(1,1); %this array can be used to get the indexes of the "cells" that display abnormal variability

        end
end

%% Trimming by Membrane constants 

% variation of Ra of more than X%
% left to be done: implement one step removal of cells using max of all 
% exclusion vectors

clear cache

Ra_variation_max=Parameter_strct(2).value;

j=1;
cc=1;

for i=1:length(Membrane_C_strct)
    
    if Membrane_C_strct(i).sample==cc
        variation_c{j,cc}=Membrane_C_strct(i).Ra;
        cc=cc+1;
    end
    
    if (i+1)<=length(Membrane_C_strct) && Membrane_C_strct(i+1).sample~=cc
        cc=1;
        j=j+1;
    end
    
end

[m,n]=size(variation_c);

for i=1:m
    
        variation_rate_cache=cell2mat(variation_c(i,1:n));
        
        for j=1:length(variation_rate_cache)-1
            variation_rate{i,j}=abs((variation_rate_cache(j)-variation_rate_cache(j+1))/variation_rate_cache(j));
        end
        
        variation_rate_m=cell2mat(variation_rate(i,:));
        
        if max(variation_rate_m)>Ra_variation_max
            exclusion_vector(i)=1;
        else
            exclusion_vector(i)=0;
        end
end

Parameter_strct(2).exclusion_vector=exclusion_vector;

%% first trimming, by minimum number of events

clear cache

min_event_num=Parameter_strct(5).value;
[m,n]=size(Datasets_all);
cache=Datasets_all;

j=1;
for i=1:n
    if numel(cache{1,i})>30
        mat_o_trim1(:,j)=cache(:,i);
        trim1_master_metadata(:,j)=master_metadata(:,i+1);
        exclusion_vector(i)=0;
        
        j=j+1;
    else
        exclusion_vector(i)=1;
    end
end
Parameter_strct(5).exclusion_vector=exclusion_vector;

j=1;
for i=1:length(exclusion_vector)

    if exclusion_vector(i)==0
        trim1_rec_meta(j)=Membrane_C_strct(i);
        j=j+1;
    end
end

% to insert trimming criteria
% parameter=parameter_strct(4).value;
% coeff_trim=poly_coeff>parameter;


%% for testing variation index visually
clear cache

[m,n]=size(mat_o_trim1);
CM=jet(n);

%axis_var=[0 70;5 45];
axis_var=[-1 2;-2 2]; %for zscore

for j=1:m
   f=figure;
   hold on
   
    for i=1:n
        cache=mat_o_trim1{j,i};
        cc=numel(cache);
        cc=(3:cc);
        o_mat_x{j,i}=cc';
        
        zeed=zscore(cache(3:end,1));
        moving_av=tsmovavg(zeed,'t',40,1);
        
        %moving_av=tsmovavg(cache,'t',40,1);
        
        scatter(cc',zeed,1,CM(i,1:3));
        
        %scatter(cc',cache,1,CM(i,1:3));
        
        plot(moving_av,'Color',CM(i,1:3),'LineWidth',0.8);
        
        %set(l(1),'color',CM(i,1:3),'LineWidth',2);
    end
    
    axis([0 600 axis_var(j,1) axis_var(j,2)])
    
    xlabel('event number')
    ylabel(variables{j})
    line([0,600],[0,0],'Color','k')
    fig=gcf;
    fig.Position=([1 1 1440 275]);
    
    hold off
end

%% Statistics


%% OUTPUT


 [n,m]=size(o);
 for j=1:m
    for i=1:n
        T=cell2table(o{i,j});
        writetable(T,strcat(char(types(j)),'--',char(variables(i)),'.csv'));
    end
 end