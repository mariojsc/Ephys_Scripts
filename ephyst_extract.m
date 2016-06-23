function [target]=ephyst_extract(struct_origin,type,variable,target_size,modifier)   

%given an original structure, the function looks for a type of file,
%variables in those files and turns them into isolated table of target_size

    target=cell(max(target_size)+10,length(struct_origin));
    
    i=1;
    
    for y=1:length(struct_origin)
        
        c2=char(struct_origin(y).type);
        
        if strcmp(c2,type)

            target{1,i}=struct_origin(y).animal;    
            target{2,i}=struct_origin(y).slice;    
            target{3,i}=struct_origin(y).type;
            target{4,i}=struct_origin(y).cell;
            target_cache=table2cell(struct_origin(y).values);
            
            if iscell(variable) %should be modified for all possibilities; this makes the function accept online str inside cells for this argument!!!
                [~,n]=find(strcmp(variable,target_cache));
                [l,~]=size(target_cache);
                k=1;
                    while k<=l
                        target{k+4,i}=str2double(target_cache{k,n}); %k+1 for the column title 
                        k=k+1;
        %             if target{k+4,i}<0
        %                target{k+4,i}=target{k+4,i}*-1;
        %             end
                    end
               
            elseif variable==0
                 n=1;
                [l,~]=size(target_cache);
                k=1;
                    while k<=l
                        target{k+4,i}=target_cache{k,n}; %k+1 for the column title 
                        k=k+1;
                    end
          
            end
            
            
            
        end
        i=i+1;
    end
    
        if modifier==0 % in order to pic single strings from files (as in the CC.txt files)
                target2=cell(size(target).*3);
                [~,m]=size(target);
                i_hold=1;
                j=1;
                i=1;
                j_hold=j;
                    while j<m
                        a=target{1,j};
                        target2{1,j_hold}=a;
                        i=1;
                            while target{1,j}==a 
                                target2{i+1,j_hold}=target{5,j};
                                i=i+1;
                                j=j+1;
                                
                                if j>m
                                    break
                                end
                                 
                            end
                        j_hold=j_hold+1;
                    end
                target=target2;
        end
    
end