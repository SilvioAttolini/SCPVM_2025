function [indexes] = correctindexes(maxorder)
counter_i = 0;
counter = 0;
for l_order = 0:maxorder
    for u_order = -l_order:l_order
        counter_i = counter_i+1;
        if mod(l_order+abs(u_order),2)==0
            counter = counter+1;
            save = counter_i;
            indexes(counter) = counter_i;
        end
        
    end
end
end