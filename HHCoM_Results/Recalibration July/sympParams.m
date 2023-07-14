% attempt 1 (July 7)

sympParams = []; 

for l = 0.0001 : 0.00002 : 0.0002
    for r = 0.014 : 0.002 : 0.02
        for d = 0.6 : 0.1 : 0.8 
            for l_r = 0.02 : 0.02 : 0.1
                for r_d = 0.025 : 0.02 : 0.1

                    sympParams = [sympParams; l r d l_r r_d]; 

                end
            end
        end 
    end 
end 


% attempt 2 (July 12)

sympParams = []; 

for l = 0.001 : 0.002 : 0.01 % changed this
    for r = 0.014 : 0.002 : 0.02
        for d = 0.6 : 0.1 : 0.8 
            for l_r = 0.06 : 0.02 : 0.14
                for r_d = 0.025 : 0.02 : 0.1

                    sympParams = [sympParams; l r d l_r r_d]; 

                end
            end
        end 
    end 
end 

size(sympParams)

% attempt 3 (July 13)

sympParams = []; 

for l = 0.001 : 0.002 : 0.01 % changed this
    for r = 0.014 : 0.002 : 0.02
        for d = 0.6 : 0.1 : 0.8 
            for l_r = 0.4 : 0.1 : 0.8
                for r_d = 0.025 : 0.02 : 0.1

                    sympParams = [sympParams; l r d l_r r_d]; 

                end
            end
        end 
    end 
end 

size(sympParams)