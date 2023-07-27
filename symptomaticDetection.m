function[dPop, ccSymp] = symptomaticDetection(pop , ...
    year , hpvScreenStartYear , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , risk , intervens, age, ...
    screenAlgs , kSymp, hystMult, udPop, udPopNoTreat, udPopTreat, udPopHyst, ccReghpvVaxInds)

% Note that this function is only run during the years before the screening year. 
% There should be no hysterectomy in these years, so I am not including any symptomatic treatment with hysterectomy. 

dPop = zeros(size(pop));
ccSymp = zeros(3, age, 3);
% for ccSymp, the last compartment represents (1) symptomatic/treated, (2) symptomatic/untreated, (3) symptomatic/treated with hysterectomy

if (year < hpvScreenStartYear)
    for d = 1 : disease 
	    for v = 1 : viral
		    for h = 1 : hpvVaxStates
			    for s = 1 : hpvNonVaxStates
				    for x = 1 : 3
					    for i = 1 : intervens
						    for a = 1 : age 
							    for r = 1 : risk
                                    if [( (x==1) && ((h==6) || (s==6)) ) || (x==2) || (x==3)] 
                             		    toSympTreat = screenAlgs.treatRetain(4); % proportion of women who move to treatment from symptoms
                                        % toSympTreatHyst = screenAlgs.colpoRetain * screenAlgs.treatRetain(4) * hystMult(x); % removed since there should be no hysterectomy
    %                                     toSympUntreat = screenAlgs.colpoRetain * (1 - screenAlgs.treatRetain(4)); 
                                        toSympUntreat = 1 - toSympTreat; % anyone who is not treated / LTFU is moved to untreated
                                        sympDetected = kSymp(x) .* pop(udPop(d,v,h,s,x,i,a,r)); % num women with symptoms
    
                                        dPop(udPopTreat(d,v,h,s,x,i,a,r)) = dPop(udPopTreat(d,v,h,s,x,i,a,r)) + toSympTreat .* sympDetected; % move to treated compartment
                                        dPop(udPopNoTreat(d,v,h,s,x,i,a,r)) = dPop(udPopNoTreat(d,v,h,s,x,i,a,r)) + toSympUntreat .* sympDetected; % move to untreated compartment
                                        % dPop(udPopHyst(d,v,h,s,x,i,a,r)) = dPop(udPopHyst(d,v,h,s,x,i,a,r)) + toSympTreatHyst .* sympDetected; % move to hyst compartment. removed since there should be no treatment by hysterectomy before the screening year. 
    
                                        dPop(udPop(d,v,h,s,x,i,a,r)) = dPop(udPop(d,v,h,s,x,i,a,r)) - sympDetected; % remove from the undetected compartment 
    
                                        % update results matrices
                                        ccSymp(x,a,1) = ccSymp(x,a,1) + sum(toSympTreat * sympDetected); %treated
                                        ccSymp(x,a,2) = ccSymp(x,a,2) + sum(toSympUntreat * sympDetected); %untreated
                                        % ccSymp(d,h,s,x,i,a,3) = ccSymp(d,h,s,x,i,a,3) + sum(toSympTreatHyst * sympDetected); %hysterectomy. removed since there should be no treatment by hysterectomy. 
								    end 
							    end 
						    end 
					    end 
				    end 
			    end 
		    end
	    end
    end 
end 

if (year >= hpvScreenStartYear)

  for d = 1 : disease 
	for v = 1 : viral
		for h = 1 : hpvVaxStates
			for s = 1 : hpvNonVaxStates
				for x = 1 : 3
					for i = 1 : intervens
						for a = 1 : age 
							for r = 1 : risk
                                if [( (x==1) && ((h==6) || (s==6)) ) || (x==2) || (x==3)] 

                                     toSympTreat = screenAlgs.treatRetain(4) * (1 - hystMult(x)); % proportion of women who move to treatment from symptoms
                                     toSympTreatHyst = screenAlgs.treatRetain(4) * hystMult(x); 
                                     toSympUntreat = 1 - (toSympTreatHyst + toSympTreat); % anyone who is not treated / LTFU is moved to untreated
                                    
%                              		toSympTreat = screenAlgs.treatRetain(4); % proportion of women who move to treatment from symptoms
%                                     % toSympTreatHyst = screenAlgs.colpoRetain * screenAlgs.treatRetain(4) * hystMult(x); % removed since there should be no hysterectomy
% %                                     toSympUntreat = screenAlgs.colpoRetain * (1 - screenAlgs.treatRetain(4)); 
%                                     toSympUntreat = 1 - toSympTreat; % anyone who is not treated / LTFU is moved to untreated
                                    sympDetected = kSymp(x) .* pop(udPop(d,v,h,s,x,i,a,r)); % num women with symptoms

                                    dPop(udPopTreat(d,v,h,s,x,i,a,r)) = dPop(udPopTreat(d,v,h,s,x,i,a,r)) + toSympTreat .* sympDetected; % move to treated compartment
                                    dPop(udPopNoTreat(d,v,h,s,x,i,a,r)) = dPop(udPopNoTreat(d,v,h,s,x,i,a,r)) + toSympUntreat .* sympDetected; % move to untreated compartment
                                    dPop(udPopHyst(d,v,h,s,x,i,a,r)) = dPop(udPopHyst(d,v,h,s,x,i,a,r)) + toSympTreatHyst .* sympDetected; % move to hyst compartment. removed since there should be no treatment by hysterectomy before the screening year. 

                                    dPop(udPop(d,v,h,s,x,i,a,r)) = dPop(udPop(d,v,h,s,x,i,a,r)) - sympDetected; % remove from the undetected compartment 

                                    % update results matrices
                                    ccSymp(x,a,1) = ccSymp(x,a,1) + sum(toSympTreat * sympDetected); %treated
                                    ccSymp(x,a,2) = ccSymp(x,a,2) + sum(toSympUntreat * sympDetected); %untreated
                                    ccSymp(x,a,3) = ccSymp(x,a,3) + sum(toSympTreatHyst * sympDetected); %hysterectomy
								end 
							end 
						end 
					end 
				end 
			end 
		end
	end
  end 
end 


dPop = sparse(dPop);