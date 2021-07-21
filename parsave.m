function parsave(fname, fivYrAgeGrpsOn , tVec , popVec , newHiv , transCD4 , ...
    hivDeaths , deaths , ccDeath , ...
    menCirc , ...
    artDist , artDistList , artTreatTracker , artDiscont , ...
    nHivDiagVec , nHivUndiagVec , nTestedNeg , nTestedUndiag , propHivDiag , ...
    currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier)

savDir = [pwd , '/HHCoM_Results/' , pathModifier, '/'];
save(fullfile(savDir , fname) , 'fivYrAgeGrpsOn' , 'tVec' ,  'popVec' , ...
    'newHiv' , 'transCD4' , ...
    'hivDeaths' , 'deaths' , 'ccDeath' , ...
    'menCirc' , ...
    'artDist' , 'artDistList' , 'artTreatTracker' , 'artDiscont' , ...
    'nHivDiagVec' , 'nHivUndiagVec' , 'nTestedNeg' , 'nTestedUndiag' , 'propHivDiag' , ...
    'currYear' , 'lastYear' , 'vaxRate' , 'vaxEff' , 'popLast' , '-v7.3') 
end
