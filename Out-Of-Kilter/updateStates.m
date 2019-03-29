function kilter_state = updateStates(kilter_state,Xij,omegas,Lij,Uij,Cij,Arches)
%% Update kilter states for all our arches according to in-kilter conditions
% João Araújo, 2019
    
for i = 1:size(kilter_state,1)

    % In-kilter condition: cij - omega(i) + omega(j) > 0
    if Cij(i) - omegas(Arches(i,1)) + omegas(Arches(i,2)) > 0
        if Xij(i) == Lij(i) 
            kilter_state(i) = true;
        end

    % In-kilter condition: cij - omega(i) + omega(j) = 0
    elseif Cij(i) - omegas(Arches(i,1)) + omegas(Arches(i,2)) == 0

        if Xij(i) == Uij(i) || ...
           (Xij(i) > Lij(i) && Xij(i) < Uij(i)) || ...
           Xij(i) == Lij(i)

           kilter_state(i) = true; 

        end

    % In-kilter condition: cij - omega(i) + omega(j) < 0
    elseif Cij(i) - omegas(Arches(i,1)) + omegas(Arches(i,2)) < 0

        if Xij(i) == Uij(i) 
            kilter_state(i) = true; 
        end

    end

end

end