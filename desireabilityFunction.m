function d = desireabilityFunction(y, m, I, effort, effortThreshold, Imax, r, type)

    if type == "smooth"
        T = sqrt(1);
        L = sqrt(.5); 
    else
        T = 1;
        L = .5;
    end
    
    if effort < effortThreshold && m <= Imax && m > 0 && I <= Imax && I > 0 && y >= L && y <= T
        d = ((y-L)/(T-L))^r; 
    else
        d = 0;
    end

end