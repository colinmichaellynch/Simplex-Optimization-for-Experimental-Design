clear

%% Load data, set constants

data = load("TrueResponseSurfaceTest.mat");

IMatrix = data.IMatrix;
IVec = IMatrix(:); 
mMatrix = data.mMatrix;
mVec = mMatrix(:); 
effortMatrix = data.effortMatrix;
effortVec = effortMatrix(:); 
accuracyMatrix = data.accuracyMatrix;
accuracyVec = accuracyMatrix(:); 

simulations = 1000; 
effortThreshold = 400; 
alpha = 1;
gamma = 2;
rho = .5;
sigma = .5; 
iterations = 20;
Imax = 40; 
Imin = 1; 
r = 2;
W = 1/3;
type = "smooth";
pointsAreCollinear = @(xy) rank(xy(2:end,:) - xy(1,:)) == 1;

xGraph = 1:Imax; 
yGraph = (effortThreshold-(W-1).*xGraph)./(xGraph + (1./W)- 1); 

%% Find point with minimal effort with BA at least at threshold

boolean = effortMatrix > effortThreshold; 
baFittedAdj = accuracyMatrix;
baFittedAdj(boolean) = NaN; 
maximum = max(max(sqrt(baFittedAdj))); 
[maxx,maxy] = find(sqrt(baFittedAdj)==maximum);

if length(maxx) > 1
    effort = maxx.*maxy+(W-1).*maxy+((1/W)-1).*maxx;
    [~, index] = min(effort); 
    maxx = maxx(index);
    maxy = maxy(index);
end

maxIFitted = IMatrix(maxx,maxy); 
maxmFitted = mMatrix(maxx,maxy); 
finalBAFitted = baFittedAdj(maxx,maxy); 
finalEffortFitted = effortMatrix(maxx,maxy); 

%% Combine response surfaces into desirability function 

tbl = table(IVec,mVec,effortVec,sqrt(accuracyVec),'VariableNames',{'I','m','Effort','BA'});
model = 'BA ~ I + m + I:m + I^2 + m^2 + I^3 + m^3 + I^4 + m^4'; 
lm = fitlm(tbl,model); 
baFitted = reshape(lm.Fitted, 40, 40); 
Coef = table2array(lm.Coefficients);
Coef = Coef(:,1);

for I = 1:Imax
    for m = 1:Imax
        effort = I*m+(W-1)*m+((1/W)-1)*I;
        y = baFitted(I, m); 
        desirability(I, m) = desireabilityFunction(y, m, I, effort, effortThreshold, Imax, r, type);
    end
end

maximum = max(max(desirability)); 
[maxx,maxy] = find(desirability==maximum);

%% Find area of effort threshold

counter = 1; 
vertexX = [];
vertexY = [];
maximumI = (effortThreshold-(W-1)*1)/(1 + (1./W)- 1); 
maximumM = (effortThreshold*W - 1 + 1*W)/(W*(1-1+W)); 
if Imax < maximumI
    vertexX(counter) = 1; 
    vertexY(counter) = Imax; 
    counter = counter+1;
    vertexX(counter) = (effortThreshold*W - Imax + 1*W)/(W*(Imax-1+W)); 
    vertexY(counter) = Imax; 
else
    vertexX(counter) = 1; 
    vertexY(counter) = maximumI; 
end
counter = counter+1;
if Imax < maximumM
    vertexX(counter) = Imax; 
    vertexY(counter) = (effortThreshold-(W-1)*Imax)/(Imax + (1/W)- 1); 
    counter = counter+1;
    vertexX(counter) = Imax; 
    vertexY(counter) = 1; 
else
    vertexX(counter) = maximumM; 
    vertexY(counter) = 1; 
end
counter = counter+1;
vertexX(counter) = 1; 
vertexY(counter) = 1; 
vertices = [vertexX', vertexY'];
[k,area] = convhull(vertices);
minAreaofTriangle = .1*area;
maxAreaofTriangle =.3*area;

%% Implement search method

for h = 1:simulations

    flag = 0;
    counterLoop = 0;

    %make sure that the initial triangle is large enough
    while flag == 0

        x = randi([Imin Imax],1,3);
        y = randi([Imin Imax],1,3);
        counter = 1;
        while counter < 4

            effort = x(counter)*y(counter) + (W-1)*x(counter) + (1/W - 1)*y(counter);
            if effort > effortThreshold
                x(counter) = randi([Imin Imax]);
                y(counter) = randi([Imin Imax]);
            else
                counter = counter+1;
            end
            vertices = [x', y'];
        end

        %if sum(diff(x)) > 0 && sum(diff(y)) > 0 
        if pointsAreCollinear(vertices)==0
            [~, area] = convhull(vertices);

            if area > minAreaofTriangle && area < maxAreaofTriangle
                flag = 1;
            end

            counterLoop = counterLoop+1;
            if counterLoop>1000
                flag = 1;
            end

        end
    end

    randX = x;
    randY = y;

    for j = 1:iterations

        %step 1: order
        for i = 1:3 
            y = mdlPred(Coef, randX(i), randY(i));
            effort = effortFunction(randX(i), randY(i), W);
            z(i) = 1-desireabilityFunction(y, randY(i), randX(i), effort, effortThreshold, Imax, r, type); %find zs
        end
        [z, rank] = sort(z, 'ascend'); %sort them in ascending order
        randX = randX(rank);
        randY = randY(rank);

        %step 2: centroid
        x0 = mean(randX(1:end-1)); 
        y0 = mean(randY(1:end-1)); 

        %step 3: reflection
        xr = x0+alpha*(x0 - randX(end));
        yr = y0+alpha*(y0 - randY(end));
        y = mdlPred(Coef, xr, yr);
        effort = effortFunction(xr, yr, W);
        zr = 1-desireabilityFunction(y, yr, xr, effort, effortThreshold, Imax, r, type); 
        boolean = zr < z; 

        if sum(boolean)==2
            randX(end) = xr;
            randY(end) = yr;
            %return to step 1

        elseif sum(boolean) == 3 

            %step 4: expansion
            xe = x0 + gamma*(xr - x0);
            ye = y0 + gamma*(yr - y0);
            y = mdlPred(Coef, xe, ye);        
            effort = effortFunction(xe, ye, W);
            ze = 1-desireabilityFunction(y, ye, xe, effort, effortThreshold, Imax, r, type); 

            if ze < zr
                randX(end) = xe;
                randY(end) = ye;
                %return to step 1
            else
                randX(end) = xr;
                randY(end) = yr;
                %return to step 1
            end

        else

            %step 5: contraction
            xc = x0 + rho*(randX(end) - x0);
            yc = y0 + rho*(randY(end) - y0);
            y = mdlPred(Coef, xc, yc);        
            effort = effortFunction(xc, yc, W);
            zc = 1-desireabilityFunction(y, yc, xc, effort, effortThreshold, Imax, r, type); 

            if zc < z(end)
                randX(end) = xc;
                randY(end) = yc;
                %return to step 1
            else

                %step 6: shrink
                randX(2) = randX(1) + sigma*(randX(2)-randX(1));
                randX(3) = randX(1) + sigma*(randX(3)-randX(1));
                randY(2) = randY(1) + sigma*(randY(2)-randY(1));
                randY(3) = randY(1) + sigma*(randY(3)-randY(1));
            end

        end

        perimeter = sqrt((randX(1)-randX(2))^2+(randY(1)-randY(2))^2)+sqrt((randX(1)-randX(3))^2+(randY(1)-randY(3))^2)+sqrt((randX(2)-randX(3))^2+(randY(2)-randY(3))^2);

        if perimeter < 2 + sqrt(2) + 1
            break
        end

    end

    finalX = round(mean(randX));
    finalY = round(mean(randY));
    finalZ = accuracyMatrix(finalX, finalY);

    numberOfComputations(h) = 2+j; 
    if (1-zc) > maximum-.025*maximum && (1-zc) < maximum+.025*maximum
        reachedPeak(h) = 1;
    else
        reachedPeak(h) = 0;
    end
    
    finalDesirability(h) = (1-zc);
    if (1-zc)>maximum
        finalDesirability(h) = maximum;
    end
end

figure(1)
tiledlayout(1,2)

nexttile
boxplot(numberOfComputations, 'Labels', {'Simplex Method'})
ylabel("Number of Computations", 'FontSize', 14)
ttl = title('A)', 'FontSize', 18);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  

nexttile
boxplot(finalDesirability, 'Labels', {'Simplex Method'});
xl =yline(maximum,'-.','Maximum');
ylim([.8, 1])
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'left';
ttl = title('B)', 'FontSize', 18);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylabel("Final Desirability", 'FontSize', 14)

sum(reachedPeak)/simulations
median(numberOfComputations)/(Imax^2-1)