function [] = MultifidBayesOpt(fun,opt)

if opt.numFidels == 2


    f2 = fun.fid1;
    f1 = fun.fid2;


    for test = 1: opt.numTests
        tic
        fprintf('Test number = %d \n', test);
        Chi = round(unscale_point(lhsdesign(opt.grid_size, opt.dims), opt.mins, opt.maxes), 4);
        grid = Chi(floor(rand(opt.numSamp, 1)*opt.grid_size),:);


        minsamples = [];
        minvalues = [];
        samples = grid;



        X1 = [];
        X2 = [];
        y1 = [];
        y2 = [];



        fprintf('Running first samples...\n');

        for j = 1:opt.numSampFid(1)
            y1(end+1,1) = f1(grid(j,:));
            X1(end+1,:) = grid(j,:);
        end

        for j = opt.numSampFid(1)+1:opt.numSampFid(1)+opt.numSampFid(2)
            y2(end+1,1) = f2(grid(j,:));
            X2(end+1,:) = grid(j,:);
        end


        [minvalues(1,1),minpos] = min(y2);
        minsamples(1,:) = X2(minpos, :);


        [ia, ~] = ismember(Chi, grid, 'rows');
        Chi(ia, :) = [];





        Iteration = 0;
        times = [];

        Budget = opt.numSampFid(end)  + opt.FidelityCost(1)*(opt.numSampFid(1));

        for i = 1:opt.max_iters
            tic


            gprMdl1 = fitrgp(X1,y1,'KernelFunction','squaredexponential');
            [ypred1,ysd1] = predict(gprMdl1, Chi);
            yvar1 = (ysd1).^2;
            [ypred2, ysd2] =  MFGP(f1, f2, X1, X2, y1, y2, Chi, opt);
            yvar2 = (ysd2).^2;


            mu = [ypred1 ypred2];
            sigma2 = [yvar1 yvar2];



            if opt.AF == 1
                [x_new, M, Chi] = MFEI(mu, sigma2, Chi, y2, opt);
            end


            if opt.AF == 2
                [x_new, M, Chi] = MFPI(mu, sigma2, Chi, y2, X1, X2, opt);
            end



            if M == 2
                Budget = Budget + 1;
                X2(end+1,:) = x_new;
                y2(end+1) = f2(x_new);
            end

            if M == 1
                Budget = Budget + opt.FidelityCost(1);
                X1(end+1,:) = x_new;
                y1(end+1) = f1(x_new);
            end



            fprintf('Iteration = %d, Budget = %d, Minimum = %f\n', i, Budget, minvalues(end));

            [minvalues(end +1,1), minposition] = min(y2);
            minsamples(end+1,:) = X2(minposition,:);
            Cost(i, 1) = Budget;
            samples(end+1,:) = x_new;
            Iteration = Iteration +1;
            save(['Test' num2str(test)])

            if Budget >= opt.Budget
                break
            end

            timeiter = toc;
            times(end+1) = timeiter;
        end

    end

end


end



