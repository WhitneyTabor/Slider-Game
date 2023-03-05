% run31 local gradient following with random position switches in continuous time
% (parameter switchtime determines how often players randomly switch)

% RESULT:  Good fraughtness and fairly good fit with inc = 0.05, 
% switchtime = 10, fsought = 0.3;  unlike the humans and APEM, this model
% shows sympathetic bowing in the Right Branch during the adverse.
%  Also fsought has to be pretty high to make the fit strong.


% Saved a in Hierarchy directory for inclusion in 2023 PLOSone
% submission

tic

addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/CellArrays';
addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/FileIO';
addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/PiecewiseLinear';
addpath '/Users/wht02001/Dropbox/Projects/Solab Dropbox/Tools/Solab Matlab Library/Fraughtness';
addpath '/Users/wht02001/Dropbox/Projects/Solab Dropbox/Tools/Solab Matlab Library/SliderGame';

trainsound = load('train');
chirpsound = load('chirp');
handelsound = load('handel');
gongsound = load('gong');

toc

%% Run model

% Set parameters

% Model params
adjacency = [1 1 1 0 0 0 0;
    1 1 0 1 1 0 0;
    1 0 1 0 0 1 1;
    0 1 0 1 0 0 0;
    0 1 0 0 1 0 0;
    0 0 1 0 0 1 0;
    0 0 1 0 0 0 1];

%inc = 0.05; switchtime = 5; maxtime = 200;
%inc = 0.05; switchtime = 10; maxtime = 200;
%inc = 0.1; switchtime = 1; maxtime = 300;
    % inc is the amount of time for the average player step.
    % So crossing time = 20 steps * 0.1 s/step = 2 s
    % maxtime = 300 s.  Based on the original parameter exploration
    % which suggested that APEM should making a switching decision
    % after every unit of travel, we're setting switchtime to 1 here
    % RESULT:  fails to find sufficiently fraught, successful runs 
    % in reasonable time

%inc = 0.1; switchtime = 50; maxtime = 300; figbase = 430; % Trying different values
%inc = 0.1; switchtime = 20; maxtime = 300; figbase = 420;
inc = 0.1; switchtime = 10; maxtime = 300; figbase = 410;
%inc = 0.1; switchtime =  5; maxtime = 300; figbase = 407; 
%inc = 0.1; switchtime = 2.5; maxtime = 300; figbase = 404;


fthresh = 0.1;
%fthresh = 0.2;
%fthresh = 0.3;
tol = 0.05;
limvals = [-1, 1];
mtype = 'C-LGRAD-RJ';
currdir = 'Hierarchy';

% Display params
%saving = 0; printing = 0;
saving = 1; printing = 1;
%figbase = 420; % See above

% Analysis params
leftbranch = [2, 4, 5];
rightbranch = [3, 6, 7];
apex = 1;
switchedorder = [1, 3, 2, 6, 7, 4, 5];
leftcolor = [0.5, 0.5, 0];
rightcolor = [0.5, 0, 0.5];
apexcolor = [0.8, 0, 0];
colorlist = [apexcolor; leftcolor; rightcolor];
%nwarp = 5000;
nwarp = 1000;
%nsought = 20;
%nsought = 200;
nsought = 1000;

% Display params


% Initialize
maxp = 1.5*maxtime/inc; % for mem alloc

nplayers = size(adjacency, 1);
darray = zeros(nwarp+1, nplayers, nsought);
aarray = zeros(nwarp+1, nplayers, nsought);
frecord = zeros(nwarp+1, nsought);

successcount = 0;
dcount = 0;
trycount = 0;
maxtries = 10000;
fprintf('\nOf %d:  ', nsought);
while (dcount < nsought) && (trycount < maxtries)
    % Initialize
    trycount = trycount + 1;

    shist = zeros(nplayers, maxp);
    thist = zeros(1, maxp);
    shist(:, 1) = zeros(nplayers, 1); % Initialize to 0 state
    success = 0;
    dirvec = 2*rpick(2, [nplayers, 1]) - 3; % Randomly pick a direction (1 or -1) for every player
    posvec = inc*dirvec; % Go one meanstep in the selected directions
    currtime = 0;
    nextind = rpick(nplayers);
    nexttime = -1; % Causes model to alter a direction on first update

    scount = 2;
    shist(:, scount) = posvec;

    % Run
    while ~success && (currtime < maxtime)
        scount = scount + 1;

        [posvec, dirvec, currtime, nexttime, nextind] = updatestatelocalgradrs(posvec, adjacency, dirvec, inc, currtime, nexttime, switchtime, nextind, limvals);

        thist(scount) = currtime;

        shist(:, scount) = posvec;

        success = nplayers - abs(sum(posvec)) < tol;
    end

    nsteps = scount;

    % Consolidate
    shist = shist(:, 1:nsteps);
    thist = thist(1:nsteps);
    totaltime = currtime;


    % Check winnning and amount of data

    if success && (scount >= 2)

        successcount = successcount + 1;

        % Time warp
        pwl = PWLF;
        inc = (totaltime-thist(1))/nwarp;
        pwl = pwl.buildf(shist', thist');
        newthist = thist(1):inc:totaltime;
        newthist = newthist';
        newshist = pwl.eval(newthist);

        % Compute fraughtness rate (this is an approximation because it applies
        % the fraghtness of each time window of the rasterized data to the whole
        % window, but with reasonably large nwarp, the value shouldn't be too far
        % off)
        nrstates = size(newshist, 1);
        fhist = zeros(1, nrstates);
        for scount = 1:nrstates
            fhist(scount) = fraught(newshist(scount, :)', adjacency);
        end
        pfraught = sum(fhist)/nrstates;

        % Selecåt relatively fraught cases

        if pfraught > fthresh
            dcount = report(dcount, 10);

            % Flip so that winning side is 1
            if sum(newshist(end, :)) < 0
                newshist = -newshist;
            end

            % If need be, interchange the branches
            if mean(mean(newshist(:, rightbranch))) <...  % Reverse the branches if need be
                    mean(mean(newshist(:, leftbranch)))   % so Left Branch is always the one that reverses
                newshist = newshist(:, switchedorder);  % Switch player assignment of position raster
                %curragrnewshist = curragrnewshist(:, switchedorder);  % Correspondingly, switch player assignment of agreement raster
            end

            % Get agreements
            agreehist = zeros(size(newshist));
            for scount = 1:size(newshist, 1)
                agreehist(scount, :) = getagreelaopt(newshist(scount, :)', adjacency);
            end

            frecord(:, dcount) = fhist';
            darray(:, :, dcount) = newshist;
            aarray(:, :, dcount) = agreehist;
            % Track subgroups

            trycount = 0;
        end
    end
end
toc

successcount

frecord = frecord(:, 1:dcount);
darray = darray(:, :, 1:dcount);
aarray = aarray(:, :, 1:dcount);

ncompleted = dcount;

dmean = mean(darray, 3);
amean = mean(aarray, 3);

%%% Plot by group
tic

addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/SummaryStats';

grouping = {1, [2, 4, 5], [3, 6, 7]};
gnames = {'Apex', 'Left Branch', 'Right Branch'};
gmeans = groupmeans(dmean, grouping);
ameans = groupmeans(amean, grouping);
% Reverse left branch
ameans(:, 2) = -ameans(:, 2);

mvals = mean(gmeans);

npoints = size(gmeans, 1);

% Mean trajectories
figure(figbase + 1);
setfig;
hold on;
xvals = (0:(npoints-1))/(npoints-1);
hh = zeros(3, 1);
for ccount = 1:3
    hh(ccount) = plot(xvals, gmeans(:, ccount), 'Color', colorlist(ccount, :));
end
for ccount = 1:3
    plot(xvals, ameans(:, ccount), 'Color', colorlist(ccount, :), 'LineWidth', 1);
end
xlabel('Normalized Time');
ylabel('Position');
tstring = [' Nruns ', num2str(ncompleted), ' Switchtime ', num2str(switchtime)];
title(mtype, tstring);
%legend(hh, gnames);

%%% Saving and printing
if saving % average data (before group means)
    pfilename = [currdir, '/', topname, '-Model-', mtype, '-nruns-', num2str(nruns), '-st-', num2str(switchtime), '-adata.mat'];
    save(pfilename, 'dmean', 'amean', 'frecord');
    fullmfilename = matlab.desktop.editor.getActiveFilename;
    namefields = splitstring(fullmfilename, '/');
    localmfilename = namefields{end};
    currdate = string(datetime('today', 'Format', 'yyyy-MM-dd'));
    destfile = mergecells([currdir, '/', topname, '-Model-', mtype, '-', currdate, '-', localmfilename]);
    copyfile(localmfilename, destfile);
end

if printing % averaging analysis
    figname = [currdir, '/', topname, '-Model-', mtype, '-nruns-', num2str(nruns), '-st-', num2str(switchtime), '-afig.jpg'];
    print('-djpeg', figname);
end

%sound(gongsound.y, gongsound.Fs);
sound(trainsound.y, trainsound.Fs);

%% Check fraughtness distribution across time

tic
findvals = find(findex);
figure(figbase + 3);
setfig;
histogram(findvals, 20);
xlabel('Trial');
ylabel('Fraughtness count');
title('Fraughtness across trials: ', tstring)
toc


toc

%sound(gongsound.y, gongsound.Fs);



%% Single run

% Simulation parameters
maxp = 3000;  % Estimated number of observation times (for memory allocation)

% Model parameters

%%% Case 1:  Hierarchy
adjacency = [1 1 1 0 0 0 0;
    1 1 0 1 1 0 0;
    1 0 1 0 0 1 1;
    0 1 0 1 0 0 0;
    0 1 0 0 1 0 0;
    0 0 1 0 0 1 0;
    0 0 1 0 0 0 1];
topname = 'Hierarchy-7';


%%% Case 2  ER 10 nodes Data 03 (made with eprob = 0.1---see below)
% ndim = 10;
% eprob = 0.1;
% currdir = 'Erdös-Renyi/Data03';
% %iteration = '01';
% %iteration = '08';
% iteration = '09'; % Classic case with first max near 0.5
% %iteration = '13'; % Mysterious case--turned up a small amount of
%                   % obdurateness at above obd threshold
% %iteration = '14';
% %iteration = '15'; % Notable:  ft correct at 3/5
% %iteration = '16'; % ot correct at 1/3; ft correct at 1/2;
% %iteration = '18'; % Notable case:  ot correct at 1/2; ft correct at 1/3
% %iteration = '20'; % Notable case:  ot correct at 0
% topname = ['Random-ER-', iteration, '-ndim-', num2str(ndim), '-eprob-', num2str(eprob)];
% tstring = [tstring, ' ', topname];
% afilename = [currdir, '/', topname, '-adj.mat'];
% %adjacency = connectedergraph(ndim, eprob, maxsteps);  %% Toggle build
% %save(afilename, 'adjacency'); %% Toggle build
% asource = load(afilename); %% Toggle load
% adjacency = asource.adjacency; %% Toggle load

maxtime = 200;

%inc = 0.05; switchtime = 5;
%inc = 0.05; switchtime = 10;
%inc = 0.1; switchtime = 10; % Fairly substantial fraughtness
inc = 0.1; switchtime = 20;
%inc = 0.1; switchtime = 2; % Low fraughtness
tol = 0.05;

% Display parameters
mtype = 'C-LGRAD-RJ';
figno = 14;

% Initialize

nplayers = size(adjacency, 1);
shist = zeros(nplayers, maxp); % Player positions
lhist = zeros(1, maxp); % Latencies to next change
thist = zeros(1, maxp); % Times of change

dhist = shist; % Direction vectors
%nhist = fhist;  % Player who changes on current step

dirvec = 2*rpick(2, [nplayers, 1]) - 3; % Randomly pick a direction (1 or -1) for every player
limvals = [-1, 1]; % Store the extremal values
posvec = inc*dirvec; % Go one step in the selected directions

scount = 2;
shist(:, scount) = posvec;

success = 0;
scount = 0;
currtime = 0;
nextind = rpick(nplayers);
nexttime = -1; % Causes model to alter a direction on first update

% Run
while ~success && (currtime < maxtime)
    scount = scount + 1;

    [posvec, dirvec, currtime, nexttime, nextind] = updatestatelocalgradrs(posvec, adjacency, dirvec, inc, currtime, nexttime, switchtime, nextind, limvals);

    thist(scount) = currtime;

    shist(:, scount) = posvec;

    success = nplayers - abs(sum(posvec)) < tol;
end

% Consolidate
shist = shist(:, 1:scount);
lhist = lhist(1:scount);
thist = thist(1:scount);
%nhist = nhist(1:scount);

% Compute time history
totaltime = sum(lhist);

% Outcome label
if success
    stroutcome = 'Success';
else
    stroutcome = 'Failure';
end

newshist = shist';
newthist = thist;

% % Time warp
% pwl = PWLF;
% inc = (totaltime-thist(1))/nwarp;
% pwl = pwl.buildf(shist', thist');
% newthist = thist(1):inc:totaltime;
% newthist = newthist';
% newshist = pwl.eval(newthist);

% Compute fraughtness rate (this is an approximation because it applies
% the fraghtness of each time window of the rasterized data to the whole
% window, but with reasonably large nwarp, the value shouldn't be too far
% off)
nrstates = size(newshist, 1);
fhist = zeros(1, nrstates);
for scount = 1:nrstates
    fhist(scount) = fraught(newshist(scount, :)', adjacency);
end
fpercent = sum(fhist)/nrstates;

% Plot
tic
figure(10);
setfig;
fillsegments(fhist, [-1, 1], [0.7, 0.7, 0.7], newthist); % bitlist, ylims, fillcolor
hh = plot(newthist, newshist);
%hh = plot(thist', shist' + 0.01*randn(size(shist')));
tstring = [mtype, ' single run: ', stroutcome, ' Fraughtness: ', num2str(fpercent)];
title(tstring);
legend(hh, num2str((1:nplayers)'));
toc
