% run22 APEM
%
% APEM as in run 19 but with continuous time exponentially sampled (as
% recommended by Iddo)
%
% Used this script to generate the APEM model figure for PLOSone 2023
% submission
%
% The beginning runs multiple runs, seeking a target number of highly
% fraught trials, then does the Averaging Analysis
% Later, there is a subsection that runs a single trial, plotting
% position vs time
% At the end, there is a section that includes the Adverse in plots

tic

addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/PiecewiseLinear';
addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/CellArrays';
addpath '/Users/wht02001/Dropbox/Projects/Solab Dropbox/Tools/Solab Matlab Library/Fraughtness';
addpath '/Users/wht02001/Dropbox/Projects/Solab Dropbox/Tools/Solab Matlab Library/SliderGame';
addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/SummaryStats';

trainsound = load('train');
chirpsound = load('chirp');
handelsound = load('handel');
gongsound = load('gong');

toc


%% Set parameters

tic

% Model params
adjacency = [1 1 1 0 0 0 0;
    1 1 0 1 1 0 0;
    1 0 1 0 0 1 1;
    0 1 0 1 0 0 0;
    0 1 0 0 1 0 0;
    0 0 1 0 0 1 0;
    0 0 1 0 0 0 1];
topname = 'Hierarchy-7';
currdir = 'Hierarchy';
limvals = [-1, 1];
agreethresh = 0.5;
posthresh = 0.5;
mtype = 'C-APEM';
%maxtime = 200; meanstep = 0.2; % Ballparked by trying runs
%maxtime = 1480; meanstep = 1/14; figbase = 500; % 
%maxtime = 300; meanstep = 1/14; figbase = 510;  % Discussion below
%maxtime = 300; meanstep = 2; figbase = 530;
%maxtime = 300; meanstep = 1/5; figbase = 540;
maxtime = 300; meanstep = 0.1; figbase = 500; 
   % Motivated by timing alignment to humans and equation of meanstep with
   % stepsize for arrow keys in lab implementation
%maxtime = 300; meanstep = 0.05; figbase = 510;
%maxtime = 300; meanstep = 0.5; figbase = 520;
%maxtime = 300; meanstep = 1; figbase = 530;
%maxtime = 300; meanstep = 1.5; figbase = 540;
%maxtime = 300; meanstep = 2; figbase = 550;
    % Motivated by comparison to random walk case


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
%nsought = 10;
%nsought = 100;
nsought = 1000;
%nsought = 4000;
maxtries = 10000;

% Display parameters
reportinterval = 50;
fthresh = 0.1;
%figbase = 150; % See Above
%saving = 1; printing = 1; printtitle = 0;
saving = 0; printing = 0; printtitle = 1;

% Note: In the human experiments, the speed of travel is 1 slider
% unit/second.  Thus, it takes 2 seconds for a player to make a
% full traverse of the slider if they go straight across.  The players 
% have 5 minutes to play (in the two lab data sets [CHECK]).  Therefore 
% they can do at most 300/2 = 150 traversals per game. 
% In C-APEM, we can think of the time units as seconds. The players
% move at rate of meanstep units per second so if we set meanstep = 1, then
% we should set maxtime = 300 to allow 60 traverses per game.  In C-APEM,
% meanstep also controls how often the players consider switching
% directions. In earlier tests (of the discrete time version of the model,
% we found tau_max = 10 and inc (the number of units moved per time step) =
% 0.1 resulted in optimal solution time and success rate, which we took
% to be what the humans were doing.
% See 
%
%    /Users/wht02001/Dropbox/Projects/Oscillation/Models/APEM/Model01hierarchy7/BatchParamExplore02/batchparamexplore02.m
%
% In other words, on average, players considered direction switches
% every tau_max/2 = 5 steps, so they each considered switching directions after
% traveling 0.5 units.  In C-APEM, posvec(i) = crop(posvec(i) + latency*dirvec(i))
% where latency is the amount of time taken for the current step and
% dirvec(i) is the rate of movement.  |dirvec(i)| is constantly 1, so the
% speed is always 1 unit/sec, the same value as for the human players.
% Thus, to match the human game, maxtime should be set to 300 seconds.
% The variable meanstep specifies the average latency.  This is thus a free
% parameter which can be adjusted to model the hypothesized average human
% attention rate.  Taking the simulation above as a guide, we note that in
% C-APEM, each player has a 1-in-7 chance of being selected to attend on
% each latency.  Therefore, to match the model result, we set 7*meanstep =
% 0.5, or meanstep = 1/14. We run this case here, but may adjust it to
% an ideal value later based on a batch parameter exploration with C-APEM 
% (see ABatch2023)
% 

% Initialize
nplayers = size(adjacency, 1);
darray = zeros(nwarp+1, nplayers, nsought);
aarray = zeros(nwarp+1, nplayers, nsought);
frecord = zeros(nwarp+1, nsought);

successcount = 0;
dcount = 0;
trycount = 0;

tic
fprintf('\nOf %d:  ', nsought);
while (dcount < nsought) && (trycount < maxtries)
    % Initialize
    trycount = trycount + 1;

    shist = zeros(nplayers, maxp);
    lhist = zeros(1, maxp);
    thist = zeros(1, maxp);
    shist(:, 1) = zeros(nplayers, 1); % Initialize to 0 state
    lhist(1) = 0;
    success = 0;
    dirvec = 2*rpick(2, [nplayers, 1]) - 3; % Randomly pick a direction (1 or -1) for every player
    posvec = meanstep*dirvec; % Go one meanstep in the selected directions
    currtime = meanstep;

    scount = 2;
    shist(:, scount) = posvec;
    lhist(scount) = meanstep;


    % Run
    while ~success && (currtime < maxtime)
        scount = scount + 1;

        [posvec, dirvec, nextind, nextlat] = updatestatecapem(posvec, dirvec, adjacency, meanstep, agreethresh, posthresh, limvals);

        lhist(scount) = nextlat;
        currtime = currtime + nextlat;
        thist(scount) = currtime;

        shist(:, scount) = posvec;

        success = abs(sum(posvec)) == nplayers;
    end
    nsteps = scount;

    % Consolidate
    shist = shist(:, 1:nsteps);
    lhist = lhist(1:nsteps); % Sequence of latencies
    thist = thist(1:nsteps);
    totaltime = currtime;


    %% Check winning

    % Initialize agreehist

    if success

        successcount = successcount + 1;

        % Time warp
        pwl = PWLF;
        inc = (totaltime-thist(1))/nwarp;
        pwl = pwl.buildf(shist', thist');
        newthist = thist(1):inc:totaltime;
        newthist = newthist';
        newshist = pwl.eval(newthist);

        % Get fraughtness rate
        nrstates = size(newshist, 1);
        fhist = zeros(1, nrstates);
        for scount = 1:nrstates
            fhist(scount) = fraught(newshist(scount, :)', adjacency);
        end
        pfraught = sum(fhist)/nrstates;

        % Selecåt relatively fraught cases

        if pfraught > fthresh
            dcount = report(dcount, reportinterval);

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

            % Store data
            frecord(:, dcount) = fhist';
            darray(:, :, dcount) = newshist;
            aarray(:, :, dcount) = agreehist;

            trycount = 0;
        end
    end
end
toc

successcount

frecord = frecord(:, 1:dcount);
darray = darray(:, :, 1:dcount);
aarray = aarray(:, :, 1:dcount);

nruns = dcount;

dmean = mean(darray, 3);
amean = mean(aarray, 3);

%% Plot by group
tic

figure(figbase + 1);
setfig;
hold on;

grouping = {1, [2, 4, 5], [3, 6, 7]};
gnames = {'Apex', 'Left Branch', 'Right Branch'};
gmeans = groupmeans(dmean, grouping);
ameans = groupmeans(amean, grouping);
% Reverse left branch
ameans(:, 2) = -ameans(:, 2);
tvals = (0:nwarp)/(nwarp+1);

%%% Find change pts

leftbranchmean = gmeans(:, 2);
ipt = findchangepts(leftbranchmean', 'MaxNumChanges', 3, 'Statistic', 'linear');

%%% Plot adverse
fcolor = [0.8, 0.8, 0.8];
leftind = tvals(ipt(2));
rightind = tvals(ipt(3));
xvals = [leftind, rightind, rightind, leftind];
%yvals = [limvals(1), limvals(1), limvals(2), limvals(2)];
yvals = [-1.5, -1.5, 1.5, 1.5]; % Overkill
fill(xvals, yvals, fcolor);

% Mean trajectories
hh = zeros(3, 1);
for ccount = 1:3
    hh(ccount) = plot(tvals, gmeans(:, ccount), 'Color', colorlist(ccount, :));
end
for ccount = 1:3
    plot(tvals, ameans(:, ccount), 'Color', colorlist(ccount, :), 'LineWidth', 1);
end

xlabel('Normalized Time');
ylabel('Position');
tstring = [mtype, ' Nruns ', num2str(nruns), ' fthresh ', num2str(fthresh), ' meanstep ', num2str(round(meanstep, 4))];
if printtitle
    title(tstring);
end
%legend(hh, gnames);
ylim(limvals);

%%% Saving and printing
if saving % average data (before group means)
    pfilename = [currdir, '/', topname, '-Model-', mtype, '-nruns-', num2str(nruns), '-fthresh-', num2str(fthresh), '-meanstep-', num2str(round(meanstep, 4)), '-adata.mat'];
    save(pfilename, 'dmean', 'amean', 'frecord');
    fullmfilename = matlab.desktop.editor.getActiveFilename;
    namefields = splitstring(fullmfilename, '/');
    localmfilename = namefields{end};
    currdate = string(datetime('today', 'Format', 'yyyy-MM-dd'));
    destfile = mergecells([currdir, '/', topname, '-Model-', mtype, '-', currdate, '-', localmfilename]);
    copyfile(localmfilename, destfile);
end

if printing % averaging analysis
    figname = [currdir, '/', topname, '-Model-', mtype, '-nruns-', num2str(nruns), '-fthresh-', num2str(fthresh), '-meanstep-', num2str(round(meanstep, 4)), '-afig.jpg'];
    print('-djpeg', figname);
end

toc

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

agreethresh = 0.44; posthresh = 0.5;
%agreethresh = 0.5; posthresh = 0.5;
%agreethresh = 0.2; posthresh = 0.5;
maxtime = 200;
meanstep = 0.1;
%meanstep = 1.5;

% Display parameters
mtype = 'C-APEM';
figno = 14;

% Initialize

nplayers = size(adjacency, 1);
shist = zeros(nplayers, maxp); % Player positions
lhist = zeros(1, maxp); % Latencies to next change
thist = zeros(1, maxp); % Times of change

dhist = shist; % Direction vectors
fhist = zeros(1, maxp); % Fraughtness at current step
%nhist = fhist;  % Player who changes on current step

dirvec = 2*rpick(2, [nplayers, 1]) - 3; % Randomly pick a direction (1 or -1) for every player
limvals = [-1, 1]; % Store the extremal values
posvec = meanstep*dirvec; % Go one step in the selected directions

scount = 2;
shist(:, scount) = posvec;

success = 0;
scount = 0;
currtime = 0;

% Run
while ~success && (currtime < maxtime)
    scount = scount + 1;

    [posvec, dirvec, nextind, nextlat] = updatestatecapem(posvec, dirvec, adjacency, meanstep, agreethresh, posthresh, limvals);

    lhist(scount) = nextlat;
    currtime = currtime + nextlat;
    thist(scount) = currtime;

    %nhist(scount) = nextind;
    fhist(scount) = fraught(posvec, adjacency);
    shist(:, scount) = posvec;
    dhist(:, scount) = dirvec;

    success = abs(sum(posvec)) == nplayers;
end

% Consolidate
shist = shist(:, 1:scount);
lhist = lhist(1:scount);
fhist = fhist(1:scount);
thist = thist(1:scount);
%nhist = nhist(1:scount);

% Compute time history
totaltime = sum(lhist);

% Compute percent fraughtness
fpercent = sum(lhist(fhist == 1))/totaltime;

% Outcome label
if success
    stroutcome = 'Success';
else
    stroutcome = 'Failure';
end

% Plot
figure(figno);
setfig;
jj = plot(thist', shist');
% Obtain colors
cmat = zeros(nplayers, 3);
for pcount = 1:nplayers
    cmat(pcount, :) = jj(pcount).Color;
end
fillsegments(fhist, [-1, 1], [0.9, 0.9, 0.9], thist'); % bitlist, ylims, fillcolor, time history
hh = plot(thist', shist');
%plot(thist', shist' + 0.01*randn(size(shist')));
tstring = [mtype, ' single run: ', stroutcome, ' Fraughtness: ', num2str(fpercent)];
title(tstring);
legend(hh, num2str((1:nplayers)'));

%% Find Adverse

%%%% NOTE: If printing is on, this loads a mat file, produces a figure from it,
%%%% adjusts the figure, and then saves to the corresponding figname, overwriting 
%%%% whatever figure was previously saved with that name

tic
% Display params
figbase = 300;
printing = 0; printtitle = 1;
%printing = 1; printtitle = 0;

% Model parameters
limvals = [-1, 1];

% Analysis params
leftbranch = [2, 4, 5];
rightbranch = [3, 6, 7];
apex = 1;
switchedorder = [1, 3, 2, 6, 7, 4, 5];
leftcolor = [0.5, 0.5, 0];
rightcolor = [0.5, 0, 0.5];
apexcolor = [0.8, 0, 0];
colorlist = [apexcolor; leftcolor; rightcolor];

%%% Load data
currdir = 'Hierarchy';
localfilename = 'Hierarchy-7-Model-C-APEM-nruns-1000-fthresh-0.1-meanstep-0.1-a';
figdata = [currdir, '/', localfilename, 'data.mat'];
tstring = 'C-APEM';

playermeans = load(figdata);

dmean = playermeans.dmean;
amean = playermeans.amean;
nwarp = size(dmean, 1) - 1;

%%% Plot by group
tic

addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/SummaryStats';

grouping = {1, [2, 4, 5], [3, 6, 7]};
gnames = {'Apex', 'Left Branch', 'Right Branch'};
gmeans = groupmeans(dmean, grouping);
ameans = groupmeans(amean, grouping);
% Reverse left branch
ameans(:, 2) = -ameans(:, 2);
tvals = (0:nwarp)/(nwarp+1);

mvals = mean(gmeans);

%%% Find change pts

leftbranchmean = gmeans(:, 2);
ipt = findchangepts(leftbranchmean,'MaxNumChanges',3,'Statistic','mean');
%ipt = findchangepts(leftbranchmean,'MaxNumChanges',3,'Statistic','rms');

%%% Plot change pts
% for changept = ipt'
%     plot([tvals(changept); tvals(changept)], [-1; 1], 'k-', 'LineWidth', 1);
% end

%%% Plot adverse
figure(figbase + 1);
setfig;
hold on;

fcolor = [0.8, 0.8, 0.8];
leftind = tvals(ipt(2));
rightind = tvals(ipt(3));
xvals = [leftind, rightind, rightind, leftind];
yvals = [limvals(1), limvals(1), limvals(2), limvals(2)];
fill(xvals, yvals, fcolor);

%%% Plot Mean trajectories

hh = zeros(3, 1);
for ccount = 1:3
    hh(ccount) = plot(tvals, gmeans(:, ccount), 'Color', colorlist(ccount, :));
end
for ccount = 1:3
    plot(tvals, ameans(:, ccount), 'Color', colorlist(ccount, :), 'LineWidth', 1);
end
xlabel('Normalized Time');
ylabel('Position');
%tstring = [mtype, ' Nruns ', num2str(nruns), ' fthresh ', num2str(fthresh), ' meanstep ', num2str(round(meanstep, 4))];
if printtitle
    title(tstring);
end
%legend(hh, gnames);

if printing % averaging analysis
    figname = [currdir, '/', localfilename, 'fig.jpg'];
    print('-djpeg', figname);
end

toc

