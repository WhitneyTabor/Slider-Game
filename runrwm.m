% run15 random walk
%
% RESULTS:  Pretty far from fitting the human data with small stepsize
% (because it's hard to achieve highly fraught runs, and even harder to
% achieve highly fraught runs that are successful)
% Using run with stepsize 2, fthresh 0.1 nsought = 1000 for
% 2023 PLOSone submission

%% Libraries

tic

addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/PiecewiseLinear';
addpath '/Users/wht02001/Dropbox/Matlab/Mfiles/SummaryStats';
addpath '/Users/wht02001/Dropbox/Projects/Solab Dropbox/Tools/Solab Matlab Library/Fraughtness';
addpath '/Users/wht02001/Dropbox/Projects/Solab Dropbox/Tools/Solab Matlab Library/SliderGame';

trainsound = load('train');
chirpsound = load('chirp');
handelsound = load('handel');
gongsound = load('gong');

toc

%% Simulation

% Simulation parameters
mtype = 'RW';
currbase = '/Users/wht02001/Dropbox/Projects/Solab Dropbox/Projects/Slider Game/Model05Compromiseexplore/WhitsRuns';
currdir = [currbase, '/Hierarchy'];
if not(isfolder(currdir))
    mkdir(currdir)
end

% Display parameters
figbase = 320;
reportinterval = 20;
%reportinterval = 100;
%saving = 1; printing = 1;
saving = 0; printing = 0;

% Model parameters
% Case: Hierarchy
topname = 'hierarchy-7';
adjacency = [1 1 1 0 0 0 0; % hierarchy
    1 1 0 1 1 0 0;
    1 0 1 0 0 1 1;
    0 1 0 1 0 0 0;
    0 1 0 0 1 0 0;
    0 0 1 0 0 1 0;
    0 0 1 0 0 0 1];
%stepsize = 0.25;
%stepsize = 0.5;
stepsize = 2; % Constant flipping
%maxp = 10000;
maxp = 2000;

% Analysis parameters
%fthresh = 0;
%fthresh = 0.001;
%fthresh = 0.02;
%fthresh = 0.01;
%fthresh = 0.1;
fthresh = 0.2;
%fthresh = 0.3;
%fthresh = 0.4;
%fthresh = 1/(2^5) + 0.01; % Somewhat arbitrary choice
leftbranch = [2, 4, 5];
rightbranch = [3, 6, 7];
apex = 1;
switchedorder = [1, 3, 2, 6, 7, 4, 5];
leftcolor = [0.5, 0.5, 0];
rightcolor = [0.5, 0, 0.5];
apexcolor = [0.8, 0, 0];
colorlist = [apexcolor; leftcolor; rightcolor];
nwarp = 1000;
%nsought = 20;
%nsought = 100;
nsought = 200;
%nsought = 500;
%nsought = 1000;
%nsought = 2000;
%nsought = 4000;

% Initialize
nplayers = size(adjacency, 1);
tstring = [topname, ' ', mtype];
darray = zeros(nwarp, nplayers, nsought);
aarray = zeros(nwarp, nplayers, nsought);
dcount = 0;
stally = 0;

% Run
tic
fprintf('\nOf %d:  ', nsought);
while dcount < nsought

    % Initialize
    ndim = size(adjacency, 1);
    currstate = zeros(ndim, 1);
    shist = zeros(ndim, maxp);
    shisst(:, 1) = zeros(ndim, 1); % Initialize to 0 state
    success = 0;
    scount = 1;

    %%% Run

    while ~success && (scount < maxp)
        %%%%%
        scount = scount + 1;
        ds = rpick(3, [ndim, 1]) - 2;
        currstate = crop(currstate + stepsize*ds, -1, 1);
        %currstate = crop(currstate + 2*ds/3, -1, 1); % Never on faces
        shist(:, scount) = currstate;

        success = abs(sum(currstate)) == ndim;
        %%%%%
    end
    nsteps = scount;

    % Consolidate, transpose
    shist = shist(:, 1:nsteps);
    shist = shist';


    % Check winning
    if success

        fvals = zeros(1, nsteps);
        for scount = 1:nsteps
            fvals(scount) = fraught(shist(scount, :)', adjacency);
        end
        pfraught = sum(fvals == 1)/nsteps;

        %         stally = stally + 1;
        %         fprintf('(%5.3f)', pfraught)
        %         if mod(stally, 10) == 0
        %             fprintf('\n');
        %         end

        % Select relatively fraught cases

        if pfraught >= fthresh
            dcount = report(dcount, reportinterval);

            % Flip so that winning side is 1
            if sum(shist(end, :)) < 0
                shist = -shist;
            end

            % Time warp
            pwl = PWLF;
            rawlength = size(shist, 1);
            inc = rawlength/nwarp;
            pwl = pwl.buildf(shist);
            newshist = pwl.eval(inc:inc:rawlength);

            % Flip so that the Right Branch is 1-biased
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

            darray(:, :, dcount) = newshist;
            aarray(:, :, dcount) = agreehist;

            % Track subgroups

        end
    end
end

nruns = dcount;


dmean = mean(darray, 3);
amean = mean(aarray, 3);

%%% Check fraughtness for aligned data

tic
fcount = 0;
findex = zeros(nwarp, 1);
for scount = 1:nwarp
    if fraught(dmean(scount, :)', adjacency)
        fcount = fcount + 1;
        findex(scount) = 1;
    end
end
mfraughtness = fcount/nwarp;
findex = logical(findex);
toc

%%% Check fraughtness distribution across time

% tic
% findvals = find(findex);
% figure(figbase + 3);
% setfig;
% histogram(findvals);
% toc

%%% Plot by group
tic

grouping = {1, [2, 4, 5], [3, 6, 7]};
gnames = {'Apex', 'Left Branch', 'Right Branch'};
gmeans = groupmeans(dmean, grouping);
ameans = groupmeans(amean, grouping);
% Reverse left branch
ameans(:, 2) = -ameans(:, 2);

% Mean trajectories
figure(figbase + 1);
hold off;
setfig;
hold on;
hh = zeros(3, 1);
for ccount = 1:3
    hh(ccount) = plot(gmeans(:, ccount), 'Color', colorlist(ccount, :));
end
for ccount = 1:3
    plot(ameans(:, ccount), 'Color', colorlist(ccount, :), 'LineWidth', 1);
end
tstring = [tstring, ' Nruns ', num2str(nruns), ' fthresh ', num2str(fthresh)];
title(tstring);
legend(hh, gnames);

%%% Saving and printing
if saving % average data (before group means)
    pfilename = [currdir, '/', topname, '-Model-', mtype, '-nruns-', num2str(nruns), '-fthresh-', num2str(fthresh), '-adata.mat'];
    save(pfilename, 'dmean', 'amean', 'frecord');
    fullmfilename = matlab.desktop.editor.getActiveFilename;
    namefields = splitstring(fullmfilename, '/');
    localmfilename = namefields{end};
    currdate = string(datetime('today', 'Format', 'yyyy-MM-dd'));
    destfile = mergecells([currdir, '/', topname, '-Model-', mtype, '-', currdate, '-', localmfilename]);
    copyfile(localmfilename, destfile);
end

if printing % averaging analysis
    figname = [currdir, '/', topname, '-Model-', mtype, '-nruns-', num2str(nruns), '-fthresh-', num2str(fthresh), '-afig.jpg'];
    print('-djpeg', figname);
end

toc

%sound(gongsound.y, gongsound.Fs);

%% Individual trial

%% Single run

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
%stepsize = 0.25;
stepsize = 2; % Constant flipping
maxp = 2000;


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

%tol = 1.5;

% Display parameters
mtype = 'Random Walk';
figbase = 24;
nwarp = 1000;

% Initialize
ndim = size(adjacency, 1);
currstate = zeros(ndim, 1);
shist = zeros(ndim, maxp);
shisst(:, 1) = zeros(ndim, 1); % Initialize to 0 state
success = 0;
scount = 1;

% Run
while ~success && (scount < maxp)
    %%%%%
    scount = scount + 1;
    ds = rpick(3, [ndim, 1]) - 2;
    currstate = crop(currstate + stepsize*ds, -1, 1);
    %currstate = crop(currstate + 2*ds/3, -1, 1); % Never on faces
    shist(:, scount) = currstate;

    success = abs(sum(currstate)) == ndim;
    %%%%%
end

% Consolidate
shist = shist(:, 1:scount);

%%% Plot

% figure(13);
% setfig
% plot(shist');
%
%
% figure(14);
% setfig;
% %plot(shist');
% plot(2.^shist(1, :));
%
% figure(15);
% setfig;
% plot(disscrephist');

%%% Outcome label
if success
    stroutcome = 'Success';
else
    stroutcome = 'Failure';
end

% Time warp
pwl = PWLF;
inc = scount/nwarp;
pwl = pwl.buildf(shist', (1:scount)');
newthist = 1:scount;
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
fpercent = sum(fhist)/nrstates;

% Plot
tic
figure(figbase + 5);
setfig;
fillsegments(fhist, [-1, 1], [0.7, 0.7, 0.7], newthist); % bitlist, ylims, fillcolor
hh = plot(newthist, newshist);
%hh = plot(thist', shist' + 0.01*randn(size(shist')));
tstring = [mtype, ' single run: ', stroutcome, ' Fraughtness: ', num2str(fpercent)];
title(tstring);
legend(hh, num2str((1:nplayers)'));
toc


