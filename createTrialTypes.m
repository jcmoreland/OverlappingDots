 function cond = createTrialTypes(cond)
% Create trial set for experiment that is fully counterbalanced

cond.blockType = cond.blockvector;  % A vector of the block types in this run
cond.blockvector = reshape(repmat(cond.blockvector,cond.nTrials,1),1,cond.nBlocks*cond.nTrials);    % Make a vector which says the block type of every trial in the experiment

%--- Trial types:
dualsameup = find(cond.blockvector==1);
dualsamedown = find(cond.blockvector==2);
dualdifup = find(cond.blockvector==3);
dualdifdown = find(cond.blockvector==4);

singleleftsameup = find(cond.blockvector==5);
singleleftsamedown = find(cond.blockvector==6);
singleleftdifup = find(cond.blockvector==7);
singleleftdifdown = find(cond.blockvector==8);

singlerightsameup = find(cond.blockvector==9);
singlerightsamedown = find(cond.blockvector==10);
singlerightdifup = find(cond.blockvector==11);
singlerightdifdown = find(cond.blockvector==12);

%--- Block Count
cond.blocknum = cond.blockvector;

%--- Main Types - more useful for analysis than anything else.
cond.dual = zeros(1,length(cond.blockvector));
cond.single = zeros(1,length(cond.blockvector));
cond.dualsame = zeros(1,length(cond.blockvector));
cond.dualdif = zeros(1,length(cond.blockvector));
cond.singlesame = zeros(1,length(cond.blockvector));
cond.singledif = zeros(1,length(cond.blockvector));

cond.dual([dualsameup,dualsamedown,dualdifup,dualdifdown]) = 1;
cond.single([singleleftsameup,singleleftsamedown,singleleftdifup,singleleftdifdown,...
    singlerightsameup,singlerightsamedown,singlerightdifup,singlerightdifdown]) = 1;
cond.dualsame([dualsameup,dualsamedown]) = 1;
cond.dualdif([dualdifup,dualdifdown]) = 1;
cond.singlesame([singleleftsameup,singleleftsamedown,singlerightsameup,singlerightsamedown]) = 1;
cond.singledif([singleleftdifup,singleleftdifdown,singlerightdifup,singlerightdifdown]) = 1;

%--- Attend direction for each side (info for the arrow)
cond.attendL = zeros(1,length(cond.blockvector));
cond.attendR = zeros(1,length(cond.blockvector));

cond.attendL([dualsameup,dualdifup,singleleftsameup,...
    singleleftdifup,singlerightsameup,singlerightdifup]) = 1;% 0 = dir1, 1 = dir2
cond.attendR([dualsameup,dualdifdown,singleleftsameup,...
    singleleftdifdown,singlerightsameup,singlerightdifdown]) = 1;% 0 = dir1, 1 = dir2
% cond.attendL and .attendR hold which field is being attended
% on each side. 0=dir1, 1=dir2.
% Encode this on 1:4 to match the for loop:
cond.attendL = cond.attendL + 1; % now 0->1,1->2
cond.attendR = cond.attendR + 3; % now 0->3,1->4
            
%--- Response side
cond.responseL = zeros(1,length(cond.blockvector));
% Single side blocks 
cond.responseL([singleleftsameup,singleleftsamedown,singleleftdifup,singleleftdifdown]) = 1; % Left side response when == 1
% Dual task blocks
cond.responseL([dualsameup(1:length(dualsameup)/2),...
    dualsamedown(1:length(dualsamedown)/2),...
    dualdifup(1:length(dualdifup)/2),...
    dualdifdown(1:length(dualdifdown)/2)]) = 1;

%--- Events
% Equal number of left and right events and events can only occur in the
% attended fields.
cbal = [0,1,1,0;0,1,0,1]; % every possible combination of target present/absent on each field
trialEvents = repmat(cbal,1,cond.nTrials*cond.nBlocks/4);
cond.trialTargetR = trialEvents(1,:); % For the attended field on the right this codes whether an event occurs.
cond.trialTargetL = trialEvents(2,:);

%--- Shuffle trials within block types
blockTypes = unique(cond.blockType);
for b = 1:length(blockTypes)
    curTrs = find(cond.blockvector==blockTypes(b));
    nT = length(curTrs);
    
    [~,index] = Shuffle(1:nT);
    tmptrialindx = find(cond.blockvector==blockTypes(b));
    
    cond.dual(tmptrialindx) = cond.dual(curTrs(index));
    cond.single(tmptrialindx) = cond.single(curTrs(index));
    cond.dualsame(tmptrialindx) = cond.dualsame(curTrs(index));
    cond.dualdif(tmptrialindx) = cond.dualdif(curTrs(index));
    cond.singlesame(tmptrialindx) = cond.singlesame(curTrs(index));
    cond.singledif(tmptrialindx) = cond.singledif(curTrs(index));
    
    cond.attendL(tmptrialindx) = cond.attendL(curTrs(index));
    cond.attendR(tmptrialindx) = cond.attendR(curTrs(index));
    
    cond.responseL(tmptrialindx) = cond.responseL(curTrs(index));
    cond.trialTargetL(tmptrialindx) = cond.trialTargetL(curTrs(index));
    cond.trialTargetR(tmptrialindx) = cond.trialTargetR(curTrs(index));

end

%--- Now randomise blocks
[cond.blockType,~] = Shuffle(cond.blockType); % Shuffle the block order
trialOrder = [];
checkType = zeros(1,length(cond.blocklabel));  % a vector to keep track of how many blcks we have done 

for i = 1:length(cond.blockType) % for each block
    blocktype = cond.blockType(i); % block type for this frst block
    % update the block count for this type
    checkType(blocktype) = checkType(blocktype) + 1;
    % find all the trials of this type
    curBlTr = find(cond.blockvector == blocktype);
    % of these take the indices of the next 1:cond.nTrials that have this block type
    curTrials = curBlTr((checkType(blocktype)-1)*cond.nTrials+ 1:checkType(blocktype)*cond.nTrials );
    % add these indices to the overall trialOrder
    trialOrder = [trialOrder, curTrials];
end
cond.dual = cond.dual(trialOrder);
cond.single = cond.single(trialOrder);
cond.dualsame = cond.dualsame(trialOrder);
cond.dualdif = cond.dualdif(trialOrder);
cond.singlesame = cond.singlesame(trialOrder);
cond.singledif = cond.singledif(trialOrder);

cond.attendL = cond.attendL(trialOrder);
cond.attendR = cond.attendR(trialOrder);

cond.responseL = cond.responseL(trialOrder);
cond.trialTargetL = cond.trialTargetL(trialOrder);
cond.trialTargetR = cond.trialTargetR(trialOrder);