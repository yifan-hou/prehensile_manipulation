% Draw multiple cones in wrench space
function drawWrenchSpace(cone_all_fix, eCone_allFix, hCone_allFix, V_control_directions, F_control_directions, ...
    cone_generators, eh_modes, number_of_modes, ...
    compatibilities, goal_id)

% blue shades
blue=[58/255,67/255,186/255];
sky=[98/255,197/255,218/255];
% navy=[11/255,17/255,113/255];
% indigo=[40/255,30/255,93/255];
teal=[72/255,170/255,173/255];
ocean=[1/255,96/255,100/255];
peacock=[1/255,45/255,54/255];
% azure=[22/255,32/255,166/255];
cerulean=[4/255,146/255,194/255];
lapis=[39/255,50/255,194/255];
spruce=[44/255,62/255,76/255];
aegean=[30/255,69/255,110/255];
% blueberry=[36/255,21/255,112/255];
% denim=[21/255,30/255,61/255];
% admiral=[6/255,16/255,148/255];
sapphire=[82/255,178/255,192/255];
artic=[130/255,237/255,253/255];
bluecolors = [blue;sky;teal;ocean;cerulean;lapis;aegean;sapphire;artic];
bluecolors = [bluecolors; bluecolors];
nbcolor = size(bluecolors, 1);

% red shades
red=[208/255,49/255,45/255];
cherry=[153/255,15/255,2/255];
rosered=[226/255,37/255,43/255];
jam=[96/255,15/255,11/255];
merlot=[84/255,31/255,27/255];
garnet=[96/255,11/255,4/255];
crimson=[184/255,15/255,10/255];
ruby=[144/255,6/255,3/255];
scarlet=[145/255,13/255,9/255];
winered=[76/255,8/255,5/255];
brick=[126/255,40/255,17/255];
apple=[169/255,27/255,13/255];
mahogany=[66/255,13/255,9/255];
blood=[113/255,12/255,4/255];
sangriared=[94/255,25/255,20/255];
berryred=[121/255,24/255,18/255];
currant=[103/255,12/255,7/255];
blushred=[188/255,84/255,73/255];
candy=[210/255,21/255,2/255];
lipstick=[156/255,16/255,3/255];
redcolors = [red;cherry;rosered;jam;merlot;garnet;crimson;ruby;scarlet;winered;brick;apple;mahogany;blood;sangriared;berryred;currant;blushred;candy;lipstick];
nrcolor = size(redcolors, 1);

rorder_b = randperm(nbcolor);
rorder_r = randperm(nrcolor);
bcount = 1;
rcount = 1;

% exam all the modes
% draw contact screws
figure(1); clf(1); hold on;

% draw the force controlled subspace
% slate=[63/255,61/255,83/255];
% drawCone(10*F_control_directions, ocean, false);
% drawCone(10*V_control_directions, blueberry, true);

% draw the origin
plot3(0,0,0,'k.','markersize',25);

% draw the whole hand cone and e cone
% drawCone(cone_all_fix, sapphire, true);
% drawCone([eCone_allFix; -hCone_allFix]', peacock, true);
drawCone(eCone_allFix', sapphire, true);
drawCone(hCone_allFix', blue, true);

axis([-0.8 0.8 -0.1 1.1 -0.4 0.25]);
set(gca, 'XTick', -0.4:0.2:0.4, 'YTick', -0.1:0.25:1.1, 'ZTick', -0.4:0.325:.25);
view(-135, 27);


% draw cones
for i = 1:number_of_modes
    if i == goal_id
        fprintf('Mode %d: (Goal)\n', i);
    else
        fprintf('Mode %d:\n', i);
    end
    texts = printModes(eh_modes(:, i), false);
    if compatibilities(i) == true
        disp(texts);
        bcolor = bluecolors(rorder_b(bcount), :);
        bcount = bcount + 1;
        drawCone(cone_generators{i}, bcolor, true, texts);

%         % projection on force plane
%         cone_projection_i = force_basis*(force_basis')*cone_generators{i};
%         drawCone(cone_projection_i, bcolor, true);
    else
        disp([texts ' Infeasible']);
%         rcolor = redcolors(randi(nrcolor), :);
        bcolor = bluecolors(rorder_b(bcount), :);
        bcount = bcount + 1;
        drawCone(cone_generators{i}, bcolor, true, texts);
    end
end


% axis([-0.8 0.8 -0.1 1.1 -0.5 0.5]);
% set(gca, 'XTick', -0.8:0.4:0.8, 'YTick', -0.1:0.25:1.1, 'ZTick', -0.5:0.25:.5);
axis([-0.4 0.4 -0.1 1.1 -0.4 0.25]);
set(gca, 'XTick', -0.4:0.2:0.4, 'YTick', -0.1:0.25:1.1, 'ZTick', -0.4:0.325:.25);
view(-135, 27);
