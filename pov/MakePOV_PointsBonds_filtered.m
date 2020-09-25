function MakePOV_PointsBonds_filtered(filenameBonds, filenamePar, outputName, ninst)

filenameBonds = 'base';

%for i = 1:nist
    i = 10000;
    fidW1 = fopen(strcat(outputName, '_', (num2str(i)), '.pov'), 'w');
    fidW2 = fopen(strcat('3plet_', (num2str(i)), '.pov'), 'w');
    filenameBonds = strcat(filenameBonds, '.', num2str(i));
    filenamePar = strcat(filenamePar, '.', num2str(i));
    neigh_info = read_demsi(filenameBonds);
    particle_info = read_demsi(filenamePar);
    [uIDS, ind, ~] = unique(sort(neigh_info(:,[1,2]),2), 'rows');
    u_neigh = neigh_info(ind,:);
    triplets = [];
    
    for i = 1:length(particle_info(:,1))
        
        start1 = particle_info(i,1);
        ids1 = [];
        ids1 = [ids1; uIDS(uIDS(:,1)==start1,:); uIDS(uIDS(:,2)==start1,end:-1:1)];
        candidates1 = ids1(:,2);
        target = start1;
        match = [];
        
        for j = 1:length(candidates1)
            start2 = candidates1(j);
            ids2 = [];
            ids2 = [ids2; uIDS(uIDS(:,1)==start2,:); uIDS(uIDS(:,2)==start2,end:-1:1)];  
            candidates3 = ids2(:,2);
            candidates3(candidates3==target) = [];
            
                    for k = 1:length(candidates3)
                        start3 = candidates3(k);
                        ids3 = [];
                        ids3 = [ids3; uIDS(uIDS(:,1)==start3,:); uIDS(uIDS(:,2)==start3,end:-1:1)];
                        if any(ids3(:,2) == target)
                            match = [match; target, candidates1(j), ids3(ids3(:,2) == target,1)];
%                            break;
                        end
                    end
        
        end
        
        if ~isempty(match)
            for w = 1:length(match(:,1))
                fprintf(fidW2, '%i %i %i\n', match(w,1), match(w,2), match(w,3));
            end
        end
        triplets = [triplets; match];
        
    end
    
    

fprintf(fidW, '#include "colors.inc"\n#include "textures.inc"\n#include "functions.inc"\n#include "glass.inc"\n');
fprintf(fidW, '\ncamera {\nlocation <0, 0, 100000>\nlook_at <200000,250000,0>\nangle 0}\n');
fprintf(fidW, '\n#macro point_cloud(A, r)\n#local _pce = dimension_size(A,1);\n#local c = 0;\n#while (c<_pce)\nsphere{A[c],r}\n#local c=c+1;\n#end\n#end\n\n#declare radiusA = 0.2;\n#declare radiusB = 0.2;\n#declare radiusC = 0.25;\n');
fprintf(fidW, '\nunion{\n');

radius = 0.3;

tline = fgetl(fidP);
timesteps = str2num(fgetl(fidP));
display(sprintf('Timestep: %f', timesteps));
tline = fgetl(fidP);
numberOfbeads = str2num(fgetl(fidP));

display(sprintf('Number of beads: %g', numberOfbeads));
tline = fgetl(fidP);
box = str2num(fgetl(fidP));

tline = fgetl(fidP);
tline = fgetl(fidP);
tline = fgetl(fidP);

tline = fgetl(fidB);
tline = fgetl(fidB);
tline = fgetl(fidB);

numberOfbonds = str2num(fgetl(fidB));
display(sprintf('Number of bonds: %g', numberOfbonds));
display(sprintf('Box: %g %g', box));

tline = fgetl(fidB);
tline = fgetl(fidB);
tline = fgetl(fidB);
tline = fgetl(fidB);
tline = fgetl(fidB);

coords = zeros(numberOfbeads,6);
bonds = zeros(numberOfbonds,3);
spacing = 15;
ind = 0;

for i = 1:numberOfbeads
    
%     ind = ind + 1;
%     if ind == spacing
%        fprintf(fidW, '\n');
%        ind = 1;
%     end
    
    local = str2num(fgetl(fidP));
%    local(2:3) = [];    
    coords(i,:) = local;
%     if i == numberOfbeads
%         fprintf(fidW, '<%g,%g,%g>', coords(i,4:end));
%     else
%         fprintf(fidW, '<%g,%g,%g>,', coords(i,4:end));
%     end
    
end

a = 0.75; b = 0.3; c = 0.3;
d = -1; e = 2.21; f = 0.23;
coords1 = coords(coords(:,4).*a + coords(:,5).*b + coords(:,6).*c -3 < 0, :);
coords2 = coords(-2*coords(:,4) - coords(:,5).*1 + coords(:,6).*5.23 + 56 > 0, :);
coordsF = union(coords1,coords2,'rows');
rem_ids = setdiff(coords(:,1), coordsF(:,1));
gel = coordsF(coordsF(:,3) == 2, :);
gel_more = coordsF(coordsF(:,3) == 7, :);
gel = [gel; gel_more];
%coords = coordsF;

for i = 1:numberOfbonds
    bonds(i,:) = str2num(fgetl(fidB));
end

for i = 1:length(rem_ids)

    [a,~] = find(rem_ids(i) == bonds(:, 2:3));
    
    if ~isempty(a)
        bonds(a,:) = [];
    end
    
end

gelBonds = bonds(bonds(:,1) == 2, :);
gelBondsOu = [];

[a1,b1] = ismember(gelBonds(:,2), gel_more(:,1)); % finding gel_more bonds
[a2,b2] = ismember(gelBonds(:,3), gel_more(:,1));
sum_avecs = a1 + a2;
logic_a = sum_avecs == 2;
gelBondsMore = gelBonds(logic_a, :); % gel_more bonds

 xxx = [gelBonds(:,2) gelBonds(:,3)];                   % get the bonds
 AAA = sparse(xxx(:,1), xxx(:,2), 1);                         % make sparse matrix from edge list
 [rows, cols] = size(AAA);
  if rows ~= cols                                             % make sure that matrix AAA is square
     if rows > cols
         column_New = zeros(rows, rows - cols);
         AAA = [AAA column_New];
     elseif cols > rows 
         rows_New = zeros(cols - rows, cols);
         AAA = [AAA; rows_New];
     end
 end
AAA = AAA + AAA.';                                           % make the graph undirected
%  h = biograph(AAA);                                           % make biograph object
%  [~, pred, ~] = traverse(h, 1, 'Directed', false);        % nodes which are not connected are outputted as NaN the column index in pred is the atom id
[~, pred, ~] = graphtraverse(AAA, 126, 'Directed', false, 'Method', 'DFS');
[~, col] = find(isnan(pred));                % find NaN nodes in the array
[~, ~, ib1] = intersect(col, xxx(:,1));           % finds bonds in first column which have to be eliminated
[~, ~, ib2] = intersect(col, xxx(:,2));           % finds bonds in second column which have to be eliminated
i = union(ib1, ib2);               % this basically tells you which rows to eliminate
%gelBonds(i, :) = [];            % eliminates the rows from the bond list
[~, ~, ib3] = intersect(col, gel(:,1));               % find the row indicies for atom list
gel(ib3,:) = [];           % remove beads from bead list
 
%logic_a(i,:) = [];
gelBonds = gelBonds(~logic_a, :);

%gel = coordsF(coordsF(:,3) == 2, :);
shell = coordsF(coordsF(:,3) == 3, :);
chains1 = coordsF(coordsF(:,3) == 4, :);
chains2 = coordsF(coordsF(:,3) == 5, :);
particles = coordsF(coordsF(:,3) == 6, :);
% gel_more = coordsF(coordsF(:,3) == 7, :);
% gel = [gel; gel_more];

aveXshell = mean(shell(:,4)); aveYshell = mean(shell(:,5)); aveZshell = mean(shell(:,6));
%distP = sqrt( (particles(:,4) - aveXshell).^2 + (particles(:,5) - aveYshell).^2 + (particles(:,6) - aveZshell).^2 ) <= 24;

%particles = particles(distP, :);
distPar = sqrt( (particles(:,4) - aveXshell).^2 + (particles(:,5) - aveYshell).^2 + (particles(:,6) - aveZshell).^2 );
particlesNearShell = (distPar < 80);
particlesNearShell = particles(particlesNearShell, :);

moleculeStart = min(particlesNearShell(:,2));
particlesNearGel = 0;
particlesInsideGel = 0;
molIDS = [];
cutoffNearGel = 3.5;
capturedP = 0;
cutoffPG = 2.5;
kurrMOL = 0;

%for ii = 1:(length(unique(particlesNearShell(:,2))) - 1)
for ii = 1:(length(particlesNearShell(:,2)))
    
        kurrMOL = particlesNearShell(ii,2);
        
        if ~ismember(kurrMOL, molIDS)
                
            data = particlesNearShell(particlesNearShell(:,2) == particlesNearShell(ii,2), :);
            dataAveX = mean(data(:,4)); dataAveY = mean(data(:,5)); dataAveZ = mean(data(:,6));

    %         for j = 1:length(data(:,1))

            distParJ = sqrt( (dataAveX - gel_more(:,4)).^2 + (dataAveY - gel_more(:,5)).^2 + (dataAveZ - gel_more(:,6)).^2 );
            indDist = (distParJ < cutoffNearGel);
            minDist = min(distParJ);

            gelNearParticle = gel_more(indDist, :);
            gelAveXnear = mean(gelNearParticle(:,4)); gelAveYnear = mean(gelNearParticle(:,5)); gelAveZnear = mean(gelNearParticle(:,6));
            distCenters = sqrt( (dataAveX - gelAveXnear).^2 + (dataAveY - gelAveYnear).^2 + (dataAveZ - gelAveZnear).^2 );

           
            if distCenters <= cutoffNearGel./2

                particlesNearShell(ii, 3)...
                    = 20;

            elseif distCenters <= cutoffNearGel.*2

%                 particlesNearShell(ii, 3)...
%                     = 30;
                  continue;
                  
            elseif minDist <= cutoffPG.*3.5

                particlesNearShell(ii, 3)...
                    = 30;

            end     
            
        end

end

distP = sqrt( (particlesNearShell(:,4) - aveXshell).^2 + (particlesNearShell(:,5) - aveYshell).^2 + (particlesNearShell(:,6) - aveZshell).^2 ) <= 24;
particlesInside = particlesNearShell(distP, :);

particlesFarShell = particlesNearShell(particlesNearShell(:,3)==30, :);
particlesNearShell = particlesNearShell(~distP, :);
particlesNearShell = particlesNearShell(particlesNearShell(:,3)==20, :);


gel = coordsF(coordsF(:,3) == 2, :);
fprintf(fidW, '#declare Geldata = array [%g]\n', length(gel));
fprintf(fidW, '{\n');

for i = 1:length(gel)
    
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end
    
    if i == length(gel)
        fprintf(fidW, '<%g,%g,%g>', gel(i,4:end));
    else
        fprintf(fidW, '<%g,%g,%g>,', gel(i,4:end));
    end
    
end

fprintf(fidW, '\n};\n');
fprintf(fidW, 'point_cloud(Geldata, radiusA)\n');

[~, m] = sort(coords(:,1));
% [~, mg] = sort(gel(:,1));
% [~, ms] = sort(shell(:,1));
% [~, mc1] = sort(chains1(:,1));
% [~, mc2] = sort(chains2(:,1));
% [~, mp] = sort(particles(:,1));
% gel = [gel(mg,1), gel(mg,2), gel(mg,3), gel(mg,4), gel(mg,5), gel(mg,6)];
% shell = [shell(ms,1), shell(ms,2), shell(ms,3), shell(ms,4), shell(ms,5), shell(ms,6)];
% chains1 = [chains1(mc1,1), chains1(mc1,2), chains1(mc1,3), chains1(mc1,4), chains1(mc1,5), chains1(mc1,6)];
% chains2 = [chains2(mc2,1), chains2(mc2,2), chains2(mc2,3), chains2(mc2,4), chains2(mc2,5), chains2(mc2,6)];
% particles = [particles(mp,1), particles(mp,2), particles(mp,3), particles(mp,4), particles(mp,5), particles(mp,6)];
coordsN = [coords(m,1), coords(m,2), coords(m,3), coords(m,4), coords(m,5), coords(m,6)];

chains1Bonds = bonds(bonds(:,1) == 3, :);
chains2Bonds = bonds(bonds(:,1) == 4, :);
[type1, type2] = scanBonds(chains2Bonds, chains2);

spacing = 5;
ind = 0;
fprintf(fidW, '\n');

for i = 1:length(gelBonds)
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end    
    if sqrt( (coordsN(gelBonds(i,2),4) - coordsN(gelBonds(i,3),4)).^2 + (coordsN(gelBonds(i,2),5) - coordsN(gelBonds(i,3),5)).^2 ...
            + (coordsN(gelBonds(i,2),6) - coordsN(gelBonds(i,3),6)).^2 ) > 5
        continue;
    else
        fprintf(fidW, ' cylinder {<%g,%g,%g>, <%g,%g,%g>, %s}', ...
        coordsN(gelBonds(i,2), 4:end), coordsN(gelBonds(i,3), 4:end), 'radiusB');
    end
end

fprintf(fidW, '\ntexture{ pigment{ rgbf < .4, .4, 1, 0.975> }}\n}');
%fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Green }  }\ntexture { pigment{ Green }  }\nfinish {\n');
%fprintf(fidW, 'ambient .1\ndiffuse .5\nreflection{ .15 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

%=========================================================================%
%=========================================================================%
%=========================================================================%

fprintf(fidW, '\nunion{\n');
fprintf(fidW, '#declare GeldataM = array [%g]\n', length(gel_more));
fprintf(fidW, '{\n');

for i = 1:length(gel_more)
    
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end
    
    if i == length(gel_more)
        fprintf(fidW, '<%g,%g,%g>', gel_more(i,4:end));
    else
        fprintf(fidW, '<%g,%g,%g>,', gel_more(i,4:end));
    end
    
end

fprintf(fidW, '\n};\n');
fprintf(fidW, 'point_cloud(GeldataM, radiusA)\n');

spacing = 5;
ind = 0;
fprintf(fidW, '\n');

for i = 1:length(gelBondsMore)
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end    
    if sqrt( (coordsN(gelBondsMore(i,2),4) - coordsN(gelBondsMore(i,3),4)).^2 + (coordsN(gelBondsMore(i,2),5) - coordsN(gelBondsMore(i,3),5)).^2 ...
            + (coordsN(gelBondsMore(i,2),6) - coordsN(gelBondsMore(i,3),6)).^2 ) > 5
        continue;
    else
        fprintf(fidW, ' cylinder {<%g,%g,%g>, <%g,%g,%g>, %s}', ...
        coordsN(gelBondsMore(i,2), 4:end), coordsN(gelBondsMore(i,3), 4:end), 'radiusB');
    end
end

fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Green }  }\ntexture { pigment{ Green }  }\nfinish {\n');
fprintf(fidW, 'ambient .1\ndiffuse .5\nreflection{ .15 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

%=========================================================================%
%=========================================================================%
%=========================================================================%

fprintf(fidW, '\nunion{\n');
fprintf(fidW, '#declare Shelldata = array [%g]\n', length(shell));
fprintf(fidW, '{\n');

for i = 1:length(shell)
    
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end
    
    if i == length(shell)
        fprintf(fidW, '<%g,%g,%g>', shell(i,4:end));
    else
        fprintf(fidW, '<%g,%g,%g>,', shell(i,4:end));
    end
    
end

fprintf(fidW, '\n};\n');
fprintf(fidW, 'point_cloud(Shelldata, radiusC)\n');


fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Blue }  }\ntexture { pigment{ Blue }  }\nfinish {\n');
fprintf(fidW, 'ambient .15\ndiffuse .6\nreflection{ .05 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fidW, '\nunion{\n');
fprintf(fidW, '#declare chains1data = array [%g]\n', length(chains1));
fprintf(fidW, '{\n');

for i = 1:length(chains1)
    
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end
    
    if i == length(chains1)
        fprintf(fidW, '<%g,%g,%g>', chains1(i,4:end));
    else
        fprintf(fidW, '<%g,%g,%g>,', chains1(i,4:end));
    end
    
end

fprintf(fidW, '\n};\n');
fprintf(fidW, 'point_cloud(chains1data, radiusA)\n');

spacing = 5;
ind = 0;
fprintf(fidW, '\n');

for i = 1:length(type1)
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end    
    if sqrt( (coordsN(type1(i,2),4) - coordsN(type1(i,3),4)).^2 + (coordsN(type1(i,2),5) - coordsN(type1(i,3),5)).^2 ...
            + (coordsN(type1(i,2),6) - coordsN(type1(i,3),6)).^2 ) > 5
        continue;
    else
        fprintf(fidW, ' cylinder {<%g,%g,%g>, <%g,%g,%g>, %s}', ...
        coordsN(type1(i,2), 4:end), coordsN(type1(i,3), 4:end), 'radiusB');
    end
end

fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Gray }  }\ntexture { pigment{ Gray }  }\nfinish {\n');
fprintf(fidW, 'ambient .15\ndiffuse .6\nreflection{ .15 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fidW, '\nunion{\n');
fprintf(fidW, '#declare chains2data = array [%g]\n', length(chains2));
fprintf(fidW, '{\n');

for i = 1:length(chains2)
    
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end
    
    if i == length(chains2)
        fprintf(fidW, '<%g,%g,%g>', chains2(i,4:end));
    else
        fprintf(fidW, '<%g,%g,%g>,', chains2(i,4:end));
    end
    
end

fprintf(fidW, '\n};\n');
fprintf(fidW, 'point_cloud(chains2data, radiusA)\n');

spacing = 5;
ind = 0;
fprintf(fidW, '\n');

for i = 1:length(type2)
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end    
    if sqrt( (coordsN(type2(i,2),4) - coordsN(type2(i,3),4)).^2 + (coordsN(type2(i,2),5) - coordsN(type2(i,3),5)).^2 ...
            + (coordsN(type2(i,2),6) - coordsN(type2(i,3),6)).^2 ) > 5
        continue;
    else
        fprintf(fidW, ' cylinder {<%g,%g,%g>, <%g,%g,%g>, %s}', ...
        coordsN(type2(i,2), 4:end), coordsN(type2(i,3), 4:end), 'radiusB');
    end
end

fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Red }  }\ntexture { pigment{ Red }  }\nfinish {\n');
fprintf(fidW, 'ambient .15\ndiffuse .6\nreflection{ .10 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

%=========================================================================%
%=========================================================================%
%=========================================================================%

fprintf(fidW, '\nunion{\n');
fprintf(fidW, '#declare Particlesdata0 = array [%g]\n', length(particlesNearShell));
fprintf(fidW, '{\n');

for i = 1:length(particlesNearShell)
    
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end
    
    if i == length(particlesNearShell)
        fprintf(fidW, '<%g,%g,%g>', particlesNearShell(i,4:end));
    else
        fprintf(fidW, '<%g,%g,%g>,', particlesNearShell(i,4:end));
    end
    
end

fprintf(fidW, '\n};\n');
fprintf(fidW, 'point_cloud(Particlesdata0, radiusC)\n');


fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Magenta }  }\ntexture { pigment{ Magenta }  }\nfinish {\n');
fprintf(fidW, 'ambient .15\ndiffuse .6\nreflection{ .05 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

%=========================================================================%
%=========================================================================%
%=========================================================================%

%=========================================================================%
%=========================================================================%
%=========================================================================%

fprintf(fidW, '\nunion{\n');
fprintf(fidW, '#declare Particlesdata1 = array [%g]\n', length(particlesFarShell));
fprintf(fidW, '{\n');

for i = 1:length(particlesFarShell)
    
    ind = ind + 1;
    if ind == spacing
       fprintf(fidW, '\n');
       ind = 1;
    end
    
    if i == length(particlesFarShell)
        fprintf(fidW, '<%g,%g,%g>', particlesFarShell(i,4:end));
    else
        fprintf(fidW, '<%g,%g,%g>,', particlesFarShell(i,4:end));
    end
    
end

fprintf(fidW, '\n};\n');
fprintf(fidW, 'point_cloud(Particlesdata1, radiusC)\n');


fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Yellow }  }\ntexture { pigment{ Yellow }  }\nfinish {\n');
fprintf(fidW, 'ambient .15\ndiffuse .6\nreflection{ .05 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

%=========================================================================%
%=========================================================================%
%=========================================================================%

%=========================================================================%
%=========================================================================%
%=========================================================================%
if ~isempty(particlesInside)
    
    fprintf(fidW, '\nunion{\n');
    fprintf(fidW, '#declare Particlesdata2 = array [%g]\n', length(particlesInside));
    fprintf(fidW, '{\n');

    for i = 1:length(particlesInside)

        ind = ind + 1;
        if ind == spacing
           fprintf(fidW, '\n');
           ind = 1;
        end

        if i == length(particlesInside)
            fprintf(fidW, '<%g,%g,%g>', particlesInside(i,4:end));
        else
            fprintf(fidW, '<%g,%g,%g>,', particlesInside(i,4:end));
        end

    end

    fprintf(fidW, '\n};\n');
    fprintf(fidW, 'point_cloud(Particlesdata2, radiusC)\n');


    fprintf(fidW, '\n#texture { Bright_Bronze scale 1 pigment{ Brown }  }\ntexture { pigment{ Brown }  }\nfinish {\n');
    fprintf(fidW, 'ambient .15\ndiffuse .6\nreflection{ .05 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

end
%=========================================================================%
%=========================================================================%
%=========================================================================%

fprintf(fidW, 'background { color White }\nlight_source {\n<75, -50, -120>\ncolor White\nfade_distance 50}');
fclose('all');