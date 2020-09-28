

filenameBonds= 'base';
filenamePar = 'particles'; 
outputName = 'newest_waves';
ninst = [0:10000:420000];

for index = ninst

    fidW1 = fopen(strcat(outputName, '_', (num2str(index)), '.pov'), 'w');
    fidW2 = fopen(strcat('3plet_', (num2str(index)), '.txt'), 'w');
    filenameBonds1 = strcat(filenameBonds, '.', num2str(index));
    filenamePar1 = strcat(filenamePar, '.', num2str(index));
    neigh_info = read_demsi(filenameBonds1);
    particle_info = read_demsi(filenamePar1);
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
    
    

    fprintf(fidW1, '#include "colors.inc"\n#include "textures.inc"\n#include "functions.inc"\n#include "glass.inc"\n#include "metals.inc"\n#include "stones.inc"\n#include "Stacked_Planes_Clouds_1.inc"\n#declare scaleP1 = 0.2;\n');
    fprintf(fidW1, 'camera {\nperspective\nlocation <-150000, 250000, 90000>\nlook_at <500000,250000,0>\nsky <0,0,1>\nangle 85\n}');
    
    %'\ncamera {\nlocation <250000, 0, 450000>\nlook_at <250000,250000,0>}\n');    
%    fprintf(fidW1, 'plane {z, -1000 pigment {color Blue }}\n');
%    fprintf(fidW1, 'union{\n');
    
%    fprintf(fidW1, '\n#macro point_cloud(A, r)\n#local _pce = dimension_size(A,1);\n#local c = 0;\n#while (c<_pce)\nsphere{A[c],r}\n#local c=c+1;\n#end\n#end\n\n#declare radiusA = 0.2;\n#declare radiusB = 0.2;\n#declare radiusC = 0.25;\n');
%    fprintf(fidW1, '\nunion{\n');

%    radius = 0.3;

    fprintf(fidW1, '//-------------------------------------------------------------------------------------//\n');
    fprintf(fidW1, 'object{ Stacked_Planes_Clouds_1 (\n');
             fprintf(fidW1, '3000, // Cloud___Base_Height,    // height of the lowest layer of clouds\n');
             fprintf(fidW1, '0.60, // Percentage_of_Blue_Sky, // 1 = no clouds, 0 = totally covered sky\n');
             fprintf(fidW1, '60,   // Number_of_Layers,       // number of planes\n');
             fprintf(fidW1, '10,   // Distance_of_Layers,     // distance between planes\n\n');

             fprintf(fidW1, '<0.356,0.35,0.41>*0.85,// Clouds_Base_Color,  // color of lower cloud parts\n');
             fprintf(fidW1, '<1,1,1>*1,           // Clouds_Top_Color,   // color of upper clouds parts\n');
             fprintf(fidW1, 'pigment{wrinkles},     // clouds pattern - i.e. pigment{granite}, pigment{agate}, pigment{wrinkles}, ...\n');
             fprintf(fidW1, '0.6 //+0.1*(0.5 - 0.5*cos(Time*2*pi*100)), // pattern turbulence,\n');
             fprintf(fidW1, '9,    // Pattern_Octaves,   // pattern modifier\n');
             fprintf(fidW1, '3,    // Pattern_Lambda,    // pattern modifier\n');
             fprintf(fidW1, '0.50, // Pattern_Omega,     // pattern modifier\n\n');

             fprintf(fidW1, '30,   // Pattern_Distance, // Moves only texture pattern up/down\n');
             fprintf(fidW1, '3000, // Clouds_Scale,     // Scaling for texture = cloud size - big: ~10000 - small: ~2500\n');
             fprintf(fidW1, '0.06, // Dimmer Factor   0< ... <1 // Dimmer factor\n');
             fprintf(fidW1, '0.20, // Pattern border, 0 ~ 1\n');
             fprintf(fidW1, '0.50, // Filter start,   0 ~ 1\n');
             fprintf(fidW1, '0.80  // Filter end,     0 ~ 1\n');
           fprintf(fidW1, ') //----------------------------------------------------------------------------//\n\n');

     fprintf(fidW1, 'hollow  // no_shadow\n');
     fprintf(fidW1, 'rotate  <0,0,100>\n'); 
     fprintf(fidW1, 'scale <-1,1,1>*30\n');
     fprintf(fidW1, 'translate< -%i, 0, 500>\n', index);
    fprintf(fidW1, '}\n\n');

    fprintf(fidW1, '// adding background blue ----------------------------------------------------\n');
    fprintf(fidW1, 'sky_sphere { pigment { gradient <0,0,1>\n');
                           fprintf(fidW1, 'color_map { [0.00 rgb <0.6,0.7,0.9>]\n');
                                       fprintf(fidW1, '[0.35 rgb <0.1,0.2,0.5>]\n');
                                       fprintf(fidW1, '[0.65 rgb <0.11,0.2,0.5>]\n');
                                       fprintf(fidW1, '[1.00 rgb <0.6,0.7,0.9>]\n');
                                     fprintf(fidW1, '}\n');
                           fprintf(fidW1, 'scale 2\n');
                         fprintf(fidW1, '} // end of pigment\n');
               fprintf(fidW1, '} //end of skysphere ---------------------------------------------\n');
    fprintf(fidW1, '// fog on the ground --------------------------------------------------------\n');
    fprintf(fidW1, 'fog { fog_type   2\n');
          fprintf(fidW1, 'distance   5000\n');
          fprintf(fidW1, 'color      White\n');
          fprintf(fidW1, 'fog_offset 0.1\n');
          fprintf(fidW1, 'fog_alt    3.5\n');
          fprintf(fidW1, 'turbulence 1.8\n');
          fprintf(fidW1, 'up<0,0,1>\n');
        fprintf(fidW1, '}\n\n');

%     fprintf(fidW1, 'plane{<0,0,1>, -100\n');
%           fprintf(fidW1, 'texture{pigment{ rgb <0.2, 0.2, 0.2> }\n');
%                   fprintf(fidW1, 'normal { bumps 0.08 scale <1,0.25,0.35>*30000 turbulence 0.6 }\n');
%                   fprintf(fidW1, 'finish { ambient 0.05 diffuse 0.55\n');
%                            fprintf(fidW1, 'brilliance 6.0 phong 0.8 phong_size 1200\n');
%                            fprintf(fidW1, 'reflection 0.6 }\n');
%                  fprintf(fidW1, '}\n\n');
%          fprintf(fidW1, '}\n');

    fprintf(fidW1, 'plane{<0,0,1>, -100\n')
    fprintf(fidW1, 'texture{pigment{ rgb <0.2, 0.2, 0.2> }\n')
    fprintf(fidW1, 'normal {\n')
      fprintf(fidW1, 'function {\n')
        fprintf(fidW1, 'f_ridged_mf(x, y, z, 0.1, 3.0, 7, 0.7, 0.7, 2)\n')
      fprintf(fidW1, '} 0.8\n')
      fprintf(fidW1, 'scale <0.13, 0.4, 0.13>*50000\n')
    fprintf(fidW1, '}\n')
    fprintf(fidW1, 'finish { ambient 0.05 diffuse 0.55\n')
    fprintf(fidW1, 'brilliance 6.0 phong 0.8 phong_size 1200\n')
    fprintf(fidW1, 'reflection 0.6 }\n')
    fprintf(fidW1, '}\n')
    fprintf(fidW1, 'interior {\n')
      fprintf(fidW1, 'ior 1.3\n')
      fprintf(fidW1, 'fade_distance 2\n')
      fprintf(fidW1, 'fade_power 500\n')
    fprintf(fidW1, '}\n')
    fprintf(fidW1, '}\n')

    display(sprintf('Timestep: %f', index));
    display(sprintf('Number of triplets: %g', length(triplets(:,1))));

    for w1 = 1:length(triplets(:,1))
        
        if w1 == 1
            
            fprintf(fidW1, '\nsmooth_triangle{<%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f> texture{pigment{ color rgb< 1, 1, 1>*0.85 } finish {ambient 0.2 diffuse 0.6 phong .75 phong_size scaleP1}}}\n', ...  %texture { White_Marble scale scaleP1 } }\n', ... 
                particle_info(particle_info(:,1) == triplets(w1,1), 3), particle_info(particle_info(:,1) == triplets(w1,1), 4), particle_info(particle_info(:,1) == triplets(w1,1), 5), ...
                0, 0, 1, particle_info(particle_info(:,1) == triplets(w1,2), 3), particle_info(particle_info(:,1) == triplets(w1,2), 4), particle_info(particle_info(:,1) == triplets(w1,2), 5), ...
                0, 0, 1, particle_info(particle_info(:,1) == triplets(w1,3), 3), particle_info(particle_info(:,1) == triplets(w1,3), 4), particle_info(particle_info(:,1) == triplets(w1,3), 5), 0, 0, 1);
            
        end
        
        fprintf(fidW1, 'smooth_triangle{<%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f> texture{pigment{ color rgb< 1, 1, 1>*0.85 } finish {ambient 0.2 diffuse 0.6 phong .75 phong_size scaleP1}}}\n', ...  %texture { White_Marble scale scaleP1 } }\n', ... 
                particle_info(particle_info(:,1) == triplets(w1,1), 3), particle_info(particle_info(:,1) == triplets(w1,1), 4), particle_info(particle_info(:,1) == triplets(w1,1), 5), ...
                0, 0, 1, particle_info(particle_info(:,1) == triplets(w1,2), 3), particle_info(particle_info(:,1) == triplets(w1,2), 4), particle_info(particle_info(:,1) == triplets(w1,2), 5), ...
                0, 0, 1, particle_info(particle_info(:,1) == triplets(w1,3), 3), particle_info(particle_info(:,1) == triplets(w1,3), 4), particle_info(particle_info(:,1) == triplets(w1,3), 5), 0, 0, 1);
        
    end
    
%    fprintf(fidW1, '\ntexture{ pigment{ rgbf < .4, .4, 1, 0.975> }}\n}');
%    fprintf(fidW1, '\n#texture { Bright_Bronze scale 1 pigment{ Green }  }\ntexture { pigment{ Green }  }\nfinish {\n');
%    fprintf(fidW1, 'ambient .1\ndiffuse .5\nreflection{ .15 1 fresnel }\nspecular 1\nconserve_energy\n}\n}\n');

    %=========================================================================%
    %=========================================================================%
    %=========================================================================%

    fprintf(fidW1, 'background { color White }\nlight_source {\n<300000, 7000000, 1000000>\ncolor White\nfade_distance 100000}');
    fclose('all');
    
end