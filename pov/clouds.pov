// POV-Ray version 3.6/3.7 scenery file "sky03a.pov"
// author: Friedrich A. Lohmueller, Nov-2013
// homepage: http://www.f-lohmueller.de
// Shadow casting planes with clouds 
// A solution of the horizont problem without fog.
// inspired by Paul Koning
//--------------------------------------------------------------------------
#version 3.7; // 3.6;                        
global_settings{ assumed_gamma 1.0 }
#default{ finish{ ambient 0.1 diffuse 0.9 }}
//--------------------------------------------------------------------------

#include "colors.inc"
#include "textures.inc"
// camera -----------------------------------------------------------
#declare Camera_0 = camera {angle 80 
                            right    x*image_width/image_height
                            location  <0.0 , 1.0 ,-3.0>
                            look_at   <0.0 , 2.0 , 0.0>}
#declare Camera_1 = camera {angle 80   // shows shadow casting clouds
                            right    x*image_width/image_height
                            location  <0.0 , 200.0 ,-3.0>
                            look_at   <0.0 , 200.0 ,50.0>}
#declare Camera_2 = camera {angle 80   // the clouds from above
                            right    x*image_width/image_height
                            location  <0.0 , 2500 ,-3.0>
                            look_at   <0.0 , 2000 ,2000.0>}
camera{Camera_0}
// sun ---------------------------------------------------------------
light_source{<1500,2500,-2500>*100 color rgb<1,1,1> }

// sky ---------------------------------------------------------------

/* // optional with ground fog
sky_sphere{ pigment{ color rgb<0.15,0.28,0.75>*0.5}   } 
// ground fog at the horizon -----------------------------------------
fog{ fog_type   2
     distance   1000
     color      rgb<1,1,1>*0.9
     fog_offset 0.1
     fog_alt    30
     turbulence 1.8
   } //---------------------------------------------------------------
*/
// without ground fog
sky_sphere{
 pigment{ gradient y
          color_map{
          [0.0 color rgb<1,1,1>             ]
          [0.3 color rgb<0.18,0.28,0.75>*0.6]
          [1.0 color rgb<0.15,0.28,0.75>*0.5] }
          scale 1.05
          translate<0,-0.05,0>
    }
}  


// spherical cloud layer --------------------------------------------
#declare R_planet = 6000000 ;
#declare R_sky    = R_planet + 2000 ;

sphere{ <0, -R_planet, 0>, R_sky  hollow
       
        texture{ pigment{ bozo turbulence 0.75
                          octaves 6  omega 0.7 lambda 2  phase 0.00 //0.15*clock
                         color_map {
                          [0.00 color rgb <0.95, 0.95, 0.95> ]
                          [0.05 color rgb <1, 1, 1>*1.25 ]
                          [0.15 color rgb <0.85, 0.85, 0.85> ]
                          [0.55 color rgbt <1, 1, 1, 1>*1 ]
                          [1.00 color rgbt <1, 1, 1, 1>*1 ]
                         } // end color_map 
                         translate< 3, 0,-1>
                         scale <0.3, 0.4, 0.2>*3
                       } // end pigment
                              
                 #if (version = 3.7 )  finish {emission 1 diffuse 0}
                 #else                 finish { ambient 1 diffuse 0}
                 #end 
                 scale 3000
               } // end interior texture
    // no_shadow 
  }

// ground ------------------------------------------------------------
sphere{ <0, -R_planet, 0>, R_planet  
 
         texture{ pigment{color rgb<0.35,0.65,0.0>*0.8}
                  normal {bumps 0.75 scale 0.015}
                } // end of texture
      } // end of plane
//--------------------------------------------------------------------

//--------------------------------------------------------------------------
//---------------------------- objects in scene ----------------------------
//--------------------------------------------------------------------------
// a mirror sphere !!!!  
sphere{ <0,0,0>,0.6 scale <1,1,1> rotate<0,0,0> translate<0,0.7,0>
        texture{ Polished_Chrome }
      }
//--------------------------------------------------------------------------
