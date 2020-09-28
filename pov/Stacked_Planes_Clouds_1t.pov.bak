// POV-Ray 3.7 Scene File "Stacked_Planes_Clouds_1.pov"
// author: Friedrich A. Lohmueller, Dec-2009 / Jan-2011
// email: Friedrich.Lohmueller_at_t-online.de
// homepage: http://www.f-lohmueller.de
//--------------------------------------------------------------------------
#version 3.6; // 3.7;
global_settings{ assumed_gamma 1.0 }
#default{ finish{ ambient 0.1 diffuse 0.9 }}
//--------------------------------------------------------------------------
#include "colors.inc"
#include "textures.inc"
#include "glass.inc"
#include "metals.inc"
#include "golds.inc"
#include "stones.inc"
#include "woods.inc"
#include "shapes.inc"
#include "shapes2.inc"
#include "functions.inc"
#include "math.inc"
#include "transforms.inc"
//-------------------------------------------------------------------------------------------------------<<<<
//------------------------------------------------------------- Camera_Position, Camera_look_at, Camera_Angle
#declare Camera_Number = 1 ;
//--------------------------------------------------------------------------------------------------------<<<<
#switch ( Camera_Number )
#case (0)
  #declare Camera_Position = < 0.00, 1.00, -7.00> ;  // front view
  #declare Camera_Look_At  = < 0.00, 1.00,  0.00> ;
  #declare Camera_Angle    =  65 ;
#break
#case (1)
  #declare Camera_Position = < 5.00, 1.00,-5.00> ;  // diagonal view up
  #declare Camera_Look_At  = < 0.00, 4.00,  0.00> ;
  #declare Camera_Angle    = 75 ;
#break
#else
  #declare Camera_Position = < 0.00, 1.00, -7.00> ;  // front view
  #declare Camera_Look_At  = < 0.00, 1.00,  0.00> ;
  #declare Camera_Angle    =  65 ;
#break
#end // of "#switch ( Camera_Number )" -----------------------------
//--------------------------------------------------------------------------------------------------------<<<<
camera{ location Camera_Position
        right    x*image_width/image_height
        angle    Camera_Angle
        look_at  Camera_Look_At
      }
//------------------------------------------------------------------------------------------------------------
// sun ---------------------------------------------------------------------
light_source{<3000,3500,-500> color White*0.8}           // sun light
light_source{ Camera_Position  color rgb<0.9,0.9,1>*0.2}  // flash light

// sky --------------------------------------------------------------
// set max_trace_level to  Number_of_Layers
global_settings { max_trace_level 65 }//(1...256) [default = 5]
//---------------------------------------------------------------------------------------
#include "sky/Stacked_Planes_Clouds_1.inc" 
//-------------------------------------------------------------------------------------// 
object{ Stacked_Planes_Clouds_1 (
         1000, // Cloud___Base_Height,    // height of the lowest layer of clouds 
         0.60, // Percentage_of_Blue_Sky, // 1 = no clouds, 0 = totally covered sky
         60,   // Number_of_Layers,       // number of planes 
         12,   // Distance_of_Layers,     // distance between planes  

         <0.356,0.35,0.41>*0.85,// Clouds_Base_Color,  // color of lower cloud parts
         <1,1,1>*1,           // Clouds_Top_Color,   // color of upper clouds parts
         pigment{wrinkles},     // clouds pattern - i.e. pigment{granite}, pigment{agate}, pigment{wrinkles}, ...
         0.6 //+0.1*(0.5 - 0.5*cos(Time*2*pi*100)), // pattern turbulence,         
         9,    // Pattern_Octaves,   // pattern modifier
         3,    // Pattern_Lambda,    // pattern modifier
         0.50, // Pattern_Omega,     // pattern modifier

         30,   // Pattern_Distance, // Moves only texture pattern up/down
         2500, // Clouds_Scale,     // Scaling for texture = cloud size - big: ~10000 - small: ~2500 
         0.06, // Dimmer Factor   0< ... <1 // Dimmer factor
         0.20, // Pattern border, 0 ~ 1 
         0.50, // Filter start,   0 ~ 1
         0.80  // Filter end,     0 ~ 1                  
       ) //----------------------------------------------------------------------------// 

 hollow  // no_shadow
 rotate  <0,100,0> 
 scale <-1,1,1>*0.3 
 translate< 0, 0, 500>
}
// adding background blue ----------------------------------------------------
sky_sphere { pigment { gradient <0,1,0>
                       color_map { [0.00 rgb <0.6,0.7,0.9>]
                                   [0.35 rgb <0.1,0.2,0.5>]
                                   [0.65 rgb <0.11,0.2,0.5>]
                                   [1.00 rgb <0.6,0.7,0.9>] 
                                 } 
                       scale 2         
                     } // end of pigment
           } //end of skysphere --------------------------------------------- 
// fog on the ground --------------------------------------------------------
fog { fog_type   2
      distance   500
      color      White  
      fog_offset 0.1
      fog_alt    3.5
      turbulence 1.8
    }
//----------------------------------------------------------------------------
// ground --------------------------------------------------------------------
plane { <0,1,0>, 0
        texture{ pigment{ color rgb<0.35,0.65,0.0>*0.72 }
	         normal { bumps 0.75 scale 0.015 }
                 finish { phong 0.1 }
               } // end of texture
      } // end of plane
//----------------------------------------------------------------------------
//---------------------------- objects in scene ------------------------------
//----------------------------------------------------------------------------

// test sphere to adjust the sun position 
sphere { <0,0,0>, 1 
         texture {  Polished_Chrome }  
         translate<0,1,0>  
       } // ----------------------------- 


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------






