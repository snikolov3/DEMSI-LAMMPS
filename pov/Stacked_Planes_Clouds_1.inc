// POV-Ray 3.6/3.7 include file "Stacked_Planes_Clouds_1.inc"                                                                                                                                                // PoVRay 3.6 include File "wheel_000.inc"
// author: Friedrich A, Lohmueller, Dec-2009
//  http://www.f-lohmueller.de/
// email: Friedrich.Lohmueller_aT_t-online.de
//------------------------------------------------------------------------------
// This include file makes clouds by stacking many transparent planes
// with a comon cloud pattern.
// Original idea for this technique from
// Stacked Sky Planes macro &quot;stsky.inc&quot;  by Jaime Vives Piqueres, 2000
// See http://www.ignorancia.org/en/index.php?page=Stacked_Planes
// -------------------------------------- --------------------------------------- 
// 
// ----------------------------------------------------------------------------- 
#ifndef( Stacked_Planes_Clouds_1_Inc_Temp)
#declare Stacked_Planes_Clouds_1_Inc_Temp = version;
#version 3.6;

//------------------------------------------------------------------------------ /////////  
#macro Stacked_Planes_Clouds_1( 
        
         Cloud_Base_Height_,  // 16000;  // height of the lowest layer of clouds
         Percentage_of_Clouds_, // 0.53, // 0 = no clouds, 1 = totally covered sky
         Number_of_Layers_,   // 60,     // number of planes or spheres
         Distance_of_Layers_, // 50,     // distance between planes or spheres 

         Clouds_Base_Color_,  // <0.356,0.35,0.41>, // color of lower cloud parts
         Clouds_Top_Color_,   // <1,1,1>*1.1,       // color of upper clouds parts
         Clouds_Pattern_,     // pigment{wrinkles}  // clouds pattern - pigment pattern like pigment{granite}
         Pattern_Turbulence_, // 0.6 +0.1*(0.5 - 0.5*cos(Time*2*pi*100)); // pattern turbulence,
         Pattern_Octaves_,    // 9,   // pattern modifier
         Pattern_Lambda_,     // 3,   // pattern modifier
         Pattern_Omega_,      // 0.5, // pattern modifier

         Texture_Distance_,   // 30,   // Distance between successive pattern textures
         Pattern_Scale_,      // 5000; // Scaling for texture = cloud size - big: ~10000 - small: ~2500  
         Dimmer_Factor_,      // 0.58; // Dimmer_Factor
         P_Border_, // pattern border, 0 ~ 1 
         F_Start_,  // filter start,   0 ~ 1
         F_End_,    // filter end,     0 ~ 1 
       ) //------------------------------------------------------------------------------------ 
//--------------------------------------------------------------------------------------------- 
//---------------------------------------------------------------------------------------------
// note : don't forget to raise max_trace_level!!!
#local Nr = 1;
#local EndNr = Number_of_Layers_;
#local C_Step_ = ((1 - Percentage_of_Clouds_ - P_Border_) / EndNr);
//---------------------------------------------------------------------------------------------
union{

#while (Nr<= EndNr) //1)//

 #local Color1= Clouds_Base_Color_*(1-(Nr/EndNr))               +Clouds_Top_Color_*(Nr/EndNr);
 #local Color2= Clouds_Base_Color_*(1-(Nr/EndNr))*Dimmer_Factor_+Clouds_Top_Color_*(Nr/EndNr)*Dimmer_Factor_;

 plane{ <0,0,1>,3000  
   texture{
     pigment{ Clouds_Pattern_
              octaves    Pattern_Octaves_ 
              lambda     Pattern_Lambda_
              omega      Pattern_Omega_
              turbulence Pattern_Turbulence_
              color_map{ //-------------------------------------------------------------------------------- 
                [ 0.00                            Clear ]
                [ Percentage_of_Clouds_           rgbf 1 ]
                [ Percentage_of_Clouds_+P_Border_ rgb Color1 filter F_Start_+(F_End_-F_Start_)*(Nr/EndNr) ]
                [ 1.00                            rgb Color2 filter F_Start_+(F_End_-F_Start_)*(Nr/EndNr) ]
              } // end of color_map ----------------------------------------------------------------------- 
              scale Pattern_Scale_    
            }// end of pigment

     finish {ambient 0.9 diffuse 0} // Problem: no shadows on Clouds ...!!! 
     translate <0, Cloud_Base_Height_ + Texture_Distance_*Nr, 0> // moves only texture pattern !!!
   } // end of texture
    
   translate <0, 0, Cloud_Base_Height_ + Distance_of_Layers_ * Nr> 
 } // end texture 
 #local Percentage_of_Clouds_= Percentage_of_Clouds_ + C_Step_*0.5;
 #local Nr=Nr+1;
 #end  // end of loop
 
} // end of union 
// ---------------------------------------------------------------------------------------------
#end// of macro -----------------------------------------------------------------// end of macro
// ---------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
// sample: 
/*
//----------------------------------------------------------------------------- textures
//---------------------------------------------------------------------------------------
#include "Stacked_Planes_Clouds_1.inc" 
//-------------------------------------------------------------------------------------// 
object{ Stacked_Planes_Clouds_1 (
         1000, // Cloud___Base_Height,    // height of the lowest layer of clouds 
         0.60, // Percentage_of_Blue_Sky, // 1 = no clouds, 0 = totally covered sky
         60,   // Number_of_Layers,       // number of planes 
         12,   // Distance_of_Layers,     // distance between planes  

         <0.356,0.35,0.41>*0.85,// Clouds_Base_Color,  // color of lower cloud parts
         <1,1,1>*1.0,           // Clouds_Top_Color,   // color of upper clouds parts
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
 scale 1 
 rotate  <0, -10,0> 
 translate<0000,0,1060>
 // scale<3,1,3>
} 
// adding background blue 
sky_sphere { pigment { gradient <0,1,0>
                       color_map { [0.00 rgb <0.6,0.7,0.9>]
                                   [0.35 rgb <0.1,0.2,0.5>]
                                   [0.65 rgb <0.11,0.2,0.5>]
                                   [1.00 rgb <0.6,0.7,0.9>] 
                                 } 
                       scale 2         
                     } // end of pigment
           } //end of skysphere --------------------------------------------- 
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
*/


#version Stacked_Planes_Clouds_1_Inc_Temp;
#end
//------------------------------------- end of include file

