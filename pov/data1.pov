// Persistence of Vision Ray Tracer Scene Description File
// File: waves_test_b
// Desc: waves.inc water wave simulation test scene
//       heightfield version
// Date: November 2001
// Auth: Christoph Hormann <chris_hormann@gmx.de>


#version 3.5;

global_settings {
  assumed_gamma 1.0
  max_trace_level 15
}

#include "colors.inc"

#default { finish { ambient 0.07 }}

camera {
  location  <-1.1,  -0.6,  0.15>
  direction y
  sky       z
  up        z
  right     (4/3)*x
  look_at   <0.0, 0.0, 0.05>
  angle 43
}


light_source {
  <-8.0,   6.2,   4.5>
  color rgb 1.3
}


#declare Scale_Z = 0.02;

#declare M_Watx =
material {
  texture {
    pigment {
      color rgbt <0.4,0.5,0.3, 0.9>
    }
    finish {
      diffuse 0.3
      ambient 0.0

      reflection {
        0.03, 1.0
        fresnel on
      }

      conserve_energy

      specular 0.4
      roughness 0.003
    }
  }
  interior {
    ior 1.33

    fade_distance 0.3
    fade_power 1001
    fade_color <0.4, 0.5, 0.47>
  }
}

#include "waves.inc"

#declare MAX_WAVES = 60;


#declare fn_Water=
function {
  Waves(pi/2, 5.2, clock*7, 700) * 0.024 + 0.5
}


#declare Siz=5;
#declare Height=0.1;

height_field {

  function 300, 300 { fn_Water(x, y, 0) }

  smooth

  rotate -90*x
  translate <-0.5,-0.5,0.5>
  scale <Siz,Siz,Height>

  material {
    M_Watx
  }

  hollow on
}

difference {

  box { < -(Siz/2+0.2), -(Siz/2+0.2),-1>, <(Siz/2+0.2), (Siz/2+0.2),0.05> }
  box { < -(Siz/2-0.1), -(Siz/2-0.1),-0.99>, <(Siz/2-0.1), (Siz/2-0.1),0.2> }

  texture {
    pigment {
      color rgb <0.9, 0.88, 0.85>
    }
  }
}

cylinder { -z, z*0.35, 0.2 translate <0.2, -0.2, 0>

  texture {
    pigment {
      color rgb <0.9, 0.88, 0.85>
    }
  }
}

plane {
  z, -1

  texture {
    pigment {
      color rgb <0.9, 0.88, 0.85>
    }
  }
  hollow on
}


fog{
   fog_type 2
   fog_alt 0.4
   fog_offset 0.4
   color rgbt <0.85,0.9,1.0, 0.0>
   distance 300
   up z
}

sphere {
  <0, 0, 0>, 1
  texture {
    pigment {
      gradient z
      color_map {
        [0.0 color rgb < 0.680, 0.735, 0.976 >]
        [0.2 color rgb < 0.350, 0.500, 0.950 >]
      }
    }
    finish {
      diffuse 0
      ambient 1
    }
  }
  scale<1000, 1000, 200>
  no_shadow
  hollow on
}



