
// To render the scene from the command line, type:
// povray -W800 -N600 +A water.pov
// This will produce water.png, an 800x600 image. For faster rendering, you can reduce the resolution and not do anti-aliasing, for example:
// povray -W320 -H240 water.pov
// To really get the look you want, check out the POV ray docs: http://www.povray.org/documentation/
// It's pretty similar to OpenGL, so have fun with it!

#include "colors.inc"
#include "textures.inc"
#include "finish.inc"

background{Black}

// If you're not doing water, or just don't like the color, change it here
#declare Water = pigment
{
	color Blue transmit 0.7
}

// Camera. Pretty straightforward. Modify as you wish.
camera {
	location <0,0,-3>
	look_at <0,0,0>
}

// Add or modify the light sources as you wish
light_source { <10, 20, -10> color White }

// This is the main "blob" part. Output your particles as spheres. Look inside for more details.
// To really get the look you want, read up here: http://www.povray.org/documentation/view/3.6.1/71/
blob
{
	// You may want to change this depending on how "sticky" you want your blobs to be. Experiment to get the look you want.
	threshold .5

		// Each sphere is this format:
		// sphere { <x,y,z>, 1, 1 pigment {Water} }
		// x,y,z specify the center. The other numbers are radius and strength (see the povray docs).
		// Here, I've defined the four centers of the blobs
	sphere { <.75,-0.5,1>, 1, 1 pigment {Water} }
	sphere { <-.75,-0.5,1>, 1, 1 pigment {Water} }
	sphere { <0,0.4,1>, 1, 1 pigment {Water} }
	sphere { <1.1,1.1,1>, 1, 1 pigment {Water} }

	// Some material info to make it transparent, etc.
	finish {
		ambient 0.0
		diffuse 0.0
		specular 0.4
		roughness 0.003
		reflection { 0.003, 1.0 fresnel on }
	}
	interior { ior 1.33 }
}

// Here, you should add whatever surrounding geometry you want.
// This is just a big marble wall to give the droplets a backdrop
box {
	<-8, -8, 2>, <8,8,4>
	texture { White_Marble scale 50 }
}
