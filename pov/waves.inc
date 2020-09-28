// Persistence of Vision Ray Tracer Include File
// File: waves.inc
// Desc: water wave simulation
//       based on technique described in
//       http://www.first.gmd.de/persons/bwalter/html/report/report.html
//       and implementation from issue 12/2001 of 'PC-Magazin'
// Date: November 2001
// Auth: Christoph Hormann <chris_hormann@gmx.de>
//       http://www.tu-bs.de/~y0013390/
//
// Version 0.3 (experimental)
// This version is updated to be compatible with POV-Ray 3.5.
// no further changes were made in comparison to version 0.2.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                                                               %
// %   Feel free to use this file for your povray work, you are    %
// %   also allowed to modify and redistribute it as long as you   %
// %   leave these comments intact and highlight any changes you   %
// %   make.                                                       %
// %                                                               %
// %   I encourage you to contact me if you have ideas for         %
// %   improvements                                                %
// %                                                               %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "gamma.inc"

#declare gravity = 9.81;
#declare MAX_WAVES = 40;
#declare xalpha = 0.008;

#declare Seed=seed(9);

#macro directionEnergy(F, Theta)
  #local Temp = cos(Theta*0.5);
  #local Temp = Temp*Temp;

  #if (F >= freqPeak)
    #local s_sm = pow(F/freqPeak, -2.5);
    #local sm = 9.77;
  #else
    #local s_sm = pow(F/freqPeak,  5.0);
    #local sm = 6.97;
  #end

  #local S = s_sm * sm;
  (gamma( S+1 ) * pow(Temp, S) / (2.0*sqrt(pi)*gamma(S+0.5)) )
#end

#macro energylDFreq(F)
  ( xalpha * gravity * gravity / ( pow(2.0*pi, 4) * pow(F, 5) ) * exp((-5.0/4.0)* pow(F/freqPeak, -4.0) ) )
#end

#macro energy2DFreqAngle(F, Theta)
  ( energylDFreq( F ) * directionEnergy( F, Theta ) )
#end

#macro waveAmplitude( F, Theta, K )
 #local Omega = 10.0;
 ( sqrt( energy2DFreqAngle( F, Theta ) * gravity * pi * pi / ( K * F * Omega ) ) )
#end


// -------------------------------------------
// Main waves macro
//
// parameters:
//    - windDirection: main wind direction (radians)
//    - windSpeed: wind speed (mainly influencing the frequency distribution)
//    - time: time control for animation
//    - Scale: total scaling factor
//
// -------------------------------------------

#macro Waves(windDirection, windSpeed, time, Scale)

  #local freqPeak = 0.13 * gravity / windSpeed;

  (

    #debug "generating function...\n"

    #local nWaves = 0;

    #while ( nWaves < MAX_WAVES )

      #local Frequency = freqPeak + rand(Seed) - 0.5;

      #if ( Frequency > 0.0 )

        #local Wave_lambda = 2*pi / ( pow(2.0 * pi * Frequency, 2) / gravity );
        #local Wave_k = 2*pi / Wave_lambda;
        #local Wave_omega = sqrt(gravity * Wave_k );
        #local Wave_freq = Wave_omega / ( 2*pi );
        #local Wave_periode = 1.0 / Wave_freq;

        #local Wave_direction = ( rand(Seed)*2 ) * pi;
        #local Wave_dirX = cos( Wave_direction ) * 0.5;
        #local Wave_dirY = sin( Wave_direction ) * 0.5;

        #local phi0 = rand(Seed)*2*pi;
        #local a0 = waveAmplitude( Wave_freq, Wave_direction - windDirection, Wave_k );
        #local Wave_amplitude = a0 * cos( phi0 );
        #local Wave_phase = a0 * sin( phi0 );

        #if ( abs(a0) >= 0.0001 )
          #local nWaves = nWaves + 1;

          // --- Function code ---

          (
          Wave_amplitude *
          sin( Wave_phase + Wave_omega*Wave_freq*( time + Wave_dirX*x*Scale + Wave_dirY*y*Scale ) )
          )

          #if ( nWaves < MAX_WAVES )
            +
          #end
        #end
      #end
    #end

  )

#end
