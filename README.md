# Golf ballistics model

Simulate golf ball flight

![ballistics_plot](https://github.com/cagrell/golfmodel/blob/master/Image_1.PNG)

Based on ODE's from [MacDonald and Hanzely (1991) "The physics of the drive in golf"
    American Journal of Physics 59(3):213-218](
    https://www.researchgate.net/publication/253220357_The_physics_of_the_drive_in_golf)

__Input:__
```
  ------------------------------------------------------
  velocity              -   initial velocity  [m/s]
  launch_angle_deg      -   launch angle      [deg]
  off_center_angle_deg  -   aim               [deg]

  spin_rpm              -   ball spin         [rpm]
  spin_angle_deg        -   spin axis angle   [deg]

  windspeed             -   wind speed        [m/s] 
  windheading_deg       -   wind angle        [deg]

  Optional input:
  ------------------------------------------------------
  mass     -   mass of ball                  [kg]
  radius   -   radius of ball                [kg]
  rho      -   air density                   [kg/m^3]
  g        -   acceleration due to gravity   [kg]
```
