//
//  fluid.hpp
//  FluidSimulation
//
//  Created by Noé on 30/06/2021.
//

#ifndef fluid_hpp
#define fluid_hpp

#include <stdio.h>
#include "../utils/utils.hpp"

int IX(int x, int y);

class Fluid {
public:
    float dt;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;
    
    float *Vx0;
    float *Vy0;
    
    Fluid(){}
    Fluid(int diffusion, int viscosity, float dt);
    void Free();
    void AddDensity(int x, int y, float amount);
    void AddVelocity(int x, int y, float amountX, float amountY);
    void Step();
    
};

#endif /* fluid_hpp */
