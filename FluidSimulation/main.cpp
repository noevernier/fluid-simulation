#define OLC_PGE_APPLICATION
#include "utils/utils.hpp"

#include "fluid/fluid.hpp"


class FluidSimulation : public olc::PixelGameEngine {
public:
    FluidSimulation() {
        sAppName = "Fluid Simulation";
    }
    
public:
    Fluid *fluid;
    int px = 0;
    int py = 0;
    bool OnUserCreate() override {
        fluid = new Fluid(0.1,0.3,0.2);
        
        return true;
    }
    
    bool OnUserUpdate(float fElapsedTime) override {
        Clear(olc::BLACK);
        
        fluid->Step();
        
        fluid->AddDensity(GetMouseX(), GetMouseY(), 1);
        fluid->AddVelocity(GetMouseX(), GetMouseY(), GetMouseX()-px, GetMouseY()-py);
        
        for (int x = 0; x < ScreenWidth(); x++){
            for (int y = 0; y < ScreenHeight(); y++){
                int c = fluid->density[x+y*N] *255;
                Draw(x, y, olc::Pixel(c % 255, c % 255, c % 255));
            }
        }
        
        px = GetMouseX();
        py = GetMouseY();
        
        return true;
    }
    
    bool OnUserDestroy() override {
        fluid->Free();
        std::cout << "memory free" << std::endl;
        return true;
    }
};


int main(int argc, char const *argv[]) {
    FluidSimulation simulation;
	if (simulation.Construct(N, N, SIZE, SIZE))
        simulation.Start();

	return 0;
}
