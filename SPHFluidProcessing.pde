// Max Omdal 2020

boolean debugMode = false;
boolean paused = false;
PShader unlitShader;
Camera cam = new Camera();
boolean[] keys;

Fluid fluid;

void setup() {
    size(1280, 960, P3D);
    cam.setPerspective();
    surface.setTitle("SPH Fluid Simulation [Max Omdal]");
    keys = new boolean[8];
    keys[0] = false; keys[1] = false; keys[2] = false; keys[3] = false;
    keys[4] = false; keys[5] = false; keys[6] = false; keys[7] = false;

    unlitShader = loadShader("unlit_frag.glsl", "unlit_vert.glsl");

    fluid = new Fluid(80);
}

void draw() {
    cam.update();
    background(210,210,220);
    lights();
    pointLight(100,100,180,0,1,0);
    if(!paused) {
        update(1.f/frameRate);
    }

    fluid.draw();
}

void update(float dt) {
    fluid.update(0.015);
}