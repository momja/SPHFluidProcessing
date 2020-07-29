// Max Omdal 2020

boolean debugMode = false;
boolean paused = false;
PShader unlitShader;
Camera cam = new Camera();
boolean[] keys;
boolean[][] bitmap;
Fluid fluid;
PGraphics painting;
PImage watercolor;
PImage canvas;
ArrayList<Particle> particlesUnderMouse = new ArrayList<Particle>();

float curSec;

void setup() {
    size(1280, 960, P3D);
    cam.setPerspective();
    surface.setTitle("SPH Fluid Simulation [Max Omdal]");
    keys = new boolean[8];
    keys[0] = false; keys[1] = false; keys[2] = false; keys[3] = false;
    keys[4] = false; keys[5] = false; keys[6] = false; keys[7] = false;

    unlitShader = loadShader("unlit_frag.glsl", "unlit_vert.glsl");

    bitmap = binarizeImage("umn_logo.jpg", 150);
    fluid = createFluidToDrawImage(bitmap);

    painting = createGraphics(1350, 1080, P2D);
    painting.imageMode(CENTER);
    watercolor = loadImage("watercolor.png");
    canvas = loadImage("canvas.jpg");
}

void draw() {
    updatePixels();
    cam.update();
    background(100,100,130);
    lights();
    pointLight(200,200,200,0,1,0);
    if(!paused) {
        update(1.f/frameRate);
        fluid.draw();
    }

    drawPainting();
}

void drawPainting() {
    float width_2 = (float)this.bitmap[0].length/20;
    float height_2 = (float)this.bitmap.length/20;

    textureMode(IMAGE);
    
    beginShape();
    noStroke();
    texture(canvas);
    vertex(-height_2,-0.0001,width_2,0,0);
    vertex(height_2,-0.0001,width_2,this.canvas.width,0);
    vertex(height_2,-0.0001,-width_2,this.canvas.width,this.canvas.height);
    vertex(-height_2,-0.0001,-width_2,0,this.canvas.height);
    endShape();

    beginShape();
    noStroke();
    texture(this.painting);
    vertex(-height_2,0,width_2,0,0);
    vertex(height_2,0,width_2,this.painting.width,0);
    vertex(height_2,0,-width_2,this.painting.width,this.painting.height);
    vertex(-height_2,0,-width_2,0,this.painting.height);
    endShape();
}

void update(float dt) {
    fluid.update(0.012);

    if (mousePressed) {
        Ray3 camToPlane = getMouseCast();
        Vec3 intersectionPt = rayPlaneCollision(camToPlane, new Vec3(0,1,0), new Vec3(0,0,0));
        if (particlesUnderMouse.size() == 0 && intersectionPt != null) {
            // Get nearest point on plane. Make sure it's less than a max distance
            // Get points near pt
            ArrayList<Particle> nearParticles = fluid.getParticlesAround(intersectionPt, 0.1);
            HashSet<Particle> checkedParticles = new HashSet<Particle>();
            for (Particle p : nearParticles) {
                if (checkedParticles.contains(p)) continue;
                checkedParticles.add(p);
                particlesUnderMouse.add(p);
            }
        }

        // Update positions of particles being pulled by mouse
        for (Particle p : particlesUnderMouse) {
            Vec3 vDir = intersectionPt.minus(p.pos);
            if (vDir.length() >= 0.5) {
                p.pos.add(vDir.times(0.3));
            }
        }
    } else {
        particlesUnderMouse = new ArrayList<Particle>();
    }

    updatePainting(dt);
}

void updatePainting(float dt) {
    // if (curSec == second()) return;
    // curSec = second();
    painting.beginDraw();

    float width_2 = (float)this.bitmap[0].length/20;
    float height_2 = (float)this.bitmap.length/20;
    // insert watercolor mark at each ball
    for (Particle p : fluid.particles) {

        float xCoord = map(p.pos.x, -height_2, height_2, 0,1350);
        float yCoord = map(p.pos.z, width_2, -width_2, 0,1080);

        painting.image(watercolor, xCoord, yCoord, 20,20);
    }
    painting.endDraw();
}

Fluid createFluidToDrawImage(boolean[][] bitmap) {
    int totalParticles = 0;
    for (int i = 0; i < bitmap[0].length; i++) {
        for (int j = 0; j < bitmap.length; j++) {
            totalParticles += bitmap[j][i] ? 1 : 0;
        }
    }

    Fluid fluid = new Fluid(floor(totalParticles*1.25));

    // Rearrange particles to be in right spots
    int partCount = 0;
    for (int i = 0; i < bitmap[0].length; i++) {
        for (int j = 0; j < bitmap.length; j++) {
            if (!bitmap[j][i]) continue;
            Particle p = fluid.particles.get(partCount);
            p.pos = new Vec3((float)j/10 - (float)bitmap.length/20, 0.01, (float)i/10 - (float)bitmap[0].length/20);
            p.pos.add(new Vec3(random(-0.01,0.01),random(0,0.01),random(-0.01,0.01)));
            p.oldPos = p.pos;

            if (partCount % 4 == 0) {
                Particle p2 = fluid.particles.get(totalParticles+partCount/4);
                p2.pos = new Vec3((float)j/10 - (float)bitmap.length/20, 0.1, (float)i/10 - (float)bitmap[0].length/20);
                p2.pos.add(new Vec3(random(-0.01,0.01),random(0,0.01),random(-0.01,0.01)));
                p2.oldPos = p2.pos;
            }

            partCount++;
        }
    }
    return fluid;
}