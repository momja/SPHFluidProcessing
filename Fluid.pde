// Max Omdal 2020

class Fluid {
    ArrayList<Particle> particles;
    ArrayList<ParticlePair> pairs; 
    Octree<Particle> octree;
    float timeStepSize = 0.001;
    int octantCapacity = 2;
    Vec3 gravity = new Vec3(0,-1000,0);
    Vec3 boundMax = new Vec3(1,2,1);
    Vec3 boundMin = new Vec3(-1,0,-1);

    // Fluid Parameters
    float K_smoothingRadius = 0.2;
    float K_stiff = 50;
    float K_stiffN = 20;
    float K_restDensity = 5;


    public Fluid(int particleCount) {
        particles = new ArrayList<Particle>();
        pairs = new ArrayList<ParticlePair>();
        Octant octPts = new Octant(new Vec3(0,1,0), new Vec3(2,2,2));
        octree = new Octree<Particle>(octPts, this.octantCapacity);
        for (int i = 0; i < particleCount; i++) {
            float randX = random(boundMin.x,boundMax.x);
            float randY = random(boundMin.y,boundMax.y);
            float randZ = random(boundMin.z,boundMax.z);
            addParticle(new Vec3(randX, randY, randZ));
        }

        // Populate Octree
        // updateOctree();
    }

    public void addParticle(Vec3 pos) {
        Particle newPart = new Particle(pos, this.K_smoothingRadius);
        particles.add(newPart);
    }

    private void createPair(Particle p1, Particle p2) {
        ParticlePair newPair = new ParticlePair(p1, p2);
        for (ParticlePair pair : pairs) {
            if (pair.equals(newPair)) {
                return;
            }
        }
        this.pairs.add(new ParticlePair(p1, p2));
    }

    private void updateOctree() {
        // Rebuild the octree
        octree.clear();
        for (Particle p : particles) {
            octree.insert(p);
        }
    }

    private void constrainToBounds(Particle p) {
        if (p.pos.x < boundMin.x) {
            // p.pos.x = boundMin.x;
            p.vel.add(new Vec3(5.1,0,0));
        }
        if (p.pos.x >= boundMax.x) {
            // p.pos.x = boundMax.x;
            p.vel.add(new Vec3(-5.1,0,0));
        }
        if (p.pos.y < boundMin.y) {
            // p.pos.y = boundMin.y;
            p.vel.add(new Vec3(0,5.1,0));
        }
        if (p.pos.y >= boundMax.y) {
            // p.pos.y = boundMax.y;
            p.vel.add(new Vec3(0,-5.1,0));
        }
        if (p.pos.z < boundMin.z) {
            // p.pos.z = boundMin.z;
            p.vel.add(new Vec3(0,0,5.1));
        }
        if (p.pos.z >= boundMax.z) {
            // p.pos.z = boundMax.z;
            p.vel.add(new Vec3(0,0,-5.1));
        }
    }

    public void update(float dt) {
        int timeSteps = ceil(dt/timeStepSize);
        float finalStepSize = dt % timeStepSize;

        this.pairs = new ArrayList<ParticlePair>();

        for (int i = 0; i < timeSteps; i++) {
            float stepSize = timeStepSize;
            if (i == timeSteps - 1) stepSize = finalStepSize;
            
            // TODO: Numerical integration goes here
            for (Particle p : particles) {
                p.vel = p.pos.minus(p.oldPos);
                constrainToBounds(p);
                p.oldPos = p.pos;
                p.vel.add(gravity.times(stepSize));
                p.pos.add(p.vel.times(stepSize));
                p.density = p.densityN = 0;
            }

            // updateOctree();

            // TODO : Use Octree for time speedup
            for (Particle p1 : particles) {
                // ArrayList<Particle> otherParticles = octree.inSameOctant(p1);
                for (Particle p2 : particles) {
                    if (p1 == p2)
                        continue;
                    if (p1.distance(p2) < K_smoothingRadius)
                        createPair(p1, p2);
                }
            }

            for (ParticlePair pair : pairs) {
                Particle A = pair.p1;
                Particle B = pair.p2;

                pair.q = 1 - A.distance(B) / K_smoothingRadius;
                pair.q2 = pow(pair.q, 2);
                pair.q3 = pow(pair.q, 3);
                
                A.density += pair.q2;
                B.density += pair.q2;
                A.densityN += pair.q3;
                B.densityN += pair.q3;
            }

            for (Particle p : particles) {
                p.pressure = K_stiff * (p.density - K_restDensity);
                p.pressureN = K_stiffN * p.densityN;
            }

            for (ParticlePair pair : pairs) {
                Particle A = pair.p1;
                Particle B = pair.p2;
                
                float pressure = A.pressure + B.pressure;
                float pressureN = A.pressureN + B.pressureN;
                float displace = (pressure*pair.q + pressureN*pair.q2) * pow(stepSize, 2);
                Vec3 a2bn = A.dirNormal(B);
                Vec3 displaceVec = a2bn.times(displace);
                A.pos.subtract(displaceVec);
                B.pos.add(displaceVec);
            }

        }
    }

    public void checkForCollisions(PShape[] rigidBodies) {
        // TODO : Implement
    }

    public void draw() {
        for (Particle particle : particles) {
            particle.draw();
        }

        if (debugMode) {
            octree.show(false);
            pushStyle();
            noFill();
            stroke(255,0,0,125);
            strokeWeight(0.5);
            for (Particle p : particles) {
                pushMatrix();
                translate(p.pos.x, p.pos.y, p.pos.z);
                sphere(K_smoothingRadius);
                popMatrix();
            }
            popStyle();
        }
    }
}