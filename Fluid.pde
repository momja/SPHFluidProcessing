// Max Omdal 2020

class Fluid {
    ArrayList<Particle> particles;
    ArrayList<ParticlePair> pairs; 
    Octree<Particle> octree;
    float timeStepSize = 0.0003;
    int octantCapacity = 20;
    Vec3 gravity = new Vec3(0,-10000,0);
    Vec3 boundMax = new Vec3(0.2,1,0.5);
    Vec3 boundMin = new Vec3(-0.2,0,-0.5);

    // Fluid Parameters
    float K_smoothingRadius = 0.03;
    float K_stiff = 5;
    float K_stiffN = 8;
    float K_restDensity = 0.5;


    public Fluid(int particleCount) {
        particles = new ArrayList<Particle>();
        pairs = new ArrayList<ParticlePair>();
        Octant octPts = new Octant(new Vec3(0,0.5,0), new Vec3(10,10,10));
        octree = new Octree<Particle>(octPts, this.octantCapacity);
        for (int i = 0; i < particleCount; i++) {
            float randX = random(boundMin.x,boundMax.x);
            float randY = random(boundMin.y,boundMax.y);
            float randZ = random(boundMin.z,boundMax.z);
            addParticle(new Vec3(randX, randY, randZ));
        }

        // Populate Octree
        updateOctree();
    }

    public void addParticle(Vec3 pos) {
        Particle newPart = new Particle(pos, this.K_smoothingRadius);
        particles.add(newPart);
    }

    private void createPair(Particle p1, Particle p2, float dist) {
        ParticlePair newPair = new ParticlePair(p1, p2);
        newPair.q = 1 - dist / (K_smoothingRadius*2);
        this.pairs.add(newPair);
    }

    private void updateOctree() {
        // Rebuild the octree
        octree.clear();
        for (Particle p : particles) {
            octree.insert(p);
        }
    }

    private void constrainToBounds(Particle p) {
        float friction = 0.2;
        // if (p.pos.x < boundMin.x) {
        //     Vec3 normal = new Vec3(1,0,0);
        //     Vec3 vNormal = normal.times(dot(p.vel, normal));
        //     Vec3 vTangent = p.vel.minus(vNormal);
        //     Vec3 impulse = vNormal.minus(vTangent.times(friction));
        //     // p.pos.x = boundMin.x;
        //     p.vel.add(impulse);
        // }
        // if (p.pos.x >= boundMax.x) {
        //     Vec3 normal = new Vec3(-1,0,0);
        //     Vec3 vNormal = normal.times(dot(p.vel, normal));
        //     Vec3 vTangent = p.vel.minus(vNormal);
        //     Vec3 impulse = vNormal.minus(vTangent.times(friction));
        //     // p.pos.x = boundMax.x;
        //     p.vel.add(impulse);            
        // }
        if (p.pos.y < boundMin.y) {
            Vec3 normal = new Vec3(0,1,0);
            Vec3 vNormal = normal.times(dot(p.vel, normal));
            Vec3 vTangent = p.vel.minus(vNormal);
            Vec3 impulse = vNormal.minus(vTangent.times(friction));
            p.pos.y = boundMin.y;
            p.vel.add(impulse);
        }
        // if (p.pos.y >= boundMax.y) {
        //     Vec3 normal = new Vec3(0,-1,0);
        //     Vec3 vNormal = normal.times(dot(p.vel, normal));
        //     Vec3 vTangent = p.vel.minus(vNormal);
        //     Vec3 impulse = vNormal.minus(vTangent.times(friction));
        //     // p.pos.y = boundMax.y;
        //     p.vel.add(impulse);
        // }
        // if (p.pos.z < boundMin.z) {
        //     Vec3 normal = new Vec3(0,0,1);
        //     Vec3 vNormal = normal.times(dot(p.vel, normal));
        //     Vec3 vTangent = p.vel.minus(vNormal);
        //     Vec3 impulse = vNormal.minus(vTangent.times(friction));
        //     // p.pos.z = boundMin.z;
        //     p.vel.add(impulse);
        // }
        // if (p.pos.z >= boundMax.z) {
        //     Vec3 normal = new Vec3(0,0,-1);
        //     Vec3 vNormal = normal.times(dot(p.vel, normal));
        //     Vec3 vTangent = p.vel.minus(vNormal);
        //     Vec3 impulse = vNormal.minus(vTangent.times(friction));
        //     // p.pos.z = boundMax.z;
        //     p.vel.add(impulse);
        // }
        p.pos.add(p.vel.times(timeStepSize));
    }

    public void update(float dt) {
        int timeSteps = ceil(dt/timeStepSize);

        this.pairs = new ArrayList<ParticlePair>();

        for (int i = 0; i < timeSteps; i++) {
            float stepSize = timeStepSize;

            for (Particle p : particles) {
                p.vel = p.pos.minus(p.oldPos).times(1.f/stepSize);
                p.oldPos = p.pos;
                p.vel.add(gravity.times(stepSize));
                p.pos.add(p.vel.times(stepSize));
                p.density = p.densityN = 0;
            }

            updatePairs();
            doubleDensityRelaxation();
            resolveCollisions();
        }
    }

    public void updatePairs() {
        updateOctree();
        for (Particle p1 : particles) {
            ArrayList<Particle> otherParticles = octree.inSameOctant(p1);
            for (Particle p2 : otherParticles) {
                if (p1 == p2) continue;

                float dist = p1.distance(p2);
                if (dist < K_smoothingRadius*2)
                    createPair(p1, p2, dist);
            }
        }
    }

    public void doubleDensityRelaxation() {
        for (ParticlePair pair : pairs) {
                Particle A = pair.p1;
                Particle B = pair.p2;

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
                float displace = (pressure*pair.q + pressureN*pair.q2) * pow(timeStepSize, 2);
                Vec3 a2bn = A.dirNormal(B);
                Vec3 displaceVec = a2bn.times(displace/2);
                A.pos.subtract(displaceVec);
                B.pos.add(displaceVec);
            }
    }

    public void resolveCollisions() {
        for (Particle p : particles) {
            constrainToBounds(p);
        }
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