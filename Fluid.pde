// Max Omdal 2020
import java.util.HashSet;

class Fluid {
    ArrayList<Particle> particles;
    ArrayList<ParticlePair> pairs; 
    Octree<Particle> octree;
    float timeStepSize = 0.006;
    int octantCapacity = 30;
    Vec3 gravity = new Vec3(0,-100,0);
    Vec3 boundMax = new Vec3(0.5,100,0.5);
    Vec3 boundMin = new Vec3(-0.5,0,-0.5);

    // Fluid Parameters
    float K_smoothingRadius = 0.1;
    float K_stiff = 30;
    float K_stiffN = 35;
    float K_restDensity = 8;
    float K_spring = 0.01;
    float K_springRestLength = 10;

    PShape particleShape;


    public Fluid(int particleCount) {
        particles = new ArrayList<Particle>();
        pairs = new ArrayList<ParticlePair>();
        Octant octPts = new Octant(new Vec3(0,2,0), new Vec3(10,20,10));
        octree = new Octree<Particle>(octPts, this.octantCapacity);
        this.particleShape = createShape(SPHERE, 0.02);
        this.particleShape.setStrokeWeight(0);
        this.particleShape.setFill(color(200,0,10));
        for (int i = 0; i < particleCount; i++) {
            float randX = random(-1, 1);
            float randY = random(-1, 1);
            float randZ = random(-1, 1);
            addParticle(new Vec3(randX, randY, randZ));
        }
        // Populate Octree
        updateOctree();
    }

    public void addParticle(Vec3 pos) {
        Particle newPart = new Particle(pos, this.K_smoothingRadius);
        newPart.particleShape = this.particleShape;
        particles.add(newPart);
    }

    public ArrayList<Particle> getParticlesAround(Vec3 pt, float radius) {
        Particle tmpParticle = new Particle(pt, radius);
        ArrayList<Particle> nearParticles = octree.inSameOctant(tmpParticle);
        return nearParticles;
    }

    private void createPair(Particle p1, Particle p2, float dist) {
        if (p1 == p2) return;
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
        // if (p.pos.x < boundMin.x) {
        //     p.pos.x = boundMin.x;
        // }
        // if (p.pos.x >= boundMax.x) {
        //     p.pos.x = boundMax.x;
        // }
        // if (p.pos.z < boundMin.z) {
        //     p.pos.z = boundMin.z;
        // }
        // if (p.pos.z >= boundMax.z) {
        //     p.pos.z = boundMax.z;
        // }
        if (p.pos.y < boundMin.y) {
            p.pos.y = boundMin.y;
        }
    }

    public void update(float dt) {
        float pairTime = 0;
        float doubleDensityTime = 0;
        float collisionTime = 0;

        int timeSteps = ceil(dt/timeStepSize);

        this.pairs = new ArrayList<ParticlePair>();

        for (int i = 0; i < timeSteps; i++) {
            float stepSize = timeStepSize;

            for (Particle p : particles) {
                // Apply gravity
                p.vel.add(gravity.times(stepSize));
            }

            // Modify Velocities with pairwise viscosity impulses
            applyViscosity();

            for (Particle p : particles) {
                p.oldPos = p.pos;
                p.pos.add(p.vel.times(stepSize));
                p.density = p.densityN = 0;
            }

            pairTime -= millis();
            updatePairs();
            pairTime += millis();
            doubleDensityTime -= millis();
            doubleDensityRelaxation();
            doubleDensityTime += millis();
            collisionTime -= millis();
            resolveCollisions();
            collisionTime += millis();

            for (Particle p : particles) {
                p.vel = p.pos.minus(p.oldPos).times(1.f/stepSize);
            }
        }
        // print("Time to make pairs (ms):\t");
        // println(pairTime);
        // print("Time to compute pos (ms):\t");
        // println(doubleDensityTime);
        // print("Time to compute coll (ms):\t");
        // println(collisionTime);
        // print("Avg number of neighbors:\t");
        // println(averageNoNbrs);
        // println();
    }

    float averageNoNbrs = 0;

    public void updatePairs() {
        averageNoNbrs = 0;
        updateOctree();
        for (Particle p1 : particles) {
            HashSet<Particle> checkedParticles = new HashSet<Particle>();
            ArrayList<Particle> otherParticles = octree.inSameOctant(p1);
            for (Particle p2 : otherParticles) {
                if (p1 == p2) continue;
                if (checkedParticles.contains(p2)) continue;
                checkedParticles.add(p2);
                averageNoNbrs++;
                float dist = p1.distance(p2);
                if (dist < K_smoothingRadius*2)
                    createPair(p1, p2, dist);
            }
        }
        averageNoNbrs /= particles.size();
    }

    public void applyViscosity() {
        for (ParticlePair pair : pairs) {
            Particle A = pair.p1;
            Particle B = pair.p2;
            float displace = timeStepSize*timeStepSize*K_spring*(1 - K_springRestLength/(K_smoothingRadius*2))*(K_springRestLength - pair.dist);
            Vec3 displaceVec = pair.btwnDir.times(displace);
            A.pos.subtract(displaceVec.times(0.5));
            B.pos.add(displaceVec.times(0.5));
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
                Vec3 displaceVec = a2bn.times(displace*0.5);
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