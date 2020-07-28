import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.HashSet; 
import java.util.HashSet; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class SPHFluidProcessing extends PApplet {

// Max Omdal 2020

boolean debugMode = false;
boolean paused = false;
PShader unlitShader;
Camera cam = new Camera();
boolean[] keys;

Fluid fluid;

public void setup() {
    
    cam.setPerspective();
    surface.setTitle("SPH Fluid Simulation [Max Omdal]");
    keys = new boolean[8];
    keys[0] = false; keys[1] = false; keys[2] = false; keys[3] = false;
    keys[4] = false; keys[5] = false; keys[6] = false; keys[7] = false;

    unlitShader = loadShader("unlit_frag.glsl", "unlit_vert.glsl");

    fluid = new Fluid(2000);
}

public void draw() {
    cam.update();
    background(210,210,220);
    lights();
    pointLight(100,100,180,0,1,0);
    if(!paused) {
        update(1.f/frameRate);
    }
    
    float drawTime = -millis();
    fluid.draw();
    drawTime += millis();
    // print("Draw Time (ms):\t");
    // println(drawTime);
}

public void update(float dt) {
    fluid.update(0.012f);
}
class Camera {
    Vec3 camLocation = new Vec3(0,0,0);
    Vec3 camLookAt = new Vec3(0,0,0);
    Vec3 camUp = new Vec3(0,-1,0);
    float radius = 3;
    int slider = 0;
    float theta = 0;
    float fov = 55;
    float nearPlaneW = 1 + 1.f/3;
    float nearPlaneH = 1;
    float nearPlaneDist = 1;
    float farPlaneDist = 1000;

    public void update() {
        if (keyPressed) {
            if (keyCode == UP) {
                slider++;
            } else if (keyCode == DOWN) {
                slider--;
            } else if (keyCode == LEFT) {
                theta -= 1.1f;
            } else if (keyCode == RIGHT) {
                theta += 1.1f;
            }
        }
        camLocation.x = cos(radians(theta))*radius;
        camLocation.y = PApplet.parseFloat(slider)/5;
        camLocation.z = sin(radians(theta))*radius;
        camera(camLocation.x, camLocation.y, camLocation.z,
               camLookAt.x,   camLookAt.y,   camLookAt.z,
               camUp.x,       camUp.y,       camUp.z);
    }

    public void setPerspective() {
         perspective(radians(fov), nearPlaneW/nearPlaneH, nearPlaneDist, farPlaneDist);
    }
}
// Max Omdal 2020


class Fluid {
    ArrayList<Particle> particles;
    ArrayList<ParticlePair> pairs; 
    Octree<Particle> octree;
    float timeStepSize = 0.003f;
    int octantCapacity = 10;
    Vec3 gravity = new Vec3(0,-1000,0);
    Vec3 boundMax = new Vec3(0.1f,1,10);
    Vec3 boundMin = new Vec3(-0.1f,0,-10);

    // Fluid Parameters
    float K_smoothingRadius = 0.04f;
    float K_stiff = 5;
    float K_stiffN = 8;
    float K_restDensity = 0.5f;

    PShape particleShape;


    public Fluid(int particleCount) {
        particles = new ArrayList<Particle>();
        pairs = new ArrayList<ParticlePair>();
        Octant octPts = new Octant(new Vec3(0,9,0), new Vec3(10,20,10));
        octree = new Octree<Particle>(octPts, this.octantCapacity);
        this.particleShape = createShape(SPHERE, 0.04f);
        this.particleShape.setStrokeWeight(0);
        this.particleShape.setFill(color(0,50,180));
        for (int i = 0; i < particleCount/10; i++) {
            for (int j = 0; j < 10; j++) {
                float randX = random(-0.01f, 0.01f);
                float randY = ((float)i)/30 + random(-0.01f, 0.01f);
                float randZ = ((float)j)/30 + random(-0.01f, 0.01f);
                addParticle(new Vec3(randX, randY, randZ));
            }
        }


        // Populate Octree
        updateOctree();
    }

    public void addParticle(Vec3 pos) {
        Particle newPart = new Particle(pos, this.K_smoothingRadius);
        newPart.particleShape = this.particleShape;
        particles.add(newPart);
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
        float friction = 0.2f;
        if (p.pos.x < boundMin.x) {
            Vec3 normal = new Vec3(1,0,0);
            Vec3 vNormal = normal.times(dot(p.vel, normal));
            Vec3 vTangent = p.vel.minus(vNormal);
            Vec3 impulse = vNormal.minus(vTangent.times(friction));
            // p.pos.x = boundMin.x;
            p.vel.add(impulse);
        }
        if (p.pos.x >= boundMax.x) {
            Vec3 normal = new Vec3(-1,0,0);
            Vec3 vNormal = normal.times(dot(p.vel, normal));
            Vec3 vTangent = p.vel.minus(vNormal);
            Vec3 impulse = vNormal.minus(vTangent.times(friction));
            // p.pos.x = boundMax.x;
            p.vel.add(impulse);            
        }
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
        if (p.pos.z < boundMin.z) {
            Vec3 normal = new Vec3(0,0,1);
            Vec3 vNormal = normal.times(dot(p.vel, normal));
            Vec3 vTangent = p.vel.minus(vNormal);
            Vec3 impulse = vNormal.minus(vTangent.times(friction));
            // p.pos.z = boundMin.z;
            p.vel.add(impulse);
        }
        if (p.pos.z >= boundMax.z) {
            Vec3 normal = new Vec3(0,0,-1);
            Vec3 vNormal = normal.times(dot(p.vel, normal));
            Vec3 vTangent = p.vel.minus(vNormal);
            Vec3 impulse = vNormal.minus(vTangent.times(friction));
            // p.pos.z = boundMax.z;
            p.vel.add(impulse);
        }
        p.pos.add(p.vel.times(timeStepSize));
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
                p.vel = p.pos.minus(p.oldPos).times(1.f/stepSize);
                p.oldPos = p.pos;
                p.vel.add(gravity.times(stepSize));
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
            strokeWeight(0.5f);
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
class Octant {
    Vec3 origin;
    Vec3 size;
    int capacity;
    float xmin, xmax, ymin, ymax, zmin, zmax;
    Vec3[] bounds = new Vec3[2];

    public Octant(Vec3 origin, Vec3 size) {
        this.origin = origin;
        this.size = size;

        this.xmin = origin.x - size.x/2;
        this.xmax = origin.x + size.x/2;
        this.ymin = origin.y - size.y/2;
        this.ymax = origin.y + size.y/2;
        this.zmin = origin.z - size.z/2;
        this.zmax = origin.z + size.z/2;
        bounds[0] = new Vec3(xmin,ymin,zmin);
        bounds[1] = new Vec3(xmax,ymax,zmax);
    }

    public boolean containsPoint(Vec3 p) {
        return (p.x >= xmin &&
                p.x < xmax &&
                p.y >= ymin &&
                p.y < ymax &&
                p.z >= zmin &&
                p.z < zmax);
    }

    public boolean rayIntersects(Ray3 ray) {
        // Credit goes to https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection

        // 1.) Check if ray is contained inside the OctantPoints
        if (containsPoint(ray.origin)) {
            return true;
        }
        // 2.) Check if ray passes through the oct
        float tmin, tmax, tymin, tymax, tzmin, tzmax;

        tmin = (bounds[ray.sign[0]].x - ray.origin.x) * ray.invDir.x;
        tmax = (bounds[1-ray.sign[0]].x - ray.origin.x) * ray.invDir.x;
        tymin = (bounds[ray.sign[1]].y - ray.origin.y) * ray.invDir.y;
        tymax = (bounds[1-ray.sign[1]].y - ray.origin.y) * ray.invDir.y;

        if ((tmin > tymax) || (tymin > tmax)) {
            return false;
        }
        if (tymin > tmin) {
            tmin = tymin;
        }
        if (tymax < tmax) {
            tmax = tymax;
        }

        tzmin = (bounds[ray.sign[2]].z - ray.origin.z) * ray.invDir.z; 
        tzmax = (bounds[1-ray.sign[2]].z - ray.origin.z) * ray.invDir.z; 
    
        if ((tmin > tzmax) || (tzmin > tmax)) {
            return false;
        }
        if (tzmin > tmin) {
            tmin = tzmin;
        }
        if (tzmax < tmax) {
            tmax = tzmax;
        }

        // if (tmax < 0.f) {
        //     return false;
        // }

        // if (ray.magnitude > 0 && ray.magnitude < tmax) {
        //     return false;
        // }

        return true; 
    }
}
// Max Omdal 2020

interface OctantInsertable {
    public Vec3 getPos();
    public boolean inOctant(Octant octant);
}
// Max Omdal 2020



class Octree<T extends OctantInsertable> {
    Octant bounds;
    int capacity;
    ArrayList<T> points;
    boolean divided = false;
    ArrayList<Octree<T>> children;
    int depth;
    int maxDepth = 6;

    public Octree(Octant bounds, int capacity) {
        this.bounds = bounds;
        this.capacity = capacity;
        this.points = new ArrayList<T>();
        this.depth = 0;
    }

    public void clear() {
        this.points = new ArrayList<T>();
        this.divided = false;
        this.children = new ArrayList<Octree<T>>();
    }

    public void insert(T p) {
        if (!p.inOctant(this.bounds)) {
            return;
        }

        if (this.depth == this.maxDepth || (!this.divided && this.points.size() < this.capacity)) {
            this.points.add(p);
        } else if (!this.divided) {
            this.points.add(p);
            subdivide();
        } else {
            for (Octree<T> child : children) {
                child.insert(p);
            }
        }
    }

    public void show(boolean withPoints) {
        if (this.divided) {
            for (Octree<T> child : children) {
                child.show(withPoints);
            }
        } else {
            pushMatrix();
            strokeWeight(0.8f);
            stroke(255,100,100);
            noFill();
            translate(bounds.origin.x, bounds.origin.y, bounds.origin.z);
            box(bounds.size.x, bounds.size.y, bounds.size.z);
            popMatrix();
        }
        if (withPoints) {
            pushStyle();
            strokeWeight(4);
            stroke(255,0,0);
            for (T p : points) {
                Vec3 pos = p.getPos();
                point(pos.x, pos.y, pos.z);
            }
            popStyle();
        }
    }

    public ArrayList<T> rayIntersectsOctants(Ray3 ray) {
        // Gets all octants that the ray is in and returns the points stored in those octants
        ArrayList<T> pointsInOctants = new ArrayList<T>();
        
        if (bounds.rayIntersects(ray)) {
            if (divided) {
                // Recursively divide
                for (Octree<T> child : children) {
                    pointsInOctants.addAll(child.rayIntersectsOctants(ray));
                }
            } else {
                // Just add the points in this oct
                pointsInOctants.addAll(this.points);
            }
        }

        return pointsInOctants;
    }

    public ArrayList<T> inSameOctant(T p) {
        ArrayList<T> pointsInOctants = new ArrayList<T>();

        if (p.inOctant(bounds)) {
            if (divided) {
                // Recursively divide
                // ArrayList<T> newPts;
                for (Octree<T> child : children) {
                    pointsInOctants.addAll(child.inSameOctant(p));
                }
            } else {
                // Just add the points in this oct
                pointsInOctants.addAll(this.points);
            }
        }

        return pointsInOctants;
    }

    private void subdivide() {
        Vec3 origin = bounds.origin;
        Vec3 size = bounds.size;
        Vec3 subdivideSize = size.times(0.5f);
        this.children = new ArrayList<Octree<T>>();

        Vec3[] octantSections = new Vec3[]{
            new Vec3(0.5f,0.5f,0.5f),
            new Vec3(-0.5f,0.5f,0.5f),
            new Vec3(-0.5f,-0.5f,0.5f),
            new Vec3(0.5f,-0.5f,0.5f),
            new Vec3(0.5f,0.5f,-0.5f),
            new Vec3(-0.5f,0.5f,-0.5f),
            new Vec3(-0.5f,-0.5f,-0.5f),
            new Vec3(0.5f,-0.5f,-0.5f)
        };
        
        for (int i = 0; i < 8; i++) {
            Octant octant = new Octant(origin.plus(subdivideSize.times(octantSections[i])), subdivideSize);
            Octree<T> child = new Octree<T>(octant, capacity);
            child.depth = this.depth + 1;
            this.children.add(child);
        }

        for (T p : points) {
            for (Octree<T> child : children) {
                child.insert(p);
            }
        }

        points = new ArrayList<T>();

        this.divided = true;
    }

}
// Max Omdal 2020

class Particle implements OctantInsertable {
    float radius = 0.04f;
    Float smoothingRadius;
    Vec3 pos;
    Vec3 oldPos;
    Vec3 vel;
    float density;
    float densityN;
    float pressure;
    float pressureN;

    PShape particleShape;


    public Particle(Vec3 pos, Float smoothingRadius) {
        this.pos = pos;
        this.oldPos = pos;
        this.vel = new Vec3();
        this.pressure = 0;
        this.pressureN = 0;
        this.density = 0;
        this.densityN = 0;
        this.smoothingRadius = smoothingRadius;
    }

    public float distance(Particle p) {
        return this.pos.distanceTo(p.pos);
    }

    public Vec3 dirNormal(Particle p) {
        Vec3 vec = p.pos.minus(this.pos);
        if (vec.length() > 0) {
            return vec.normalized();
        }
        return new Vec3();
    }

    public void draw() {
        push();
        if (this.particleShape != null) {
            translate(this.pos.x, this.pos.y, this.pos.z);
            shape(this.particleShape);
        } else {
            stroke(0,50,180);
            strokeWeight(20);
            point(this.pos.x, this.pos.y, this.pos.z);
        }
        pop();
    }

    /*
    OctantInsertable
    ----------------
    */

    @Override
    public boolean inOctant(Octant octant) {
        // Check to see if any part of particle's radius is in octant
        if (octant.containsPoint(this.pos)) {
            return true;
        }

        // float distSquared = this.smoothingRadius*this.smoothingRadius;

        // if (pos.x < octant.xmin) distSquared -= pow(pos.x - octant.xmin, 2);
        // else if (pos.x >= octant.xmax) distSquared -= pow(pos.x - octant.xmax, 2);
        // if (pos.y < octant.ymin) distSquared -= pow(pos.y - octant.ymin, 2);
        // else if (pos.y >= octant.ymax) distSquared -= pow(pos.y - octant.ymax, 2);
        // if (pos.z < octant.zmin) distSquared -= pow(pos.z - octant.zmin, 2);
        // else if (pos.z >= octant.zmax) distSquared -= pow(pos.z - octant.zmax, 2);

        // return distSquared > 0;

        float x = max(octant.xmin, min(pos.x, octant.xmax));
        float y = max(octant.ymin, min(pos.y, octant.ymax));
        float z = max(octant.zmin, min(pos.z, octant.zmax));

        Vec3 closestPt = new Vec3(x,y,z);
        float distance = pos.distanceTo(closestPt);

        return distance < this.smoothingRadius;
    }

    @Override
    public Vec3 getPos() {
        return this.pos;
    } 
}

class ParticlePair {
    Particle p1;
    Particle p2;

    float q;
    float q2;
    float q3;

    public ParticlePair(Particle p1, Particle p2) {
        this.p1 = p1;
        this.p2 = p2;
    }

    public boolean equals(ParticlePair pair) {
        if ((this.p1 == pair.p1 && this.p2 == pair.p2) ||
            (this.p2 == pair.p1 && this.p1 == pair.p2)) {
            return true;
        }
        return false;
    }
}
class Ray3 {
    Vec3 origin;
    Vec3 direction;
    float magnitude = -1;
    Vec3 invDir;
    int[] sign = new int[3];
    
    public Ray3(Vec3 origin, Vec3 dir, float magnitude) {
        this.origin = origin;
        this.direction = dir.normalized();
        this.magnitude = magnitude;

        invDir = new Vec3(1.f/dir.x, 1.f/dir.y, 1.f/dir.z);
        sign[0] = PApplet.parseInt(invDir.x < 0);
        sign[1] = PApplet.parseInt(invDir.y < 0);
        sign[2] = PApplet.parseInt(invDir.z < 0);
    }

    public Ray3(Vec3 origin, Vec3 dir) {
        this.origin = origin;
        this.direction = dir.normalized();

        invDir = new Vec3(1.f/dir.x, 1.f/dir.y, 1.f/dir.z);
        sign[0] = PApplet.parseInt(invDir.x < 0);
        sign[1] = PApplet.parseInt(invDir.y < 0);
        sign[2] = PApplet.parseInt(invDir.z < 0);
    }

    public void debugDraw() {
        pushStyle();
        stroke(255,255,0);
        strokeWeight(1);
        if (magnitude > 0) line(origin.x, origin.y, origin.z, origin.x+direction.x*magnitude, origin.y+direction.y*magnitude, origin.z+direction.z*magnitude);
        else line(origin.x, origin.y, origin.z, origin.x+direction.x*10, origin.y+direction.y*10, origin.z+direction.z*10);
        popStyle();
    }

    public Vec3 pointAtTime(float t) {
        return origin.plus(direction.times(t));
    }
}
public Ray3 getMouseCast() {
  Vec3 w = cam.camLookAt.minus(cam.camLocation).normalized();
  Vec3 u = cross(w, cam.camUp).normalized();
  Vec3 v = cross(u, w).normalized();

  w.mul(-1);

  float m3dx = map(mouseX, 0, width, -cam.nearPlaneW/2, cam.nearPlaneW/2);
  float m3dy = map(mouseY, 0, height, -cam.nearPlaneH/2, cam.nearPlaneH/2);
  float m3dz = -1;

  float m3dx_world = m3dx*u.x + m3dy*v.x + m3dz*w.x + cam.camLocation.x;
  float m3dy_world = m3dx*u.y + m3dy*v.y + m3dz*w.y + cam.camLocation.y;
  float m3dz_world = m3dx*u.z + m3dy*v.z + m3dz*w.z + cam.camLocation.z;

  Vec3 m_world = new Vec3(m3dx_world, m3dy_world, m3dz_world);
  Vec3 rayDir = m_world.minus(cam.camLocation);
  rayDir.normalize();
  return new Ray3(cam.camLocation, rayDir);
}

public void keyPressed() {
    if (key == ' ')
        paused = !paused;
    if (key == 'p')
        debugMode = !debugMode;
    if (key=='w')
        keys[0]=true;
    if (key=='a')
        keys[1]=true;
    if (key=='s')
        keys[2]=true;
    if (key=='d')
        keys[3]=true;
    if (keyCode==LEFT)
        keys[4]=true;
    if (keyCode==RIGHT)
        keys[5]=true;
    if (keyCode==UP)
        keys[6]=true;
    if (keyCode==DOWN)
        keys[7]=true;
}

public void keyReleased() {
    if (key=='w')
        keys[0]=false;
    if (key=='a')
        keys[1]=false;
    if (key=='s')
        keys[2]=false;
    if (key=='d')
        keys[3]=false;
    if (keyCode==LEFT)
        keys[4]=false;
    if (keyCode==RIGHT)
        keys[5]=false;
    if (keyCode==UP)
        keys[6]=false;
    if (keyCode==DOWN)
        keys[7]=false;
}

public void mouseClicked() {
}
//Vector Library [2D]
//CSCI 5611 Vector 3 Library [Incomplete]

//Instructions: Add 3D versions of all of the 2D vector functions
//              Vec3 must also support the cross product.
public class Vec2 {
  public float x, y;

  public Vec2() {
    this.x = 0;
    this.y = 0;
  }
  
  public Vec2(float x, float y){
    this.x = x;
    this.y = y;
  }

  public Vec2(Vec2 v) {
    this.x = v.x;
    this.y = v.y;
  }
  
  public String toString(){
    return "(" + x + ", " + y + ")";
  }
  
  public float length(){
    return sqrt(x*x + y*y);
  }

  public float lengthSqr() {
    return (x*x + y*y);
  }
  
  public Vec2 plus(Vec2 rhs){
    return new Vec2(rhs.x+this.x, rhs.y+this.y);
  }
  
  public void add(Vec2 rhs){
    this.x += rhs.x;
    this.y += rhs.y;
  }
  
  public Vec2 minus(Vec2 rhs){
    return new Vec2(this.x-rhs.x, this.y-rhs.y);
  }
  
  public void subtract(Vec2 rhs){
    this.x -= rhs.x;
    this.y -= rhs.y;
  }
  
  public Vec2 times(float rhs){
    return new Vec2(this.x*rhs, this.y*rhs);
  }
  
  public void mul(float rhs){
    this.x *= rhs;
    this.y *= rhs;
  }
  
  public void normalize(){
    float magnitude = this.length();
    this.x /= magnitude;
    this.y /= magnitude;
  }
  
  public Vec2 normalized(){
    float magnitude = this.length();
    return new Vec2(this.x/magnitude, this.y/magnitude);
  }
  
  public float distanceTo(Vec2 rhs){
    return this.minus(rhs).length();
  }
}

public Vec2 interpolate(Vec2 a, Vec2 b, float t){
  return a.plus((b.minus(a)).times(t));
}

public float interpolate(float a, float b, float t) {
   return a + (b - a)*t; 
}

public float dot(Vec2 a, Vec2 b){
  return a.x*b.x + a.y*b.y;
}

public Vec2 projAB(Vec2 a, Vec2 b){
  return b.times(dot(a, b));
}
//Vector Library [2D]
//CSCI 5611 Vector 3 Library [Incomplete]

//Instructions: Add 3D versions of all of the 2D vector functions
//              Vec3 must also support the cross product.
public class Vec3 {
  public float x, y, z;
  
  public Vec3(float x, float y, float z){
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public Vec3(Vec3 copyVec) {
    this.x = copyVec.x;
    this.y = copyVec.y;
    this.z = copyVec.z;
  }

  public Vec3() {
    this.x = 0;
    this.y = 0;
    this.z = 0;
  }

  public Vec3(PVector v) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
  }

  public boolean equals(Vec3 v) {
    return this.x == v.x && this.y == v.y && this.z == v.z;
  }
  
  public String toString(){
    return "(" + x + ", " + y + ", " + z + ")";
  }
  
  public float length(){
    if (x == 0 && y == 0 && z == 0) return 0.f;
    return sqrt(x*x + y*y + z*z);
  }
  
  public Vec3 plus(Vec3 rhs){
    return new Vec3(rhs.x+this.x, rhs.y+this.y, rhs.z+this.z);
  }
  
  public void add(Vec3 rhs){
    this.x += rhs.x;
    this.y += rhs.y;
    this.z += rhs.z;
  }
  
  public Vec3 minus(Vec3 rhs){
    return new Vec3(this.x-rhs.x, this.y-rhs.y, this.z-rhs.z);
  }
  
  public void subtract(Vec3 rhs){
    this.x -= rhs.x;
    this.y -= rhs.y;
    this.z -= rhs.z;
  }
  
  public Vec3 times(float rhs){
    return new Vec3(this.x*rhs, this.y*rhs, this.z*rhs);
  }

  public Vec3 times(Vec3 v) {
    return new Vec3(this.x*v.x, this.y*v.y, this.z*v.z);
  }
  
  public void mul(float rhs){
    this.x *= rhs;
    this.y *= rhs;
    this.z *= rhs;
  }
  
  public void normalize(){
    float magnitude = this.length();
    this.x /= magnitude;
    this.y /= magnitude;
    this.z /= magnitude;
  }
  
  public Vec3 normalized(){
    float magnitude = this.length();
    assert magnitude > 0 : "zero magnitude";
    return new Vec3(this.x/magnitude, this.y/magnitude, this.z/magnitude);
  }
  
  public float distanceTo(Vec3 rhs){
    return this.minus(rhs).length();
  }

  public void clamp(float minLength, float maxLength) {
    float curLength = length();
    if (curLength > maxLength) {
      this.normalize();
      this.mul(maxLength);
    } else if (curLength < minLength) {
      this.normalize();
      this.mul(minLength);
    }
  }
}

public Vec3 interpolate(Vec3 a, Vec3 b, float t){
  return a.plus((b.minus(a)).times(t));
}

public float dot(Vec3 a, Vec3 b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

public Vec3 cross(Vec3 a, Vec3 b){
  return new Vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

public Vec3 projAB(Vec3 a, Vec3 b){
  return b.times(dot(a, b));
}

public Vec3 reflect(Vec3 d, Vec3 n) {
  Vec3 r = d.minus(n.times(dot(d,n.normalized())*2));
  return r;
}

public boolean pointLiesOnTriangle(Vec3 point, Vec3 vert1, Vec3 vert2, Vec3 vert3, Vec3 e1, Vec3 e2) {
  // See if point on a plane is within a triangle
  // Vec3 ep = point.minus(vert1);
  // float d11 = dot(e1,e1);
  // float d12 = dot(e1,e2);
  // float d22 = dot(e2,e2);
  // float dp1 = dot(ep,e1);
  // float dp2 = dot(ep,e2);
  // float D       = d11*d22 - d12*d12;
  // float D_beta  = d22*dp1 - d12*dp2;
  // float D_gamma = d11*dp2 - d12*dp1;
  // float beta = D_beta/D;
  // float gamma = D_gamma/D;
  // float alpha = 1 - beta + gamma;
  // print(alpha);
  // print(" ");
  // print(beta);
  // print(" ");
  // print(gamma);
  // println();
  // return (alpha > 0.0000001 && alpha < 1.0000001);

  // Source inspired by Scratchapixel:
  // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates

  Vec3 surfaceNormal = cross(e1, e2);
  Vec3 C;

  Vec3 edge0 = vert2.minus(vert1);
  Vec3 vp0 = point.minus(vert1);
  C = cross(edge0, vp0);
  if (dot(surfaceNormal, C) < 0) { return false; }

  Vec3 edge1 = vert3.minus(vert2);
  Vec3 vp1 = point.minus(vert2);
  C = cross(edge1, vp1);
  if (dot(surfaceNormal, C) < 0) { return false; }

  Vec3 edge2 = vert1.minus(vert3);
  Vec3 vp2 = point.minus(vert3);
  C = cross(edge2, vp2);
  if (dot(surfaceNormal, C) < 0) { return false; }

  return true;
}
  public void settings() {  size(1280, 960, P3D); }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "SPHFluidProcessing" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
