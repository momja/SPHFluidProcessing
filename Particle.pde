// Max Omdal 2020

class Particle implements OctantInsertable {
    float radius = 0.1;
    Float smoothingRadius;
    Vec3 pos;
    Vec3 oldPos;
    Vec3 vel;
    float density;
    float densityN;
    float pressure;
    float pressureN;


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
        noStroke();
        fill(0,50,180);
        translate(this.pos.x, this.pos.y, this.pos.z);
        sphere(this.radius);
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
