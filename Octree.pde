// Max Omdal 2020

class Octree<T extends OctantInsertable> {
    Octant bounds;
    int capacity;
    ArrayList<T> points;
    boolean divided = false;

    Octree<T> q111,q011,q001,q101,q110,q010,q000,q100;

    public Octree(Octant bounds, int capacity) {
        this.bounds = bounds;
        this.capacity = capacity;
        this.points = new ArrayList<T>();
    }

    public void clear() {
        this.points = new ArrayList<T>();
        this.divided = false;

        q111 = q011 = q001 = q101 = q110 = q010 = q000 = q100 = null;
    }

    public void insert(T p) {
        if (!p.inOctant(this.bounds)) {
            return;
        }

        if (!this.divided && this.points.size() < this.capacity) {
            this.points.add(p);
        } else {
            if (!this.divided) {
                subdivide();
            }

            q111.insert(p);
            q011.insert(p);
            q001.insert(p);
            q101.insert(p);
            q110.insert(p);
            q010.insert(p);
            q000.insert(p);
            q100.insert(p);
        }
    }

    public void show(boolean withPoints) {
        if (this.divided) {
            q111.show(withPoints);
            q011.show(withPoints);
            q001.show(withPoints);
            q101.show(withPoints);
            q110.show(withPoints);
            q010.show(withPoints);
            q000.show(withPoints);
            q100.show(withPoints);
        } else {
            pushMatrix();
            strokeWeight(0.8);
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
                pointsInOctants.addAll(q111.rayIntersectsOctants(ray));
                pointsInOctants.addAll(q011.rayIntersectsOctants(ray));
                pointsInOctants.addAll(q001.rayIntersectsOctants(ray));
                pointsInOctants.addAll(q101.rayIntersectsOctants(ray));
                pointsInOctants.addAll(q110.rayIntersectsOctants(ray));
                pointsInOctants.addAll(q010.rayIntersectsOctants(ray));
                pointsInOctants.addAll(q000.rayIntersectsOctants(ray));
                pointsInOctants.addAll(q100.rayIntersectsOctants(ray));
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
                pointsInOctants.addAll(q111.inSameOctant(p));
                pointsInOctants.addAll(q011.inSameOctant(p));
                pointsInOctants.addAll(q001.inSameOctant(p));
                pointsInOctants.addAll(q101.inSameOctant(p));
                pointsInOctants.addAll(q110.inSameOctant(p));
                pointsInOctants.addAll(q010.inSameOctant(p));
                pointsInOctants.addAll(q000.inSameOctant(p));
                pointsInOctants.addAll(q100.inSameOctant(p));
            } else {
                // Just add the points in this oct
            }
                pointsInOctants.addAll(this.points);
        }

        return pointsInOctants;
    }

    private void subdivide() {
        Vec3 origin = bounds.origin;
        Vec3 size = bounds.size;
        Vec3 subdivideSize = size.times(0.5);

        Octant q111_b = new Octant(origin.plus(subdivideSize.times(new Vec3(0.5,0.5,0.5))), subdivideSize);
        Octant q011_b = new Octant(origin.plus(subdivideSize.times(new Vec3(-0.5,0.5,0.5))), subdivideSize);
        Octant q001_b = new Octant(origin.plus(subdivideSize.times(new Vec3(-0.5,-0.5,0.5))), subdivideSize);
        Octant q101_b = new Octant(origin.plus(subdivideSize.times(new Vec3(0.5,-0.5,0.5))), subdivideSize);
        Octant q110_b = new Octant(origin.plus(subdivideSize.times(new Vec3(0.5,0.5,-0.5))), subdivideSize);
        Octant q010_b = new Octant(origin.plus(subdivideSize.times(new Vec3(-0.5,0.5,-0.5))), subdivideSize);
        Octant q000_b = new Octant(origin.plus(subdivideSize.times(new Vec3(-0.5,-0.5,-0.5))), subdivideSize);
        Octant q100_b = new Octant(origin.plus(subdivideSize.times(new Vec3(0.5,-0.5,-0.5))), subdivideSize);
        q111 = new Octree<T>(q111_b, capacity);
        q011 = new Octree<T>(q011_b, capacity);
        q001 = new Octree<T>(q001_b, capacity);
        q101 = new Octree<T>(q101_b, capacity);
        q110 = new Octree<T>(q110_b, capacity);
        q010 = new Octree<T>(q010_b, capacity);
        q000 = new Octree<T>(q000_b, capacity);
        q100 = new Octree<T>(q100_b, capacity);

        // for (T p : points) {
        //     q111.insert(p);
        //     q011.insert(p);
        //     q001.insert(p);
        //     q101.insert(p);
        //     q110.insert(p);
        //     q010.insert(p);
        //     q000.insert(p);
        //     q100.insert(p);
        // }

        // points = new ArrayList<T>();

        this.divided = true;
    }

}