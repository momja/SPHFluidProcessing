// Max Omdal 2020
import java.util.HashSet;


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
        Vec3 subdivideSize = size.times(0.5);
        this.children = new ArrayList<Octree<T>>();

        Vec3[] octantSections = new Vec3[]{
            new Vec3(0.5,0.5,0.5),
            new Vec3(-0.5,0.5,0.5),
            new Vec3(-0.5,-0.5,0.5),
            new Vec3(0.5,-0.5,0.5),
            new Vec3(0.5,0.5,-0.5),
            new Vec3(-0.5,0.5,-0.5),
            new Vec3(-0.5,-0.5,-0.5),
            new Vec3(0.5,-0.5,-0.5)
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