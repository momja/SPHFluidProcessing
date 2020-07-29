public Vec3 rayPlaneCollision(Ray3 ray, Vec3 planeNormal, Vec3 planeOrigin) {
    float denominator = dot(planeNormal, ray.direction);
    if (abs(denominator) <= 0.000001) {
        // No ray plane intersection exists
        return null;
    }

    float D = dot(planeOrigin, planeNormal);

    float numerator = -(dot(planeNormal, ray.origin) - D);

    float t = numerator/denominator;

    if (t < 0) {
        // Haven't hit yet
        return null;
    }
    
    Vec3 p = ray.origin.plus(ray.direction.times(t));

    return p;
}