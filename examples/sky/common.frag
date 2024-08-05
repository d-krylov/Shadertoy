#ifndef COMMON_FRAG
#define COMMON_FRAG

#define PLANE_ID  0.0
#define SKY_ID    1.0

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Hit {
  vec3 position;
  vec3 normal;
  float t;
  float id;
};

vec2 MIN(vec2 a, vec2 b) {
  return (a.x < b.x) ? a : b;
}

// Intersection

bool sphereIntersect(in Ray ray, in vec4 sphere, out vec2 t) {
  vec3 origin_to_center = ray.origin - sphere.xyz;
  float b = dot(origin_to_center, ray.direction);
  float c = dot(origin_to_center, origin_to_center) - sphere.w * sphere.w;
  float D = b * b - c;
  if (D < 0.0) return false;
  t = vec2(-b - sqrt(D), -b + sqrt(D));
  return true;
}

#endif // COMMON_FRAG