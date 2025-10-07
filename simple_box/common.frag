#ifndef COMMON_FRAG
#define COMMON_FRAG

#define PI       (3.141592653589793)
#define INFINITY (1e23)

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Hit {
  vec3 normal;
  vec3 position;
  bool inside;
  float t;
  int id;
};

#define DIFFUSE  (0)
#define EMISSIVE (1)
#define SPECULAR (2)

struct Material {
  vec3 albedo;
  int type;  
};

// https://iquilezles.org/articles/intersectors/

vec2 sphere_intersect(Ray ray, vec3 center, float radius) {
  vec3 center_to_origin = ray.origin - center;
  float b = dot(center_to_origin, ray.direction);
  float c = dot(center_to_origin, center_to_origin) - radius * radius;
  float D = b * b - c;
  return (D < 0.0) ? vec2(-1.0) : vec2(-b - sqrt(D), -b + sqrt(D));
}

Hit sphere_intersect(Ray ray, vec3 center, float radius, int id) {
  vec2 t = sphere_intersect(ray, center, radius);
  Hit hit; hit.id = -1; hit.t = INFINITY;
  if (t == vec2(-1.0)) return hit;
  hit.inside = t.x < 0.0;
  hit.t = hit.inside ? t.y : t.x;
  hit.position = ray.origin + ray.direction * hit.t;
  hit.normal = normalize(hit.position - center) * (hit.inside ? -1.0 : 1.0);
  hit.id = id;
  return hit;
}

vec4 triangle_intersect(Ray ray, vec3 v0, vec3 v1, vec3 v2) {
  vec3 v0_v1 = v1 - v0;
  vec3 v0_v2 = v2 - v0;
  vec3 v0_origin = ray.origin - v0;
  vec3 normal = cross(v0_v1, v0_v2);
  vec3 q = cross(v0_origin, ray.direction);
  float d = 1.0 / dot(ray.direction, normal);
  float u = d * dot(-q, v0_v2);
  float v = d * dot(+q, v0_v1);
  float t = d * dot(-normal, v0_origin);
  if (u < 0.0 || v < 0.0 || (u + v) > 1.0) t = INFINITY;
  return vec4(t, normalize(normal));
}

void set_box(out vec3 [30] triangles, float x, float y, float z) {
  triangles = vec3 [] (
    vec3(-x, -y, -z), vec3(+x, -y, -z), vec3(+x, +y, -z), vec3(+x, +y, -z), vec3(-x, +y, -z), vec3(-x, -y, -z),
    vec3(-x, -y, +z), vec3(-x, -y, -z), vec3(-x, +y, -z), vec3(-x, +y, -z), vec3(-x, +y, +z), vec3(-x, -y, +z),
    vec3(+x, -y, -z), vec3(+x, -y, +z), vec3(+x, +y, +z), vec3(+x, +y, +z), vec3(+x, +y, -z), vec3(+x, -y, -z),
    vec3(-x, +y, -z), vec3(+x, +y, -z), vec3(+x, +y, +z), vec3(+x, +y, +z), vec3(-x, +y, +z), vec3(-x, +y, -z),
    vec3(-x, -y, +z), vec3(+x, -y, +z), vec3(+x, -y, -z), vec3(+x, -y, -z), vec3(-x, -y, -z), vec3(-x, -y, +z)
  );
}

uint wang_hash(inout uint seed) {
  seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
  seed *= uint(9);
  seed = seed ^ (seed >> 4);
  seed *= uint(0x27d4eb2d);
  seed = seed ^ (seed >> 15);
  return seed;
}

float random(inout uint state) {
  return float(wang_hash(state)) / 4294967296.0;
}

vec3 random_sphere_unit_vector(inout uint seed) {
  float z = random(seed) * 2.0 - 1.0;
  float a = random(seed) * 2.0 * PI;
  float r = sqrt(1.0 - z * z);
  float x = r * cos(a);
  float y = r * sin(a);
  return vec3(x, y, z);
}

vec2 random_disk_unit_vector(inout uint seed) {
  float u = random(seed);
  float f = random(seed) * 2.0 * PI; 
  float r = sqrt(u);
  return r * vec2(sin(f), cos(f));
}

vec3 random_hemisphere_unit_vector(vec3 normal, inout uint seed) {
  vec3 rv = random_sphere_unit_vector(seed);
	vec3 uu = normalize(cross(normal, vec3(0.0, 1.0, 1.0)));
	vec3 vv = cross(uu, normal);
	vec3 rr = vec3(rv.x * uu + rv.y * vv + rv.z * normal);  
  return normalize(rr);
}

#endif // COMMON_FRAG