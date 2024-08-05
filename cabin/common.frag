#ifndef COMMON_FRAG
#define COMMON_FRAG

#define WALL_ID 0.0
#define FLOOR_ID 1.0

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Hit {
  float t;
  float id;
  vec3 position;
  vec3 normal;
};

struct Material {
  vec3 albedo;
};

mat2 rotate(float a) {
  return mat2(cos(a), sin(a), -sin(a), cos(a));
}

// SDF
float sdBox(vec3 p, vec3 b) {
  vec3 q = abs(p) - b;
  return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float sdCylinder(in vec3 p, in vec2 rh) {
  vec2 d = vec2(length(p.xz), abs(p.y)) - rh;
  return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}

vec3 ACES(const vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

vec3 getTexture(sampler2D s, vec3 p, vec3 n) {
  return texture(s, p.xy).xyz * abs(n.z) +
         texture(s, p.xz).xyz * abs(n.y) +
         texture(s, p.zy).xyz * abs(n.x); 
}

#endif // COMMON_FRAG