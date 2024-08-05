#ifndef COMMON_FRAG
#define COMMON_FRAG

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

#define OCTAVES 4

mat2 rotate(float a) {
  return mat2(cos(a), sin(a), -sin(a), cos(a));
}

float hash(in vec2 p) {
  float q = dot(p, vec2(12.9898, 78.233));
  return fract(sin(q) * 43758.5453123);
}

// FBM
float noise(in vec2 x) {
  vec2 i = floor(x);
  vec2 f = fract(x);
  vec2 u = f * f * (3.0 - 2.0 * f);
    
  float a = hash(i + vec2(0.0, 0.0));
  float b = hash(i + vec2(1.0, 0.0));
  float c = hash(i + vec2(0.0, 1.0));
  float d = hash(i + vec2(1.0, 1.0));
    
  return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float fbm(in vec2 x, in float H) {    
  float G = pow(2.0, -H);
  float f = 1.0;
  float a = 1.0;
  float t = 0.0;
  for (int i = 0; i < OCTAVES; i++) {
    t += a * noise(f * x);
    f *= 2.0;
    a *= G;
  }
  return t;
}

vec3 ACES(const vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

#endif // COMMON_FRAG