#ifndef COMMON_FRAG
#define COMMON_FRAG

#define PI (3.14159265358979323)

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

vec2 MIN(vec2 p, vec2 q) {
  return (p.x < q.x) ? p : q;
}

float smin_root(float a, float b, float k) {
  float x = b - a;
  return 0.5 * (a + b - sqrt(x * x + 4.0 * k * k));
}

float sd_capsule(vec3 p, float r, float h) {
  p.y -= clamp(p.y, 0.0, h);
  return length(p) - r;
}

float repeat_angle(vec2 p, float n) {
  float sp = 2.0 * PI / n;
  float an = atan(p.y, p.x);
  float id = floor(an / sp);
  return sp * id;
}

mat2 rotate(float a) {
  return mat2(cos(a), sin(a), -sin(a), cos(a));
}

float random1f(vec2 p) {
  return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453123);
}

float random1f(vec3 p) {
  return fract(sin(dot(p, vec3(12.9898, 78.233, 37.119))) * 43758.5453123);
}

vec2 random2f(vec2 p)  {
  vec2 k = vec2(0.3183099, 0.3678794);
  float n = dot(vec2(111.0, 113.0), p);
  return fract(n * fract(k * n));
}

float noise(vec3 x) {
  vec3 i = floor(x);
  vec3 f = fract(x);

  f = f * f * (3.0 - 2.0 * f);

  float A = random1f(i + vec3(0.0, 0.0, 0.0));
  float B = random1f(i + vec3(1.0, 0.0, 0.0));
  float C = random1f(i + vec3(0.0, 1.0, 0.0));
  float D = random1f(i + vec3(1.0, 1.0, 0.0));
  float E = random1f(i + vec3(0.0, 0.0, 1.0));
  float F = random1f(i + vec3(1.0, 0.0, 1.0));
  float G = random1f(i + vec3(0.0, 1.0, 1.0));
  float H = random1f(i + vec3(1.0, 1.0, 1.0));

  return mix(mix(mix(A, B, f.x), mix(C, D, f.x), f.y),
             mix(mix(E, F, f.x), mix(G, H, f.x), f.y), f.z);
}

float fbm(vec3 x, float H, int octaves) {    
  float G = pow(2.0, -H);
  float f = 1.0;
  float a = 1.0;
  float t = 0.0;
  for (int i = 0; i < octaves; i++) {
    t += a * noise(f * x);
    f *= 2.0;
    a *= G;
  }
  return t;
}

#endif // COMMON_FRAG