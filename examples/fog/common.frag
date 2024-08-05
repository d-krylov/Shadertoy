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

// https://www.shadertoy.com/view/XsXfRH
float hash(vec3 p) {
  p = 50.0 * fract(p * 0.3183099 + vec3(0.71, 0.113, 0.419));
  return -1.0 + 2.0 * fract(p.x * p.y * p.z * (p.x + p.y + p.z));
}

mat2 rotate(float a) {
  return mat2(cos(a), sin(a), -sin(a), cos(a));
}

float noise(in vec3 x) {
  vec3 p = floor(x);
  vec3 w = fract(x);

  vec3 u = w * w * w * (w * (w * 6.0 - 15.0) + 10.0);

  float a = hash(p + vec3(0.0, 0.0, 0.0));
  float b = hash(p + vec3(1.0, 0.0, 0.0));
  float c = hash(p + vec3(0.0, 1.0, 0.0));
  float d = hash(p + vec3(1.0, 1.0, 0.0));
  float e = hash(p + vec3(0.0, 0.0, 1.0));
  float f = hash(p + vec3(1.0, 0.0, 1.0));
  float g = hash(p + vec3(0.0, 1.0, 1.0));
  float h = hash(p + vec3(1.0, 1.0, 1.0));

  float ab = mix(a, b, u.x);
  float cd = mix(c, d, u.x);
  float ef = mix(e, f, u.x);
  float gh = mix(g, h, u.x);

  return mix(mix(ab, cd, u.y), mix(ef, gh, u.y), u.z);
}

#define OCTAVES 3

float fbm(in vec3 x, in float H) {    
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

float BeerLambert(float absorption, float distance) {
  return exp(-absorption * distance);
}

#endif // COMMON_FRAG