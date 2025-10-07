#ifndef COMMON_FRAG
#define COMMON_FRAG

#define PI                (3.141592653589793)
#define KM                (1000.0)
#define EARTH_RADIUS      (6360.0 * KM)
#define ATMOSPHERE_RADIUS (6380.0 * KM)
#define HR                (8.0 * KM)
#define HM                (1.2 * KM)
#define G                 (0.76)
#define BM                (vec3(21e-6, 21e-6, 21e-6))
#define BR                (vec3(5.8e-6, 13.5e-6, 33.1e-6))
#define C                 (vec3(0.0, -EARTH_RADIUS, 0.0))

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

mat2 rotate(float a) {
  return mat2(cos(a), sin(a), -sin(a), cos(a));
}

vec2 random2f(vec2 p) {
  float a = dot(p, vec2(13.0, 57.0));
  float b = dot(p, vec2(57.0, 113.0));
  return fract(sin(vec2(a, b)) * 43758.5453);
}

float random1f(vec2 p) {
  float a = dot(p, vec2(13.0, 57.0));
  return fract(sin(a) * 43758.5453);
}

// https://iquilezles.org/articles/distfunctions2d/

float sd_circle(vec2 p, float r) {
  return length(p) - r;
}

// https://iquilezles.org/articles/distfunctions/

float sd_sphere(vec3 p, float radius) {
  return length(p) - radius;
}

// https://iquilezles.org/articles/intersectors/

float sphere_intersect(Ray ray, vec3 center, float R) {
  vec3 center_to_origin = ray.origin - center;
  float b = dot(center_to_origin, ray.direction);
  float c = dot(center_to_origin, center_to_origin) - R * R;
  float D = b * b - c;
  if (D < 0.0) return -1.0;
  float t1 = -b - sqrt(D), t2 = -b + sqrt(D);
  return (t1 >= 0.0) ? t1 : t2;
}

// https://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds/simulating-sky/simulating-colors-of-the-sky.html

float rayleigh_phase_function(float mu) {
  return 3.0 * (1.0 + mu * mu) / (16.0 * PI);
}

float mie_phase_function(float mu, float g) {
  float g2 = g * g;
  float numerator = (1.0 + mu * mu) * (1.0 - g2);
  float denominator = (2.0 + g2) * pow(1.0 + g2 - 2.0 * g * mu, 1.5);
  return 3.0 * numerator / (8.0 * PI * denominator);
}

float noise(vec2 x, sampler2D noise_texture) {
  vec2 i = floor(x);
  vec2 f = fract(x);
  vec2 u = f * f * (3.0 - 2.0 * f);
    
  ivec2 p = ivec2(i);
  float a = texelFetch(noise_texture, (p + ivec2(0, 0)) & 255, 0 ).x;
	float b = texelFetch(noise_texture, (p + ivec2(1, 0)) & 255, 0 ).x;
	float c = texelFetch(noise_texture, (p + ivec2(0, 1)) & 255, 0 ).x;
	float d = texelFetch(noise_texture, (p + ivec2(1, 1)) & 255, 0 ).x;
    
  return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float fbm(vec2 x, float H, int octaves, sampler2D noise_texture) {    
  float g = pow(2.0, -H);
  float f = 1.0;
  float a = 1.0;
  float t = 0.0;
  for (int i = 0; i < octaves; i++) {
    t += a * noise(f * x, noise_texture);
    f *= 2.0;
    a *= g;
  }
  return t;
}

// https://www.shadertoy.com/view/XslGRr

float noise(vec3 p, sampler2D noise_texture) {
  vec3 i = floor(p);
  vec3 f = fract(p); 
  f = f * f * (3.0 - 2.0 * f);
	vec2 uv = (i.xy + vec2(37.0, 239.0) * i.z) + f.xy;
  vec2 rg = textureLod(noise_texture, (uv + 0.5) / 256.0, 0.0).yx;
	return mix(rg.x, rg.y, f.z) * 2.0 - 1.0;
}

float fbm(vec3 x, sampler2D noise_texture, float H, int octaves) {    
  float g = pow(2.0, -H);
  float f = 1.0;
  float a = 1.0;
  float t = 0.0;
  for (int i = 0; i < octaves; i++) {
    t += a * noise(f * x, noise_texture);
    f *= 2.0;
    a *= g;
  }
  return t;
}

#endif // COMMON_FRAG