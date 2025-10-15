#ifndef COMMON_FRAG
#define COMMON_FRAG

#define PI (3.141592653589793)

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

vec2 MIN(vec2 p, vec2 q) { return (p.x < q.x) ? p : q; }
vec2 MAX(vec2 p, vec2 q) { return (p.x > q.x) ? p : q; }

vec2 SUB(vec2 p, vec2 q) { return MAX(vec2(-p.x, p.y), q); }

mat2 rotate(float a) {
  return mat2(cos(a), sin(a), -sin(a), cos(a));
}

float repeat_angle(vec2 p, float n) {
  float sp = 2.0 * PI / n;
  float an = atan(p.y, p.x);
  float id = round(an / sp);
  return sp * id;
}

float smin_quadratic(float a, float b, float k) {
  float s = 4.0 * k;
  float h = max(s - abs(a - b), 0.0) / s;
  return min(a, b) - h * h * k;
}

float random1f(vec2 p) {
  float a = dot(p, vec2(13.0, 57.0));
  return fract(sin(a) * 43758.5453);
}

// SDF

float sd_capsule(vec3 p, float r, float h) {
  p.y -= clamp(p.y, 0.0, h);
  return length(p) - r;
}

float sd_triangle_isosceles(vec2 p, vec2 size) {
  p.x = abs(p.x);
	vec2 a = p - size * clamp(dot(p, size) / dot(size, size), 0.0, 1.0);
  vec2 b = p - size * vec2(clamp(p.x / size.x, 0.0, 1.0), 1.0);
  float k = sign(size.y);
  float d = min(dot(a, a), dot(b, b));
  float s = max(k * (p.x * size.y - p.y * size.x), k * (p.y - size.y));
	return sqrt(d) * sign(s);
}

float sd_triangular_prism(vec3 p, vec2 size, float width) {
  float d = sd_triangle_isosceles(p.xy, size * vec2(1.0, -1.0));
  vec2 w = vec2(d, abs(p.z) - width);
  return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

float sd_box(vec3 p, vec3 size) {
  vec3 q = abs(p) - size;
  return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float sd_torus(vec3 p, float R, float r) {
  vec2 q = vec2(length(p.xz) - R, p.y);
  return length(q) - r;
}

float sd_circle(vec2 p, float r) {
  return length(p) - r;
}

float sd_octahedron(vec3 p, float s) {
  p = abs(p);
  return (p.x + p.y + p.z - s) * sqrt(3.0) / 3.0;
}

// TEXTURE

vec4 texture_map(sampler2D texture_sampler, vec3 p, vec3 n, float k) {
  vec4 x = texture(texture_sampler, p.yz);
  vec4 y = texture(texture_sampler, p.zx);
  vec4 z = texture(texture_sampler, p.xy);
  vec3 w = pow(abs(n), vec3(k));  
  return (x * w.x + y * w.y + z * w.z) / (w.x + w.y + w.z);
}

mat3 look_at(vec3 origin, vec3 target) {
  vec3 wu = vec3(0.0, 1.0, 0.0);
  vec3 f = normalize(origin - target);
  vec3 r = normalize(cross(wu, f));
  vec3 u = normalize(cross(f, r));
  return mat3(r, u, f);
}

// https://www.shadertoy.com/view/XslGRr

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

float fbm(vec2 x, sampler2D noise_texture, float H, int octaves) {    
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