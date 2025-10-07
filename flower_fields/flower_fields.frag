/*
TODO:
  1. Clouds on ShaderToy are bad.
*/

#include "common.frag"
#iChannel0 "assets/noise.png"

#define EPSILON          (0.001)
#define TERRAIN_ID       (0.0)
#define GRASS_ID         (1.0)
#define AA               (1)
#define GRASS_FREQUENCY  (0.5)


float get_height(vec2 p) {
  vec2 m1 = vec2(+600.0, -1200.0);
  vec2 m2 = vec2(-400.0, -400.0);
  float h1 = 100.0 * smoothstep(800.0, 100.0, length(p - m1));
  float h2 = 50.0 * smoothstep(500.0, 100.0, length(p - m2));
  return h1 + h2;
}

vec2 sd_terrain(vec3 p) {
  float h = get_height(p.xz);
  return vec2(p.y - h, TERRAIN_ID);
}

// https://www.shadertoy.com/view/dd2cWh

float sd_grass_blade(vec2 p, float r) {
  float bias = 0.1;
  float ret = sd_circle(p - vec2(r, 0.0), r);
  float sub = sd_circle(p - vec2(r, bias), r - bias);
  ret = max(ret, -sub);
  ret = max(ret, -p.y);
  ret = max(ret, p.x - 1.0);
  return ret;
}

float sd_grass_blade(vec3 p, float thickness) {
  float d = max(0.0, sd_grass_blade(p.xy, 2.0));
  float f = sqrt(d * d + p.z * p.z) - thickness;
  return f;
}

vec2 sd_grass(vec3 p, float n) {
  p.y -= get_height(p.xz);
  vec2 id = round(p.xz / n);
  vec2 offset = sign(p.xz - n * id);
  float d = 1e20;
  const int N = 2;
  for (int j = 0; j < N; j++ ) {
    for (int i = 0; i < N; i++) {
      vec2 rid = id + vec2(i, j) * offset;      
      vec3 r = p; 
      r.xz -= n * rid;
      r.xz *= rotate(2.0 * PI * random1f(rid));
      d = min(d, sd_grass_blade(r, 0.01));
    }
  }
  return vec2(d, GRASS_ID);
}

vec2 map(vec3 p) {
  vec2 ret = sd_terrain(p);
  vec2 grass = sd_grass(p, GRASS_FREQUENCY);
  ret = MIN(ret, grass);
  return ret;
}

// https://iquilezles.org/articles/normalsSDF/

vec3 get_normal(vec3 p) {
  const float h = 0.01; 
  const vec2 k = vec2(1.0, -1.0);
  return normalize(k.xyy * map(p + k.xyy * h).x + 
                   k.yyx * map(p + k.yyx * h).x + 
                   k.yxy * map(p + k.yxy * h).x + 
                   k.xxx * map(p + k.xxx * h).x);
}

Hit march(Ray ray, float near, float far, float step_size, int step_count) {
  Hit hit; 
  hit.id = -1.0;
  float t = near;
  for (int i = 0; i < step_count && t < far; i++) {
    vec3 p = ray.origin + ray.direction * t;
    vec2 d = map(p);
    if (d.x < EPSILON) {
      hit.t = t;
      hit.id = d.y; 
      hit.position = p; 
      hit.normal = get_normal(hit.position);
      break;
    }
    t += step_size * d.x;
  }
  return hit;
}

float get_clouds(vec3 p) {
  float cloud = fbm(0.002 * p, iChannel0, 0.5, 3);
  cloud = smoothstep(0.1, 0.8, cloud);
  cloud *= cloud * 10.0;
  return cloud;
}

vec2 densities(vec3 p) {
  vec2 ret;
  float h = length(p - C) - EARTH_RADIUS;
  float cloud = 0.0;
  if (4000.0 < h && h < 10000.0) {
    cloud = get_clouds(p + vec3(0.0, 0.0, -100.0 * iTime));
  }
  ret.x = exp(-h / HR);
  ret.y = exp(-h / HM) + cloud;
  return ret;
}

vec3 compute_sun_transmittance(vec3 p, vec3 sun_direction, float sun_step_count, vec2 depth) {
  Ray sun_ray = Ray(p, sun_direction);
  vec3 A = vec3(0.0, 0.0, 0.0);
  float Ls = sphere_intersect(sun_ray, C, ATMOSPHERE_RADIUS);
  if (Ls > 0.0) {
    float sun_step_size = Ls / sun_step_count;
    vec2 sun_depth = vec2(0.0, 0.0);
    for (float j = 0.0; j < sun_step_count; j++) {
      vec3 p_sun = p + sun_direction * sun_step_size * j;
      vec2 d_sun = densities(p_sun);
      sun_depth += d_sun * sun_step_size;
    }
    A = exp(-(BR * (sun_depth.x + depth.x) + BM * (sun_depth.y + depth.y)));
  }
  return A;
}

vec3 scatter(Ray ray, vec3 sun_direction, float step_count, float sun_step_count) {
  float L = sphere_intersect(ray, C, ATMOSPHERE_RADIUS);
  float mu = dot(sun_direction, ray.direction);
  float phase_r = rayleigh_phase_function(mu);
  float phase_m = mie_phase_function(mu, G);
  float step_size = L / step_count;
  vec3 R = vec3(0.0), M = vec3(0.0);
  vec2 depth = vec2(0.0, 0.0);
  for (float i = 0.0; i < step_count; i++) {
    vec3 p = ray.origin + ray.direction * step_size * i;
    vec2 d = densities(p);
    d *= step_size; depth += d;
    vec3 A = compute_sun_transmittance(p, sun_direction, sun_step_count, depth);
    if (A != vec3(0.0)) {
      R += A * d.x;
      M += A * d.y;
    } else {
      return vec3(0.0, 0.0, 0.0);
    }
  }
  return 6.0 * (R * BR * phase_r + M * BM * phase_m);
}

vec3 apply_fog(vec3 color, float b, float t, vec3 ray_direction, vec3 sun_direction) {
  float fog_amount = 1.0 - exp(-t * b);
  float sun_amount = max(dot(ray_direction, sun_direction), 0.0);
  vec3 fog_color = vec3(0.5, 0.6, 0.7);
  vec3 sun_color = vec3(1.0, 0.6, 0.6);
  vec3 ret_color = mix(fog_color, sun_color, pow(sun_amount, 8.0));
  return mix(color, ret_color, fog_amount);
}

vec3 get_sun_direction() {
  vec3 sun_direction = vec3(0.8, 0.35, -1.0);
  return normalize(sun_direction);
}

// https://iquilezles.org/articles/rmshadows/

float soft_shadow(Ray ray, float near, float far, int steps) {
  float ret = 1.0;
  float t = near;
  for (int i = 0; i < steps && t < far; i++) {
    vec3 p = ray.origin + ray.direction * t;
    float d = map(p).x;
    ret = min(ret, 10.0 * d / t);
    if (ret < EPSILON) break;
    t += d;
  }
  ret = clamp(ret, 0.0, 1.0);
  return ret * ret * (3.0 - 2.0 * ret);
}

vec3 get_material(Hit hit) {
  vec3 color = vec3(0.0);
  if (hit.id == TERRAIN_ID) {
    color = vec3(0.05, 0.05, 0.05);
  } else if (hit.id == GRASS_ID) {
    color = vec3(0.1, 0.5, 0.1);
    color -= 0.25 * smoothstep(2.0, 0.0, hit.position.y);
  }
  return color;
}

float ambient_occlusion(vec3 p, vec3 n) {
  float ao_step = 0.1;
  float ao = 0.0;
  for (float i = 1.0; i < 8.0; i += 1.0) {
    float d = ao_step * i;
    ao += max(0.0, (d - map(p + n * d).x) / d);
  }
  return clamp(1.0 - ao * 0.1, 0.0, 1.0);
}

// https://iquilezles.org/articles/outdoorslighting/

vec3 get_light(Hit hit) {
  vec3 ambient = get_material(hit);
  vec3 sun_direction = get_sun_direction();
  Ray shadow_ray = Ray(hit.position, sun_direction);
  float shadow = soft_shadow(shadow_ray, 0.1, 100.0, 32);
  float ao = ambient_occlusion(hit.position, hit.normal);
  float sun = clamp(dot(sun_direction, hit.normal), 0.0, 1.0);
  float sky = clamp(0.5 + 0.5 * hit.normal.y, 0.0, 1.0);
  float indirect = clamp(dot(hit.normal, normalize(sun_direction * vec3(-1.0, 0.0, -1.0))), 0.0, 1.0);
  vec3 ret = sun * vec3(1.64, 1.0, 0.99) * pow(vec3(shadow), vec3(1.0,1.2,1.5));
  ret += sky * vec3(0.16, 0.20, 0.28) * ao;
  ret += indirect * vec3(0.40, 0.28, 0.20) * ao;
  return ret * ambient;
}

vec3 render(Ray ray) {
  Hit hit = march(ray, 0.0, 5000.0, 0.5, 8000);
  if (hit.id != -1.0) {
    vec3 color = get_light(hit);
    color = apply_fog(color, 0.0002, hit.t, ray.direction, get_sun_direction());
    return color;
  }

  return scatter(ray, get_sun_direction(), 64.0, 8.0);
}

void mainImage(out vec4 out_color, in vec2 in_position) {
  int Z = min(iFrame, 0);
  vec3 total = vec3(0.0);

  for (int m = Z; m < AA; m++) {
    for (int n = Z; n < AA; n++) {
      vec2 offset = vec2(float(m), float(n)) / float(AA) - 0.5;
      vec2 uv = (2.0 * (in_position + offset) - iResolution.xy) / iResolution.y;

      Ray ray;
      ray.origin = vec3(0.0, 8.0, 8.0);
      ray.direction = normalize(vec3(uv, -1.0)); 

      vec3 color = render(ray);

      // gamma correction
      color = pow(color, vec3(0.4545));
      total += color;
    }
  }

  total /= float(AA * AA);
  out_color = vec4(total, 1.0);
}
