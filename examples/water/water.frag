#include "common.frag"
#iChannel0 "assets/pebbles.png"

#define STEP_SIZE        0.5
#define MAXIMUM_STEPS    5000
#define FAR              100.0
#define EPSILON          0.001

bool is_water = true;

vec2 scene(vec3 p) {
  float h = texture(iChannel0, 0.1 * p.xz).r;
  float d = p.y - 1.5 * h;
  float ret = d;

  if (is_water) {
    float water = p.y - 2.0 + 0.1 * sin(5.0 * h + iTime);

    ret = min(water, ret);
  }

  float id = float(ret == d);
  return vec2(ret, id);
}

vec3 getNormal(in vec3 p) {
  const float h = 0.01; 
  const vec2 k = vec2(1.0, -1.0);
  return normalize(k.xyy * scene(p + k.xyy * h).x + 
                   k.yyx * scene(p + k.yyx * h).x + 
                   k.yxy * scene(p + k.yxy * h).x + 
                   k.xxx * scene(p + k.xxx * h).x);
}

Hit trace(in Ray ray, float near) {
  Hit hit;
  hit.id = -1.0;
  float t = near;
  for (int i = 0; i < MAXIMUM_STEPS; i++) {
    vec3 p = ray.origin + ray.direction * t;
    vec2 d = scene(p);
    if (d.x < EPSILON) {
      hit.id = d.y;
      hit.position = p;
      hit.normal = getNormal(hit.position);
      break;
    }
    t += STEP_SIZE * d.x;
    if (t > FAR) { break; }
  }
  return hit;
}


vec3 rayMarch(Ray ray) {
  vec3 color = vec3(0.0);
  Hit hit = trace(ray, 0.0);
  vec3 LD = normalize(vec3(1.0, 1.0, 0.0));
  if (hit.id != -1.0) {
    if (hit.id == 0.0) {
      color = vec3(0.1, 0.5, 0.5);
      is_water = false;
      Ray refracted;
      refracted.origin = hit.position;
      refracted.direction = refract(ray.direction, hit.normal, AIR_IOR / WATER_IOR);
      Hit water_hit = trace(refracted, 0.1);
      if (water_hit.id == 1.0) {
        color *= vec3(0.5, 0.3, 0.2);
        color *= dot(LD, water_hit.normal);
      }
    }
  }
  return color;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 uv = (2.0 * fragCoord - iResolution.xy) / iResolution.xy;
  uv.y *= iResolution.y / iResolution.x;

  vec3 color = vec3(0.0);

  Ray ray;
  ray.origin = vec3(0.0, 5.0, 10.0);
  ray.direction = normalize(vec3(uv, -1.0));

  color = rayMarch(ray);


  fragColor = vec4(color, 0.0);
}