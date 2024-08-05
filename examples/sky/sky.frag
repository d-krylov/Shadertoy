#include "common.frag"
#iChannel0 "assets/pebbles.png"

#define STEP_SIZE        0.5
#define MAXIMUM_STEPS    5000
#define FAR              100.0
#define EPSILON          0.001

vec2 scene(vec3 p) {
  float plane = p.y;
  return vec2(plane, PLANE_ID);
}

vec3 sky(in Ray ray) {
  float R0 = 1000.0;
  float R1 = 1500.0;
  vec2 t0; vec2 t1;
  bool b0 = sphereIntersect(ray, vec4(0.0, -500.0, 0.0, R0), t0);
  bool b1 = sphereIntersect(ray, vec4(0.0, -500.0, 0.0, R1), t1);
  float ray_begin = max(t0.x, t0.y);
  float ray_end   = max(t1.x, t1.y);
  ray.origin += ray.direction * ray_begin;
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

vec3 getChess(vec2 p, vec3 x, vec3 y) {
  vec2 q = floor(p);
  float m = mod(q.x + q.y, 2.0);
  return mix(x, y, m);
}

vec3 getMaterial(Hit hit) {
  if (hit.id == PLANE_ID) {
    return getChess(hit.position.xz, vec3(0.0, 1.0, 0.0), vec3(0.0, 0.0, 1.0));
  } else if (hit.id == SKY_ID) {
    return vec3(0.0, 0.0, 1.0);
  }
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 uv = (2.0 * fragCoord - iResolution.xy) / iResolution.xy;
  uv.y *= iResolution.y / iResolution.x;

  vec3 color = vec3(0.0);

  Ray ray;
  ray.origin = vec3(0.0, 2.0, 0.0);
  ray.direction = normalize(vec3(uv, -1.0));

  Hit hit = trace(ray, 0.0);

  vec3 LD = normalize(vec3(1.0, 2.0, -1.0));

  if (hit.id != -1.0) {
    color = getMaterial(hit);
    color *= dot(LD, hit.normal);
  } else {
    color = sky(ray);
  }



  fragColor = vec4(color, 0.0);
}