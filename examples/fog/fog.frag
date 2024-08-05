#include "common.frag"

#define STEP_SIZE        0.5
#define MAXIMUM_STEPS    5000
#define FAR              100.0
#define EPSILON          0.001

vec2 scene(vec3 p) {
  float terrain = p.y;
  return vec2(terrain, 0.0);
}

vec3 getNormal(in vec3 p) {
  const float h = 0.01; 
  const vec2 k = vec2(1.0, -1.0);
  return normalize(k.xyy * scene(p + k.xyy * h).x + 
                   k.yyx * scene(p + k.yyx * h).x + 
                   k.yxy * scene(p + k.yxy * h).x + 
                   k.xxx * scene(p + k.xxx * h).x);
}

Hit trace(in Ray ray, float near, out float fog) {
  Hit hit;
  hit.id = -1.0;
  float t = near;
  for (int i = 0; i < MAXIMUM_STEPS; i++) {
    vec3 p = ray.origin + ray.direction * t;
    vec2 d = scene(p);
    fog += 0.05 * fbm(p - vec3(0.0, 0.0, iTime), 0.5);
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

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 uv = (2.0 * fragCoord - iResolution.xy) / iResolution.xy;
  uv.y *= iResolution.y / iResolution.x;

  Ray ray;
  ray.origin = vec3(0.0, 4.0, 9.0);
  ray.direction = normalize(vec3(uv, -1.0));

  vec3 color = vec3(0.0);

  float fog = 0.0;
  Hit hit = trace(ray, 0.0, fog);

  vec3 LD = normalize(vec3(1.0, 1.0, 0.0));

  if (hit.id != -1.0) {
    color = getChess(hit.position.xz, vec3(0.0, 1.0, 0.0), vec3(0.0, 0.0, 1.0));
    color *= dot(hit.normal, LD); 
  }

  fragColor = vec4(color, 0.0);
}