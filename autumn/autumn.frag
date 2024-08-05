#include "common.frag"
#iChannel0 "assets/organic.jpg"

#define STEP_SIZE        0.5
#define MAXIMUM_STEPS    5000
#define FAR              100.0
#define EPSILON          0.001

vec2 scene(vec3 p) {
  float h = 10.0 * fbm(0.02 * p.xz, 0.5);

  float water = p.y - 6.0;

  float ret = min(water, p.y - h);

  float id = float(ret == water);

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

Material GetMaterial(in Hit hit) {
  Material m;
  if (hit.id == 0.0) {
    m.albedo = texture(iChannel0, 0.1 * hit.position.xz).rgb;
    m.albedo.y += 0.2 * hit.normal.y;
  } else if (hit.id == 1.0) {
    m.albedo = vec3(0.1, 0.4, 0.5);
  }
  return m;
}

vec3 rayMarch(Ray ray) {
  vec3 color = vec3(1.0);
  for (int i = 0; i < 4; i++) {
    Hit hit = trace(ray, 0.1);
    if (hit.id == 1.0) {
      ray.origin = hit.position;
      ray.direction = normalize(refract(ray.direction, hit.normal, 0.5));
      Material m = GetMaterial(hit);
      color *= m.albedo;
    } else if (hit.id == 0.0) {
      Material m = GetMaterial(hit);
      color *= m.albedo;
      break;
    }
  }
  return color;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 uv = (2.0 * fragCoord - iResolution.xy) / iResolution.xy;
  
  uv.y *= iResolution.y / iResolution.x;

  vec3 color = vec3(0.0);

  Ray ray;
  ray.origin = vec3(-15.0, 10.0, 30.0);
  ray.direction = normalize(vec3(uv, -1.0));

  Hit hit = trace(ray, 0.0);

  vec3 LD = normalize(vec3(1.0, 1.0, 0.0));

  color = rayMarch(ray);

  color *= dot(hit.normal, LD);

  color = ACES(color);

  fragColor = vec4(color, 0.0);
}