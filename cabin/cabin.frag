#include "common.frag"
#iChannel0 "assets/tree.jpg"

#define STEP_SIZE        0.5
#define MAXIMUM_STEPS    5000
#define FAR              100.0
#define EPSILON          0.001

#define LOG_R            0.5            
#define ROOM_SIZE        5.0

vec2 home(vec3 p) {
  vec3 q0 = p; q0.xz = abs(q0.xz); vec3 q1 = q0 - vec3(0.0, LOG_R, 0.0);
  q0.y -= clamp(2.0 * LOG_R * round(0.5 * q0.y / LOG_R),  0.0, ROOM_SIZE);
  q1.y -= clamp(2.0 * LOG_R * round(0.5 * q1.y / LOG_R), -1.0, ROOM_SIZE);
  float x_wall = sdCylinder(q0.xzy - vec3(ROOM_SIZE, 0.0, 0.0), vec2(LOG_R, ROOM_SIZE));
  float z_wall = sdCylinder(q1.yxz - vec3(0.0, 0.0, ROOM_SIZE), vec2(LOG_R, ROOM_SIZE));
  q0 = p; q0.z -= clamp(0.5 * round(2.0 * q0.z), -ROOM_SIZE, ROOM_SIZE);
  q0.y = abs(q0.y - 0.5 * ROOM_SIZE);
  float floor = sdBox(q0 - vec3(0.0, ROOM_SIZE / 2.0, 0.0), vec3(ROOM_SIZE, 0.1, 0.235));
  float window = sdBox(p - vec3(0.0, ROOM_SIZE / 2.0, -ROOM_SIZE), vec3(1.0, 1.0, LOG_R + 1.0));
  float ret = max(-window, min(floor, min(x_wall, z_wall)));
  float id = FLOOR_ID * float(ret == floor) + WALL_ID * float(ret != floor);
  return vec2(ret, id);
}

vec2 scene(vec3 p) {
  vec2 ret = home(p);

  return ret;
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

Material getMaterial(Hit hit) {
  Material m;
  if (hit.id == WALL_ID) {
    vec2 uv = vec2(hit.position.x + hit.position.z, atan(hit.normal.x + hit.normal.z, hit.normal.y));
    m.albedo = texture(iChannel0, 0.2 * uv).rgb;
  } else if (hit.id == FLOOR_ID) {
    m.albedo = getTexture(iChannel0, hit.position, hit.normal);
  }
  return m;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 uv = (2.0 * fragCoord - iResolution.xy) / iResolution.xy;
  uv.y *= iResolution.y / iResolution.x;

  vec3 color = vec3(0.0);

  Ray ray;
  ray.origin = vec3(0.0, 2.0, 0.0);
  ray.direction = normalize(vec3(uv, -1.0));

  ray.direction.xz *= rotate(iTime);

  Hit hit = trace(ray, 0.0);
  
  vec3 LD = normalize(vec3(1.0, 2.0, 1.0));

  if (hit.id != -1.0) {
    color = getMaterial(hit).albedo;
    color *= 0.5 * dot(hit.normal, LD);
  }

  color = ACES(color);
  color = pow(color, vec3(1.0 / 2.2));

  fragColor = vec4(color, 0.0);
}