#include "common.frag"
#iChannel0 "assets/ground.jpg"
#iChannel1 "assets/brick.jpg"
#iChannel2 "assets/noise.png"

const float EPSILON          = 0.001;
const float TERRAIN_ID       = 0.0;
const float TREE_ID          = 1.0;
const float TEMPLE_ID        = 2.0;
const float TEMPLE_INSIDE_ID = 3.0;
const float TEMPLE_WINDOW_ID = 4.0;
const float GRASS_ID         = 5.0;

float get_height(vec2 p) {
  float d = smoothstep(0.0, 150.0, length(p));
  float terrain = fbm(0.2 * p, iChannel2, 1.0, 4);
  float h = terrain * d;
  return h;
}

vec2 sd_terrain(vec3 p) {
  float h = get_height(p.xz);
  return vec2(p.y - h, TERRAIN_ID);
}

vec2 sd_tree(vec3 p, float deep, float height, float radius, float number) {
  float tree = sd_capsule(p, radius, height);
  for (float i = 0.0; i < deep; i++) {
    p.xz *= rotate(PI / 6.0);
    float repeat = repeat_angle(p.xz, number);
    p.y -= height;
    p.xz *= rotate(repeat);
    p.xy *= rotate(PI / 4.0); 
    p = p.yxz;
    radius /= sqrt(number);
    height /= 1.5;
    float branch = sd_capsule(p, radius, height);
    tree = smin_quadratic(tree, branch, 0.1);
  }
  return vec2(tree, TREE_ID);
}

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

vec2 sd_temple(vec3 p) {
  vec3 b1_size = vec3(35.0, 20.0, 30.0);
  vec3 b2_size = vec3(0.35, 1.5, 1.1) * b1_size;
  vec3 b3_size = vec3(1.0, 1.5, 0.0) * b2_size + vec3(0.0, 0.0, b2_size.x);
  vec3 b4_size = vec3(1.0, 0.4, 0.0) * b2_size + vec3(0.0, 0.0, 3.0);
  vec3 w1_size = vec3(3.0, 6.0, 5.0);
  vec3 q = p; q.xz *= rotate(PI / 4.0); // For rotating octahedron on tower
  vec3 r = p; r.x = abs(r.x);
  float r1_size = 10.0, r2_size = 20.0, r3_size = 15.0;
  float offset1 = b2_size.z - b1_size.z;
  float offset2 = b2_size.z + b4_size.z;
  vec2 w1 = vec2(sd_box(r - vec3(0.7 * b1_size.x, b1_size.yz), w1_size), TEMPLE_WINDOW_ID);
  vec2 i1 = vec2(sd_box(p - vec3(0.0, b1_size.y, 0.0), b1_size), TEMPLE_INSIDE_ID);
  float t1 = sd_torus(p.xzy - vec3(0.0, 1.8 * b2_size.y, b2_size.z + offset1).xzy, 6.0, 1.0);
  float b1 = sd_box(p - vec3(0.0, b1_size.y, 0.0), b1_size);
  float b2 = sd_box(p - vec3(0.0, b2_size.y, offset1), b2_size);
  float b3 = sd_box(p - vec3(0.0, b3_size.y, 0.0), b3_size);
  float b4 = sd_box(p - vec3(0.0, b4_size.y, offset2), b4_size);
  float r1 = sd_triangular_prism(p - vec3(0.0, 2.0 * b2_size.y + r1_size, offset1), vec2(b2_size.x, r1_size), b2_size.z);
  float r2 = sd_triangular_prism(p - vec3(0.0, 2.0 * b1_size.y + r2_size, 0.0), vec2(b1_size.x, r2_size), b1_size.z);
  float r3 = sd_triangular_prism(p - vec3(0.0, 2.0 * b4_size.y + r3_size, offset2), vec2(b4_size.x, r3_size), b4_size.z);
  float r4 = sd_octahedron(q - vec3(0.0, 2.0 * b3_size.y, 0.0), sqrt(2.0) * b3_size.z);
  vec2 temple = vec2(min(t1, min(r3, min(b4, min(r4, min(b3, min(r2, min(r1, min(b1, b2)))))))), TEMPLE_ID);
  vec2 ret = SUB(w1, temple);
  ret = SUB(i1, ret);
  return ret;
}

vec3 pretransform(vec3 p, float kx) {
  p.y -= 0.5 * sin(p.x) + 0.4 * cos(p.z);
  p.x -= kx * sin(p.y);
  return p;
}

vec2 sd_trees(vec3 p) {
  vec2 ret = vec2(1e9, -1.0);
  vec3 p1 = vec3(70.0, 0.0, 0.0);
  vec3 p2 = vec3(-60.0, 0.0, 50.0);
  vec3 p3 = vec3(50.0, 0.0, 60.0);
  vec3 p4 = vec3(20.0, 0.0, 65.0);
  if (length(p.xz - p1.xz) < 40.0) {
    vec3 q = pretransform(p, 0.2); 
    ret = MIN(ret, sd_tree(q - p1, 4.0, 20.0, 1.5, 4.0));
  } else if (length(p.xz - p2.xz) < 40.0) {
    vec3 q = pretransform(p, 0.0); 
    ret = MIN(ret, sd_tree(q - p2, 4.0, 25.0, 1.5, 5.0));
  } else if (length(p.xz - p3.xz) < 20.0) {
    vec3 q = pretransform(p, 0.4);    
    ret = MIN(ret, sd_tree(q - p3, 5.0, 8.0, 0.5, 4.0));
  } else if (length(p.xz - p4.xz) < 20.0) {
    vec3 q = pretransform(p, 0.4);    
    ret = MIN(ret, sd_tree(q - p4, 6.0, 9.0, 0.5, 3.0));
  }
  return ret;
}

vec2 map(vec3 p) {
  //p.xz *= rotate(iTime);
  vec2 ret = sd_temple(p);
  vec2 terrain = sd_terrain(p);
  vec2 grass = sd_grass(p, 0.5);
  ret = MIN(ret, terrain);
  ret = MIN(ret, sd_trees(p));
  ret = MIN(ret, grass);
  return ret;
}

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

const vec3 MOON_DIRECTION = normalize(vec3(1.0, 1.0, -1.0));
const vec3 SKY_COLOR = pow(vec3(0.16, 0.22, 0.33), vec3(2.2));
const vec3 MOON_COLOR = pow(vec3(0.66, 0.76, 0.88), vec3(2.2));
const vec3 EVIL_POSITION = vec3(0.0, 10.0, 0.0);
const vec3 EVIL_COLOR = 5.0 * vec3(1.0, 0.15, 0.0);

vec3 get_sky_color(Ray ray) {
  vec3 clouds = vec3(0.0);
  float s = 0.25;
  for (int i = 0; i < 2; ++i) {
    vec2 p = ray.direction.xz / (ray.direction.y) - s * iTime;
    clouds += fbm(p, iChannel2, 1.0, 3);
    s *= 1.35;
  }
  vec3 color = SKY_COLOR + 0.05 * clouds;
  vec2 moon_position = ray.direction.xy / ray.direction.z + vec2(0.0, 0.4);
  color = mix(color, 0.2 * MOON_COLOR + clouds, smoothstep(0.2, 0.1, length(moon_position)));
  return color;
}

vec3 get_fog_color(Ray ray) {
  vec3 color = SKY_COLOR;
  color += MOON_COLOR * pow(max(dot(MOON_DIRECTION, ray.direction), 0.0), 16.0) * max(0.0, -ray.direction.z);
  return color / (color + 1.0);
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
  vec3 color = vec3(0.1, 0.1, 0.1);
  if (hit.id == TERRAIN_ID) {
    color = vec3(0.2 + hit.position.y, hit.position.y, 0.2);
  } else if (hit.id == TEMPLE_ID) {
    vec3 c1 = texture_map(iChannel0, 0.01 * hit.position, hit.normal, 1.0).rgb;
    vec3 c2 = texture_map(iChannel1, 0.08 * hit.position, hit.normal, 1.0).rgb;
    color = mix(c1, c2, 0.8);
  } else if (hit.id == TEMPLE_WINDOW_ID) {
    color = vec3(0.1, 0.0, 0.0);
  } else if (hit.id == TEMPLE_INSIDE_ID) {
    color = vec3(0.5, 0.5, 0.5);
  }
  return color;
}

vec3 render(Ray ray) {
  vec3 color = vec3(0.0);
  Hit hit = march(ray, 0.0, 500.0, 0.8, 1000);

  if (hit.id != -1.0) {
    vec3 evil_direction = normalize(EVIL_POSITION - hit.position);
    color = 0.2 * get_material(hit);
    Ray moon_shadow_ray = Ray(hit.position, MOON_DIRECTION);
    Ray evil_shadow_ray = Ray(hit.position, -evil_direction);
    float fog = 1.0 - exp(-0.03 * hit.t);
    float moon_shadow = soft_shadow(moon_shadow_ray, 0.2, 500.0, 64);
    float evil_shadow = soft_shadow(evil_shadow_ray, 0.2, 500.0, 64);
    float NdotL = clamp(dot(MOON_DIRECTION, hit.normal), 0.0, 1.0);
    float NdotE = clamp(dot(evil_direction, hit.normal), 0.0, 1.0);
    color += NdotL * MOON_COLOR * moon_shadow;
    color += NdotE * EVIL_COLOR * evil_shadow;
    color = mix(color, get_fog_color(ray), fog);
  } else {
    color = get_sky_color(ray);
  }

  return color;
}

void mainImage(out vec4 out_color, in vec2 in_position) {
  vec2 uv = (2.0 * in_position - iResolution.xy) / iResolution.x;

  vec3 origin = vec3(60.0, 10.0, 100.0);
  vec3 target = vec3(0.0, 30.0, 0.0);

  Ray ray;
  ray.origin = origin;
  ray.direction = normalize(look_at(origin, target) * vec3(uv, -1.0)); 

  vec3 color = render(ray);
  color = pow(color, vec3(1.0 / 2.2));

  out_color = vec4(color, 1.0);
}