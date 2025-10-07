#include "common.frag"
#iChannel0 "simple_box.frag"

#define NUMBER_OF_BOUNCES (4)
#define NUMBER_OF_SAMPLES (4)

void add_hit(inout Hit current, in Hit next) {
  if ((next.t > 0.0 && next.t < current.t) || (current.t < 0.0)) current = next;
}

// Front, Back, Left, Right, Top, Bottom
Hit box_intersect(Ray ray, float x, float y, float z, int id) {
  vec3 [30] triangles; set_box(triangles, x, y, z);
  Hit hit; hit.id = -1; hit.t = INFINITY;
  for (int i = 0; i < triangles.length() / 3; i++) {
    vec3 p0 = triangles[3 * i + 0];
    vec3 p1 = triangles[3 * i + 1];
    vec3 p2 = triangles[3 * i + 2];
    vec4 t_and_normal = triangle_intersect(ray, p0, p1, p2);
    if (t_and_normal.x != INFINITY) {
      hit.t = min(hit.t, t_and_normal.x);
      hit.id = id + i / 2;
      hit.normal = t_and_normal.yzw;
      hit.position = ray.origin + hit.t * ray.direction;
    }
  }
  return hit;
}


Material get_material(Hit hit) {
  Material material;
  if (hit.id == 0)       material = Material(vec3(0.4, 0.5, 0.5), DIFFUSE);
  else if (hit.id == 1)  material = Material(vec3(0.4, 0.2, 0.1), DIFFUSE);
  else if (hit.id == 2)  material = Material(vec3(3.5, 1.2, 3.7), EMISSIVE);
  else if (hit.id == 3)  material = Material(vec3(0.8, 0.8, 0.8), DIFFUSE);
  else if (hit.id == 4)  material = Material(vec3(0.1, 0.6, 0.6), DIFFUSE);
  else if (hit.id == 5)  material = Material(vec3(4.5, 4.1, 4.2), EMISSIVE);
  else if (hit.id == 6)  material = Material(vec3(0.5, 0.2, 0.4), DIFFUSE);
  else if (hit.id == 7)  material = Material(vec3(0.1, 0.2, 0.6), DIFFUSE);
  else if (hit.id == 8)  material = Material(vec3(0.5, 0.1, 0.2), DIFFUSE);
  else if (hit.id == 9)  material = Material(vec3(0.2, 0.7, 0.3), DIFFUSE);
  else if (hit.id == 10) material = Material(vec3(0.2, 0.1, 0.2), DIFFUSE);
  else if (hit.id == 11) material = Material(vec3(0.2, 0.4, 0.4), DIFFUSE);
  return material;
}


Hit intersect_scene(Ray ray) {
  Hit hit = sphere_intersect(ray, vec3(+7.0, +2.5, +4.0), 2.5, 0);
  add_hit(hit, sphere_intersect(ray, vec3(-7.0, -2.0, +0.0), 3.0, 1));
  add_hit(hit, sphere_intersect(ray, vec3(-6.0, -8.0, +8.0), 2.0, 2));
  add_hit(hit, sphere_intersect(ray, vec3(+6.0, -7.0, +9.0), 1.5, 3));
  add_hit(hit, sphere_intersect(ray, vec3(+3.0, -4.0, -7.0), 3.0, 4));
  add_hit(hit, sphere_intersect(ray, vec3(+0.0, +7.0, +0.0), 3.0, 5));
  add_hit(hit, sphere_intersect(ray, vec3(+0.0, -7.5, +2.0), 2.5, 6));
  add_hit(hit, box_intersect(ray, 10.0, 10.0, 10.0, 7));
  return hit;
}

vec3 trace(Ray ray, int bounces, inout uint seed) {
  vec3 total = vec3(0.0);
  vec3 throughput = vec3(1.0);

  for (int i = 0; i < bounces; i++) {
    Hit hit = intersect_scene(ray);

    if (hit.id == -1) {
      total += vec3(0.1);
      break;
    }
    
    Material material = get_material(hit);

    ray.origin = hit.position + hit.normal * 0.01;
    ray.direction = random_hemisphere_unit_vector(hit.normal, seed);

    if (material.type == EMISSIVE) {
      total += throughput * material.albedo;
    }

    throughput *= material.albedo;
  }

  return total;
}

void mainImage(out vec4 out_color, in vec2 in_position) {
  vec2 uv = (2.0 * in_position - iResolution.xy) / iResolution.y;

  uint seed = (uint(in_position.x) * 1973u + uint(in_position.y) * 9277u + uint(iFrame) * 26699u) | 1u;

  Ray ray;
  ray.origin = vec3(0.0, 0.0, 20.0);
  ray.direction = normalize(vec3(uv, -1.0));
  vec3 total = vec3(0.0);

  for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
    vec2 origin_offset = (vec2(random(seed), random(seed)) - vec2(0.5)) / iResolution.y;
    ray.direction = normalize(vec3(uv + origin_offset, -1.0));
    total += trace(ray, NUMBER_OF_BOUNCES, seed) / float(NUMBER_OF_SAMPLES); 
  }

  vec3 current_color = texture(iChannel0, in_position / iResolution.xy).rgb;

  total = mix(current_color, total, 1.0f / float(iFrame + 1));

  out_color = vec4(total, 1.0);
}