#ifndef COMMON_FRAG
#define COMMON_FRAG

#define AIR_IOR   1.0
#define WATER_IOR 1.33

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



#endif // COMMON_FRAG