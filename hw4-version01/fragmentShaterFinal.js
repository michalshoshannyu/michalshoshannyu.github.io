let fragmentShader = `
precision mediump float;

uniform float uTime, uFL;
uniform vec3 uCursor;
uniform vec2 uR;
uniform mat4 uA[` + NQ + `];  // Quadric coefficient matrices for object A
uniform mat4 uB[` + NQ + `];  // Quadric coefficient matrices for object B
uniform mat4 uC[` + NQ + `];  // Quadric coefficient matrices for object C


varying vec3 vPos;

float noise(vec3 point) { float r = 0.; for (int i=0;i<16;i++) {
  vec3 D, p = point + mod(vec3(i,i/4,i/8) , vec3(4.0,2.0,2.0)) +
       1.7*sin(vec3(i,5*i,8*i)), C=floor(p), P=p-C-.5, A=abs(P);
  C += mod(C.x+C.y+C.z,2.) * step(max(A.yzx,A.zxy),A) * sign(P);
  D=34.*sin(987.*float(i)+876.*C+76.*C.yzx+765.*C.zxy);P=p-C-.5;
  r+=sin(6.3*dot(P,fract(D)-.5))*pow(max(0.,1.-2.*dot(P,P)),4.);
} return .5 * sin(r); }

float turbulence(vec3 P) {
   float f = 0., s = 1.;
   for (int i = 0 ; i < 9 ; i++) {
      f += abs(noise(s * P)) / s;
      s *= 2.;
      P = vec3(.866*P.x + .5*P.z, P.y + 100., -.5*P.x + .866*P.z);
   }
   return f;
}

vec3 marble(vec3 pos) {
   float v = turbulence(pos);
   float s = sqrt(.5 + .5 * sin(20. * pos.y + 8. * v));
   return vec3(.8,.7,.5) * vec3(s,s*s,s*s*s);
}

mat3 refraction(vec3 V, vec3 W, vec3 N, float n1, float n2) {
  vec3 C3 = dot(N, W) * N;
  vec3 S3 = W - C3;
  vec3 S4 = S3 * n2 / n1;
  vec3 C4 = - sqrt(1. - S4 * S4) * N;
  vec3 W3 = C4 + S4;

  vec3 P = V + 0.00001 * W3;
  return mat3(P, W3, 0, 0, 0);
}

mat4 coefficients(mat4 Q) {
    float a = Q[0].x; 
    float b = Q[1].y; 
    float c = Q[2].z; 
    float d = Q[1].z + Q[2].y;
    float e = Q[0].z + Q[2].x; 
    float f = Q[0].y + Q[1].x;
    float g = Q[0].w + Q[3].x; 
    float h = Q[1].w + Q[3].y; 
    float i = Q[2].w + Q[3].z; 
    float j = Q[3].w;

    return mat4(a, 0, 0, 0,
                f, b, 0, 0,
                e, d, c, 0,
                g, h, i, j);
}

vec2 rayIntersection(vec3 V, vec3 W, mat4 Q) {
    mat4 quadricQ = coefficients(Q);
    vec4 V4 = vec4(V, 1.);
    vec4 W4 = vec4(W, 0.);

    float a = quadricQ[0].x;
    float f = quadricQ[1].x;
    float b = quadricQ[1].y;
    float e = quadricQ[2].x;
    float d = quadricQ[2].y;
    float c = quadricQ[2].z;
    float g = quadricQ[3].x;
    float h = quadricQ[3].y;
    float i = quadricQ[3].z;
    float j = quadricQ[3].w;

    float A = a * W4.x * W4.x + b * W4.y * W4.y + c * W4.z * W4.z +
              d * W4.y * W4.z + e * W4.z * W4.x + f * W4.x * W4.y;

    float B = 2.0 * (a * V4.x * W4.x + b * V4.y * W4.y + c * V4.z * W4.z) +
              d * (V4.y * W4.z + V4.z * W4.y) +
              e * (V4.z * W4.x + V4.x * W4.z) +
              f * (V4.x * W4.y + V4.y * W4.x);

    float C = a * V4.x * V4.x + b * V4.y * V4.y + c * V4.z * V4.z +
              d * V4.y * V4.z + e * V4.z * V4.x + f * V4.x * V4.y +
              g * V4.x + h * V4.y + i * V4.z + j;

    float discriminant = B * B - 4. * A * C;
    if (discriminant < 0.0) {
        return vec2(-1.0, -1.0);  // No intersection
    }

    float sqrtDiscriminant = sqrt(discriminant);
    float t0 = (-B - sqrtDiscriminant) / (2.0 * A);
    float t1 = (-B + sqrtDiscriminant) / (2.0 * A);
    
    return vec2(t0, t1);
}

vec3 norm(mat4 Q, vec3 P){
    mat4 quadricQ = coefficients(Q); 

    float a = quadricQ[0].x;
    float f = quadricQ[1].x;
    float b = quadricQ[1].y;
    float e = quadricQ[2].x;
    float d = quadricQ[2].y;
    float c = quadricQ[2].z;
    float g = quadricQ[3].x;
    float h = quadricQ[3].y;
    float i = quadricQ[3].z;
    float j = quadricQ[3].w;

    float x = P.x;
    float y = P.y;
    float z = P.z;

    vec3 N = normalize(vec3(
        2.0*a*x + e*z + f*y + g,
        2.0*b*y + d*z + f*x + h,
        2.0*c*z + d*y + e*x + i
    ));
    return N;
}

mat3 rayQuadric(vec3 V, vec3 W, mat4 quadricA, mat4 quadricB, mat4 quadricC, int enterExit) {

    vec2 rayAIntersection = rayIntersection(V, W, quadricA);
    vec2 rayBIntersection = rayIntersection(V, W, quadricB);
    vec2 rayCIntersection = rayIntersection(V, W, quadricC);

    float rayAEnter = rayAIntersection.x;
    float rayBEnter = rayBIntersection.x;
    float rayCEnter = rayCIntersection.x;

    float rayAExit = rayAIntersection.y;
    float rayBExit = rayBIntersection.y;
    float rayCExit = rayCIntersection.y;

    float tExit = min(min(rayAExit, rayBExit), rayCExit);
    mat4 exitMatrix = quadricA;

    if (tExit == rayBExit) {
        exitMatrix = quadricB;
    }
    if (tExit < rayCExit) {
        exitMatrix = quadricC;
    }

    float tEnter = max(max(rayAEnter, rayBEnter), rayCEnter);
    mat4 enterMatrix = quadricA;

    if (tEnter == rayBEnter) {
        enterMatrix = quadricB;
    }
    if (tEnter == rayCEnter) {
        enterMatrix = quadricC;
    }

    if (enterExit == 1) {
        vec3 P = V + tExit * W;
        vec3 N = -norm(exitMatrix, P);
        return mat3(tEnter, tExit, 0, P, N);
    } else {
        vec3 P = V + tEnter * W;
        vec3 N = norm(enterMatrix, P);
        return mat3(tEnter, tExit, 0, P, N);
    }
}

vec3 outRefraction(vec3 color, vec3 V, vec3 W) {
    float closestDist = 200.;

    vec3 oColor = color;
    for (int i = 0; i < `+NQ+`; i++) {

        mat4 quadricA = uA[i];
        mat4 quadricB = uB[i];
        mat4 quadricC = uC[i];

        mat3 rayQ = rayQuadric(V, W, quadricA, quadricB, quadricC, -1);

        float rayExit = rayQ[0].y;
        float rayEnter = rayQ[0].x;
        vec3 P = rayQ[1];
        vec3 N = rayQ[2];

        if (rayExit > rayEnter && rayEnter > 0. && rayEnter < closestDist) {
            closestDist = rayEnter;

            vec3 reflectDir = W - 2. * N * dot(N, W);
            vec3 lightDir = normalize(vec3(-10.0, 3.0,1.0));
            vec3 lightColor = oColor * max(0., dot(N, lightDir));
            vec3 lightSpecColor = oColor * pow(max(0., dot(reflectDir, lightDir)), 0.8);
            oColor += oColor * lightColor;
            oColor += lightSpecColor;
        }
    }
    return oColor;
}

vec3 shaderRefraction(vec3 color, vec3 W, vec3 P, vec3 N, mat4 quadricA, mat4 quadricB, mat4 quadricC) {

    vec3 oColor = color;
    vec3 reflectDir = W - 2. * N * dot(N, W);
    vec3 lightDir = normalize(vec3(-10.0, 3.0,1.0));
    vec3 lightColor = oColor * max(0., dot(N, lightDir));
    vec3 lightSpecColor = oColor * pow(max(0., dot(reflectDir, lightDir)), 0.3);
    oColor += oColor * lightColor;
    oColor += lightSpecColor;

    mat3 enterRef = refraction(P, W, N, 1.2, 1.5);
    vec3 enterRefDir = enterRef[1];
    vec3 enterP = enterRef[0];

    mat3 rayQ = rayQuadric(enterP, enterRefDir, quadricA, quadricB, quadricC, 1);
    float refracExit = rayQ[0].y;

    vec3 exitP = rayQ[1];
    vec3 exitN = rayQ[2];

    mat3 exitRef = refraction(exitP, enterRefDir, exitN, 1.2, 1.5);
    vec3 exitRefDir = exitRef[1];
    vec3 exitedP = exitRef[0];

    vec3 refractedColor = outRefraction(oColor, exitedP, exitRefDir);
    oColor = mix(refractedColor, color, 0.5);

    return oColor;
}

vec3 shadersRefraction(vec3 color, vec3 P, vec3 W) {
    float closestDist = 200.;

    vec3 oColor = color;
    for (int i = 0; i < `+NQ+`; i++) {

        mat4 quadricA = uA[i];
        mat4 quadricB = uB[i];
        mat4 quadricC = uC[i];

        mat3 rayQ = rayQuadric(P, W, quadricA, quadricB, quadricC, -1);

        float rayeExit = rayQ[0].y;
        float rayEnter = rayQ[0].x;
        vec3 newP = rayQ[1];
        vec3 N = rayQ[2];

        if (rayeExit > rayEnter && rayEnter > 0. && rayEnter < closestDist) {
            closestDist = rayEnter;
            oColor = shaderRefraction(color, W, newP, N, quadricA, quadricB, quadricC);
        }
    }
    return oColor;
}

vec3 quadricShader(vec3 color, vec3 W, vec3 P, vec3 N, mat4 quadricA, mat4 quadricB, mat4 quadricC, vec3 L) {

    vec3 oColor = vec3(1.0, 0.5, 0.0);
    
    vec3 reflectDir = W - 2. * N * dot(N, W);
    vec3 lightDir = vec3(1.0, 0.0, 1.0);
    vec3 lightColor = oColor * max(0., dot(N, lightDir));
    vec3 lightSpecColor = oColor * pow(max(0., dot(reflectDir, lightDir)), 0.1);
    float noi = abs(noise(lightSpecColor * reflectDir)) ;
    float na = noise(vPos * 100.0) ;
    oColor += oColor * lightColor;
    oColor += lightSpecColor + na * noi;

    
    mat3 entereRef = refraction(P, W, N, 2., 6.);
    vec3 enterRefDir = entereRef[1];
    vec3 enterP = entereRef[0];

    mat3 rayQ = rayQuadric(enterP, enterRefDir, quadricA, quadricB, quadricC, 1);
    float refracExit = rayQ[0].y;
    float refracEntert = rayQ[0].x;

    vec3 ExitP = rayQ[1];
    vec3 exitN = rayQ[2];

    mat3 exitRef = refraction(ExitP, enterRefDir, exitN, 1.2, 1.5);
    vec3 exitRefDir = exitRef[1];
    vec3 exitedP = exitRef[0];

    vec3 refractedColor = shadersRefraction(oColor, exitedP, exitRefDir);
    oColor = mix(refractedColor, oColor, 0.8)+noi;

    if (refracEntert > 0.0) {
        vec3 reP = P + .001 * reflectDir;
        vec3 reflectedColor = shadersRefraction(color, reP, reflectDir) + na * noi;
        oColor = mix(color, reflectedColor, 0.5);
    }


    return oColor;
}


vec3 quadricsShader(vec3 color, vec3 V, vec3 W, vec3 L) {
  float closestDist = 200.;

  vec3 oColors = color;
  for (int i = 0; i < `+NQ+`; i++) {
    mat4 quadricA = uA[i];
    mat4 quadricB = uB[i];
    mat4 quadricC = uC[i];

    mat3 rayQ = rayQuadric(V, W, quadricA, quadricB, quadricC, -1);

    float rayExit = rayQ[0].y;
    float rayEnter = rayQ[0].x;

    if (rayExit > rayEnter && rayEnter > 0. && rayEnter < closestDist) {
      closestDist = rayEnter;

      vec3 P = rayQ[1];
      vec3 N = rayQ[2];
      oColors = quadricShader(color, W, P, N, quadricA, quadricB, quadricC, L);
    }
  }
  return oColors;
}

void main(void) {
    float na = noise(vPos * 100.0) ;
    vec3 litColor = vec3(1.0, 0.9, 0.8);
    vec3 origin = vec3(.0, .0, 10.);
    vec3 color = vec3(.0, .4, .2+na);

    
    vec4 pos = vec4(vPos, 1.0);
    vec4 transformedPosition = uA[0] * pos;  
    transformedPosition.xy *= uR / min(uR.x, uR.y);
    transformedPosition.xy += origin.xy;
    vec3 W = normalize(vec3(transformedPosition.xy, -uFL));
    vec3 lightDir = normalize(vec3(-5.0, 10.0, 5.0)); 

    color = quadricsShader(color, origin, W, lightDir);
    color = clamp(color, 0.0, 1.0);

    gl_FragColor = vec4(sqrt(color), 1.0);
}

`;