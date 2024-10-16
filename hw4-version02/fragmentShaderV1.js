let fragmentShader = `
precision mediump float;

uniform float uTime, uFL;
uniform vec3 uCursor;
varying vec3 vPos;
uniform vec2 uR;


uniform mat4 uA[` + NQ + `];  // Quadric coefficient matrices for object A
uniform mat4 uB[` + NQ + `];  // Quadric coefficient matrices for object B
uniform mat4 uC[` + NQ + `];  // Quadric coefficient matrices for object C


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


// Function to calculate ray-quadratic intersection
vec2 rayQ(vec3 V, vec3 W, float a, float b, float c, float d, float e, float f, float g, float h, float i, float j) {
   float A = a * W.x * W.x + b * W.y * W.y + c * W.z * W.z +
              d * W.y * W.z + e * W.z * W.x + f * W.x * W.y;

    float B = 2.0 * (a * V.x * W.x + b * V.y * W.y + c * V.z * W.z) +
              d * (V.y * W.z + V.z * W.y) +
              e * (V.z * W.x + V.x * W.z) +
              f * (V.x * W.y + V.y * W.x);

    float C = a * V.x * V.x + b * V.y * V.y + c * V.z * V.z +
              d * V.y * V.z + e * V.z * V.x + f * V.x * V.y +
              g * V.x + h * V.y + i * V.z + j;

    float discriminant = B * B - 4.0 * A * C;
    if (discriminant < 0.0) {
        return vec2(-1.0, -1.0);  // No intersection
    }
    
    float sqrtDiscriminant = sqrt(discriminant);
    float t0 = (-B - sqrtDiscriminant) / (2.0 * A);
    float t1 = (-B + sqrtDiscriminant) / (2.0 * A);
    
    return vec2(t0, t1);  // Return entry and exit points
}

// Function to calculate intersection with three quadrics and return t values, point, and normal
mat3 rayQ3(vec3 V, vec3 W, mat4 A, mat4 B, mat4 C, int inOut) {

       
        float a = A[0].x; 
        float b = A[1].y; 
        float c = A[2].z; 
        float d = A[1].z + A[2].y;
        float e = A[0].z + A[2].x; 
        float f = A[0].y + A[1].x;
        float g = A[0].w + A[3].x; 
        float h = A[1].w + A[3].y; 
        float i = A[2].w + A[3].z; 
        float j = A[3].w;

        float bA = B[0].x;
        float bB = B[1].y;
        float bC = B[2].z;
        float bD = B[1].z + B[2].y;
        float bE = B[0].z + B[2].x;
        float bF = B[0].y + B[1].x;
        float bG = B[0].w + B[3].x;
        float bH = B[1].w + B[3].y;
        float bI = B[2].w + B[3].z;
        float bJ = B[3].w;

        float cA = C[0].x;
        float cB = C[1].y;
        float cC = C[2].z;
        float cD = C[1].z + C[2].y;
        float cE = C[0].z + C[2].x;
        float cF = C[0].y + C[1].x;
        float cG = C[0].w + C[3].x;
        float cH = C[1].w + C[3].y;
        float cI = C[2].w + C[3].z;
        float cJ = C[3].w;

        

    vec2 tA = rayQ(V, W,  a, b, c, d, e, f, g, h, i, j);
    vec2 tB = rayQ(V, W,  bA, bB, bC, bD, bE, bF, bG, bH, bI,bJ);
    vec2 tC = rayQ(V, W, cA, cB, cC, cD, cE, cF, cG, cH, cI, cJ);

    // Take the maximum entry point and the minimum exit point
    float tIn = max(max(tA.x, tB.x), tC.x);
    float tOut = min(min(tA.y, tB.y), tC.y);
    
    mat3 result;
    result[0] = vec3(tIn, tOut, 0.0);  // t values
    
    if (tIn > 0.0 && tIn < tOut) {
        // Calculate intersection point P
        vec3 P = V + tIn * W;
        
        mat4 Q;
        if (inOut == 0) {
            Q = A; 
        } else if (inOut == 1) {
            Q = C;  
        } else if (inOut == 2) {
            Q = B; 
        }

        vec3 N = normalize(vec3(
            2.0 * Q[0].x * P.x + Q[0].y * P.y + Q[0].z * P.z + Q[0].w,
            2.0 * Q[1].y * P.y + Q[1].x * P.x + Q[1].z * P.z + Q[1].w,
            2.0 * Q[2].z * P.z + Q[2].x * P.x + Q[2].y * P.y + Q[2].w
        ));
        
        result[1] = P;  // Surface point
        result[2] = N;  // Surface normal

        } else {
            result[1] = vec3(0.0);  // No intersection
            result[2] = vec3(0.0);  // No normal
        }
        
        return result;
}


vec3 refractRay(vec3 W1, vec3 N, float n1, float n2) {
    float eta = n1 / n2;
    float cosThetaI = -dot(N, W1);
    float sin2ThetaT = eta * eta * (1.0 - cosThetaI * cosThetaI);
    
    if (sin2ThetaT > 1.0) {
        return vec3(0.0);  
    }
    
    float cosThetaT = sqrt(1.0 - sin2ThetaT);
    return eta * W1 + (eta * cosThetaI - cosThetaT) * N;
}

const vec3 bgColor = vec3(0.0, 0.0, .5); // Blue background

void main() {
    vec3 rayOrigin = vec3(0.0, 0.0,7.0); 
    //vec3 rayOrigin = vec3(vPos.xyz); 
    vec3 rayDir = normalize(vec3(vPos.xy, -uFL)); 

    // vec4 pos = vec4(vPos, 1.0);
    // vec4 transformedPosition = uA[0] * pos; 
    // transformedPosition.xy *= uR / min(uR.x, uR.y);
    // transformedPosition.xy += rayOrigin.xy;
    // vec3 W = normalize(vec3(transformedPosition.xy, -uFL));

    //debuging
    // gl_FragColor = vec4((rayDir + 1.0) / 2.0, 1.0);
    // return; 

    vec3 lightDir = normalize(vec3(-10.0, 3.0,1.0)); // Light direction
    vec3 color = bgColor; 

    // Iterate over all shapes
    for (int k = 0; k < ` + NQ + `; k++) {
        mat4 A = uA[k];
        mat4 B = uB[k];
        mat4 C = uC[k];

        mat3 R = rayQ3(rayOrigin, rayDir, A, B, C, 0); 
        
        float tIn = R[0].x;
        float tOut = R[0].y;

        //debuging
    //    gl_FragColor = vec4(tIn, tOut, 0.0, 1.0);
    //     return; 
        
        if (tIn > 0. && tIn < tOut) {
            vec3 P = R[1];  // Surface point
            vec3 N = R[2];  // Surface normal


            // Reflection direction
            vec3 reflectionDir = rayDir - 2.0 * dot(N, rayDir) * N;
            vec3 reflectionOrigin = P + N * 0.001;  // Offset - to avoid self-intersection

            //compute color contribution
            vec3 reflectionColor = vec3(0.0); 
           
            bool inShadow = false;
       
            for (int m = 0; m < ` + NQ + `; m++) {
                if (m == k) continue;  
                vec3 viewDir = normalize(-rayDir);  // Ray direction points toward camera
                
                // Handle interactions with other shapes
                mat3 reflectedHit = rayQ3(P , lightDir, uA[m], uB[m], uC[m], 1);  // inOut = 1 for exit
                if (reflectedHit[0].x > 0.0 && reflectedHit[0].x < tOut) {
                    float reflectedTIn = reflectedHit[0].x;
                    float reflectedTOut = reflectedHit[0].y;

                    if (reflectedTIn > 0. && reflectedTIn < reflectedTOut) {
                    vec3 reflectedP = reflectedHit[1];  // New intersection point
                    vec3 reflectedN = reflectedHit[2];  // New surface normal

                     // Add shading for reflected point
                    // vec3 reflectionDir = reflect(-lightDir, reflectedN);
                    // float specularStrength = 0.8;  // Control how intense the specular highlight is
                    // float shininess = 32.0;        // Control how sharp the specular highlight is
                    // vec3 viewDir = normalize(-rayDir);  // Ray direction points toward camera
                    // float spec = pow(max(dot(reflectionDir, viewDir), 0.0), shininess);
                    // vec3 specular = specularStrength * spec * vec3(1.0);

                    float reflectedDiffuse = max(dot(reflectedN, lightDir), 0.0);
                    // reflectionColor = vec3(1.0) * (reflectedDiffuse + specular);  // Add specular reflection

                     reflectionColor = vec3(1.0) * reflectedDiffuse;  // Simple shading for now
                  }
                   color = mix(color, reflectionColor, 0.5);  // Blend between object color and reflection
                    inShadow = true; 
                }
            }
           
            if (inShadow) {
             inShadow = true; 
            } else {
            float diffuse = max(dot(N, lightDir), 0.0);
             vec3 R = rayDir - 2. * N * dot(N, rayDir);
             vec3 s = color * pow(max(dot(R, lightDir),0.), 8.);
             float na = noise(vPos * 200.0) ;
             float noi = abs(noise(s * R)) ;
             color = vec3(2.0, 1.0, 0.0) + na + diffuse + s + noi ; 
            }
            
         }
    }
  
    gl_FragColor = vec4(color, 1.0);  // Output the final color
}

`;
