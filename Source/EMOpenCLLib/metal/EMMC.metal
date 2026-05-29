//
// EMMC.metal
//
// Metal Shading Language (MSL) translation of opencl/EMMC.cl — the Monte Carlo
// BSE electron-trajectory kernel.  This is a direct, behaviour-preserving port:
// the LFSR113 RNG, the Lambert projection and the MC physics are reproduced
// exactly so that, fed identical prime seeds, the output matches the OpenCL
// kernel.  See MetalMigrationPlan.md (Phase 1).
//
// Argument buffer indices match the OpenCL clSetKernelArg indices used by the
// host (mod_MCOpenCL / mod_SEMCLwrappers / mod_EBSDFull) so the same host code
// drives either backend:
//   0 Lamx  1 Lamy  2 E  3 count  4 z  5 rho  6 A  7 num_el
//   8 seeds 9 sig  10 omega  11 depth  12 energy  13 steps
//
// Scalars are received as `constant T&` (bound via setBytes); array arguments
// as `device T*` (bound via setBuffer).  Apple GPUs have no fp64, but this
// kernel is entirely single precision, so the port is exact.
//

#include <metal_stdlib>
using namespace metal;

#define RAND_MAX  2147483647.0f
#define PI        3.14159f

struct LambertStruct {
    float x;
    float y;
};

struct lfsrret {
    int z1;
    int z2;
    int z3;
    int z4;
    int rand;
};

//--------------------------------------------------------------------------
// lfsr113 combined Tausworthe RNG (identical bit operations to EMMC.cl)
//--------------------------------------------------------------------------
static lfsrret lfsr113_Bits(int z11, int z22, int z33, int z44)
{
    lfsrret ret;
    int z1 = z11, z2 = z22, z3 = z33, z4 = z44;
    int b;
    b  = ((z1 << 6) ^ z1) >> 13;
    z1 = ((z1 & 4294967294U) << 18) ^ b;
    b  = ((z2 << 2) ^ z2) >> 27;
    z2 = ((z2 & 4294967288U) << 2) ^ b;
    b  = ((z3 << 13) ^ z3) >> 21;
    z3 = ((z3 & 4294967280U) << 7) ^ b;
    b  = ((z4 << 3) ^ z4) >> 12;
    z4 = ((z4 & 4294967168U) << 13) ^ b;
    ret.z1 = z1;
    ret.z2 = z2;
    ret.z3 = z3;
    ret.z4 = z4;
    ret.rand = (z1 ^ z2 ^ z3 ^ z4);
    return ret;
}

//--------------------------------------------------------------------------
// Lambert projection of a point of the unit sphere (port of LambertSphereToPlane)
// Takes the direction cosines as a float3 (the OpenCL version took float xyz[3]).
//--------------------------------------------------------------------------
static LambertStruct LambertSphereToPlane(float3 normxyz)
{
    float q;
    float LPssPi2  = 0.886226925452758f;
    float LPssPio2 = 1.253314137315500f;
    LambertStruct ret;
    ret.x = 0.0f;
    ret.y = 0.0f;
    float3 xyz;
    float mag = sqrt(normxyz.x*normxyz.x + normxyz.y*normxyz.y + normxyz.z*normxyz.z);
    xyz.x = normxyz.x/mag;
    xyz.y = normxyz.y/mag;
    xyz.z = normxyz.z/mag;
    if (fabs(xyz.z) == 1.0f) {
        ret.x = 0.0f;
        ret.y = 0.0f;
    }
    else {
        if (fabs(xyz.y) <= fabs(xyz.x) && (xyz.x != 0.0f)) {
            q = (fabs(xyz.x)/xyz.x) * sqrt(2.0f*(1.0f + xyz.z));
            ret.x = q * LPssPi2;
            ret.y = q * atan(xyz.y/xyz.x)/LPssPi2;
        }
        else {
            if (xyz.y != 0.0f) {
                q = (fabs(xyz.y)/xyz.y) * sqrt(2.0f*(1.0f + xyz.z));
                ret.x = q * atan(xyz.x/xyz.y)/LPssPi2;
                ret.y = q * LPssPi2;
            }
        }
    }
    ret.x = ret.x/LPssPio2;
    ret.y = ret.y/LPssPio2;
    return ret;
}

//--------------------------------------------------------------------------
// MC : Monte Carlo BSE electron scattering kernel (port of __kernel void MC)
//--------------------------------------------------------------------------
kernel void MC(device   float* Lamx   [[buffer(0)]],
               device   float* Lamy   [[buffer(1)]],
               constant float& E      [[buffer(2)]],
               constant int&   count  [[buffer(3)]],
               constant float& z      [[buffer(4)]],
               constant float& rho    [[buffer(5)]],
               constant float& A      [[buffer(6)]],
               constant int&   num_el [[buffer(7)]],
               device   int*   seeds  [[buffer(8)]],
               constant float& sig    [[buffer(9)]],
               constant float& omega  [[buffer(10)]],
               device   float* depth  [[buffer(11)]],
               device   float* energy [[buffer(12)]],
               constant int&   steps  [[buffer(13)]],
               uint2 gid [[thread_position_in_grid]])
{
    int tx = int(gid.x);
    int ty = int(gid.y);
    float dir_cos[3];
    LambertStruct ret;

    int id = count*ty + tx;
    float rand;
    int z11, z22, z33, z44;
    lfsrret retrnd;

    int counter1, counter2;
    float4 c_new, r_new;
    float E_new, alpha, de_ds, phi, psi, mfp, sig_eNA, step, dsq, dsqi, absc0z;

    float J;    // Joy, Monte Carlo simulation for Electron Microscopy and Microanalysis
    J = (9.76f*z + 58.5f*powr(z,-0.19f))*1E-3f;

    float4 r0 = float4(0.0f, 0.0f, 0.0f, 0.0f);
    float4 c0 = float4(cos(omega)*sin(sig), sin(omega)*sin(sig), cos(sig), 0.0f);
    float escape_depth;

// Setting all values to -10. Any value other than -10 will denote a backscattered
// electron with the x and y component of the Lambert projection.
    for (int i = 0; i < num_el; ++i){
        Lamx[num_el*id + i]   = -10.0f;
        Lamy[num_el*id + i]   = -10.0f;
        depth[num_el*id + i]  = 10.0f;
        energy[num_el*id + i] = 0.0f;
    }

    for (int i = 0; i < num_el; ++i){
        z11 = seeds[4*id];
        z22 = seeds[4*id + 1];
        z33 = seeds[4*id + 2];
        z44 = seeds[4*id + 3];
        retrnd = lfsr113_Bits(z11,z22,z33,z44);
        seeds[4*id]     = retrnd.z1;
        seeds[4*id + 1] = retrnd.z2;
        seeds[4*id + 2] = retrnd.z3;
        seeds[4*id + 3] = retrnd.z4;
        rand = fabs(retrnd.rand/RAND_MAX);
        r0 = float4(0.0f, 0.0f, 0.0f, 0.0f);
        c0 = float4(cos(omega)*sin(sig), sin(omega)*sin(sig), cos(sig), 0.0f);
        E_new = E;
        c_new = c0;
        escape_depth = 0.0f;
        alpha = (3.4E-3f)*powr(z,0.66667f)/E_new;
        sig_eNA = (5.21f * 602.2f)*z*z/E_new/E_new*4.0f*PI/alpha/(1.0f+alpha)*pow(E_new+511.0f,2.0f)/pow(E_new+1024.0f,2.0f);

        mfp = A * 1.0e7f/(rho*sig_eNA);
        step = -mfp * log(rand);
        r_new = r0 + step*c_new;
        r0 = r_new;
        de_ds = -0.00785f*(z/(A*E_new)) * log(1.166f*E_new/J + 0.9911f);
        E_new += step*rho*de_ds;

        counter1 = 0;
        counter2 = 0;

        while (counter1 < steps){
            alpha = (3.4e-3f)*powr(z,0.66667f)/E_new;
            sig_eNA = (5.21f * 602.2f)*z*z/E_new/E_new*4.0f*PI/alpha/(1.0f+alpha)*pow(E_new+511.0f,2.0f)/pow(E_new+1024.0f,2.0f);
            mfp = A * 1.0e7f/(rho*sig_eNA);

            z11 = seeds[4*id];
            z22 = seeds[4*id + 1];
            z33 = seeds[4*id + 2];
            z44 = seeds[4*id + 3];
            retrnd = lfsr113_Bits(z11,z22,z33,z44);
            seeds[4*id]     = retrnd.z1;
            seeds[4*id + 1] = retrnd.z2;
            seeds[4*id + 2] = retrnd.z3;
            seeds[4*id + 3] = retrnd.z4;
            rand = fabs(retrnd.rand/RAND_MAX);
            step = -mfp * log(rand);

            de_ds = -0.00785f*(z/(A*E_new)) * log(1.166f*E_new/J + 0.9911f);

            z11 = seeds[4*id];
            z22 = seeds[4*id + 1];
            z33 = seeds[4*id + 2];
            z44 = seeds[4*id + 3];
            retrnd = lfsr113_Bits(z11,z22,z33,z44);
            seeds[4*id]     = retrnd.z1;
            seeds[4*id + 1] = retrnd.z2;
            seeds[4*id + 2] = retrnd.z3;
            seeds[4*id + 3] = retrnd.z4;
            rand = fabs(retrnd.rand/RAND_MAX);
            phi = acos(1.0f - ((2.0f*alpha*rand)/(1.0f + alpha - rand)));

            z11 = seeds[4*id];
            z22 = seeds[4*id + 1];
            z33 = seeds[4*id + 2];
            z44 = seeds[4*id + 3];
            retrnd = lfsr113_Bits(z11,z22,z33,z44);
            seeds[4*id]     = retrnd.z1;
            seeds[4*id + 1] = retrnd.z2;
            seeds[4*id + 2] = retrnd.z3;
            seeds[4*id + 3] = retrnd.z4;
            rand = fabs(retrnd.rand/RAND_MAX);
            psi = 2.0f*PI*rand;

// new direction cosines of the electrons after scattering event
            if ((c0.z >= 0.99999f) || (c0.z <= -0.99999f) ){
                absc0z = fabs(c0.z);
                c_new = float4(sin(phi) * cos(psi), sin(phi) * sin(psi), (c0.z/absc0z)*cos(phi), 0.0f);
            }
            else {
                dsq = sqrt(1.0f-c0.z*c0.z);
                dsqi = 1.0f/dsq;
                c_new = float4(sin(phi)*(c0.x*c0.z*cos(psi) - c0.y*sin(psi))*dsqi + c0.x*cos(phi), sin(phi) * (c0.y * c0.z * cos(psi) + c0.x * sin(psi)) * dsqi + c0.y * cos(phi), -sin(phi) * cos(psi) * dsq + c0.z * cos(phi), 0.0f);
            }

            if (fabs(c_new.z) > 1.0E-5f){
                escape_depth = r_new.z/c_new.z;
            }

            r_new = r0 + step*c_new;

            r0 = r_new;
            c0 = c_new;
            E_new += step*rho*de_ds;
            if (r0.z <= 0 && counter2 == 0){
                dir_cos[0] = c0.x;
                dir_cos[1] = c0.y;
                dir_cos[2] = c0.z;
                if(dir_cos[0] != 0.0f && dir_cos[1] != 0.0f && dir_cos[2] != 0.0f){
                    ret = LambertSphereToPlane(float3(dir_cos[0], dir_cos[1], dir_cos[2]));
                }
                Lamx[num_el*id + i]   = ret.x;
                Lamy[num_el*id + i]   = ret.y;
                depth[num_el*id + i]  = escape_depth;
                energy[num_el*id + i] = E_new;
                counter2 = 1;
            }

            counter1++ ;
        }
    }
}
