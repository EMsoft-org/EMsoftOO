//
// EMMCxyz.metal
//
// MSL translation of opencl/EMMCxyz.cl — the Monte Carlo kernel used by
// mod_MCOpenCL's 'Ivol' (interaction-volume) mode.  Same LFSR113 RNG and MC
// physics as EMMC.metal, but the output is the (x,y,z) exit position of each
// backscattered electron rather than a Lambert projection (no Lambert here).
//
// Argument buffer indices match the OpenCL clSetKernelArg order in mod_MCOpenCL
// (Ivol branch):
//   0 Lamx  1 Lamy  2 Lamz  3 E  4 count  5 z  6 rho  7 A
//   8 num_el 9 seeds 10 sig 11 omega 12 steps
//
// See MetalMigrationPlan.md (Phase 3).
//

#include <metal_stdlib>
using namespace metal;

#define RAND_MAX  2147483647.0f
#define PI        3.14159f

struct lfsrret {
    int z1;
    int z2;
    int z3;
    int z4;
    int rand;
};

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

kernel void MCxyz(device   float* Lamx   [[buffer(0)]],
                  device   float* Lamy   [[buffer(1)]],
                  device   float* Lamz   [[buffer(2)]],
                  constant float& E      [[buffer(3)]],
                  constant int&   count  [[buffer(4)]],
                  constant float& z      [[buffer(5)]],
                  constant float& rho    [[buffer(6)]],
                  constant float& A      [[buffer(7)]],
                  constant int&   num_el [[buffer(8)]],
                  device   int*   seeds  [[buffer(9)]],
                  constant float& sig    [[buffer(10)]],
                  constant float& omega  [[buffer(11)]],
                  constant int&   steps  [[buffer(12)]],
                  uint2 gid [[thread_position_in_grid]])
{
    int tx = int(gid.x);
    int ty = int(gid.y);

    int id = count*ty + tx;
    float rand;
    int z11, z22, z33, z44;
    lfsrret retrnd;

    int counter1, counter2;
    float4 c_new, r_new;
    float E_new, alpha, de_ds, phi, psi, mfp, sig_eNA, step, dsq, dsqi, absc0z;

    float J;
    J = (9.76f*z + 58.5f*powr(z,-0.19f))*1E-3f;

    float4 r0 = float4(0.0f, 0.0f, 0.0f, 0.0f);
    float4 c0 = float4(cos(omega)*sin(sig), sin(omega)*sin(sig), cos(sig), 0.0f);
    float escape_depth;

    for (int i = 0; i < num_el; ++i){
        Lamx[num_el*id + i] = -100000.0f;
        Lamy[num_el*id + i] = -100000.0f;
        Lamz[num_el*id + i] = -100000.0f;
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

        alpha = (3.4E-3f)*powr(z,0.67f)/E_new;
        sig_eNA = (5.21f * 602.3f)*((z*z)/(E_new*E_new))*((4.0f*PI)/(alpha*(1+alpha)))*((E_new + 511.0f)*(E_new + 511.0f)/((E_new + 1024.0f)*(E_new + 1024.0f)));
        mfp = A/(rho*sig_eNA);
        step = -mfp * log(rand);
        r_new = r0 + step*c_new*1.0e7f;
        r0 = r_new;
        counter1 = 0;
        counter2 = 0;

        while (counter1 < steps){
            alpha = (3.4E-3f)*powr(z,0.67f)/E_new;
            sig_eNA = (5.21f * 602.3f)*((z*z)/(E_new*E_new))*((4*PI)/(alpha*(1+alpha)))*((E_new + 511.0f)*(E_new + 511.0f)/((E_new + 1024.0f)*(E_new + 1024.0f)));
            mfp = A/(rho*sig_eNA);

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

            de_ds = -78500.0f*(z/(A*E_new)) * log(1.166f*E_new/J + 0.9911f);

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
            phi = acos(1 - ((2*alpha*rand)/(1 + alpha - rand)));

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
            psi = 2*PI*rand;

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
            r_new = r0 + step*c_new*1.0e7f;
            c0 = c_new;
            E_new += step*rho*de_ds;
            if (r_new.z <= 0 && counter2 == 0){
                Lamx[num_el*id + i] = r0.x;
                Lamy[num_el*id + i] = r0.y;
                Lamz[num_el*id + i] = r0.z;
                counter2 = 1;
            }
            r0 = r_new;

            counter1++ ;
        }
    }
}
