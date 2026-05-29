//
// EMMCfoil.metal
//
// MSL translation of opencl/EMMCfoil.cl — the Monte Carlo kernel used by
// mod_MCOpenCL's 'foil' mode.  Same LFSR113 RNG, Lambert projection and MC
// physics as EMMC.metal, but for a foil geometry: electrons are tracked through
// a slab of the given `thickness` and those transmitted out the far (z >=
// thickness) side are accumulated in the southern hemisphere (LamxSH/LamySH).
//
// Argument buffer indices match the OpenCL clSetKernelArg order in mod_MCOpenCL
// (foil branch):
//   0 E  1 count  2 z  3 rho  4 A  5 num_el  6 seeds  7 sig  8 omega
//   9 depth  10 energy  11 steps  12 thickness  13 LamxSH  14 LamySH
//
// See MetalMigrationPlan.md (Phase 3).
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

kernel void MC(constant float& E         [[buffer(0)]],
               constant int&   count     [[buffer(1)]],
               constant float& z         [[buffer(2)]],
               constant float& rho       [[buffer(3)]],
               constant float& A         [[buffer(4)]],
               constant int&   num_el    [[buffer(5)]],
               device   int*   seeds     [[buffer(6)]],
               constant float& sig       [[buffer(7)]],
               constant float& omega     [[buffer(8)]],
               device   float* depth     [[buffer(9)]],
               device   float* energy    [[buffer(10)]],
               constant int&   steps     [[buffer(11)]],
               constant float& thickness [[buffer(12)]],
               device   float* LamxSH    [[buffer(13)]],
               device   float* LamySH    [[buffer(14)]],
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
    float E_new, alpha, de_ds, phi, psi, mfp, sig_eNA, step, dsq, dsqi, absc0z, th;
    th = thickness;

    float J;
    J = (9.76f*z + 58.5f*powr(z,-0.19f))*1E-3f;
    float xinit, zinit;
    xinit = 0.0f;
    zinit = 0.0f;

    float4 r0 = float4(xinit, 0.0f, -zinit, 0.0f);
    float4 c0 = float4(cos(omega)*sin(sig), sin(omega)*sin(sig), cos(sig), 0.0f);
    float escape_depth;

    for (int i = 0; i < num_el; ++i){
        LamxSH[num_el*id + i] = -10.0f;
        LamySH[num_el*id + i] = -10.0f;
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
        r0 = float4(xinit, 0.0f, -zinit, 0.0f);
// change the sign of the incidence angle since we're accumulating in the Southern hemisphere.
        c0 = float4(-cos(omega)*sin(sig), -sin(omega)*sin(sig), cos(sig), 0.0f);
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
                escape_depth = (thickness - r_new.z)/c_new.z;
            }

            r_new = r0 + step*c_new*1.0e7f;
            r0 = r_new;
            c0 = c_new;
            E_new += step*rho*de_ds;

// check whether the electron has escaped; depends on the sample (foil) shape.
// BSEs (escaped on the +z / front side) are eliminated; transmitted electrons
// (z >= thickness) are accumulated in the southern hemisphere.
            if ( r0.z <= 0  && counter2 == 0) {
                counter2 = 1;
            }
            if ( r0.z >= thickness && counter2 == 0 ) {
                dir_cos[0] = c0.x;
                dir_cos[1] = c0.y;
                dir_cos[2] = c0.z;
                if (dir_cos[2] < 0.0f){
                    dir_cos[2] = -dir_cos[2];
                }
                if(dir_cos[0] != 0.0f && dir_cos[1] != 0.0f && dir_cos[2] != 0.0f){
                    ret = LambertSphereToPlane(float3(dir_cos[0], dir_cos[1], dir_cos[2]));
                    LamxSH[num_el*id + i] = ret.x;
                    LamySH[num_el*id + i] = ret.y;
                    depth[num_el*id + i]  = escape_depth;
                    energy[num_el*id + i] = E_new;
                    counter2 = 1;
                }
            }

            counter1++ ;
        }
    }
}
