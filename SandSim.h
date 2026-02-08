#include <stdint.h>
#include <stdio.h>
#include <string.h>

// fixed-point configuration (Q16.16, 32-bit only)
typedef int32_t fixed_t;
#define FIXED_SHIFT 16
#define FIXED_ONE ((fixed_t)1 << FIXED_SHIFT)
#define FIXED_FROM_INT(x) ((fixed_t)((x) << FIXED_SHIFT))
#define FIXED_FROM_RATIO(num, den) ((fixed_t)((((fixed_t)(num)) * FIXED_ONE) / (den)))

static inline fixed_t fixed_mul(fixed_t a, fixed_t b) {
    int32_t a_hi = a >> 16;
    int32_t a_lo = a & 0xFFFF;
    int32_t b_hi = b >> 16;
    int32_t b_lo = b & 0xFFFF;

    int32_t hi = a_hi * b_hi;
    int32_t mid = a_hi * b_lo + a_lo * b_hi;
    uint32_t lo = (uint32_t)a_lo * (uint32_t)b_lo;

    return (hi << 16) + mid + (int32_t)(lo >> 16);
}

static inline fixed_t fixed_div(fixed_t numerator, fixed_t denominator) {
    if (denominator == 0) {
        return 0;
    }
    int sign = 1;
    uint32_t num = (uint32_t)numerator;
    uint32_t den = (uint32_t)denominator;
    if (numerator < 0) {
        num = (uint32_t)(-numerator);
        sign = -sign;
    }
    if (denominator < 0) {
        den = (uint32_t)(-denominator);
        sign = -sign;
    }

    uint32_t integer_part = num / den;
    uint32_t remainder = num % den;
    uint32_t quotient = integer_part << FIXED_SHIFT;

    for (int i = 0; i < FIXED_SHIFT; i++) {
        remainder <<= 1;
        if (remainder >= den) {
            remainder -= den;
            quotient |= (1U << (FIXED_SHIFT - 1 - i));
        }
    }

    int32_t result = (int32_t)quotient;
    return sign < 0 ? -result : result;
}

static inline fixed_t fixed_mul_int(fixed_t value, int32_t multiplier) {
    return (fixed_t)(value * multiplier);
}

#define FIXED_MUL(a, b) fixed_mul((a), (b))
#define FIXED_DIV(a, b) fixed_div((a), (b))
#define FIXED_CLAMP(x, lo, hi) ((x) < (lo) ? (lo) : ((x) > (hi) ? (hi) : (x)))

//模拟所用的参数
#define NumberOfParticles 128U //粒子数量
#define ParticleRadius FIXED_FROM_RATIO(5, 100) //粒子半径，这同时也是网格的间距
#define Spacing FIXED_FROM_RATIO(1, 10) //网格间距
#define CellNumX 16U //x轴方向的网格数量
#define CellNumY 16U //y轴方向的网格数量
#define CellCount CellNumX*CellNumY //网格总数
#define dt FIXED_FROM_RATIO(1, 100) //时间步长
#define BOUNCYNESS FIXED_FROM_RATIO(-1, 5) //墙壁的弹性。

#define overRelaxiation FIXED_FROM_RATIO(19, 10)
#define stiffnessCoefficient FIXED_FROM_RATIO(3, 5)

//用于访问粒子位置的一些方便的宏
#define XID(n) 2*(n)    //访问x位置
#define YID(n) 2*(n)+1  //访问y位置
#define INDEX(x,y) ((x)*CellNumY+y) //根据网格的xy坐标访问对应的网格



void ParticleIntegrate(fixed_t xAcceleration, fixed_t yAcceleration);  //根据所提供的加速度仿真粒子的行动，注意墙壁碰撞应该在这里进行。
void PushParticlesApart(unsigned int nIters);                     //将粒子相互推开

void density_update(void);
void particles_to_grid(void);
void compute_grid_forces(unsigned int nIters);
void grid_to_particles(void);

void InitParticles(void);//初始化粒子，还没完成
void visualize_grid(void); //显示
