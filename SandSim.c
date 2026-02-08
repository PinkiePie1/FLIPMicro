#include "SandSim.h"
#include "main.h"

//加快运算速度，一些除法可以先算。
static fixed_t invertSpacing;
static fixed_t particleRestDensity = 0;
static fixed_t halfSpacing;
static fixed_t minPos;
static fixed_t maxXPos;
static fixed_t maxYPos;
static fixed_t maxXBoundary;
static fixed_t maxYBoundary;
const fixed_t nudge = FIXED_FROM_RATIO(1, 100);
const fixed_t flipBlend = FIXED_FROM_RATIO(1, 5);
const fixed_t half = FIXED_FROM_RATIO(1, 2);

//各种数据的缓存
fixed_t particlePos[NumberOfParticles*2]; //粒子的位置，x位置为2*n，y位置为2n+1
fixed_t particleVel[NumberOfParticles*2]; //粒子的速度，x速度为2*n，y速度为2n+1

fixed_t uVel[(CellNumX + 1) * CellNumY];   //u分量 (x方向)
fixed_t vVel[CellNumX * (CellNumY + 1)];   //v分量 (y方向)
fixed_t uPrev[(CellNumX + 1) * CellNumY];  //上一帧u
fixed_t vPrev[CellNumX * (CellNumY + 1)];  //上一帧v
fixed_t uWeights[(CellNumX + 1) * CellNumY];
fixed_t vWeights[CellNumX * (CellNumY + 1)];
unsigned int gridType[CellCount]; //网格是气体还是液体

unsigned int Count[CellCount+1U]; //计算粒子碰撞所用的缓存，其中每个entry对应了应该去哪里找粒子。
fixed_t particleDensity[CellCount];//框内粒子情况，可用此作为显示。
unsigned int particlePosId[NumberOfParticles]; //经过hashgrid排序之后的粒子位置
unsigned int baseGridType[CellCount];
void printLocation(unsigned int n);

static int clamp_index(int value, int min_value, int max_value) {
    if (value < min_value) {
        return min_value;
    }
    if (value > max_value) {
        return max_value;
    }
    return value;
}

static fixed_t clamp_fixed(fixed_t value, fixed_t min_value, fixed_t max_value) {
    if (value < min_value) {
        return min_value;
    }
    if (value > max_value) {
        return max_value;
    }
    return value;
}

static inline int fixed_to_int(fixed_t value) {
    return (int)(value >> FIXED_SHIFT);
}

#define U_INDEX(x,y) ((x) * CellNumY + (y))
#define V_INDEX(x,y) ((x) * (CellNumY + 1) + (y))

static uint32_t isqrt32(uint32_t value) {
    uint32_t op = value;
    uint32_t res = 0;
    uint32_t one = 1U << 30;

    while (one > op) {
        one >>= 2;
    }

    while (one != 0) {
        if (op >= res + one) {
            op -= res + one;
            res = (res >> 1) + one;
        } else {
            res >>= 1;
        }
        one >>= 2;
    }

    return res;
}

static void accumulate_u(fixed_t x, fixed_t y, fixed_t value) {
    fixed_t fx = FIXED_MUL(x, invertSpacing);
    fixed_t fy = FIXED_MUL(y, invertSpacing) - half;
    int x0 = (int)(fx >> FIXED_SHIFT);
    int y0 = (int)(fy >> FIXED_SHIFT);
    x0 = clamp_index(x0, 0, CellNumX - 2);
    y0 = clamp_index(y0, 0, CellNumY - 2);
    fixed_t tx = fx - FIXED_FROM_INT(x0);
    fixed_t ty = fy - FIXED_FROM_INT(y0);
    tx = FIXED_CLAMP(tx, 0, FIXED_ONE);
    ty = FIXED_CLAMP(ty, 0, FIXED_ONE);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    fixed_t w00 = FIXED_MUL((FIXED_ONE - tx), (FIXED_ONE - ty));
    fixed_t w10 = FIXED_MUL(tx, (FIXED_ONE - ty));
    fixed_t w01 = FIXED_MUL((FIXED_ONE - tx), ty);
    fixed_t w11 = FIXED_MUL(tx, ty);

    uVel[U_INDEX(x0, y0)] += FIXED_MUL(value, w00);
    uWeights[U_INDEX(x0, y0)] += w00;
    uVel[U_INDEX(x1, y0)] += FIXED_MUL(value, w10);
    uWeights[U_INDEX(x1, y0)] += w10;
    uVel[U_INDEX(x0, y1)] += FIXED_MUL(value, w01);
    uWeights[U_INDEX(x0, y1)] += w01;
    uVel[U_INDEX(x1, y1)] += FIXED_MUL(value, w11);
    uWeights[U_INDEX(x1, y1)] += w11;
}

static void accumulate_v(fixed_t x, fixed_t y, fixed_t value) {
    fixed_t fx = FIXED_MUL(x, invertSpacing) - half;
    fixed_t fy = FIXED_MUL(y, invertSpacing);
    int x0 = (int)(fx >> FIXED_SHIFT);
    int y0 = (int)(fy >> FIXED_SHIFT);
    x0 = clamp_index(x0, 0, CellNumX - 2);
    y0 = clamp_index(y0, 0, CellNumY - 2);
    fixed_t tx = fx - FIXED_FROM_INT(x0);
    fixed_t ty = fy - FIXED_FROM_INT(y0);
    tx = FIXED_CLAMP(tx, 0, FIXED_ONE);
    ty = FIXED_CLAMP(ty, 0, FIXED_ONE);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    fixed_t w00 = FIXED_MUL((FIXED_ONE - tx), (FIXED_ONE - ty));
    fixed_t w10 = FIXED_MUL(tx, (FIXED_ONE - ty));
    fixed_t w01 = FIXED_MUL((FIXED_ONE - tx), ty);
    fixed_t w11 = FIXED_MUL(tx, ty);

    vVel[V_INDEX(x0, y0)] += FIXED_MUL(value, w00);
    vWeights[V_INDEX(x0, y0)] += w00;
    vVel[V_INDEX(x1, y0)] += FIXED_MUL(value, w10);
    vWeights[V_INDEX(x1, y0)] += w10;
    vVel[V_INDEX(x0, y1)] += FIXED_MUL(value, w01);
    vWeights[V_INDEX(x0, y1)] += w01;
    vVel[V_INDEX(x1, y1)] += FIXED_MUL(value, w11);
    vWeights[V_INDEX(x1, y1)] += w11;
}

static fixed_t sample_u_filtered(fixed_t x, fixed_t y, const fixed_t *grid, fixed_t *weight_out) {
    fixed_t fx = FIXED_MUL(x, invertSpacing);
    fixed_t fy = FIXED_MUL(y, invertSpacing) - half;
    int x0 = (int)(fx >> FIXED_SHIFT);
    int y0 = (int)(fy >> FIXED_SHIFT);
    x0 = clamp_index(x0, 0, CellNumX - 2);
    y0 = clamp_index(y0, 0, CellNumY - 2);
    fixed_t tx = fx - FIXED_FROM_INT(x0);
    fixed_t ty = fy - FIXED_FROM_INT(y0);
    tx = FIXED_CLAMP(tx, 0, FIXED_ONE);
    ty = FIXED_CLAMP(ty, 0, FIXED_ONE);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    fixed_t w00 = FIXED_MUL((FIXED_ONE - tx), (FIXED_ONE - ty));
    fixed_t w10 = FIXED_MUL(tx, (FIXED_ONE - ty));
    fixed_t w01 = FIXED_MUL((FIXED_ONE - tx), ty);
    fixed_t w11 = FIXED_MUL(tx, ty);

    int n = CellNumY;
    int nr0 = x0 * n + y0;
    int nr1 = x1 * n + y0;
    int nr2 = x1 * n + y1;
    int nr3 = x0 * n + y1;
    int offset = n;

    fixed_t valid0 = (gridType[nr0] != AIR_CELL || gridType[nr0 - offset] != AIR_CELL) ? FIXED_ONE : 0;
    fixed_t valid1 = (gridType[nr1] != AIR_CELL || gridType[nr1 - offset] != AIR_CELL) ? FIXED_ONE : 0;
    fixed_t valid2 = (gridType[nr2] != AIR_CELL || gridType[nr2 - offset] != AIR_CELL) ? FIXED_ONE : 0;
    fixed_t valid3 = (gridType[nr3] != AIR_CELL || gridType[nr3 - offset] != AIR_CELL) ? FIXED_ONE : 0;

    fixed_t d = FIXED_MUL(valid0, w00) + FIXED_MUL(valid1, w10) +
        FIXED_MUL(valid2, w11) + FIXED_MUL(valid3, w01);
    if (weight_out != NULL) {
        *weight_out = d;
    }
    if (d == 0) {
        return 0;
    }

    fixed_t num = 0;
    if (valid0) {
        num += FIXED_MUL(grid[U_INDEX(x0, y0)], w00);
    }
    if (valid1) {
        num += FIXED_MUL(grid[U_INDEX(x1, y0)], w10);
    }
    if (valid2) {
        num += FIXED_MUL(grid[U_INDEX(x1, y1)], w11);
    }
    if (valid3) {
        num += FIXED_MUL(grid[U_INDEX(x0, y1)], w01);
    }

    return FIXED_DIV(num, d);
}

static fixed_t sample_v_filtered(fixed_t x, fixed_t y, const fixed_t *grid, fixed_t *weight_out) {
    fixed_t fx = FIXED_MUL(x, invertSpacing) - half;
    fixed_t fy = FIXED_MUL(y, invertSpacing);
    int x0 = (int)(fx >> FIXED_SHIFT);
    int y0 = (int)(fy >> FIXED_SHIFT);
    x0 = clamp_index(x0, 0, CellNumX - 2);
    y0 = clamp_index(y0, 0, CellNumY - 2);
    fixed_t tx = fx - FIXED_FROM_INT(x0);
    fixed_t ty = fy - FIXED_FROM_INT(y0);
    tx = FIXED_CLAMP(tx, 0, FIXED_ONE);
    ty = FIXED_CLAMP(ty, 0, FIXED_ONE);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    fixed_t w00 = FIXED_MUL((FIXED_ONE - tx), (FIXED_ONE - ty));
    fixed_t w10 = FIXED_MUL(tx, (FIXED_ONE - ty));
    fixed_t w01 = FIXED_MUL((FIXED_ONE - tx), ty);
    fixed_t w11 = FIXED_MUL(tx, ty);

    int n = CellNumY;
    int nr0 = x0 * n + y0;
    int nr1 = x1 * n + y0;
    int nr2 = x1 * n + y1;
    int nr3 = x0 * n + y1;
    int offset = 1;

    fixed_t valid0 = (gridType[nr0] != AIR_CELL || gridType[nr0 - offset] != AIR_CELL) ? FIXED_ONE : 0;
    fixed_t valid1 = (gridType[nr1] != AIR_CELL || gridType[nr1 - offset] != AIR_CELL) ? FIXED_ONE : 0;
    fixed_t valid2 = (gridType[nr2] != AIR_CELL || gridType[nr2 - offset] != AIR_CELL) ? FIXED_ONE : 0;
    fixed_t valid3 = (gridType[nr3] != AIR_CELL || gridType[nr3 - offset] != AIR_CELL) ? FIXED_ONE : 0;

    fixed_t d = FIXED_MUL(valid0, w00) + FIXED_MUL(valid1, w10) +
        FIXED_MUL(valid2, w11) + FIXED_MUL(valid3, w01);
    if (weight_out != NULL) {
        *weight_out = d;
    }
    if (d == 0) {
        return 0;
    }

    fixed_t num = 0;
    if (valid0) {
        num += FIXED_MUL(grid[V_INDEX(x0, y0)], w00);
    }
    if (valid1) {
        num += FIXED_MUL(grid[V_INDEX(x1, y0)], w10);
    }
    if (valid2) {
        num += FIXED_MUL(grid[V_INDEX(x1, y1)], w11);
    }
    if (valid3) {
        num += FIXED_MUL(grid[V_INDEX(x0, y1)], w01);
    }

    return FIXED_DIV(num, d);
}

void InitParticles(){
    //在每个格子里生成一个。如果全填满了则不再生成。
    int p_num = 0;
    invertSpacing = FIXED_DIV(FIXED_ONE, Spacing);
    halfSpacing = FIXED_MUL(Spacing, half);
    minPos = Spacing;
    maxXPos = fixed_mul_int(Spacing, CellNumX - 1);
    maxYPos = fixed_mul_int(Spacing, CellNumY - 1);
    maxXBoundary = fixed_mul_int(Spacing, CellNumX) - ParticleRadius;
    maxYBoundary = fixed_mul_int(Spacing, CellNumY) - ParticleRadius;

    for (int x = 0; x < CellNumX; x++) {
        for (int y = 0; y < CellNumY; y++) {
            if (x == 0 || y == 0 || x == CellNumX - 1 || y == CellNumY - 1) {
                baseGridType[INDEX(x, y)] = SOLID_CELL;
            } else {
                baseGridType[INDEX(x, y)] = AIR_CELL;
            }
        }
    }

    for (int i = 1; i < CellNumY ; i++){
        for(int j = 1; j < CellNumX; j++){

            int particleIndex = p_num;
            particlePos[XID(particleIndex)] = fixed_mul_int(Spacing, j) + halfSpacing;
            particlePos[YID(particleIndex)] = fixed_mul_int(Spacing, i) + halfSpacing;
            particleVel[XID(particleIndex)] = 0;
            particleVel[YID(particleIndex)] = 0;
            if(p_num++ >= NumberOfParticles){return;} 
            PRINT("generated a particle %d at %ld, %ld. \n",particleIndex,(long)particlePos[XID(particleIndex)],(long)particlePos[YID(particleIndex)]);
            
        }
    }

    return;
}

/*
 * @brief 根据所提供的加速度仿真粒子的行动
 *
 * @param xAccleration,yAcceleration xy方向的加速度
 * @return none
 *
 * @note 也会仿真墙壁。 
 */ 
void ParticleIntegrate(fixed_t xAcceleration, fixed_t yAcceleration){
    for(unsigned int i = 0; i < NumberOfParticles ; i++){
        //先算速度
        particleVel[XID(i)] += FIXED_MUL(xAcceleration, dt);
        particleVel[YID(i)] += FIXED_MUL(yAcceleration, dt);
        //再算位置
        particlePos[XID(i)] += FIXED_MUL(particleVel[XID(i)], dt);
        particlePos[YID(i)] += FIXED_MUL(particleVel[YID(i)], dt); 

        // 边界。因为网格数据结构的原因，第一行和第一列是没有速度量的，所以得把粒子挤出去，第一行和第一列不放粒子。
        // X方向边界
        if(particlePos[XID(i)] < minPos){
            particlePos[XID(i)] = minPos + nudge;
            particleVel[XID(i)] = FIXED_MUL(particleVel[XID(i)], BOUNCYNESS); // 反弹阻尼
        }
        if(particlePos[XID(i)] >= maxXBoundary){
            particlePos[XID(i)] = maxXBoundary - nudge;
            particleVel[XID(i)] = FIXED_MUL(particleVel[XID(i)], BOUNCYNESS);
        }
        
        // Y方向边界
        if(particlePos[YID(i)] < minPos){
            particlePos[YID(i)] = minPos + nudge;
            particleVel[YID(i)] = FIXED_MUL(particleVel[YID(i)], BOUNCYNESS);
        }
        if(particlePos[YID(i)] >= maxYBoundary){
            particlePos[YID(i)] = maxYBoundary - nudge;
            particleVel[YID(i)] = FIXED_MUL(particleVel[YID(i)], BOUNCYNESS);
        }

    }
    return;
}

/*
 * @brief 相互推开粒子，模拟粒子之间的碰撞
 *
 * @param nIters 仿真重复的次数
 * @return none
 *
 * @note 也会仿真墙壁。 
 */ 
void PushParticlesApart(unsigned int nIters){

    //清零缓存
    memset(Count,0,sizeof(Count));
    memset(particlePosId,0,sizeof(particlePosId));

    //数数
    for (unsigned int i=0;i<NumberOfParticles;i++){
        fixed_t x = particlePos[XID(i)];
        fixed_t y = particlePos[YID(i)];

        unsigned int xi = (unsigned int)fixed_to_int(FIXED_MUL(x, invertSpacing)); //x和y的网格坐标，即这个粒子在网格中所处的位置
        unsigned int yi = (unsigned int)fixed_to_int(FIXED_MUL(y, invertSpacing)); 

        Count[INDEX(xi,yi)] ++;
    }

    unsigned int first = 0;
    
    //partial sum
    for (unsigned int i=0;i<CellCount;i++){
        first += Count[i];
        Count[i] = first;
    } 
    Count[CellCount] = first; //guard
    
    //将排序好的粒子塞入对应的位置缓存中
    for (unsigned int i=0;i<NumberOfParticles;i++){
        fixed_t x = particlePos[XID(i)];//取出xy坐标
        fixed_t y = particlePos[YID(i)];

        unsigned int xi = (unsigned int)fixed_to_int(FIXED_MUL(x, invertSpacing)); //x和y的网格坐标，即这个粒子在网格中所处的位置
        unsigned int yi = (unsigned int)fixed_to_int(FIXED_MUL(y, invertSpacing)); 

        int gridindex = --Count[INDEX(xi,yi)]; //对应位置先减1，再放入排序好了的目标缓存
        particlePosId[gridindex] = i; //放入排序的缓存
        PRINT("index is %d %d",xi,yi);
    }
    
//开始把粒子推开

fixed_t minDist = ParticleRadius * 2;
fixed_t minDist_scaled = minDist >> 8;
uint32_t minDist2 = (uint32_t)minDist_scaled * (uint32_t)minDist_scaled;
do {
    for (unsigned int i = 0; i < NumberOfParticles; i++) {
        fixed_t px = particlePos[XID(i)];
        fixed_t py = particlePos[YID(i)];

        unsigned int pxi = (unsigned int)fixed_to_int(FIXED_MUL(px, invertSpacing));
        unsigned int pyi = (unsigned int)fixed_to_int(FIXED_MUL(py, invertSpacing));
        unsigned int x0 = pxi > 0 ? pxi - 1 : 0;
        unsigned int y0 = pyi > 0 ? pyi - 1 : 0;
        unsigned int x1 = (pxi + 1 < CellNumX) ? pxi + 1 : CellNumX - 1;
        unsigned int y1 = (pyi + 1 < CellNumY) ? pyi + 1 : CellNumY - 1;

        for (int xi = x0; xi <= x1; xi++) {
            for (int yi = y0; yi <= y1; yi++) {
                unsigned int cellNr = INDEX(xi,yi);
                unsigned int first = Count[cellNr];
                unsigned int last = Count[cellNr + 1];
                for (int j = first; j < last; j++) { //遍历xi，yi格子中的所有粒子，并将其和第i个粒子对比，看是否有碰撞。
                    int id = particlePosId[j];
                    if (id == i)
                        continue;
                    fixed_t qx = particlePos[XID(id)];
                    fixed_t qy = particlePos[YID(id)];

                    fixed_t dx = qx - px;
                    fixed_t dy = qy - py;
                    int32_t dx_scaled = dx >> 8;
                    int32_t dy_scaled = dy >> 8;
                    uint32_t d2 = (uint32_t)(dx_scaled * dx_scaled) +
                        (uint32_t)(dy_scaled * dy_scaled);
                    if (d2 > minDist2)
                        continue;
                    else if(d2 == 0){
                        PRINT("kiss detected on %d,%d",i,id);
                        continue;
                    }
                    PRINT("bounce detected on %d,%d",i,id);
                    fixed_t d = (fixed_t)isqrt32(d2) << 8;
                    fixed_t s = FIXED_DIV((minDist - d) >> 1, d);
                    dx = FIXED_MUL(dx, s);
                    dy = FIXED_MUL(dy, s);
                    particlePos[XID(i)] -= dx;
                    particlePos[YID(i)] -= dy;
                    particlePos[XID(id)] += dx;
                    particlePos[YID(id)] += dy;
                }

        }
    }
    

    }

    nIters--; 
    } while ( nIters > 0 );
    
    //推完了，再处理一次边界
    for(int i = 0; i < NumberOfParticles; i++){
        if(particlePos[XID(i)] < minPos){
            particlePos[XID(i)] = minPos + nudge;
        }
        if(particlePos[XID(i)] >= maxXBoundary){
            particlePos[XID(i)] = maxXBoundary - nudge;
        }
        
        // Y方向边界
        if(particlePos[YID(i)] < minPos){
            particlePos[YID(i)] = minPos + nudge;
        }
        if(particlePos[YID(i)] >= maxYBoundary){
            particlePos[YID(i)] = maxYBoundary - nudge;
        }
    }
    return;
}

/*
 * @brief 寻找有粒子的框。
 *
 * @param none
 * @return none
 *
 * @note 
 */ 
 void density_update(){

    memset(particleDensity,0,sizeof(particleDensity));
    for (unsigned int i=0;i<NumberOfParticles;i++){
        fixed_t x = clamp_fixed(particlePos[XID(i)], minPos, maxXPos);
        fixed_t y = clamp_fixed(particlePos[YID(i)], minPos, maxYPos);

        fixed_t x_shifted = x - halfSpacing;
        fixed_t y_shifted = y - halfSpacing;

        int x0 = fixed_to_int(FIXED_MUL(x_shifted, invertSpacing));
        int y0 = fixed_to_int(FIXED_MUL(y_shifted, invertSpacing));
        x0 = clamp_index(x0, 0, CellNumX - 2);
        y0 = clamp_index(y0, 0, CellNumY - 2);
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        fixed_t tx = FIXED_MUL(x_shifted - fixed_mul_int(Spacing, x0), invertSpacing);
        fixed_t ty = FIXED_MUL(y_shifted - fixed_mul_int(Spacing, y0), invertSpacing);
        tx = FIXED_CLAMP(tx, 0, FIXED_ONE);
        ty = FIXED_CLAMP(ty, 0, FIXED_ONE);

        fixed_t sx = FIXED_ONE - tx;
        fixed_t sy = FIXED_ONE - ty;

        particleDensity[INDEX(x0,y0)] += FIXED_MUL(sx, sy);
        particleDensity[INDEX(x1,y0)] += FIXED_MUL(tx, sy);
        particleDensity[INDEX(x1,y1)] += FIXED_MUL(tx, ty);
        particleDensity[INDEX(x0,y1)] += FIXED_MUL(sx, ty);
    }

    if (particleRestDensity == 0) {
        fixed_t sum = 0;
        int numFluidCells = 0;

        for (int i = 0; i < CellCount; i++) {
            if (gridType[i] == FLUID_CELL) {
                sum += particleDensity[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0) {
            particleRestDensity = FIXED_DIV(sum, FIXED_FROM_INT(numFluidCells));
        }
    }

    return;
 }

 /*
 * @brief 粒子速度换算到网格。
 *
 * @param none
 * @return none
 *
 * @note 
 */ 
void particles_to_grid() {
    //保存上一帧投影后的网格速度供FLIP使用
    memcpy(uPrev, uVel, sizeof(uVel));
    memcpy(vPrev, vVel, sizeof(vVel));

    //清零密度，网格类型和网格速度
    memcpy(gridType, baseGridType, sizeof(gridType));
    memset(uVel,0,sizeof(uVel));
    memset(vVel,0,sizeof(vVel));
    memset(uWeights,0,sizeof(uWeights));
    memset(vWeights,0,sizeof(vWeights));

    for(int p=0; p<NumberOfParticles; p++){

        fixed_t x = particlePos[XID(p)];
        fixed_t y = particlePos[YID(p)];
        //计算粒子对应的网格位置和delta x y
        int xcell = fixed_to_int(FIXED_MUL(x, invertSpacing));
        int ycell = fixed_to_int(FIXED_MUL(y, invertSpacing));
        xcell = clamp_index(xcell, 0, CellNumX - 1);
        ycell = clamp_index(ycell, 0, CellNumY - 1);

        PRINT("calculate particle %d at %ld,%ld. xcell: %d, ycell: %d.\n",p,(long)x,(long)y,xcell,ycell);
        
        int gridIndex = INDEX(xcell,ycell);//对应的网格在缓存中的位置。
        if (gridType[gridIndex] == AIR_CELL) {
            gridType[gridIndex] = FLUID_CELL;//对应网格是液体网格
        }

        accumulate_u(x, y, particleVel[XID(p)]);
        accumulate_v(x, y, particleVel[YID(p)]);

    }
    //传递完了之后应该除掉总权重。
    for (int i=0;i<CellNumX + 1;i++){
       for (int j=0;j<CellNumY;j++){
        if(uWeights[U_INDEX(i,j)]){ //显然如果这里是空气那就不用算了，不应该处以0
            uVel[U_INDEX(i,j)] = FIXED_DIV(uVel[U_INDEX(i,j)], uWeights[U_INDEX(i,j)]);
         }
       }
    }
    for (int i=0;i<CellNumX;i++){
       for (int j=0;j<CellNumY + 1;j++){
        if(vWeights[V_INDEX(i,j)]){ //显然如果这里是空气那就不用算了，不应该处以0
            vVel[V_INDEX(i,j)] = FIXED_DIV(vVel[V_INDEX(i,j)], vWeights[V_INDEX(i,j)]);
        }
       }
    }
    for (int y = 0; y < CellNumY; y++) {
        uVel[U_INDEX(0, y)] = 0;
        uVel[U_INDEX(CellNumX, y)] = 0;
    }
    for (int x = 0; x < CellNumX; x++) {
        vVel[V_INDEX(x, 0)] = 0;
        vVel[V_INDEX(x, CellNumY)] = 0;
    }
    return;
}

 /*
 * @brief 计算不可压缩性
 *
 * @param none
 * @return none
 *
 * @note 
 */ 
void compute_grid_forces(unsigned int nIters) {
    //根据网格速度，计算无法压缩的液体,nIters为迭代数量，需要算n次。
    do{
        //对每个grid进行迭代
        for(int xcell = 1; xcell < CellNumX - 1; xcell++){
            for(int ycell = 1; ycell < CellNumY - 1; ycell++){
                if(gridType[INDEX(xcell,ycell)] != FLUID_CELL){
                    continue;//如果这个grid是空气，就跳过。
                } else {
                //如果grid是液体，可以计算divergence。首先求出这个grid的xy坐标。

                //计算divergence
                fixed_t d = uVel[U_INDEX(xcell + 1, ycell)]
                    - uVel[U_INDEX(xcell, ycell)]
                    + vVel[V_INDEX(xcell, ycell + 1)]
                    - vVel[V_INDEX(xcell, ycell)];//计算divergence
                if (particleRestDensity > 0) {
                    fixed_t compression = particleDensity[INDEX(xcell,ycell)] - particleRestDensity;
                    if (compression > 0) {
                        compression = FIXED_MUL(compression, stiffnessCoefficient);
                        d -= compression;
                    }
                }

                int s1 = gridType[INDEX(xcell - 1, ycell)] != SOLID_CELL ? 1 : 0;
                int s2 = gridType[INDEX(xcell, ycell - 1)] != SOLID_CELL ? 1 : 0;
                int s3 = gridType[INDEX(xcell + 1, ycell)] != SOLID_CELL ? 1 : 0;
                int s4 = gridType[INDEX(xcell, ycell + 1)] != SOLID_CELL ? 1 : 0;

                int s_sum = s1 + s2 + s3 + s4;
                if (s_sum == 0) {
                    continue;
                }

                fixed_t p = FIXED_DIV(-d, FIXED_FROM_INT(s_sum));
                p = FIXED_MUL(p, overRelaxiation);

                if (s1) {
                    uVel[U_INDEX(xcell, ycell)] -= p;
                }
                if (s3) {
                    uVel[U_INDEX(xcell + 1, ycell)] += p;
                }
                if (s2) {
                    vVel[V_INDEX(xcell, ycell)] -= p;
                }
                if (s4) {
                    vVel[V_INDEX(xcell, ycell + 1)] += p;
                }



            }
            }
        }
        nIters--;       
    } while (nIters>0);

    //墙壁边界条件：边界速度为0
    for (int y = 0; y < CellNumY; y++) {
        uVel[U_INDEX(0, y)] = 0;
        uVel[U_INDEX(CellNumX, y)] = 0;
    }
    for (int x = 0; x < CellNumX; x++) {
        vVel[V_INDEX(x, 0)] = 0;
        vVel[V_INDEX(x, CellNumY)] = 0;
    }
    return;
}
 /*
 * @brief 把网格速度还给粒子
 *
 * @param none
 * @return none
 *
 * @note 
 */ 

void grid_to_particles() {
    for(int p=0; p<NumberOfParticles; p++){
        fixed_t x = clamp_fixed(particlePos[XID(p)], minPos, maxXPos);
        fixed_t y = clamp_fixed(particlePos[YID(p)], minPos, maxYPos);

        fixed_t weight = 0;
        fixed_t pic_vx = sample_u_filtered(x, y, uVel, &weight);
        if (weight > 0) {
            fixed_t prev_vx = sample_u_filtered(x, y, uPrev, NULL);
            fixed_t flip_vx = particleVel[XID(p)] + (pic_vx - prev_vx);
            particleVel[XID(p)] = FIXED_MUL(pic_vx, FIXED_ONE - flipBlend) + FIXED_MUL(flip_vx, flipBlend);
        }

        weight = 0;
        fixed_t pic_vy = sample_v_filtered(x, y, vVel, &weight);
        if (weight > 0) {
            fixed_t prev_vy = sample_v_filtered(x, y, vPrev, NULL);
            fixed_t flip_vy = particleVel[YID(p)] + (pic_vy - prev_vy);
            particleVel[YID(p)] = FIXED_MUL(pic_vy, FIXED_ONE - flipBlend) + FIXED_MUL(flip_vy, flipBlend);
        }
    }

    return;
}

void visualize_grid() {
    //显示
    // 初始化缓冲区为全'-'
    char visual_buffer[CellNumY][CellNumX+1];
    memset(visual_buffer, '-', sizeof(visual_buffer));

    printf("\e[1;1H\e[2J");//清空屏幕；。
    
    // 标记粒子位置
    for(int p=0; p<NumberOfParticles; p++){
        int x = fixed_to_int(FIXED_MUL(particlePos[XID(p)], invertSpacing));
        int y = fixed_to_int(FIXED_MUL(particlePos[YID(p)], invertSpacing));

        visual_buffer[y][x] = 'x';
    }

    // 添加字符串终止符
    for(int j=0; j<CellNumY; j++){
        visual_buffer[j][CellNumX] = '\0';
    }

    
    printf("PIC Simulation (X: %d, Y: %d, Particles: %d)\n", 
           CellNumX, CellNumY, NumberOfParticles);
    for(int j=0; j<CellNumY; j++){
        printf("%s\n", visual_buffer[j]);
    }
    printLocation(0);
    printLocation(1);
    fflush(stdout); 
    
}

static void print_fixed_value(fixed_t value) {
    if (value < 0) {
        printf("-");
        value = -value;
    }
    int32_t integer = value >> FIXED_SHIFT;
    uint32_t frac = (uint32_t)(value & (FIXED_ONE - 1));
    uint32_t frac_scaled = (frac * 100U) >> FIXED_SHIFT;
    printf("%ld.%02u", (long)integer, frac_scaled);
}

void printLocation(unsigned int n){
    printf("Particle %d location:",n);
    print_fixed_value(particlePos[XID(n)]);
    printf(",");
    print_fixed_value(particlePos[YID(n)]);
    printf(", speed is ");
    print_fixed_value(particleVel[XID(n)]);
    printf(",");
    print_fixed_value(particleVel[YID(n)]);
    printf("\n");
}
