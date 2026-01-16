#include <random>
#include "SFMT.h"

using namespace std;

sfmt_t sfmt;
static uint32_t g_seed;

void SFMTinitialize() {
    sfmt_init_gen_rand(&sfmt, unsigned(time(0)));
}

unsigned int drand() {
    return sfmt_genrand_uint64(&sfmt);
}

double prand() {
    return sfmt_genrand_res53(&sfmt);
}

void fastSrand(){
    g_seed = time(NULL);
}

uint32_t fastRand() {
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}