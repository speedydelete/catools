/*
Copyright 2025 speedydelete

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

searches for orthogonal NRSS (nonadjustible reduced speed of ship) that face east in square-grid range-1 Moore neighborhood cellular automata (without B0)
see https://conwaylife.com/forums/viewtopic.php?f=11&t=6352&p=218310 for more informatio
to compile: gcc -Wall -Werror -Ofast -o nrss nrss.c
to use: nrss <engine-count> <max-x-seperation> <max-period> <randomize-soups-1-or-0> <state-file>
when randomization is off it will try every possible combination of engines
*/

#include <stdbool.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <time.h>

// parameters

#define RULESTR "B2-ak3ce4eikqrz5-iknq6-ek8/S1c2aek3aekn4eiknry5eiky6-ei7c8"

// maximum height and width, these are the base-2 logarithms
// should not be higher than 16
#define HEIGHT 8
#define WIDTH 12

// don't change
typedef uint_fast16_t uint16;
typedef uint_fast32_t uint32;
typedef uint_fast64_t uint64;
#define HEIGHTVALUE (1 << HEIGHT)
#define WIDTHVALUE (1 << WIDTH)
#define SIZE (HEIGHT + WIDTH)
#define SIZEVALUE (1 << SIZE)

// starting position
#define STARTX 64
#define STARTY ((1 << HEIGHT) / 2 - 64)

// min and max y seperation between engines
#define MINY 7
#define MAXY 12

// put an engine in the data array at the specified position
// i is the array index of the (x, y) coordinates of where to put the engine
// the default one is 2o$o$2o!
static inline void put_engine(uint8_t data[], uint32 i) {
    data[i] = 1;
    data[i + 1] = 1;
    data[i + WIDTHVALUE] = 1;
    data[i + WIDTHVALUE + WIDTHVALUE] = 1;
    data[i + WIDTHVALUE + WIDTHVALUE + 1] = 1;
}

// engine dimensions
#define ENGINEHEIGHT 3
#define ENGINEWIDTH 2

// the number of phases of the engine
#define ENGINEPHASES 128

// whether to skip oscillators
#define SKIPOSCILLATORS 1

// whether to skip duplicate speeds
#define SKIPDUPLICATES 1

// maximum number of ships
#define MAXSHIPS 4096

// generations between phase checks
#define CHECKINTERVAL 64

// whether to reduce the period to lowest terms
#define REDUCEPERIOD 1

// minimum period
#define MINPERIOD 3

// uncomment to make it work in stupid online C compilers
// #define BRUH

/*
hardcoded transition table
indexed by 0b(abcdefghi) where the neighborhood is:
adg
beh
cfi
to make a new one, go to https://speedydelete.com/int_tools, open up the console, and put this in:
(() => {
let trs = parseRule('B2-ak3ce4eikqrz5-iknq6-ek8/S1c2aek3aekn4eiknry5eiky6-ei7c8');
let out = new Uint8Array(512);
for (let i = 0; i < 512; i++) {
    let j = (i & 273) | ((i & 32) << 2) | ((i & 4) << 4) | ((i & 128) >> 2) | ((i & 2) << 2) | ((i & 64) >> 4) | ((i & 8) >> 2);
    out[i] = trs[j];
}
return '{' + out.join(', ') + '}';
})()
*/
const uint8_t transitions[512] = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1};

// DEBUG 1 logs every generation and bigger step
// DEBUG 2 logs the state of the pattern as well
// DEBUG 3 logs like everything
#define DEBUG 0

// end parameters


#ifndef BRUH
#include <unistd.h>
#include <fcntl.h>
#endif

#if DEBUG > 2
#define malloc(size) ({ \
    printf("Allocating %lu bytes", (size)); \
    void* out = malloc((size)); \
    printf(": %p\n", out); \
    out; \
})
#define free(ptr) ({ \
    printf("Freeing %p\n", (ptr)); \
    free((ptr)); \
})
#endif

uint32 engines;
uint16 max_x_sep;
uint16 max_period;
bool use_random_soups;
char* state_file;


uint8_t data[SIZEVALUE] = {0};
uint16 top = 0;
uint16 bottom = 0;
uint16 left = 0;
uint16 right = 0;

void clear() {
    uint32 stop = ((uint32)bottom << WIDTH) + left;
    uint32 max = ((uint32)top << WIDTH) + right;
    for (uint32 start = ((uint32)top << WIDTH) + left; start < stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            data[i] = 0;
        }
        max += WIDTHVALUE;
    }
}


uint8_t temp_data[SIZEVALUE];

static inline bool run_row(uint32 i, uint16* lowX, uint16* highX) {
    uint16 tr = 0;
    tr |= (uint16)data[i - WIDTHVALUE - 1] << 5;
    tr |= (uint16)data[i - 1] << 4;
    tr |= (uint16)data[i + WIDTHVALUE - 1] << 3;
    tr |= (uint16)data[i - WIDTHVALUE] << 2;
    tr |= (uint16)data[i] << 1;
    tr |= (uint16)data[i + WIDTHVALUE];
    uint32 max = i - left + right + 1;
    uint8_t value;
    bool any_changed = false;
    uint16 x = left - 1;
    for (; i <= max; i++) {
        tr = (tr << 3) & 511;
        tr |= (uint16)data[i - WIDTHVALUE + 1] << 2;
        tr |= (uint16)data[i + 1] << 1;
        tr |= (uint16)data[i + WIDTHVALUE + 1];
        #if DEBUG > 2
        printf("transition: %"PRIu8" %"PRIuFAST16" %"PRIu8"\n", data[i], tr, transitions[tr]);
        #endif
        value = transitions[tr];
        if (value) {
            any_changed = true;
            if (x < *lowX) {
                *lowX = x;
            }
            if (x > *highX) {
                *highX = x;
            }
        }
        temp_data[i] = value;
        x++;
    }
    return any_changed;
}

bool run_generation() {
    uint16 lowX = WIDTHVALUE;
    uint16 highX = 0;
    uint16 lowY = HEIGHTVALUE;
    uint16 highY = 0;
    bool any_changed;
    uint32 i = ((top - 1) << WIDTH) + left - 1;
    for (uint16 y = top - 1; y <= bottom; y++) {
        #if DEBUG > 2
        printf("i: %"PRIuFAST32"\n", i);
        #endif
        any_changed = run_row(i, &lowX, &highX);
        if (any_changed) {
            if (y < lowY) {
                lowY = y;
            }
            if (y > highY) {
                highY = y;
            }
        }
        i += WIDTHVALUE;
    }
    // idk why this is needed, someone should figure out why because this probably slows it down
    clear();
    top = lowY;
    bottom = highY + 1;
    left = lowX;
    right = highX + 1;
    uint32 stop = ((uint32)bottom << WIDTH) + left;
    uint32 max = ((uint32)top << WIDTH) + right;
    for (uint32 start = (top << WIDTH) + left; start < stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            data[i] = temp_data[i];
        }
        max += WIDTHVALUE;
    }
    return lowX < highX;
}


typedef struct engine_phase {
    uint32 height;
    uint32 width;
    uint8_t data[];
} engine_phase;

engine_phase* engine_phases[ENGINEPHASES];

void generate_phases() {
    clear();
    put_engine(data, (STARTY << WIDTH) + STARTX);
    top = STARTY;
    bottom = STARTY + ENGINEHEIGHT;
    left = STARTX;
    right = STARTX + ENGINEWIDTH;
    for (uint16 i = 0; i < ENGINEPHASES; i++) {
        #if DEBUG > 0
        printf("Generating phase %"PRIuFAST16"\n", i);
        #endif
        uint16 height = bottom - top;
        uint16 width = right - left;
        uint32 size = height * width;
        // printf("Size: %"PRIuFAST32", Total size: %lu\n", size, sizeof(engine_phase) + size);
        engine_phase* phase = malloc(sizeof(engine_phase) + size);
        phase->height = height;
        phase->width = width;
        uint32 j = 0;
        uint32 stop = ((uint32)bottom << WIDTH) + left;
        uint32 max = ((uint32)top << WIDTH) + right;
        for (uint32 start = (top << WIDTH) + left; start < stop; start += WIDTHVALUE) {
            for (uint32 i = start; i < max; i++) {
                // printf("%"PRIuFAST32"\n", j);
                phase->data[j] = data[i];
                j++;
            }
            max += WIDTHVALUE;
        }
        // printf("Placing phase %"PRIuFAST16"\n", i);
        engine_phases[i] = phase;
        run_generation();
    }
    #if DEBUG > 0
    printf("Phases generated\n");
    #endif
}


static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

uint64_t rng_state[4];

uint64_t rng() {
    #define s rng_state
	const uint64_t result = rotl(s[1] * 5, 7) * 9;
	const uint64_t t = s[1] << 17;
	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = rotl(s[3], 45);
    #undef s
	return result;
}

#ifndef BRUH
void init_rng() {
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) {
        perror("Error opening /dev/urandom");
        exit(1);
    }
    ssize_t size = read(fd, rng_state, 4 * sizeof(uint64_t));
    if (size < (ssize_t)(4 * sizeof(uint64_t))) {
        perror("Error reading /dev/urandom");
        close(fd);
        exit(1);
    }
    close(fd);
}
#else
void init_rng() {
    printf("Enter random seed: ");
    scanf("%"PRIu64, rng_state);
    rng_state[1] = 0x0123456789ABCDEF;
    rng_state[2] = 0x5555555555555555;
    rng_state[3] = 0x0F1E2D3C4B5A6978;
    for (int i = 0; i < 128; i++) {
        rng();
    }
}
#endif


uint16 randint(uint16 range) {
    if (range == 0) {
        return 0;
    }
    uint32_t value = rng();
    uint32_t max = ((((uint64_t)1 << 32) - 1) / (uint64_t)range) * (uint64_t)range;
    while (value > max) {
        value = rng();
    }
    uint16 out = value % range;
    if (out >= range) {
        return randint(range);
    } else {
        return out;
    }
}

uint8_t initial_pattern[SIZEVALUE];
uint16 ip_bottom = 0;
uint16 ip_right = 0;

typedef struct engine_info {
    uint16 x;
    uint16 y;
    uint16 phase;
} engine_info;

engine_info* global_engines;

void create_soup() {
    clear();
    uint32 stop = ((uint32)ip_bottom << WIDTH) + STARTX;
    uint32 max = (STARTY << WIDTH) + ip_right;
    for (uint32 start = (STARTY << WIDTH) + STARTX; start < stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            data[i] = 0;
        }
        max += WIDTHVALUE;
    }
    uint16 y;
    uint16 x;
    engine_phase* phase;
    // clock_t start = clock();
    if (use_random_soups) {
        top = STARTY;
        bottom = 0;
        left = STARTX;
        right = 0;
        x = STARTY;
        y = STARTX;
        for (uint16 i = 0; i < engines; i++) {
            phase = engine_phases[randint(ENGINEPHASES)];
            // idk why this is needed, someone should figure out why because this probably slows it down
            if (phase == NULL) {
                generate_phases();
                create_soup();
                return;
            }
            uint16 j = 0;
            for (uint16 cy = 0; cy < phase->height; cy++) {
                uint32 i = ((uint32)(y + cy) << WIDTH) + x;
                for (uint16 cx = 0; cx < phase->width; cx++) {
                    initial_pattern[i] = phase->data[j];
                    i++;
                    j++;
                }
            }
            if (x + phase->width > right) {
                right = x + phase->width;
            }
            if (y + phase->height > bottom) {
                bottom = y + phase->height;
            }
            x = STARTX + randint(max_x_sep);
            y += MINY + randint(MAXY - MINY + 1);
        }
    } else {
        top = STARTY;
        bottom = 0;
        left = STARTX;
        right = 0;
        y = STARTY;
        for (uint16 i = 0; i < engines; i++) {
            engine_info engine = global_engines[i];
            y += engine.y;
            phase = engine_phases[engine.phase];
            // #if DEBUG > 0
            // printf("Placing phase %"PRIuFAST16" at x = %"PRIuFAST16", y = %"PRIuFAST16"\n", engine.phase, engine.x + STARTX, engine.y + y);
            // #endif
            for (uint16 cy = 0; cy < phase->height; cy++) {
                uint32 i = ((y + cy) << WIDTH) + STARTX + engine.x;
                uint16 j = cy * phase->width;
                uint16 end = j + phase->height;
                for (; j < end; j++) {
                    initial_pattern[i] = phase->data[j];
                    i++;
                }
            }
            if (STARTX + engine.x + phase->width > right) {
                right = STARTX + engine.x + phase->width;
            }
            if (STARTY + engine.y + phase->height > bottom) {
                bottom = STARTY + engine.y + phase->height;
            }
        }
        for (uint16 i = engines - 1; i >= 0; i++) {
            engine_info engine = global_engines[i];
            if (i == 0) {
                engine.phase = (engine.phase + 1) % ENGINEPHASES;
            } else {
                if (engine.y != MAXY) {
                    engine.y++;
                    break;
                } else {
                    engine.y = MINY;
                    if (engine.x != STARTX + max_x_sep) {
                        engine.x++;
                        break;
                    } else {
                        engine.x = STARTX;
                        engine.phase = (engine.phase + 1) % ENGINEPHASES;
                    }
                }
            }
        }
    }
    stop = ((uint32)bottom << WIDTH) + left;
    max = ((uint32)top << WIDTH) + right;
    for (uint32 start = (top << WIDTH) + left; start < stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            data[i] = initial_pattern[i];
        }
        max += WIDTHVALUE;
    }
    ip_bottom = bottom;
    ip_right = right;
    // printf("Generated soup in %.3f seconds\n", ((double)clock() - (double)start) / CLOCKS_PER_SEC);
}


typedef struct pattern_data {
    uint16 top;
    uint16 left;
    uint16 height;
    uint16 width;
    uint32 size;
    uint32 population;
    uint64_t hash;
    uint32_t data[];
} pattern_data;

uint32 phase_number;
pattern_data** phase_cache;

void cache_phase() {
    uint16 height = bottom - top;
    uint16 width = right - left;
    uint32 size = height * width;
    uint32 population = 0;
    uint32 data_length = size >> 5;
    if (size % 32 != 0) {
        data_length++;
    }
    if (data_length == 0) {
        data_length = 1;
    }
    // printf("height = %"PRIuFAST16", width = %"PRIuFAST16", data_length = %"PRIuFAST32"\n", height, width, data_length);
    pattern_data* out = malloc(sizeof(pattern_data) + data_length * sizeof(uint32_t));
    // printf("%"PRIuFAST16" %"PRIuFAST16" %"PRIuFAST32"\n", ip_bottom, ip_right, data_length);
    out->top = top;
    out->left = left;
    out->height = height;
    out->width = width;
    out->size = size;
    for (uint32 i = 0; i < data_length; i++) {
        out->data[i] = 0;
    }
    uint32 data_i = 0;
    uint16 bit_number = 0;
    uint32 stop = ((uint32)bottom << WIDTH) + left;
    uint32 max = ((uint32)top << WIDTH) + right;
    uint32 num_cells = 0;
    for (uint32 start = (top << WIDTH) + left; start < stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            // if (data_i > data_length) {
            //     printf("you are probably writing past the array bounds ");
            // }
            if (data[i]) {
                population++;
            }
            out->data[data_i] |= data[i] << bit_number;
            bit_number++;
            if (bit_number == 32) {
                bit_number = 0;
                data_i++;
            }
            num_cells++;
            // printf("%"PRIuFAST32" ", num_cells);
        }
        max += WIDTHVALUE;
    }
    // printf("\n");
    uint64_t hash = 0;
    for (uint32 i = 0; i < data_length; i += 4) {
        hash += (uint64_t)out->data[i] << 32;
        if (data_length - i > 4) {
            hash += (uint64_t)out->data[i + 1];
            hash ^= (uint64_t)out->data[i + 2] << 32;
            hash ^= (uint64_t)out->data[i + 3];
        } else if (i + 1 < data_length) {
            hash += (uint64_t)out->data[i + 1];
            if (i + 2 < data_length) {
                hash ^= (uint64_t)out->data[i + 2] << 32;
                if (i + 3 < data_length) {
                    hash ^= (uint64_t)out->data[i + 3];
                }
            }
        }
        hash = (hash << 16) | (hash >> 48);
    }
    out->hash = hash;
    out->population = population;
    phase_cache[phase_number] = out;
}


#if REDUCEPERIOD > 0
uint16 gcd(uint16 a, uint16 b) {
    while (b != 0) {
        uint16 temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}
#endif

uint64_t check_for_spaceship() {
    pattern_data* data;
    pattern_data* current = phase_cache[phase_number];
    uint16 period = CHECKINTERVAL;
    // printf("Checking\n");
    for (uint16 i = phase_number - 1; i < UINT_FAST16_MAX; i--) {
        data = phase_cache[i];
        // printf("Checking generation %"PRIuFAST16" %"PRIuFAST32" %"PRIuFAST32"\n", i, current->population, data->population);
        if (current->hash == data->hash && current->population == data->population && current->height == data->height && current->width == data->width) {
            // #if DEBUG > 0
            // printf("Match found\n");
            // #endif
            for (uint16 i = 0; i < current->size >> 5; i++) {
                if (current->data[i] != data->data[i]) {
                    goto not_found;
                }
            }
            int32_t dx = current->left - data->left;
            int32_t dy = current->top - data->top;
            if (dx < 0) {
                dx = -dx;
            }
            if (dy < 0) {
                dy = -dy;
            }
            if (dx == 0 || dy != 0) {
                return period << 32;
            } else {
                #if REDUCEPERIOD > 0
                uint16 num = gcd(dx, period);
                dx /= num;
                period /= num;
                #endif
                return ((uint64_t)period << 32) | (uint64_t)dx;
            }
            not_found:;
        }
        period += CHECKINTERVAL;
    }
    return 0;
}


uint32 ships = 0;
uint64_t speeds[MAXSHIPS];
char* rles;

uint64_t parse_speed(char* data, uint32 i, uint32 end) {
    uint64_t out = atoi(data + i);
    uint32 begin_period = 0;
    for (; i < end; i++) {
        if (data[i] == '/') {
            begin_period = i + 1;
            break;
        }
    }
    if (begin_period == 0) {
        printf("Invalid speed\n");
        exit(1);
    }
    out |= (uint64_t)atoi(data + begin_period) << 32;
    return out;
}

void read_state() {
    #ifndef BRUH
    FILE* f = fopen(state_file, "r");
    if (f == 0) {
        perror("Error opening state file");
        exit(1);
    }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* data = malloc(size);
    size_t bytes = fread(data, 1, size, f);
    if (bytes < (size_t)size) {
        perror("Error reading state file");
        fclose(f);
        exit(1);
    }
    ships = atoi(data);
    uint32 i;
    for (i = 0; i < size; i++) {
        if (data[i] == '\n') {
            goto parse_speeds;
        }
    }
    printf("Invalid state file: could not find second line\n");
    fclose(f);
    exit(1);
    parse_speeds:;
    i++;
    uint32 start = i;
    uint16 current_speed = 0;
    for (; i < size; i++) {
        if (data[i] == ' ') {
            speeds[current_speed] = parse_speed(data, start, i);
            current_speed++;
            start = i;
        } else if (data[i] == '\n') {
            goto parse_rles;
        }
    }
    printf("Invalid state file: could not find third line\n");
    fclose(f);
    exit(1);
    parse_rles:;
    i++;
    uint32 rle_size = size - i;
    if (rles != NULL) {
        free(rles);
    }
    if (rle_size != 0) {
        rles = malloc(rle_size);
        memcpy(rles, data + i, rle_size - 1);
        rles[rle_size - 1] = '\0';
    } else {
        rles = NULL;
    }
    free(data);
    fclose(f);
    #endif
}

void add_ship(uint64_t speed) {
    #if SKIPDUPLICATES > 0
    for (uint16 i = 0; i < ships; i++) {
        if (speeds[i] == speed) {
            // printf("Duplicate %"PRIu64"c/%"PRIu64" found\n", speed & 65535, speed >> 32);
            return;
        }
    }
    #endif
    speeds[ships] = speed;
    ships++;
    printf("%"PRIu64"c/%"PRIu64" found! (%"PRIuFAST16" NRSS total)\n", speed & 65535, speed >> 32, ships);
    #ifndef BRUH
    FILE* f = fopen(state_file, "w");
    if (f == 0) {
        perror("Error opening state file");
        exit(1);
    }
    #else
    printf("=== begin file data ===\n");
    #define fprintf(x, y, ...) printf(y, ## __VA_ARGS__)
    #define fputc(char, f) printf("%c", char)
    #endif
    fprintf(f, "%"PRIuFAST32" NRSS\n", ships);
    for (uint16 i = 0; i < ships; i++) {
        fprintf(f, "%"PRIu64"c/%"PRIu64" ", speeds[i] & 65535, speeds[i] >> 32);
    }
    // printf("ships: %"PRIuFAST16"\n", ships);
    #ifndef BRUH
    if (ships != 1 && rles != NULL) {
        fprintf(f, "\n%s", rles);
    }
    #endif
    fprintf(f, "\n# %"PRIu64"c/%"PRIu64" ", speed & 65535, speed >> 32);
    uint16 height = ip_bottom - STARTY;
    uint16 width = ip_right - STARTX;
    fprintf(f, "\nx = %"PRIuFAST16", y = %"PRIuFAST16", rule = "RULESTR"\n", width, height);
    uint32 stop = ((uint32)ip_bottom << WIDTH) + STARTX;
    uint32 max = (STARTY << WIDTH) + ip_right;
    char prev = '\0';
    uint16 count = 0;
    #define addchar(c) { \
        if (prev == (c)) { \
            count++; \
        } else { \
            for (uint16 i = 0; i < count; i++) { \
                fputc(prev, f); \
            } \
            prev = (c); \
            count = 1; \
        } \
    }
    for (uint32 start = (STARTY << WIDTH) + STARTX; start < stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            addchar(initial_pattern[i] ? 'o' : 'b');
        }
        max += WIDTHVALUE;
        addchar('$');
    }
    #undef addchar
    fprintf(f, "!\n");
    #ifndef BRUH
    fclose(f);
    read_state();
    #else
    printf("=== end file data ===\n");
    #endif
}


uint64 soup_count;

void run_soup() {
    #if DEBUG > 0
    printf("Creating soup... ");
    #endif
    create_soup();
    #if DEBUG > 0
    printf("complete\n");
    #endif
    phase_number = 0;
    uint64_t speed;
    // #define topm1 (top - 1)
    // #define bottomp1 (bottom + 1)
    // #define leftm1 (left - 1)
    // #define rightp1 (right + 1)
    // #define heightp2 (bottomp1 - topm1)
    // #define widthp2 (rightp1 - leftm1)
    // printf("top: %"PRIuFAST16", bottom: %"PRIuFAST16", left: %"PRIuFAST16", right: %"PRIuFAST16, top, bottom, left, right);
    // printf("\nx = %"PRIuFAST16", y = %"PRIuFAST16"\n", widthp2, heightp2);
    // char* row = malloc((widthp2 + 1) * sizeof(char));
    // row[widthp2] = '\0';
    // for (uint16 y = topm1; y < bottomp1; y++) {
    //     uint32 i = ((uint32)y << WIDTH) + leftm1;
    //     for (uint16 x = 0; x < widthp2; x++) {
    //         row[x] = data[i] ? 'o' : 'b';
    //         i++;
    //     }
    //     printf("%s$", row);
    // }
    // free(row);
    uint16 i;
    // clock_t start_time = clock();
    for (i = 0; i < max_period; i++) {
        uint16 pop = 0;
        uint32 stop = ((uint32)bottom << WIDTH) + left;
        uint32 max = ((uint32)top << WIDTH) + right;
        for (uint32 start = (top << WIDTH) + left; start < stop; start += WIDTHVALUE) {
            for (uint32 i = start; i < max; i++) {
                if (data[i]) {
                    pop++;
                }
            }
            max += WIDTHVALUE;
        }
        #if DEBUG > 0
        #define topm1 (top - 1)
        #define bottomp1 (bottom + 1)
        #define leftm1 (left - 1)
        #define rightp1 (right + 1)
        #define heightp2 (bottomp1 - topm1)
        #define widthp2 (rightp1 - leftm1)
        printf("Running generation %"PRIuFAST16" (population %"PRIuFAST32")\n", i, pop);
        #if DEBUG > 1
        printf("x = %"PRIuFAST16", y = %"PRIuFAST16"\n", widthp2, heightp2);
        char* row = malloc((widthp2 + 1) * sizeof(char));
        row[widthp2] = '\0';
        for (uint16 y = topm1; y < bottomp1; y++) {
            uint32 i = ((uint32)y << WIDTH) + leftm1;
            for (uint16 x = 0; x < widthp2; x++) {
                row[x] = data[i] ? '1' : '0';
                i++;
            }
            printf("%s\n", row);
        }
        free(row);
        #endif
        #endif
        if (top < 2 || bottom > HEIGHTVALUE - 2 || left < 2 || right > WIDTHVALUE - 2) {
            break;
        }
        if (!run_generation()) {
            break;
        }
        if (i % CHECKINTERVAL == 0) {
            #if DEBUG > 0
            printf("Checking for spaceship... ");
            #endif
            cache_phase();
            if ((speed = check_for_spaceship()) != 0) {
                if (speed < MINPERIOD) {
                    #if DEBUG > 0
                    printf("complete, less than min period\n");
                    #endif
                    break;
                }
                #if SKIPOSCILLATORS > 0
                if ((speed & 65535) == 0) {
                    #if DEBUG > 0
                    printf("complete, skipped oscillator\n");
                    #endif
                    break;
                }
                #endif
                #if DEBUG > 0
                printf("complete, true\n");
                #endif
                add_ship(speed);
                break;
            }
            phase_number++;
            #if DEBUG > 0
            printf("complete, false\n");
            #endif
        }
    }
    for (uint16 i = 0; i < phase_number; i++) {
        free(phase_cache[i]);
    }
    // printf("Soup stabilized after %"PRIuFAST16" generations (%.3f seconds)\n", i, (double)(clock() - start_time) / (double)CLOCKS_PER_SEC);
    soup_count++;
}


clock_t start_clock;
clock_t prev_clock;
uint64_t soups;
uint64_t prev_soups;
int64_t max_soups;

void show_status_force(clock_t current) {
    if (max_soups < 0) {
        printf("%"PRIu64" soups completed (%.3f soups/second current, %.3f overall)\n", soups, (double)(soups - prev_soups) / ((double)(current - prev_clock) / CLOCKS_PER_SEC), (double)soups / ((double)(current - start_clock) / CLOCKS_PER_SEC));
    } else {
        printf("%"PRIu64" soups completed (%.3f%%, %.3f soups/second current, %.3f overall)\n", soups, (double)soups / (double)max_soups * 100, (double)(soups - prev_soups) / ((double)(current - prev_clock) / CLOCKS_PER_SEC), (double)soups / ((double)(current - start_clock) / CLOCKS_PER_SEC));
    }
}

void show_status() {
    clock_t current = clock();
    if ((((double)current - (double)prev_clock) / (double)CLOCKS_PER_SEC) >= 10) {
        show_status_force(current);
        prev_clock = current;
        prev_soups = soups;
    }
}

void cleanup() {
    show_status_force(clock());
    free(phase_cache);
    for (uint16 i = 0; i < ENGINEPHASES; i++) {
        free(engine_phases[i]);
    }
    free(rles);
    free(global_engines);
}

void on_sigint(int something) {
    printf("\n");
    cleanup();
    exit(1);
}

int main(int argc, char** argv) {
    #ifndef BRUH
    if (argc != 6) {
        printf("Usage: nrss <engine-count> <max-x-seperation> <max-period> <randomize-soups-1-or-0> <state-file>\n");
        return 1;
    }
    #else
    if (argc != 5) {
        printf("Usage: nrss <engine-count> <max-x-seperation> <max-period> <randomize-soups-1-or-0>\n");
        return 1;
    }
    #endif
    engines = atoi(argv[1]);
    max_x_sep = atoi(argv[2]);
    max_period = atoi(argv[3]);
    use_random_soups = (bool)atoi(argv[4]);
    #ifndef BRUH
    state_file = argv[5];
    #endif
    phase_cache = malloc((max_period / CHECKINTERVAL + 1) * sizeof(pattern_data*));
    generate_phases();
    read_state();
    global_engines = malloc(sizeof(engine_info) * engines);
    for (uint16 i = 0; i < engines; i++) {
        global_engines[i].x = 0;
        global_engines[i].y = i == 0 ? 0 : MINY;
        global_engines[i].phase = 0;
    }
    start_clock = clock();
    prev_clock = start_clock;
    prev_soups = 0;
    signal(SIGINT, on_sigint);
    if (use_random_soups) {
        init_rng();
        max_soups = -1;
        printf("Starting search\n");
        while (true) {
            run_soup();
            soups++;
            show_status();
        }
    } else {
        max_soups = ENGINEPHASES;
        for (uint16 i = 0; i < engines - 1; i++) {
            max_soups *= (int64_t)ENGINEPHASES * (int64_t)(MAXY - MINY) * ((int64_t)max_x_sep + 1);
        }
        printf("Searching %"PRIu64" soups\n", max_soups);
        for (uint64_t i = 0; i < max_soups; i++) {
            run_soup();
            soups++;
            show_status();
        }
    }
    cleanup();
    return 0;
}
