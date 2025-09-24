
/*
Copyright 2025 speedydelete

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// searches for orthogonal NRSS (nonadjustible reduced speed of ship) that face east in square-grid range-1 Moore neighborhood cellular automata
// see https://conwaylife.com/forums/viewtopic.php?f=11&t=6352&p=218310 for more informatio
// to compile: gcc -Wall -Werror -Ofast -o nrss nrss.c
// to use: nrss <engine-count> <max-x-seperation> <max-period> <randomize-soups-1-or-0> <state-file>
// when randomization is off it will try every possible combination of engines

#include <stdbool.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <time.h>

// parameters

#define RULESTR "B2-ak3ce4eikqrz5-iknq6-ek8/S1c2aek3aekn4eiknry5eiky6-ei7c8"

// maximum height and width, these are the base-2 logarithms
// if a pattern exceeds these, it is reported as potential infinite growth
// should not be higher than 16
#define HEIGHT 12
#define WIDTH 12

// starting position
#define STARTX 64
#define STARTY ((1 << HEIGHT) / 2 - 64)

// min and max y seperation between engines
#define MINY 7
#define MAXY 12

// put an engine in the data array at the specified position
// i is the array index of the (x, y) coordinates of where to put the engine
// the default one is 2o$o$2o!
static inline void put_engine(uint8_t data[], uint_fast32_t i) {
    data[i] = 1;
    data[i + 1] = 1;
    data[i + WIDTH] = 1;
    data[i + WIDTH + WIDTH] = 1;
    data[i + WIDTH + WIDTH + 1] = 1;
}

// engine dimensions
#define ENGINEHEIGHT 3
#define ENGINEWIDTH 2

// the number of phases of the engine
#define ENGINEPHASES 100

// whether to skip oscillators, comment out to keep oscillators
#define SKIPOSCILLATORS

// whether to skip potential linear growth, uncomment to skip it
// #define SKIPLINEARGROWTH

// maximum number of ships
#define MAXSHIPS 4096

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


typedef uint_fast16_t uint16;
typedef uint_fast32_t uint32;
typedef uint_fast64_t uint64;

#define HEIGHTVALUE (1 << HEIGHT)
#define WIDTHVALUE (1 << WIDTH)
#define SIZE (HEIGHT + WIDTH)
#define SIZEVALUE (1 << SIZE)

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
    for (uint32 start = (top << WIDTH) + left; start <= stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            data[i] = 0;
        }
        max += WIDTHVALUE;
    }
}


uint8_t temp_data[SIZEVALUE];

static inline bool run_row(uint32 i, uint16* lowX, uint16* highX) {
    uint16 tr = 0;
    tr |= (uint16)data[i - WIDTH - 1] << 5;
    tr |= (uint16)data[i - 1] << 4;
    tr |= (uint16)data[i + WIDTH - 1] << 3;
    tr |= (uint16)data[i - WIDTH] << 2;
    tr |= (uint16)data[i] << 1;
    tr |= (uint16)data[i + WIDTH];
    uint32 max = i - left + right + 2;
    uint8_t value;
    bool any_changed = false;
    for (i++; i <= max; i++) {
        tr = (tr << 3) & 511;
        tr |= (uint16)data[i - WIDTH + 1] << 2;
        tr |= (uint16)data[i + 1] << 1;
        tr |= (uint16)data[i + WIDTH + 1];
        value = transitions[tr];
        if (value) {
            any_changed = true;
            if (i < *lowX) {
                *lowX = i;
            }
            if (i > *highX) {
                *highX = i;
            }
        }
        temp_data[i] = value;
    }
    return any_changed;
}

bool run_generation() {
    uint16 lowX = WIDTH;
    uint16 highX = 0;
    uint16 lowY = HEIGHT;
    uint16 highY = 0;
    bool any_changed;
    uint32 i = ((top - 1) << WIDTH) + left - 1;
    for (uint16 y = top - 1; y <= bottom + 1; y++) {
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
    top = lowY;
    bottom = highY;
    left = lowX;
    right = highX;
    uint32 stop = ((uint32)bottom << WIDTH) + left;
    uint32 max = ((uint32)top << WIDTH) + right;
    for (uint32 start = (top << WIDTH) + left; start <= stop; start += WIDTHVALUE) {
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
    put_engine(data, (STARTY << WIDTH) + STARTX);
    top = STARTY;
    bottom = STARTY + ENGINEHEIGHT;
    left = STARTX;
    right = STARTX + ENGINEWIDTH;
    printf("%"PRIuFAST16" %"PRIuFAST16" %"PRIuFAST16" %"PRIuFAST16"\n", top, bottom, left, right);
    for (int i = 0; i < ENGINEPHASES; i++) {
        uint16 height = bottom - top;
        uint16 width = right - left;
        uint32 size = height * width;
        printf("%"PRIuFAST16" %"PRIuFAST16" %"PRIuFAST16" %"PRIuFAST16"\n", top, bottom, left, right);
        engine_phase* phase = malloc(sizeof(engine_phase) + size);
        phase->height = height;
        phase->width = width;
        int j = 0;
        uint32 stop = ((uint32)bottom << WIDTH) + left;
        uint32 max = ((uint32)top << WIDTH) + right;
        for (uint32 start = (top << WIDTH) + left; start <= stop; start += WIDTHVALUE) {
            for (uint32 i = start; i < max; i++) {
                phase->data[j++] = data[i];
            }
            max += WIDTHVALUE;
        }
        engine_phases[i] = phase;
        run_generation();
    }
}


static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

uint64_t rng_state[4];

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

uint16 randint(uint16 range) {
    uint32_t value = rng();
    uint32_t max = (((uint64_t)1 << 32) - 1) / range;
    while (value > max) {
        value = rng();
    }
    return (uint16)(value % range);
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
    uint32 stop = ((uint32)ip_bottom << WIDTH) + STARTX;
    uint32 max = (STARTY << WIDTH) + ip_right;
    for (uint32 start = (STARTY << WIDTH) + STARTX; start <= stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            data[i] = 0;
        }
        max += WIDTHVALUE;
    }
    put_engine(initial_pattern, (STARTY << WIDTH) + STARTX);
    uint16 y = STARTY;
    uint16 x;
    engine_phase* phase;
    if (use_random_soups) {
        for (uint16 i = 0; i < engines - 1; i++) {
            x = STARTX + randint(max_x_sep);
            y += MINY + randint(MAXY - MINY + 1);
            phase = engine_phases[randint(ENGINEPHASES)];
            for (uint16 cy = 0; cy < phase->height; cy++) {
                uint32 i = ((y + cy) << WIDTH) + x;
                for (uint16 cx = 0; cx < phase->width; cx++) {
                    initial_pattern[i] = phase->data[i];
                    i++;
                }
            }
        }
    } else {
        for (int i = engines - 2; i >= 0; i++) {
            engine_info engine = global_engines[i];
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
        for (uint16 i = 0; i < engines - 1; i++) {
            engine_info engine = global_engines[i];
            y += engine.y;
            phase = engine_phases[engine.phase];
            for (uint16 cy = 0; cy < phase->height; cy++) {
                uint32 i = (y + cy) << WIDTH;
                for (uint16 cx = 0; cx < phase->width; cx++) {
                    initial_pattern[i] = phase->data[i];
                    i++;
                }
            }
        }
    }
}


typedef struct pattern_data {
    uint64_t hash;
    uint32 population;
    uint16 top;
    uint16 left;
    uint16 height;
    uint16 width;
    uint32 size;
    uint32_t data[];
} pattern_data;

uint32 generation;
pattern_data** previous_generation_cache;

void cache_pattern_data() {
    uint16 height = top - bottom;
    uint16 width = right - left;
    uint32 size = height * width;
    uint32 population = 0;
    uint32 data_length = size >> 5;
    pattern_data* out = malloc(sizeof(pattern_data) + data_length * sizeof(uint32_t));
    out->top = top;
    out->left = left;
    out->height = height;
    out->width = width;
    out->size = size;
    for (int i = 0; i < size; i++) {
        out->data[i] = data[i];
    }
    uint32 data_i = 0;
    uint16 bit_number = 0;
    uint32 stop = ((uint32)bottom << WIDTH) + left;
    uint32 max = ((uint32)top << WIDTH) + right;
    for (uint32 start = (top << WIDTH) + left; start <= stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            if (data[i]) {
                population++;
            }
            out->data[data_i] |= data[i] << bit_number;
            bit_number++;
            if (bit_number == 32) {
                bit_number = 0;
                data_i++;
            }
        }
        max += WIDTHVALUE;
    }
    uint64_t hash = 0;
    for (int i = 0; i < data_length; i += 4) {
        hash += (uint64_t)out->data[i] << 32;
        hash += (uint64_t)out->data[i + 1];
        hash ^= (uint64_t)out->data[i + 2] << 32;
        hash ^= (uint64_t)out->data[i + 3];
        hash = (hash << 16) | (hash >> 48);
    }
    out->hash = hash;
    previous_generation_cache[generation] = out;
}


uint64_t check_for_spaceship() {
    pattern_data* data;
    pattern_data* current = previous_generation_cache[generation];
    uint16 period = 1;
    for (uint16 i = generation - 1; i >= 0; i++) {
        data = previous_generation_cache[i];
        if (current->hash == data->hash && current->population == data->population && current->height == data->height && current->width == data->width) {
            for (int i = 0; i < current->size; i++) {
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
            if (dy != 0 || (dx == 0 && dy != 0)) {
                return period << 32;
            } else {
                return ((uint64_t)period << 32) | ((uint64_t)dy << 16) | (uint64_t)dx;
            }
            not_found:;
        }
    }
    return 0;
}


uint32 ships = 0;
uint64_t speeds[MAXSHIPS];
char* rles;

uint64_t parse_speed(char* data, uint32 i, uint32 end) {
    uint32 begin_period = 0;
    for (; i < end; i++) {
        if (data[i] == '/') {
            begin_period = i + 1;
        }
    }
    if (begin_period == 0) {
        printf("Invalid speed\n");
        exit(1);
    }
    uint64_t out = atoi(data);
    out |= (uint64_t)atoi(data + begin_period) << 32;
    return out;
}

void read_state() {
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
    uint32 rle_size = size - i - 1;
    if (rles != NULL) {
        free(rles);
    }
    rles = malloc(rle_size);
    memcpy(data + i, rles, rle_size);
    free(data);
}

void add_ship(uint64_t speed) {
    printf("%"PRIu64"c/%"PRIu64" found (%"PRIuFAST16" NRSS total)!\n", speed & 65535, speed >> 32, ships);
    FILE* f = fopen(state_file, "w");
    if (f == 0) {
        perror("Error opening state file");
        exit(1);
    }
    ships++;
    speeds[ships] = speed;
    fprintf(f, "%"PRIuFAST32" NRSS", ships);
    for (uint16 i = 0; i < ships; i++) {
        fprintf(f, "%"PRIu64"c/%"PRIu64" ", speeds[i] & 65535, speeds[i] >> 32);
    }
    fprintf(f, "%s", rles);
    fprintf(f, "# %"PRIu64"c/%"PRIu64" ", speed & 65535, speed >> 32);
    uint16 height = ip_bottom - STARTY;
    uint16 width = ip_right - STARTX;
    fprintf(f, "x = %"PRIuFAST16", y = %"PRIuFAST16", rule = "RULESTR"\n", width, height);
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
    for (uint32 start = (STARTY << WIDTH) + STARTX; start <= stop; start += WIDTHVALUE) {
        for (uint32 i = start; i < max; i++) {
            addchar(data[i] ? 'o' : 'b');
            data[i] = 0;
        }
        max += WIDTHVALUE;
        addchar('$');
    }
    #undef addchar
    fprintf(f, "!\n");
    fclose(f);
    read_state();
}


uint64 soup_count;

void run_soup() {
    create_soup();
    generation = 0;
    uint64_t data;
    for (int i = 0; i < max_period; i++) {
        if (!run_generation()) {
            break;
        }
        if (top == 0 || bottom == HEIGHT || left == 0 || right == WIDTH) {
            #ifndef SKIPLINEARGROWTH
            #endif
            break;
        }
        cache_pattern_data();
        if ((data = check_for_spaceship()) != 0) {
            #ifdef SKIPOSCILLATORS
            if (data == (uint64_t)1 << 33) {
                break;
            }
            #endif
            for (uint16 speed = 0; speed < ships; speed++) {
                if (speed == data) {
                    goto duplicate;
                }
            }
            add_ship(data);
            duplicate:;
            break;
        }
    }
    for (int i = 0; i < generation; i++) {
        free(previous_generation_cache[i]);
    }
    soup_count++;
}


clock_t start_clock;
clock_t prev_clock;
uint64_t soups;
int64_t max_soups;

void show_status_force(clock_t current) {
    if (max_soups < 0) {
        printf("%"PRIu64" soups completed (%.3f soups/second current, %.3f overall)\n", soups, (double)(current - prev_clock) / CLOCKS_PER_SEC * 100, (double)(current - prev_clock) / CLOCKS_PER_SEC * 100);
    } else {
        printf("%"PRIu64" soups completed (%.3f%%, %.3f soups/second current, %.3f overall)\n", soups, (double)soups / (double)max_soups * 100, (double)(current - prev_clock) / CLOCKS_PER_SEC * 100, (double)(current - prev_clock) / CLOCKS_PER_SEC * 100);
    }
}

void show_status() {
    clock_t current = clock();
    if ((((double)current - (double)prev_clock) / (double)CLOCKS_PER_SEC) >= 10) {
        show_status_force(current);
        prev_clock = current;
    }
}

void cleanup() {
    show_status_force(clock());
    free(previous_generation_cache);
    for (int i = 0; i < ENGINEPHASES; i++) {
        free(engine_phases[i]);
    }
    free(rles);
}

void on_sigint(int something) {
    cleanup();
    exit(1);
}

int main(int argc, char** argv) {
    if (argc != 6) {
        printf("Usage: nrss <engine-count> <max-x-seperation> <max-period> <randomize-soups-1-or-0> <state-file>\n");
        return 1;
    }
    engines = atoi(argv[1]);
    max_x_sep = atoi(argv[2]);
    max_period = atoi(argv[3]);
    use_random_soups = (bool)atoi(argv[4]);
    state_file = argv[5];
    previous_generation_cache = malloc(max_period * sizeof(pattern_data*));
    generate_phases();
    read_state();
    prev_clock = clock();
    signal(SIGINT, on_sigint);
    if (use_random_soups) {
        max_soups = -1;
        while (true) {
            run_soup();
            soups++;
            show_status();
        }
    } else {
        max_soups = ENGINEPHASES;
        for (int i = 0; i < engines - 1; i++) {
            max_soups *= (int64_t)ENGINEPHASES * (int64_t)(MAXY - MINY) * (int64_t)max_x_sep;
        }
        for (uint64_t i = 0; i < max_soups; i++) {
            run_soup();
            soups++;
            show_status();
        }
    }
    // stuff
    cleanup();
    return 0;
}


