/* 
 *  author: Suhas Vittal
 *  date:   10 July 2023
 * */

#ifndef UF_FAST_COMMON_h
#define UF_FAST_COMMON_h

//
// Temporary defines: (place in their own header in practice)
//

#define UFF_N_DETECTORS 16 
#define UFF_BOUNDARY    UFF_N_DETECTORS
#define MAX_DEFECTS 32

//
// Actual definitions:
//

typedef long long int   int64_t;
typedef int             int32_t;
typedef char            int8_t;

typedef unsigned long long int  uint64_t;
typedef unsigned int            uint32_t;
typedef unsigned char           uint8_t;

typedef uint32_t   event_t;

#define AM_ADJACENT_MASK    0x01
#define AM_FRAME_MASK       0x3e

#define AM_HALF_1_OFFSET    6
#define AM_HALF_2_OFFSET    7

extern uint8_t     adjacency_matrix[UFF_N_DETECTORS][UFF_N_DETECTORS];
                            // Each adjacency matrix entry has 8 bits such that
                            //  bit 0 = is adjacent (is there even an edge??)
                            //  bit 1-5 = any frame changes
                            //  bit 6-7 = half edges grown (if fully grown, =0b11)
                            // Diagonal elements store boundary connection data.
#define AL_WIDTH        24
#define AL_VALID_MASK   (1 << 31)

extern uint32_t    adjacency_lists[UFF_N_DETECTORS][AL_WIDTH];
                            // bits 0-30 = other event (surely you don't need 31 bits!)
                            // bit 31 = valid

struct cluster_t {    // 40/72 bits.
    uint8_t             data;   // bits 0-5 is size, bit 6 is parity, bit 7 is frozen
    struct cluster_t*   root;   // Either 32-bit or 64-bit, depending on addressability.
};

// We're using the same FIFO structure for vertices and edges.
// Simply pack the endpoints together when using edges.

#define FIFO_SIZE   256
#define FIFO_MASK   ((FIFO_SIZE)-1)

struct fifo_t {
    event_t     data[256];
    uint8_t        size;
    uint8_t        head;
};

//
//  HELPER FUNCTIONS:
// 

uint8_t*   get_adj_data(event_t, event_t);

uint8_t    get_size(struct cluster_t*);
uint8_t    can_expand(struct cluster_t*);

void    fifo_init(struct fifo_t*);
void    fifo_push(struct fifo_t*, event_t);
event_t fifo_pop(struct fifo_t*);

struct cluster_t*   uf_find(struct cluster_t*);
void                uf_merge(struct cluster_t*, struct cluster_t*);

//
// MAIN FUNCTION:
//

#define PACKED_INPUT_DATA_ARRAY_SIZE    (((UFF_N_DETECTORS) >> 6)+1)

int32_t decode(uint64_t*, uint64_t*);

#endif  // UF_FAST_COMMON_h

 
