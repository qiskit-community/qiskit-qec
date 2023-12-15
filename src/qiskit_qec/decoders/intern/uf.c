/*
 *  author: Suhas Vittal
 *  date:   9 July 2023
 *
 *  This is a standalone version of UF in C.
 *  No libraries. Nothing.
 *
 *  I found that uniform UF is better than
 *  non-uniform UF (maybe because I didn't
 *  implement it correctly?). So this
 *  is uniform UF.
 *
 *  Nothing should be allocated on the heap.
 * */

// An external header file should define the following
//
//  UFF_N_DETECTORS -- number of events in the decoding graph
//  BOUNDARY    -- event id corresponding to boundary
//  MAX_DEFECTS
//

#include "uf_common.h" // Replace with data header later.

#define NULL    0

//#define DEBUG

#ifdef DEBUG
#include <stdio.h>
#endif

//
// Declarations of internal structures:
//

struct cluster_t    cluster_table[MAX_DEFECTS];
uint8_t                clusters_in_use;

struct cluster_t*   vertex_to_cluster[UFF_N_DETECTORS+1];

//
// DECODER DEFINITION:
//  Returns the pauli frame update. Right now, it is only setup to do
//  one stabilizer type instead of both.
//

uint64_t
decode(uint64_t* xin, uint64_t* corr_p) {
    struct fifo_t   boundary_vertices;
    struct fifo_t   forest;
    uint8_t            forest_incidence_ctr[UFF_N_DETECTORS+1];

    fifo_init(&boundary_vertices);
    fifo_init(&forest);

    // Reset some structures.
    for (uint32_t i = 0; i < UFF_N_DETECTORS; ++i) {
        vertex_to_cluster[i] = NULL;
        forest_incidence_ctr[i] = 0;
        for (uint8_t j = 0; j < AL_WIDTH; ++j) {
            event_t w = adjacency_lists[i][j];
            if (!(w & AL_VALID_MASK))   continue;
            w &= ~AL_VALID_MASK;
            uint8_t* adj_entry = get_adj_data(i, w);
            *adj_entry &= ~(0x3 << AM_HALF_1_OFFSET);
        }
    } 
    vertex_to_cluster[UFF_BOUNDARY] = NULL;
    forest_incidence_ctr[UFF_BOUNDARY] = 0;

    clusters_in_use = 0;

    event_t defects[MAX_DEFECTS];
    uint8_t hw = 0;     // Hamming weight (or number of defects)

    // Identify which defects are active.
#ifdef DEBUG
    printf("===========================\nDetectors:");
#endif

    for (uint32_t i = 0; i < PACKED_INPUT_DATA_ARRAY_SIZE; ++i) {
        if (xin[i] == 0)   continue;
        uint8_t* blk8 = (uint8_t*)(xin + i);
        for (uint8_t j = 0; j < 8; ++j) {
            if (blk8[j] == 0)   continue;
            for (uint8_t k = 0; k < 8; ++k) {
                if (blk8[j] & (1 << k)) {
                    event_t ee = (i << 6) | (j << 3) | k;
#ifdef DEBUG
                    printf(" %u", ee);
#endif
                    defects[hw++] = ee;
                    // Create cluster for defect.
                    struct cluster_t* cl = cluster_table + (clusters_in_use++);
                    cl->data = 1 << 6;
                    cl->root = cl;
                    vertex_to_cluster[ee] = cl;
                    
                    fifo_push(&boundary_vertices, ee);

                    if (hw > MAX_DEFECTS)   return 1;
                }
            }
        }
    }
    if (hw & 0x1) {
        // Add cluster for the boundary.
        defects[hw++] = UFF_BOUNDARY;
#ifdef DEBUG
        printf(" %u", UFF_BOUNDARY);
#endif
        struct cluster_t* cl = cluster_table + (clusters_in_use++);
        cl->data = 0x3 << 6;
        cl->root = cl;
        vertex_to_cluster[UFF_BOUNDARY] = cl;

        fifo_push(&boundary_vertices, UFF_BOUNDARY);
        // Also flip bit in events
        xin[UFF_BOUNDARY >> 6] |= 1L << (UFF_BOUNDARY & 0x3f);
    }
#ifdef DEBUG
    printf("\n");
#endif
    struct fifo_t next_boundary;

    while (boundary_vertices.size) {
#ifdef DEBUG
        printf("----------------------------\n");
#endif
        fifo_init(&next_boundary);
        // Go through the boundary vertices and grow the clusters.
        for (uint8_t i = 0; i < boundary_vertices.size; ++i) {
            event_t v = boundary_vertices.data[(boundary_vertices.head + i) & FIFO_MASK];
            if (v == UFF_BOUNDARY) {
                for (uint8_t j = 0; j < UFF_N_DETECTORS; ++j) {
                    uint8_t* e_data_ptr = &adjacency_matrix[j][j];
                    if (!(*e_data_ptr & 0x1))   continue;   // Not connected to the boundary.
                    *e_data_ptr = (*e_data_ptr & (AM_ADJACENT_MASK | AM_FRAME_MASK))
                                | ((*e_data_ptr & (1 << AM_HALF_1_OFFSET)) << 1)
                                | (1 << AM_HALF_1_OFFSET);
                }
            } else {
                for (uint8_t j = 0; j < AL_WIDTH; ++j) {
                    event_t w = adjacency_lists[v][j];
                    if (!(w & AL_VALID_MASK))  continue;
                    w &= ~AL_VALID_MASK;
                    uint8_t* e_data_ptr = get_adj_data(v, w);
                    *e_data_ptr = (*e_data_ptr & (AM_ADJACENT_MASK | AM_FRAME_MASK))
                                | ((*e_data_ptr & (1 << AM_HALF_1_OFFSET)) << 1)
                                | (1 << AM_HALF_1_OFFSET);
                }
            }
        }

        // Now go through and merge any clusters.
        while (boundary_vertices.size) {
            event_t v = fifo_pop(&boundary_vertices);
            if (v == UFF_BOUNDARY)  continue;
            struct cluster_t* clv = vertex_to_cluster[v];
            struct cluster_t* clvrt = uf_find(clv);
            // Merge the cluster with other new vertices or with other clusters.
            uint8_t neighbors_in_cluster = 0;
            uint8_t degree = 0;
#ifdef DEBUG
            printf("Growing cluster %x (r = %x) from %u:\n", clv, clvrt, v);
#endif
            for (uint8_t i = 0; i < AL_WIDTH; ++i) {
                event_t w = adjacency_lists[v][i];
                if (!(w & AL_VALID_MASK))  continue;
                w &= ~AL_VALID_MASK;
                ++degree;

                uint8_t e_data = *get_adj_data(v, w);
                if ((e_data >> AM_HALF_1_OFFSET) == 0x3) {
                    ++neighbors_in_cluster;
                    // Either w is part of a cluster or not. If it
                    // is part of a cluster, we must merge.
                    if (vertex_to_cluster[w] == NULL) {
                        // Not part of a cluster.
                        vertex_to_cluster[w] = clv;
                        clvrt->data |= (w == UFF_BOUNDARY) << 7;
                        fifo_push(&boundary_vertices, w);
                    } else {
                        struct cluster_t* clw = vertex_to_cluster[w];
                        struct cluster_t* clwrt = uf_find(clw);

                        if (clvrt == clwrt) continue;

                        uf_merge(clvrt, clwrt);
                        clvrt = uf_find(clvrt);   // May have been modified.
                    }
#ifdef DEBUG
                    printf("\tMerging with %u\n", w);
#endif
                    fifo_push(&forest, v);
                    fifo_push(&forest, w);
                    ++forest_incidence_ctr[v];
                    ++forest_incidence_ctr[w];
                }
            }
            
            if (neighbors_in_cluster < degree && can_expand(clvrt)) {
                fifo_push(&next_boundary, v);
            }
        }
        // Transfer the data between next_boundary and boundary_vertices.
        fifo_init(&boundary_vertices);
        while (next_boundary.size) {
            fifo_push(&boundary_vertices, fifo_pop(&next_boundary));
        }
    }

#ifdef DEBUG
    printf("Peeling:\n");
#endif

    // Peel the forest to get the correction.
    while (forest.size) {
        event_t v = fifo_pop(&forest);
        event_t w = fifo_pop(&forest);
        // Figure out which vertex is the pendant vertex.
        if (v == w) {
            forest_incidence_ctr[v] -= 2;
            continue;
        } else if (forest_incidence_ctr[v] == 0 || forest_incidence_ctr[w] == 0) {
            continue;
        } else if (forest_incidence_ctr[v] > 1 && forest_incidence_ctr[w] > 1) {
            // Not a leaf, just cycle it back through the fifo.
            fifo_push(&forest, v);
            fifo_push(&forest, w);
#ifdef DEBUG
            printf("\tcould not apply %u, %u (ctrs = %u, %u)\n",
                    v, w, forest_incidence_ctr[v], forest_incidence_ctr[w]);
            printf("\t\tfifo = (size = %u, head = %u)\n", forest.size, forest.head);
#endif
            continue;
        } else if (forest_incidence_ctr[v] > 1) {
            event_t tmp = w;
            w = v;
            v = tmp;
        }
        --forest_incidence_ctr[v];
        --forest_incidence_ctr[w];

        if (xin[v >> 6] & (1L << (v & 0x3f))) {
#ifdef DEBUG
            printf("\tapplied correction across events %u and %u (ctrs = %u, %u\n",
                    v, w, forest_incidence_ctr[v], forest_incidence_ctr[w]);
#endif
            // Update correction.
            *corr_p ^= (*get_adj_data(v, w) & AM_FRAME_MASK) >> 1;
            // Clear pendant vertex bit and flip other vertex bit in input.
            xin[v >> 6] &= ~(1L << (v & 0x3f));
            xin[w >> 6] ^= 1L << (w & 0x3f);
        } else {
#ifdef DEBUG
            printf("\tdid not apply correction across events %u and %u (ctrs = %u, %u\n",
                    v, w, forest_incidence_ctr[v], forest_incidence_ctr[w]);
#endif
        }
    }
    return 0;
}


//
// HELPER FUNCTION DEFINITIONS:
//

inline uint8_t*
get_adj_data(event_t x, event_t y) {
    if (x == UFF_BOUNDARY)      return &adjacency_matrix[y][y];
    else if (y == UFF_BOUNDARY) return &adjacency_matrix[x][x];
    else                        return x < y ? &adjacency_matrix[x][y] 
                                                : &adjacency_matrix[y][x];
}

inline uint8_t
get_size(struct cluster_t* cl) {
    return cl->data & 0x3f;
}

inline uint8_t
can_expand(struct cluster_t* cl) {
    return (cl->data & 0xc0) == (1 << 6);
}

inline void
fifo_init(struct fifo_t* f_p) {
    f_p->size = 0;
    f_p->head = 0;
}

inline void
fifo_push(struct fifo_t* f_p, event_t e) {
    f_p->data[(f_p->head + f_p->size++) & FIFO_MASK] = e;
}

event_t
fifo_pop(struct fifo_t* f_p) {
    --f_p->size;
    event_t out = f_p->data[f_p->head++];
    if (f_p->head == FIFO_SIZE) f_p->head = 0;
    return out;
}

struct cluster_t*   
uf_find(struct cluster_t* cl) {
    struct cluster_t* curr = cl;
    while (curr->root != curr) {
        curr = curr->root;
    }
    return curr;
}

void
uf_merge(struct cluster_t* rta, struct cluster_t* rtb) {
    uint8_t asz = get_size(rta);
    uint8_t bsz = get_size(rtb);
    if (asz > bsz) {
        uf_merge(rtb, rta);
        return;
    }

    rta->root = rtb;
    rta->data = ((asz + bsz) & 0x3f)
                    | ((rta->data ^ rtb->data) & (1 << 6))
                    | ((rta->data | rtb->data) & (1 << 7));
    rtb->data = rta->data;
}

