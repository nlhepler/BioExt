
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int min_overlap;
    int min_reads;
    int tol_gaps;
    int tol_ambigs;
} opts_t;

typedef struct {
    int col;
    int ins;
    char nuc;
} triple_t;

typedef struct {
    triple_t * data;
    int len;
    int lpos;
    int rpos;
    int ncontrib;
} aligned_t;

char nuc2bits( const char nuc );

char bits2nuc( const char bits );

aligned_t * merge__(
    const int nreads,
    const aligned_t * const reads,
    const opts_t * const opts,
    int * nclusters
    );

void aligned_destroy( aligned_t * const read );

#ifdef __cplusplus
}
#endif
