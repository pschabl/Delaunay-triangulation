/**
	Implementation of Delaunay triangulation using the quad-edge data structure
	and a "divide and conquer" algorithm as proposed by Guibas & Stolfi.
	The program can also compute the convex hull using the "Graham scan" method.
	On success, the program will output requested image files in TGA format.

	Notes:

	. Currently point coordinates are of type signed short.

	Compiling:

	. debug build:
		gcc -Wall -O2 -g -DDEBUG_VERSION guibas_stolfi.c

	. release build:
		gcc -Wall -O2 -DNDEBUG guibas_stolfi.c

	Developed on a x86_64 GNU/Linux machine.

	Platform dependend/linux-specific code:

	. int clock_gettime( clockid_t clk_id, struct timespec *tp )
		Used for timing.
	. char *getenv( const char *name )
		Get environment variable. Here: 'LANG'.
	. char *setlocale( int category, const char *locale )
		Set/query the programs current locale.
		Here for printing numbers with thousand separator.

	TODO

	. Compute and draw the Voronoi diagram.
	. Read points from a file, e.g. PNT
	. Command line argument evaluation is incomplete and not robust.
	. During debugging all images get drawn, regardless of image options.
	. Create postscript plot. Use pslib?
	. Allow for ccordinates in floating or fixed point,
	  and if so, be aware of precision and rounding issues!).
	. Store results in a file, rather than creating a diagram.
	. Profiling and optimization.

	FIXME

	. There's a bug (in the triangulation?) if the image size > 704 x 704.
*/

#include <assert.h>
#include <locale.h> // setlocale
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h> // getenv
#include <string.h>
#include <time.h>

//------------------------------------------------------------------------------

/** Image area/dimensions. */
#define MIN_WIDTH			32
#define DEFAULT_WIDTH		256
#define MAX_WIDTH			4096
#define MIN_HEIGHT			32
#define DEFAULT_HEIGHT		256
#define MAX_HEIGHT			4096
#define BORDER_PERCENT		2
/** Range for number of points. */
#define INDEX_BITS			16LLU
#define MAX_POINTS			((int)(((1LLU<<INDEX_BITS)+6LLU)/3LLU))
#define MIN_POINTS			3
#define DEFAULT_POINTS		42
/** Option bits. */
#define OPT_MANDATORY		0x01 /** Option must be given. */
#define OPT_ARGS			0x02 /** Option needs an argument. */
#define OPT_BREV			0x04 /** Can use abbreviation. */
#define OPT_DEFAULT			0x08 /** Option has a default value. */
#define OPT_STRING			0x10 /** Expecting a string. */
#define OPT_DECIMAL			0x20 /** Expecting a decimal value. */
#define OPT_FLOAT			0x40 /** Expecting a floating point value. */
/** Float comparison. */
#define SLOPE_INF			(1E+38f)
#define EPSILON_TO_ZERO		(1E-7f)
/** Random number generation. */
#define PCG_MULTIPLIER_64	6364136223846793005ULL
#define PCG_INITIALIZER_64	0x4d595df4d0f33173ULL
#define RNG_SANITY_BREAK	128
/** In-circle test results. */
#define OUTSIDE				(-1)	/*! Point is outside the circle. */
#define COCIRCULAR			0		/*! Points are on the same circle. */
#define INSIDE				1		/*! Point is inside the circle. */
/** Triangle winding. */
#define CW					(-1)	/*! Clock-wise. */
#define COLINEAR			0		/*! Points are co-linear. */
#define CCW					1		/*! Counter-clock-wise. */
/** Comparison. */
#define SMALLER				(-1)
#define EQUAL				0
#define GREATER				1
/** Edge flags. */
#define EDGE_USED			0x01 /*! 1 = edge in use. */
#define EDGE_VERTEX			0x02 /*! 1 = edge is a vertex, 0 = face. */
#define EDGE_VISITED		0x04 /*! 1 = edge visited. */
/** Image option flags. */
#define DRAW_VORONOI		0x01 /*! Draw Voronoi diagram in the image. */
#define DRAW_DELAUNAY		0x02 /*! Draw Delaunay diagram in the image. */
#define DRAW_CONVEX_HULL	0x04 /*! Draw convex hull in the image. */
/** Targa file format. */
#define TGA_UNCOMPRESSED_RGB	2
#define TGA_UNCOMPRESSED_BW		3

//------------------------------------------------------------------------------

#define max(x,y)			((x)>(y)?(x):(y))
#define min(x,y)			((x)<(y)?(x):(y))
#define clamp(x,lo,hi)		max((lo),min((x),(hi)))
#define array_size(x)		(sizeof((x))/sizeof((x)[0]))

//------------------------------------------------------------------------------

enum {
	/** Command line options. */
	OPT_INPUT_FILE,		// --input-file		-i file-name
	OPT_WIDTH,			// --width			-W number
	OPT_HEIGHT,			// --height			-H number
	OPT_NUM_POINTS,		// --points			-P number
	OPT_DRAW_VORONOI,	// --draw-voronoi 	-V
	OPT_DRAW_DELAUNAY,	// --draw-delaunay 	-D
	OPT_DRAW_HULL,		// --draw-convex 	-C
	OPT_VERBOSE,		// --verbose 		-v
	OPT_HELP,			// --help			-h
	OPT_MAX
};

enum {
	/** Timers. */
	T_DELAUNAY,
	T_DRAW_DELAUNAY,
	T_GRAHAM,
	T_DRAW_HULL,
	T_MAX
};

//------------------------------------------------------------------------------

typedef unsigned char u8;
typedef signed short int s16;
typedef unsigned short int u16;
typedef unsigned int u32;
typedef unsigned long int lui;
#if __WORD_SIZE == 64
typedef unsigned long int u64;
#else
typedef unsigned long long int u64;
#endif

typedef s16 coord_t;

//------------------------------------------------------------------------------

typedef union {
	int decimal;
	float floating_point;
	const char *string;
} Value; // 8 B

typedef struct {
	u16 seen;
	u16 id;
	u16 argc;
	u16 flags;
	Value value;
} Option; // 16 B

typedef struct {
	coord_t x;
	coord_t y;
} Point; // 4 B

typedef struct {
	Point from;
	Point to;
} Line_Segment; // 8 B

typedef struct {
	Point a, b, c;
} Triangle; // 12 B

typedef struct {
	int size;
	int capacity;
	Point *points;
} Point_Stack; // 16 B

typedef union {
	int face;
	Point vertex;
} Edge_Data; // 4 B

typedef struct Edge Edge;

struct Edge { // Could be expanded with i.e. color or normal.
	u16 id;			/*! Index in the parent quad-edge. */
	u16 rack;		/*! For storing bits. */
	Edge_Data data;	/*! Vertex or face. */
	Edge *next;		/*! Next CCW quad-edge with the same origin. */
}; // 16 B

typedef struct {
	/*! Edges: primal, ccw rotated, symmetric and cw rotated. */
	Edge e[ 4 ];
} Quad_Edge; // 64 B, fits into a cache line.

typedef struct {
	Edge *left;
	Edge *right;
} Edge_Pair; // 16 B

typedef struct {
	int point_count;				/*! Number of points to triangulate. */
	int last_free_index;			/*! Track highest index in use. */
	int free_index_count;			/*! Count free indices. */
	int edge_cap;					/*! Computed max. #no of edges. */
	int segment_count;				/*! Count line segments. */
	int pad2;
	Quad_Edge *quad_edges;			/*! List of quad-edges. */
	int *free_indices;				/*! Free list of quad-edge indices. */
	Point *points;					/*! List of points. */
	Point_Stack *convex_hull;		/*! Points along the convex hull. */
	Line_Segment *line_segments;	/*! Line sgments of the triangulation. */
} Delaunay; // 64 B so it fits in a cache line.

typedef struct {
	u8 r, g, b;
} RGB; // 3 B

typedef struct {
	char id_length;
	char color_map_type;
	char data_type_code;
	s16 color_map_origin_x;
	s16 color_map_origin_y;
	char color_map_depth;
	s16 x_origin;
	s16 y_origin;
	u16 width;
	u16 height;
	char bits_per_pixel;
	char img_descriptor;
} __attribute__((__packed__)) TGA_Header; // 18 B

typedef struct {
	int width;
	int height;
	int flags;
	int pad0;
	Point *line;
	RGB *data;
} Image; // 32 B

//------------------------------------------------------------------------------

/** Diagram colors. */
static const RGB POINT_COLOR		= { .r = 250U, .g = 255U, .b = 255U };
static const RGB DELAUNAY_COLOR		= { .r =  15U, .g = 255U, .b =  15U };
//static const RGB VORONOI_COLOR		= { .r =  15U, .g =  15U, .b = 255U };
static const RGB CONVEX_HULL_COLOR  = { .r = 127U, .g = 127U, .b = 127U };

/** Tables for accessing members of a quad edge. */
static const int tbl_sym[ ] = {  2,  2, -2, -2 };
static const int tbl_rot[ ] = {  1,  1,  1, -3 };
static const int tbl_inv[ ] = {  3, -1, -1, -1 };

/** Long option format. */
static const char S_OPT_LONG[ ][ 16 ] = { // Must not be empty!
	[ 0 ] = "--input-file",
	[ 1 ] = "--width",
	[ 2 ] = "--height",
	[ 3 ] = "--points",
	[ 4 ] = "--draw-voronoi",
	[ 5 ] = "--draw-delaunay",
	[ 6 ] = "--draw-convex",
	[ 7 ] = "--verbose",
	[ 8 ] = "--help"
};

/** Abbreviated option format. */
static const char S_OPT_BREV[ ][ 4 ] = {
	[ 0 ] = "-i",
	[ 1 ] = "-w",
	[ 2 ] = "-h",
	[ 3 ] = "-p",
	[ 4 ] = "-V",
	[ 5 ] = "-D",
	[ 6 ] = "-C",
	[ 7 ] = "-v",
	[ 8 ] = ""
};

/** Option descriptions. */
static const char S_OPT_DESCR[ ][ 32 ] = {
	[ 0 ] = "load data from file (not impl.)",
	[ 1 ] = "image width in pixels",
	[ 2 ] = "image height in pixels",
	[ 3 ] = "point count",
	[ 4 ] = "draw Voronoi outlines",
	[ 5 ] = "draw triangulation",
	[ 6 ] = "draw convex hull",
	[ 7 ] = "be verbose",
	[ 8 ] = "display this help and exit"
};

/** Helper functions for printing debug info. */

#ifdef DEBUG_VERSION

static Point qe_org( Edge *e );
static Point qe_dst( Edge *e );

void print_edge( const char *title, Edge *e ) {
	Point b = qe_dst( e );
	Point a = qe_org( e );

	printf( "\t%-16s (%3d,%3d) to (%3d,%3d)\n",
	title, a.x, a.y, b.x, b.y );
}

void print_edge_pair( Edge_Pair ep ) {
	if( ep.left ) {
		print_edge( "ep left", ep.left );
	}
	if( ep.right ) {
		print_edge( "ep right", ep.right );
	}
}

void print_point( Point p ) {
	printf( "\t  point (%3d,%3d)\n", p.x, p.y );
}

void print_triangle( const Triangle *t ) {
	printf( "\t  (%3d,%3d), (%3d,%3d), (%3d,%3d)\n",
		t->a.x, t->a.y, t->b.x, t->b.y, t->c.x, t->c.y );
}

void print_qe( const Quad_Edge *qe ) {
	const Edge *e;

	e = &qe->e[ 0 ];
	printf( "\t  # %d org (%3d,%3d), next %p\n",
		e->id, e->data.vertex.x, e->data.vertex.y, e->next );
	e = &qe->e[ 2 ];
	printf( "\t  # %d sym (%3d,%3d), next %p\n",
		e->id, e->data.vertex.x, e->data.vertex.y, e->next );
}

void print_points( Point *points, int num_points ) {
	int i;
	Point p;
	printf( "\n\tListing %d points:\n", num_points );

	for( i = 0; i < num_points; i++ ) {
		p = points[ i ];
		printf( "\t  (%3d,%3d)\n", p.x, p.y );
	}
}

void debug_sizes( void ) {
	printf( "\n\tStructure size in bytes:\n" );
	printf( "\t  TGA_Header:   \t %3lu\n", sizeof( TGA_Header ) );
	printf( "\t  Image:        \t %3lu\n", sizeof( Image ) );
	printf( "\t  Option:       \t %3lu\n", sizeof( Option ) );
	printf( "\t  Point:        \t %3lu\n", sizeof( Point ) );
	printf( "\t  Edge:         \t %3lu\n", sizeof( Edge ) );
	printf( "\t  Edge_Pair:    \t %3lu\n", sizeof( Edge_Pair ) );
	printf( "\t  Quad_Edge:    \t %3lu\n", sizeof( Quad_Edge ) );
	printf( "\t  Line_Segment: \t %3lu\n", sizeof( Line_Segment ) );
	printf( "\t  Point_Stack:  \t %3lu\n", sizeof( Point_Stack ) );
	printf( "\t  Delaunay:     \t %3lu\n", sizeof( Delaunay ) );
}

#endif /* DEBUG_VERSION */

//------------------------------------------------------------------------------

static double milliseconds( void ) {
	struct timespec time;
	clock_gettime( CLOCK_MONOTONIC_RAW, &time );
	return ( double )( time.tv_sec * 1000 )
	   + ( double ) time.tv_nsec / 1000000.0;
}

inline static void zero_buffer( void *buffer, lui size ) {
	u8 *dst = ( u8 * ) buffer;
	u8 const *end = dst + size;

	while( dst < end ) {
		*dst++ = 0;
	}
}

inline static int slen( const char *s ) {
	// Returns the string length.
	const char *cur = s;

	while( cur && ( 0 != *cur ) ) {
		cur++;
	}
	return ( int )( cur - s );
}

inline static int sncmp( const char *left, const char *right, int len ) {
	// Tells wether two strings are equal (0) or not (-1).
	int i;
	char cl, cr;

	for( i = 0; i < len; i++ ) {
		cl = *left;
		cr = *right;

		if( ( cl != cr ) || ( 0 == cl ) || ( 0 == cr ) ) {
			return -1;
		}
		left++;
		right++;
	}
	return 0;
}

inline static
float remap( float x, float t1, float t2, float s1, float s2 ) {
	// Re-map x from the range [t1, t2] to the range [s1, s2].
	return ( ( x - t1 ) / ( t2 - t1 ) ) * ( s2 - s1 ) + s1;
}

inline static float absf( float g ) {
	return ( g < 0.0f ) ? -g : g;
}

inline static int chk_fp_zero( float g ) {
	return ( absf( g ) >= EPSILON_TO_ZERO ) ? 0 : 1;
}

inline static int chk_fp_equal( float f, float g ) {
	int a = ( f < 0.0f );
	int b = ( g < 0.0f );

	if( a ^ b ) { // Values of opposite signs.
		return 0;
	} else {
		return chk_fp_zero( f - g );
	}
}
//------------------------------------------------------------------------------

inline static u32 pcg_rotate_cw_32( u32 value, u32 rotation ) {
	return ( value >> rotation ) | ( value << ( ( -rotation ) & 31 ) );
}

inline static u32 pcg_random( u64 *state ) {
	u64 old_state = *state;
	*state *= PCG_MULTIPLIER_64;
	return pcg_rotate_cw_32(
		( u32 )( ( old_state ^ ( old_state >> 18U ) ) >> 27U ),
		( u32 )( old_state >> 59U ) );
}

static Point random_point( u64 *state, int width, int height ) {
	Point result;
	result.x = ( coord_t )( pcg_random( state ) % width );
	result.y = ( coord_t )( pcg_random( state ) % height );
	return result;
}
//------------------------------------------------------------------------------

inline static int delaunay_triangles( int n, int k ) {
	// Returns #of triangles for set of n points,
	// with k points on the convex hull.
	return 2 * n - 2 - k;
}

inline static int delaunay_edges( int n, int k ) {
	// Returns #of edges for set of n points,
	// with k points on the convex hull.
	return 3 * n - 3 - k;
}

//------------------------------------------------------------------------------

inline static int point_flatten( Point p, int dim ) {
	// Returns an index into a flat array derived from 2d point.
	return( p.y * dim + p.x );
}

inline static Point point_unflatten( int index, int dim ) {
	// Returns a 2d point from a flat array index.
	return ( Point ) { .x = index % dim, .y = index / dim };
}

inline static int same_point( Point a, Point b ) {
	return( ( a.x == b.x ) && ( a.y == b.y ) );
}

inline static Point point_sub( Point a, Point b ) {
	// Returns a - b, that is, a vector from b to a.
	return ( Point ){ .x = a.x - b.x, .y = a.y - b.y };
}

inline static int sign( int n ) {
	// Returns the sign of a given integer.
	return ( n < 0 ) ? -1 : 1;
}

inline static int cross_product( Point a, Point b ) {
	return( a.x * b.y - a.y * b.x );
}

inline static int inner_product( Point a, Point b ) {
	return a.x * b.x + a.y * b.y;
}

inline static int det( Point a, Point b, Point c ) {
	// Computes twice the area of the oriented triangle (a, b, c).
	return( ( b.x - a.x ) * ( c.y - a.y ) - ( b.y - a.y ) * ( c.x - a.x ) );
}

inline static int distance_squared( Point a, Point b ) {
	Point c = point_sub( a, b );
	return( c.x * c.x + c.y * c.y );
}

inline static int length_squared( Point a ) {
	return( a.x * a.x + a.y * a.y );
}

inline static int is_ccw( Point p1, Point p2, Point p3 ) {
	// Computes the winding of a triangle.
	int d = det( p1, p2, p3 );
	return( ( d > 0 ) ? CCW : ( d < 0 ) ? CW : COLINEAR );
}

inline static int in_circle( Point a, Point b, Point c, Point d ) {
	// Checks the relation of point d to the circum-circle of
	// the triangle (a, b, c) via projection onto a paraboloid.
	int result =
		length_squared( a ) * det( b, c, d ) -
		length_squared( b ) * det( a, c, d ) +
		length_squared( c ) * det( a, b, d ) -
		length_squared( d ) * det( a, b, c );
//	printf( "\tin_circle: % d\n", result );
	return( ( result < 0 ) ? OUTSIDE : ( result > 0 ) ? INSIDE : COCIRCULAR );
}

/**
	Quad-edge methods
______________________________________________________________________________*/

inline static Quad_Edge *qe_primal( Edge *e ) {
	return ( Quad_Edge * )( e - e->id );
}

inline static Edge *qe_sym( Edge *e ) {
	return e + tbl_sym[ e->id ];
}

inline static Edge *qe_rot( Edge *e ) {
	return e + tbl_rot[ e->id ];
}

inline static Edge *qe_inv( Edge *e ) {
	return e + tbl_inv[ e->id ];
}

inline static Edge *qe_onext( Edge *e ) {
	return e->next;
}

inline static Edge *qe_lnext( Edge *e ) {
	return qe_rot( qe_onext( qe_inv( e ) ) );
}

inline static Edge *qe_rnext( Edge *e ) {
	return qe_inv( qe_onext( qe_rot( e ) ) );
}

inline static Edge *qe_dnext( Edge *e ) {
	return qe_sym( qe_onext( qe_sym( e ) ) );
}

inline static Edge *qe_oprev( Edge *e ) {
	return qe_rot( qe_onext( qe_rot( e ) ) );
}

inline static Edge *qe_lprev( Edge *e ) {
	return qe_sym( qe_onext( e ) );
}

inline static Edge *qe_rprev( Edge *e ) {
	return qe_onext( qe_sym( e ) );
}

inline static Edge *qe_dprev( Edge *e ) {
	return qe_inv( qe_onext( qe_inv( e ) ) );
}

inline static Point qe_org( Edge *e ) {
	return e->data.vertex;
}

inline static Point qe_dst( Edge *e ) {
	return qe_sym( e )->data.vertex;
}

inline static void qe_set_org( Edge *e, Point p ) {
	e->data.vertex = p;
}

inline static void qe_set_dst( Edge *e, Point p ) {
	qe_sym( e )->data.vertex = p;
}

inline static int qe_right_of( Point p, Edge *e ) {
	return( CCW == is_ccw( p, qe_dst( e ), qe_org( e ) ) );
}

inline static int qe_left_of( Point p, Edge *e ) {
	return( CCW == is_ccw( p, qe_org( e ), qe_dst( e ) ) );
}

inline static int qe_valid( Edge *e, Edge *base ) {
	return qe_right_of( qe_dst( e ), base );
}

inline static void qe_init( Quad_Edge *qe ) {
	// Initialize the members of a quad-edge.
	qe->e[ 0 ].id = 0; qe->e[ 0 ].next = &qe->e[ 0 ];
	qe->e[ 0 ].rack = EDGE_USED | EDGE_VERTEX;
	qe->e[ 1 ].id = 1; qe->e[ 1 ].next = &qe->e[ 3 ];
	qe->e[ 1 ].rack = EDGE_USED;
	qe->e[ 2 ].id = 2; qe->e[ 2 ].next = &qe->e[ 2 ];
	qe->e[ 2 ].rack = EDGE_USED | EDGE_VERTEX;
	qe->e[ 3 ].id = 3; qe->e[ 3 ].next = &qe->e[ 1 ];
	qe->e[ 3 ].rack = EDGE_USED;
}

inline static Edge *qe_make_edge( Delaunay *delaunay, Point a, Point b ) {
	// Returns a pointer to the first edge of the next free quad-edge.
	int index;

	if( delaunay->free_index_count > 0 ) {
		index = delaunay->free_indices[ --delaunay->free_index_count ];
	} else if( delaunay->last_free_index < delaunay->edge_cap ) {
		index = delaunay->last_free_index++;
	} else {
		return 0;
	}
	Quad_Edge *result = &delaunay->quad_edges[ index ];
	assert( ! ( EDGE_USED & result->e[ 0 ].rack ) );

	qe_init( result );
	qe_set_org( result->e, a );
	qe_set_dst( result->e, b );
	return result->e;
}

inline static void qe_splice( Edge *a, Edge * b ) {
	// Separates or joins two rings at the given edges.
	// The function is its own inverse.
	Edge *alpha = qe_rot( qe_onext( a ) );
	Edge *beta = qe_rot( qe_onext( b ) );
	Edge *temp = qe_onext( alpha );
	alpha->next = qe_onext( beta );
	beta->next = temp;
	temp = qe_onext( a );
	a->next = qe_onext( b );
	b->next = temp;
}

inline static void qe_swap( Edge *e ) {
	// Swap an edge in a quadrilateral. */
	Edge *a = qe_oprev( e );
	Edge *b = qe_oprev( qe_sym ( e ) );
	// The first two splice operations disconnect the edge.
	qe_splice( e, a );
	qe_splice( qe_sym( e ), b );
	// The next two splice operations connect the edge again.
	qe_splice( e, qe_lnext( a ) );
	qe_splice( qe_sym( e ), qe_lnext( b ) );
	qe_set_org( e, qe_dst( a ) );
	qe_set_dst( e, qe_dst( b ) );
}

inline static void qe_remove_edge( Delaunay *delaunay, Edge *e ) {
	// Remove a quad-edge.
	int index;
	qe_splice( e, qe_oprev( e ) );
	qe_splice( qe_sym( e ), qe_oprev( qe_sym( e ) ) );
	Quad_Edge *qe = qe_primal( e );

	index = qe - delaunay->quad_edges;
	assert( delaunay->last_free_index > index );
	// FIXME Setting the rack to zero ought to be sufficient,
	//		 but keep zero-ing until code is bug free, because
	//		 a point coordinates of (0,0) will be spotted.
	qe->e[ 0 ].rack = 0;
//	memset( qe, 0, sizeof( Quad_Edge ) );

	if( 1 == ( delaunay->last_free_index - index ) ) {
		// If the last element is removed, then there's
		// no need to add that index to the free list.
		delaunay->last_free_index--;
	} else {
		assert( delaunay->free_index_count < delaunay->edge_cap );
		delaunay->free_indices[ delaunay->free_index_count++ ] = index;
	}
}

inline static Edge *qe_connect( Delaunay *delaunay, Edge *a, Edge *b ) {
	// Connect two edges.
	Edge *e = qe_make_edge( delaunay, qe_dst( a ), qe_org( b ) );
	assert( e );
	qe_splice( e, qe_lnext( a ) );
	qe_splice( qe_sym( e ), b );
	return e;
}

inline static void make_pair( Edge_Pair *ep, Edge *l, Edge *r ) {
	// Combine two edges into a pair.
	ep->left = l;
	ep->right = r;
}

static Edge_Pair triangulate( Delaunay *delaunay, int start, int end ) {
	// Assumes points come sorted by X-first.
	Point org, dst, p, *points;
	Edge_Pair result, L, R;
	Edge *a, *b, *c, *ldo, *rdo, *ldi, *rdi, *base, *lcand, *rcand, *t;
	int middle, test, len;

	len = end - start + 1;
	assert( len > 1 );
	points = delaunay->points;

	if( 2 == len ) {
		a = qe_make_edge( delaunay, points[ start ], points[ start + 1 ] );
		make_pair( &result, a, qe_sym( a ) );
	} else if( 3 == len ) {
		a = qe_make_edge( delaunay, points[ start     ], points[ start + 1 ] );
		b = qe_make_edge( delaunay, points[ start + 1 ], points[ start + 2 ] );
		qe_splice( qe_sym( a ), b );
		test = is_ccw( points[ start ], points[ start + 1 ], points[ start + 2 ] );

		if( COLINEAR == test ) {
			make_pair( &result, a, qe_sym( b ) );
		} else {
			c = qe_connect( delaunay, b, a );

			if( CCW == test ) {
				make_pair( &result, a, qe_sym( b ) );
			} else {
				make_pair( &result, qe_sym( c ), c );
			}
		}
	} else { // len >= 4
		middle = ( start + end ) / 2;
		L = triangulate( delaunay, start, middle );
		R = triangulate( delaunay, middle + 1, end );
		ldo = L.left;	ldi = L.right; // Left outer-inner
		rdi = R.left;	rdo = R.right; // Right inner-outer

		while( 1 ) { // Compute the lower common tangent of L and R.
			if( qe_left_of( qe_org( rdi ), ldi ) ) {
				ldi = qe_lnext( ldi );
			} else if( qe_right_of( qe_org( ldi ), rdi ) ) {
				rdi = qe_rprev( rdi );
			} else {
				break;
			}
		}
		base = qe_connect( delaunay, qe_sym( rdi ), ldi );

		if( same_point( ldi->data.vertex, ldo->data.vertex ) ) {
			ldo = qe_sym( base );
		}
		if( same_point( rdi->data.vertex, rdo->data.vertex ) ) {
			rdo = base;
		}
		while( 1 ) { // Merge loop.
			org = qe_org( base );
			dst = qe_dst( base );
			lcand = qe_onext( qe_sym( base ) );

			if( qe_valid( lcand, base ) ) {
				while( 1 ) {
					p = qe_dst( qe_onext( lcand ) );

					if( same_point( org, p ) ||
						( INSIDE != in_circle( dst, org, qe_dst( lcand ), p ) ) )
					{
						break;
					}
					t = qe_onext( lcand );
					qe_remove_edge( delaunay, lcand );
					lcand = t;
				}
			}
			rcand = qe_oprev( base );

			if( qe_valid( rcand, base ) ) {
				while( 1 ) {
					p = qe_dst( qe_oprev( rcand ) );

					if( same_point( dst, p ) ||
						( INSIDE != in_circle( dst, org, qe_dst( rcand ), p ) ) )
					{
						break;
					}
					t = qe_oprev( rcand );
					qe_remove_edge( delaunay, rcand );
					rcand = t;
				}
			}
			if( ! qe_valid( lcand, base ) && ! qe_valid( rcand, base ) ) {
				break;
			}
			if( ! qe_valid( lcand, base ) ||
				( qe_valid( rcand, base ) &&
				( INSIDE == in_circle(
					qe_dst( lcand ),
					qe_org( lcand ),
					qe_org( rcand ),
					qe_dst( rcand ) ) ) ) )
			{
				base = qe_connect( delaunay, rcand, qe_sym( base ) );
			} else {
				base = qe_connect( delaunay, qe_sym( base ), qe_sym( lcand ) );
			}
		}
		make_pair( &result, ldo, rdo );
	}
	return result;
}
//------------------------------------------------------------------------------

inline static int bresenham( Point *line, Point a, Point b ) {
	// Computes points for a line from a to b and returns the point count.
	int dx_err = 0;
	int dy_err = 0;
	int result = 0;
	int dx = b.x - a.x;
	int dy = b.y - a.y;
	int two_dx = dx + dx;
	int two_dy = dy + dy;
	int cur_x = a.x; // Start at a and move towards b.
	int cur_y = a.y;
	int x_inc = 1;
	int y_inc = 1;

	if( dx < 0 ) {
		x_inc = -1;
		dx = -dx;
		two_dx = -two_dx;
	}
	if( dy < 0 ) {
		y_inc = -1;
		dy = -dy;
		two_dy = -two_dy;
	}
	line[ result++ ] = a; // Add first point.

	if( ( dx == 0 ) && ( 0 == dy ) ) { // No more points on the line.
		return result;
	}
	if( dy <= dx ) { // Slope <= 1
		while( cur_x != b.x ) {
			cur_x += x_inc;
			dx_err += two_dy;

			if( dx_err > dx ) {
				cur_y += y_inc;
				dx_err -= two_dx;
			}
			line[ result++ ] = ( Point ){ .x = cur_x, .y = cur_y };
		}
	} else { // Slope is large, reverse roles of X and Y.
		while( cur_y != b.y ) {
			cur_y += y_inc;
			dy_err += two_dx;

			if( dy_err > dy ) {
				cur_x += x_inc;
				dy_err -= two_dy;
			}
			line[ result++ ] = ( Point ){ .x = cur_x, .y = cur_y };
		}
	}
	return result;
}

static int draw_line( Image *image, Point a, Point b, RGB color ) {
	int i, j;
	int count = bresenham( image->line, a, b );

	for( i = 0; i < count; i++ ) {
		j = point_flatten( image->line[ i ], image->width );
		// Sanity test for an out of bounds line.
		assert( j < image->width * image->height );
		image->data[ j ] = color;
	}
	return count;
}

static void draw_convex_hull( Delaunay *delaunay, Image *image ) {
	int i, j, size;
	Point a, b;
	Point_Stack *stack = delaunay->convex_hull;
	size = stack->size;

	for( i = 0; i < size; i++ ) {
		a = stack->points[ i ];

		if( 1 == ( size - i ) ) {
			b = stack->points[ 0 ];
		} else {
			b = stack->points[ i + 1 ];
		}
		draw_line( image, a, b, CONVEX_HULL_COLOR );
	}
	// Points / vertices.
	for( i = 0; i < delaunay->point_count; i++ ) {
		a = delaunay->points[ i ];
		j = point_flatten( a, image->width );
		image->data[ j ] = POINT_COLOR;

		if( a.x > 0 ) {
			image->data[ j - 1 ] = POINT_COLOR;
		}
		if( a.x < ( ( image->width - 1 ) ) ) {
			image->data[ j + 1 ] = POINT_COLOR;
		}
		if( a.y > 0 ) {
			image->data[ j - image->width ] = POINT_COLOR;
		}
		if( a.y < ( image->height - 1 ) ) {
			image->data[ j + image->width ] = POINT_COLOR;
		}
	}
}

static void delaunay_line_segments( Delaunay *delaunay ) {
	int i, j, dupl_count, qe_count, found, s1, s2, cap;
	Point q, r;
	Quad_Edge *qe;
	Line_Segment *ls;
	dupl_count = qe_count = 0;
	cap = delaunay->edge_cap;

	for( i = 0; i < cap; i++ ) {
		qe = &delaunay->quad_edges[ i ];

		if( EDGE_USED & qe->e[ 0 ].rack ) {
			q = qe_org( &qe->e[ 0 ] );
			r = qe_dst( &qe->e[ 0 ] );
//			assert( ! same_point( q, r ) );

			if( same_point( q, r ) ) {
				continue;
			}
			qe_count++;
			found = 0;

			for( j = 0; j < delaunay->segment_count; j++ ) {
				ls = &delaunay->line_segments[ j ];

				s1 = ( same_point( q, ls->from ) && same_point( r, ls->to ) );
				s2 = ( same_point( r, ls->from ) && same_point( q, ls->to ) );

				if( s1 || s2 ) {
					dupl_count++;
					found = 1;
					break;
				}
			}
			if( ! found ) {
				ls = &delaunay->line_segments[ delaunay->segment_count++ ];
				ls->from = q;
				ls->to = r;
			}
		}
	}
	printf( "\n\tVisited %'d of %'d quad-edges, got %d duplicate segments.\n",
		qe_count, cap, dupl_count );
}

inline static void add_edges(
	Point *points, int *num_points,
	Edge **edges, int *num_edges,
	Edge *e )
{
	int edge_count, point_count;
	Point p;
	Edge *cur, *sym;

	cur = e;
	edge_count = *num_edges;
	point_count = *num_points;

	do {
		cur->rack |= EDGE_VISITED;
		p = qe_org( cur );
		sym = qe_sym( cur );
		edges[ edge_count++ ] = sym;
		points[ point_count++ ] = p;
		cur = qe_lnext( cur );
	} while( cur != e );

	*num_points += ( point_count - *num_points );
	*num_edges += ( edge_count - *num_edges );
}

static void create_triangles( Edge_Pair *ep ) {
	int i, cross, edge_count, point_count;
	Point p1, p2, p3;
	Edge *cur, *sym, *e;
	Edge *edges[ 512 ];
	Point points[ 512 ];
	Triangle tri;
	i = edge_count = point_count = 0;
	e = ( Edge * ) qe_primal( ep->left );

	while( 1 ) {
		p3 = qe_dst( qe_onext( e ) );
		p1 = point_sub( qe_dst( e ), p3 );
		p2 = point_sub( qe_org( e ), p3 );
		cross = cross_product( p1, p2 );

		if( cross >= 0) {
			break;
		}
		e = qe_onext( e );
	}
	add_edges( points, &point_count, edges, &edge_count, e );
	point_count = 0; // Clear points.

	while( i < edge_count ) {
		e = edges[ i++ ];

		if( ! ( e->rack & EDGE_VISITED ) ) {
			add_edges( points, &point_count, edges, &edge_count, e );
		}
	}
	printf( "\n\tCreating triangles: found %'d points, yielding %'d triangles\n",
		point_count, point_count / 3 );

	for( i = 0; i < point_count; i += 3 ) {
		tri.a = points[ i ];
		tri.b = points[ i + 1 ];
		tri.c = points[ i + 2 ];
	}
}

static void draw_delaunay( Delaunay *delaunay, Image *image ) {
	int i, j;
	Point a, b;
	Line_Segment *ls;

	delaunay_line_segments( delaunay );

	for( j = 0; j < delaunay->segment_count; j++ ) {
		ls = &delaunay->line_segments[ j ];
		a = ls->from;
		b = ls->to;
		draw_line( image, a, b, DELAUNAY_COLOR );
	}
	// Points / vertices.
	for( i = 0; i < delaunay->point_count; i++ ) {
		a = delaunay->points[ i ];
		j = point_flatten( a, image->width );
		image->data[ j ] = POINT_COLOR;

		if( a.x > 0 ) {
			image->data[ j - 1 ] = POINT_COLOR;
		}
		if( a.x < ( ( image->width - 1 ) ) ) {
			image->data[ j + 1 ] = POINT_COLOR;
		}
		if( a.y > 0 ) {
			image->data[ j - image->width ] = POINT_COLOR;
		}
		if( a.y < ( image->height - 1 ) ) {
			image->data[ j + image->width ] = POINT_COLOR;
		}
	}
}
//------------------------------------------------------------------------------

inline static void swap_points( Point *points, int a, int b ) {
	Point tmp = points[ a ];
	points[ a ] = points[ b ];
	points[ b ] = tmp;
}

inline static int x_first( Point a, Point b ) {
	// Compares two points, first by X, then by Y.
	// Returns -1 if a < b, 0 if a == b, 1 if a > b.
	int result = 0;

	if( a.x < b.x ) {
		result =  -1;
	} else if( a.x > b.x ) {
		result =  1;
	} else if( a.y < b.y ) {
		result =  -1;
	} else if( a.y > b.y ) {
		result =  1;
	}
	return result;
}

inline static int y_first( Point a, Point b ) {
	// Compares two points, first by Y, then by X.
	// Returns -1 if a < b, 0 if a == b, 1 if a > b.
	int result = 0;

	if( a.y < b.y ) {
		result =  -1;
	} else if( a.y > b.y ) {
		result =  1;
	} else if( a.x < b.x ) {
		result =  -1;
	} else if( a.x > b.x ) {
		result =  1;
	}
	return result;
}

typedef int ( *axis_cmp )( Point a, Point b );

void quick_sort_points( const axis_cmp cmp, Point *buffer, int len ) {
	if( len < 2 ) return;
	int i, j;
	// For simplicity the choice of the pivot is not random.
	Point pivot = buffer[ len >> 1 ];

	for( i = 0, j = len - 1; ; i++, j-- ) {
		while( -1 == cmp( buffer[ i ], pivot ) ) {
			i++;
		}
		while( 1 == cmp( buffer[ j ], pivot ) ) {
			j--;
		}
		if( i >= j ) {
			break;
		}
		swap_points( buffer, i, j );
	}
	quick_sort_points( cmp, buffer, i );
	quick_sort_points( cmp, buffer + i, len - i );
}

inline static void stack_push( Point_Stack *stack, Point p ) {
	assert( stack->size < stack->capacity );
	stack->points[ stack->size++ ] = p;
}

inline static void stack_pop( Point_Stack *stack ) {
	assert( stack->size > 0 );

	if( stack->size > 0 ) {
		stack->size--;
	}
}

inline static Point stack_top( Point_Stack *stack ) {
	// Caller responsible for checking that the stack size is > 0!
	return stack->points[ stack->size - 1 ];
}

inline static Point stack_next( Point_Stack *stack ) {
	// Caller responsible for checking that the stack size is > 1!
	return stack->points[ stack->size - 2 ];
}

inline static float slope( Point a, Point b ) {
	// Returns the slope of a line segment.
	Point c = point_sub( a, b );

	if( 0 == c.x ) {
		return SLOPE_INF;
	}
	return ( float ) c.y / ( float ) c.x;
}

inline static
int slope_cmp( Point p, Point anchor, float m, int dist ) {
	// Notiz: this is about 2 (sometimes up to 4) times faster
	// than comparing polar angles, which uses atan2f().
	int result, d, s_neg, m_neg, s_low;
	float s = slope( p, anchor );

	if( chk_fp_equal( s, m ) ) { // Prefer the one with greater distance.
		d = distance_squared( p, anchor );
		result = ( d < dist ) ? GREATER : ( d > dist ) ? SMALLER : EQUAL;
	} else { // Different in value and/or sign.
		s_neg = ( s < 0.0f );
		m_neg = ( m < 0.0f );
		s_low = ( s < m );

		if( SLOPE_INF == m ) { // s != m, so slope could be undefined.
			result = s_neg ? SMALLER : GREATER;
		} else if( s_neg ) {
			if( m_neg ) {
				result = s_low ? GREATER : SMALLER;
			} else {
				result = SMALLER;
			}
		} else if( m_neg ) { // ! s_neg
			result = GREATER;
		} else {
			result = s_low ? GREATER : SMALLER;
		}
	}
	return result;
}

static void slope_sort( Point *buffer, Point p, int len ) {
	if( len < 2 ) return;
	int i, j;
	Point pivot = buffer[ len >> 1 ]; // Notiz: not random for simplicity.
	int dist = distance_squared( pivot, p );
	float m = slope( pivot, p );

	for( i = 0, j = len - 1; ; i++, j-- ) {
		while( SMALLER == slope_cmp( buffer[ i ], p, m, dist ) ) {
			i++;
		}
		while( GREATER == slope_cmp( buffer[ j ], p, m, dist ) ) {
			j--;
		}
		if( i >= j ) {
			break;
		}
		swap_points( buffer, i, j );
	}
	slope_sort( buffer, p, i );
	slope_sort( buffer + i, p, len - i );
}

static int graham_scan( Delaunay *delaunay ) {
	// Compute the convex hull for a set of points and return the point count.
	int i, smallest;
	Point p, anchor;
	Point *points = delaunay->points;
	Point_Stack *stack = delaunay->convex_hull;
	const int capacity = delaunay->point_count;
	// Step 1: find the point(s) with the smallest Y,
	// of those chose the one with the smallest X.
	smallest = 0;
	anchor = points[ 0 ];

	for( i = 1; i < capacity; i++ ) {
		p = points[ i ];

		if( SMALLER == y_first( p, anchor ) ) { // p < anchor
			anchor = p;
			smallest = i;
		}
	}
	points[ smallest ] = points[ 0 ];
	points[ 0 ] = anchor;
	// Step 2: sort rest of points by the slope of the vector p - anchor.
	slope_sort( points + 1, anchor, capacity - 1 );
	// Step 3: compute hull.
	for( i = 0; i < capacity; i++ ) {
		p = points[ i ];

		while( ( stack->size > 1 ) &&
			( CCW == is_ccw( stack_next( stack ), stack_top( stack ), p ) ) )
		{
			stack_pop( stack );
		}
		stack_push( stack, p );
	}
	return stack->size;
}

static int create_ps_plot( Delaunay *delaunay, char *ps_data, int cap ) {
	int i, num_points, result;
	Point p, *points;
	const float lo = 0.0f;
	const float hi = 1024.0f;
	float x, y;

	result = 0;
	num_points = delaunay->point_count;
	points = delaunay->points;

	result += sprintf( ps_data, "%%!\n" );
	result += sprintf( ps_data + result, "/seg {moveto lineto stroke} def\n" );
	result += sprintf( ps_data + result, "/bullet {0.02 0 360 arc fill} def\n" );
	result += sprintf( ps_data + result, "3 72 mul dup scale\n" );
	result += sprintf( ps_data + result, "1.0 1.5 translate\n" );
	result += sprintf( ps_data + result, "0 setlinecap\n2 setlinejoin\n" );

	for( i = 0; i < num_points; i++ ) {
		p = points[ i ];
		x = remap( ( float ) p.x, lo, hi, 0.0f, 1.0f );
		y = remap( ( float ) p.y, lo, hi, 0.0f, 1.0f );

		result += sprintf( ps_data + result, "/p%d {%f %f} def\n",
			i + 1, x, y );
	}
	result += sprintf( ps_data + result,
	"0.7 setgray\n0.03 setlinewidth\n[ ] 0 setdash\n" );

	/** Draw Voronoi. */
#if 0
	result += sprintf( ps_data + result,
	"0 setgray\n1 2 72 mul div setlinewidth\n[ ] 0 setdash\n" );
#endif
	for( i = 0; i < num_points; i++ ) {
	   result += sprintf( ps_data + result, "p%d bullet\n", i + 1 );
	}
	result += sprintf( ps_data + result, "showpage\n%%%%EndProlog" );
	*( ps_data + result ) = 0;
	return result;
}

//------------------------------------------------------------------------------

static int write_file( const char *file, char *data, int size ) {
	int written;
	FILE *fh = fopen( file, "w" );

	if( ! fh ) {
		fprintf( stderr,
			"\n\tError: could not open file %s for writing!\n", file );
		return -1;
	}
	written = fwrite( data, 1, size, fh );
	fclose( fh );

	if( written != size ) {
		fprintf( stderr, "\n\tError: output file %s corrupt!\n", file );
		return -1;
	} else {
		printf( "\tWrote %'d bytes to %s\n",
		size, file );
	}
	return written;
}

static void write_tga( const char *file, Image *image, char code, char bps ) {
	TGA_Header tga = { };
	tga.data_type_code = code;
	tga.width = image->width;
	tga.height = image->height;
	tga.bits_per_pixel = bps;
	int written;
	const int data_size = image->width * image->height * sizeof( RGB );
	const int header_size = sizeof( TGA_Header );

	FILE *texfh = fopen( file, "w" );

	if( ! texfh ) {
		fprintf( stderr, "\n\tError: could not open %s for writing!", file );
	} else {
		written = fwrite( &tga, 1, header_size, texfh );
		written += fwrite( image->data, 1, data_size, texfh );
		fclose( texfh );

		if( written != ( data_size + header_size ) ) {
			fprintf( stderr, "\n\tError: output file %s corrupt!\n", file );
		} else {
			printf( "\tWrote %'d bytes to %s\n",
			data_size + header_size, file );
		}
	}
}
//------------------------------------------------------------------------------

static void usage( const char *exe, int brief ) {
	int i, offset;
	char buffer[ 1024 ]; // Currently needed bytes: 562

	if( brief ) {
		sprintf( buffer, "Try: %s --help for more information\n", exe );
	} else {
		offset = sprintf( buffer, "Usage: %s [OPTION]...\n", exe );
		offset += sprintf( buffer + offset, "\nOptions are:\n\n" );

		for( i = 0; i < OPT_MAX; i++ ) {
			if( strlen( S_OPT_BREV[ i ] ) > 0 ) {
				offset += sprintf( buffer + offset, "   %s, %-16s   %s\n",
					S_OPT_BREV[ i ], S_OPT_LONG[ i ], S_OPT_DESCR[ i ] );
			} else {
				offset += sprintf( buffer + offset, "       %-16s   %s\n",
					S_OPT_LONG[ i ], S_OPT_DESCR[ i ] );
			}
		}
		offset += sprintf( buffer + offset,
			"\nUnspecified options default to the following:\n" );
		sprintf( buffer + offset,
			"    --width %d --height %d --points %d --draw-delaunay\n",
				DEFAULT_WIDTH, DEFAULT_HEIGHT, DEFAULT_POINTS );
	}
	fprintf( stderr, buffer );
}

inline static
int scan_decimal( char *msg, const char *arg, Value *value, int option ) {
	int result = sscanf( arg, "%d", &value->decimal );

	if( 1 != result ) {
		sprintf( msg,
			"Error: type mismatch for option %s, expected a decimal, got '%s' !",
			S_OPT_LONG[ option ], arg );
	}
	return result;
}

inline static
int scan_float( char *msg, const char *arg, Value *value, int option ) {
	int result = sscanf( arg, "%f", &value->floating_point );

	if( 1 != result ) {
		sprintf( msg,
			"Error: type mismatch for option %s, expected a float, got '%s' !",
			S_OPT_LONG[ option ], arg );
	}
	return result;
}

static int get_options( int argc, char *argv[ ], Option *options ) {
	u16 id;
	int i, j, count, done, exp, found, len, len_long, len_brev;
	int mandatory, exp_args, error;
	Option *opt;
	const char *arg;
	char msg[ 128 ];
	count = exp = done = 0;

	for( i = 0; i < OPT_MAX; i++ ) {
		opt = &options[ i ];
		mandatory = ( opt->flags & OPT_MANDATORY );

		if( mandatory ) {
			exp++;
		}
	}
	i = 1;

	while( i < argc ) {
		arg = argv[ i ];
		len = strlen( arg );
		found = 0;

		for( j = 0; j < OPT_MAX; j++ ) {
			opt = &options[ j ];
			id = opt->id;

			len_long = strlen( S_OPT_LONG[ id ] );
			len_brev = strlen( S_OPT_BREV[ id ] );

			if( len_long == len ) {
				if( 0 == strncmp( S_OPT_LONG[ id ], arg, len ) ) {
					found = 1;
					break;
				}
			} else if( len_brev == len ) {
				if( 0 == strncmp( S_OPT_BREV[ id ], arg, len ) ) {
					found = 1;
					break;
				}
			}
		}
		if( found ) {
			if( 1 == opt->seen ) {
				sprintf( msg, "Error: duplicate option %s!", S_OPT_LONG[ id ] );
				goto lbl_err;
			} else {
				opt->seen = 1;
			}
			exp_args = ( opt->flags & OPT_ARGS );

			if( exp_args ) {
				if( 1 == ( argc - i ) ) { // Error. No more args.
					sprintf( msg, "Error: option %s requires a value!", arg );
					goto lbl_err;
				}
				i++;

				if( opt->flags & OPT_DECIMAL ) {
					error = scan_decimal( msg, argv[ i ], &opt->value, id );
				} else if( opt->flags & OPT_FLOAT ) {
					error = scan_float( msg, argv[ i ], &opt->value, id );
				} else if( opt->flags & OPT_STRING ) {
					opt->value.string = argv[ i ];
					error = 1;
				}
				if( error < 1 ) {
					goto lbl_err;
				}
			}
			if( opt->flags & OPT_MANDATORY ) {
				count++;
			}
		} else {
			sprintf( msg, "Error: option %s not recognized!", arg );
			goto lbl_err;
		}
		i++;
	}
	if( exp != count ) {
		sprintf( msg, "Error: expected %d arguments!", exp );
		goto lbl_err;
	}
	return 0;

lbl_err:
	fprintf( stderr, "%s\n", msg );
	return -1;
}

static void make_option( Option *options, int id, int argc, int flags ) {
	Option *opt = &options[ id ];

	opt->seen = 0;
	opt->id = id;
	opt->argc = argc;
	opt->value.string = 0;

	if( strlen( S_OPT_BREV[ opt->id ] ) > 0 ) {
		flags |= OPT_BREV;
	}
	opt->flags = flags;
}

//------------------------------------------------------------------------------

#ifdef DEBUG_VERSION

void test_functions( void ) {
	const int dim = 4096;
	int idx, i, j;
	Point p, q;

	for( i = 0; i < dim; i++ ) {
		for( j = 0; j < dim; j++ ) {
			p = ( Point ){ .x = i, .y = j };
			idx = point_flatten( p, dim );
			q = point_unflatten( idx, dim );
			assert( same_point( p, q ) );
		}
	}
}

void test_predicates( void ) {
	/** Verify predicates in_circle and is_ccw. */
#define N 4

	const int results[ ] = { INSIDE, COCIRCULAR, OUTSIDE };
	int i, j, k, d, test;
	Point p;
	Triangle *t;
	Triangle ccw[ ] = {
		{ { 0, 0 }, { 8, 0 }, { 0, 8 } },
		{ { 8, 0 }, { 0, 8 }, { 0, 0 } },
		{ { 0, 8 }, { 0, 0 }, { 8, 0 } },
	};
	Point points[ ][ N ] = {
		{
			{ 1, 1 }, { 1, 2 }, { 2, 1 }, { 2, 2 }
		}, {
			{ 0, 0 }, { 8, 0 }, { 0, 8 }, { 8, 8 }
		}, {
			{ 9, 0 }, { 0, 9 }, { 9, 9 }, { -1, -1 }
		}
	};
	t = &ccw[ 0 ];

	for( i = 0; i < array_size( ccw ); i++ ) {
		d = det( t->a, t->b, t->c );
		assert( d > 0 );

		d = det( t->a, t->c, t->b );
		assert( d < 0 );

		for( k = 0; k < array_size( results ); k++ ) {
			for( j = 0; j < N; j++ ) {
			p = points[ k ][ j ];
			test = in_circle( t->a, t->b, t->c, p );
	//		printf( "\n\ttest: i = %d, k = %d, j = %d\n", i, k, j );
			assert( results[ k ] == test );
			}
		}
		t++;
	}
	Edge *e;
	Quad_Edge qe;
	Point above = { .x = 1, .y = 3 };
	Point below = { .x = 1, .y = 1 };
	Point org = { .x = 1, .y = 2 };
	Point dst = { .x = 2, .y = 2 };

	qe_init( &qe );
	e = qe.e;
	qe_set_org( e, org );
	qe_set_dst( e, dst );

	test = qe_left_of( above, e );
	assert( test );
	test = qe_left_of( below, e );
	assert( ! test );
	test = qe_right_of( above, e );
	assert( ! test );
	test = qe_right_of( below, e );
	assert( test );
#undef N
}
#endif /* DEBUG_VERSION */

int main( int argc, char *argv[ ] ) {
	coord_t bx, by;
	lui sz, accum;
	int i, j, k, error, count, line_len;
	int width, height, num_points, tri_cap, edge_cap;
	double t_start, t_end, timings[ T_MAX ];
	Image image;
	u8 *storage;
	Point p;
	Edge_Pair ep;
	Delaunay delaunay;
	u64 entropy;
	char *nl;
	const char *exe, *input_file;
	Option options[ OPT_MAX ];

	width = DEFAULT_WIDTH;
	height = DEFAULT_HEIGHT;
	num_points = DEFAULT_POINTS;
	storage = 0;
	exe = argv[ 0 ];
	input_file = 0;
	entropy = PCG_INITIALIZER_64;
	memset( &image, 0, sizeof( Image ) );
	memset( timings, 0, sizeof( timings ) );

	if( argc > 1 ) {
		make_option( options, OPT_INPUT_FILE,	 1, OPT_ARGS | OPT_STRING );
		make_option( options, OPT_WIDTH,		 1, OPT_ARGS | OPT_DECIMAL );
		make_option( options, OPT_HEIGHT,		 1, OPT_ARGS | OPT_DECIMAL );
		make_option( options, OPT_NUM_POINTS,	 1, OPT_ARGS | OPT_DECIMAL );
		make_option( options, OPT_DRAW_VORONOI,	 0, 0 );
		make_option( options, OPT_DRAW_DELAUNAY, 0, 0 );
		make_option( options, OPT_DRAW_HULL,	 0, 0 );
		make_option( options, OPT_HELP,			 0, 0 );
		make_option( options, OPT_VERBOSE,		 0, 0 );

		error = get_options( argc, argv, options );

		if( 0 != error ) {
			goto lbl_usage_brief;
		}
		if( options[ OPT_HELP ].seen ) {
			goto lbl_usage_long;
		}
		if( options[ OPT_WIDTH ].seen ) {
			width = clamp( options[ OPT_WIDTH ].value.decimal,
			MIN_WIDTH, MAX_WIDTH );
		}
		if( options[ OPT_HEIGHT ].seen ) {
			height = clamp( options[ OPT_HEIGHT ].value.decimal,
			MIN_HEIGHT, MAX_HEIGHT );
		}
		if( options[ OPT_NUM_POINTS ].seen ) {
			num_points = clamp( options[ OPT_NUM_POINTS ].value.decimal,
			MIN_POINTS, ( width * height ) / 2 );
		}
		if( options[ OPT_DRAW_VORONOI ].seen ) {
			image.flags |= DRAW_VORONOI;
		}
		if( options[ OPT_DRAW_DELAUNAY ].seen ) {
			image.flags |= DRAW_DELAUNAY;
		}
		if( options[ OPT_DRAW_HULL ].seen ) {
			image.flags |= DRAW_CONVEX_HULL;
		}
		if( options[ OPT_INPUT_FILE ].seen ) {
			input_file = options[ OPT_INPUT_FILE ].value.string;
		}
	}
	if( 0 == image.flags ) {
		// FIXME To compute and draw the Dealunay will be the default for release.
//		image.flags = DRAW_DELAUNAY;
		image.flags = DRAW_DELAUNAY | DRAW_CONVEX_HULL;
	}
	nl = getenv( "LANG" );

	if( nl ) {
		setlocale( LC_ALL, nl );
	} else {
		setlocale( LC_ALL, "" );
	}
	// TODO When implementing reading points from a file, load it here.
	//		Then the memory allocation scheme needs to be revisited.

	// Notiz: Use the upper bounds for the #of triangles/edges
	//        in order to allocate enough space.
	// In any valid triangulation the convex hull has at least 3 points.
	tri_cap = delaunay_triangles( num_points, 3 );
	edge_cap = delaunay_edges( num_points, 3 );
	line_len = width + height;
	sz = ( lui )( line_len + 2 * num_points ) * sizeof( Point )
		+ ( lui )( width * height ) * sizeof( RGB )
		+ ( lui ) edge_cap * sizeof( Quad_Edge )
		+ ( lui ) edge_cap * sizeof( int )
		+ ( lui ) edge_cap * sizeof( Line_Segment )
		+ sizeof( Point_Stack );

	storage = calloc( 1, sz );

	if( ! storage ) {
		fprintf( stderr,
			"\n\tError: alloction of %'lu bytes of memory failed!\n", sz );
		goto lbl_end;
	}
	// Initialize data structures and pointers.
	accum = 0;

	delaunay.point_count = num_points;
	delaunay.last_free_index = 0;
	delaunay.free_index_count = 0;
	delaunay.edge_cap = edge_cap;
	delaunay.segment_count = 0;
	delaunay.quad_edges = ( Quad_Edge * ) storage;
	accum += ( lui ) edge_cap * sizeof( Quad_Edge );
	delaunay.free_indices = ( int * )( storage + accum );
	accum += ( lui ) edge_cap * sizeof( int );
	delaunay.points = ( Point * )( storage + accum );
	accum += ( lui ) num_points * sizeof( Point );
	delaunay.convex_hull = ( Point_Stack * )( storage + accum );
	delaunay.convex_hull->size = 0;
	delaunay.convex_hull->capacity = num_points;
	accum += sizeof( Point_Stack );
	delaunay.convex_hull->points = ( Point * )( storage + accum );
	accum += ( lui ) num_points * sizeof( Point );
	delaunay.line_segments = ( Line_Segment * )( storage + accum );
	accum += ( lui ) edge_cap * sizeof( Line_Segment );

	image.width = width;
	image.height = height;
	image.line = ( Point * )( storage + accum );
	accum += ( lui ) line_len * sizeof( Point );
	image.data = ( RGB * )( storage + accum );
	accum += ( lui )( width * height ) * sizeof( RGB );

	printf( "\n\tComputed memory size %'10lu, accumulated size: %'10lu\n", sz, accum );

	assert( accum == sz );

	printf( "\n\tMemory usage in bytes: \t %'10lu\n", sz );
	printf( "\t  quad edges:    \t %'10lu\n", ( lui ) edge_cap * sizeof( Quad_Edge ) );
	printf( "\t  free list:     \t %'10lu\n", ( lui ) edge_cap * sizeof( int ) );
	printf( "\t  random points: \t %'10lu\n", ( lui ) num_points * sizeof( Point ) );
	printf( "\t  hull stack:    \t %'10lu\n", sizeof( Point_Stack ) );
	printf( "\t  hull points:   \t %'10lu\n", ( lui ) num_points * sizeof( Point ) );
	printf( "\t  line segments: \t %'10lu\n", ( lui ) edge_cap * sizeof( Line_Segment ) );
	printf( "\t  output line:   \t %'10lu\n", ( lui ) line_len * sizeof( Point ) );
	printf( "\t  output image:  \t %'10lu\n", ( lui )( width * height ) * sizeof( RGB ) );

#ifdef DEBUG_VERSION
	test_functions( );
	test_predicates( );
	debug_sizes( );
#endif /* DEBUG_VERSION */

	// Compute a border of the image area. All points will be situated inside.
	bx = ( coord_t )( ( ( float ) width  / 100.0f ) * ( float ) BORDER_PERCENT );
	by = ( coord_t )( ( ( float ) height / 100.0f ) * ( float ) BORDER_PERCENT );
	accum = 0;

	if( input_file ) {
		// stat file, open for reading
		// read first line, parse options
		// alloc space
		// read points
		// close file
	} else { // Create random points.
		for( i = 0; i < num_points; i++ ) {
			count = 0;

			while( 1 ) {
				count++;
				p = random_point( &entropy, width - 2 * bx, height - 2 * by );
				p.x += bx;
				p.y += by;

				for( j = 0; j < i; j++ ) {
					if( same_point( p, delaunay.points[ j ] ) ) {
						break;
					}
				}
				if( i == j ) { // New point is unique.
					break;
				}
				if( count > RNG_SANITY_BREAK ) {
					fprintf( stderr,
					"\n\tError: random point generation corrupt!\n" );
					goto lbl_end;
				}
			}
			accum += ( lui ) count;
			delaunay.points[ i ] = p;
		}
	}
	printf( "\n\tSetup: %'d points, image size: %d x %d\n",
		num_points, width, height );

	if( DRAW_CONVEX_HULL & image.flags ) { // Draw triangulation.
		t_start = milliseconds( );
		k = graham_scan( &delaunay );
		t_end = milliseconds( );
		timings[ T_GRAHAM ] = t_end - t_start;

		if( k < 0 ) {
			goto lbl_end;
		}
		t_start = milliseconds( );
		draw_convex_hull( &delaunay, &image );
		t_end = milliseconds( );
		timings[ T_DRAW_HULL ] = t_end - t_start;
		write_tga( "delaunay_convex_hull.tga", &image, TGA_UNCOMPRESSED_RGB, 24 );
		memset( image.data, 0, ( lui )( image.width * image.height ) * sizeof( RGB ) );

#ifdef DEBUG_VERSION
		tri_cap = delaunay_triangles( num_points, k );
		edge_cap = delaunay_edges( num_points, k );
		printf( "\t%d points on the convex hull => %d triangle(s), %d edges.\n",
			k, tri_cap, edge_cap );
#endif /* DEBUG_VERSION */

	}
	if( DRAW_DELAUNAY & image.flags ) { // Draw triangulation.
		t_start = milliseconds( );
		quick_sort_points( x_first, delaunay.points, delaunay.point_count );
//		print_points( delaunay.points, delaunay.point_count );
		ep = triangulate( &delaunay, 0, num_points - 1 );
		t_end = milliseconds( );
		timings[ T_DELAUNAY ] = t_end - t_start;

//		create_triangles( &ep );

		t_start = milliseconds( );
		draw_delaunay( &delaunay, &image );
		t_end = milliseconds( );
		timings[ T_DRAW_DELAUNAY ] = t_end - t_start;

//		create_triangles( &ep );
#ifdef DEBUG_VERSION
		printf( "\tLine segment count: %d\n", delaunay.segment_count );
#endif /* DEBUG_VERSION */
		write_tga( "delaunay_triangulation.tga", &image, TGA_UNCOMPRESSED_RGB, 24 );
	}
#if 0
	sz = 1024 * 1024 * 8; // 8 Mb
	char *ps_data;
	ps_data = calloc( 1, sz );

	if( ps_data ) {
	sz = create_ps_plot( &delaunay, ps_data, sz );
	write_file( "delaunay.ps", ps_data, sz );
	free( ps_data );
	} else {
	fprintf( stderr,
		"\n\tError: alloction of %d bytes of memory failed!\n", sz );
	}
#endif

	printf( "\n\tTimings in ms:\n" );

	if( DRAW_CONVEX_HULL & image.flags ) { // Draw convex hull.
		printf( "\t  Graham scan:   %.6lf, output: %.6lf\n",
			timings[ T_GRAHAM ], timings[ T_DRAW_HULL ] );
	}
	if( DRAW_DELAUNAY & image.flags ) { // Draw triangulation.
		printf( "\t  Triangulation: %.6lf, output: %.6lf\n",
			timings[ T_DELAUNAY ], timings[ T_DRAW_DELAUNAY ] );
	}
	goto lbl_end;

lbl_usage_brief:
	usage( exe, 1 );
	goto lbl_end;

lbl_usage_long:
	usage( exe, 0 );

lbl_end:
	if( storage ) {
		free( storage );
	}

	return 0;
}
/**
	End of program. Below is stuff in testing.
*/

#if 0
inline static
Point circum_center( Point a, Point b, Point c ) {
	Point result;
	float d;
	int la, lb, lc, aby, bcy, cay, acx, bax, cbx;

	la = length_squared( a );
	lb = length_squared( b );
	lc = length_squared( c );

	aby = a.y - b.y;
	bcy = b.y - c.y;
	cay = c.y - a.y;

	acx = a.x - c.x;
	bax = b.x - a.x;
	cbx = c.x - b.x;

	d = ( float )( 2 * ( a.x * bcy + b.x * cay + c.x * aby ) );
	result.x = ( int ) roundf( ( float )( la * bcy + lb * cay + lc * aby ) / d );
	result.y = ( int ) roundf( ( float )( la * cbx + lb * acx + lc * bax ) / d );
	return result;
}
#endif
#if 0
inline static float polar_angle( Point a, Point b ) {
	Point c = point_sub( a, b );
	return atan2f( c.y, c.x );
}

inline static
int polar_cmp( Point p, Point anchor, float angle, int dist ) {
	int d, result;
	float f = polar_angle( p, anchor );

	if( chk_fp_equal( f, angle ) ) {
	d = distance_squared( p, anchor );
		result = ( d < dist ) ? GREATER : ( d > dist ) ? SMALLER : EQUAL;
	} else {
	result = ( f < angle ) ? GREATER : SMALLER;
	}
	return result;
}

static void polar_sort( Point *buffer, Point p, int len ) {
	if( len < 2 ) return;
	int i, j;
	Point pivot = buffer[ len >> 1 ]; // Notiz: not random for simplicity.
	int dist = distance_squared( pivot, p );
	float angle = polar_angle( pivot, p );

	for( i = 0, j = len - 1; ; i++, j-- ) {
	while( SMALLER == polar_cmp( buffer[ i ], p, angle, dist ) ) {
		i++;
	}
	while( GREATER == polar_cmp( buffer[ j ], p, angle, dist ) ) {
		j--;
	}
	if( i >= j ) {
		break;
	}
	swap_points( buffer, i, j );
	}
	polar_sort( buffer, p, i );
	polar_sort( buffer + i, p, len - i );
}
#endif
