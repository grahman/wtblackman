//
//  wtblackman~.c
//  simplemsp~
//
//  Created by Graham Barab on 1/11/16.
//
//

#include "wtblackman~.h"
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "ext.h"			// standard Max include, always required (except in Jitter)
#include "ext_obex.h"		// required for "new" style objects
#include "z_dsp.h"			// required for MSP objects

#define SAMPLES (2048)
#define FIX_ZERO(x) (x ? x : 0.000001)

// struct to represent the object's state
typedef struct _wtblackman {
    t_pxobject          ob;			// the object itself (t_pxobject in MSP instead of t_object)
    double R1;                      // Current sample rate
    unsigned long long N;           // Number of samples of original waveform.
    double f;                       // Frequency of oscillator
    double f_1;                     // Previous value of f
    double index;                   // Phase accumulator
} t_wtblackman;

static double *wavetable;
static unsigned long long instances;

// method prototypes
void *t_wtblackman_new(t_symbol *s, long argc, t_atom *argv);
void t_wtblackman_free(t_wtblackman *x);
void t_wtblackman_assist(t_wtblackman *x, void *b, long m, long a, char *s);
void t_wtblackman_float(t_wtblackman *x, double f);
void t_wtblackman_dsp64(t_wtblackman *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void t_wtblackman_perform64(t_wtblackman *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *t_wtblackman_class = NULL;


//***********************************************************************************************


/* ---------------- HELPER FUNCTIONS ---------------- */

int approx_zero(double x)
{
    double thresh = 5E-15;
    return (x <= 0 + thresh && x >= 0 - thresh);
}

/* Taken from http://stackoverflow.com/questions/12276675/modulus-with-negative-numbers-in-c */
static long mod(long a, long b)
{ return (a%b+b)%b; }

void init_wavetable(double *table, unsigned long long N /*,double sample_rate */)
{
    unsigned int i;
    
//    post("init wavetable called");
    for (i = 0; i < N; ++i)
        table[i] = 0.42 - (0.5 * cos(2 * M_PI * i / (double)N))+ (0.08 * cos(4.0 * M_PI * i / (double)N));
    
    /* Cos wavetable oscillator for tests */
//    for (i = 0; i < N; ++i) {
//        table[i] = cos(2.0 * M_PI * i / (double)N);
//    }
}

/* Function interpolate4:
    4-pt Lagrange interpolation implementation */
static inline double interpolate4(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double index)
{
    int i;
    double L[4] = {0};
    double out = 0;
    
    for (i = 0; i < 4; ++i) {
        switch(i) {
            case 0: {
                L[i] = ((index - x1) * (index - x2) * (index - x3)) / ((x0 - x1) * (x0 - x2) * (x0 - x3));
                break;
            }
            case 1: {
                L[i] = ((index - x0) * (index - x2) * (index - x3)) / ((x1 - x0) * (x1 - x2) * (x1 - x3));
                break;
            }
            case 2: {
                L[i] = ((index - x0) * (index - x1) * (index - x3)) / ((x2 - x0) * (x2 - x1) * (x2 - x3));
                break;
            }
            case 3: {
                L[i] = ((index - x0) * (index - x1) * (index - x2)) / ((x3 - x0) * (x3 - x1) * (x3 - x2));
                break;
            }
        }
    }
    
    for (i = 0; i < 4; ++i) {
        switch(i) {
            case 0: {
                out += (y0 * L[i]);
                break;
            }
            case 1: {
                out += (y1 * L[i]);
                break;
            }
            case 2: {
                out += (y2 * L[i]);
                break;
            }
            case 3: {
                out += (y3 * L[i]);
                break;
            }
        }
    }
    return out;
}

/* ------------------End Helper Functions -------------- */
void ext_main(void *r)
{
    // object initialization, note the use of dsp_free for the freemethod, which is required
    // unless you need to free allocated memory, in which case you should call dsp_free from
    // your custom free function.
    
    t_class *c = class_new("wtblackman~", (method)t_wtblackman_new, (method)t_wtblackman_free, (long)sizeof(t_wtblackman), 0L, A_GIMME, 0);
    
    class_addmethod(c, (method)t_wtblackman_float,		"float",	A_FLOAT, 0);
    class_addmethod(c, (method)t_wtblackman_dsp64,		"dsp64",	A_CANT, 0);
    class_addmethod(c, (method)t_wtblackman_assist,	"assist",	A_CANT, 0);
    
    class_dspinit(c);
    class_register(CLASS_BOX, c);
    t_wtblackman_class = c;
}


void *t_wtblackman_new(t_symbol *s, long argc, t_atom *argv)
{
    t_wtblackman *x = (t_wtblackman *)object_alloc(t_wtblackman_class);
    
    if (x) {
        if (!wavetable) {
            wavetable = calloc(SAMPLES, sizeof(double));
            if (!wavetable) {
                /* Handle Error? */
            }
            init_wavetable(wavetable, SAMPLES);
        }
        x->N = SAMPLES;
        dsp_setup((t_pxobject *)x, 2);	// MSP inlets: arg is # of inlets and is REQUIRED!
        // use 0 if you don't need inlets
        outlet_new(x, "signal"); 		// signal outlet (note "signal" rather than NULL)
        ++instances;    // Global static variable
    }
    return (x);
}

void t_wtblackman_free(t_wtblackman *x)
{
    if (!--instances) {
        free(wavetable);
        wavetable = NULL;
    }
    dsp_free;
}


void t_wtblackman_assist(t_wtblackman *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_INLET) { //inlet
        switch (a) {
            case 0: {
                sprintf(s, "Frequency in Hz.");
                break;
            }
            case 1: {
                sprintf(s, "Sync signal - positive edge resets phasor to beginning of cycle (+ phase offset)");
                break;
            }
            default:
                break;
        }
//        sprintf(s, "I am inlet %ld", a);
    }
    else {	// outlet
        sprintf(s, "Output");
    }
}


void t_wtblackman_float(t_wtblackman *x, double f)
{
//    x->f = f;
}


// registers a function for the signal chain in Max
void t_wtblackman_dsp64(t_wtblackman *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
//    post("my sample rate is: %f", samplerate);
    
    // instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
    // the arguments passed are:
    // 1: the dsp64 object passed-in by the calling function
    // 2: the symbol of the "dsp_add64" message we are sending
    // 3: a pointer to your object
    // 4: a pointer to your 64-bit perform method
    // 5: flags to alter how the signal chain handles your object -- just pass 0
    // 6: a generic pointer that you can use to pass any additional data to your perform method

    x->index = 0;
    x->R1 = samplerate;
    x->N = SAMPLES;
    object_method(dsp64, gensym("dsp_add64"), x, t_wtblackman_perform64, 0, NULL);
}

// this is the 64-bit perform method audio vectors
void t_wtblackman_perform64(t_wtblackman *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
    
    double a;       // Coefficient: how much faster target frequency is from orig frequency
    double f0;
    long n = sampleframes;
    double *out = *outs;
    double *frqs =ins[0];
    
    f0 = x->R1 / FIX_ZERO((double)SAMPLES);
    
    x->f = *frqs++;
    a = x->f / f0;

    while (n--) {
        x->index += a;
        if (x->index > x->N - 1)
            x->index = fmod(x->index, x->N);
        *out++ = interpolate4((((long)x->index) - 1),
                              wavetable[mod((((long)x->index) - 1), x->N)],
                              ((long)x->index),
                              wavetable[mod(((long)x->index), x->N)],
                              (((long)x->index) + 1),
                              wavetable[mod((((long)x->index) + 1), x->N)],
                              (((long)x->index) + 2),
                              wavetable[mod((((long)x->index) + 2), x->N)],
                              x->index);
        x->f_1 = x->f;
        x->f = *frqs++;
        
        if (x->f != x->f_1) {
            /* Specified frequency has changed */
            a = x->f / f0;
        }
    }
}

