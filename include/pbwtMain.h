#ifndef MAIN
#define MAIN

#define _POSIX_C_SOURCE 200809L
#include <string.h>
#include "pbwt.h"
#include "version.h"
#include "math.h"

static PBWT *playGround (PBWT *p);
static void prettyPlot (PBWT *p, FILE *fp, int K);
static void exportSiteInfo (PBWT *p, FILE *fp, int f1, int f2);
static void siteFrequencySpectrum (PBWT *p);
static void recordCommandLine (int argc, char *argv[]);
int pbwtMain (int argc, char *argv[]);

#endif