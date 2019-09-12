#ifndef _TEXTURE_H
#define _TEXTURE_H
#include "math.h"
#include <iostream>
#include <iomanip>
using namespace std;

float mean(float *p, int N);
float moment(float *p, int N, int k);
float contrast(float *p, int N);
float entropy(float *p, int N);
float energy(float *p, int N);
float moment_k(float *p, int N, int k);
float dissimilarity_cooc(float *p, int N);
float contraste_cooc(float *p, int N);

#endif