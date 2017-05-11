#ifndef __H_MAIN__
#define __H_MAIN__

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

const char sssp_openmp_file[] = "./result/sssp_openmp.result";
const char sssp_cpu_file[] = "./result/sssp_cpu.result";
const char apsp_cpu_file[] = "./result/apsp_cpu.result";
const char apsp_openmp_file[]  = "./result/apsp_openmp.result";
const char cc_cpu_file[] = "./result/cc_cpu.result";
const char cc_openmp_file[] = "./result/cc_openmp.result";
const char bc_cpu_file[] = "./result/bc_cpu.result";
const char bc_openmp_file[] = "./result/bc_openmp.result";
const char runtime_and_seedup_file[] = "./result/runtime_and_seedup.result";

const int INF = 0x7fffffff;
const int V = 218;

double btm1_sum = 0;
double btm2_sum = 0;
double btm3_sum = 0;
double btm4_sum = 0;
double btm5_sum = 0;
double btm6_sum = 0;
double btm7_sum = 0;
double btm8_sum = 0;

path_t sa[V];
path_t dist[V*V];//apsp
path_t global_cc[V];//cc

path_t dist_bc[V];//
index_t sp_count[V];
path_t bc[V];

inline void print_sssp_openmp()
{
    FILE * fp = fopen(sssp_openmp_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%g\n", sa[i]);
    }
    fclose(fp);
}

inline void print_sssp_cpu()
{
FILE * fp = fopen(sssp_cpu_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%g\n", sa[i]);
    }
    fclose(fp);
}

inline void print_apsp_cpu()
{
    FILE * fp = fopen(apsp_cpu_file, "w");
    for(int i=0; i<V; ++i)
    {
        for(int j=0; j<V; ++j)
            fprintf(fp, "%g ", dist[i*V+j]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

inline void print_apsp_openmp()
{
    FILE * fp = fopen(apsp_openmp_file, "w");
    for(int i=0; i<V; ++i)
    {
        for(int j=0; j<V; ++j)
            fprintf(fp, "%g ", dist[i*V+j]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

inline void print_cc_cpu()
{
    FILE * fp = fopen(cc_cpu_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%g\n", global_cc[i]);
    }
    fclose(fp);
}
inline void print_cc_openmp()
{
    FILE * fp = fopen(cc_openmp_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%g\n", global_cc[i]);
    }
    fclose(fp);
}

inline void print_bc_cpu()
{
    FILE * fp = fopen(bc_cpu_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%d %g\n", i, bc[i]);
        //fprintf(fp, "%g\n", dist[i]);
    }
    fclose(fp);
}

inline void print_bc_openmp()
{
    FILE * fp = fopen(bc_openmp_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%d %g\n", i, bc[i]);
        //fprintf(fp, "%g\n", dist[i]);
    }
    fclose(fp);
}

inline void print_runtime_and_seedup()
{
    FILE * fp = fopen(runtime_and_seedup_file, "w");
    	fprintf(fp, "%s\n"," avg runting time (30times):");
        fprintf(fp, "%s %g\n", "sssp_cpu:", btm1_sum/30.0);
	fprintf(fp, "%s %g\n", "apsp_cpu:", btm2_sum/30.0);
	fprintf(fp, "%s %g\n", "  cc_cpu:", btm3_sum/30.0);
	fprintf(fp, "%s %g\n", "  bc_cpu:", btm4_sum/30.0);
	fprintf(fp, "%s %g\n", "sssp_openmp:", btm5_sum/30.0);
	fprintf(fp, "%s %g\n", "apsp_openmp:", btm6_sum/30.0);
	fprintf(fp, "%s %g\n", "  cc_openmp:", btm7_sum/30.0);
	fprintf(fp, "%s %g\n", "  bc_openmp:", btm8_sum/30.0);
	fprintf(fp, "%s\n","avg speedup:");
	fprintf(fp, "%s %g\n", "sssp_speedup:", btm1_sum/btm5_sum);
	fprintf(fp, "%s %g\n", "apsp_speedup:", btm2_sum/btm6_sum);
	fprintf(fp, "%s %g\n", "cc_speedup:", btm3_sum/btm7_sum);
	fprintf(fp, "%s %g\n", "bc_speedup:", btm4_sum/btm8_sum);
	
    
    fclose(fp);
}
#endif
