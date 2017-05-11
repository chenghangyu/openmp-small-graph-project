
#include "main.h"

using namespace std;

void sssp_openmp(int root, graph *g, int thd_num)
{	
    index_t v = g->vert_count;
    index_t e = g->edge_count;

    sa[root] = 0;
    for(index_t i=1; i<v; ++i)
    {
        sa[i] = INF;
    }
    
    
    int level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=1; i<v; ++i)
        {	
        	#pragma omp parallel for num_threads(thd_num)
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(sa[g->csr[j]] < INF)
                {
                    if(sa[i] > sa[g->csr[j]] + g->weight[j])
                    {
                        sa[i] = sa[g->csr[j]] + g->weight[j];
                        if(!flag)
                            flag = true;
                    }
                }
            }
        }
        
    }

}

void sssp_cpu(int root, graph *g)
{	
    index_t v = g->vert_count;
    index_t e = g->edge_count;
   
    sa[root] = 0;
    for(index_t i=1; i<v; ++i)
    {
        sa[i] = INF;
    }
    
    int level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=1; i<v; ++i)
        {
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(sa[g->csr[j]] < INF)
                {
                    if(sa[i] > sa[g->csr[j]] + g->weight[j])
                    {
                        sa[i] = sa[g->csr[j]] + g->weight[j];
                        if(!flag)
                            flag = true;
                    }
                }
            }
        }
    }
    
}

void sssp_cpu_for_apsp(int root, graph *g)
{	
    index_t v = g->vert_count;
    index_t e = g->edge_count;


    for(index_t i=0; i<v; ++i)
    {
        sa[i] = INF;
    }
    sa[root] = 0;
    
    int level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=0; i<v; ++i)
        {
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(sa[g->csr[j]] < INF)
                {
                    if(sa[i] > sa[g->csr[j]] + g->weight[j])
                    {
                        sa[i] = sa[g->csr[j]] + g->weight[j];
                        flag = true;
                    }
                }
            }
        }
    }
    for(index_t i=0; i<v; ++i)
        dist[root*V+i] = sa[i];
}

void sssp_openmp_for_apsp(int root, graph *g)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;

	path_t sa[V];
    for(index_t i=0; i<v; ++i)
    {
        sa[i] = INF;
    }
    sa[root] = 0;
    
    int level = 0;
    bool flag = true;
    
    while(flag)
    {
        flag = false;    
        for(index_t i=0; i<v; ++i)
        {
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(sa[g->csr[j]] < INF)
                {
                    if(sa[i] > sa[g->csr[j]] + g->weight[j])
                    {
                        sa[i] = sa[g->csr[j]] + g->weight[j];
                        flag = true;
                    }
                }
            }
        }
    }
    
    for(index_t i=0; i<v; ++i)
        dist[root*V+i] = sa[i];
}

void cc(int root, graph *g)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;

    for(index_t i=0; i<v; ++i)
    {
        sa[i] = INF;
    }
    sa[root] = 0;
    
    int level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=0; i<v; ++i)
        {
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(sa[g->csr[j]] < INF)
                {
                    if(sa[i] > sa[g->csr[j]] + g->weight[j])
                    {
                        sa[i] = sa[g->csr[j]] + g->weight[j];
                        flag = true;
                    }
                }
            }
        }
    }
    path_t sum = 0;
    for(index_t i=0; i<v; ++i)
    {
        //if a vertex is unreachable, the distance is 0
       
        if(sa[i] < INF && sa[i] > 0)
            sum += 1.0/sa[i];
    }
    
    global_cc[root] = sum;
}

void cc_openmp(int root, graph *g)
{
    index_t v = g->vert_count;
    index_t e = g->edge_count;
	path_t sa[V];

    for(index_t i=0; i<v; ++i)
    {
        sa[i] = INF;
    }
    sa[root] = 0;
    
    int level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=0; i<v; ++i)
        {
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(sa[g->csr[j]] < INF)
                {
                    if(sa[i] > sa[g->csr[j]] + g->weight[j])
                    {
                        sa[i] = sa[g->csr[j]] + g->weight[j];
                        flag = true;
                    }
                }
            }
        }
    }
    path_t sum = 0;
    for(index_t i=0; i<v; ++i)
    {
        //if a vertex is unreachable, the distance is 0
        if(sa[i] < INF && sa[i] > 0)
            sum += 1.0/sa[i];
    }
    //printf("%g\n", sum);
    global_cc[root] = sum;
}

index_t sssp_for_bc_cpu(index_t root, graph *g)
{
    memset(sp_count, 0, sizeof(index_t) * V);
    sp_count[root] = 1;
    for(index_t i=0; i<V; ++i)
    {
        dist_bc[i] = INF;
        sa[i] = INF;
    }
    dist_bc[root] = 0;
    sa[root] = 0;

    index_t level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=0; i<V; ++i)
        {
            bool flag_one = false;
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(dist_bc[g->csr[j]] < INF)
                {
                    if(dist_bc[i] > dist_bc[g->csr[j]] + g->weight[j])
                    {
                        dist_bc[i] = dist_bc[g->csr[j]] + g->weight[j];
                        sp_count[i] = 0;
                        sa[i] = sa[g->csr[j]] + 1;
                        if(sa[i] > level)
                            level = sa[i];
                        if(!flag_one)
                            flag_one = true;
                    }

                }
            }
            if(flag_one)
            {
                flag = true;
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                    if(dist_bc[i] == dist_bc[g->csr[j]] + g->weight[j])
                    {
                        sp_count[i] += sp_count[g->csr[j]];
                    }
            }

        }
    }
    return level;
}

void bc_one(index_t root, graph *g, index_t level)
{ 
    path_t bc_tmp[V];
    memset(bc_tmp, 0, sizeof(path_t)*V);
    for(index_t cur=level; cur>=0; --cur)
    {
        for(index_t i=0; i<V; ++i)
        {
            if(sa[i] == cur)
            {
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                {
                    index_t w = g->csr[j];
                    if(dist_bc[w] == dist_bc[i] + g->weight[j])
                    {
                        //if(sp_count[w] != 0)
                        bc_tmp[i] += sp_count[i]*1.0*(1+bc_tmp[w])/sp_count[w];
                    }

                }

            }
        }
    }
    bc_tmp[root] = 0;
    for(index_t i=0; i<V; ++i)
    {
        bc[i] += bc_tmp[i];
    }
}

void bc_all(graph *g)
{
    memset(bc, 0, sizeof(path_t)*V);
    for(index_t i=0; i<g->vert_count; ++i)
    {
        index_t level = sssp_for_bc_cpu(i, g);
        bc_one(i, g, level);
    }
}

void bc_one_openmp(index_t root, graph *g)
{   

path_t dist[V];
index_t sa[V];
index_t sp_count[V];


//sssp begin
    memset(sp_count, 0, sizeof(index_t) * V);
    sp_count[root] = 1;
    for(index_t i=0; i<V; ++i)
    {
        dist[i] = INF;
        sa[i] = INF;
    }
    dist[root] = 0;
    sa[root] = 0;

    index_t level = 0;
    bool flag = true;
    while(flag)
    {
        flag = false;    
        for(index_t i=0; i<V; ++i)
        {
            bool flag_one = false;
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(dist[g->csr[j]] < INF)
                {
                    if(dist[i] > dist[g->csr[j]] + g->weight[j])
                    {
                        dist[i] = dist[g->csr[j]] + g->weight[j];
                        sp_count[i] = 0;
                        sa[i] = sa[g->csr[j]] + 1;
                        if(sa[i] > level)
                            level = sa[i];
                        if(!flag_one)
                            flag_one = true;
                    }

                }
            }
            if(flag_one)
            {
                flag = true;
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                    if(dist[i] == dist[g->csr[j]] + g->weight[j])
                    {
                        sp_count[i] += sp_count[g->csr[j]];
                    }
            }

        }
    }//sssp end    
   
    
    path_t bc_tmp[V];
    memset(bc_tmp, 0, sizeof(path_t)*V);
    for(index_t cur=level; cur>=0; --cur)
    {
        for(index_t i=0; i<V; ++i)
        {
            if(sa[i] == cur)
            {
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                {
                    index_t w = g->csr[j];
                    if(dist[w] == dist[i] + g->weight[j])
                    {
                        bc_tmp[i] += sp_count[i]*1.0*(1+bc_tmp[w])/sp_count[w];
                    }

                }

            }
        }
    }
    
    
    
    bc_tmp[root] = 0;
    for(index_t i=0; i<V; ++i)
    {
    #pragma omp atomic
        bc[i] += bc_tmp[i];
    }
    
}

void bc_all_openmp(graph *g,int thd_num)
{
    memset(bc, 0, sizeof(path_t)*V);
    
    
    #pragma omp parallel for num_threads(thd_num)
    for(index_t i=0; i<g->vert_count; ++i)
    {
        bc_one_openmp(i, g);
    }

}

//main
int main(int args, char ** argv)
{
    printf("Input: ./smallgraph_openmp_test ./dataset/begin.bin ./dataset/adjacent.bin ./dataset/weight.bin tread_num\n");

    if(args != 5)
        exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
    const char *weight_filename = argv[3];    
    const int thd_num = atoi(argv[4]);
    graph *g = new graph(beg_filename, csr_filename, weight_filename);
    
    
//sssp cpu    
   
    for(int k=0;k<30;k++){
	    double btm1_start=wtime();
	    sssp_cpu(0, g);
	    btm1_sum += wtime()-btm1_start;
	    }
	    
    std::cout<<"\nSSSP   cpu   linear  success(30times), AVG runtime is: "<<btm1_sum/30.0<<" second(s)\n";
    print_sssp_cpu();

//apsp cpu
    memset(dist, 0, V*V*sizeof(path_t));
    
    for(int k=0;k<30;k++){
    	    double btm2_start=wtime();
	    for(int source=0; source<V; ++source)
	    {
		sssp_cpu_for_apsp(source, g);// g->vert_count, g->edge_count);
	    }
	    btm2_sum=wtime()-btm2_start;
	    }
    std::cout<<"APSP   cpu   linear  success(30times), AVG runtime is: "<<btm2_sum/30.0<<" seconds(s)\n";
    print_apsp_cpu();
 
//cc cpu
memset(global_cc, 0, V*sizeof(path_t));

    
    for(int k=0; k<30; k++)
    {
        double btm3_start = wtime();
        for(int source=0; source<V; ++source)
        {
            cc(source, g);// g->vert_count, g->edge_count);
        }
        btm3_sum += wtime() - btm3_start;
    }
    
    std::cout<<"  CC   cpu   linear  success(30times), AVG runtime is: "<<btm3_sum/30.0<<" seconds(s)\n";
    print_cc_cpu();
//bc cpu

    for(int k=0; k<30; k++)
    {
    double btm4_start=wtime();
    bc_all(g);
    btm4_sum += wtime() - btm4_start;
    }
    std::cout<<"  BC   cpu   linear  success(30times), AVG runtime is: "<<btm4_sum/30.0<<" seconds(s)\n";
    print_bc_cpu();
    
    std::cout<<"\nCPU OpenMp threads: "<<thd_num<<"\n";
//sssp openmp    

    for(int k=0; k<30; k++)
    {
        double btm5_start = wtime();
    sssp_openmp(0, g,thd_num);// g->vert_count, g->edge_count);
    btm5_sum += wtime() - btm5_start;
    }
       std::cout<<"SSSP OpenMp parallel success(30times), AVG runtime is: "<<btm5_sum/30.0<<" seconds(s)\n";
    print_sssp_openmp();
    
//apsp openmp    
    memset(dist, 0, V*V*sizeof(path_t));
    for(int k=0;k<30;k++)
    {
    double btm6_start=wtime();
    #pragma omp parallel for num_threads(thd_num)
    for(int source=0; source<V; ++source)
    {
        sssp_openmp_for_apsp(source, g);// g->vert_count, g->edge_count);
    }
      btm6_sum += wtime() - btm6_start;
    }
    
    std::cout<<"APSP OpenMp parallel success(30times), AVG runtime is: "<<btm6_sum/30.0<<" seconds(s)\n";
    print_apsp_openmp();
    
//cc openmp
    memset(global_cc, 0, V*sizeof(path_t));
  
    
    for(int k=0;k<30;k++)
    {
    double btm7_start = wtime();
    #pragma omp parallel for num_threads(thd_num)
    for(int source=0; source<V; ++source)
    {
        cc_openmp(source, g);// g->vert_count, g->edge_count);
    }
    btm7_sum += wtime()-btm7_start;
    }
     
    std::cout<<"  CC OpenMp parallel success(30times), AVG runtime is: "<<btm7_sum/30.0<<" seconds(s)\n";
    print_cc_openmp();
    
//bc openmp
    
    for(int k=0;k<30;k++)
    {  
    double btm8_start = wtime();
    bc_all_openmp(g,thd_num);  
     btm8_sum +=wtime()-btm8_start;
    }
    print_bc_openmp();
    
//output
    print_runtime_and_seedup();
    
    std::cout<<"  BC OpenMp parallel success(30times), AVG runtime is: "<<btm8_sum/30.0<<" seconds(s)\n"; 
    std::cout<<"\nSSSP speed up is: "<<btm1_sum/btm5_sum<<"\n";
    std::cout<<"APSP speed up is: "<<btm2_sum/btm6_sum<<"\n";
    std::cout<<"  CC speed up is: "<<btm3_sum/btm7_sum<<"\n";
    std::cout<<"  BC speed up is: "<<btm4_sum/btm8_sum<<"\n";
    return 0;
}
