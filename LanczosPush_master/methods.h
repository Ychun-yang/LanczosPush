#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <fstream>
#include <future>
#include <filesystem>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <numeric>
#include <thread>
#include <sys/time.h> 
#include <time.h>
#include <cstring>
#include<iomanip>
//#include <sys/stat.h>
#include "graph.h"
#include "random.h"
#include <cstdlib>
#include <set>
#include <Eigen/Dense>
//#include <boost/math/special_functions.hpp>
//#include <boost/filesystem.hpp>

using Vector = std::vector<double>;
using Matrix = std::vector<vector<double>>;

double get_current_time_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

void generateQueryNode(int num_query_nodes, int n, ofstream& outf) {
    for(int i=0; i<num_query_nodes; i++) {
        int node = drand() % n;
        outf << node << endl;
    }
}

bool create_directory(const std::string& path) {
    std::string command = "mkdir -p " + path;
    int status = std::system(command.c_str());
    return status == 0;
}


//Algorithm

// Power Method

void groundtruth(graph& g, int s, int t, const double eps, string filename)
{
    double start = get_current_time_sec();

    int n = g.n, itr = 0;
    vector<double> pi(n, 0);
    vector<double> residual(n, 0);
    
    residual[s] = 1;
    residual[t] = -1;
    //double r_sum = 1.;
    int L=10000;
    double cur=0;
    while (itr < L)
    {
        vector<double> new_residual(n, 0);

        for (int i = 0; i < n; i++)
        {
            pi[i] += residual[i];
            
            double increment = residual[i] / g.weightdegree[i];
            for (int j = 0; j < g.degree[i]; j++)
            {
                int ne = g.AdjList[i][j];
                new_residual[ne] += increment;
            }
        }

        residual.swap(new_residual);
        if (itr%10 == 0)
        {
            double cur_new = pi[s]/g.weightdegree[s]-pi[t]/g.weightdegree[t];
            if (abs(cur_new-cur)<eps){break;}
            cur = cur_new;
        }

        //r_sum *= 1. - alpha;
        itr ++;
    }
    
    double end = get_current_time_sec();

    //double time = end-begin;
    //cout << "l1-error: " << r_sum << endl;
    //cout << "iter_num: " << itr << endl;

    double gt = pi[s]/g.weightdegree[s]-pi[t]/g.weightdegree[t];

    string directory = "result/gt/" + filename;
    create_directory(directory);
    ofstream fout("result/gt/" + filename +"/"+to_string(s)+".txt");
    fout<<gt<<endl;
    fout<<start-end<<endl;
    fout<<itr<<endl;
    //fout<<s<<" "<<t<<" "<<gt<<" "<<time<<endl;
    fout.close();
}


// Lanczos Iteration
void lanczos(graph& g, int s, int t, int k, string filename)
{
    int n = g.n;
    //int k = 50;
    Matrix Q = Matrix(k+2, Vector(n,0.0));
    //Matrix T = Matrix(k, Vector(k, 0.0));
    vector<double> beta(k+2,0.0);
    vector<double> alpha(k+1,0.0);
    Q[1][s]=1./sqrt(g.weightdegree[s]) / sqrt( 1./g.weightdegree[s] + 1./g.weightdegree[t] );
    Q[1][t]=-1./sqrt(g.weightdegree[t]) / sqrt( 1./g.weightdegree[s] + 1./g.weightdegree[t] );

    double start = get_current_time_sec();
    for (int itr=1; itr<=k; itr++)
    {
        //compute q_i+1 = A * q_i - beta_i * q_i-1
        
        for (int idx=0; idx<n; idx++)
        {
            
            for (int adj_idx=0; adj_idx<g.degree[idx]; adj_idx++)
            {
                int ne = g.AdjList[idx][adj_idx];
                Q[itr+1][ne]+=Q[itr][idx] / sqrt( (double) g.weightdegree[idx]*g.weightdegree[ne]);
            }

            Q[itr+1][idx]-=beta[itr]*Q[itr-1][idx];

        }

        //compute alpha_i = <q_i+1,q_i>
        alpha[itr] = inner_product(Q[itr+1].begin(), Q[itr+1].end(), Q[itr].begin(), 0.0);

        //compute q_i+1 = q_i+1 - alpha_i * q_i
        for (int idx=0; idx<n; idx++)
        {
            Q[itr+1][idx] = Q[itr+1][idx] - alpha[itr] * Q[itr][idx];
        }

        //compute beta_i+1 = ||q_i+1||
        beta[itr+1] = sqrt(inner_product(Q[itr+1].begin(), Q[itr+1].end(), Q[itr+1].begin(), 0.0));

        //compute q_i+1 = q_i+1 / beta_i+1
        for (int idx=0; idx<n; idx++)
        {
            Q[itr+1][idx] = Q[itr+1][idx] / beta[itr+1];
        }
    }

    //compute f(T_k)e_1 using power method

    //int max_itr = (int) 1/Alpha* 10*log (10);
    int max_itr=k*k;
    vector<double> fT(k,0);
    vector<double> residual(k, 0);
    residual[0]=1;
    for (int itr=0; itr<max_itr; itr++)
    {
        vector<double> new_residual(k, 0);

        for (int i = 0; i < k; i++)
        {
            fT[i] += residual[i];
            
            //compute r_new = (1-Alpha) * T_k * r_cur
            if (i == 0)
            {
                new_residual[0] += alpha[1]* residual[0];
                new_residual[0] += beta[2]* residual[1];
            }
            else if (i == k-1)
            {
                new_residual[k-1] += alpha[k]* residual[k-1];
                new_residual[k-1] += beta[k]* residual[k-1];
            }
            else if (i >=1 && i<=k-2)
            {
                new_residual[i] += alpha[i+1]* residual[i];
                new_residual[i] += beta[i+1]* residual[i-1];
                new_residual[i] += beta[i+2]* residual[i+1];
            }

        }
        residual.swap(new_residual);
    }

    //compute pi = Q * f(T) * e_1
    
    double end = get_current_time_sec();

    double value = fT[0] * ( 1./(double)g.weightdegree[s] + 1./(double)g.weightdegree[t] );
    string directory = "result/lanczos/" + filename +"/"+to_string(s);
    create_directory(directory);
    ofstream fout("result/lanczos/" + filename +"/"+to_string(s) + "/" +to_string(k)+".txt");
    fout<<value<<endl;
    fout<<end-start<<endl;
    fout.close();

    //ifstream gt_f("result/" + filename + "/gt/" + to_string(s) + to_string(t) + ".txt");
    //double gt_value;
    //gt_f >> gt_value;
    //gt_f.close();
    //double error = abs(gt_value-value);

}


// Lanczos Push
void lanczos_Local(graph& g, const int s, const int t, const double eps, int k, string filename)
{
    int n = g.n;
    //const vector<vector<int>> &G = graph.a;
    //int k = 50;
    Matrix Q = Matrix(k+2, Vector(n,0.0));
    vector<vector<int>> candidate_set(k+2, vector<int>());
    vector<vector<int>> supp(k+2, vector<int>());

    vector<int> candidate_count(k+2,0);
    vector<int> supp_count(k+2,0);
    vector<vector<bool>> Is_Insupp(k+2,vector<bool>(n, false));

    vector<double> beta(k+2,0.0);
    vector<double> alpha(k+1,0.0);
    Q[1][s]=1./sqrt(g.weightdegree[s]) / sqrt( 1./g.weightdegree[s] + 1./g.weightdegree[t] );
    Q[1][t]=-1./sqrt(g.weightdegree[t]) / sqrt( 1./g.weightdegree[s] + 1./g.weightdegree[t] );
    candidate_set[1].push_back(s);
    candidate_set[1].push_back(t);
    candidate_count[1]=2;
    supp[1].push_back(s);
    supp[1].push_back(t);
    supp_count[1]=2;
    double total_count=0;
    //bool flags = 0;

    double start = get_current_time_sec();

    for (int itr=1; itr<=k; itr++)
    {

        //if(flags==0){
        //compute q_i+1 = A * q_i - beta_i * q_i-1
        total_count+=supp_count[itr];
        for (int idx=0; idx<candidate_count[itr]; idx++)
        {
            int cur=candidate_set[itr][idx];
            if (Is_Insupp[itr+1][cur] == 0)
            {
                supp[itr+1].push_back(cur);
                supp_count[itr+1]++;
                Is_Insupp[itr+1][cur] = 1;
            }

            for (int adj_idx=0; adj_idx<g.degree[cur]; adj_idx++)
            {
                int ne = g.AdjList[cur][adj_idx];
                Q[itr+1][ne]+=Q[itr][cur] / sqrt( (double) g.weightdegree[cur]*g.weightdegree[ne]);
                if (Is_Insupp[itr+1][ne] == 0)
                {
                    supp[itr+1].push_back(ne);
                    supp_count[itr+1]++;
                    Is_Insupp[itr+1][ne] = 1;
                }
            }

        }

        for (int idx=0; idx<candidate_count[itr-1]; idx++)
        {
            int cur=candidate_set[itr-1][idx];
            Q[itr+1][cur]-=beta[itr]*Q[itr-1][cur];
            if (Is_Insupp[itr+1][cur] == 0)
            {
                supp[itr+1].push_back(cur);
                supp_count[itr+1]++;
                Is_Insupp[itr+1][cur] = 1;
            }
            
        }

        //compute alpha_i = <q_i+1,q_i>
        for (int idx=0; idx<supp_count[itr]; idx++)
        {
            int cur=supp[itr][idx];
            alpha[itr] += Q[itr+1][cur]*Q[itr][cur];
        }

        //compute q_i+1 = q_i+1 - alpha_i * q_i
        for (int idx=0; idx<candidate_count[itr]; idx++)
        {
            int cur=candidate_set[itr][idx];
            Q[itr+1][cur] = Q[itr+1][cur] - alpha[itr] * Q[itr][cur];
        }

        //compute beta_i+1 = ||q_i+1||
        for (int idx=0; idx<supp_count[itr+1]; idx++)
        {
            int cur=supp[itr+1][idx];
            beta[itr+1] += Q[itr+1][cur]*Q[itr+1][cur];
        }
        beta[itr+1]=sqrt(beta[itr+1]);

        //compute q_i+1 = q_i+1 / beta_i+1
        for (int idx=0; idx<supp_count[itr+1]; idx++)
        {
            int cur=supp[itr+1][idx];
            Q[itr+1][cur] = Q[itr+1][cur] / beta[itr+1];
            if (abs(Q[itr+1][cur] / sqrt((double) g.degree[cur])) > eps)
            {
                candidate_set[itr+1].push_back(cur);
                candidate_count[itr+1]++;
            }
        }
        //}

        //if(total_count>n/5){flags=1;}
    }
    //cout<<"flags"<<flags<<endl;
    //compute f(T_k)e_1 using power method

    //int max_itr = (int) 1/Alpha* 10*log (10);
    int max_itr=k*k;
    vector<double> fT(k,0);
    vector<double> residual(k, 0);
    residual[0]=1;
    for (int itr=0; itr<max_itr; itr++)
    {
        vector<double> new_residual(k, 0);

        for (int i = 0; i < k; i++)
        {
            fT[i] += residual[i];
            
            //compute r_new = (1-Alpha) * T_k * r_cur
            if (i == 0)
            {
                new_residual[0] += alpha[1]* residual[0];
                new_residual[0] += beta[2]* residual[1];
            }
            else if (i == k-1)
            {
                new_residual[k-1] += alpha[k]* residual[k-1];
                new_residual[k-1] += beta[k]* residual[k-1];
            }
            else if (i >=1 && i<=k-2)
            {
                new_residual[i] += alpha[i+1]* residual[i];
                new_residual[i] += beta[i+1]* residual[i-1];
                new_residual[i] += beta[i+2]* residual[i+1];
            }

        }
        residual.swap(new_residual);
    }

    double end = get_current_time_sec();


    //cout << "Lanczos local time: " << (double)(end - begin) / CLOCKS_PER_SEC << " s" << endl;

    //compute pi = Q * f(T) * e_1

    double value = fT[0] * ( 1./(double)g.weightdegree[s] + 1./(double)g.weightdegree[t] );
    string directory = "result/lanczos_local/" + filename+"/"+to_string(s);
    create_directory(directory);
    ofstream fout("result/lanczos_local/" + filename+"/"+ to_string(s)+"/" + to_string(eps)+".txt");
    fout<<value<<endl;
    fout<<end-start<<endl;
    fout.close();

}

void compare_results_prefix_plot(graph& g, int s, int t, string filename, string method, double time_result[][7], double error_result[][7],double deg_err_result[][7],double nnz[][7],double conductance[][7], int &method_pointer) {
    cout << method << endl;
    int n = g.n;
    double method_time;
    double* gt = new double[n];
    double* method_result = new double[n];
    for(int i=0; i<5; i++) {
        double eps = time_result[i][0];
        ifstream gt_f("result/" + filename + "/gt/" + to_string(s) + to_string(t) + ".txt");
        ifstream method_f("result/" + filename + "/" + method + "/" + to_string(eps)+"/"+ to_string(s) + to_string(t) +".txt");
        //cout<<"result/" + filename + "/" + method + "/" +to_string(alpha)+"/" + to_string(eps)+"/"+ to_string(s) + ".txt"<<endl;
        for(int j=0; j<n; j++) {
        gt_f >> gt[j];
        method_f >> method_result[j];
        }
        method_f >> method_time;
        //cout<<"method time: "<<method_time<<endl;
        gt_f.close();
        method_f.close();
        double l1_error = 0.;
        double deg_error = 0.;
        uint nnz_entries = 0 ;
        for(int j=0; j<n; j++) {
        //double err=0;
        //if(g.weightdegree[i]>0){err = abs(method_result[i] - gt[i])/g.weightdegree[i];}
        double err = abs(method_result[j] - gt[j])/g.weightdegree[j];
        double err_abs = abs(method_result[j] - gt[j]);
        //if(deg_error < err) deg_error = err;
        if(deg_error < err){deg_error = err;}
        //l1_error += err*g.weightdegree[i];
        l1_error+=err_abs;
        if(method_result[j]!=0){nnz_entries++;}
        }
        //double cur_conductance=sweepCut(g, s, method_result);
        //cout << "deg_normalize_error: " << deg_error << endl;
        //cout << "l1_error: " << l1_error << endl;
        //cout << "nnz entries: " << nnz << endl;
        //cout << "time: " << method_time << endl;
        time_result[i][method_pointer] += method_time;
        error_result[i][method_pointer] += l1_error;
        deg_err_result[i][method_pointer] += deg_error;
        nnz[i][method_pointer] += nnz_entries;
        //conductance[i][method_pointer] +=cur_conductance;
    }
    
    delete[] gt;
    delete[] method_result;

    method_pointer++;
}

void plot_results(graph& g, int num_query_nodes, string filename) {

    double time_result[6][2];
    double error_result[6][2];
    for(int i=0; i<6; i++) {
        for(int j=0; j<2; j++) {
            time_result[i][j] = 0;
            error_result[i][j] = 0;
        }
    }


    int method_pointer;

    ifstream infile("data/" + filename + ".query");
    for(int i=0; i<num_query_nodes; i++) {
        int s,t;
        infile >> s;
        infile >> t;
        vector<double> eps={0.1,0.05, 0.01, 0.005,0.001,0.0005,0.0001};
        vector<int> k={5,10,15,20,25,30};
        //double eps=0.005;int k=20;
        for(int i=0; i<6; i++) {
            //double eps=pow(10,-i-1);
            //int k=10*(i+1);
            //double eps2=0.5-0.1*i;
            ifstream gt_f("result/gt/" + filename +"/"+to_string(s)+".txt");
            ifstream method_f1("result/lanczos/" + filename+"/"+ to_string(s) +"/" + to_string(k[i])+".txt");
            ifstream method_f2("result/lanczos_local/" + filename+"/" + to_string(s)+"/" + to_string(eps[i])+".txt");

            double gt_value, v1,v2;
            double t1,t2;
            gt_f >> gt_value;
            gt_f.close();
            method_f1 >> v1;
            method_f1 >> t1;
            method_f1.close();
            error_result[i][0]+= abs(gt_value-v1);
            time_result[i][0]+= t1;

            method_f2 >> v2;
            method_f2 >> t2;
            method_f2.close();
            error_result[i][1]+= abs(gt_value-v2);
            time_result[i][1]+= t2;
        }
    }
    infile.close();
    for(int i=0; i<6; i++) {
        for(int j=1; j<2; j++) {
            time_result[i][j] /= (double)num_query_nodes;
            error_result[i][j] /= (double)num_query_nodes;
        }
    }

    string directory = "plot_result/" + filename;
    create_directory(directory);

    ofstream fout1("plot_result/" + filename + "/time.txt");
    for(int i=0; i<6; i++) {
        for(int j=0; j<2; j++) {
            fout1 << time_result[i][j] << " ";
        }
        fout1 << endl;
    }
    fout1.close();

    ofstream fout2("plot_result/" + filename + "/error.txt");
    for(int i=0; i<6; i++) {
        for(int j=0; j<2; j++) {
            fout2 << error_result[i][j] << " ";
        }
        fout2 << endl;
    }
    fout2.close();
}