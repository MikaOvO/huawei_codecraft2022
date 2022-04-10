#include <iostream>
#include <map>
#include <string>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>
#include <queue>
#include <assert.h>
#include <unistd.h>
#include <ctime>
#include <set>
#include <unordered_set>
#include <omp.h>
#include <unistd.h>
#include <signal.h>
#include <functional>

using namespace std;

typedef pair<int, int> P;
typedef pair<P, int> LP;
typedef pair<int, P> RP;
typedef pair<P, P> PP;

const int maxn = 1e5+10;
const int inf = 1e9;

ofstream *debug_file = nullptr;
ofstream *result_file = nullptr;

string output_dir = "/output";
string data_dir = "/data";
string info = "";

const int MAXT = 8928 + 1;
const int MAXM = 135 + 1;
const int MAXN = 135 + 1; 
const int MAXK = 100 + 1;
const int BASE = 10000;

int has_time_cost = 1;

int this_ab_index = 0;

// params
const int consumer_nd_block = 10;
const int producer_block = 1;
const int dff_producer_block = 1;
double Cik = 0.0;
int simT = 100;

clock_t begin_time;
void sig_handler(int num);
void CheckTime() {
    double time = (double)(clock() - begin_time) / CLOCKS_PER_SEC;
    if (time > 275.0) {
        sig_handler(0);
        exit(0);
    }
}

double Rand() {
    int range = 10000;
    return (double)(rand()%range)/range;
}


// AB test
int is_ab;
map<string, vector<int> > ab_candidate_map;
map<string, int> ab_map;

// utils
int StringToInt(string &str);
void Strip(string& str, const char& chr);
vector<string> Split(const string& str, const string& delim);

// IO
void ReadIn();
void Output();

// name map
map<string, int> producers_name_id_map;
map<string, int> consumers_name_id_map;
map<string, int> stream_name_id_map[MAXT];

// number
int times, producer_number, consumer_number;

int cost_max_index, can_full_use_time, cost_max_index90, can_full_use_time90;

int pt[MAXT];

int base_cost;

// func

void WorkTimeBaseline(int time, int time_index, int nd_write);
void PreWork(int nd_write, int left_weight, int right_weight);
void PreFirst(int nd_write, int left_weight, int right_weight);

void Reset(int reset_has_cost, int reset_has_use_full_time, int reset_write);
void Reset(int reset_has_cost, int reset_has_use_full_time, int reset_write, int time);
int EndWork(int sort_type, int loop_type);
int EndSwapWork(int sort_type, int loop_type);
int EndTryFullUse();
int EndWorkWithCost(int sort_type, int loop_type);
double EndWorkWithCostVec(int sort_type, int loop_type);
void EndTryDiscard();

// 暂时不考虑为0，也就是扔掉边缘的情况
double CalCost(int Ci, int Wi) {
    if (Wi <= base_cost) return base_cost;
    double tmp = (double) (Wi - base_cost) * (double) (Wi - base_cost) / (double) Ci  + (double) Wi;
    return tmp;
}

// 计算此次放置会多出多少耗费
// double CalCost(int Ci, int Wi, int cur, int ci) {
//     return CalCost(Ci, max(Wi, cur + ci)) - CalCost(Ci, Wi);
// }

// 考虑增长率
double CalCost(int Ci, int Wi, int cur, int ci) {
    double x1 = CalCost(Ci, max(Wi, cur + ci)) - CalCost(Ci, Wi);
    // double x2 = (double)(max(Wi, cur + ci) - base_cost) / (double)Ci;
    // x1 /= (double((max(Wi, cur + ci) - Wi)) + 0.00000001);
    x1 += Cik * Ci;
    return x1;
}

struct Producer {
    string name;
    int pagerank;
    double pr,nxt_pr;
    int bandwidth;
    int has_cost;
    int has_use_full_time;
    int can_visit_point[MAXM];
    vector<int> can_visit_point_vec;
    int need_time_cost_sum[MAXT];
    int is_full_use_time[MAXT]; // 是否为5%使用日
    int time_cost[MAXT];
    int ini_need_time_cost_sum[MAXT]; // 便于未来reset
    int lst_ab_has_cost;
    int cost90;
    int discard; // 是否扔掉这个边缘，暂时全不扔
    unordered_set<int> consumer_set[MAXT];
    Producer() {
        for (int i = 0; i < MAXT; ++i) {
            is_full_use_time[i] = 0;
            need_time_cost_sum[i] = 0;
            time_cost[i] = 0;
            consumer_set[i].clear();
        }
        cost90 = 0;
        can_visit_point_vec.clear();
        has_cost = 0;
        has_use_full_time = 0;
        discard = 0;
        lst_ab_has_cost = 0;
    }
    void Reset(int reset_has_cost, int reset_has_use_full_time) {
        if(reset_has_cost) has_cost = base_cost;
        if(reset_has_use_full_time) has_use_full_time = 0;
        for (int i = 1; i <= times; ++i) {
            if(reset_has_use_full_time) is_full_use_time[i] = 0;
            need_time_cost_sum[i] = ini_need_time_cost_sum[i];
            time_cost[i] = 0;
        }
    }
    void Reset(int reset_has_cost, int reset_has_use_full_time, int time) {
        if(reset_has_cost) has_cost = base_cost;
        if(reset_has_use_full_time) has_use_full_time = 0;
        if(reset_has_use_full_time) is_full_use_time[time] = 0;
        need_time_cost_sum[time] = ini_need_time_cost_sum[time];
        time_cost[time] = 0;
    }
    
    bool operator < (const Producer& producer) const {
        return bandwidth < producer.bandwidth;
    }
    int TimeRemain(int time) {
        return bandwidth - time_cost[time];
    }
    int HasCostRemain(int time) {
        // 得到不花钱的部分， has cost可能超过bd变成V
        return min(has_cost, bandwidth) - time_cost[time];
    }
    int GetScore() {
        int *p = new int[MAXT];
        for (int time = 1; time <= times; ++time) {
            p[time] = time_cost[time];
        }
        sort(p + 1, p + 1 + times); 
        // int ret = 0;
        // for (int time = cost_max_index; time >= 2; --time) {
        //     ret += (p[time] - p[time - 1]) * time * 5; 
        // }        
        // ?
        int ret = p[cost_max_index];
        delete[] p;
        return ret;
    }
    double GetAnswer(int update_has_cost, int update_full_use_time) {
        int *p = new int[MAXT];
        for (int time = 1; time <= times; ++time) {
            p[time] = time_cost[time];
        }
        sort(p + 1, p + 1 + times); 
        double ret = CalCost(bandwidth, p[cost_max_index]);
        if (cost90) ret = CalCost(bandwidth, p[cost_max_index90]);
        if (p[times] == 0) ret = 0.0;
        // 更新容量
        if (update_has_cost) {
            has_cost = max(base_cost, p[cost_max_index]);
            if (cost90) has_cost = max(base_cost, p[cost_max_index90]);
        }
        if (update_full_use_time) {
            for (int i = 1; i <= times; ++i) {
                is_full_use_time[i] = 0;
            }
            int num = 0;
            for (int i = 1; i <= times; ++i) {
                if (num == can_full_use_time) break;
                if (cost90 == 0 && time_cost[i] <= p[cost_max_index]) continue;
                if (cost90 == 1 && time_cost[i] <= p[cost_max_index90]) continue;
                ++num;
                is_full_use_time[i] = 1;
            }
            int *indexs = new int[MAXT];
            // shuffle ?
            for (int time = 1; time <= times; ++time) indexs[time] = time;
            for (int i = 1; i <= times; ++i) {
                int time = indexs[i];
                if (num == can_full_use_time) break;
                if (cost90 == 0 && time_cost[i] < p[cost_max_index]) continue;
                if (cost90 == 1 && time_cost[i] < p[cost_max_index90]) continue;
                ++num;
                is_full_use_time[time] = 1;
            }
            delete[] indexs;
        }
        delete[] p;
        return ret;
    }
} producers[MAXN];

struct TimeNode {
    int sum_cost;
    int ini_sum_cost;
    int sum_value;
    TimeNode() {
        sum_cost = 0;
        sum_value = 0;
    }
    void Reset() {
        sum_cost = ini_sum_cost;
        sum_value = ini_sum_cost;
    }
} time_node[MAXT];

struct Stream {
    string name;
    Stream() {
        
    }
} streams[MAXT][MAXK];

struct Consumer {
    string name;
    int pagerank;
    double pr,nxt_pr;
    vector<int> can_visit_point_vec;
    int can_visit_point[MAXN];
    int time_need[MAXT][MAXK];
    int ini_time_need[MAXT][MAXK];
    int degrees;
    Consumer() {
        can_visit_point_vec.clear();
        degrees = 0;
    }
    void Reset() {
        for (int i = 1; i <= times; ++i) {
            for (int j = 1; j <= pt[i]; ++j)
                time_need[i][j] = ini_time_need[i][j];
        }
    }
    void Reset(int i) {
        for (int j = 1; j <= pt[i]; ++j) {
            time_need[i][j] = ini_time_need[i][j];
        }
    }
    void NodifyProducer(int time, int cost_bandwidth) {
        for (auto& producer_id : can_visit_point_vec) {
            producers[producer_id].need_time_cost_sum[time] -= cost_bandwidth;
        }
    }
    void NodifyTimeNode(int time, int cost_bandwidth) {
        time_node[time].sum_cost -= cost_bandwidth;
        time_node[time].sum_value -= cost_bandwidth;
    }
} consumers[MAXM];

// time, consumer_id, stream_id -> producer_id
int info_bandwidth[MAXT][36][MAXK]; 
int best_answer = 2000000000;
int best_info_bandwidth[MAXT][36][MAXK];

int GetAnswer() {
    double ret = 0.0;
    for (int i = 1; i <= producer_number; ++i) ret += producers[i].GetAnswer(0, 0);
    ret += 0.5;
    return (int)ret;
}

vector<PP> time_consumer_vec[MAXT];
vector<PP> time_consumer_vec_degree[MAXT];
vector<PP> all_consumer_vec;

bool consumer_cmp(const PP& p1, const PP& p2) {
    if (p1.first.first / consumer_nd_block == p2.first.first / consumer_nd_block) {
        return p1.first.second < p2.first.second;
    }
    return p1.first.first < p2.first.first;
}


void Init() {
    //srand(298758457);
    srand((unsigned)('a'+'s'+'o'+'u'+'l'));
    //srand((unsigned)time(NULL));
    ::cost_max_index = (times * 95 + 99) / 100;
    ::cost_max_index90 = (times * 90 + 99) / 100;
    
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "cost_max_index: " << cost_max_index 
                      << endl;
    }
    ::can_full_use_time = times - cost_max_index;
    ::can_full_use_time90 = times - cost_max_index90;
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            if (producers[producer_id].can_visit_point[consumer_id]) {
                ++consumers[consumer_id].degrees;
                producers[producer_id].can_visit_point_vec.emplace_back(consumer_id);
                consumers[consumer_id].can_visit_point_vec.emplace_back(producer_id);
                for (int time = 1; time <= times; ++time) {
                    for (int stream_id = 1; stream_id <= pt[time]; ++stream_id)
                        producers[producer_id].need_time_cost_sum[time] += consumers[consumer_id].time_need[time][stream_id];
                }
            }
        }
    }

    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     producers[producer_id].pr = 1.0 / producer_number;
    // }
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     consumers[consumer_id].pr = 1.0 / consumer_number;
    // }

    // for (int k = 0; k < 100; ++k) {
    //     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //         producers[producer_id].nxt_pr = 0.0;
    //     }
    //     for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //         consumers[consumer_id].nxt_pr = 0.0;
    //     }    

    //     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //         for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //             if (producers[producer_id].can_visit_point[consumer_id] == 0) continue;
    //             producers[producer_id].nxt_pr += consumers[consumer_id].pr / consumers[consumer_id].can_visit_point_vec.size();
    //             consumers[consumer_id].nxt_pr += producers[producer_id].pr / producers[producer_id].can_visit_point_vec.size();
    //         }
    //     }

    //     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //         producers[producer_id].pr = producers[producer_id].nxt_pr;
    //     }
    //     for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //         consumers[consumer_id].pr = consumers[consumer_id].nxt_pr;
    //     }
    // }

    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     producers[producer_id].pagerank = producers[producer_id].pr * 10000000;
    //     if (debug_file != nullptr) {
    //         (*debug_file) << "producer_pr, "
    //                       << " producer_id: " <<producer_id
    //                       << " degree: " << producers[producer_id].can_visit_point_vec.size()
    //                       << " pr: " << producers[producer_id].pagerank
    //                       << endl; 
    //     }
    // }
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     consumers[consumer_id].pagerank = consumers[consumer_id].pr * 10000000;
    //     if (debug_file != nullptr) {
    //         (*debug_file) << "consumer_pr, "
    //                       << " consumer_id: " << consumer_id
    //                       << " degree: " << consumers[consumer_id].can_visit_point_vec.size()
    //                       << " pr: " << consumers[consumer_id].pagerank
    //                       << endl; 
    //     }
    // }

    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int time = 1; time <= times; ++time) {
            for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
                time_node[time].sum_cost += consumers[consumer_id].time_need[time][stream_id];
                time_node[time].sum_value += consumers[consumer_id].time_need[time][stream_id];
            }
        }
    }

    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int time = 1; time <= times; ++time) {
            for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
                // 写死每天的情况
                if (consumers[consumer_id].time_need[time][stream_id] > 0) {
                    // all_consumer_vec.emplace_back(PP(P(consumers[consumer_id].time_need[time][stream_id], 1), 
                    //                                  P(time, consumer_id * BASE + stream_id)));
                    // time_consumer_vec[time].emplace_back(PP(P(consumers[consumer_id].time_need[time][stream_id], -(int)consumers[consumer_id].can_visit_point_vec.size()), 
                    //                                         P(consumer_id, stream_id)));
                    time_consumer_vec[time].emplace_back(PP(P(consumers[consumer_id].time_need[time][stream_id], -(int)consumers[consumer_id].can_visit_point_vec.size()), 
                                                            P(consumer_id, stream_id)));
                    time_consumer_vec_degree[time].emplace_back(PP(P(-consumers[consumer_id].can_visit_point_vec.size(), consumers[consumer_id].time_need[time][stream_id]), 
                                                            P(consumer_id, stream_id)));
                }                    
                consumers[consumer_id].ini_time_need[time][stream_id] = consumers[consumer_id].time_need[time][stream_id];
            }
        }
    }

    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int time = 1; time <= times; ++time) {
            producers[producer_id].ini_need_time_cost_sum[time] = producers[producer_id].need_time_cost_sum[time];
        }
    }
    
    // sort(all_consumer_vec.begin(), all_consumer_vec.end());
    // reverse(all_consumer_vec.begin(), all_consumer_vec.end());

    for (int time = 1; time <= times; ++time) {
        sort(time_consumer_vec[time].begin(), time_consumer_vec[time].end(), consumer_cmp);
        reverse(time_consumer_vec[time].begin(), time_consumer_vec[time].end());
        sort(time_consumer_vec_degree[time].begin(), time_consumer_vec_degree[time].end());
        reverse(time_consumer_vec_degree[time].begin(), time_consumer_vec_degree[time].end());
        time_node[time].ini_sum_cost = time_node[time].sum_cost;
    }



    // ab_candidate_map["pre_producer_sort"] = {0, 1, 2};
    // ab_candidate_map["pre_consumer_sort"] = {0, 1, 2};
}

void MoveSomeBandWidth(int time, int from_producer_id, int to_producer_id, int consumer_id, int stream_id, int bandwidth, int update_use_cost, int nd_delete) {
    if (debug_file != nullptr && is_ab == 0) {
        if ( time <= 5 && consumer_id <= 3 && bandwidth > 0 ) {
            (*debug_file) << "MoveSomeBandWidth time: " << time  
                          << " from_producer_id: " << from_producer_id
                          << " to_producer_id: " << to_producer_id 
                          << " consumer_id: " << consumer_id
                          << " stream_id: " << stream_id
                          << " bandwidth: "<< bandwidth 
                          << endl;
        }
    }
    producers[to_producer_id].time_cost[time] += bandwidth;
    producers[from_producer_id].time_cost[time] -= bandwidth;

    if (update_use_cost && producers[to_producer_id].is_full_use_time[time] == 0) {
        producers[to_producer_id].has_cost = max(producers[to_producer_id].has_cost, producers[to_producer_id].time_cost[time]);
    }

    // not nodify

    info_bandwidth[time][consumer_id][stream_id] = to_producer_id;
    // 延迟删除
    if (nd_delete == 1) {
        producers[from_producer_id].consumer_set[time].erase(consumer_id * BASE + stream_id);
    }
    producers[to_producer_id].consumer_set[time].insert(consumer_id * BASE + stream_id);
}

void AddSomeBandWidth(int time, int producer_id, int consumer_id, int stream_id, int bandwidth, int update_use_cost=0, int nd_write=1) {
    if (bandwidth == 0) return ;
    if (debug_file != nullptr && is_ab == 0) {
        if (time <= 3 && producer_id <= 3 && consumer_id <= 3 && stream_id <= 3) {
            (*debug_file) << "AddSomeBandWidth time: " << time 
                          << " producer_id: " << producer_id
                          << " consumer_id: " << consumer_id 
                          << " stream_id: " << stream_id
                          << " bandwidth: "<< bandwidth 
                          << endl;
        }
    }
    // 把耗费放到当天里，注意策略里不需要写这一块了
    producers[producer_id].time_cost[time] += bandwidth;
    consumers[consumer_id].time_need[time][stream_id] -= bandwidth;
    consumers[consumer_id].NodifyProducer(time, bandwidth);
    consumers[consumer_id].NodifyTimeNode(time, bandwidth);

    // 更新一下理论花费（实际上比这要大，因为预先分配的不一定真的那最高的5%
    if (update_use_cost && producers[producer_id].is_full_use_time[time] == 0) {
            // 更新容量
            producers[producer_id].has_cost = max(producers[producer_id].has_cost, producers[producer_id].time_cost[time]);
    }

    if (nd_write == 0) return ;

    info_bandwidth[time][consumer_id][stream_id] = producer_id;
    producers[producer_id].consumer_set[time].insert(consumer_id * BASE + stream_id);
}

// 降序首次适应
bool DFF(int time, vector<P>& producer_vec, int update_use_cost=0, int nd_write=1) {
    int *lst = new int[MAXN];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) lst[producer_id] = -1;
    for (auto &pp : producer_vec) {
        lst[pp.second] = pp.first;
    }
    int consumer_id, producer_id, stream_id, bandwidth;
    int best_producer_id, best_value, best_degree;
    for (auto& cp : time_consumer_vec[time]) {
        consumer_id = cp.second.first;
        stream_id = cp.second.second;
        bandwidth = consumers[consumer_id].time_need[time][stream_id];
        if (bandwidth == 0) continue;
        best_producer_id = -1;
        best_value = -1;
        best_degree = 100000;
        for (int i = 0 ; i < consumers[consumer_id].can_visit_point_vec.size(); ++i) {
            int producer_id = consumers[consumer_id].can_visit_point_vec[i];
            // 找能放得下他最大的那一个
            if (lst[producer_id] >= bandwidth) {
                if (lst[producer_id] / dff_producer_block > best_value / dff_producer_block ||
                    lst[producer_id] / dff_producer_block == best_value / dff_producer_block && producers[producer_id].can_visit_point_vec.size() < best_degree) {
                    best_value = lst[producer_id];
                    best_producer_id = producer_id;
                    best_degree = producers[producer_id].can_visit_point_vec.size();
                }
            }
        }
        if (best_producer_id >= 1) {
            AddSomeBandWidth(time, best_producer_id, consumer_id, stream_id, bandwidth, update_use_cost, nd_write);    
            lst[best_producer_id] -= bandwidth;    
        }
    }
    delete[] lst;
    return true;
}

void PreWork(int nd_write, int left_weight, int right_weight) {
    int *ban = new int[MAXT];

    for (int i = 1; i <= times; ++i) ban[i] = 0;

    for (int k = 1; k <= can_full_use_time * producer_number; ++k) {
        int best_producer_id = -1, best_time = -1, can_min_last = -1, time_sum_use = -1, node_degree = 100000;
        for (int time = 1; time <= times; ++time) {
            if (ban[time] == 1) continue;
            if (time_sum_use < time_node[time].sum_value) {
                //(time_sum_use == time_node[time].sum_cost && can_last > can_min_last) ||  
                time_sum_use = time_node[time].sum_value;
                best_time = time;
            }
        }

        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].has_use_full_time == can_full_use_time) continue;
            if (producers[producer_id].is_full_use_time[best_time] == 1) continue;
            if (producers[producer_id].discard == 1) continue;
            int can_last = min(producers[producer_id].need_time_cost_sum[best_time], producers[producer_id].TimeRemain(best_time));
            // 按value来排(can_cost)
            // 不先考虑放不满一半的
            if (left_weight != -1 && can_last * left_weight < producers[producer_id].TimeRemain(best_time) * right_weight) continue;
            if (can_last / producer_block > can_min_last / producer_block ||
                can_last / producer_block == can_min_last / producer_block && producers[producer_id].can_visit_point_vec.size() < node_degree) {
                can_min_last = can_last;
                best_producer_id = producer_id;
                node_degree = producers[producer_id].can_visit_point_vec.size();
            }
        }

        if (time_sum_use <= 0) {
            break ;
        }
        // 仍有可以降低的天
        if (best_producer_id <= 0) {
            --k;
            ban[best_time] = 1;
            continue ;
        }

        // todo
        for (auto& cp : time_consumer_vec[best_time]) {
            int consumer_id = cp.second.first;
            int stream_id = cp.second.second;
            if (consumers[consumer_id].time_need[best_time][stream_id] == 0) continue;
            if (producers[best_producer_id].can_visit_point[consumer_id] == 0) continue;
            int bandwidth = producers[best_producer_id].TimeRemain(best_time);
            if (bandwidth >= consumers[consumer_id].time_need[best_time][stream_id]) {
                bandwidth = consumers[consumer_id].time_need[best_time][stream_id];
                AddSomeBandWidth(best_time, best_producer_id, consumer_id, stream_id, bandwidth, 0, 0);
            }
        }

        // update time_node sum value
        time_node[best_time].sum_value += producers[best_producer_id].lst_ab_has_cost;
        
        // update producer
        producers[best_producer_id].is_full_use_time[best_time] = 1;
        ++producers[best_producer_id].has_use_full_time;

        if (debug_file != nullptr && k <= 100 && is_ab == 0) {
            (*debug_file) << "prework full_use, producer_id: " << best_producer_id
                          << " time: " << best_time 
                          << " his number: " << producers[best_producer_id].has_use_full_time 
                          << " last: " << can_min_last
                          << " degree: " << node_degree
                          << " time_sum_cost: " << time_sum_use 
                          << endl;
        }
    }
    delete[] ban;
}

void PreFirst(int nd_write) {
    vector<pair<pair<int,int>,int> > can_cost_time;
    vector<LP > tmp;  
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        // tmp.emplace_back(LP(P(producers[producer_id].bandwidth / producer_block, -(int)producers[producer_id].can_visit_point_vec.size()), producer_id));
        tmp.emplace_back(LP(P(producers[producer_id].bandwidth / producer_block, 1), producer_id));
    }
    // 边缘从度数小的开始
    sort(tmp.begin(), tmp.end());
    reverse(tmp.begin(), tmp.end());
    int up;
    for (auto& p : tmp) {
        int producer_id = p.second;
        can_cost_time.clear();
        for (int time = 1; time <= times; ++time) {
            int can_cost = producers[producer_id].need_time_cost_sum[time];
            can_cost = min(can_cost, producers[producer_id].bandwidth);
            // -3 -3, -3 -2, -2 -100
            can_cost_time.emplace_back(make_pair(make_pair(-can_cost, -time_node[time].sum_value), time));
        }
        sort(can_cost_time.begin(), can_cost_time.end());
        // reverse(can_cost_time.begin(), can_cost_time.end());
        up = can_full_use_time;
        if (producers[producer_id].cost90) up = can_full_use_time90;
        for (int i = 0; i < up; ++i) {
            int time = can_cost_time[i].second;

            // update time_node sum value
            time_node[time].sum_value += producers[producer_id].lst_ab_has_cost;

            producers[producer_id].is_full_use_time[time] = 1;
            ++producers[producer_id].has_use_full_time;

            for (auto& cp : time_consumer_vec[time]) {
                int consumer_id = cp.second.first;
                int stream_id = cp.second.second;
                if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
                    continue;
                }
                int cost_bandwidth = producers[producer_id].TimeRemain(time);
                if (cost_bandwidth < consumers[consumer_id].time_need[time][stream_id]) continue;
                if (consumers[consumer_id].time_need[time][stream_id] == 0) continue;
                AddSomeBandWidth(time, producer_id, consumer_id, stream_id, consumers[consumer_id].time_need[time][stream_id], 0, nd_write);              
            }
        }
    }
}

void PreWork(int nd_write) {
    clock_t func_begin_time = clock();
    
    PreWork(nd_write, 3, 1);
    PreWork(nd_write, -1, -1);

    // Reset(0, 0);
    // #pragma omp parallel for num_threads(4)  
    // for(int time = 1; time <= times; ++time) {
    //     vector<P> producer_vec;
    //     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //         if (producers[producer_id].is_full_use_time[time]) {
    //             producer_vec.emplace_back(P(producers[producer_id].bandwidth, producer_id));
    //         }
    //     }
    //     DFF(time, producer_vec, 0, nd_write);
    // }

    if (debug_file != nullptr && is_ab == 0) {
        int *tmp_time_cost = new int[MAXT];
        int num = 0;
        (*debug_file) <<"prework5%:"<< endl;
        for (int i = 1; i <= min(producer_number, 10); ++i) {
            num = 0;
            for (int time = 1; time <= times; ++time) tmp_time_cost[time] = producers[i].time_cost[time];
            sort(tmp_time_cost + 1, tmp_time_cost + 1 + times);
            (*debug_file) << " producer_id: " << i << " bandwidth: " << producers[i].bandwidth <<endl;
            for (int time = cost_max_index ; time <= times; ++time) {
                (*debug_file) << " " << tmp_time_cost[time];
                if (num == 20) {
                    num = 0;
                    (*debug_file) << endl;
                }
                ++num;
            } 
            (*debug_file) << endl;
        }
        delete[] tmp_time_cost;
        (*debug_file) << endl;
    }
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "prework_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }
}

void WorkPre5(int time, int time_index, int nd_write) {
    vector<P> f_producer_vec;
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        if (producers[producer_id].is_full_use_time[time] == 1) {
            f_producer_vec.emplace_back(P(producers[producer_id].bandwidth, producer_id));    
        }
    }
    DFF(time, f_producer_vec, 0, nd_write);
}

void WorkPre5AndNoCost(int time, int time_index, int nd_write) {
    vector<P> f_producer_vec;
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        if (producers[producer_id].is_full_use_time[time] == 1) {
            f_producer_vec.emplace_back(P(producers[producer_id].bandwidth, producer_id));    
        } else {
            f_producer_vec.emplace_back(P(producers[producer_id].HasCostRemain(time), producer_id));        
        }
    }
    DFF(time, f_producer_vec, 0, nd_write);
}

void WorkNoCost(int time, int time_index, int nd_write) {
    vector<P> f_producer_vec;
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        if (producers[producer_id].is_full_use_time[time] == 0) {
            f_producer_vec.emplace_back(P(producers[producer_id].HasCostRemain(time), producer_id));        
        }
    }
    DFF(time, f_producer_vec, 0, nd_write);
}

int vis_time[MAXT];
void WorkTimeBaseline(int time, int time_index, int nd_write) {
    clock_t func_begin_time = clock();

    // todo
    vector<LP> producer_vec;
    for (int i = 1; i <= producer_number; ++i) {
        if (producers[i].discard == 1) continue;
        producer_vec.emplace_back(LP(P(producers[i].can_visit_point_vec.size(), producers[i].bandwidth), i));
    }
    sort(producer_vec.begin(), producer_vec.end());
    reverse(producer_vec.begin(), producer_vec.end());
    
    // vector<P> consumer_vec;
    // for (int i = 1; i <= consumer_number; ++i) {
    //     consumer_vec.emplace_back(P(consumers[i].can_visit_point_vec.size(), i));
    // }
    // sort(consumer_vec.begin(), consumer_vec.end());
    
    // 先放不花钱的
    // vector<P> f_producer_vec;
    // vector<RP> f_consumer_vec;
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     if (producers[producer_id].is_full_use_time[time] == 0) {
    //         f_producer_vec.emplace_back(P(producers[producer_id].has_cost, producer_id));
    //     }
    // }
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
    //         if (consumers[consumer_id].time_need[time][stream_id] > 0)
    //             f_consumer_vec.emplace_back(RP(consumers[consumer_id].time_need[time][stream_id], P(consumer_id, stream_id)));
    //     }
    // }
    // DFF(time, f_consumer_vec, f_producer_vec, 0, nd_write);
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
    //         int nd_bandwidth = consumers[consumer_id].time_need[time][stream_id];
    //         for(auto& pp : producer_vec) {
    //             int producer_id = pp.second;
    //             if (producers[producer_id].can_visit_point[consumer_id] == 0) continue;
    //             if (producers[producer_id].HasCostRemain(time) >= nd_bandwidth) {
    //                 AddSomeBandWidth(time, producer_id, consumer_id, stream_id, nd_bandwidth, 1, 1);
    //                 nd_bandwidth = 0;
    //                 break ;
    //             }
    //         }
    //     }
    // }
    for (auto& cp : time_consumer_vec[time]) {
        int consumer_id = cp.second.first;
        int stream_id = cp.second.second;
        int nd_bandwidth = consumers[consumer_id].time_need[time][stream_id];
        if (nd_bandwidth == 0) continue ;
        int best_producer_id = -1;
        double min_cost = 1e12;            
        for(auto& pp : producer_vec) {
            int producer_id = pp.second;
            if (producers[producer_id].discard == 1) continue;
            if (producers[producer_id].can_visit_point[consumer_id] == 0) continue;
            if (producers[producer_id].TimeRemain(time) < nd_bandwidth) continue;
            double tmp = (int)CalCost(producers[producer_id].bandwidth, producers[producer_id].has_cost, producers[producer_id].time_cost[time], nd_bandwidth);
            // tmp = -(int)producers[producer_id].bandwidth;
            // 同时刻先找的最大的
            if (tmp < min_cost) {
                min_cost = tmp;
                best_producer_id = producer_id;
            }
        }
        if (best_producer_id == -1) continue;
        AddSomeBandWidth(time, best_producer_id, consumer_id, stream_id, nd_bandwidth, 1, nd_write);
    }
    if (debug_file != nullptr && time_index <= 10 && is_ab == 0) {
        (*debug_file) << "workbaseline_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }

}

void WorkAllTime(int nd_write) {
    vector<LP> producer_vec;
    for (int i = 1; i <= producer_number; ++i) {
        if (producers[i].discard == 1) continue;
        producer_vec.emplace_back(LP(P(producers[i].bandwidth, 1), i));
    }
    sort(producer_vec.begin(), producer_vec.end());
    reverse(producer_vec.begin(), producer_vec.end());

    for (auto& cp : all_consumer_vec) {
        int time = cp.second.first;
        int consumer_id = cp.second.second / BASE;
        int stream_id = cp.second.second % BASE;
        int nd_bandwidth = consumers[consumer_id].time_need[time][stream_id];
        if (nd_bandwidth == 0) continue ;
        int best_producer_id = -1;
        double min_cost = 1e12;            
        for(auto& pp : producer_vec) {
            int producer_id = pp.second;
            if (producers[producer_id].discard == 1) continue;
            if (producers[producer_id].can_visit_point[consumer_id] == 0) continue;
            if (producers[producer_id].TimeRemain(time) < nd_bandwidth) continue;
            double tmp = (int)CalCost(producers[producer_id].bandwidth, producers[producer_id].has_cost, producers[producer_id].time_cost[time], nd_bandwidth);
            // tmp = -(int)producers[producer_id].bandwidth;
            // 同耗费先找的 最大的
            if (tmp < min_cost) {
                min_cost = tmp;
                best_producer_id = producer_id;
            }
        }
        if (best_producer_id == -1) continue;
        AddSomeBandWidth(time, best_producer_id, consumer_id, stream_id, nd_bandwidth, 1, nd_write);
    }     
}



void Work() {
    // discard
    // int d_number = 5, bd_number = 5;
    // vector<P> producer_vec;
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     producer_vec.emplace_back(P(producers[producer_id].bandwidth, producer_id));
    // }
    // sort(producer_vec.begin(), producer_vec.end());
    // for (auto& pp : producer_vec) {
    //     if (bd_number == 0) break;
    //     --bd_number;
    //     int producer_id = pp.second;
    //     int flag = 1;
    //     for (auto& consumer_id : producers[producer_id].can_visit_point_vec) {
    //         if (consumers[consumer_id].degrees <= 3) {
    //             flag = 0;
    //             break;
    //         }
    //     }
    //     if (flag) {
    //         producers[producer_id].discard = 1;
    //         for (auto& consumer_id : producers[producer_id].can_visit_point_vec) {
    //             --consumers[consumer_id].degrees;
    //         }
    //     }
    // }
    // producer_vec.clear();
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     producer_vec.emplace_back(P(producers[producer_id].can_visit_point_vec.size(), producer_id));
    // }
    // sort(producer_vec.begin(), producer_vec.end());
    // for (auto& pp : producer_vec) {
    //     if (d_number == 0) break;
    //     --d_number;
    //     int producer_id = pp.second;
    //     if (producers[producer_id].discard == 1) continue;
    //     int flag = 1;
    //     for (auto& consumer_id : producers[producer_id].can_visit_point_vec) {
    //         if (consumers[consumer_id].degrees <= 3) {
    //             flag = 0;
    //             break;
    //         }
    //     }
    //     if (flag) {
    //         producers[producer_id].discard = 1;
    //         for (auto& consumer_id : producers[producer_id].can_visit_point_vec) {
    //             --consumers[consumer_id].degrees;
    //         }
    //     }
    // }

    // PreWork(0);
    // if (result_file != nullptr) {
    //     long long all_has_cost = 0;
    //     int max_need = 0;
    //     for (int i = 1; i <= producer_number; ++i) {
    //         for (int time = 1; time <= times; ++time) {
    //             max_need = max(max_need, time_node[time].sum_cost);
    //             all_has_cost += producers[i].time_cost[time];
    //         }
    //     }
    //     (*result_file) << "after 5%, "
    //                     << " all_has_cost: " << all_has_cost
    //                     << " max_need: " << max_need
    //                     << endl;
    // }
    // todo
    // initToK();
    // 设置容量
    for (int i = 1; i <= producer_number; ++i) producers[i].has_cost = base_cost;

    vector<P> time_vec;
    for (int time = 1; time <= times; ++time) time_vec.emplace_back(P(time_node[time].sum_cost, time));
    // sort(time_vec.begin(), time_vec.end());
    // reverse(time_vec.begin(), time_vec.end());
    // random_shuffle(time_vec.begin(), time_vec.end());
    // for (int index = 0; index < times; ++index) {
    //     int time = time_vec[index].second;
    //     WorkTimeBaseline(time, index + 1, 0, 0);
    // }
    int upi = 30;
    for (int i = 1; i <= upi; ++i) {
        CheckTime();
        this_ab_index = i;
        if (i > 15) {
            random_shuffle(time_vec.begin(), time_vec.end());
            //  for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            //     int tmp = max(0, producers[producer_id].has_cost - base_cost);
            //     tmp = tmp * 5 / 100; // -10%
            //     producers[producer_id].has_cost = tmp; 
            //     producers[producer_id].has_cost = min(producers[producer_id].bandwidth, producers[producer_id].has_cost);
            // }
        }
        if (i == upi) is_ab = 0;
        // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        //     int tmp = max(0, producers[producer_id].has_cost - base_cost);
        //     tmp = tmp * 5 / 100; // -10%
        //     producers[producer_id].has_cost += tmp; 
        //     producers[producer_id].has_cost = min(producers[producer_id].bandwidth, producers[producer_id].has_cost);
        // }
        
        Reset(0, 1, 1);
        
        PreFirst(0);
        
        if (result_file != nullptr && is_ab == 1 && i == 1) {
            long long all_has_cost = 0;
            int max_need = 0;
            for (int i = 1; i <= producer_number; ++i) {
                for (int time = 1; time <= times; ++time) {
                    max_need = max(max_need, time_node[time].sum_cost);
                    all_has_cost += producers[i].time_cost[time];
                }
            }
            (*result_file) << "after 5%, "
                            << " all_has_cost: " << all_has_cost
                            << " max_need: " << max_need
                            << endl;
        }
        // Reset(0, 0, 0);
        // for (int index = 0; index < times; ++index) {
        //     int time = time_vec[index].second;
        //     WorkPre5(time, index + 1, 0);
        // }
        // WorkAllTime(0);
        Reset(0, 0, 0);

        for (int index = 0; index < times; ++index) {
            int time = time_vec[index].second;
            WorkPre5AndNoCost(time, index + 1, 1);
            WorkTimeBaseline(time, index + 1, 1);
        }
        // Reset(0, 0, 0);
        
        // for (int index = 0; index < times; ++index) {
        //     int time = time_vec[index].second;
        //     WorkPre5AndNoCost(time, index + 1, 1);
        //     WorkTimeBaseline(time, index + 1, 1);
        // }

        Cik = 0.0;

        int tmp = 0;
        EndWork(1, 0);
        EndWorkWithCostVec(1, 0);
        for (int j = 1; j <= 7; ++j) {
            if (EndWorkWithCostVec(0, 0) == 0) {
                ++tmp;
            } else {
                tmp = 0;
            }
            // if(j % 3 == 0) EndTryFullUse();
            if (tmp == 2) break;
        }
        EndSwapWork(1, 0);
        EndSwapWork(0, 0);
        
        if (debug_file != nullptr) {
            (*debug_file) << "AB_answer: " << GetAnswer()
                          << endl;
        }
        if (GetAnswer() < best_answer) {
            best_answer = GetAnswer();
            for (int time = 1; time <= times; ++time) {
                for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
                    for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
                        best_info_bandwidth[time][consumer_id][stream_id] = info_bandwidth[time][consumer_id][stream_id];
                    }
                }
            }
        }
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) { 
            producers[producer_id].lst_ab_has_cost = producers[producer_id].has_cost;
        }
        // if (i == 1) {
        //     vector<pair<long long, int> > discard_vec;
        //     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        //         long long sum = 0;
        //         for (int time = 1; time <= times; ++time) sum += producers[producer_id].time_cost[time];
        //         discard_vec.emplace_back(P(sum, producer_id));
        //     }
        //     sort(discard_vec.begin(), discard_vec.end());
        //     int discard_num = producer_number / 30;
        //     for (int k = 0; k < discard_num; ++k) {
        //         int producer_id = discard_vec[k].second;
        //         int flag = 1;
        //         for (auto& consumer_id : producers[producer_id].can_visit_point_vec) {
        //             if (consumers[consumer_id].degrees <= 3) {
        //                 flag = 0;
        //                 break;
        //             }
        //         }
        //         if (flag) {
        //             producers[producer_id].discard = 1;
        //             for (auto& consumer_id : producers[producer_id].can_visit_point_vec) {
        //                 --consumers[consumer_id].degrees;
        //             }
        //         }
        //     }
        // }
    }
    if (result_file != nullptr) {
        long long all_has_cost = 0;
        for (int i = 1; i <= producer_number; ++i) {
            for (int time = 1; time <= times; ++time) {
                if (producers[i].is_full_use_time[time])
                    all_has_cost += producers[i].time_cost[time];
            }
        }
        (*result_file) << "End 5%, "
                       << " all_has_cost: " << all_has_cost
                       << endl;
    }
}


unordered_set<int>::iterator it[MAXT];
unordered_set<int> tmp[MAXT];
int EndWorkWithCost(int sort_type, int loop_type) {
    clock_t func_begin_time = clock();
    
    // todo
    vector<P> producer_vec;
    int *pre_has_cost = new int[MAXN];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        producers[producer_id].GetAnswer(1, 1);
        producer_vec.emplace_back(1, producer_id);
        pre_has_cost[producer_id] = producers[producer_id].has_cost;
    }

    if(sort_type == 0) random_shuffle(producer_vec.begin(), producer_vec.end());
    if(sort_type == 1) {
        producer_vec.clear();
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].discard == 1) continue;
            producer_vec.emplace_back(producers[producer_id].GetScore(), producer_id);
        }
        sort(producer_vec.begin(), producer_vec.end());
        reverse(producer_vec.begin(), producer_vec.end());
    }

    int print_number = 0;
    int all_win = 0;
    double cost_all_win = 0.0;
    priority_queue<P> max_que;
    P max_p, cmax_p;
    int nxt;
    int from_producer_id, to_producer_id, consumer_id, stream_id, bandwidth, nd_bandwidth;
    double best_down = 1.0, best_to_producer_id; 
    double this_down ;

    for (int i = 0; i < producer_vec.size(); ++i) {
        from_producer_id = producer_vec[i].second;
        for (int time = 1; time <= times; ++time) {
            tmp[time].clear();
        }
        while (max_que.size()) max_que.pop();
        int Max = base_cost; // 放到这个就没必要了
        for (int time = 1; time <= times; ++time) {
            if (producers[from_producer_id].is_full_use_time[time] == 1) continue;
            max_que.push(P(producers[from_producer_id].time_cost[time], time));
            it[time] = producers[from_producer_id].consumer_set[time].begin();
            // if (debug_file != nullptr && i<= 3) {
            //     (*debug_file) << "time set,"
            //                   << " from_producer_id: " << from_producer_id 
            //                   << " time: " << time
            //                   << " it: " << (producers[from_producer_id].consumer_set[time].end() == it[time])
            //                   << " size: " << producers[from_producer_id].consumer_set[time].size()
            //                   << endl;
            // }
        }
        while (true) {
            max_p = max_que.top(); 
            max_que.pop();
            if (max_p.first <= base_cost) break;
            if (max_que.size()) nxt = max_que.top().first;
            else nxt = 0;
            int time = max_p.second;
            // if (debug_file != nullptr && i<= 3) { 
            //     (*debug_file) << "EndQueue:" 
            //                   << " max_p: " << max_p.first << " " << max_p.second
            //                   << " pre_has_cost: " << producers[from_producer_id].time_cost[time]
            //                   << " set_size: " << producers[from_producer_id].consumer_set[time].size()
            //                   << endl;
            // }
            while (it[time] != producers[from_producer_id].consumer_set[time].end()) {
                auto& cp = *it[time];
                consumer_id = cp / BASE;
                stream_id = cp % BASE;
                if (consumers[consumer_id].ini_time_need[time][stream_id] == 0) {
                    continue ;
                }
                int down;
                if (loop_type == 0) down = 0;
                if (loop_type == 1) down = i;
                best_to_producer_id = -1;
                best_down = 0.0;
                nd_bandwidth = consumers[consumer_id].ini_time_need[time][stream_id];
                for (int j = producer_vec.size() - 1; j >= down; --j) {
                    if (j == i) continue;
                    to_producer_id = producer_vec[j].second;
                    if (producers[to_producer_id].can_visit_point[consumer_id] == 0) continue;
                    if (producers[to_producer_id].is_full_use_time[time] == 1) {
                        bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]; 
                    } else {
                        bandwidth = producers[to_producer_id].has_cost - producers[to_producer_id].time_cost[time];     
                    }
                    if (bandwidth >= nd_bandwidth) {
                        best_to_producer_id = to_producer_id;
                        break ;
                    }
                    bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time];
                    if (bandwidth >= nd_bandwidth) {
                        // todo
                        // - xx
                        double tmp = CalCost(producers[from_producer_id].bandwidth, max(producers[from_producer_id].time_cost[time] - nd_bandwidth, nxt))
                                    -CalCost(producers[from_producer_id].bandwidth, producers[from_producer_id].time_cost[time]);
                        // + xx
                        tmp += CalCost(producers[to_producer_id].bandwidth, producers[to_producer_id].has_cost, producers[to_producer_id].time_cost[time], nd_bandwidth);
                        
                        if (tmp < best_down) {
                            best_down = tmp;
                            best_to_producer_id = to_producer_id;
                        }
                    } 
                }
                if (best_to_producer_id >= 1) {
                    if (debug_file != nullptr && i<= 3 && print_number <= 100) {
                        (*debug_file) << "EndWorkWithCost," 
                                      << " time: " << time
                                      << " from_producer_id: " << from_producer_id
                                      << " best_to_producer_id: " << best_to_producer_id
                                      << " consumer_id: " << consumer_id
                                      << " stream_id: " << stream_id
                                      << " nd_bandwidth: " << nd_bandwidth
                                      << " pre_cost: " << producers[from_producer_id].time_cost[time]
                                      << " best_down: " <<best_down
                                      << endl;
                        ++print_number;
                    }
                    MoveSomeBandWidth(time, from_producer_id, best_to_producer_id, consumer_id, stream_id, consumers[consumer_id].ini_time_need[time][stream_id], 1, 0);
                    tmp[time].insert(consumer_id * BASE + stream_id);
                    break ;
                } else {
                    // 枚举下一个需求
                    it[time]++; 
                }
            }
            // if (debug_file != nullptr && i<= 3) { 
            //     (*debug_file) << "find?"
            //                   << " " << (it[time] != producers[from_producer_id].consumer_set[time].end())
            //                   << endl;
            // }
            if (it[time] != producers[from_producer_id].consumer_set[time].end()) {
                max_que.push(P(producers[from_producer_id].time_cost[time], time));
                // 找到了 相当于下一次从下一个
                it[time]++;
            } else {
                break ;
            }
        }
        for (int time = 1; time <= times; ++time) {
            for (auto& consumer_info : tmp[time]) {
                producers[from_producer_id].consumer_set[time].erase(consumer_info);
            }
        }
        producers[from_producer_id].GetAnswer(1, 1);
        all_win += pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost;
        cost_all_win += CalCost(producers[from_producer_id].bandwidth, pre_has_cost[from_producer_id]) - 
                        CalCost(producers[from_producer_id].bandwidth, producers[from_producer_id].has_cost);
        if (debug_file != nullptr && is_ab == 0) {
            (*debug_file) << "EndWork, producer_id: " << from_producer_id
                          << " pre_has_cost: " << pre_has_cost[from_producer_id]
                          << " win: " << pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost
                          << endl;
        }
    }
    delete[] pre_has_cost;
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "endwork_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }
    if (debug_file != nullptr) {
        (*debug_file) << "endwork_all_win: " << all_win
                      << " endwork_cost_all_win: " << cost_all_win
                      << endl;
    }
    // 一直循环吧
    return 1;
}

int EndWork(int sort_type, int loop_type) {
    clock_t func_begin_time = clock();
    
    // todo
    vector<P> time_vec;
    vector<P> producer_vec;
    int *pre_has_cost = new int[MAXN];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        producers[producer_id].GetAnswer(1, 1);
        producer_vec.emplace_back(1, producer_id);
        pre_has_cost[producer_id] = producers[producer_id].has_cost;
    }

    if(sort_type == 0) random_shuffle(producer_vec.begin(), producer_vec.end());
    if(sort_type == 1) {
        producer_vec.clear();
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].discard == 1) continue;
            producer_vec.emplace_back(producers[producer_id].GetScore(), producer_id);
        }
        sort(producer_vec.begin(), producer_vec.end());
        reverse(producer_vec.begin(), producer_vec.end());
    }

    int all_win = 0;
    double cost_all_win = 0.0;
    unordered_set<int> tmp;
    int from_producer_id, to_producer_id, consumer_id, stream_id, bandwidth;
    for (int i = 0; i < producer_vec.size(); ++i) {
        from_producer_id = producer_vec[i].second;
        time_vec.clear();
        tmp.clear();
        int Max = base_cost; // 放到这个就没必要了
        for (int time = 1; time <= times; ++time) {
            if (producers[from_producer_id].is_full_use_time[time] == 1) continue;
            time_vec.emplace_back(P(producers[from_producer_id].time_cost[time], time));
        }
        sort(time_vec.begin(), time_vec.end());
        reverse(time_vec.begin(), time_vec.end());
        for (auto& tp : time_vec) {
            int time = tp.second;
            int cur = tp.first;
            if (tp.first <= Max) break;
            for (auto& cp : producers[from_producer_id].consumer_set[time]) {
                if (cur <= Max) break;
                consumer_id = cp / BASE;
                stream_id = cp % BASE;
                if (consumers[consumer_id].ini_time_need[time][stream_id] == 0) {
                    continue ;
                }
                int down;
                if (loop_type == 0) down = 0;
                if (loop_type == 1) down = i;
                for (int j = producer_vec.size() - 1; j >= down; --j) {
                    if (j == i) continue;
                    to_producer_id = producer_vec[j].second;
                    if (producers[to_producer_id].can_visit_point[consumer_id] == 0) continue;
                    if (producers[to_producer_id].is_full_use_time[time] == 1) {
                        bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]; 
                    } else {
                        bandwidth = producers[to_producer_id].has_cost - producers[to_producer_id].time_cost[time];     
                    }
                    if (bandwidth >= consumers[consumer_id].ini_time_need[time][stream_id]) {
                        MoveSomeBandWidth(time, from_producer_id, to_producer_id, consumer_id, stream_id, consumers[consumer_id].ini_time_need[time][stream_id], 0, 0);
                        cur -= consumers[consumer_id].ini_time_need[time][stream_id];
                        tmp.insert(consumer_id * BASE + stream_id);
                        break ;
                    }
                }
            }
            for (auto& consumer_info : tmp) {
                producers[from_producer_id].consumer_set[time].erase(consumer_info);
            }
            Max = max(Max, cur);
        }
        producers[from_producer_id].GetAnswer(1, 1);
        all_win += pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost;
        cost_all_win += CalCost(producers[from_producer_id].bandwidth, pre_has_cost[from_producer_id]) - 
                        CalCost(producers[from_producer_id].bandwidth, producers[from_producer_id].has_cost);
        if (debug_file != nullptr && is_ab == 0) {
            (*debug_file) << "EndWork, producer_id: " << from_producer_id
                          << " pre_has_cost: " << pre_has_cost[from_producer_id]
                          << " win: " << pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost
                          << endl;
        }
    }

    delete[] pre_has_cost;
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "endwork_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }
    if (debug_file != nullptr) {
        (*debug_file) << "endwork_all_win: " << all_win
                      << " endwork_cost_all_win: " << cost_all_win
                      << endl;
    }
    return all_win;
}

double EndWorkWithCostVec(int sort_type, int loop_type) {
    clock_t func_begin_time = clock();
    
    // todo
    vector<P> time_vec;
    vector<P> producer_vec;
    int *pre_has_cost = new int[MAXN];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        producers[producer_id].GetAnswer(1, 1);
        producer_vec.emplace_back(1, producer_id);
        pre_has_cost[producer_id] = producers[producer_id].has_cost;
    }

    if(sort_type == 0) random_shuffle(producer_vec.begin(), producer_vec.end());
    if(sort_type == 1) {
        producer_vec.clear();
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].discard == 1) continue;
            producer_vec.emplace_back(producers[producer_id].GetScore(), producer_id);
        }
        sort(producer_vec.begin(), producer_vec.end());
        reverse(producer_vec.begin(), producer_vec.end());
    }

    int all_win = 0;
    double cost_all_win = 0.0;
    unordered_set<int> tmp;
    int from_producer_id, to_producer_id, consumer_id, stream_id, bandwidth, nd_bandwidth;
    int best_to_producer_id;
    double best_down;
    for (int i = 0; i < producer_vec.size(); ++i) {
        from_producer_id = producer_vec[i].second;
        time_vec.clear();
        tmp.clear();
        int Max = base_cost; // 放到这个就没必要了
        for (int time = 1; time <= times; ++time) {
            if (producers[from_producer_id].is_full_use_time[time] == 1) continue;
            time_vec.emplace_back(P(producers[from_producer_id].time_cost[time], time));
        }
        sort(time_vec.begin(), time_vec.end());
        reverse(time_vec.begin(), time_vec.end());
        for (int tpi = 0; tpi < time_vec.size(); ++tpi) {
            auto& tp = time_vec[tpi];
            int nxt = base_cost;
            if (tpi < (int)time_vec.size() - 1) nxt = time_vec[tpi + 1].first;
            int time = tp.second;
            int cur = tp.first;
            if (tp.first <= Max) break;
            for (auto& cp : producers[from_producer_id].consumer_set[time]) {
                if (cur <= Max) break;
                consumer_id = cp / BASE;
                stream_id = cp % BASE;
                if (consumers[consumer_id].ini_time_need[time][stream_id] == 0) {
                    continue ;
                }
                int down;
                if (loop_type == 0) down = 0;
                if (loop_type == 1) down = i;
                nd_bandwidth = consumers[consumer_id].ini_time_need[time][stream_id];
                best_down = 0.0;
                best_to_producer_id = -1;
                for (int j = producer_vec.size() - 1; j >= down; --j) {
                    if (j == i) continue;
                    to_producer_id = producer_vec[j].second;
                    if (producers[to_producer_id].can_visit_point[consumer_id] == 0) continue;
                    if (producers[to_producer_id].is_full_use_time[time] == 1) {
                        bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]; 
                    } else {
                        bandwidth = producers[to_producer_id].has_cost - producers[to_producer_id].time_cost[time];     
                    }
                    if (bandwidth >= nd_bandwidth) {
                        best_to_producer_id = to_producer_id;
                        break ;
                    }
                    bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]; 
                    if (bandwidth >= nd_bandwidth) {
                        // todo
                        // - xx
                        double tmp = CalCost(producers[from_producer_id].bandwidth, max(producers[from_producer_id].time_cost[time] - nd_bandwidth, nxt))
                                    -CalCost(producers[from_producer_id].bandwidth, max(producers[from_producer_id].time_cost[time], nxt));
                        // + xx
                        tmp += CalCost(producers[to_producer_id].bandwidth, producers[to_producer_id].has_cost, producers[to_producer_id].time_cost[time], nd_bandwidth);
                        
                        if (tmp < best_down) {
                            best_down = tmp;
                            best_to_producer_id = to_producer_id;
                        }
                    } 
                }
                if (debug_file != nullptr && best_down < 0.0 && i <= 3) {
                    (*debug_file) << "down! best_down: " << best_down
                                  << " bandwidth: " << nd_bandwidth
                                  << endl;
                }
                if (best_to_producer_id >= 1) {
                    MoveSomeBandWidth(time, from_producer_id, best_to_producer_id, consumer_id, stream_id, consumers[consumer_id].ini_time_need[time][stream_id], 1, 0);
                    tmp.insert(consumer_id * BASE + stream_id);
                    cur -= consumers[consumer_id].ini_time_need[time][stream_id]; 
                }
            }
            for (auto& consumer_info : tmp) {
                producers[from_producer_id].consumer_set[time].erase(consumer_info);
            }
            Max = max(Max, cur);
        }
        producers[from_producer_id].GetAnswer(1, 1);
        all_win += pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost;
        cost_all_win += CalCost(producers[from_producer_id].bandwidth, pre_has_cost[from_producer_id]) - 
                        CalCost(producers[from_producer_id].bandwidth, producers[from_producer_id].has_cost);
        if (debug_file != nullptr && is_ab == 0) {
            (*debug_file) << "EndWork, producer_id: " << from_producer_id
                          << " pre_has_cost: " << pre_has_cost[from_producer_id]
                          << " win: " << pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost
                          << endl;
        }
    }

    delete[] pre_has_cost;
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "endwork_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }
    if (debug_file != nullptr) {
        (*debug_file) << "endwork_all_win: " << all_win
                      << " endwork_cost_all_win: " << cost_all_win
                      << endl;
    }
    return cost_all_win;
}

int EndSwapWork(int sort_type, int loop_type) {
    clock_t func_begin_time = clock();
    
    // todo
    vector<P> time_vec;
    vector<P> producer_vec;
    int *pre_has_cost = new int[MAXN];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        producers[producer_id].GetAnswer(1, 1);
        producer_vec.emplace_back(1, producer_id);
        pre_has_cost[producer_id] = producers[producer_id].has_cost;
    }

    if(sort_type == 0) random_shuffle(producer_vec.begin(), producer_vec.end());
    if(sort_type == 1) {
        producer_vec.clear();
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].discard == 1) continue;
            producer_vec.emplace_back(producers[producer_id].GetScore(), producer_id);
        }
        sort(producer_vec.begin(), producer_vec.end());
        reverse(producer_vec.begin(), producer_vec.end());
    }

    int all_win = 0;
    double cost_all_win = 0.0;
    unordered_set<int> tmp;
    int from_producer_id, to_producer_id, consumer_id, stream_id, bandwidth;
    int consumer_id2, stream_id2, pc;
    int find_flag;
    for (int i = 0; i < producer_vec.size(); ++i) {
        from_producer_id = producer_vec[i].second;
        time_vec.clear();
        tmp.clear();
        int Max = base_cost; // 放到这个就没必要了
        for (int time = 1; time <= times; ++time) {
            if (producers[from_producer_id].is_full_use_time[time] == 1) continue;
            time_vec.emplace_back(P(producers[from_producer_id].time_cost[time], time));
        }
        sort(time_vec.begin(), time_vec.end());
        reverse(time_vec.begin(), time_vec.end());
        for (auto& tp : time_vec) {
            int time = tp.second;
            int cur = tp.first;
            if (tp.first <= Max) break;
            for (auto& cp : producers[from_producer_id].consumer_set[time]) {
                if (cur <= Max) break;
                consumer_id = cp / BASE;
                stream_id = cp % BASE;
                find_flag = 0;
                if (consumers[consumer_id].ini_time_need[time][stream_id] == 0) {
                    continue ;
                }
                int down;
                if (loop_type == 0) down = 0;
                if (loop_type == 1) down = i;
                for (int j = producer_vec.size() - 1; j >= down; --j) {
                    if (j == i) continue;
                    to_producer_id = producer_vec[j].second;
                    if (producers[to_producer_id].can_visit_point[consumer_id] == 0) continue;
                    if (producers[to_producer_id].is_full_use_time[time] == 1) {
                        bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]; 
                    } else {
                        bandwidth = producers[to_producer_id].has_cost - producers[to_producer_id].time_cost[time];     
                    }
                    if (bandwidth >= consumers[consumer_id].ini_time_need[time][stream_id]) {
                        MoveSomeBandWidth(time, from_producer_id, to_producer_id, consumer_id, stream_id, consumers[consumer_id].ini_time_need[time][stream_id], 0, 0);
                        cur -= consumers[consumer_id].ini_time_need[time][stream_id];
                        tmp.insert(consumer_id * BASE + stream_id);
                        find_flag = 1;
                        break ;
                    }
                }
                if (find_flag == 1) continue;
                for (int j = producer_vec.size() - 1; j >= down; --j) {
                    if (j == i) continue;
                    to_producer_id = producer_vec[j].second;
                    if (producers[to_producer_id].can_visit_point[consumer_id] == 0) continue;
                    if (producers[to_producer_id].is_full_use_time[time] == 1) {
                        bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]; 
                    } else {
                        bandwidth = producers[to_producer_id].has_cost - producers[to_producer_id].time_cost[time];     
                    }
                    for (auto& to_cp : producers[to_producer_id].consumer_set[time]) {
                        consumer_id2 = to_cp / BASE;
                        stream_id2 = to_cp % BASE;
                        if (producers[from_producer_id].can_visit_point[consumer_id2] == 0) continue;
                        if (consumers[consumer_id2].ini_time_need[time][stream_id2] >= consumers[consumer_id].ini_time_need[time][stream_id]) continue;
                        pc = consumers[consumer_id].ini_time_need[time][stream_id] - consumers[consumer_id2].ini_time_need[time][stream_id2];
                        if (bandwidth >= pc) {
                            MoveSomeBandWidth(time, from_producer_id, to_producer_id, consumer_id, stream_id, consumers[consumer_id].ini_time_need[time][stream_id], 0, 0);
                            MoveSomeBandWidth(time, to_producer_id, from_producer_id, consumer_id2, stream_id2, consumers[consumer_id2].ini_time_need[time][stream_id2], 0, 0);
                            cur -= consumers[consumer_id].ini_time_need[time][stream_id];
                            cur += consumers[consumer_id2].ini_time_need[time][stream_id2];
                            tmp.insert(consumer_id * BASE + stream_id);
                            find_flag = 1;
                            break ;
                        }
                    }
                    if (find_flag) break;
                }
            }
            for (auto& consumer_info : tmp) {
                producers[from_producer_id].consumer_set[time].erase(consumer_info);
            }
            Max = max(Max, cur);
        }
        producers[from_producer_id].GetAnswer(1, 1);
        all_win += pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost;
        cost_all_win += CalCost(producers[from_producer_id].bandwidth, pre_has_cost[from_producer_id]) - 
                        CalCost(producers[from_producer_id].bandwidth, producers[from_producer_id].has_cost);
        if (debug_file != nullptr && is_ab == 0) {
            (*debug_file) << "EndSwapWork, producer_id: " << from_producer_id
                          << " pre_has_cost: " << pre_has_cost[from_producer_id]
                          << " win: " << pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost
                          << endl;
        }
    }

    delete[] pre_has_cost;
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "endswapwork_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }
    if (debug_file != nullptr) {
        (*debug_file) << "endswapwork_all_win: " << all_win
                      << " endswapwork_cost_all_win: " << cost_all_win
                      << endl;
    }
    return all_win;
}

int EndTryFullUse() {
    clock_t func_begin_time = clock();
    
    // todo
    vector<P> time_vec;
    vector<P> producer_vec;
    int *pre_has_cost = new int[MAXN];
    int *producer_vec_id = new int[MAXN];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].discard == 1) continue;
        producers[producer_id].GetAnswer(1, 1);
        producer_vec.emplace_back(1, producer_id);
        pre_has_cost[producer_id] = producers[producer_id].has_cost;
    }

    random_shuffle(producer_vec.begin(), producer_vec.end());

    for (int i = 0 ; i < producer_vec.size(); ++i) {
        producer_vec_id[producer_vec[i].second] = i;
    }

    int all_win = 0;
    unordered_set<int> tmp;
    int from_producer_id, to_producer_id, consumer_id, stream_id, bandwidth;
    for (int i = 0; i < producer_vec.size(); ++i) {
        to_producer_id = producer_vec[i].second;
        time_vec.clear();
        tmp.clear();
        for (int time = 1; time <= times; ++time) {
            // 5%理论已经没办法放了
            if (producers[to_producer_id].is_full_use_time[time] == 1) continue;
            time_vec.emplace_back(P(producers[to_producer_id].time_cost[time], time));
        }
        for (auto& tp : time_vec) {
            int time = tp.second;
            int time_id = -1, lb = 0, mid, ub = (int)time_consumer_vec[time].size() - 1;
            int lst = producers[to_producer_id].has_cost - producers[to_producer_id].time_cost[time];
            while (lb <= ub) {
                mid = (lb + ub) / 2;
                if (time_consumer_vec[time][mid].first.first <= lst) {
                    time_id = mid;
                    ub = mid - 1;
                } else {
                    lb = mid + 1;
                }
            }
            for (int j = time_id; j < time_consumer_vec[time].size(); ++j) {
                consumer_id = time_consumer_vec[time][j].second.first;
                stream_id = time_consumer_vec[time][j].second.second;
                from_producer_id = info_bandwidth[time][consumer_id][stream_id];
                bandwidth = time_consumer_vec[time][j].first.first;
                if (producers[to_producer_id].can_visit_point[consumer_id] == 0) continue;
                if (producer_vec_id[to_producer_id] <= producer_vec_id[from_producer_id] ) continue;
                if (producers[from_producer_id].is_full_use_time[time] == 1) continue;
                if (lst < bandwidth) continue;
                MoveSomeBandWidth(time, from_producer_id, to_producer_id, consumer_id, stream_id, bandwidth, 0, 1);
                lst -= bandwidth;
            }
        }
    }
    return 0;
}

// int nxt_hascost[MAXN];
// int tmp_up[MAXN][MAXT];
// void EndTryDiscard() {
//     clock_t func_begin_time = clock();
    
//     // todo
//     vector<P> time_vec;
//     vector<P> producer_vec;
//     int *pre_has_cost = new int[MAXN];

//     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
//         producers[producer_id].GetAnswer(1, 1);
//     }

//     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
//         if (producers[producer_id].discard == 1) continue;
//         int tmp = 0;
//         for (int time = 1; time <= times; ++time) tmp = max(tmp, producers[producer_id].time_cost[time]);
//         producer_vec.emplace_back(tmp, producer_id);
//     }
//     sort(producer_vec.begin(), producer_vec.end());
    
//     unordered_set<int> tmp;
//     int from_producer_id, to_producer_id, consumer_id, stream_id, bandwidth;
//     int nd_bandwidth;
//     for (int i = 0; i < producer_vec.size() / 2; ++i) {
//         from_producer_id = producer_vec[i].second;
//         time_vec.clear();
//         tmp.clear();
//         for (int time = 1; time <= times; ++time) {
//             // if (producers[from_producer_id].is_full_use_time[time] == 1) continue;
//             time_vec.emplace_back(P(producers[from_producer_id].time_cost[time], time));
//         }
//         sort(time_vec.begin(), time_vec.end());
//         reverse(time_vec.begin(), time_vec.end());
//         for (int j = 1; j <= producer_number; ++j) nxt_hascost[j] = producers[j].has_cost;
//         for (int j = 1; j <= producer_number; ++j) {
//             for (int time = 1; time <= times; ++time) {
//                 tmp_up[j][time] = 0;
//             }
//         }
//         for (auto& tp : time_vec) {
//             int time = tp.second;
//             int cur = tp.first;
//             for (auto& cp : producers[from_producer_id].consumer_set[time]) {
//                 consumer_id = cp / BASE;
//                 stream_id = cp % BASE;
//                 if (consumers[consumer_id].ini_time_need[time][stream_id] == 0) {
//                     continue ;
//                 }
//                 int down;
//                 for (int j = producer_vec.size() - 1; j >= 0; --j) {
//                     if (j == i) continue;
//                     to_producer_id = producer_vec[j].second;
//                     if (producers[to_producer_id].discard == 1) continue;
//                     if (producers[to_producer_id].can_visit_point[consumer_id] == 0) continue;
//                     bandwidth = producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]; 
                   
//                     if (bandwidth >= consumers[consumer_id].ini_time_need[time][stream_id]) {
//                         // MoveSomeBandWidth(time, from_producer_id, to_producer_id, consumer_id, stream_id, consumers[consumer_id].ini_time_need[time][stream_id], 0, 0);
                        
//                     }
//                 }
//             }
//         }
//         producers[from_producer_id].GetAnswer(1, 1);
//         all_win += pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost;
//         cost_all_win += CalCost(producers[from_producer_id].bandwidth, pre_has_cost[from_producer_id]) - 
//                         CalCost(producers[from_producer_id].bandwidth, producers[from_producer_id].has_cost);
//         if (debug_file != nullptr && is_ab == 0) {
//             (*debug_file) << "EndWork, producer_id: " << from_producer_id
//                           << " pre_has_cost: " << pre_has_cost[from_producer_id]
//                           << " win: " << pre_has_cost[from_producer_id] - producers[from_producer_id].has_cost
//                           << endl;
//         }
//     }

//     delete[] pre_has_cost;
//     if (debug_file != nullptr && is_ab == 0) {
//         (*debug_file) << "endwork_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
//                       << endl;
//     }
//     if (debug_file != nullptr) {
//         (*debug_file) << "endwork_all_win: " << all_win
//                       << " endwork_cost_all_win: " << cost_all_win
//                       << endl;
//     }
//     return all_win;   
// }

void Reset(int reset_has_cost, int reset_has_use_full_time, int reset_write) {
    for (int i = 1; i <= times; ++i) {
        time_node[i].Reset();
    }
    for (int i = 1; i <= producer_number; ++i) {
        producers[i].Reset(reset_has_cost, reset_has_use_full_time);
    }
    for (int i = 1; i <= consumer_number; ++i) {
        consumers[i].Reset();
    }
    if (reset_write) {
        for (int i = 1; i <= producer_number; ++i) {
            for (int time = 1; time <= times; ++time)
                producers[i].consumer_set[time].clear();
        }
    }
}

void Reset(int reset_has_cost, int reset_has_use_full_time, int reset_write, int time) {
    for (int i = 1; i <= times; ++i) {
        time_node[i].Reset();
    }
    for (int i = 1; i <= producer_number; ++i) {
        producers[i].Reset(reset_has_cost, reset_has_use_full_time, time);
    }
    for (int i = 1; i <= consumer_number; ++i) {
        consumers[i].Reset(time);
    }
    if (reset_write) {
        for (int i = 1; i <= producer_number; ++i) {
            producers[i].consumer_set[time].clear();
        }
    }
}

void sig_handler(int num) {
    Output();
    if (result_file != nullptr && info[0] != '!') {
        int sum_cost = 0;
        int check_flag = 1;
        for (int i = 1; i <= consumer_number; ++i) {
            for (int time = 1; time <= times; ++time) {
                for (int stream_id = 1; stream_id <= pt[time]; ++stream_id)
                    if (consumers[i].time_need[time][stream_id]) 
                    {
                        check_flag = 0;
                    }
            }
        }
        for (int i = 1; i <= producer_number; ++i) {
            for (int time = 1; time <= times; ++time) {
                if (producers[i].time_cost[time] < 0) 
                {
                    check_flag = 0;
                }
            }
        }
        if (check_flag == 0) {
            (*result_file) << "check fail!" << endl;
        }
        (*result_file) << "data_dir: " << data_dir <<endl;
        (*result_file) << "sum_cost: " << GetAnswer() 
                       << " best_answer: " << best_answer 
                       << " cost_time: " << (double)(clock() - begin_time) / CLOCKS_PER_SEC << "s" 
                       << endl;
        (*result_file) << endl;
    }    
    if (result_file != nullptr && debug_file != nullptr && is_ab == 0) {
        for (int i = 1; i <= producer_number; ++i) {
            sort(producers[i].time_cost + 1, producers[i].time_cost + 1 + times);
        }
        int num = 0;
        (*debug_file) << "********end********" << endl;
        for (int i = 1; i <= producer_number; ++i) {
            num = 0;
            (*debug_file) << "producer_id: " << i << " bandwidth: " << producers[i].bandwidth << endl;
            // after sort !
            (*debug_file) << "sample pre95%: " <<producers[i].time_cost[1];
            if (times <= 100) {
                for (int time = 2; time <= cost_max_index; ++time) {
                    (*debug_file) << " " << producers[i].time_cost[time];
                    if (num == 20) {
                        num = 0;
                        (*debug_file) << endl;
                    }
                    ++num;
                }
            } else {
                for (int time = 2; time <= cost_max_index; time += times * 2 / 100 ) {
                    (*debug_file) << " " << producers[i].time_cost[time];
                    if (num == 20) {
                        num = 0;
                        (*debug_file) << endl;
                    }
                    ++num;
                }    
            }
            num = 0;
            (*debug_file) << endl;
            (*debug_file) << "cost: " << producers[i].GetAnswer(0, 0) << endl;
            (*debug_file) << "last5%: ";
            for (int time = cost_max_index + 1; time <= times; ++time) {
                (*debug_file) << " " << producers[i].time_cost[time];
                if (num == 20) {
                    num = 0;
                    (*debug_file) << endl;
                }
                ++num;
            } 
            (*debug_file) << endl;
        }
    }
    if (debug_file != nullptr) debug_file->close();
    if (result_file != nullptr) result_file->close();
}

int main(int argc, char *argv[]) {
    // signal(SIGALRM, sig_handler);
    // alarm(2);
    begin_time = clock();
    // local or oj
    ::info = "xxx";
    for (int i = 1; i < argc; ++i) {
        int len = strlen(argv[i]);
        string tmp = "";
        for (int j = 0; j < len; ++j) tmp += argv[i][j];
        Strip(tmp, '\r');
        vector<string> str_vec = Split(tmp, "=");
        if (str_vec[0] == "info") {
            info = str_vec[1];
        } 
        if (str_vec[0] == "debug_file") {
            debug_file = new ofstream();
            debug_file->open(str_vec[1], ios::out);
        } 
        if (str_vec[0] == "result_file") {
            result_file = new ofstream();
            result_file->open(str_vec[1], ofstream::app);
        } 
        if (str_vec[0] == "output_dir") {
            output_dir = str_vec[1];
        } 
        if (str_vec[0] == "data_dir") {
            data_dir = str_vec[1];
        }
    }
    if (result_file != nullptr && info[0] != '!') {
        time_t now = time(0);
        (*result_file) << "info: " << endl;
        (*result_file) << "time: " << ctime(&now);
    }
    ReadIn();
    Init();
    is_ab = 1;
    // for (int k = 1; k <= 10; ++k) Work();
    // ABWork();
    // is_ab = 0;
    Work();
    // WorkBinarySearch();
    Output();

    if (result_file != nullptr && info[0] != '!') {
        int sum_cost = 0;
        int check_flag = 1;
        for (int i = 1; i <= consumer_number; ++i) {
            for (int time = 1; time <= times; ++time) {
                for (int stream_id = 1; stream_id <= pt[time]; ++stream_id)
                    if (consumers[i].time_need[time][stream_id]) 
                    {
                        check_flag = 0;
                    }
            }
        }
        if (check_flag == 0) {
            (*result_file) << "check fail! consumer != 0" << endl;
        }
        check_flag = 1;
        for (int i = 1; i <= producer_number; ++i) {
            for (int time = 1; time <= times; ++time) {
                if (producers[i].time_cost[time] < 0) 
                {
                    check_flag = 0;
                }
            }
        }
        if (check_flag == 0) {
            (*result_file) << "check fail! producer < 0" << endl;
        }
        check_flag = 1;
        for (int i = 1; i <= producer_number; ++i) {
            for (int time = 1; time <= times; ++time) {
                if (producers[i].time_cost[time] > producers[i].bandwidth) 
                {
                    check_flag = 0;
                }
            }
        }
        if (check_flag == 0) {
            (*result_file) << "check fail! producer > bd" << endl;
        }
        check_flag = 1;
        for (int i = 1; i <= producer_number; ++i) {
            if (producers[i].discard == 0) continue;
            for (int time = 1; time <= times; ++time) {
                if (producers[i].time_cost[time] > 0) 
                {
                    check_flag = 0;
                }
            }
        }
        if (check_flag == 0) {
            (*result_file) << "check fail! discard" << endl;
        }
        check_flag = 1;
        for (int time = 1; time <= times; ++time) {
            for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
                for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
                    if (consumers[consumer_id].ini_time_need[time][stream_id] == 0) continue;
                    if (consumers[consumer_id].can_visit_point[info_bandwidth[time][consumer_id][stream_id]] == 0)
                        check_flag = 0;
                }
            }
        }
        if (check_flag == 0) {
            (*result_file) << "check fail! edge" << endl;
        }
        (*result_file) << "data_dir: " << data_dir <<endl;
        (*result_file) << "sum_cost: " << GetAnswer() 
                       << " best_answer: " << best_answer 
                       << " cost_time: " << (double)(clock() - begin_time) / CLOCKS_PER_SEC << "s" 
                       << endl;
        (*result_file) << endl;
    }    
    if (result_file != nullptr && debug_file != nullptr && is_ab == 0) {
        for (int i = 1; i <= producer_number; ++i) {
            sort(producers[i].time_cost + 1, producers[i].time_cost + 1 + times);
        }
        int num = 0;
        (*debug_file) << "********end********" << endl;
        for (int i = 1; i <= producer_number; ++i) {
            num = 0;
            (*debug_file) << "producer_id: " << i << " bandwidth: " << producers[i].bandwidth << endl;
            // after sort !
            (*debug_file) << "sample pre95%: " <<producers[i].time_cost[1];
            if (times <= 100) {
                for (int time = 2; time <= cost_max_index; ++time) {
                    (*debug_file) << " " << producers[i].time_cost[time];
                    if (num == 20) {
                        num = 0;
                        (*debug_file) << endl;
                    }
                    ++num;
                }
            } else {
                for (int time = 2; time <= cost_max_index; time += times * 2 / 100 ) {
                    (*debug_file) << " " << producers[i].time_cost[time];
                    if (num == 20) {
                        num = 0;
                        (*debug_file) << endl;
                    }
                    ++num;
                }    
            }
            num = 0;
            (*debug_file) << endl;
            (*debug_file) << "cost: " << producers[i].GetAnswer(0, 0) << endl;
            (*debug_file) << "last5%: ";
            for (int time = cost_max_index + 1; time <= times; ++time) {
                (*debug_file) << " " << producers[i].time_cost[time];
                if (num == 20) {
                    num = 0;
                    (*debug_file) << endl;
                }
                ++num;
            } 
            (*debug_file) << endl;
        }
    }
    if (debug_file != nullptr) debug_file->close();
    if (result_file != nullptr) result_file->close();
    return 0;
}

void ReadIn() {
    // read consumers
    ifstream read_file;
    read_file.open(data_dir + "/demand.csv", ios::in);
    string data;
    int lines = 0;
    vector<string> vec_data;
    int cur_time = 0, cur_stream = 0;
    string lst_time = "????";
    while (read_file >> data) {
        ++lines;
        Strip(data, '\r');
        vec_data = Split(data, ",");    
        if (lines == 1) {
            for (int i = 2; i < vec_data.size(); ++i) {
                consumers_name_id_map[vec_data[i]] = i - 1;
                consumers[i - 1].name = vec_data[i];
            }
            consumer_number = vec_data.size() - 2;
            continue;
        }
        if (vec_data[0] != lst_time) {
            pt[cur_time] = cur_stream;
            ++cur_time;
            cur_stream = 0;
            lst_time = vec_data[0];
        }
        ++cur_stream;
        stream_name_id_map[cur_time][vec_data[1]] = cur_stream;
        streams[cur_time][cur_stream].name = vec_data[1];
        for (int i = 2; i < vec_data.size(); ++i) {
            consumers[i - 1].time_need[cur_time][cur_stream] = StringToInt(vec_data[i]);
            if (debug_file != nullptr && i <= 5 && lines - 1 <= 5) {
                (*debug_file) << "consumers: " << i - 1
                              << " time: " << cur_time 
                              << " stream: " << cur_stream
                              << " name: " << consumers[i - 1].name 
                              << " bandwidth: " << consumers[i - 1].time_need[cur_time][cur_stream] 
                              << endl;
            }
        }
    }
    times = cur_time;
    pt[cur_time] = cur_stream;
    
    read_file.close();
    lines = 0;
    
    // read producers
    read_file.open(data_dir + "/site_bandwidth.csv", ios::in);
    while (read_file >> data) {
        ++lines;
        Strip(data, '\r');
        vec_data = Split(data, ",");  
        if (lines == 1) continue;
        producers_name_id_map[vec_data[0]] = lines - 1;
        producers[lines - 1].name = vec_data[0];
        producers[lines - 1].bandwidth = StringToInt(vec_data[1]);
        if (debug_file != nullptr && lines - 1 <= 3) {
            (*debug_file) << "producers: " << lines - 1 << " bandwidth: " << producers[lines - 1].bandwidth 
            << " name: " << producers[lines - 1].name << endl;
        }
    }
    producer_number = lines - 1;

    // sort(producers + 1, producers + 1 + producer_number);
    // for (int i = 1; i <= producer_number; ++i) producers_name_id_map[producers[i].name] = i;

    read_file.close();   
    lines = 0;

    // read config

    int qos_down = -1;

    read_file.open(data_dir + "/config.ini", ios::in);
    while (read_file >> data) {
        ++lines;
        Strip(data, '\r');
        vec_data = Split(data, "=");        
        if (lines == 1) continue;
        if (lines == 2) {
            qos_down = StringToInt(vec_data[1]);
            if (debug_file != nullptr) {
                (*debug_file) << "qos_down: " << qos_down << endl;
            }
        }
        if (lines == 3) {
            base_cost = StringToInt(vec_data[1]);
            if (debug_file != nullptr) {
                (*debug_file) << "base_cost: " << base_cost << endl;
            }
        }
    }
    read_file.close();   
    lines = 0;

    // read qos map
    read_file.open(data_dir + "/qos.csv", ios::in);
    vector<int> consumer_ids;
    consumer_ids.emplace_back(-1); // index: 0 
    while (read_file >> data) {
        ++lines;
        Strip(data, '\r');
        vec_data = Split(data, ",");        
        if (lines == 1) {
            for (int i = 1; i < vec_data.size(); ++i) {
                consumer_ids.emplace_back(consumers_name_id_map[vec_data[i]]);
            }
            continue;
        }
        int producer_id = producers_name_id_map[vec_data[0]];
        for (int i = 1; i < vec_data.size(); ++i) {
            if (StringToInt(vec_data[i]) < qos_down) {
                int consumer_id = consumer_ids[i];
                producers[producer_id].can_visit_point[consumer_id] = 1;
                consumers[consumer_id].can_visit_point[producer_id] = 1;
                if (debug_file != nullptr && producer_id <= 3 && consumer_id <= 3) {
                    (*debug_file) << "producers " << producer_id << " consumers: "
                         << consumer_id << " has edge." <<endl;
                }
            } else {
                int consumer_id = consumer_ids[i];
                producers[producer_id].can_visit_point[consumer_id] = 0;
                consumers[consumer_id].can_visit_point[producer_id] = 0;                
            }
        }
    }
    read_file.close();   
    lines = 0;   

    if (debug_file != nullptr) {
        (*debug_file) << "read cost " << (double)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << endl;
    }
}

void Output() {
    if (output_dir == "NULL") return ; // local and neednot
    ofstream file;
    file.open(output_dir + "/solution.txt", ios::out);
    vector<int> consumer_get_info[MAXM][MAXM]; // consumer_id, producer_id -> stream_id
    for (int time = 1; time <= times; ++time) {
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
                if (best_info_bandwidth[time][consumer_id][stream_id] <= 0) {
                    continue ;
                }
                int producer_id = best_info_bandwidth[time][consumer_id][stream_id];
                consumer_get_info[consumer_id][producer_id].emplace_back(
                    stream_id
                );
            }
        }
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {   
            int begin_flag = 1;
            file<<consumers[consumer_id].name<<":";
            for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
                if (consumer_get_info[consumer_id][producer_id].size() == 0) continue;
                if (begin_flag == 1) {
                    begin_flag = 0;
                } else {
                    file<<",";
                }
                if (consumer_get_info[consumer_id][producer_id].size() > 0) {
                    file<<"<"<<producers[producer_id].name;
                }
                for(auto &p : consumer_get_info[consumer_id][producer_id]) {
                    file<<',';
                    file<<streams[time][p].name;
                }
                if (consumer_get_info[consumer_id][producer_id].size() > 0) {
                    file<<">";
                }
            }
            if (time != times || consumer_id != consumer_number) file<<'\n';
        }
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            for (int producer_id = 1; producer_id <= producer_number; ++producer_id) 
                consumer_get_info[consumer_id][producer_id].clear();
        }
    }
    file.close();
}

int StringToInt(string &str) {
    int ret = 0;
    for (int i = 0; i < str.length(); ++i) {
        ret = ret * 10 + str[i] - '0';
    }
    return ret;
}

void Strip(string& str, const char& chr) {
    while (str.length() && str[str.length() - 1] == chr) {
        str.pop_back();
    }
}

vector<string> Split(const string& str, const string& delim) {
    vector<string> ret;
    string tmp = "";
    for (int i = 0; i < str.length(); ++i) {
        int is_split = 1;
        for (int j = 0; j < delim.length(); ++j) {
            if (str[i] != delim[j]) {
                is_split = 0;
                break;
            }
        }
        if (is_split) {
            ret.emplace_back(tmp);
            tmp = "";
            i += delim.length() - 1;
        } else {
            tmp += str[i];
        }
    }
    if (tmp != "") {
        ret.emplace_back(tmp);
    }
    return ret;    
}