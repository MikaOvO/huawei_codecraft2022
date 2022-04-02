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

using namespace std;

typedef pair<int, int> P;
typedef pair<P, int> LP;
typedef pair<int, P> RP;
typedef pair<P, P> PP;

ofstream *debug_file = nullptr;
ofstream *result_file = nullptr;

string output_dir = "/output";
string data_dir = "/data";
string info = "";

const int MAXT = 8928 + 1;
const int MAXM = 135 + 1;
const int MAXN = 135 + 1; 
const int MAXK = 100 + 1;

int has_time_cost = 1;

clock_t begin_time;

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

int cost_max_index, can_full_use_time;

int pt[MAXT];

int base_cost;

double CalCost(int Ci, int Wi) {
    if (Wi == 0) return 0;
    if (Wi <= base_cost) return base_cost;
    double tmp = (double) (Wi - base_cost) * (double) (Wi - base_cost) / (double) Ci  + (double) Wi;
    return tmp;
}

// 计算此次放置会多出多少耗费
double CalCost(int Ci, int Wi, int cur, int ci) {
    return CalCost(Ci, max(Wi, cur + ci)) - CalCost(Ci, Wi);
}

struct Producer {
    string name;
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
    Producer() {
        for (int i = 0; i < MAXT; ++i) {
            is_full_use_time[i] = 0;
            need_time_cost_sum[i] = 0;
            time_cost[i] = 0;
        }
        can_visit_point_vec.clear();
        has_cost = 0;
        has_use_full_time = 0;
        lst_ab_has_cost = 0;
    }
    void Reset(int reset_has_cost, int reset_has_use_full_time) {
        if(reset_has_cost) has_cost = 0;
        if(reset_has_use_full_time) has_use_full_time = 0;
        for (int i = 1; i <= times; ++i) {
            if(reset_has_use_full_time) is_full_use_time[i] = 0;
            need_time_cost_sum[i] = ini_need_time_cost_sum[i];
            time_cost[i] = 0;
        }
    }
    bool operator < (const Producer& producer) const {
        return bandwidth < producer.bandwidth;
    }
    int TimeRemain(int time) {
        return bandwidth - time_cost[time];
    }
    int HasCostRemain(int time) {
        return has_cost - time_cost[time];
    }
    double GetAnswer(int update_has_cost, int update_full_use_time) {
        int *p = new int[MAXT];
        for (int time = 1; time <= times; ++time) {
            p[time] = time_cost[time];
        }
        sort(p + 1, p + 1 + times); 
        double ret = CalCost(bandwidth, p[cost_max_index]);
        if (update_has_cost) {
            has_cost = p[cost_max_index];
        }
        if (update_full_use_time) {
            for (int i = 1; i <= times; ++i) {
                is_full_use_time[i] = 0;
            }
            int num = 0;
            for (int i = 1; i <= times; ++i) {
                if (num == can_full_use_time) break;
                if (time_cost[i] <= p[cost_max_index]) continue;
                ++num;
                is_full_use_time[i] = 1;
            }
            int *indexs = new int[MAXT];
            // shuffle ?
            for (int time = 1; time <= times; ++time) indexs[time] = time;
            for (int i = 1; i <= times; ++i) {
                int time = indexs[i];
                if (num == can_full_use_time) break;
                if (time_cost[time] < p[cost_max_index]) continue;
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
    vector<int> can_visit_point_vec;
    int can_visit_point[MAXN];
    int time_need[MAXT][MAXK];
    int ini_time_need[MAXT][MAXK];
    Consumer() {
        can_visit_point_vec.clear();
    }
    void Reset() {
        for (int i = 1; i <= times; ++i) {
            for (int j = 1; j <= pt[i]; ++j)
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

// func
void Reset(int reset_has_cost, int reset_has_use_full_time);

// time, consumer_id, stream_id -> producer_id
int info_bandwidth[MAXT][36][MAXK]; 
int best_answer = 2000000000;

int GetAnswer() {
    double ret = 0.0;
    for (int i = 1; i <= producer_number; ++i) ret += producers[i].GetAnswer(0, 0);
    ret += 0.5;
    return (int)ret;
}

void Init() {
    srand(298758457);
    //srand((unsigned)time(NULL));
    ::cost_max_index = (times * 95 + 99) / 100;
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "cost_max_index: " << cost_max_index 
                      << endl;
    }
    ::can_full_use_time = times - cost_max_index;

    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            if (producers[producer_id].can_visit_point[consumer_id]) {
                producers[producer_id].can_visit_point_vec.emplace_back(consumer_id);
                consumers[consumer_id].can_visit_point_vec.emplace_back(producer_id);
                for (int time = 1; time <= times; ++time) {
                    for (int stream_id = 1; stream_id <= pt[time]; ++stream_id)
                        producers[producer_id].need_time_cost_sum[time] += consumers[consumer_id].time_need[time][stream_id];
                }
            }
        }
    }

    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int time = 1; time <= times; ++time) {
            for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
                time_node[time].sum_cost += consumers[consumer_id].time_need[time][stream_id];
                time_node[time].sum_value += consumers[consumer_id].time_need[time][stream_id];
            }
        }
    }

    // 暂时不需要初始化多次跑取最优

    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int time = 1; time <= times; ++time) {
            for (int stream_id = 1; stream_id <= pt[time]; ++stream_id)
                consumers[consumer_id].ini_time_need[time][stream_id] = consumers[consumer_id].time_need[time][stream_id];
        }
    }

    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int time = 1; time <= times; ++time) {
            producers[producer_id].ini_need_time_cost_sum[time] = producers[producer_id].need_time_cost_sum[time];
        }
    }
    
    for (int time = 1; time <= times; ++time) {
        time_node[time].ini_sum_cost = time_node[time].sum_cost;
    }

    // ab_candidate_map["pre_producer_sort"] = {0, 1, 2};
    // ab_candidate_map["pre_consumer_sort"] = {0, 1, 2};
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
        if (producers[producer_id].has_cost == 0) {
            producers[producer_id].has_cost = max(base_cost, producers[producer_id].time_cost[time]);    
        } 
        producers[producer_id].has_cost = max(producers[producer_id].has_cost, producers[producer_id].time_cost[time]);
    }

    if (nd_write == 0) return ;

    info_bandwidth[time][consumer_id][stream_id] = producer_id;
}

// 降序首次适应
bool DFF(int time, vector<RP>& consumer_vec, vector<P>& producer_vec, int update_use_cost=0, int nd_write=1) {
    sort(consumer_vec.begin(), consumer_vec.end());
    reverse(consumer_vec.begin(), consumer_vec.end());
    int *lst = new int[MAXN];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) lst[producer_id] = -1;
    for (auto &pp : producer_vec) {
        lst[pp.second] = pp.first;
    }
    int consumer_id, producer_id, stream_id, bandwidth;
    int best_producer_id, best_value;
    for (auto& cp : consumer_vec) {
        consumer_id = cp.second.first;
        stream_id = cp.second.second;
        bandwidth = consumers[consumer_id].time_need[time][stream_id];
        best_producer_id = -1;
        best_value = -1;
        for (auto& producer_id : consumers[consumer_id].can_visit_point_vec) {
            if (lst[producer_id] >= bandwidth) {
                if (lst[producer_id] > best_value) {
                    best_value = lst[producer_id];
                    best_producer_id = producer_id;
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


// // 降序最佳适应，优先考虑度数（先处理不太好放的
// bool DDFF(int time, vector<PP>& consumer_vec, vector<P>& producer_vec, int update_use_cost=0, int nd_write=1) {
//     sort(consumer_vec.begin(), consumer_vec.end());
    
//     int consumer_id, producer_id, stream_id, bandwidth;
//     int best_producer_id, best_value, best_degree;
//     for (auto& cp : consumer_vec) {
//         consumer_id = cp.second.first;
//         stream_id = cp.second.second;
//         bandwidth = consumers[consumer_id].time_need[time][stream_id];
//         best_producer_id = -1;
//         best_value = -1;
//         best_degree = 1000000;
//         for (auto& pp : producer_vec) {
//             producer_id = pp.second;
//             if (producers[producer_id].can_visit_point[consumer_id] == 0) continue;
//             if (producers[producer_id].TimeRemain(time) >= bandwidth) {
//                 if ((int)producers[producer_id].can_visit_point_vec.size() < best_degree ||
//                     (int)producers[producer_id].can_visit_point_vec.size() == best_degree && producers[producer_id].TimeRemain(time) > best_value) {
//                     best_value = producers[producer_id].TimeRemain(time);
//                     best_degree = (int)producers[producer_id].can_visit_point_vec.size();
//                     best_producer_id = producer_id;
//                 }
//             }
//         }
//         if (best_producer_id >= 1)
//             AddSomeBandWidth(time, best_producer_id, consumer_id, stream_id, bandwidth, update_use_cost, nd_write);    
//     }
//     return true;
// }
// // 降序最佳适应
// bool DBF(int time, vector<RP>& consumer_vec, vector<P>& producer_vec, int update_use_cost=0, int nd_write=1) {
//     sort(consumer_vec.begin(), consumer_vec.end());
//     reverse(consumer_vec.begin(), consumer_vec.end());

//     int consumer_id, producer_id, stream_id, bandwidth;
//     int best_producer_id, best_value;
//     for (auto& cp : consumer_vec) {
//         consumer_id = cp.second.first;
//         stream_id = cp.second.second;
//         bandwidth = consumers[consumer_id].time_need[time][stream_id];
//         best_producer_id = -1;
//         best_value = 1000000000;
//         for (auto& pp : producer_vec) {
//             producer_id = pp.second;
//             if (producers[producer_id].can_visit_point[consumer_id] == 0) continue;
//             if (producers[producer_id].TimeRemain(time) >= bandwidth) {
//                 if (producers[producer_id].TimeRemain(time) < best_value) {
//                     best_value = producers[producer_id].TimeRemain(time);
//                     best_producer_id = producer_id;
//                 }
//             }
//         }
//         if (best_producer_id >= 1)
//             AddSomeBandWidth(time, best_producer_id, consumer_id, stream_id, bandwidth, update_use_cost, nd_write);    
//     }
//     return true;
// }

void PreWork() {
    clock_t func_begin_time = clock();
    
    vector<pair<int,int> > consumer_tmp;  
    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        consumer_tmp.emplace_back(make_pair(consumers[consumer_id].can_visit_point_vec.size(), consumer_id));
    }
    // 客户从度数小的先开始
    sort(consumer_tmp.begin(), consumer_tmp.end());
    int *ban = new int[MAXT];
    for (int i = 1; i <= times; ++i) ban[i] = 0;
    for (int k = 1; k <= can_full_use_time * producer_number; ++k) {
        int best_producer_id = -1, best_time = -1, can_min_last = -1, time_sum_use = -1, node_degree = -1;
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
            int can_last = min(producers[producer_id].need_time_cost_sum[best_time], producers[producer_id].TimeRemain(best_time));
            // 按value来排(can_cost)
            if (can_last > can_min_last ) {
                can_min_last = can_last;
                best_producer_id = producer_id;
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
        for (auto& cp : consumer_tmp) {
            int consumer_id = cp.second;
            if (producers[best_producer_id].can_visit_point[consumer_id] == 0) continue;
            for (int stream_id = 1; stream_id <= pt[best_time]; ++stream_id) {
                int bandwidth = producers[best_producer_id].TimeRemain(best_time);
                if (bandwidth >= consumers[consumer_id].time_need[best_time][stream_id]) {
                    bandwidth = consumers[consumer_id].time_need[best_time][stream_id];
                    AddSomeBandWidth(best_time, best_producer_id, consumer_id, stream_id, bandwidth, 0, 0);
                }
            }
        }

        // update time_node sum value
        time_node[best_time].sum_value += producers[best_producer_id].lst_ab_has_cost;
        
        // update producer
        producers[best_producer_id].is_full_use_time[best_time] = 1;
        ++producers[best_producer_id].has_use_full_time;

        if (debug_file != nullptr && k <= 100) {
            (*debug_file) << "prework full_use, producer_id: " << best_producer_id
                          << " time: " << best_time 
                          << " his number: " << producers[best_producer_id].has_use_full_time 
                          << " last: " << can_min_last
                          << " degree: " << node_degree
                          << " time_sum_cost: " << time_sum_use <<endl;
        }
    }

    Reset(0, 0);
    vector<P> producer_vec;
    vector<RP> consumer_vec;
    for(int time = 1; time <= times; ++time) {
        producer_vec.clear();
        consumer_vec.clear();
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].is_full_use_time[time]) {
                producer_vec.emplace_back(P(producers[producer_id].bandwidth, producer_id));
            }
        }
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            for (int stream_id =1; stream_id <= pt[time]; ++stream_id) {
                consumer_vec.emplace_back(RP(consumers[consumer_id].time_need[time][stream_id], P(consumer_id, stream_id)));
            }
        }
        DFF(time, consumer_vec, producer_vec, 0, 1);
    }

    if (debug_file != nullptr) {
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
    delete[] ban;
    if (debug_file != nullptr) {
        (*debug_file) << "prework_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }
}




int vis_time[MAXT];
void WorkTimeBaseline(int time, int time_index) {
    clock_t func_begin_time = clock();
    // todo
    vector<P> producer_vec;
    for (int i = 1; i <= producer_number; ++i) {
        producer_vec.emplace_back(P(producers[i].bandwidth, i));
    }
    sort(producer_vec.begin(), producer_vec.end());
    reverse(producer_vec.begin(), producer_vec.end());
    
    // vector<P> consumer_vec;
    // for (int i = 1; i <= consumer_number; ++i) {
    //     consumer_vec.emplace_back(P(consumers[i].can_visit_point_vec.size(), i));
    // }
    // sort(consumer_vec.begin(), consumer_vec.end());
    
    // 先放不花钱的
    vector<P> f_producer_vec;
    vector<RP> f_consumer_vec;
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].is_full_use_time[time] == 0) {
            f_producer_vec.emplace_back(P(producers[producer_id].has_cost, producer_id));
        }
    }
    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int stream_id =1; stream_id <= pt[time]; ++stream_id) {
            if (consumers[consumer_id].time_need[time][stream_id] > 0)
                f_consumer_vec.emplace_back(RP(consumers[consumer_id].time_need[time][stream_id], P(consumer_id, stream_id)));
        }
    }
    DFF(time, f_consumer_vec, f_producer_vec, 0, 1);
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

    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int stream_id = 1; stream_id <= pt[time]; ++stream_id) {
            int nd_bandwidth = consumers[consumer_id].time_need[time][stream_id];
            if (nd_bandwidth == 0) continue ;
            int best_producer_id = -1;
            double min_cost = 1e12;            
            for(auto& pp : producer_vec) {
                int producer_id = pp.second;
                if (producers[producer_id].can_visit_point[consumer_id] == 0) continue;
                if (producers[producer_id].TimeRemain(time) < nd_bandwidth) continue;
                double tmp = (int)CalCost(producers[producer_id].bandwidth, producers[producer_id].has_cost, producers[producer_id].time_cost[time], nd_bandwidth);
                if (tmp < min_cost) {
                    min_cost = tmp;
                    best_producer_id = producer_id;
                }
            }
            AddSomeBandWidth(time, best_producer_id, consumer_id, stream_id, nd_bandwidth, 1, 1);
        }
    }

    if (debug_file != nullptr && time_index <= 10) {
        (*debug_file) << "workbaseline_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }

}

void EndWork();

void Work() {
    PreWork();
    // todo
    vector<P> time_vec;
    for (int time = 1; time <= times; ++time) time_vec.emplace_back(P(time_node[time].sum_cost, time));
    // sort(time_vec.begin(), time_vec.end());
    // reverse(time_vec.begin(), time_vec.end());
    for (int index = 0; index < times; ++index) {
        int time = time_vec[index].second;
        WorkTimeBaseline(time, index + 1);
    }
}

bool Check(int nd_write) {
    vector<P> f_producer_vec;
    vector<RP> f_consumer_vec;  
    for (int time = 1; time <= times; ++time) {
        f_consumer_vec.clear();
        f_producer_vec.clear();
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].is_full_use_time[time] == 0) {
                f_producer_vec.emplace_back(P(producers[producer_id].has_cost, producer_id));
            } else {
                f_producer_vec.emplace_back(P(producers[producer_id].bandwidth, producer_id));
            }
        }
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            for (int stream_id =1; stream_id <= pt[time]; ++stream_id) {
                if (consumers[consumer_id].time_need[time][stream_id] > 0)
                    f_consumer_vec.emplace_back(RP(consumers[consumer_id].time_need[time][stream_id], P(consumer_id, stream_id)));
            }
        }
        DFF(time, f_consumer_vec, f_producer_vec, 0, nd_write);
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            for (int stream_id =1; stream_id <= pt[time]; ++stream_id) {
                if (consumers[consumer_id].time_need[time][stream_id] > 0)
                    return false;
            }
        }
    }
    return true;
}

void WorkBinarySearch() {
    PreWork();

    double lb = 0.0, ub = 100.0, mid, best;
    for (int k = 1; k <= 15; ++k) {
        mid = (lb + ub) / 2.0;
        Reset(1, 0);
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (base_cost > producers[producer_id].bandwidth) {
                producers[producer_id].has_cost = producers[producer_id].bandwidth;
            } else {
                producers[producer_id].has_cost = min(producers[producer_id].bandwidth, 
                                                      int(mid * (double)producers[producer_id].bandwidth + (1.0 - mid) * (double)base_cost)); 
            }
        }
        if (Check(0)) {
            best = mid; 
            ub = mid;
        } else {
            lb = mid;
        }
        if (debug_file != nullptr) {
            (*debug_file) << "WorkBinarySearch, lb: " <<lb
                          << " ub: " <<ub
                          << endl;
        }
    }
    Reset(1, 0);
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (base_cost > producers[producer_id].bandwidth) {
            producers[producer_id].has_cost = producers[producer_id].bandwidth;
        } else {
            producers[producer_id].has_cost = min(producers[producer_id].bandwidth, 
                                                    int(best * (double)producers[producer_id].bandwidth + (1.0 - best) * (double)base_cost)); 
        }
    }
    Check(1);
}

void EndWork() {
    clock_t func_begin_time = clock();
    
    // todo
    
    if (debug_file != nullptr) {
        (*debug_file) << "endwork_time_cost: " << (double)(clock() - func_begin_time) / CLOCKS_PER_SEC << "s" 
                      << endl;
    }
}

void Reset(int reset_has_cost, int reset_has_use_full_time) {
    for (int i = 1; i <= times; ++i) {
        time_node[i].Reset();
    }
    for (int i = 1; i <= producer_number; ++i) {
        producers[i].Reset(reset_has_cost, reset_has_use_full_time);
    }
    for (int i = 1; i <= consumer_number; ++i) {
        consumers[i].Reset();
    }
}

void ABWork() {
    for (auto& p : ab_candidate_map) ab_map[p.first] = p.second[0];
    for (auto& p : ab_candidate_map) {
        best_answer = 2000000000;
        int best_v = -1;
        for (auto& v : p.second)
        {
            ab_map[p.first] = v;
            Work();
            int answer = GetAnswer();
            if (answer < best_answer) {
                best_answer = answer;
                best_v = v;
            }
        }
        ab_map[p.first] = best_v;
    }    
}

int main(int argc, char *argv[]) {
    begin_time = clock();
    // local or oj
    ::info = "";
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
    // is_ab = 1;
    // for (int k = 1; k <= 10; ++k) Work();
    // ABWork();
    is_ab = 0;
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
        (*result_file) << "sum_cost: " << GetAnswer() << " cost_time: " << (double)(clock() - begin_time) / CLOCKS_PER_SEC << "s" << endl;
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
                if (info_bandwidth[time][consumer_id][stream_id] <= 0) {
                    continue ;
                }
                int producer_id = info_bandwidth[time][consumer_id][stream_id];
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