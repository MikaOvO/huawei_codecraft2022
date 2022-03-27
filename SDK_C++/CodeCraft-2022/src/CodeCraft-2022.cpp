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

ofstream *debug_file = nullptr;
ofstream *result_file = nullptr;

string output_dir = "/output";
string data_dir = "/data";
string info = "";

const int MAXT = 8928 + 1;
const int MAXM = 135 + 1;
const int MAXN = 135 + 1; 

clock_t begin_time;

// algorithm
// 最小费用流
struct MinCostFlow {
    const static int MAXN = 135 * 2;
    const static int MAXM = 135 * 35 * 4;
    const static int INF = 1e9;
    struct Edge
    {
        int from, to, cap, flow, cost, nxt;
    };
    Edge edge[MAXM];
    int head[MAXN], edgenum;
    int pre[MAXN];//记录增广路径上 到达点i的边的编号
    int dist[MAXN];
    bool vis[MAXN];
    int source, sink;//超级源点 超级汇点
    void init()
    {
        edgenum = 0;
        memset(head, -1, sizeof(head));
    }
    void addEdge(int u, int v, int w, int c)
    {
        Edge E1 = {u, v, w, 0, c, head[u]};
        edge[edgenum] = E1;
        head[u] = edgenum++;
        Edge E2 = {v, u, 0, 0, -c, head[v]};
        edge[edgenum] = E2;
        head[v] = edgenum++;
    }
    bool SPFA(int s, int t)//寻找花销最少的路径
    {
        //跑一遍SPFA 找s--t的最少花销路径 且该路径上每一条边不能满流
        //若存在 说明可以继续增广，反之不能
        queue<int> Q;
        for(int i = 0; i < MAXN; ++i) dist[i] = INF;
        memset(vis, false, sizeof(vis));
        memset(pre, -1, sizeof(pre));
        dist[s] = 0;
        vis[s] = true;
        Q.push(s);
        while(!Q.empty())
        {
            int u = Q.front();
            Q.pop();
            vis[u] = false;
            for(int i = head[u]; i != -1; i = edge[i].nxt)
            {
                Edge E = edge[i];
                if(dist[E.to] > dist[u] + E.cost && E.cap > E.flow)//可以松弛 且 没有满流
                {
                    dist[E.to] = dist[u] + E.cost;
                    pre[E.to] = i;//记录前驱边 的编号
                    if(!vis[E.to])
                    {
                        vis[E.to] = true;
                        Q.push(E.to);
                    }
                }
            }
        }
        return pre[t] != -1;//可达返回true
    }
    void MCMF(int s, int t, int &cost, int &flow)
    {
        flow = 0;//总流量
        cost = 0;//总费用
        while(SPFA(s, t))//每次寻找花销最小的路径
        {
            int Min = INF;
            //通过反向弧 在源点到汇点的最少花费路径 找最小增广流
            for(int i = pre[t]; i != -1; i = pre[edge[i^1].to])
            {
                Edge E = edge[i];
                Min = min(Min, E.cap - E.flow);
            }
            //增广
            for(int i = pre[t]; i != -1; i = pre[edge[i^1].to])
            {
                edge[i].flow += Min;
                edge[i^1].flow -= Min;
                cost += edge[i].cost * Min;//增广流的花销
            }
            flow += Min;//总流量累加
        }
    }
};

// 最大流
struct MaxFlow {
    const static int maxm = 135 * 35 * 4;
    const static int maxn = 135 * 2;
    const static int INF = 1e9;
    struct Edge {
        int u,v,f,nxt;
    };
    Edge edges[maxm];
    int head[maxn],tot;
    void clear() {
        tot = 0;
        for(int i=0;i<maxn;i++) head[i]=-1;
    }
    void addEdge(int u,int v,int f) {
        edges[tot].u=u;edges[tot].v=v;edges[tot].f=f;
        edges[tot].nxt=head[u];head[u]=tot++;
        edges[tot].u=v;edges[tot].v=u;edges[tot].f=0;
        edges[tot].nxt=head[v];head[v]=tot++;
    }

    int depth[maxn],cur[maxn];
    bool bfs(int S,int T) // 分层图
    {
        queue<int> Q;
        memset(depth,0,sizeof(depth));
        depth[S]=1;
        Q.push(S);
        while(!Q.empty())
        {
            int u=Q.front();Q.pop();
            for(int i=head[u];~i;i=edges[i].nxt)
            {
                if(edges[i].f>0 && depth[edges[i].v]==0)
                {
                    depth[edges[i].v]=depth[u]+1;
                    Q.push(edges[i].v);
                    if(edges[i].v==T) return true;
                }
            }
        }
        return false;
    }

    int dfs(int u,int T,int flow) // 找增广路
    {
        if(u==T) return flow;
        for(int &i=cur[u];~i;i=edges[i].nxt) // &符号，i改变同时也改变cur[u]的值，达到记录当前弧的目的
        {
            if(depth[edges[i].v]==depth[u]+1 && edges[i].f!=0)
            {
                int sub=dfs(edges[i].v,T,min(flow,edges[i].f)); // 向下增广
                if(sub>0)
                {
                    edges[i].f-=sub;
                    edges[i^1].f+=sub;
                    return sub;
                }
            }
        }
        return 0;
    }

    int Dinic(int n,int S,int T)
    {
        int ans=0;
        while(bfs(S,T))
        {
            for(int i=0;i<=n;i++) // //每一次建立完分层图后都要把cur置为每一个点的第一条边
                cur[i]=head[i];
            while(int f=dfs(S,T,INF))
                ans+=f;
        }
        return ans;
    }
    //n最多的点数，S为源点，T为汇点
    /*
    tot=0;
    memset(head,-1,sizeof(head));
    addEdge(u,v,c);
    */
};


/*
 * 使用前调用函数vector<int> ProduceHasCost(vector<int> times)即可
 */

// AB
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

// number
int times, producer_number, consumer_number;

int cost_max_index, can_full_use_time;

struct Producer {
    string name;
    int bandwidth;
    int has_cost;
    int has_use_full_time;
    long long waste;
    int can_visit_point[MAXM];
    vector<int> can_visit_point_vec;
    int need_time_cost_sum[MAXT];
    int is_full_use_time[MAXT];
    int time_cost[MAXT];
    int ini_need_time_cost_sum[MAXT];
    int flow_is_full_use_time[MAXT];
    priority_queue<int, vector<int>, greater<int> > cost_que;
    Producer() {
        for (int i = 0; i < MAXT; ++i) {
            is_full_use_time[i] = 0;
            flow_is_full_use_time[i] = 0;
            need_time_cost_sum[i] = 0;
            time_cost[i] = 0;
        }
        can_visit_point_vec.clear();
        has_cost = 0;
        has_use_full_time = 0;
    }
    void Reset() {
        has_cost = 0;
        has_use_full_time = 0;
        for (int i = 1; i <= times; ++i) {
            is_full_use_time[i] = 0;
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
    // 1e6 * 1e4
    void CalFlowNeed() {
        for (int i = 1; i <= times; ++i) {
            flow_is_full_use_time[i] = 0;
        }
        int *p = new int[MAXT];
        for (int time = 1; time <= times; ++time) {
            p[time] = time_cost[time];
        }
        sort(p + 1, p + 1 + times);
        has_cost = p[cost_max_index];
        int num = 0;
        // QAQ
        for (int time = 1; time <= times; ++time) {
            if (num == can_full_use_time) break;
            if (time_cost[time] <= p[cost_max_index]) continue;
            ++num;
            flow_is_full_use_time[time] = 1;
        }
        int *indexs = new int[MAXT];
        for (int time = 1; time <= times; ++time) indexs[time] = time;
        for (int i = 1; i <= times; ++i) {
            int time = indexs[i];
            if (num == can_full_use_time) break;
            if (time_cost[time] < p[cost_max_index]) continue;
            ++num;
            flow_is_full_use_time[time] = 1;
        }
        delete[] p;
    }
    long long GetWaste(int update_use_cost=0) {
        // 重新计算is_full_use_time
        for (int i = 1; i <= times; ++i) {
            is_full_use_time[i] = 0;
        }
        int *p = new int[MAXT];
        for (int time = 1; time <= times; ++time) {
            p[time] = time_cost[time];
        }
        sort(p + 1, p + 1 + times);
        long long ret = 0;
        for (int time = 1; time <= cost_max_index; ++time) {
            ret += p[time];
        }
        if (update_use_cost) {
            has_cost = p[cost_max_index];
        }
        int num = 0;
        // QAQ
        for (int time = 1; time <= times; ++time) {
            if (num == can_full_use_time) break;
            if (time_cost[time] <= p[cost_max_index]) continue;
            ++num;
            is_full_use_time[time] = 1;
        }
        for (int time = 1; time <= times; ++time) {
            if (num == can_full_use_time) break;
            if (time_cost[time] < p[cost_max_index]) continue;
            ++num;
            is_full_use_time[time] = 1;
        }
        waste = (long long)cost_max_index * (long long)p[cost_max_index] - ret;
        delete[] p;
        return waste;
    }
} producers[MAXN];

struct TimeNode {
    int sum_cost;
    int ini_sum_cost;
    TimeNode() {
        sum_cost = 0;
    }
    void Reset() {
        sum_cost = ini_sum_cost;
    }
} time_node[MAXT];

struct Consumer {
    string name;
    vector<int> can_visit_point_vec;
    int can_visit_point[MAXN];
    int time_need[MAXT];
    int ini_time_need[MAXT];
    Consumer() {
        can_visit_point_vec.clear();
        for (int i = 1; i <= times; ++i) {
            time_need[i] = 0;
        }
    }
    void Reset() {
        for (int i = 1; i <= times; ++i) {
            time_need[i] = ini_time_need[i];
        }
    }
    void NodifyProducer(int time, int cost_bandwidth) {
        for (auto& producer_id : can_visit_point_vec) {
            producers[producer_id].need_time_cost_sum[time] -= cost_bandwidth;
        }
    }
    void NodifyTimeNode(int time, int cost_bandwidth) {
        time_node[time].sum_cost -= cost_bandwidth;
    }
} consumers[MAXM];

// key-(producer_id, consumer_id) value-bandwidth
int info_bandwidth[MAXT][MAXN][36]; 
int best_answer = 2000000000;


void Init() {
    srand(131);
    // srand((unsigned)time(NULL));
    ::cost_max_index = (times * 95 + 99) / 100;
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "cost_max_index: " << cost_max_index <<endl;
    }
    ::can_full_use_time = times - cost_max_index;

    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            if (producers[producer_id].can_visit_point[consumer_id]) {
                producers[producer_id].can_visit_point_vec.emplace_back(consumer_id);
                consumers[consumer_id].can_visit_point_vec.emplace_back(producer_id);
                for (int time = 1; time <= times; ++time) {
                    producers[producer_id].need_time_cost_sum[time] += consumers[consumer_id].time_need[time];
                }
            }
        }
    }

    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int time = 1; time <= times; ++time) {
            time_node[time].sum_cost += consumers[consumer_id].time_need[time];
        }
    }

    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int time = 1; time <= times; ++time) {
            consumers[consumer_id].ini_time_need[time] = consumers[consumer_id].time_need[time];
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

    ab_candidate_map["pre_producer_sort"] = {0, 1, 2};
    ab_candidate_map["pre_consumer_sort"] = {0, 1, 2};
    ab_candidate_map["worktime_producer_sort"] = {0, 1, 2};
    ab_candidate_map["worktime_consumer_sort"] = {0, 1, 2};
    ab_candidate_map["work_time_sort"] = {0, 1, 2};
}

void MoveSomeBandWidth(int time, int from_producer_id, int to_producer_id, int consumer_id, int bandwidth, int update_use_cost=0) {
    if (debug_file != nullptr && is_ab == 0) {
        if ( time <= 2 && bandwidth > 0 ) {
            (*debug_file) << "MoveSomeBandWidth time: " << time << " from_producer_id: " << from_producer_id
                          << " to_producer_id: " << to_producer_id << " bandwidth: "<< bandwidth <<endl;
        }
    }
    producers[to_producer_id].time_cost[time] += bandwidth;
    producers[from_producer_id].time_cost[time] -= bandwidth;

    // not nodify

    info_bandwidth[time][from_producer_id][consumer_id] -= bandwidth;
    info_bandwidth[time][to_producer_id][consumer_id] += bandwidth;
}

void AddSomeBandWidth(int time, int producer_id, int consumer_id, int bandwidth, int update_use_cost=0, int nd_write=1) {
    if (bandwidth == 0) return ;
    if (debug_file != nullptr && is_ab == 0) {
        if (time <= 5 && producer_id <= 5 && consumer_id <= 5) {
            (*debug_file) << "AddSomeBandWidth time: " << time << " producer_id: " << producer_id
                          << " consumer_id: " << consumer_id << " bandwidth: "<< bandwidth <<endl;
        }
    }
    // 把耗费放到当天里，注意策略里不需要写这一块了
    producers[producer_id].time_cost[time] += bandwidth;
    consumers[consumer_id].time_need[time] -= bandwidth;
    consumers[consumer_id].NodifyProducer(time, bandwidth);
    consumers[consumer_id].NodifyTimeNode(time, bandwidth);

    // 更新一下理论花费（实际上比这要大，因为预先分配的不一定真的那最高的5%
    if (update_use_cost && producers[producer_id].is_full_use_time[time] == 0) {
        producers[producer_id].has_cost = max(producers[producer_id].has_cost, producers[producer_id].time_cost[time]);
    }

    if (nd_write == 0) return ;

    info_bandwidth[time][producer_id][consumer_id] += bandwidth;
}

int GetAnswer() {
    int sum_cost = 0;
    for (int i = 1; i <= producer_number; ++i) {
        sort(producers[i].time_cost + 1, producers[i].time_cost + 1 + times);
        sum_cost += producers[i].time_cost[cost_max_index];
    }
    return sum_cost;
}

void TryFullUse() {
    vector<pair<int,int> > tmp;  
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        tmp.emplace_back(make_pair(producers[producer_id].can_visit_point_vec.size(), producer_id));
    }

    vector<pair<int,int> > consumer_tmp;  
    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        consumer_tmp.emplace_back(make_pair(consumers[consumer_id].can_visit_point_vec.size(), consumer_id));
    }
    // 客户从度数小的先开始
    sort(consumer_tmp.begin(), consumer_tmp.end());
    // 边缘从度数小的开始
    sort(tmp.begin(), tmp.end());
    int can_all_time_full_use_flag = 0;
    for (auto &p : tmp) {
        int producer_id = p.second;
        can_all_time_full_use_flag = 1;
        for (int time = 1; time <= times; ++time) {
            if (producers[producer_id].need_time_cost_sum[time] < producers[producer_id].TimeRemain(time)) {
                can_all_time_full_use_flag = 0;
                break ;
            }
        }
        if (can_all_time_full_use_flag == 0) {
            continue ;
        }

        int block = 5000;
        for (int time = 1; time <= times; ++time) {
            for (int k = 1; k <= 1000; ++k) {
                if (producers[producer_id].TimeRemain(time) == 0) {
                    break ;
                }
                for (auto& p : consumer_tmp) {
                    int consumer_id = p.second;
                    if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
                        continue;
                    }
                    int cost_bandwidth = min(producers[producer_id].TimeRemain(time), consumers[consumer_id].time_need[time]);
                    // cost_bandwidth = min(cost_bandwidth, block);
                    if (cost_bandwidth == 0) continue;
                    AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth, 1, 0);   
                }
            }
        }

        if (debug_file != nullptr && is_ab == 0) {
            (*debug_file) << "all_time_full_use: producer_id: " << producer_id << endl;
        }    
    }
}

void PreWork() {
    // prework
    // vector<pair<pair<int,int>,int> > can_cost_time;
    // vector<pair<int,int> > tmp;  
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     tmp.emplace_back(make_pair(producers[producer_id].can_visit_point_vec.size(), producer_id));
    // }
    // vector<pair<int,int> > consumer_tmp;  
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     consumer_tmp.emplace_back(make_pair(consumers[consumer_id].can_visit_point_vec.size(), consumer_id));
    // }
    // // 客户从度数小的先开始
    // sort(consumer_tmp.begin(), consumer_tmp.end());
    // if (ab_map["pre_consumer_sort"] == 1) reverse(consumer_tmp.begin(), consumer_tmp.end());
    // if (ab_map["pre_consumer_sort"] == 2) random_shuffle(consumer_tmp.begin(), consumer_tmp.end());
    // // 边缘从度数小的开始
    // sort(tmp.begin(), tmp.end());
    // if (ab_map["pre_producer_sort"] == 1) reverse(tmp.begin(), tmp.end());
    // if (ab_map["pre_producer_sort"] == 2) random_shuffle(tmp.begin(), tmp.end());
    
    // for (auto& p : tmp) {
    //     int producer_id = p.second;
    //     can_cost_time.clear();
    //     for (int time = 1; time <= times; ++time) {
    //         int can_cost = 0;
    //         int all_cost = 0;
    //         for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //             // 注释掉以后就是统计当天所有花费了，先把花费高的搞一搞
    //             if (consumers[consumer_id].can_visit_point[producer_id] == 1) {
    //                 can_cost += consumers[consumer_id].time_need[time];
    //             }
    //             all_cost += consumers[consumer_id].time_need[time];
    //         }
    //         can_cost = min(can_cost, producers[producer_id].bandwidth);
    //         can_cost_time.emplace_back(make_pair(make_pair(-can_cost, -all_cost), time));
    //     }
    //     sort(can_cost_time.begin(), can_cost_time.end());
    //     // reverse(can_cost_time.begin(), can_cost_time.end());
    //     for (int i = 0; i < can_full_use_time; ++i) {
    //         int time = can_cost_time[i].second;
    //         producers[producer_id].is_full_use_time[time] = 1;
    //         // 也许这里平均分也会好一点？
    //         int block = max(1, (producers[producer_id].bandwidth - producers[producer_id].time_cost[time]) / 1000);
    //         block = 5000;
    //         for (int k = 1; k <= 1; ++k) {
    //             for (auto& cp : consumer_tmp) {
    //                 int consumer_id = cp.second;
    //                 if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
    //                     continue;
    //                 }
    //                 int cost_bandwidth = min(producers[producer_id].bandwidth - producers[producer_id].time_cost[time],
    //                                     consumers[consumer_id].time_need[time]);
    //                 // cost_bandwidth = min(cost_bandwidth, block);
    //                 if (cost_bandwidth == 0) continue;
    //                 AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth);               
    //             }
    //         }
    //     }
    // }

    // 一种感觉上更优的方法去处理5%
    // 10000^2 / 20 = 5e6, 5e6 * 30 * 100 ~ 1e10 不可以接受
    // 优化了一下，复杂度除以30
    vector<pair<int,int> > consumer_tmp;  
    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        consumer_tmp.emplace_back(make_pair(consumers[consumer_id].can_visit_point_vec.size(), consumer_id));
    }
    // 客户从度数小的先开始
    sort(consumer_tmp.begin(), consumer_tmp.end());

    for (int k = 1; k <= can_full_use_time * producer_number; ++k) {
        int best_producer_id = -1, best_time = -1, can_min_last = -1, time_sum_use = -1, node_degree = -1;
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (producers[producer_id].has_use_full_time == can_full_use_time) {
                continue ;
            }
            for (int time = 1; time <= times; ++time) {
                if (producers[producer_id].bandwidth != producers[producer_id].TimeRemain(time)) {
                    continue ;
                }
                int can_last = min(producers[producer_id].need_time_cost_sum[time], producers[producer_id].TimeRemain(time));
                if (time_sum_use < time_node[time].sum_cost || 
                   //(time_sum_use == time_node[time].sum_cost && can_last > can_min_last) ||  
                   (time_sum_use == time_node[time].sum_cost && can_last > can_min_last) ) {
                    time_sum_use = time_node[time].sum_cost;
                    can_min_last = can_last;
                    best_producer_id = producer_id;
                    best_time = time;
                }
            }
        }
        if (time_sum_use <= 0) {
            break ;
        }
        producers[best_producer_id].is_full_use_time[best_time] = 1;
        ++producers[best_producer_id].has_use_full_time;
        if (debug_file != nullptr && k <= 1000) {
            (*debug_file) << "prework full_use, producer_id: " << best_producer_id
                          << " time: " << best_time 
                          << " his number: " << producers[best_producer_id].has_use_full_time 
                          << " last: " << can_min_last
                          << " degree: " << node_degree
                          << " time_sum_cost: " << time_sum_use <<endl;
        }
        int block = 2000;
        for (int k = 1; k <= 1000; ++k) {
            if (producers[best_producer_id].TimeRemain(best_time) == 0) {
                break ;
            }
            for (auto& p : consumer_tmp) {
                int consumer_id = p.second;
                if (consumers[consumer_id].can_visit_point[best_producer_id] == 0) {
                    continue;
                }
                int cost_bandwidth = min(producers[best_producer_id].TimeRemain(best_time), consumers[consumer_id].time_need[best_time]);
                // cost_bandwidth = min(cost_bandwidth, block);
                if (cost_bandwidth == 0) continue;
                AddSomeBandWidth(best_time, best_producer_id, consumer_id, cost_bandwidth, 0, 1);   
            }
        }
    }
    if (debug_file != nullptr) {
        int *tmp_time_cost = new int[MAXT];
        int num = 0;
        (*debug_file) <<"prework5%:"<< endl;
        for (int i = 1; i <= producer_number; ++i) {
            num = 0;
            for (int time = 1; time <= times; ++time) tmp_time_cost[time] = producers[i].time_cost[time];
            sort(tmp_time_cost + 1, tmp_time_cost + 1 + times);
            (*debug_file) << "producer_id: " << i << "bandwidth: " << producers[i].bandwidth <<endl;
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


    // 尝试更平均一点？
    // for (int k = 1; k <= 10; ++k) {
    //     int max_time_id, max_time = -1;
    //     for (int i = 1; i <= times; ++i) {
    //         if (time_node[i].sum_cost > max_time) {
    //             max_time = time_node[i].sum_cost;
    //             max_time_id = i;
    //         }
    //     }
    //     for (int time = 1; time <= times; ++time) {
    //         if (time == max_time_id) continue;
    //         for (int producer_id = 1; producer_id <= producer_number; ++producer_number) {
    //             if ()
    //         }
    //     }
    // }
}

// has fix
void WorkTimeMaxFlow(int time,int only_use_has_cost=0) {
    MaxFlow mf;
    mf.clear();
    int source = 0, target = consumer_number + producer_number + 1;
    // run full use (not cost)
    // for (int i = 1; i <= producer_number; ++i) {
    //     if (producers[i].is_full_use_time[time]) {
    //         mf.addEdge(source, i, producers[i].TimeRemain(time));
    //     }
    // }
    // for (int i = 1; i <= consumer_number; ++i) {
    //    mf.addEdge(producer_number + i, target, consumers[i].time_need[time]);
    // }
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //         if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
    //             continue;
    //         }
    //         mf.addEdge(producer_id, producer_number + consumer_id, mf.INF);
    //     }
    // }
    // mf.Dinic(target, source, target);
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     for (int i = mf.head[producer_number + consumer_id]; i != -1; i = mf.edges[i].nxt) {
    //         if (1 <= mf.edges[i].v && mf.edges[i].v <= producer_number && mf.edges[i].f > 0) {
    //             AddSomeBandWidth(time, mf.edges[i].v, consumer_id, mf.edges[i].f);
    //         }
    //     }
    // }

    // run other
    mf.clear();
    for (int i = 1; i <= producer_number; ++i) {
        if (only_use_has_cost == 0)
            mf.addEdge(source, i, producers[i].TimeRemain(time));
        else
            mf.addEdge(source, i, min(producers[i].TimeRemain(time), producers[i].has_cost));
    }
    for (int i = 1; i <= consumer_number; ++i) {
       mf.addEdge(producer_number + i, target, consumers[i].time_need[time]);
    }
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
                continue;
            }
            mf.addEdge(producer_id, producer_number + consumer_id, mf.INF);
        }
    }
    mf.Dinic(target, source, target);
    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        for (int i = mf.head[producer_number + consumer_id]; i != -1; i = mf.edges[i].nxt) {
            if (1 <= mf.edges[i].v && mf.edges[i].v <= producer_number && mf.edges[i].f > 0) {
                AddSomeBandWidth(time, mf.edges[i].v, consumer_id, mf.edges[i].f);
            }
        }
    }
}

void AverageAdd(int time, int consumer_id, int is_only_use_has_cost) {
    vector<pair<int,int> > vec; // can_cost, producer_id
    int tmp;
    int cost_bandwidth;
    for (auto& producer_id : consumers[consumer_id].can_visit_point_vec) {
        vec.emplace_back(make_pair(producers[producer_id].can_visit_point_vec.size(), producer_id));
    }   
    // 这里排序很玄学
    // sort(vec.begin(), vec.end());
    // reverse(vec.begin(), vec.end());
    int producer_id;
    for (int i = 0; i < vec.size(); ++i) {
        producer_id = vec[i].second;
        if (is_only_use_has_cost)
            tmp = max(0, producers[producer_id].has_cost - producers[producer_id].time_cost[time]);
        else
            tmp = producers[producer_id].bandwidth - producers[producer_id].time_cost[time];
        cost_bandwidth = min(tmp, (consumers[consumer_id].time_need[time] + ((int)vec.size() - i - 1)) / ((int)vec.size() - i));
        int producer_id = vec[i].second;
        if (is_only_use_has_cost)
            AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth, 0, 1);
        else
            AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth, 1, 1); 
    } 
    reverse(vec.begin(), vec.end());
    for (int i = 0; i < vec.size(); ++i) {
        if (consumers[consumer_id].time_need[time] == 0) break;
        producer_id = vec[i].second;
        if (is_only_use_has_cost)
            tmp = max(0, producers[producer_id].has_cost - producers[producer_id].time_cost[time]);
        else
            tmp = producers[producer_id].bandwidth - producers[producer_id].time_cost[time];
        cost_bandwidth = min(tmp, consumers[consumer_id].time_need[time]);
        int producer_id = vec[i].second;
        if (is_only_use_has_cost)
            AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth, 0, 1);
        else
            AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth, 1, 1); 
    } 
}

int vis_time[MAXT];
void WorkTimeBaseline(int time, int time_index) {
    for (int i = 1; i <= times; ++i) {
        vis_time[i] = 0;
    }
    if (debug_file != nullptr && is_ab == 0 && time_index <= 1) {
        (*debug_file) << "time: " << time << " " << time_index << endl; 
        int sum_cost = 0;
        for (int i = 1; i <= consumer_number; ++i) sum_cost += consumers[i].time_need[time];
        (*debug_file) << "need cost: " << sum_cost << " " << time_node[time].sum_cost << endl; 
    }
    vector<pair<int,int> > can_cost_time;
    vector<pair<pair<int,int>,int> > tmp;  
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        int min_sum_cost = 1000000000;
        for (int time = 1; time <= times; ++time) {
            if(vis_time[time] == 0)
                min_sum_cost = min(min_sum_cost, producers[producer_id].need_time_cost_sum[time]);
        }
        tmp.emplace_back(make_pair(make_pair(producers[producer_id].can_visit_point_vec.size(), 0)
                                    , producer_id));
    }
    vector<pair<int,int> > consumer_tmp;  
    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        consumer_tmp.emplace_back(make_pair(consumers[consumer_id].can_visit_point_vec.size(), consumer_id));
    }
    // 客户从度数小的先开始
    sort(consumer_tmp.begin(), consumer_tmp.end());
    if (ab_map["worktime_consumer_sort"] == 1) reverse(consumer_tmp.begin(), consumer_tmp.end());
    if (ab_map["worktime_consumer_sort"] == 2) random_shuffle(consumer_tmp.begin(), consumer_tmp.end());
    // 边缘从度数大的开始
    sort(tmp.begin(), tmp.end());
    reverse(tmp.begin(), tmp.end());
    if (ab_map["worktime_producer_sort"] == 1) reverse(tmp.begin(), tmp.end());
    if (ab_map["worktime_producer_sort"] == 2) random_shuffle(tmp.begin(), tmp.end());

    for (auto &cp : consumer_tmp) {
        int consumer_id = cp.second;
        // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        //     if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
        //         continue;
        //     }
        //     int cost_bandwidth = min(max(0, producers[producer_id].has_cost - 
        //                                     (producers[producer_id].bandwidth - producers[producer_id].time_remain_bandwidth)), 
        //                              consumers[consumer_id].time_need);
        //     if (cost_bandwidth == 0) continue;
        //     AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth);
        // }
        int block = max(1, consumers[consumer_id].time_need[time] / producer_number);
        block = 2000000;
        // 这里没必要135。。随便设置的，最差情况应该是135？
        // 先分配不花钱的（理论不花钱，实际上？
        if (time_index <= 10000) {
            AverageAdd(time, consumer_id, 1);
            // for (int k = 1; k <= 300; ++k) {
            //     if (consumers[consumer_id].time_need[time] == 0) {
            //         break ;
            //     }
            //     for(int i = tmp.size() - 1; i >= 0 ; --i) {
            //         auto& p = tmp[i];
            //         int producer_id = p.second;
            //         if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
            //             continue;
            //         }
            //         int cost_bandwidth = min(max(0, producers[producer_id].has_cost - producers[producer_id].time_cost[time]),
            //                                 consumers[consumer_id].time_need[time]);
            //         cost_bandwidth = min(cost_bandwidth, block);
            //         if (cost_bandwidth == 0) continue;
            //         AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth, 1);
            //     }
            // }
        } else {
            // 使用网络流分配试一试?
            WorkTimeMaxFlow(time, 1);
        }
        AverageAdd(time, consumer_id, 0);
        // for (int k = 1; k <= 300; ++k) {
        //     if (consumers[consumer_id].time_need[time] == 0) {
        //         break ;
        //     }
        //     for(int i = 0; i < tmp.size() ; ++i) {
        //         auto& p = tmp[i];
        //         int producer_id = p.second;
        //         if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
        //             continue;
        //         }
        //         int cost_bandwidth = min(max(0, producers[producer_id].bandwidth - producers[producer_id].time_cost[time]),
        //                                 consumers[consumer_id].time_need[time]);
        //         cost_bandwidth = min(cost_bandwidth, block);
        //         if (cost_bandwidth == 0) continue;
        //         AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth, 1);
        //     }
        // }
    }
}


void WorkTimeMinCostFlow(int time) {
    // MaxFlow mf;
    // mf.clear();
    // // run full use (not cost)
    // int source = 0, target = consumer_number + producer_number + 1;
    // for (int i = 1; i <= producer_number; ++i) {
    //     if (producers[i].is_full_use_time[time]) {
    //         mf.addEdge(source, i, producers[i].TimeRemain(time));
    //     }
    // }
    // for (int i = 1; i <= consumer_number; ++i) {
    //    mf.addEdge(producer_number + i, target, consumers[i].time_need[time]);
    // }
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //         if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
    //             continue;
    //         }
    //         mf.addEdge(producer_id, producer_number + consumer_id, mf.INF);
    //     }
    // }
    // mf.Dinic(target, source, target);
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     for (int i = mf.head[producer_number + consumer_id]; i != -1; i = mf.edges[i].nxt) {
    //         if (1 <= mf.edges[i].v && mf.edges[i].v <= producer_number && mf.edges[i].f > 0) {
    //             AddSomeBandWidth(time, mf.edges[i].v, consumer_id, mf.edges[i].f);
    //         }
    //     }
    // }

    // 代码直接用会有编译问题，部分变量名修改了，可以参考逻辑写费用流代码
    // MinCostFlow mcf;
    // mcf.init();
    // source = 0, target = consumer_number + producer_number + 1;
    // for (int i = 1; i <= producer_number; ++i) {
    //     if (producers[i].is_full_use_time[time]) {
    //         continue ;
    //     }
    //     mcf.addEdge(source, i, producers[i].time_remain_bandwidth - producers[i].has_cost, 1);
    //     mcf.addEdge(source, i, producers[i].has_cost, 0);
    // }
    // for (int i = 1; i <= consumer_number; ++i) {
    //    mcf.addEdge(producer_number + i, target, consumers[i].time_need, 0);
    // }
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //         if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
    //             continue;
    //         }
    //         mcf.addEdge(producer_id, producer_number + consumer_id, mcf.INF, 0);
    //     }
    // }
    // int cost=0, flow;
    // mcf.MCMF(source, target, cost, flow);
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     for (int i = mcf.head[producer_id]; i != -1; i = mcf.edge[i].nxt) {
    //         if (producer_number + 1 <= mcf.edge[i].to && mcf.edge[i].to <= producer_number + consumer_number && mcf.edge[i].flow > 0) {
    //             AddSomeBandWidth(time, producer_id, mcf.edge[i].to - producer_number, mcf.edge[i].flow);
    //         }
    //     }
    // }
    // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //     if (producers[producer_id].is_full_use_time[time]) {
    //         continue ;
    //     }
    //     producers[producer_id].has_cost = max(producers[producer_id].has_cost, 
    //                                           producers[producer_id].bandwidth - producers[producer_id].time_remain_bandwidth);
    // }
}
void EndWork();
void Reset();
void EndFlow(int nd_write);
void Work() {
    Reset();
    // 提前处理5%
    PreWork();
    // 尝试把每一天全能装满的全装满
    // TryFullUse();
    
    // vector<pair<int,int> > simplex_vec; // time_sum_cost, time
    // for (int time = 1; time <= times; ++time) {
    //     int sum_cost = 0;
    //     for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //         sum_cost += consumers[consumer_id].time_need[time];
    //     }
    //     simplex_vec.emplace_back(make_pair(sum_cost, time));
    // }
    // sort(simplex_vec.begin(), simplex_vec.end());
    // reverse(simplex_vec.begin(), simplex_vec.end());

    // while (simplex_vec.size() > 5) simplex_vec.pop_back();
    // vector<int> simplex_time_vec;
    // for (auto& p : simplex_vec) simplex_time_vec.emplace_back(p.second);
    // vector<int> simplex_has_cost;
    // simplex_has_cost = simplex.ProduceHasCost(simplex_time_vec);

    // // for (int i = 1; i <= producer_number; ++i) {
    // //     producers[i].has_cost = simplex_has_cost[i - 1];
    // // }

    // if (debug_file != nullptr && is_ab == 0) {
    //     for (int i = 1; i <= producer_number; ++i) {
    //         (*debug_file) << "simplex_has_cost, producer_id: " << i 
    //                       << " has_cost: " << simplex_has_cost[i - 1]
    //                       << endl;
    //     }
    // }

    // 求完全图下的比较优解 解法有问题
    // #define ll long long
    // vector<ll> vec; ll sum=0;
    // priority_queue<ll> cque;
    // priority_queue<pair<ll,int>> pque;
    // for(int i=1;i<=producer_number;++i) pque.push(make_pair(producers[i].bandwidth, can_full_use_time));
    // for(int i=1;i<=times;i++)
    // {
    //     ll tmp = 0;
    //     for(int j=1;j<=consumer_number;j++) tmp+=consumers[j].time_need[i];
    //     cque.push(tmp);
    // }
    // ll ans; pair<ll,int> p;
    // while(pque.size())
    // {
    //     p=pque.top(); pque.pop();
    //     ll tmp=cque.top(); cque.pop();
    //     tmp-=p.first;
    //     cque.push(tmp);
    //     p.second--;
    //     if(p.second) pque.push(make_pair(p.first,p.second));
    // }
    // ans=cque.top();
    // cout<<output_dir<<" "<<ans<<endl;
    
    // 处理5%后的理论最小耗费
    int max_sum_cost = 0;
    for (int i = 1; i <= times; ++i) {
        max_sum_cost = max(max_sum_cost, time_node[i].sum_cost);
    }
    if (result_file != nullptr && info[0] != '!' && is_ab == 0) {
        (*result_file) << "after 5%, max_sum_cost: " << max_sum_cost << endl;
    }
    int *p = new int[MAXT];
    if (debug_file != nullptr) {
        (*debug_file) << "after 5%, time_sum_cost: "; 
        int num = 0;
        for (int i = 1; i <= times; ++i) p[i] = time_node[i].sum_cost;
        sort(p + 1, p + 1 + times);
        for (int i = 1; i <= times; i += 100) {
            num++;
            (*debug_file) << p[i] << " ";
            if (num % 20 == 0) (*debug_file) << endl;
        }
        (*debug_file) << endl;
    }
    delete[] p;

    vector<pair<int,int> > time_cost_vec; // time_sum_cost, time
    for (int time = 1; time <= times; ++time) {
        int sum_cost = 0;
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            sum_cost += consumers[consumer_id].time_need[time];
        }
        time_cost_vec.emplace_back(make_pair(sum_cost, time));
    }
    int index = 0;
    sort(time_cost_vec.begin(), time_cost_vec.end());
    reverse(time_cost_vec.begin(), time_cost_vec.end());

    // ab一次太慢了
    while (is_ab == 1 && time_cost_vec.size() > 2000) time_cost_vec.pop_back();

    if (ab_map["work_time_sort"] == 1) reverse(time_cost_vec.begin(), time_cost_vec.end());
    if (ab_map["work_time_sort"] == 2) random_shuffle(time_cost_vec.begin(), time_cost_vec.end());
    for (auto& p : time_cost_vec) {
        int time = p.second;
        ++index;
        // 调用不同策略
        WorkTimeBaseline(time, index);
        // WorkTimeMaxFlow(time);
        
        vis_time[time] = 1;

        if ( index % 1000 == 0) EndWork();

        // 本地做一下检查
        if (debug_file != nullptr && is_ab == 0 || result_file != nullptr) {
            for (int i = 1; i <= producer_number; ++i) {
                // assert(producers[i].TimeRemain(time) >= 0);
                // producers[i].cost[time] = producers[i].bandwidth - producers[i].time_remain_bandwidth;
            }
            for (int i = 1; i <= consumer_number; ++i) {
                // assert(consumers[i].time_need[time] == 0);
            }
        }
    }
    
    memset(info_bandwidth, 0, sizeof(info_bandwidth));
    for (int k = 1; k <= 100; ++k) EndFlow(0);
    EndFlow(1);
    EndWork();


    if (debug_file != nullptr) {
        for (auto& p : ab_map) {
            (*debug_file) << "ab_params: "
                          << " key: " << p.first
                          << " value: " << p.second 
                          << endl;  
        }
        (*debug_file) << "ab_answer: " << GetAnswer() << endl;
    }
}

vector<pair<int, int> > has_cost_vec[MAXT]; // consumer, bandwidth
void EndWork() {
    for (int i = 1; i <= times; ++i) {
        has_cost_vec[i].clear();
    }
    // 重新计算has cost以后不去改变has cost，因此答案不会变劣
    int* tmp_cost = new int[MAXM]; //记录一下每个其他边缘暂时花了多少
    vector<pair<long long, int> > producer_vec; // waste, producer 
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        producers[producer_id].GetWaste(1);
        producer_vec.emplace_back(make_pair(producers[producer_id].waste, producer_id));
    }
    sort(producer_vec.begin(), producer_vec.end());
    reverse(producer_vec.begin(), producer_vec.end());
    for (int k = 1; k <= 10; ++k) {   
        int continue_flag = 0;
        for (int i = 0; i < producer_number; ++i) {
            int from_producer_id = producer_vec[i].second;
            for (int time = 1; time <= times; ++time) {
                has_cost_vec[time].clear();
                for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
                    if (info_bandwidth[time][from_producer_id][consumer_id]) {
                        has_cost_vec[time].emplace_back(make_pair(consumer_id, info_bandwidth[time][from_producer_id][consumer_id]));
                    }
                }
            }
            int min_can_has_cost = 0; // 能把耗费优化到多少，因为只考虑最大的，所以迁移一堆不一定合适，迁移到这个点就可以了
            for (int time = 1; time <= times; ++time) {
                if (producers[from_producer_id].is_full_use_time[time]) {
                    continue ;
                }
                for (int j = 1; j <= producer_number; ++j) {
                    tmp_cost[j] = 0;
                }
                int time_win = 0;
                for (auto& p : has_cost_vec[time]) {
                    int bandwidth = p.second;
                    for (auto &to_producer_id : consumers[p.first].can_visit_point_vec) {
                        // 从花费大的往小的迁移，理论更好，可以换排序实验，换了以后这里也要换
                        if (k == 1 && producers[to_producer_id].waste > producers[from_producer_id].waste && producers[to_producer_id].is_full_use_time[time] == 0) {
                            continue ;
                        }
                        if (to_producer_id == from_producer_id) {
                            continue ;
                        }
                        int tmp = max(0, producers[to_producer_id].has_cost - (producers[to_producer_id].time_cost[time] + tmp_cost[to_producer_id]));
                        // 5%, 不占花费
                        if (producers[to_producer_id].is_full_use_time[time]) {
                            tmp = max(0, producers[to_producer_id].bandwidth - (producers[to_producer_id].time_cost[time] + tmp_cost[to_producer_id]));
                        }
                        tmp = min(tmp, bandwidth);
                        bandwidth -= tmp;
                        time_win += tmp;
                        tmp_cost[to_producer_id] += tmp;
                    }
                }
                min_can_has_cost = max(min_can_has_cost, producers[from_producer_id].time_cost[time] - time_win);
            }
            for (int time = 1; time <= times; ++time) {
                if (producers[from_producer_id].is_full_use_time[time]) {
                    continue ;
                }
                int time_need_win = max(0, producers[from_producer_id].time_cost[time] - min_can_has_cost);
                for (auto& p : has_cost_vec[time]) {
                    if (time_need_win == 0) {
                        break ;
                    }
                    int bandwidth = p.second;
                    for (auto &to_producer_id : consumers[p.first].can_visit_point_vec) {
                        // 从花费大的往小的迁移，理论更好，可以换排序实验，换了以后这里也要换
                        if (k == 1 && producers[to_producer_id].waste > producers[from_producer_id].waste && producers[to_producer_id].is_full_use_time[time] == 0) {
                            continue ;
                        }
                        if (to_producer_id == from_producer_id) {
                            continue ;
                        }
                        int tmp = max(0, producers[to_producer_id].has_cost - producers[to_producer_id].time_cost[time]);
                        if (producers[to_producer_id].is_full_use_time[time]) {
                            tmp = max(0, producers[to_producer_id].bandwidth - producers[to_producer_id].time_cost[time]);
                        }
                        tmp = min(tmp, bandwidth);
                        tmp = min(tmp, time_need_win);
                        bandwidth -= tmp;
                        time_need_win -= tmp;
                        MoveSomeBandWidth(time, from_producer_id, to_producer_id, p.first, tmp);
                    }
                }
            }
            if (debug_file != nullptr && is_ab == 0) {
                (*debug_file) << "end_work, producer_id: " << from_producer_id << " can win: " 
                            << producers[from_producer_id].has_cost - min_can_has_cost << " bandwidth." <<endl;
            }
            if (producers[from_producer_id].has_cost > min_can_has_cost) {
                continue_flag = 1;
            }
            producers[from_producer_id].has_cost = min_can_has_cost;
        }
        if (continue_flag == 0) {
            break ;
        }
    }
    delete[] tmp_cost;
}

void EndFlow(int nd_write) {
    int* p = new int[MAXN];
    int* cp = new int[MAXN];
    int* c = new int[MAXN];
    int* cc = new int[MAXN];
    
    for (int i = 1; i <= producer_number; ++i) p[i] = i;
    random_shuffle(p + 1, p + 1 + producer_number);
    for (int i = 1; i <= consumer_number; ++i) c[i] = i;
    random_shuffle(c + 1, c + 1 + consumer_number);
    
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        producers[producer_id].CalFlowNeed();
    }
    int pre_answer = 0;
    for (int i = 1; i <= producer_number; ++i) pre_answer += producers[i].has_cost;
    // init
    int* pre_has_cost = new int[MAXN];
    for (int i = 1; i <= producer_number; ++i) {
        pre_has_cost[i] = producers[i].has_cost;
    }
    Reset();

    MaxFlow mf;
    mf.clear();

    int source = 0, target = consumer_number + producer_number + 1;
    for (int time = 1; time <= times; ++time) {
        mf.clear();
        for (int i = 1; i <= producer_number; ++i) {
            if (producers[p[i]].flow_is_full_use_time[time] == 1)
                mf.addEdge(source, i, producers[p[i]].bandwidth);
            // else
            //     mf.addEdge(source, i, pre_has_cost[p[i]]);
        }
        for (int i = 1; i <= consumer_number; ++i) {
            mf.addEdge(producer_number + i, target, consumers[c[i]].time_need[time]);
        }
        for (int i = 1; i <= producer_number; ++i) {
            for (int j = 1; j <= consumer_number; ++j) {
                if (producers[p[i]].can_visit_point[c[j]] == 0) {
                    continue;
                }
                mf.addEdge(i, producer_number + j, mf.INF);
            }
        }
        mf.Dinic(target, source, target);
        for (int i = 1; i <= producer_number; ++i) {
            if (producers[p[i]].flow_is_full_use_time[time] == 0)
                mf.addEdge(source, i, pre_has_cost[p[i]]);
        }
        mf.Dinic(target, source, target);
        for (int j = 1; j <= consumer_number; ++j) {
            for (int i = mf.head[producer_number + j]; i != -1; i = mf.edges[i].nxt) {
                if (1 <= mf.edges[i].v && mf.edges[i].v <= producer_number && mf.edges[i].f > 0) {
                    AddSomeBandWidth(time, p[mf.edges[i].v], c[j], mf.edges[i].f, 1, nd_write);
                }
            }
        }
    }    
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        producers[producer_id].CalFlowNeed();
    }
    int lst_answer = 0;
    for (int i = 1; i <= producer_number; ++i) lst_answer += producers[i].has_cost;
    
    if (debug_file != nullptr && is_ab == 0) {
        (*debug_file) << "end_flow: " 
                      << " win: " << pre_answer - lst_answer
                      << endl;
        (*debug_file) << "end_flow: has_cost: "; 
        for (int i = 1; i <= producer_number; ++i) (*debug_file) << producers[i].has_cost << " ";
        (*debug_file) << endl;
    }

    delete[] p;
    delete[] cp;
    delete[] c;
    delete[] cc;
    delete[] pre_has_cost;
}

void Reset() {
    for (int i = 1; i <= times; ++i) {
        time_node[i].Reset();
    }
    for (int i = 1; i <= producer_number; ++i) {
        producers[i].Reset();
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
    is_ab = 1;
    // ABWork();
    is_ab = 0;
    Work();
    Output();
    if (result_file != nullptr && info[0] != '!') {
        int sum_cost = 0;
        int check_flag = 1;
        for (int i = 1; i <= consumer_number; ++i) {
            for (int time = 1; time <= times; ++time) {
                if (consumers[i].time_need[time]) check_flag = 0;
            }
        }
        if (check_flag == 0) {
            (*result_file) << "check fail!" << endl;
        }
        for (int i = 1; i <= producer_number; ++i) {
            sort(producers[i].time_cost + 1, producers[i].time_cost + 1 + times);
            sum_cost += producers[i].time_cost[cost_max_index];
        }
        (*result_file) << "data_dir: " << data_dir <<endl;
        (*result_file) << "sum_cost: " << sum_cost << " cost_time: " << (double)(clock() - begin_time) / CLOCKS_PER_SEC << "s" << endl;
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
            (*debug_file) << "cost: " << producers[i].time_cost[cost_max_index] << endl;
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
    while (read_file >> data) {
        ++lines;
        Strip(data, '\r');
        vec_data = Split(data, ",");    
        if (lines == 1) {
            for (int i = 1; i < vec_data.size(); ++i) {
                consumers_name_id_map[vec_data[i]] = i;
                consumers[i].name = vec_data[i];
            }
            consumer_number = vec_data.size() - 1;
            continue;
        }
        for (int i = 1; i < vec_data.size(); ++i) {
            consumers[i].time_need[lines - 1] = StringToInt(vec_data[i]);
            if (debug_file != nullptr && i <= 3 && lines - 1 <= 3) {
                (*debug_file) << "consumers: " << i << " day: " << lines - 1 << " name: " << consumers[i].name 
                << " bandwidth: " << consumers[i].time_need[lines - 1] << endl;
            }
        }
    }
    times = lines - 1;
    
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
        qos_down = StringToInt(vec_data[1]);
        if (debug_file != nullptr) {
            (*debug_file) << "qos_down: " << qos_down << endl;
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
    vector<pair<int,int> > consumer_get_info[MAXM]; // producer_id, band_width
    for (int time = 1; time <= times; ++time) {
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
                if (info_bandwidth[time][producer_id][consumer_id] <= 0) {
                    continue ;
                }
                consumer_get_info[consumer_id].emplace_back(
                    make_pair(producer_id, info_bandwidth[time][producer_id][consumer_id])
                );
            }
        }
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {   
            int begin_flag = 1;
            file<<consumers[consumer_id].name<<":";
            for(auto &p : consumer_get_info[consumer_id]) {
                if (begin_flag) {
                    begin_flag = 0;
                } else {
                    file << ',';
                }
                file<<'<';
                file<<producers[p.first].name;
                file<<',';
                file<<p.second;
                file<<'>';
            }
            if (time != times || consumer_id != consumer_number) file<<'\n';
        }
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            consumer_get_info[consumer_id].clear();
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
