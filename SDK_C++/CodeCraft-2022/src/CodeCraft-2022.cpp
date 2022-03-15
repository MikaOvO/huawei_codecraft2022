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

using namespace std;

ofstream *debug_file = nullptr;
ofstream *result_file = nullptr;

string output_dir = "/output";
string data_dir = "/data";

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
    int time_remain_bandwidth;
    int can_visit_point[MAXM];
    int is_full_use_time[MAXT];
    int cost[MAXT];
    priority_queue<int, vector<int>, greater<int> > cost_que;
    Producer() {
        for (int i = 0; i < MAXT; ++i) is_full_use_time[i] = 0;
        has_cost = 0;
        has_use_full_time = 0;
    }
    bool operator < (const Producer& producer) const {
        return bandwidth < producer.bandwidth;
    }
} producers[MAXN];

struct Consumer {
    string name;
    int time_need_bandwidth;
    int can_visit_point[MAXN];
    int need_bandwidth[MAXT];
    Consumer() {
    }
} consumers[MAXM];

// On (producer_id, consumer_id)-info Give value-bandwidth
map<pair<int,int>, int> info_bandwidth[MAXT]; 

void Init() {
    srand(11);
    ::cost_max_index = (times * 95 + 99) / 100;
    if (debug_file != nullptr) {
        (*debug_file) << "cost_max_index: " << cost_max_index <<endl;
    }
    ::can_full_use_time = times - cost_max_index;
    int *p = new int[MAXT];
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        // for (int i = 1; i <= cost_max_index; ++i) producers[producer_id].cost_que.push(0);
        for (int i = 1; i <= times; ++i) p[i] = i;
        random_shuffle(p + 1, p + 1 + times);
        for (int i = 1; i <= can_full_use_time; ++i) {
            producers[producer_id].is_full_use_time[p[i]] = 1;
            if (debug_file != nullptr) {
                (*debug_file) << "producer: " << producer_id << " on time: " << p[i] << " full use." << endl;
            }
        }
    }
    delete[] p;
}

void AddSomeBandWidth(int time, int producer_id, int consumer_id, int bandwidth) {
    // 把耗费放到当天里，注意策略里不需要写这一块了
    producers[producer_id].time_remain_bandwidth -= bandwidth;
    consumers[consumer_id].time_need_bandwidth -= bandwidth;

    if (info_bandwidth[time].find(make_pair(producer_id, consumer_id)) == info_bandwidth[time].end()) 
        info_bandwidth[time][make_pair(producer_id, consumer_id)] = 0;
    info_bandwidth[time][make_pair(producer_id, consumer_id)] += bandwidth;
}

void WorkTimeBaseline(int time) {
    for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
        // for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        //     if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
        //         continue;
        //     }
        //     int cost_bandwidth = min(max(0, producers[producer_id].has_cost - 
        //                                     (producers[producer_id].bandwidth - producers[producer_id].time_remain_bandwidth)), 
        //                              consumers[consumer_id].time_need_bandwidth);
        //     if (cost_bandwidth == 0) continue;
        //     AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth);
        // }
        for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
            if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
                continue;
            }
            int cost_bandwidth = min(producers[producer_id].time_remain_bandwidth, consumers[consumer_id].time_need_bandwidth);
            if (cost_bandwidth == 0) continue;
            AddSomeBandWidth(time, producer_id, consumer_id, cost_bandwidth);
        }
    }
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        producers[producer_id].has_cost = max(producers[producer_id].has_cost, 
                                              producers[producer_id].bandwidth - producers[producer_id].time_remain_bandwidth);
    }
}

void WorkTimeMaxFlow(int time) {
    MaxFlow mf;
    mf.clear();
    int source = 0, target = consumer_number + producer_number + 1;
    // run full use (not cost)
    for (int i = 1; i <= producer_number; ++i) {
        if (producers[i].is_full_use_time[time]) {
            mf.addEdge(source, i, producers[i].time_remain_bandwidth);
        }
    }
    for (int i = 1; i <= consumer_number; ++i) {
       mf.addEdge(producer_number + i, target, consumers[i].time_need_bandwidth);
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

    // run other
    mf.clear();
    for (int i = 1; i <= producer_number; ++i) {
        mf.addEdge(source, i, producers[i].time_remain_bandwidth);
    }
    for (int i = 1; i <= consumer_number; ++i) {
       mf.addEdge(producer_number + i, target, consumers[i].time_need_bandwidth);
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

void WorkTimeMinCostFlow(int time) {
    MaxFlow mf;
    mf.clear();
    // run full use (not cost)
    int source = 0, target = consumer_number + producer_number + 1;
    for (int i = 1; i <= producer_number; ++i) {
        if (producers[i].is_full_use_time[time]) {
            mf.addEdge(source, i, producers[i].time_remain_bandwidth);
        }
    }
    for (int i = 1; i <= consumer_number; ++i) {
       mf.addEdge(producer_number + i, target, consumers[i].time_need_bandwidth);
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

    MinCostFlow mcf;
    mcf.init();
    source = 0, target = consumer_number + producer_number + 1;
    for (int i = 1; i <= producer_number; ++i) {
        if (producers[i].is_full_use_time[time]) {
            continue ;
        }
        mcf.addEdge(source, i, producers[i].time_remain_bandwidth - producers[i].has_cost, 1);
        mcf.addEdge(source, i, producers[i].has_cost, 0);
    }
    for (int i = 1; i <= consumer_number; ++i) {
       mcf.addEdge(producer_number + i, target, consumers[i].time_need_bandwidth, 0);
    }
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
            if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
                continue;
            }
            mcf.addEdge(producer_id, producer_number + consumer_id, mcf.INF, 0);
        }
    }
    int cost=0, flow;
    mcf.MCMF(source, target, cost, flow);
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        for (int i = mcf.head[producer_id]; i != -1; i = mcf.edge[i].nxt) {
            if (producer_number + 1 <= mcf.edge[i].to && mcf.edge[i].to <= producer_number + consumer_number && mcf.edge[i].flow > 0) {
                AddSomeBandWidth(time, producer_id, mcf.edge[i].to - producer_number, mcf.edge[i].flow);
            }
        }
    }
    for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
        if (producers[producer_id].is_full_use_time[time]) {
            continue ;
        }
        producers[producer_id].has_cost = max(producers[producer_id].has_cost, 
                                              producers[producer_id].bandwidth - producers[producer_id].time_remain_bandwidth);
    }
}

void Work() {
    for (int time = 1; time <= times; ++time) {
        for (int i = 1; i <= producer_number; ++i) {
            producers[i].time_remain_bandwidth = producers[i].bandwidth;
        }
        for (int i = 1; i <= consumer_number; ++i) {
            consumers[i].time_need_bandwidth = consumers[i].need_bandwidth[time];
        }
        // 调用不同策略
        WorkTimeMaxFlow(time);
        // 本地做一下检查
        if (debug_file != nullptr || result_file != nullptr) {
            for (int i = 1; i <= producer_number; ++i) {
                assert(producers[i].time_remain_bandwidth >= 0);
                producers[i].cost[time] = producers[i].bandwidth - producers[i].time_remain_bandwidth;
            }
            for (int i = 1; i <= consumer_number; ++i) {
                assert(consumers[i].time_need_bandwidth == 0);
            }
        }
    }
}

int main(int argc, char *argv[]) {
    begin_time = clock();
    // local or oj
    for (int i = 1; i < argc; ++i) {
        int len = strlen(argv[i]);
        string tmp = "";
        for (int j = 0; j < len; ++j) tmp += argv[i][j];
        Strip(tmp, '\r');
        vector<string> str_vec = Split(tmp, "=");
        if (str_vec[0] == "debug_file") {
            debug_file = new ofstream();
            debug_file->open(str_vec[1], ios::out);
        } 
        if (str_vec[0] == "result_file") {
            result_file = new ofstream();
            result_file->open(str_vec[1], ios::out);
        } 
        if (str_vec[0] == "output_dir") {
            output_dir = str_vec[1];
        } 
        if (str_vec[0] == "data_dir") {
            data_dir = str_vec[1];
        }
    }
    ReadIn();
    Init();
    Work();
    Output();
    if (result_file != nullptr) {
        int sum_cost = 0;
        for (int i = 1; i <= producer_number; ++i) {
            sort(producers[i].cost + 1, producers[i].cost + 1 + times);
            sum_cost += producers[i].cost[cost_max_index];
        }
        (*result_file) << "sum_cost: " << sum_cost <<endl;
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
            consumers[i].need_bandwidth[lines - 1] = StringToInt(vec_data[i]);
            if (debug_file != nullptr) {
                (*debug_file) << "consumers: " << i << " day: " << lines - 1 << " name: " << consumers[i].name 
                << " bandwidth: " << consumers[i].need_bandwidth[lines - 1] << endl;
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
        if (debug_file != nullptr) {
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
                if (debug_file != nullptr) {
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
        for (auto &pair : info_bandwidth[time]) {
            consumer_get_info[pair.first.second].emplace_back(
                make_pair(pair.first.first, pair.second)
            );
            if (debug_file != nullptr) {
                consumers[pair.first.second].need_bandwidth[time] -= pair.second;
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
