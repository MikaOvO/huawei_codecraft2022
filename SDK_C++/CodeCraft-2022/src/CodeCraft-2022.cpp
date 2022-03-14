#include <iostream>
#include <map>
#include <string>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>

using namespace std;

ofstream *debug_file = nullptr;
ofstream *result_file = nullptr;

string output_dir = "/output";
string data_dir = "/data";

const int MAXT = 8928 + 1;
const int MAXM = 35 + 1;
const int MAXN = 135 + 1; 

clock_t begin_time;

// utils
int StringToInt(string &str);
void Strip(string& str, const char& chr);
vector<string> Split(const string& str, const string& delim);

// IO
void ReadIn();
void Output();

map<string, int> producers_name_id_map;
map<string, int> consumers_name_id_map;

int times, producer_number, consumer_number;

struct Producer {
    string name;
    int bandwidth;
    int time_remain_bandwidth;
    int can_visit_point[MAXM];
    Producer() {}
} producers[MAXN];

struct Consumer {
    string name;
    int can_visit_point[MAXN];
    int time_need_bandwidth;
    int need_bandwidth[MAXT];
    Consumer() {
    }
} consumers[MAXM];

// On (producer_id, consumer_id)-info Give value-bandwidth
map<pair<int,int>, int> info_bandwidth[MAXT]; 

void WorkTime(int time) {
    //  QAQ
    // for (int consumer_id = 1; consumer_id <= consumer_number; ++consumer_id) {
    //     for (int producer_id = 1; producer_id <= producer_number; ++producer_id) {
    //         if (consumers[consumer_id].can_visit_point[producer_id] == 0) {
    //             continue;
    //         }
    //         int cost_bandwidth = min(producers[producer_id].time_remain_bandwidth, consumers[consumer_id].time_need_bandwidth);
    //         if (cost_bandwidth == 0) continue;
    //         producers[producer_id].time_remain_bandwidth -= cost_bandwidth;
    //         consumers[consumer_id].time_need_bandwidth -= cost_bandwidth;
    //         if (info_bandwidth[time].find(make_pair(producer_id, consumer_id)) == info_bandwidth[time].end()) 
    //             info_bandwidth[time][make_pair(producer_id, consumer_id)] = 0;
    //         info_bandwidth[time][make_pair(producer_id, consumer_id)] += cost_bandwidth;
    //     }
    // }
}

void Work() {
    for (int time = 1; time <= times; ++time) {
        for (int i = 1; i <= producer_number; ++i) {
            producers[i].time_remain_bandwidth = producers[i].bandwidth;
        }
        for (int i = 1; i <= consumer_number; ++i) {
            consumers[i].time_need_bandwidth = consumers[i].need_bandwidth[time];
        }
        WorkTime(time);
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
    Work();
    Output();
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
            if (StringToInt(vec_data[i]) >= qos_down) {
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
    if (debug_file != nullptr) {
        for(int time = 1; time <= times; ++time)
            for (int i = 1; i <= consumer_number; ++i) {
                if (consumers[i].need_bandwidth[time]) {
                    (*debug_file) << "check fail !!!" <<endl;
                    return ;
                }
            }
    }
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
