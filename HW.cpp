#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>


std::ostream &qDebug();
int qRound(double);

int sin_table[371];
int cos_table[371];
const double PI = 3.141592653589;
const double scale_factor = PI / 360.0;
const double scale_factor_2 = 360.0 / PI;
int exp_times = 0;
int global_file_num = 1;

#define user_path "/home/mini/Desktop/Data/%1"
#define exp_times_path "" // QString(user_path).arg("exp_times.txt").toStdString().c_str()
#define file(file_name)                                                                  \
    "" // QString(user_path).arg(QString("%1_%2_%3.csv").arg(exp_times).arg(file_name).arg(file_num)).toStdString().c_str()

struct init_helper
{
    init_helper()
    {
        std::ifstream ifs(exp_times_path);
        if (ifs.good())
            ifs >> exp_times;
        else
            exp_times = 0;

        exp_times++;

        std::ofstream ofs(exp_times_path);
        ofs << exp_times;

        for (std::size_t i = 0; i <= 360; i++)
        {
            sin_table[i] = std::sin(i * scale_factor) * 256.0;
            cos_table[i] = std::cos(i * scale_factor) * 256.0;
        }
    }
} helper;

struct imu_data_type
{
    int timestamp;
    double orientation;
    double odometry;
    imu_data_type(double _ori, double _odometry, int _timestamp)
    {
        orientation = _ori;
        odometry = _odometry;
        timestamp = _timestamp;
    }
};

struct point_type
{

    int x;
    int y;
    point_type() {}
    point_type(int _x, int _y)
    {
        x = _x;
        y = _y;
    }
};

struct line_type
{

    int theta;

    int r;
    int weight;
    line_type() {}
    line_type(int _theta, int _r, int _weight)
    {
        theta = _theta;
        r = _r;
        weight = _weight;
    }
};

int fast_sin(int angle) { return sin_table[angle]; }
int fast_cos(int angle) { return cos_table[angle]; }

int point_distance_sqr(int p1_x, int p1_y, int p2_x, int p2_y)
{
    return std::pow(p1_x - p2_x, 2) + std::pow(p1_y - p2_y, 2);
}

std::vector<line_type> line_cluster(const std::vector<line_type> &lines)
{
    std::size_t n_lines = lines.size();
    std::vector<unsigned char> bool_vec(n_lines, true);
    std::vector<std::vector<line_type> > groups;
    for (std::size_t i = 0; i < n_lines; i++)
    {
        line_type line_i = lines[i];
        if (bool_vec[i])
        {
            groups.push_back(std::vector<line_type>(1, line_i));
            bool_vec[i] = false;
            for (std::size_t j = i + 1; j < n_lines; j++)
            {
                line_type line_j = lines[j];
                if (bool_vec[j] && point_distance_sqr(line_i.theta, line_i.r,
                                                      line_j.theta, line_j.r) <= 900)
                {
                    bool_vec[j] = false;
                    groups.back().push_back(line_j);
                }
            }
        }
    }

    std::size_t n_groups = groups.size();

    std::vector<line_type> possible_return_value;
    for (std::size_t i = 0; i < n_groups; i++)
    {
        std::size_t group_size = groups[i].size();
        int avg_theta = 0, avg_r = 0, sum_weight = 0;
        for (std::size_t j = 0; j < group_size; j++)
        {
            avg_theta += groups[i][j].theta * groups[i][j].weight;
            avg_r += groups[i][j].r * groups[i][j].weight;
            sum_weight += groups[i][j].weight;
        }
        avg_theta /= sum_weight;
        avg_r /= sum_weight;
        possible_return_value.push_back(line_type(avg_theta, avg_r, sum_weight));
    }

    return possible_return_value;
}

int hough_matrix[370][411];
point_type laser_points[371];

void do_hough(int x, int y)
{
    for (std::size_t alpha = 0; alpha < 360; alpha++)
    {
        int r = (x * fast_sin(alpha) - y * fast_cos(alpha)) >> 16;
        if (r >= -200 && r <= 200)
            hough_matrix[alpha][r + 200]++;
    }
}

int num_distance(int a, int b) { return std::abs(a - b); }

void get_laser_points(const short *lidar_data, std::size_t start_angle,
                      std::size_t end_angle)
{
    for (std::size_t i = start_angle; i <= end_angle; i++)
    {
        laser_points[i] =
            point_type(lidar_data[i] * fast_cos(i), lidar_data[i] * fast_sin(i));
    }
}

std::vector<line_type> find_lines(const short *lidar_data, std::size_t start_angle = 0,
                                  std::size_t end_angle = 360)
{

    std::memset(hough_matrix, 0, sizeof(hough_matrix));

    get_laser_points(lidar_data, start_angle, end_angle);

    for (std::size_t i = start_angle; i <= end_angle; i++)
    {
        if (laser_points[i].x != 0 && laser_points[i].y != 0)
            do_hough(laser_points[i].x, laser_points[i].y);
    }

    std::vector<line_type> lines;
    for (std::size_t alpha = 0; alpha < 360; alpha++)
    {
        for (std::size_t r = 0; r < 401; r++)
        {
            int weight = hough_matrix[alpha][r];
            if (weight > 35)
                lines.push_back(line_type(alpha, r, weight));
        }
    }

    lines = line_cluster(lines);
    lines = line_cluster(lines);

    std::size_t n_lines = lines.size();
    for (std::size_t i = 0; i < n_lines; i++)
        lines[i].r -= 200;

    return lines;
}

std::vector<line_type> last_believable_result;

double zero_orientation;

double process_ori(double ori)
{
    while (ori < 0)
        ori += 2 * PI;
    while (ori >= 2 * PI)
        ori -= 2 * PI;
    return ori;
}

template <typename T> void save_value(const char *file_name_prefix, int file_num, T value)
{
    std::ofstream fs(file(file_name_prefix));
    fs << value;
}

void save_lidar_data(const short *lidar_data, const char *file_name_prefix, int file_num)
{
    std::ofstream ofs(file(file_name_prefix));
    for (std::size_t i = 0; i < 361; i++)
        ofs << lidar_data[i] << '\n';
}

void save_imu_data(imu_data_type imu_data, const char *file_name_prefix, int file_num)
{
    std::ofstream ofs(file(file_name_prefix));
    ofs << imu_data.timestamp << ',' << imu_data.orientation << ',' << imu_data.odometry;
}

void save_lines(const std::vector<line_type> &lines, const char *file_name_prefix,
                int file_num)
{
    std::ofstream ofs(file(file_name_prefix));
    std::size_t n = lines.size();

    for (std::size_t i = 0; i < n; i++)
        ofs << lines[i].theta << ',' << lines[i].r << '\n';
}

std::vector<line_type> process_1_line(const std::vector<line_type> &lines)
{
    std::vector<line_type> ret_lines = lines;
    int wall_dis = num_distance(last_believable_result[0].r, last_believable_result[1].r);
    if (ret_lines[0].r < 0)
        ret_lines.push_back(line_type(ret_lines[0].theta, ret_lines[0].r + wall_dis, 0));
    else
        ret_lines.push_back(line_type(ret_lines[0].theta, ret_lines[0].r - wall_dis, 0));
    return ret_lines;
}

std::vector<line_type> process_lines(const std::vector<line_type> &lines,
                                     const short *lidar_data, imu_data_type imu_data)
{
    int current_heading =
        qRound(process_ori(zero_orientation - imu_data.orientation) * scale_factor_2);

    const int line_heading = 360 - current_heading;

    std::vector<line_type> ret_lines;
    std::size_t n_lines = lines.size();
    for (std::size_t i = 0; i < n_lines; i++)
    {
        if (num_distance(lines[i].theta, line_heading) < 30)
            ret_lines.push_back(lines[i]);
    }

    std::size_t n_ret_lines = ret_lines.size();

    if (n_ret_lines == 2 && std::max(ret_lines[0].r, ret_lines[1].r) > 0 &&
        std::min(ret_lines[0].r, ret_lines[1].r) < 0)
    {
        if (num_distance(ret_lines[0].r, ret_lines[1].r) < 350)
        {
            last_believable_result = ret_lines;
            qDebug() << "[INFO] `process_lines`: 2 lines found. (" << ret_lines[0].theta
                     << ',' << ret_lines[0].r << ") (" << ret_lines[1].theta << ','
                     << ret_lines[1].r << ')';
            return ret_lines;
        }
        else
        {
            if (std::abs(ret_lines[0].r) > std::abs(ret_lines[1].r))
                ret_lines.pop_back();
            else
                ret_lines.erase(ret_lines.begin());

            ret_lines = process_1_line(ret_lines);
            qDebug() << "[WARN] `process_lines`: 2 line found, but they have something "
                        "wrong. Try to fix it. ("
                     << ret_lines[0].theta << ',' << ret_lines[0].r << ") ("
                     << ret_lines[1].theta << ',' << ret_lines[1].r << ')';
            return ret_lines;
        }
    }

    if (n_ret_lines == 1)
    {
        ret_lines = process_1_line(ret_lines);
        qDebug() << "[WARN] `process_lines`: 1 line found. Try to fix it. ("
                 << ret_lines[0].theta << ',' << ret_lines[0].r << ") ("
                 << ret_lines[1].theta << ',' << ret_lines[1].r << ')';
        last_believable_result = ret_lines;
        return ret_lines;
    }

    if (n_ret_lines == 0)
    {
        qDebug() << "[ERRO] `process_lines`: 0 line found. Data stored. Cannot fix it.";
        save_lidar_data(lidar_data, "Err points", imu_data.timestamp);
        ret_lines.clear();
        return ret_lines;
    }

    save_lidar_data(lidar_data, "Err points", imu_data.timestamp);
    save_lines(ret_lines, "Err lines", imu_data.timestamp);

    std::vector<line_type> left_lines;
    std::vector<line_type> right_lines;

    for (std::size_t i = 0; i < n_ret_lines; i++)
    {
        if (ret_lines[i].r < 0)
            left_lines.push_back(ret_lines[i]);
        else
            right_lines.push_back(ret_lines[i]);
    }

    ret_lines.clear();

    std::size_t n_left_lines = left_lines.size();
    if (n_left_lines > 0)
    {
        int avg_left_theta = 0, avg_left_r = 0;
        for (std::size_t i = 0; i < n_left_lines; i++)
        {
            avg_left_theta += left_lines[i].theta;
            avg_left_r += left_lines[i].r;
        }
        avg_left_theta /= n_left_lines;
        avg_left_r /= n_left_lines;
        ret_lines.push_back(line_type(avg_left_theta, avg_left_r, 0));
    }
    std::size_t n_right_lines = right_lines.size();
    if (n_right_lines > 0)
    {
        int avg_right_theta = 0, avg_right_r = 0;
        for (std::size_t i = 0; i < n_right_lines; i++)
        {
            avg_right_theta += right_lines[i].theta;
            avg_right_r += right_lines[i].r;
        }
        avg_right_theta /= n_right_lines;
        avg_right_r /= n_right_lines;
        ret_lines.push_back(line_type(avg_right_theta, avg_right_r, 0));
    }

    if (ret_lines.size() == 1)
        ret_lines = process_1_line(ret_lines);
    qDebug() << "[ERRO] `process_lines` found wrong lines. Data stored. Try to fix it. ("
             << ret_lines[0].theta << ',' << ret_lines[0].r << ") (" << ret_lines[1].theta
             << ',' << ret_lines[1].r << ')';
    return ret_lines;
}

bool need_init = true;

bool init(const short *lidar_data, imu_data_type imu_data)
{
    if (!need_init)
        return true;

    std::vector<line_type> lines = find_lines(lidar_data);
    if (lines.size() == 2 && num_distance(lines[0].theta, lines[1].theta) < 8 &&
        std::max(lines[0].r, lines[1].r) > 0 && std::min(lines[0].r, lines[1].r) < 0)
    {
        double init_orientation = process_ori(imu_data.orientation);
        double init_heading = 360.0 - (lines[0].theta + lines[1].theta) * 0.5;
        zero_orientation = process_ori(init_orientation + init_heading * scale_factor);

        last_believable_result = lines;
        qDebug() << "[INFO] `init`: 2 lines found. (" << lines[0].theta << ','
                 << lines[0].r << ") (" << lines[1].theta << ',' << lines[1].r << ')';
        need_init = false;
        return true;
    }

    qDebug() << "[ERRO] `init`: Cannot find valid lines.";
    return false;
}

struct output_type
{
    int steer, speed;
    output_type(int _steer, int _speed)
    {
        steer = _steer;
        speed = _speed;
    }
};

template <typename T> T clamp(T x, T down, T up)
{
    if (x > up)
        return up;
    if (x < down)
        return down;
    return x;
}

double process_error(double err)
{
    return 90.0 - 90.0 * std::exp(-std::abs(err) * 0.025);
}

const double PID_P = 6.0;

const int straight_speed = 160;
const int turn_speed = 100;

output_type go_straight(const short *lidar_data, imu_data_type imu_data,
                        bool use_dis_to_left_arg = false, int dis_to_left_arg = 0,
                        bool use_dis_to_right_arg = false, int dis_to_right_arg = 0,
                        bool use_error = false,double Error = 0,
                        bool use_speed = false,int input_speed = 0)
{
    static output_type last_output(0, 0);
    const output_type default_output(0, 150);
    if (!init(lidar_data, imu_data))
    {
        qDebug() << "[ERRO] `go_straight`: `init` returns `false`. Output ( 0 , 100 ) by "
                    "default.";
        zero_orientation = process_ori(imu_data.orientation + PI * 0.5);
        return default_output;
    }
    std::vector<line_type> lines =
        process_lines(find_lines(lidar_data), lidar_data, imu_data);

    save_lines(lines, "Lines_Data", global_file_num);
    if (lines.size() == 0)
    {
        qDebug() << "[ERRO] `go_straight`: `process_lines` returns 0 lines. Output ("
                 << (last_output.steer / 1) << ',' << last_output.speed
                 << ") by default.";
        output_type out(last_output.steer / 2, 100);
        last_output = out;
        return out;
    }

    double current_heading = 360.0 - (lines[0].theta + lines[1].theta) * 0.5;
    zero_orientation = process_ori(imu_data.orientation + current_heading * scale_factor);

    int left_r = std::min(lines[0].r, lines[1].r);
    int right_r = std::max(lines[0].r, lines[1].r);
    double pid_p = PID_P;
    double error;
    if (use_dis_to_left_arg)
    {
        error = -left_r - dis_to_left_arg;
        pid_p = 2.0 * PID_P;
    }
    else if (use_dis_to_right_arg)
    {
        error = dis_to_right_arg - right_r;
        pid_p = 2.0 * PID_P;
    }
    else
        error = -(left_r + right_r) * 0.5;

    double expect_heading = 180.0;

    if(use_error){ // 直接使用已经算好的error
        error = Error;
    }

    double correction_item = process_error(error);
    if (error < 0)
        expect_heading -= correction_item;
    else
        expect_heading += correction_item;

    qDebug() << "[INFO] `go_straight`: Current heading and error are (" << current_heading
             << ',' << error << ')';
    qDebug() << "[INFO] `go_straight`: Expect heading is" << expect_heading;

    int output_steer = qRound(pid_p * (current_heading - expect_heading));

    output_steer = clamp(output_steer, -400, 400);

    qDebug() << "[INFO] `go_straight`: Output (" << output_steer << ',' << straight_speed
             << ')';
    
    int output_speed = 0;
    if(use_speed) output_speed = input_speed;
    else output_speed = straight_speed;

    output_type out(output_steer, output_speed);
    last_output = out;
    return out;
}

enum turn
{
    turn_left,
    turn_right,
    turn_none
};

bool turn_right_flag = false;
bool turn_left_flag = false;
bool has_turned_left = false;
int turn_right_cont = 0;

void judge_turn(const short *lidar_data)
{
    bool front_barrier = false;
    int counter = 0;

    for (std::size_t i = 150; i <= 210; i++)
    {
        if (lidar_data[i] > 0 && lidar_data[i] < 220)
            counter++;
    }

    if (counter > 50)
        front_barrier = true;

    int threshold = 600;

    if (front_barrier)
        threshold = 300;

    counter = 0;
    for (std::size_t i = 0; i < 30; i++)
    {
        if (lidar_data[i] == 0 || lidar_data[i] > threshold)
            counter++;
    }

    if (counter > 20)
    {
        turn_right_flag = true;
        turn_right_cont++;
        qDebug() << "[INFO] `judge_turn`: Start turning right.";
        return;
    }

    counter = 0;
    for (std::size_t i = 330; i <= 360; i++)
    {
        if (lidar_data[i] == 0 || lidar_data[i] > threshold)
            counter++;
    }

    if (counter > 20 && !has_turned_left)
    {
        turn_left_flag = true;
        has_turned_left = true;
        qDebug() << "[INFO] `judge_turn`: Start turning left.";
        return;
    }
}

output_type turn_right_func(imu_data_type imu_data)
{
    double difference = process_ori(zero_orientation - imu_data.orientation);

    if (difference > PI)
        difference -= 2.0 * PI;

    int output_steer = 400;
    qDebug() << "[INFO] `turn_right_func`: Output (" << output_steer << ',' << turn_speed
             << ')';
    if (std::abs(difference * scale_factor_2) < 10)
    {
        turn_right_flag = false;
        need_init = true;
        qDebug() << "[INFO] `turn_right_func`: End turning right.";
        output_steer = 0;
    }
    return output_type(output_steer, turn_speed);
}

output_type turn_left_func(imu_data_type imu_data)
{
    double difference = process_ori(imu_data.orientation - zero_orientation + PI);

    if (difference > PI)
        difference -= 2.0 * PI;

    int output_steer = -400;
    qDebug() << "[INFO] `turn_left_func`: Output (" << output_steer << ',' << turn_speed
             << ')';
    if (std::abs(difference * scale_factor_2) < 10)
    {
        turn_left_flag = false;
        need_init = true;
        qDebug() << "[INFO] `turn_left_func`: End turning left.";
        output_steer = 0;
    }
    return output_type(output_steer, turn_speed);
}

bool circled = false;
bool circle_flag = false;

void judge_circle(const short *lidar_data, imu_data_type imu_data)
{
    if (turn_right_cont != 4 || circled)
        return;
    int barrier_start = 0, barrier_end = 0;
    int head_laser =
        360 - process_ori(zero_orientation - imu_data.orientation) * scale_factor_2;

    for (std::size_t i = head_laser - 15; i < head_laser + 15; i++)
    {
        if (i < 360 && i >= 0)
        {
            if (barrier_start == 0 && (lidar_data[i] == 0 || lidar_data[i] > 500) &&
                lidar_data[i + 1] < 200)
                barrier_start = i + 1;
            if (barrier_end == 0 && (lidar_data[i + 1] == 0 || lidar_data[i + 1] > 500) &&
                lidar_data[i] < 200)
                barrier_end = i;
        }
    }

    int barrier_width_sqr = std::pow(lidar_data[barrier_start], 2) +
                            std::pow(lidar_data[barrier_end], 2) -
                            2 * lidar_data[barrier_start] * lidar_data[barrier_end] *
                                std::cos(barrier_end - barrier_start);
    if (barrier_start != 0 && barrier_end != 0 && barrier_width_sqr > 196)
    {
        circle_flag = true;
        qDebug() << "[INFO] `judge_circle`: Circle start!";
    }
}

enum circle_state
{
    circle_left,
    circle_right
};

output_type circle(const short *lidar_data, imu_data_type imu_data)
{
    static circle_state state;
    static int circle_count = 0;

    static int last_distance_to_right = 0, last_distance_to_left = 0;
    static bool is_init = false;

    if (!is_init)
    {
        is_init = true;
        state = circle_right;
        last_distance_to_right = lidar_data[0];
        last_distance_to_left = lidar_data[360];
        qDebug() << "[INFO] `circle`: Switch to right.";
    }

    switch (state)
    {
    case circle_left:
        if (lidar_data[0] - last_distance_to_right > 50)
        {
            state = circle_right;
            qDebug() << "[INFO] `circle`: Switch to right.";
            circle_count++;
        }
        break;
    case circle_right:
        if (lidar_data[360] - last_distance_to_left > 50)
        {
            state = circle_left;
            qDebug() << "[INFO] `circle`: Switch to left.";
            circle_count++;
        }
        break;
    }
    last_distance_to_right = lidar_data[0];
    last_distance_to_left = lidar_data[360];
    if (circle_count == 4)
    {
        circled = true;
        circle_flag = false;
        qDebug() << "[INFO] `circle`: Circle end.";
        return go_straight(lidar_data, imu_data);
    }

    output_type out(0, 0);
    switch (state)
    {
    case circle_left:
        out = go_straight(lidar_data, imu_data, true, 40);
        break;
    case circle_right:
        out = go_straight(lidar_data, imu_data, false, 0, true, 40);
        break;
    }

    qDebug() << "[INFO] `circle`: Output (" << out.steer << ',' << out.speed << ')';
    return out;
}

output_type selector(const short *lidar_data, imu_data_type imu_data)
{
    static bool is_init = false;
    static int start_stamp = 0;
    if (!is_init)
    {
        start_stamp = imu_data.timestamp - 1;
        is_init = true;
    }
    qDebug() << "\n[INFO] `selector`: timestamp:" << imu_data.timestamp - start_stamp;
    save_lidar_data(lidar_data, "LiDAR_Data", global_file_num);
    save_imu_data(imu_data, "IMU_Data", global_file_num);
    global_file_num++;
    if (circle_flag)
        return circle(lidar_data, imu_data);
    if (turn_right_flag)
        return turn_right_func(imu_data);
    if (turn_left_flag)
        return turn_left_func(imu_data);
    judge_turn(lidar_data);
    judge_circle(lidar_data, imu_data);
    return go_straight(lidar_data, imu_data);
}

//期末部分

bool IsCar = false;//是否检测到小车的全局变量
double threshold_distance1 = 500;
double threshold_distance2 = 100;
int error_cnt = 0;
int max_error_cnt = 5;

// *检测小车*
bool Detect_The_Car(const short *lidar_data){
    //检测激光雷达突变

    int start_detect = 70;
    int end_detetct = 290;

    int threshold = 250;
    int cnt = 0;

    for(int i = start_detect;i <= end_detetct;i++){
        if(int(lidar_data[i+1] - lidar_data[i]) > threshold || int(lidar_data[i] - lidar_data[i+1]) > threshold) cnt ++;
    }
    qDebug()<<"Mutation Times : "<<cnt;

    if(cnt >= 2){
        // *cnt > 2*一般不会出现
        qDebug()<<"Detect The Car";
        return true;
    }

    else if(cnt == 1){
        qDebug()<<"Detect The Wall";
        return false;
    }

    qDebug()<<"Detect Nothing";
    return false;
}

//*查找小车所在的突变*
std::pair<int,int> Find_Mutation(const short *lidar_data){
    if(!IsCar){
        qDebug()<<"There Is No Car";
        return std::pair<int,int>(-1,-1);
    }

    int Mutation[360];
    for(int i = 0;i < 360; i++){
        Mutation[i] = lidar_data[i+1] - lidar_data[i];
    }

    int Max_Distance = -100000;
    int Min_Distance = 100000;
    int start,end;//start : 小车纸板右边界对应的度数

    for(int i = 0 ; i < 360 ; i++){
        if(Mutation[i] < Min_Distance){
            start = i + 1;
            Min_Distance = Mutation[i];
        }
        if(Mutation[i] > Max_Distance){
            end = i;
            Max_Distance = Mutation[i];
        }
    }

    bool IF_PRINT = true;
    if(IF_PRINT){
        qDebug()<<"First Mutation : "<<"1. "<<lidar_data[start - 1]<<"2. "<<lidar_data[start]<<"3. "<<lidar_data[start + 1];
        qDebug()<<"Second Mutation : "<<"1. "<<lidar_data[end - 1]<<"2. "<<lidar_data[end]<<"3. "<<lidar_data[end + 1];
    }

    return std::pair<int,int>(start,end);
}

//*计算error*
std::pair<double,double> Detect_Rascal(const short *lidar_data){
    //(error,distance)
    double MiniLength = 5.0;
    double MaxLength = 40.0;

    std::pair<int,int> result = Find_Mutation(lidar_data);
    //计算夹角
    double angle = (result.second - result.first)*M_PI/360.0;  
    //计算激光雷达到纸板两端的距离
    int start = lidar_data[result.first];
    int end = lidar_data[result.second];
    //余弦定理计算纸板长度
    double PaperBoardLength = 0.0;
    PaperBoardLength = std::sqrt(start*start + end*end + 2*start*end*std::cos(angle));
    //在保证小车纸板有足够长度的前提下，这意味着前面检测小车的函数出现了失误
    if(PaperBoardLength < MiniLength || PaperBoardLength > MaxLength){
        qDebug()<<"PaperBoardLength : "<<PaperBoardLength<<" Too small or Too long , Maybe Something Wrong";
        return std::pair<double,double>(0.0,-1);
    }

    double angle1 = 0.0; // start与纸板的夹角
    angle1 = std::acos((start*start + PaperBoardLength*PaperBoardLength - end*end)/(2*start*PaperBoardLength));
    //计算激光雷达到纸板中点的距离
    double middle = 0.0;
    middle = std::sqrt(start*start + PaperBoardLength*PaperBoardLength*0.25 - start*PaperBoardLength*std::cos(angle1));

    if(middle > threshold_distance1){
        qDebug()<<"distance to the car too long";
        return std::pair<double,double>(0.0,middle);
    }
    else if(middle < threshold_distance2){
        qDebug()<<"distance too short";
        return std::pair<double,double>(0.0,middle);
    }
    else{
        double error = middle*cos(result.first + angle/2);
        qDebug()<<"Caculate error: "<<error;
        return std::pair<double,double>(error,middle);
    }

    //以防没有返回值
    return std::pair<double,double>(0.0,0.0);
}

//*寻找小车*
output_type Find_Car(const short *lidar_data,imu_data_type imu_data){
    qDebug()<<"Find Car ";
    //为转弯作准备
    bool front_barrier = false;
    int counter = 0;
    for (int i = 150; i <= 210; i++){
        if (lidar_data[i] > 0 && lidar_data[i] < 220) counter++;
    }
    if (counter > 50 ) front_barrier = true;
    int threshold = 600;
    if (front_barrier) threshold = 300;

    counter = 0;
    //优先向一个方向转弯
    for (int i = 0; i < 30; i++){
        if (lidar_data[i] == 0 || lidar_data[i] > threshold) counter++;
    }
    if(counter >= 20){
        qDebug()<<"Turn Right";
        return turn_right_func(imu_data);
    } 

    counter = 0;
    for (int i = 360; i >= 330; i--){
        if (lidar_data[i] == 0 || lidar_data[i] > threshold) counter++;
    }
    if(counter >= 20){
        qDebug()<<"Turn Left";
        return turn_left_func(imu_data);
    } 

    //既不左转也不右转
    return go_straight(lidar_data,imu_data);
}

//*跟随小车*
output_type Follow_Car(const short *lidar_data,imu_data_type imu_data){
    qDebug()<<"Follow Car ";

    if(error_cnt == max_error_cnt){//发现纸板长度异常多次，可能把墙的凸起识别为了小车
        qDebug()<<"Confirm Detect Wrong thing , not the car";
        error_cnt = 0;
        return Find_Car(lidar_data,imu_data);
    }

    std::pair<double,double> result = Detect_Rascal(lidar_data);
    if(result.second < 0){//正常直行
        error_cnt ++;
        return go_straight(lidar_data,imu_data);
    }
    else if(result.second > threshold_distance1){
        return go_straight(lidar_data,imu_data);
    }
    else if(result.second < threshold_distance2){
        //距离过近，暂时减速或者停车
        int input_speed = 0;
        return go_straight(lidar_data,imu_data,false,0,false,0,false,0,true,input_speed);
    }
    else{
        return go_straight(lidar_data,imu_data,false,0,false,0,true,result.first);
    }
}

// *行为选择*
output_type selector_1(const short *lidar_data, imu_data_type imu_data){
    static bool is_init = false;
    static int start_stamp = 0;
    if (!is_init)
    {
        start_stamp = imu_data.timestamp - 1;
        is_init = true;
    }
    qDebug() << "\n[INFO] `selector_1`: timestamp:" << imu_data.timestamp - start_stamp;

    IsCar = Detect_The_Car(lidar_data);

    if(IsCar){
        return Follow_Car(lidar_data,imu_data);
    }
    else{
        return Find_Car(lidar_data,imu_data);
    }
}
