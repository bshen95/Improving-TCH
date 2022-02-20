//
// Created by Bojie Shen on 10/8/20.
//


#include <chrono>
#include "my_timer.h"

my_timer::my_timer() = default;


double my_timer::get_time_nano()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();
}

void my_timer::start()
{
    start_time = std::chrono::steady_clock::now();
    stop_time = start_time;
}

void my_timer::stop()
{
    stop_time =  std::chrono::steady_clock::now();
}

double my_timer::current_time_nano(std::chrono::time_point<std::chrono::steady_clock> current_time){
    return std::chrono::duration_cast<std::chrono::nanoseconds>(current_time - start_time).count();
}
double my_timer::elapsed_time_nano()
{

    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();
}

void my_timer::reset()
{
}

double my_timer::elapsed_time_micro()
{
    return elapsed_time_nano() / 1000.0;
}
double my_timer::elapsed_time_mins()
{
    return elapsed_time_nano() / 1000.0/ 1000000.0 / 60.0;;
}
