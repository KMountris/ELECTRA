/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/utilities/timer.hpp"


namespace ELECTRA {

Timer::Timer() : beg_(timer_clock_::now())
{}


Timer::~Timer()
{}


void Timer::Reset()
{
    this->beg_ = timer_clock_::now();
}


double Timer::ElapsedMilliSecs() const
{
    return std::chrono::duration_cast<timer_millisec_>
            (timer_clock_::now() - beg_).count();
}


double Timer::ElapsedSecs() const
{
    return std::chrono::duration_cast<timer_second_>
            (timer_clock_::now() - beg_).count();
}


double Timer::ElapsedMinutes() const
{
    return std::chrono::duration_cast<timer_second_>
            (timer_clock_::now() - beg_).count() / 60.;
}


double Timer::ElapsedHours() const
{
    return std::chrono::duration_cast<timer_second_>
            (timer_clock_::now() - beg_).count() / 3600.;
}


std::string Timer::PrintElapsedTime() const
{
    if (this->ElapsedMilliSecs() < 1000.) { return std::to_string(this->ElapsedMilliSecs()) + " ms"; }
    else if (this->ElapsedMilliSecs() > 1000. &&
             this->ElapsedSecs() < 60.) { return std::to_string(this->ElapsedSecs()) + " s"; }
    else if (this->ElapsedSecs() > 60. &&
             this->ElapsedSecs() < 3600.) { return std::to_string(this->ElapsedMinutes()) + " mins"; }
    else { return std::to_string(this->ElapsedHours()) + " hours"; }
}


}  //end of namespace ELECTRA
