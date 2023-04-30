#ifndef READMAPPING_LOGGING_H
#define READMAPPING_LOGGING_H

#include <string>
#include <iostream>
#include <mutex>

enum class LogLevel {
    DEBUG,
    INFO,
    WARN,
    ERROR,
    INFO_ONLY
};


class Logger {
public:
    Logger() : log_level(LogLevel::INFO) {}

    explicit Logger(LogLevel level) : log_level(level) {}

    inline void
    log_msg(LogLevel level, const std::string &func, const std::string &msg,
            std::basic_ostream<char> &stream = std::cout) {
        if ((log_level != LogLevel::INFO_ONLY && level < log_level)
            || (log_level == LogLevel::INFO_ONLY && level != LogLevel::INFO)) {
            return;
        }


        std::string level_str, color_str, color_end = "\033[0m";
        switch (level) {
            case LogLevel::DEBUG:
                level_str = "DEBUG";
                color_str = "\033[0;37m";
                break;
            case LogLevel::INFO:
                level_str = "INFO";
                color_str = "";
                color_end = "";
                break;
            case LogLevel::WARN:
                level_str = "WARN";
                color_str = "\033[0;33m";
                break;
            case LogLevel::ERROR:
                level_str = "ERROR";
                color_str = "\033[0;31m";
                break;
            case LogLevel::INFO_ONLY:
                throw std::runtime_error("LogLevel::INFO_ONLY should not be used in log_msg");
        }

        std::lock_guard<std::mutex> lock(log_mutex);
        stream << "\r" << std::string(70 + 7, ' ') << "\r" << std::flush;
        stream << color_str << "[" << level_str << "] " << func << ": " << msg << color_end << std::endl;

        if (progress_bar_enabled) {
            display_progress_bar(false);
        }
    }

    inline void
    display_progress_bar(bool need_lock = false) {
        int percent = progress_current * 100 / progress_total;
        std::string bar;
        for (int i = 0; i < 50; i++) {
            if (i < (percent / 2)) {
                bar.replace(i, 1, "=");
            } else if (i == (percent / 2)) {
                bar.replace(i, 1, ">");
            } else {
                bar.replace(i, 1, " ");
            }
        }

        if (need_lock) std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "\r";
        std::cout << progress_current << "/" << progress_total << " ";
        std::cout << "[" << bar << "] ";
        std::cout.width(3);
        std::cout << percent << "%     " << std::flush;
    }

    inline void
    create_progress_bar(int total) {
        progress_total = total;
        progress_current = 0;
        progress_bar_enabled = true;
        display_progress_bar();
    }

    inline void
    update_progress_bar(int incr = 1) {
        std::lock_guard<std::mutex> lock(incr_mutex);
        progress_current += incr;
        display_progress_bar();
    }

    inline void
    finish_progress_bar() {
        progress_current = progress_total;
        display_progress_bar();
        progress_bar_enabled = false;
        std::cout << std::endl;
    }

private:
    LogLevel log_level;

    std::mutex log_mutex;
    std::mutex incr_mutex;
    int progress_total = 0;
    int progress_current = 0;
    bool progress_bar_enabled = false;
};

#endif //READMAPPING_LOGGING_H
