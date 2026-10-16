#pragma once
#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

// Benchmark host scheduler. Thread lifetime is outside measured physics steps.
class TaskPool {
  public:
    explicit TaskPool(unsigned count) {
        for (unsigned i=0;i<count;++i) workers_.emplace_back([this] {
            std::size_t seen = 0;
            std::unique_lock lock(mutex_);
            for (;;) {
                ready_.wait(lock, [&] { return stopping_ || epoch_ != seen; });
                if (stopping_) return;
                seen = epoch_;
                lock.unlock(); consume(); lock.lock();
                if (--remaining_ == 0) finished_.notify_one();
            }
        });
    }
    ~TaskPool() {
        { std::lock_guard lock(mutex_); stopping_ = true; }
        ready_.notify_all();
        for (auto &worker : workers_) worker.join();
    }
    void run(std::size_t count, const std::function<void(std::size_t)> &job) {
        {
            std::lock_guard lock(mutex_);
            count_ = count; next_ = 0; job_ = &job; error_ = {};
            remaining_ = workers_.size(); ++epoch_;
        }
        ready_.notify_all();
        consume();
        std::unique_lock lock(mutex_);
        finished_.wait(lock, [&] { return remaining_ == 0; });
        if (error_) std::rethrow_exception(error_);
    }
  private:
    void consume() {
        for (;;) {
            const auto i = next_.fetch_add(1, std::memory_order_relaxed);
            if (i >= count_) return;
            try { (*job_)(i); }
            catch (...) { std::lock_guard lock(mutex_); if (!error_) error_ = std::current_exception(); }
        }
    }
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable ready_, finished_;
    std::atomic<std::size_t> next_{0};
    std::size_t epoch_{0}, count_{0}, remaining_{0};
    const std::function<void(std::size_t)> *job_{};
    std::exception_ptr error_;
    bool stopping_{false};
};
