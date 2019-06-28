#include <random>
#include <chrono>

class random {
    std::mt19937_64 engine;
    std::uniform_real_distribution<double> dist;
    
    public:
    random(const uint64_t offset) :
    engine(static_cast<uint64_t> (std::chrono::system_clock::to_time_t(
                                      std::chrono::system_clock::now())) +
                                  offset),
    dist(0.0, 1.0) {};
    void set_dist(const double low, const double high) {
        dist = std::uniform_real_distribution<double> (low, high);
    }
    double get_number() {
        return dist(engine);
    }
};

extern "C" {
    void *random_construct(const int seed) {
        return new class random(static_cast<uint64_t> (seed));
    }

    void random_set_dist(void *r,
                         const double low,
                         const double high) {
        static_cast<class random *> (r)->set_dist(low, high);
    }

    double random_get_number(void *r) {
        return static_cast<class random *> (r)->get_number();
    }
    
    void random_destroy(void *r) {
        delete  static_cast<class random *> (r);
    }
}
