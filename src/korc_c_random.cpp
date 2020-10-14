#include <random>

class random_U {
    std::mt19937_64 engine_U;
    std::uniform_real_distribution<double> dist_U;
    
    public:
    random_U(uint64_t seed_U) : engine_U(seed_U), dist_U(0.0, 1.0) {};
    double get_number_U() {
        return dist_U(engine_U);
    }
};

class random_N {
    std::mt19937_64 engine_N;
    std::normal_distribution<double> dist_N;
    
    public:
    random_N(uint64_t seed_N) : engine_N(seed_N), dist_N(0.0, 1.0) {};
    double get_number_N() {
        return dist_N(engine_N);
    }
};

extern "C" {
    void *random_construct_U(int seed_N) {
        return new class random_U(static_cast<uint64_t> (seed_N));
    }

    void *random_construct_N(int seed_N) {
        return new class random_N(static_cast<uint64_t> (seed_N));
    }
  
    double random_get_number_U(void *r) {
        return static_cast<class random_U *> (r)->get_number_U();
    }

    double random_get_number_N(void *r) {
        return static_cast<class random_N *> (r)->get_number_N();
    }
  
    void random_N_destroy(void *r) {
        delete  static_cast<class random_N *> (r);
    }
  
    void random_U_destroy(void *r) {
        delete  static_cast<class random_U *> (r);
    }
}
