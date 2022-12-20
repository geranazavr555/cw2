#define CILK_ENABLED
#ifdef CILK_ENABLED
#include <cilk/cilk.h>
#else
#define cilk_for for
#define cilk_spawn
#define cilk_sync
#endif
#include <iostream>
#include <cstdint>
#include <vector>
#include <queue>
#include <limits>
#include <functional>
#include <chrono>
#include <atomic>
#include <iomanip>

struct graph_type {
    using vertex_id_type = uint32_t;
    using adj_list_type = std::vector<vertex_id_type>;
    using adj_lists_type = std::vector<adj_list_type>;

    vertex_id_type vertex_count;
    adj_lists_type adj_lists;

    explicit graph_type(vertex_id_type n) : vertex_count(n), adj_lists(n) {}
};

using vtype = graph_type::vertex_id_type;
using dtype = vtype;
constexpr dtype INF = std::numeric_limits<dtype>::max();
constexpr vtype VINF = std::numeric_limits<vtype>::max();

graph_type *generate_cube(vtype n) {
    auto ***tmp = new vtype **[n];
    vtype vid = 0;
    for (vtype i = 0; i < n; ++i) {
        tmp[i] = new vtype *[n];
        for (vtype j = 0; j < n; ++j) {
            tmp[i][j] = new vtype[n];
            for (vtype k = 0; k < n; ++k) {
                tmp[i][j][k] = vid++;
            }
        }
    }

    auto *graph = new graph_type(vid);
    for (vtype i = 0; i < n; ++i) {
        for (vtype j = 0; j < n; ++j) {
            for (vtype k = 0; k < n; ++k) {
                vtype cur = tmp[i][j][k];
                if (k > 0)
                    graph->adj_lists[cur].push_back(tmp[i][j][k - 1]);
                if (k < n - 1)
                    graph->adj_lists[cur].push_back(tmp[i][j][k + 1]);
                if (j > 0)
                    graph->adj_lists[cur].push_back(tmp[i][j - 1][k]);
                if (j < n - 1)
                    graph->adj_lists[cur].push_back(tmp[i][j + 1][k]);
                if (i > 0)
                    graph->adj_lists[cur].push_back(tmp[i - 1][j][k]);
                if (i < n - 1)
                    graph->adj_lists[cur].push_back(tmp[i + 1][j][k]);
            }
        }
    }

    for (vtype i = 0; i < n; ++i) {
        for (vtype j = 0; j < n; ++j)
            delete[] tmp[i][j];
        delete[] tmp[i];
    }
    delete[] tmp;

    return graph;
}

inline std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

std::pair<uint64_t, dtype *>
bench_micros(graph_type const *graph, const std::function<dtype *(graph_type const *, vtype)> &bfs) {
    auto start = now();
    dtype *dist = bfs(graph, 0);
    auto finish = now();

    return {
            std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count(),
            dist
    };
}

dtype *bfs_sequential(graph_type const *graph, vtype start) {
    std::queue<vtype> queue;
    auto n = graph->vertex_count;
    auto *dist = new dtype[graph->vertex_count];
    std::fill(dist, dist + n, INF);
    queue.push(start);
    dist[queue.front()] = 0;

    while (!queue.empty()) {
        vtype v = queue.front();
        queue.pop();

        for (auto to: graph->adj_lists[v]) {
            if (dist[to] == INF) {
                dist[to] = dist[v] + 1;
                queue.push(to);
            }
        }
    }

    return dist;
}

using atype = vtype;

#undef NOT_DELETE

#undef PARALLEL_FILTER

#ifdef PARALLEL_FILTER

const atype BLOCK_SIZE = 1000;

template<typename T>
void filter_impl(T const *a, atype l, atype r, T*& result, atype& result_size) {
    if (l + BLOCK_SIZE >= r) {
        result = new T[r - l];
        result_size = 0;
        for (atype i = l; i < r; ++i) {
            if (a[i] < VINF)
                result[result_size++] = a[i];
        }
    } else {
        atype m = (l + r) / 2;
        T* left_result;
        T* right_result;
        atype left_size, right_size;
        cilk_spawn filter_impl(a, l, m, left_result, left_size);
        filter_impl(a, m, r,right_result, right_size);
        cilk_sync;

        result = new T[left_size + right_size];
        cilk_for (int i = 0; i < left_size + right_size; ++i) {
            if (i < left_size)
                result[i] = left_result[i];
            else
                result[i] = right_result[i - left_size];
        }

#ifndef NOT_DELETE
        delete[] left_result;
        delete[] right_result;
#endif

        result_size = left_size + right_size;
    }
}

#endif

template<typename T>
std::vector<T> *filter(T const *a, atype n) {
#ifdef PARALLEL_FILTER
    T* result;
    atype result_size;
    filter_impl(a, 0, n, result, result_size);
    auto *vec = new std::vector<T>(result_size);
    cilk_for (int i = 0; i < result_size; ++i)
        (*vec)[i] = result[i];
#ifndef NOT_DELETE
    delete[] result;
#endif
    return vec;
#else
    auto *result = new std::vector<T>();
    for (atype i = 0; i < n; ++i) {
        if (a[i] < VINF)
            result->push_back(a[i]);
    }

    return result;
#endif
}

#undef PARALLEL_SCAN

#ifdef PARALLEL_SCAN
inline atype mylog(atype x) {
    atype cnt = 0;
    while (x > 0) {
        ++cnt;
        x >>= 1;
    }
    return cnt;
}
#endif

template<typename T>
T *scan(T const *a, atype n) {
#ifdef PARALLEL_SCAN
    T* step0 = new T[n];
    T* step1 = new T[n];
    T* steps[] = {step0, step1};

    cilk_for (atype i = 0; i < n; ++i) {
        step0[i] = a[i];
    }

    auto steps_cnt = mylog(n - 1);

    for (atype i = 0; i < steps_cnt; ++i) {
        cilk_for (atype j = 0; j < n; ++j) {
            if (j < (1u << i)) {
                steps[(i + 1) & 1u][j] = steps[i & 1u][j];
            } else {
                steps[(i + 1) & 1u][j] = steps[i & 1u][j] + steps[i & 1u][j - (1u << i)];
            }
        }
    }

#ifndef NOT_DELETE
    delete[] steps[(steps_cnt + 1) & 1];
#endif
    return steps[steps_cnt & 1];

#else
    T *pref = new T[n];
    pref[0] = 0;
    for (int i = 1; i < n; ++i)
        pref[i] = pref[i - 1] + a[i];
    return pref;
#endif
}

#define BENCHMARK

dtype *bfs_parallel(graph_type const *graph, vtype start) {
    auto *front = new std::vector<vtype>;
    front->push_back(start);
    auto *used = new std::atomic_bool[graph->vertex_count];
    used[start] = true;

    auto *dist = new dtype[graph->vertex_count];
    cilk_for (int i = 0; i < graph->vertex_count; ++i) {
        dist[i] = INF;
    }
    dist[start] = 0;

#ifdef BENCHMARK
    uint64_t sum_micros1 = 0;
    uint64_t sum_micros2 = 0;
    uint64_t sum_micros3 = 0;
    uint64_t sum_micros4 = 0;
    uint64_t sum_micros5 = 0;
    uint64_t sum_micros6 = 0;
#endif

    while (!front->empty()) {
#ifdef BENCHMARK
        auto start1 = now();
#endif
        auto front_size = front->size();
        auto *deg = new dtype[front_size];
        cilk_for (int i = 0; i < front_size; ++i) {
            vtype v = (*front)[i];
            deg[i] = graph->adj_lists[v].size();
        }

#ifdef BENCHMARK
        auto finish1 = now();
        sum_micros1 += std::chrono::duration_cast<std::chrono::microseconds>(finish1 - start1).count();
#endif

#ifdef BENCHMARK
        auto start2 = now();
#endif
        auto *deg_pref = scan(deg, front_size);
#ifdef BENCHMARK
        auto finish2 = now();
        sum_micros2 += std::chrono::duration_cast<std::chrono::microseconds>(finish2 - start2).count();
#endif

#ifdef BENCHMARK
        auto start3 = now();
#endif
        auto* vui = deg;

        cilk_for (int i = 0; i < front_size; ++i) {
            vui[i] = deg_pref[i];
        }
#ifdef BENCHMARK
        auto finish3 = now();
        sum_micros3 += std::chrono::duration_cast<std::chrono::microseconds>(finish3 - start3).count();
#endif

        dtype new_front_size = deg_pref[front_size - 1] + graph->adj_lists[front->back()].size();
        auto *new_front_raw = new vtype[new_front_size];

#ifdef BENCHMARK
        auto start4 = now();
#endif
        cilk_for (int i = 0; i < new_front_size; ++i) {
            new_front_raw[i] = VINF;
        }
#ifdef BENCHMARK
        auto finish4 = now();
        sum_micros4 += std::chrono::duration_cast<std::chrono::microseconds>(finish4 - start4).count();
#endif

#ifdef BENCHMARK
        auto start5 = now();
#endif
        cilk_for (int i = 0; i < front_size; ++i) {
            vtype v = (*front)[i];
            auto const &to_list = graph->adj_lists[v];
            cilk_for (int j = 0; j < to_list.size(); ++j) { // NOLINT(modernize-loop-convert)
                vtype to = to_list[j];
                bool expected = false;
                if (used[to].compare_exchange_strong(expected, true)) {
                    new_front_raw[vui[i]++] = to;
                    dist[to] = dist[v] + 1;
                }
            }
        }

#ifdef BENCHMARK
        auto finish5 = now();
        sum_micros5 += std::chrono::duration_cast<std::chrono::microseconds>(finish5 - start5).count();
#endif

#ifndef NOT_DELETE
        delete[] deg_pref;
        delete[] vui;
        delete front;
#endif
#ifdef BENCHMARK
        auto start6 = now();
#endif
        front = filter(new_front_raw, new_front_size);
#ifdef BENCHMARK
        auto finish6 = now();
        sum_micros6 += std::chrono::duration_cast<std::chrono::microseconds>(finish6 - start6).count();
#endif
#ifndef NOT_DELETE
        delete[] new_front_raw;
#endif
    }

#ifdef BENCHMARK
    std::cerr << "deg build            " << sum_micros1 << std::endl;
    std::cerr << "scans                " << sum_micros2 << std::endl;
    std::cerr << "vui build            " << sum_micros3 << std::endl;
    std::cerr << "new_front_raw init   " << sum_micros4 << std::endl;
    std::cerr << "main loop            " << sum_micros5 << std::endl;
    std::cerr << "filters              " << sum_micros6 << std::endl;
    std::cerr << "------------------------------------" << std::endl;
#endif

#ifndef NOT_DELETE
    delete[] used;
    delete front;
#endif

    return dist;
}

const int N[] = {200};
const int ATTEMPTS = 5;

int main() {
    for (auto n: N) {
        std::cout << "Current size: " << n << std::endl;
        uint64_t par_micros_sum = 0;
        uint64_t seq_micros_sum = 0;
        auto const *graph = generate_cube(n);

        for (int j = 0; j < ATTEMPTS; ++j) {
            auto [seq_micros, seq_array] = bench_micros(graph, bfs_sequential);
            seq_micros_sum += seq_micros;

            auto [par_micros, par_array] = bench_micros(graph, bfs_parallel);
            par_micros_sum += par_micros;

            for (size_t k = 0; k < n * n * n - 1; ++k)
                if (par_array[k] != seq_array[k]) {
                    std::cout << "ERROR par in position " << k << std::endl;
                    exit(1);
                }

            delete[] par_array;
            delete[] seq_array;
        }

        std::cout << "sequential micros: " << (seq_micros_sum / ATTEMPTS) << std::endl;
        std::cout << "  parallel micros: " << (par_micros_sum / ATTEMPTS) << std::endl;
        std::cout << "Boost: " << std::setprecision(6)
                  << static_cast<double>(seq_micros_sum) / static_cast<double>(par_micros_sum) << std::endl;
        delete graph;
    }

    return 0;
}
