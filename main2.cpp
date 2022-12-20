#include <iostream>
#include <cstdint>
#include <vector>
#include <queue>
#include <limits>
#include <functional>
#include <chrono>
#include <atomic>
#include <iomanip>

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/parallel.h>

struct graph_type {
    using vertex_id_type = uint32_t;
    using adj_list_type = parlay::sequence<vertex_id_type>;
    using adj_lists_type = parlay::sequence<adj_list_type>;

    vertex_id_type vertex_count;
    vertex_id_type n;

    explicit graph_type(vertex_id_type n) : vertex_count(n * n * n), n(n) {}

    [[nodiscard]] adj_list_type adj_lists(vertex_id_type v) const {
        vertex_id_type i, j, k;
        i = v / (n * n);
        j = (v / n) % n;
        k = v % n;

        parlay::sequence<vertex_id_type> result;
        if (k > 0)
            result.push_back(i * n * n + j * n + k - 1);
        if (k < n - 1)
            result.push_back(i * n * n + j * n + k + 1);
        if (j > 0)
            result.push_back(i * n * n + (j - 1) * n + k);
        if (j < n - 1)
            result.push_back(i * n * n + (j + 1) * n + k);
        if (i > 0)
            result.push_back((i - 1) * n * n + j * n + k);
        if (i < n - 1)
            result.push_back((i + 1) * n * n + j * n + k);
        return result;
    }

    [[nodiscard]] uint32_t adj_lists_size(vertex_id_type v) const {
        vertex_id_type i, j, k;
        i = v / (n * n);
        j = (v / n) % n;
        k = v % n;

        uint32_t result = 0;
        if (k > 0)
            result++;
        if (k < n - 1)
            result++;
        if (j > 0)
            result++;
        if (j < n - 1)
            result++;
        if (i > 0)
            result++;
        if (i < n - 1)
            result++;
        return result;
    }
};

using vtype = graph_type::vertex_id_type;
using dtype = vtype;
constexpr dtype INF = std::numeric_limits<dtype>::max();
constexpr vtype VINF = std::numeric_limits<vtype>::max();

inline graph_type generate_cube(vtype n) {
    return graph_type(n);
}

inline std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

std::pair<uint64_t, parlay::sequence<dtype>>
bench_micros(graph_type const &graph, const std::function<parlay::sequence<dtype>(graph_type const&, vtype)> &bfs) {
    auto start = now();
    parlay::sequence<dtype> dist = bfs(graph, 0);
    auto finish = now();

    return {
            std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count(),
            dist
    };
}

parlay::sequence<dtype> bfs_sequential(graph_type const &graph, vtype start) {
    std::queue<vtype> queue;
    parlay::sequence<dtype> dist(static_cast<size_t>(graph.vertex_count), INF);
    queue.push(start);
    dist[queue.front()] = 0;

    while (!queue.empty()) {
        vtype v = queue.front();
        queue.pop();

        for (auto to: graph.adj_lists(v)) {
            if (dist[to] == INF) {
                dist[to] = dist[v] + 1;
                queue.push(to);
            }
        }
    }

    return dist;
}

#undef EXAMPLE
#undef PIZDEC
#undef INNER_PARALLEL

parlay::sequence<dtype> bfs_parallel(graph_type const &graph, vtype start) {
#ifdef EXAMPLE
    auto visited = parlay::tabulate<std::atomic<bool>>(graph.vertex_count, [&] (long i) {
        return (i == start); });

    parlay::sequence<dtype> dist(static_cast<size_t>(graph.vertex_count), INF);
    dist[start] = 0;

    parlay::sequence<vtype> frontier(1,start);
    dtype cur_dist = 1;
    while (!frontier.empty()) {
        // get out edges of the frontier and flatten
        auto out = parlay::flatten(parlay::map(frontier, [&] (vtype u) {return graph.adj_lists[u];}));

        // keep the v that succeed in setting the visited array
        frontier = filter(out, [&] (auto&& v) {
            bool expected = false;
            if ((!visited[v]) && visited[v].compare_exchange_strong(expected, true)) {
                dist[v] = cur_dist;
                return true;
            } else
                return false;
        });

        cur_dist++;
    }

    return dist;
#else
#ifdef PIZDEC
    parlay::sequence<vtype> front;
    front.push_back(start);

    std::vector<std::atomic_bool> used(graph.vertex_count);
    used[start] = true;

    parlay::sequence<dtype> dist(static_cast<size_t>(graph.vertex_count), INF);
    dist[start] = 0;

    while (!front.empty()) {
        front = parlay::filter(std::move(parlay::flatten(std::move(parlay::map(std::move(front), [&graph, &used, &dist](vtype v) {
            auto const &to_list = graph.adj_lists[v];
            return parlay::tabulate(to_list.size(), [&to_list, &used, &dist, v] (long i) {
                vtype to = to_list[i];
                bool expected = false;
                if (used[to].compare_exchange_strong(expected, true)) {
                    dist[to] = dist[v] + 1;
                    return static_cast<vtype>(to);
                } else
                    return VINF;
            });
        })))), [](vtype v) {
            return v < VINF;
        });
    }

    return dist;
#else
    parlay::sequence<vtype> front;
    front.push_back(start);

    std::vector<std::atomic_bool> used(graph.vertex_count);
    used[start] = true;

    parlay::sequence<dtype> dist(static_cast<size_t>(graph.vertex_count), INF);
    dist[start] = 0;

    dtype cur_dist = 1;
    while (!front.empty()) {
        auto front_size = front.size();

        auto deg = parlay::map(front, [&graph](vtype v) {
            return static_cast<vtype>(graph.adj_lists_size(v));
        });

        vtype new_front_size = parlay::scan_inplace(deg, parlay::plus<vtype>());

        parlay::sequence<vtype> vui(deg);

        parlay::sequence<vtype> new_front_raw(static_cast<size_t>(new_front_size), VINF);

        parlay::parallel_for(0, front_size, [&front, &graph, &used, &new_front_raw, &dist, &vui, cur_dist](size_t i) {
            vtype v = front[i];
            auto to_list = graph.adj_lists(v);
#ifdef INNER_PARALLEL
            parlay::parallel_for(0, to_list.size(), [&i, &to_list, &used, &new_front_raw, &dist, &vui, cur_dist](size_t j) {
#else
            for (int j = 0; j < to_list.size(); ++j) { // NOLINT(modernize-loop-convert)
#endif
                vtype to = to_list[j];
                bool expected = false;
                if (used[to].compare_exchange_strong(expected, true)) {
                    new_front_raw[vui[i]++] = to;
                    dist[to] = cur_dist;
                }
#ifdef INNER_PARALLEL
            });
#else
            }
#endif
        });

        cur_dist++;

        front = parlay::filter(new_front_raw, [](vtype v){return v < VINF;});
    }

    return dist;
#endif
#endif
}

const int N[] = {50, 100, 250, 500};
const int ATTEMPTS = 5;

int main() {
    for (auto n: N) {
        std::cout << "Current size: " << n << std::endl;
        uint64_t par_micros_sum = 0;
        uint64_t seq_micros_sum = 0;
        auto graph = generate_cube(n);

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
        }

        std::cout << "sequential micros: " << (seq_micros_sum / ATTEMPTS) << std::endl;
        std::cout << "  parallel micros: " << (par_micros_sum / ATTEMPTS) << std::endl;
        std::cout << "Boost: " << std::setprecision(6)
                  << static_cast<double>(seq_micros_sum) / static_cast<double>(par_micros_sum) << std::endl;
    }

    return 0;
}
