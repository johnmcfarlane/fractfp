#include <vector>

#include <sg14/auxiliary/boost.simd.h>

// better to ensure vertical strips and have a single y value?

namespace mandelbrot {
    constexpr auto pack_shift = 2;
    constexpr auto pack_n = 1 << pack_shift;
    constexpr auto pack_mask = pack_n - 1;

    constexpr auto point_to_pack_index(auto point) noexcept { return point >> pack_shift; }
    constexpr auto point_to_pack_element(auto point) noexcept { return point & pack_mask; }

    template<typename T>
    void nuke_vector(std::vector<T>& v) noexcept {
        std::memset(v.data(), 0, sizeof(T) * v.size());
    }

    template<typename T, std::size_t N = pack_n>
    using pack = boost::simd::pack<T, N>;

    template<typename T>
    using pack_vector = std::vector<pack<T>>;

    template<typename T>
    using vec2 = pack<T, 2>;

    template<typename Rep, typename Counter = int>
    pack<Counter> generate(pack<Rep> const c_x, pack<Rep> const c_y, Counter limit) noexcept {
        pack<Counter> counters{0};

        auto x = pack<Rep>{pack<Rep>(c_x)};
        auto y = pack<Rep>{pack<Rep>(c_y)};
        auto confined = pack<boost::simd::logical<Rep>>{true};

        while (limit--) {
            auto xx = x*x;
            auto yy = y*y;

            auto currently_confined = xx + yy <= 4;
            if (boost::simd::none(currently_confined)) {
                break;
            }
            confined = currently_confined & confined;
            counters = boost::simd::if_inc(confined, counters);

            auto x_temp = xx - yy + c_x;
            y = 2*x*y + c_y;
            x = x_temp;
        }
        return counters;
    };

    template<typename Rep, typename Counter = int>
    std::vector<Counter> generate(
            vec2<Rep> const origin, vec2<Rep> const x_extent, vec2<Rep> const y_extent,
            vec2<int> const resolution, Counter const limit) {
        // dimensions
        auto num_points = static_cast<std::size_t>(resolution[0]) * resolution[1];
        auto num_packs = (num_points + pack_mask) / pack_n;

        // results
        auto subset = std::vector<Counter>(num_points);

        // c
        auto c_x = pack_vector<Rep>{num_packs};
        auto c_y = pack_vector<Rep>{num_packs};

        auto dx = vec2<Rep>{ x_extent[0] / resolution[0], x_extent[1] / resolution[1] };
        auto dy = vec2<Rep>{ y_extent[0] / resolution[0], y_extent[1] / resolution[1] };
        for (auto row = 0; row != resolution[1]; ++row) {
            auto line_index = row * resolution[0];
            for (auto column = 0; column != resolution[0]; ++column) {
                auto point_index = line_index + column;
                auto pack_index = point_to_pack_index(point_index);
                auto pack_element = point_to_pack_element(point_index);

                c_x[pack_index][pack_element] = origin[0] + column * dx[0] + row * dy[0];
                c_y[pack_index][pack_element] = origin[1] + column * dx[1] + row * dy[1];
            }
        }

        auto counters = pack_vector<Counter>{num_packs, pack<Counter>{0}};
        for (auto i = 0; i != num_packs; ++ i) {
            counters[i] = generate(c_x[i], c_y[i], limit);
        }

        // extract results
        for (auto out_index = 0; out_index != num_points; ++ out_index) {
            subset[out_index] = counters[out_index >> pack_shift][out_index & pack_mask];
        }

        return subset;
    };
}

int main() {
    auto resolution = mandelbrot::vec2<int>{128, 64};
    auto subset = mandelbrot::generate<float>({-2, -2}, {4, 0}, {0, 4}, resolution, 65536);

    for (auto y = 0; y != resolution[1]; ++ y) {
        auto line_index = y * resolution[0];
        for (auto x = 0; x != resolution[0]; ++ x) {
            auto index = line_index + x;
            auto iterations = subset[index];
            std::putchar((iterations>=127)?' ':' ' + iterations);
        }
        std::putchar(' '+y);
        std::putchar('\n');
    }

    return 0;
}