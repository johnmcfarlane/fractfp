#include "mandelbrot.h"

#include <sg14/auxiliary/boost.simd.h>

namespace {
    // zero = no override
    constexpr auto override_pack_size = std::size_t{};

    constexpr auto use_fixed_point = true;
    constexpr auto print_stats = true;
}

namespace {
    template<typename Scalar>
    void print_info(mandelbrot::Geometry<Scalar> geometry, int max_iterations) {
        if (!print_stats) {
            return;
        }

        using pack_traits = mandelbrot::pack_traits<Scalar>;
        std::printf(
                "set[%d][%d] block[%d][%d] pack[%d] iterations:%d\n",
                geometry.resolution[1],
                geometry.resolution[0],
                pack_traits::block_height,
                pack_traits::block_width,
                pack_traits::pack_size,
                max_iterations);
    }

    void display(int const *subset_first, mandelbrot::vec2<int> resolution, int max_iterations) {
        auto const max_displayed_iterations = std::min(max_iterations, 96);
        for (auto y = 0; y != resolution[1]; ++y) {
            for (auto x = 0; x != resolution[0]; ++subset_first, ++x) {
                auto iterations = *subset_first;
                std::putchar((iterations >= max_displayed_iterations) ? ' ' : ' ' + iterations);
            }
            std::putchar('\n');
        }
    }
}

namespace mandelbrot {
    template<typename Rep, int Exponent>
    struct scalar_traits<sg14::fixed_point<Rep, Exponent>> {
        using scalar = sg14::fixed_point<Rep, Exponent>;

        using _fixed_point_rep = std::conditional_t<
                override_pack_size,
                boost::simd::pack<Rep, override_pack_size>,
                boost::simd::pack<Rep>>;

        using scalar_pack = sg14::fixed_point<_fixed_point_rep, Exponent>;

        static constexpr auto pack_size = scalar_pack::rep::static_size;

        static void set(scalar_pack &fpp, int index, scalar const &value) noexcept {
            fpp.data()[index] = value.data();
        }
    };
}

int main() {
    using namespace mandelbrot;

    using scalar = std::conditional_t<
            use_fixed_point,
            sg14::fixed_point<std::int64_t, -28>,
            float>;

    auto const resolution = vec2<int>{80, 40};

    using coordinate = typename pack_traits<scalar>::coordinate;
    auto geometry = Geometry<scalar>{
            coordinate{-2, -2},
            coordinate{4, 0}, coordinate{0, 4},
            resolution
    };
    auto const max_iterations = 1000000;

    print_info<>(geometry, max_iterations);
    auto subset = generate<>(geometry, max_iterations);

    display(subset.get(), geometry.resolution, max_iterations);

    return 0;
}
