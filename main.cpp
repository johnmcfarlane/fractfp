#include "mandelbrot.h"

#include <sg14/auxiliary/boost.simd.h>

#define USE_FIXED_POINT

namespace mandelbrot {
    template<typename Rep, int Exponent>
    struct scalar_traits<sg14::fixed_point<Rep, Exponent>> {
        using rep = Rep;
        using scalar = sg14::fixed_point<Rep, Exponent>;
        using scalar_pack = sg14::fixed_point<boost::simd::pack<rep>, Exponent>;

        static void set(scalar_pack &fpp, int index, scalar const &value) noexcept {
            fpp.data()[index] = value.data();
        }
    };
}

namespace {
    void display(int const* subset_first, mandelbrot::vec2<int> resolution, int max_iterations) {
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

int main() {
    using namespace mandelbrot;

#if defined(USE_FIXED_POINT)
    using scalar = sg14::fixed_point<std::int64_t, -28>;
#else
    using scalar = float;
#endif

    auto const resolution = vec2<int>{80, 40};

    using coordinate = typename pack_traits<scalar>::coordinate;
    auto geometry = Geometry<scalar>{
            coordinate{-2, -2},
            coordinate{4, 0}, coordinate{0, 4},
            resolution
    };
    auto const max_iterations = 1000000;

    auto subset = generate<>(geometry, max_iterations);

    display(subset.get(), geometry.resolution, max_iterations);

    return 0;
}
