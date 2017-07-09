#include <memory>

#include <sg14/auxiliary/boost.simd.h>

#define USE_FIXED_POINT

#if defined(USE_FIXED_POINT)
using Rep = std::int64_t;
constexpr auto Exponent = -28;
using Scalar = sg14::fixed_point<Rep, Exponent>;
constexpr auto pack_size = boost::simd::pack<Rep>::static_size;
#else
using Scalar = float;
constexpr auto pack_size = boost::simd::pack<Scalar>::static_size;
#endif

namespace {
    namespace mandelbrot {
        // pack constants
        template<int N>
        constexpr auto log2() {
            static_assert(N > 0, "input must be positive");
            static_assert(!(N&(N-1)), "input must be power of two");

            auto result = -1;
            for (auto n = N; n; ++result, n >>= 1) {
            }
            return result;
        }

        constexpr auto pack_shift = log2<pack_size>();

        // block constants
        // (a block is a rectangle of points within the Mandelbrot set, stored in a pack)
        constexpr auto block_width_shift = pack_shift / 2;
        constexpr auto block_height_shift = pack_shift - block_width_shift;

        constexpr auto block_width = 1 << block_width_shift;
        constexpr auto block_height = 1 << block_height_shift;

        static_assert(block_width * block_height == pack_size, "block does not fit snugly into pack");

        constexpr auto block_width_mask = block_width - 1;
        constexpr auto block_height_mask = block_height - 1;

        ////////////////////////////////////////////////////////////////////////////////
        // types

        template<typename T>
        using pack = boost::simd::pack<T, pack_size>;

        template<typename T>
        using vec2 = boost::simd::pack<T, 2>;

        using sg14::fixed_point;

        template<typename Rep, int Exponent = 0>
        using fixed_point_pack = fixed_point<pack<Rep>, Exponent>;

        ////////////////////////////////////////////////////////////////////////////////
        // utility fns

        template<typename Scalar>
        void set(pack<Scalar> &p, int index, Scalar value) noexcept {
            p[index] = value;
        }

        template<typename PackRep, int PackExponent, typename ElementRep, int ElementExponent>
        void set(fixed_point<pack<PackRep>, PackExponent> &fpp, int index,
                 fixed_point<ElementRep, ElementExponent> const &value) noexcept {
            set(fpp.data(), index, value.data());
        }

        constexpr auto num_steps(int value, int step_bits) {
            auto step = 1 << step_bits;
            auto mask = step - 1;
            return (value + mask) >> step_bits;
        }

        ////////////////////////////////////////////////////////////////////////////////
        // set generation functions

        template<typename Scalar>
        struct Geometry {
            vec2<Scalar> origin;
            vec2<Scalar> x_extent;
            vec2<Scalar> y_extent;
            vec2<int> resolution;
        };

        template<typename ScalarPack>
        pack<int> calculate(ScalarPack const &c_x, ScalarPack const &c_y, int limit) noexcept {
            pack<int> counters{0};

            auto x = c_x;
            auto y = c_y;
            auto confined = decltype(x < y){true};
            auto radius = ScalarPack{4};

            while (limit--) {
                auto xx = x * x;
                auto yy = y * y;

                auto currently_confined = xx + yy <= radius;
                confined = currently_confined & confined;
                if (boost::simd::none(confined)) {
                    break;
                }
                counters = boost::simd::if_inc(confined, counters);

                auto x_temp = xx - yy + c_x;

                y = 2 * x * y + c_y;
                x = x_temp;
            }

            return counters;
        }

        template<typename Scalar, typename ScalarPack>
        void generate_c(Geometry<Scalar> const &geometry,
                        vec2<std::unique_ptr<ScalarPack[]>> &c) {
            auto block_resolution = vec2<int>{
                    num_steps(geometry.resolution[0], block_width_shift),
                    num_steps(geometry.resolution[1], block_height_shift)};

            // populate c
            using Vec2 = vec2<Scalar>;
            auto dx = Vec2{geometry.x_extent[0] / geometry.resolution[0],
                           geometry.x_extent[1] / geometry.resolution[1]};
            auto dy = Vec2{geometry.y_extent[0] / geometry.resolution[0],
                           geometry.y_extent[1] / geometry.resolution[1]};

            for (auto block_row = 0; block_row != block_resolution[1]; ++block_row) {
                for (auto block_column = 0; block_column != block_resolution[0]; ++block_column) {
                    auto block_index = block_row * block_resolution[0] + block_column;
                    auto &c_x_pack = c[0][block_index];
                    auto &c_y_pack = c[1][block_index];

                    for (auto element_index = 0, point_row = 0; point_row != block_height; ++point_row) {
                        auto absolute_point_row = point_row + (block_row << block_height_shift);
                        for (auto point_column = 0; point_column != block_width; ++element_index, ++point_column) {
                            auto absolute_point_column = point_column + (block_column << block_width_shift);

                            set(c_x_pack,
                                element_index,
                                geometry.origin[0] + absolute_point_column * dx[0] + absolute_point_row * dy[0]);
                            set(c_y_pack,
                                element_index,
                                geometry.origin[1] + absolute_point_column * dx[1] + absolute_point_row * dy[1]);
                        }
                    }
                }
            }
        }

        auto extract_results(
                std::unique_ptr<pack<int>[]> &counters,
                vec2<int> block_resolution, vec2<int> const resolution) {
            auto num_points = static_cast<std::size_t>(resolution[0]) * resolution[1];
            auto subset = std::unique_ptr<int[]>(new int[num_points]);
            for (auto point_row = 0; point_row != resolution[1]; ++point_row) {
                auto block_row = point_row >> block_height_shift;

                for (auto point_column = 0; point_column != resolution[0]; ++point_column) {
                    auto block_column = point_column >> block_width_shift;

                    auto point_index = point_row * resolution[0] + point_column;
                    auto block_index = block_row * block_resolution[0] + block_column;

                    auto point_relative_row = point_row & block_height_mask;
                    auto point_relative_column = point_column & block_width_mask;
                    auto point_relative_index = point_relative_row * block_width + point_relative_column;

                    subset[point_index] = counters[block_index][point_relative_index];
                }
            }

            return subset;
        }

        template<typename Scalar, typename ScalarPack>
        auto generate(Geometry<Scalar> const &geometry, int const limit) {
            // dimensions
            auto block_resolution = vec2<int> {
                    num_steps(geometry.resolution[0], block_width_shift),
                    num_steps(geometry.resolution[1], block_height_shift)};
            auto num_blocks = block_resolution[0] * block_resolution[1];

            // for each point for each block
            vec2<std::unique_ptr<ScalarPack[]>> c;

            // allocate c
            c[0].reset(new ScalarPack[num_blocks]);
            c[1].reset(new ScalarPack[num_blocks]);

            // calculate c
            generate_c<Scalar, ScalarPack>(geometry, c);

            // generate set
            auto counters = std::unique_ptr<pack<int>[]>{new pack<int>[num_blocks]};
            for (auto i = 0; i != num_blocks; ++i) {
                counters[i] = calculate<ScalarPack>(c[0][i], c[1][i], limit);
            }

            // return iterations
            return extract_results(counters, block_resolution, geometry.resolution);
        }
    }
}

int main() {

#if defined(USE_FIXED_POINT)
    using ScalarPack = mandelbrot::fixed_point_pack<Rep, Exponent>;
#else
    using ScalarPack = boost::simd::pack<Scalar>;
#endif

    using Coordinate = mandelbrot::vec2<Scalar>;
    using Resolution = mandelbrot::vec2<int>;

    auto const resolution = Resolution{80, 40};
    auto const max_iterations = 1000000;

    auto subset = mandelbrot::generate<Scalar, ScalarPack>(
            mandelbrot::Geometry < Scalar > {
                    Coordinate{-2, -2},
                    Coordinate{4, 0}, Coordinate{0, 4},
                    resolution},
            max_iterations);

    auto const max_displayed_iterations = std::min(max_iterations, 96);
    for (auto y = 0; y != resolution[1]; ++y) {
        auto line_index = y * resolution[0];
        for (auto x = 0; x != resolution[0]; ++x) {
            auto index = line_index + x;
            auto iterations = subset[index];
            std::putchar((iterations >= max_displayed_iterations) ? ' ' : ' ' + iterations);
        }
        std::putchar('\n');
    }

    return 0;
}
