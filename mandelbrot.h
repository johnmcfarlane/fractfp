#pragma once

#include <boost/simd/pack.hpp>

#include <memory>

namespace mandelbrot {
    template<typename Scalar>
    struct scalar_traits {
        using rep = Scalar;
        using scalar = Scalar;
        using scalar_pack = boost::simd::pack<scalar>;

        static void set(scalar_pack &p, int index, scalar value) noexcept {
            p[index] = value;
        }
    };

    // pack constants
    template<int N>
    constexpr auto log2() {
        static_assert(N > 0, "input must be positive");
        static_assert(!(N & (N - 1)), "input must be power of two");

        auto result = -1;
        for (auto n = N; n; ++result, n >>= 1) {
        }
        return result;
    }

    template<typename Scalar>
    struct pack_traits {
        using scalar_traits = mandelbrot::scalar_traits<Scalar>;
        using scalar = typename scalar_traits::scalar;
        using scalar_pack = typename scalar_traits::scalar_pack;
        using rep = typename scalar_traits::rep;

        using coordinate = boost::simd::pack<scalar, 2>;

        constexpr static auto pack_size = boost::simd::pack<rep>::static_size;
        constexpr static auto pack_shift = log2<pack_size>();

        // block constants
        // (a block is a rectangle of points within the Mandelbrot set, stored in a pack)
        constexpr static auto block_height_shift = pack_shift / 2;
        constexpr static auto block_width_shift = pack_shift - block_height_shift;

        constexpr static auto block_width = 1 << block_width_shift;
        constexpr static auto block_height = 1 << block_height_shift;

        static_assert(block_width * block_height == pack_size, "block does not fit snugly into pack");

        constexpr static auto block_width_mask = block_width - 1;
        constexpr static auto block_height_mask = block_height - 1;

        template<typename T>
        using pack = boost::simd::pack<T, pack_size>;

        using integer_pack = pack<int>;
    };

    ////////////////////////////////////////////////////////////////////////////////
    // types

    template<typename T>
    using vec2 = boost::simd::pack<T, 2>;

    ////////////////////////////////////////////////////////////////////////////////
    // utility fns

    constexpr auto num_steps(int const value, int const step_bits) {
        auto const step = 1 << step_bits;
        auto const mask = step - 1;
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

    template<typename Scalar>
    auto calculate(
            typename scalar_traits<Scalar>::scalar_pack const &c_x,
            typename scalar_traits<Scalar>::scalar_pack const &c_y,
            int limit) noexcept {
        using pack_traits = mandelbrot::pack_traits<Scalar>;
        typename pack_traits::integer_pack counters{0};

        auto x = c_x;
        auto y = c_y;
        auto const radius = typename pack_traits::scalar_pack{2};
        auto const radius_sq = radius*radius;
        auto confined = decltype(radius_sq <= radius_sq){true};

        while (limit--) {
            auto const xx = x * x;
            auto const yy = y * y;

            auto const currently_confined = xx + yy <= radius_sq;
            confined = currently_confined & confined;
            if (boost::simd::none(confined)) {
                break;
            }

            counters = boost::simd::if_inc(confined, counters);

            auto const x_temp = xx - yy + c_x;

            y = 2 * x * y + c_y;
            x = x_temp;
        }

        return counters;
    }

    template<typename Scalar>
    void generate_c(Geometry<Scalar> const &geometry,
                    vec2<std::unique_ptr<typename scalar_traits<Scalar>::scalar_pack[]> > &c) {
        using scalar = Scalar;
        using pack_traits = mandelbrot::pack_traits<scalar>;
        auto const block_resolution = vec2<int>{
                num_steps(geometry.resolution[0], pack_traits::block_width_shift),
                num_steps(geometry.resolution[1], pack_traits::block_height_shift)};

        // populate c
        using Vec2 = vec2<scalar>;
        auto const dx = Vec2{geometry.x_extent[0] / geometry.resolution[0],
                       geometry.x_extent[1] / geometry.resolution[1]};
        auto const dy = Vec2{geometry.y_extent[0] / geometry.resolution[0],
                       geometry.y_extent[1] / geometry.resolution[1]};

        for (auto block_row = 0; block_row != block_resolution[1]; ++block_row) {
            for (auto block_column = 0; block_column != block_resolution[0]; ++block_column) {
                auto const block_index = block_row * block_resolution[0] + block_column;
                auto &c_x_pack = c[0][block_index];
                auto &c_y_pack = c[1][block_index];

                for (auto element_index = 0, point_row = 0; point_row != pack_traits::block_height; ++point_row) {
                    auto const absolute_point_row = point_row + (block_row << pack_traits::block_height_shift);
                    for (auto point_column = 0;
                         point_column != pack_traits::block_width; ++element_index, ++point_column) {
                        auto absolute_point_column = point_column + (block_column << pack_traits::block_width_shift);

                        pack_traits::scalar_traits::set(
                                c_x_pack,
                                element_index,
                                geometry.origin[0] + absolute_point_column * dx[0] + absolute_point_row * dy[0]);
                        pack_traits::scalar_traits::set(
                                c_y_pack,
                                element_index,
                                geometry.origin[1] + absolute_point_column * dx[1] + absolute_point_row * dy[1]);
                    }
                }
            }
        }
    }

    template<typename Scalar>
    auto extract_results(
            std::unique_ptr<typename pack_traits<Scalar>::integer_pack[]> &counters,
            vec2<int> const block_resolution, vec2<int> const resolution) {
        using pack_traits = mandelbrot::pack_traits<Scalar>;
        auto const num_points = static_cast<std::size_t>(resolution[0]) * resolution[1];
        auto subset = std::unique_ptr<int[]>(new int[num_points]);
        for (auto point_row = 0; point_row != resolution[1]; ++point_row) {
            auto block_row = point_row >> pack_traits::block_height_shift;

            for (auto point_column = 0; point_column != resolution[0]; ++point_column) {
                auto block_column = point_column >> pack_traits::block_width_shift;

                auto point_index = point_row * resolution[0] + point_column;
                auto block_index = block_row * block_resolution[0] + block_column;

                auto point_relative_row = point_row & pack_traits::block_height_mask;
                auto point_relative_column = point_column & pack_traits::block_width_mask;
                auto point_relative_index = point_relative_row * pack_traits::block_width + point_relative_column;

                subset[point_index] = counters[block_index][point_relative_index];
            }
        }

        return subset;
    }

    template<typename Scalar>
    auto generate(Geometry<Scalar> const &geometry, int const limit) {
        using pack_traits = mandelbrot::pack_traits<Scalar>;
        using scalar_pack = typename pack_traits::scalar_pack;
        using integer_pack = typename pack_traits::integer_pack;

        // dimensions
        auto const block_resolution = vec2<int> {
                num_steps(geometry.resolution[0], pack_traits::block_width_shift),
                num_steps(geometry.resolution[1], pack_traits::block_height_shift)};
        auto const num_blocks = block_resolution[0] * block_resolution[1];

        // for each point for each block
        vec2<std::unique_ptr<scalar_pack[]>> c;

        // allocate c
        c[0].reset(new scalar_pack[num_blocks]);
        c[1].reset(new scalar_pack[num_blocks]);

        // calculate c
        generate_c<Scalar>(geometry, c);

        // generate set
        auto counters = std::unique_ptr<integer_pack[]>{new integer_pack[num_blocks]};
        for (auto i = 0; i != num_blocks; ++i) {
            counters[i] = calculate<Scalar>(c[0][i], c[1][i], limit);
        }

        // return iterations
        return extract_results<Scalar>(counters, block_resolution, geometry.resolution);
    }
}
