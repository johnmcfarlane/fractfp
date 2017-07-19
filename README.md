# mandelbrot

Simple Mandelbrot generator using 
[fixed_point](https://github.com/johnmcfarlane/fixed_point) and
[Boost.SIMD](https://github.com/numscale/boost.simd/) libraries.

Requires:
* cmake
* C++14 compiler
* Boost 1.60

Tested on:
* Debian 9 (w. Clang 3.8 and GCC 7.1)

Incantation:
```sh
git clone https://github.com/johnmcfarlane/mandelbrot.git
cd mandelbrot
cmake -DCMAKE_BUILD_TYPE=Release .
make
./mandelbrot
```

[feedback welcome](https://github.com/johnmcfarlane)
