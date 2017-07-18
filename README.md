# mandelbrot
Simple Mandelbrot generator using 
[fixed_point](https://github.com/johnmcfarlane/fixed_point) and
[Boost.SIMD](https://github.com/numscale/boost.simd/) libraries.

```sh
git clone git@github.com:johnmcfarlane/mandelbrot.git
cd mandelbrot
cmake -DCMAKE_BUILD_TYPE=Release .
make
./mandelbrot
```