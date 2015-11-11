/**
 * GPU COMPUTING 2014
 *
 * Device-based highly accurate timer (resolution approx 0.5 microseconds)
 * that also automatically synchronises
 *
 * Usage examples: for instance in 05-vecadd-benchmark.cu
 *
 * Explanation of the underlying CUDA Event API: in the lecture on synchronisation
 *
 * (c) 2014 dominik.goeddeke@math.tu-dortmund.de
 */

#ifndef GPUTIMER_GUARD
#define GPUTIMER_GUARD


// includes
#include <iostream>
#include <cstdlib>
#include <cuda_runtime.h>


class GPUtimer
{
  private:

  // stream that this timer works in
  cudaStream_t _stream;

  // start and stop events
  cudaEvent_t _start, _stop;

  public:

  // CTOR, uses default stream 0
  GPUtimer () : _stream(0)
  {
    cudaError_t status = cudaEventCreate(&_start);
    if (status != cudaSuccess)
    {
      std::cerr << "CUDA Error (in GPUtimer CTOR): " << cudaGetErrorString(status) << std::endl;
      exit(1);
    }
    status = cudaEventCreate(&_stop);
    if (status != cudaSuccess)
    {
      std::cerr << "CUDA Error (in GPUtimer CTOR): " << cudaGetErrorString(status) << std::endl;
      exit(1);
    }
   // std::cout << "GPUtimer initialised to stream " << _stream << std::endl;
  }

  // starts the GPU timer
  void start ()
  {
    cudaEventRecord(_start, _stream);
  }

  // stops the timer and returns the elapsed time in seconds
  double stop ()
  {
    cudaEventRecord(_stop, _stream);
    cudaEventSynchronize(_stop);
    float elapsed;
    cudaEventElapsedTime(&elapsed, _start, _stop);
    return elapsed/1000.0;
  }

};
  
#endif // GPUTIMER_GUARD

