NRM Benchmarks: Benchmarks for progress mesurement
==================================================

These benchmarks are designed to offer various node-level workload patterns,
based on basic computational science kernels, along with application-side
reports on how fast the applications are progressing towards their figure of
merit.

## Requirements:

* autoconf
* automake
* pkg-config
* argo nrm downstream libraries

## Installation

```
sh autogen.sh
./configure
make -j install
```

# General Architecture

All the benchmarks rely on OpenMP to distribute their workload across a node
topology. As much as possible, these workloads should be designed to reach the
anticipated bottleneck.

All benchmarks should report their configuration at the end of the program, as
well as basic performance metrics (total runtime, figure-of-merit). This output
should allow one to verify that the benchmark is stressing the right part of
the system.

The `src` directory is organized to into subdirectories based on the type of
reporting that the benchmark does, and other similar criterias. The name of
each benchmark reflects all its categories.
